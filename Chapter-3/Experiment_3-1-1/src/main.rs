use rayon::prelude::*; // Provides parallel iterators
use rand::prelude::*;  // Provides random utilities (trait imports)
use rand::rng;         // Newer function instead of the deprecated `thread_rng`
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::sync::atomic::{AtomicBool, Ordering};

/// A parallel & thread-safe Bloom Filter.
///
/// # Fields
/// - `bits`: A vector of atomic booleans used to represent bits in the filter.
/// - `k`: The number of hash functions (seed variations).
/// - `m`: The total size of the bit array (number of bits).
pub struct BloomFilter {
    /// Vector of AtomicBool to allow concurrent writes without locks.
    bits: Vec<AtomicBool>,
    /// Number of hash functions to use.
    k: usize,
    /// Size of the bit array.
    m: usize,
}

impl BloomFilter {
    /// Creates a new Bloom Filter with `m` bits and `k` hash functions.
    ///
    /// # Arguments
    ///
    /// * `m` - The number of bits in the filter.
    /// * `k` - The number of distinct hash functions to use.
    ///
    /// # Returns
    ///
    /// * `BloomFilter` - An instance of the Bloom Filter.
    pub fn new(m: usize, k: usize) -> Self {
        // Initialize `m` atomic booleans to `false`.
        let bits = (0..m).map(|_| AtomicBool::new(false)).collect();
        
        Self { bits, k, m }
    }

    /// Inserts an item into the Bloom Filter by setting the relevant bits to `true`.
    ///
    /// # Type Parameters
    /// * `T` - A type that can be hashed. Must also be `Sync` for parallel iteration.
    ///
    /// # Arguments
    ///
    /// * `item` - A reference to the item to insert.
    ///
    /// This method uses parallel iteration over the `k` hash seeds.
    pub fn insert<T: Hash + Sync>(&self, item: &T) {
        // For each of the k hash functions, compute a position and set that bit to `true`.
        (0..self.k).into_par_iter().for_each(|seed| {
            let pos = self.hash_with_seed(item, seed) % self.m;
            // Store `true` in the bit array at index `pos`.
            self.bits[pos].store(true, Ordering::Relaxed);
        });
    }

    /// Checks whether an item may be in the Bloom Filter.
    ///
    /// # Type Parameters
    /// * `T` - A type that can be hashed. Must also be `Sync` for parallel iteration.
    ///
    /// # Arguments
    /// * `item` - A reference to the item to check.
    ///
    /// # Returns
    /// * `bool` - `true` if the item *might* be in the filter, `false` if definitely not.
    ///
    /// This method uses parallel iteration over the k hash seeds. 
    /// If *any* required bit is `false`, we immediately know the item isn't in the filter.
    pub fn contains<T: Hash + Sync>(&self, item: &T) -> bool {
        // Check all k positions; if any are `false`, return false.
        (0..self.k).into_par_iter().all(|seed| {
            let pos = self.hash_with_seed(item, seed) % self.m;
            self.bits[pos].load(Ordering::Relaxed)
        })
    }

    /// Computes a hash for an item combined with a seed.
    ///
    /// # Type Parameters
    /// * `T` - A type that can be hashed.
    ///
    /// # Arguments
    /// * `item` - The item to be hashed.
    /// * `seed` - The seed value to differentiate the hash functions.
    ///
    /// # Returns
    /// * `usize` - A hashed value used to index into the bit vector.
    #[inline]
    fn hash_with_seed<T: Hash>(&self, item: &T, seed: usize) -> usize {
        let mut hasher = DefaultHasher::new();
        item.hash(&mut hasher); // Hash the item
        seed.hash(&mut hasher); // Mix in the seed
        hasher.finish() as usize
    }
}

/// Generates a random DNA sequence of specified length using the bases A, C, G, and T.
///
/// # Arguments
///
/// * `len` - The length of the DNA sequence to generate.
///
/// # Returns
///
/// * `String` - A random DNA sequence (e.g. "ACGTG...").
fn generate_random_dna(len: usize) -> String {
    // Use the new `rng()` function for thread-local RNG.
    let mut r = rng();
    let bases = ['A', 'C', 'G', 'T'];

    (0..len)
        // For each position, pick one of the 4 bases at random.
        .map(|_| *bases.choose(&mut r).unwrap())
        .collect()
}

/// Generates a vector of random DNA sequences.
///
/// # Arguments
///
/// * `num_sequences` - The number of DNA sequences to generate.
/// * `seq_length` - The length of each DNA sequence.
///
/// # Returns
///
/// * `Vec<String>` - A vector containing the generated DNA sequences.
fn generate_synthetic_genomic_data(num_sequences: usize, seq_length: usize) -> Vec<String> {
    (0..num_sequences)
        .map(|_| generate_random_dna(seq_length))
        .collect()
}

fn main() {
    // ======= Bloom Filter Parameters =======
    // `bloom_size` (m): total number of bits in the Bloom Filter.
    // `num_hashes` (k): number of distinct hash functions to use.
    let bloom_size = 50_000;
    let num_hashes = 5;

    // ======= Create the Bloom Filter =======
    let bloom = BloomFilter::new(bloom_size, num_hashes);

    // ======= Generate Synthetic Genomic Data =======
    // For demonstration, we'll generate 5,000 DNA sequences, each 20 bases long.
    let num_sequences = 5_000;
    let seq_length = 20;
    let dataset = generate_synthetic_genomic_data(num_sequences, seq_length);

    // ======= Insert Sequences into the Bloom Filter =======
    // We use a parallel iterator to insert each sequence in parallel.
    dataset.par_iter().for_each(|seq| {
        bloom.insert(seq);
    });

    // ======= Test Membership =======
    // We'll test membership for:
    // 1) Two sequences we know are in the dataset (should return true, unless an unexpected false negative occurs, which is theoretically impossible with a proper Bloom Filter).
    // 2) Two random sequences likely not in the dataset (may still return true if there's a false positive).
    let test_sequences = vec![
        dataset[0].clone(),      // definitely in
        dataset[1].clone(),      // definitely in
        generate_random_dna(seq_length), // likely not in
        generate_random_dna(seq_length), // likely not in
    ];

    // Print out test results
    for seq in &test_sequences {
        let result = bloom.contains(seq);
        println!("Sequence: {} => in Bloom Filter? {}", seq, result);
    }
}