use rayon::prelude::*;
use rand::prelude::*;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use rand::rng; 
/// Represents the resulting MinHash signature of a dataset.
#[derive(Debug)]
pub struct MinHashSignature {
    pub signature: Vec<u64>,
}
/// A MinHasher computes MinHash signatures for sets of items.
///
/// # Fields
/// * `num_hashes` - The number of distinct hash functions to apply (often called "bands" in MinHash).
pub struct MinHasher {
    pub num_hashes: usize,
}
impl MinHasher {
    /// Creates a new `MinHasher`.
    ///
    /// # Arguments
    /// * `num_hashes` - Number of hash functions (seeds) to use for MinHash.
    pub fn new(num_hashes: usize) -> Self {
        MinHasher { num_hashes }
    }
    /// Computes the MinHash signature for a collection of hashable items.
    ///
    /// This version parallelizes over each hash function and iterates over items in parallel,
    /// taking the minimum hash value for each seed.
    ///
    /// # Arguments
    /// * `items` - A slice of items (e.g., strings) that implement `Hash`.
    ///
    /// # Returns
    /// * `MinHashSignature` containing the minimum hash value for each seed.
    pub fn compute_signature<T: Hash + Sync>(&self, items: &[T]) -> MinHashSignature {
        // For each seed in [0..num_hashes], compute the min hash among all items in parallel.
        let signature: Vec<u64> = (0..self.num_hashes)
            .into_par_iter()
            .map(|seed| {
                // For each item, compute hash_with_seed and take the minimum over all items.
                items
                    .par_iter()
                    .map(|item| self.hash_with_seed(item, seed))
                    .min()
                    .unwrap_or(u64::MAX)
            })
            .collect();
        MinHashSignature { signature }
    }
    /// Computes the approximate Jaccard similarity between two MinHash signatures.
    ///
    /// # Arguments
    /// * `sig_a` - First MinHash signature.
    /// * `sig_b` - Second MinHash signature.
    ///
    /// # Returns
    /// * A floating-point value representing the fraction of matching hash values
    ///   across all seeds (0.0 to 1.0).
    pub fn similarity(&self, sig_a: &MinHashSignature, sig_b: &MinHashSignature) -> f64 {
        // Count matches in parallel, then sum them up.
        let matches = (0..self.num_hashes)
            .into_par_iter()
            .map(|i| if sig_a.signature[i] == sig_b.signature[i] { 1 } else { 0 })
            .sum::<usize>();
        matches as f64 / self.num_hashes as f64
    }
    /// Internal helper to hash an item with a specified seed. We mix both the item's hash
    /// and the seed into a `DefaultHasher`.
    fn hash_with_seed<T: Hash>(&self, item: &T, seed: usize) -> u64 {
        let mut hasher = DefaultHasher::new();
        item.hash(&mut hasher);
        seed.hash(&mut hasher);
        hasher.finish()
    }
}
/// Generates a random DNA sequence (consisting of A, C, G, T) of a specified length.
///
/// # Arguments
/// * `len` - Length of the DNA sequence.
///
/// # Returns
/// * A `String` representing a random DNA sequence.
fn generate_random_dna_sequence(len: usize) -> String {
    let mut rng = rng();
    let bases = ['A', 'C', 'G', 'T'];
    (0..len)
        .map(|_| *bases.choose(&mut rng).unwrap())
        .collect()
}
/// Generates `num_sequences` random DNA strings, each of length `seq_len`.
fn generate_synthetic_genomic_data(num_sequences: usize, seq_len: usize) -> Vec<String> {
    (0..num_sequences)
        .map(|_| generate_random_dna_sequence(seq_len))
        .collect()
}
fn main() {
    // Create two synthetic "genomic" datasets.
    // In practice, these could be real DNA reads or gene fragments.
    let set_size = 5_000;
    let seq_length = 20;
    let set1 = generate_synthetic_genomic_data(set_size, seq_length);
    let set2 = generate_synthetic_genomic_data(set_size, seq_length);
    // Create a MinHasher with 100 hash functions (seeds).
    let hasher = MinHasher::new(100);
    // Compute MinHash signatures in parallel.
    let sig1 = hasher.compute_signature(&set1);
    let sig2 = hasher.compute_signature(&set2);
    // Print the first few signature values just for illustration.
    println!("MinHash signature (first few) for set1: {:?}", &sig1.signature[..5]);
    println!("MinHash signature (first few) for set2: {:?}", &sig2.signature[..5]);
    // Compute approximate Jaccard similarity.
    let similarity = hasher.similarity(&sig1, &sig2);
    println!("Approximate Jaccard similarity = {:.3}", similarity);
}
