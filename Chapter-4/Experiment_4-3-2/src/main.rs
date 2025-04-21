use rand::{rng, Rng}; // Use the newer rng() if that is what's required in your environment.
// Removed `std::collections::HashMap` because it's not used.

/// A basic struct representing a motif with a position weight matrix (PWM).
/// For simplicity, we store probabilities in a 2D vector, motif_length x 4.
#[derive(Debug)]
struct MotifModel {
    pwm: Vec<[f64; 4]>, // e.g., [p(A), p(C), p(G), p(T)] for each position
}

impl MotifModel {
    pub fn new_random(motif_length: usize) -> Self {
        // Initialize each position with random probabilities
        let mut rng = rng(); // replaced `rand::thread_rng()` with `rng()`
        let mut pwm = vec![[0.25; 4]; motif_length];
        for i in 0..motif_length {
            let mut sum = 0.0;
            for nt in 0..4 {
                // replaced `rng.gen_range(0.0..1.0)` with `rng.random_range(0.0..1.0)`
                let r = rng.random_range(0.0..1.0);
                pwm[i][nt] = r;
                sum += r;
            }
            // Normalize each row so that p(A)+p(C)+p(G)+p(T) = 1
            for nt in 0..4 {
                pwm[i][nt] /= sum;
            }
        }
        MotifModel { pwm }
    }

    /// Returns the probability of a k-mer under the current PWM model.
    pub fn score_kmer(&self, kmer: &[u8]) -> f64 {
        let mut prob = 1.0;
        for (i, &nt_byte) in kmer.iter().enumerate() {
            let col_index = match nt_byte {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 4, // invalid
            };
            if col_index < 4 {
                prob *= self.pwm[i][col_index];
            } else {
                prob *= 1e-9;
            }
        }
        prob
    }

    /// M-step: updates the PWM entries from the partial count matrix
    pub fn update_from_counts(&mut self, partial_counts: &Vec<[f64; 4]>) {
        let motif_length = self.pwm.len();
        for i in 0..motif_length {
            let mut sum = partial_counts[i].iter().sum::<f64>();
            if sum < 1e-9 {
                sum = 1.0;
            }
            for nt in 0..4 {
                self.pwm[i][nt] = partial_counts[i][nt] / sum;
            }
        }
    }
}

/// A minimal function simulating one E-step + partial M-step in HPC context, returning local counts.
fn em_iteration(model: &MotifModel, sequences: &[Vec<u8>]) -> Vec<[f64; 4]> {
    let motif_length = model.pwm.len();

    // We accumulate partial counts for each position x nucleotide
    let mut local_counts = vec![[0.0; 4]; motif_length];

    // For each read, we attempt to align it in a simplified manner (slide the motif and pick best position).
    // We'll treat the best position as fully “responsible” for this read.
    for seq in sequences {
        if seq.len() >= motif_length {
            let mut best_prob = -1.0;
            let mut best_start = 0;
            // Slide motif along the read, pick position with highest PWM probability
            for start in 0..=(seq.len() - motif_length) {
                let kmer = &seq[start..start+motif_length];
                let prob = model.score_kmer(kmer);
                if prob > best_prob {
                    best_prob = prob;
                    best_start = start;
                }
            }
            // Update partial counts for that best alignment
            let best_kmer = &seq[best_start..best_start+motif_length];
            for (i, &nt_byte) in best_kmer.iter().enumerate() {
                let col_index = match nt_byte {
                    b'A' | b'a' => 0,
                    b'C' | b'c' => 1,
                    b'G' | b'g' => 2,
                    b'T' | b't' => 3,
                    _ => continue,
                };
                local_counts[i][col_index] += 1.0;
            }
        }
    }

    local_counts
}

fn main() {
    // Suppose HPC ephemeral containers parse different data segments. For demonstration:
    let sequences = vec![
        b"ACGTACGTAC".to_vec(),
        b"GTTACACAC".to_vec(),
        b"ATATGGGG".to_vec(),
    ];

    // Initialize random motif model of length 4
    let mut model = MotifModel::new_random(4);

    println!("Initial motif model: {:?}", model.pwm);

    // E-step + partial M-step in local HPC container
    let partial_counts = em_iteration(&model, &sequences);

    // In HPC, partial_counts would be merged with other containers via sum. 
    model.update_from_counts(&partial_counts);

    println!("Updated motif model: {:?}", model.pwm);
}