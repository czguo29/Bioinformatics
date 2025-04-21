use std::collections::HashMap;

/// A simple HMM with two states: MOTIF and NONMOTIF. Each state has transition probabilities 
/// to itself and to the other state. The emission model is a probability distribution over nucleotides.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum HMMState {
    Motif,
    NonMotif,
}

/// Defines an HMM for a basic motif-finding scenario:
/// transition_probs[s1][s2] = P(S_{i} = s2 | S_{i-1} = s1)
/// emission_probs[state][nt] = P(nt | state)
struct MotifHMM {
    transition_probs: HashMap<(HMMState, HMMState), f64>,
    emission_probs: HashMap<(HMMState, char), f64>,
}

impl MotifHMM {
    /// Initialize a minimal HMM with user-defined transition and emission probabilities.
    fn new(
        transition_probs: HashMap<(HMMState, HMMState), f64>,
        emission_probs: HashMap<(HMMState, char), f64>,
    ) -> Self {
        MotifHMM {
            transition_probs,
            emission_probs,
        }
    }

    /// Runs the Viterbi algorithm on a DNA sequence. Returns the most-likely path of states.
    fn viterbi(&self, seq: &str) -> Vec<HMMState> {
        let seq_chars: Vec<char> = seq.chars().collect();
        let n = seq_chars.len();

        // States: [Motif, NonMotif]
        let states = [HMMState::Motif, HMMState::NonMotif];

        // dp[i][st] = maximum log probability of any path that ends in state st at position i
        let mut dp = vec![HashMap::new(); n];
        // backpointer[i][st] = the state that yields the best path to st at position i
        let mut backpointer = vec![HashMap::new(); n];

        // Initialization: i=0
        for &st in &states {
            // Define initial_prob as f64 to avoid ambiguity
            let initial_prob: f64 = 0.5;
            let emit = *self.emission_probs.get(&(st, seq_chars[0])).unwrap_or(&1e-9);
            // Remove unnecessary parentheses
            dp[0].insert(st, initial_prob.ln() + emit.ln());
        }

        // Recurrence
        for i in 1..n {
            for &st_curr in &states {
                let emit_prob = *self
                    .emission_probs
                    .get(&(st_curr, seq_chars[i]))
                    .unwrap_or(&1e-9);
                let mut best_score = f64::NEG_INFINITY;
                let mut best_prev = st_curr;

                for &st_prev in &states {
                    let prev_score = dp[i - 1].get(&st_prev).unwrap_or(&f64::NEG_INFINITY);
                    let trans = *self
                        .transition_probs
                        .get(&(st_prev, st_curr))
                        .unwrap_or(&1e-9);
                    let candidate = prev_score + trans.ln() + emit_prob.ln();

                    if candidate > best_score {
                        best_score = candidate;
                        best_prev = st_prev;
                    }
                }
                dp[i].insert(st_curr, best_score);
                backpointer[i].insert(st_curr, best_prev);
            }
        }

        // Termination: pick the best final state
        let mut best_final_state = HMMState::NonMotif;
        let mut best_final_score = f64::NEG_INFINITY;
        for &st in &states {
            let score = dp[n - 1].get(&st).unwrap_or(&f64::NEG_INFINITY);
            if *score > best_final_score {
                best_final_score = *score;
                best_final_state = st;
            }
        }

        // Traceback
        let mut path = vec![best_final_state; n];
        let mut current_state = best_final_state;
        for i in (1..n).rev() {
            current_state = backpointer[i][&current_state];
            path[i - 1] = current_state;
        }
        path
    }
}

fn main() {
    // Example transition probabilities
    let mut t = HashMap::new();
    t.insert((HMMState::Motif, HMMState::Motif), 0.8);
    t.insert((HMMState::Motif, HMMState::NonMotif), 0.2);
    t.insert((HMMState::NonMotif, HMMState::Motif), 0.05);
    t.insert((HMMState::NonMotif, HMMState::NonMotif), 0.95);

    // Example emission probabilities
    let mut e = HashMap::new();
    for &c in &['A','C','G','T'] {
        e.insert((HMMState::NonMotif, c), 0.25);
    }
    e.insert((HMMState::Motif, 'A'), 0.4);
    e.insert((HMMState::Motif, 'C'), 0.2);
    e.insert((HMMState::Motif, 'G'), 0.2);
    e.insert((HMMState::Motif, 'T'), 0.2);

    let hmm = MotifHMM::new(t, e);

    // Sample DNA sequence
    let dna_sequence = "ACGGAATACACGG";

    let path = hmm.viterbi(dna_sequence);
    println!("Sequence:     {}", dna_sequence);
    println!("Most-likely:  {:?}", path);
}