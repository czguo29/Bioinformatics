use bio::io::fasta;
use rayon::prelude::*;
use std::error::Error;

/// A simple, naive function to count how many times `motif` appears in `seq`.
/// Overlapping occurrences are allowed by advancing the search one character at a time.
fn count_occurrences(seq: &str, motif: &str) -> usize {
    let mut count = 0;
    let mut start = 0;

    while let Some(pos) = seq[start..].find(motif) {
        count += 1;
        // Move ahead by one to find potential overlapping matches.
        start += pos + 1;
    }

    count
}

fn main() -> Result<(), Box<dyn Error>> {
    // 1) Specify the motif we want to search for
    let motif = "GATTACA";

    // 2) Open the FASTA file (assume 'example.fasta' for demonstration).
    //    In a real pipeline, you might parse this from command-line arguments.
    let reader = fasta::Reader::from_file(r"C:\Users\Christyane Zabdi\Downloads\CZ-UI S2\Sem 2\BGVR\Chapter 1\Experiment_1-1-3\src\example.fasta")?;

    // 3) Collect all sequences into a vector of owned Strings.
    //    This is done sequentially, but is typically fast enough for moderate files.
    //    For massive data, chunk-based or streaming approaches might be required.
    let sequences: Vec<String> = reader
        .records()
        .filter_map(|rec_res| {
            match rec_res {
                Ok(rec) => Some(String::from_utf8_lossy(rec.seq()).to_string()),
                Err(_) => None  // skip or handle malformed records
            }
        })
        .collect();

    // 4) Use Rayon’s parallel iterator to split the workload across available threads.
    //    Each sequence is processed by `count_occurrences`, tallying how many times
    //    the motif appears. The final `.sum()` aggregates the total matches.
    let total_matches: usize = sequences
        .par_iter()
        .map(|seq| count_occurrences(seq, motif))
        .sum();

    // 5) Print the result
    println!("Total occurrences of motif '{}' across all sequences: {}", motif, total_matches);

    Ok(())
}
