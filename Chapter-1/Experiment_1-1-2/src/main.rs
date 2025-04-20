use bio::io::fasta;
use std::error::Error;

/// Demonstrates Rust’s zero-cost iterator abstractions by reading a FASTA file and
/// calculating the total GC content of sequences that pass a certain length filter.
/// Despite using map/filter/sum in a chained functional style, the compiler optimizes
/// these calls to efficient loops under the hood.
fn main() -> Result<(), Box<dyn Error>> {
    // 1) Open the FASTA file with the `bio::io::fasta` crate
    let reader = fasta::Reader::from_file(r"C:\Users\Christyane Zabdi\Downloads\CZ-UI S2\Sem 2\BGVR\Chapter 1\Experiment_1-1-2\src\example.fasta")?;

    // 2) Chain iterator adapters to parse, filter, and process sequences
    //    - `map` extracts the sequence and converts it to String
    //    - `filter` excludes sequences shorter than 50 nucleotides
    //    - another `map` computes the GC count
    //    - `sum()` accumulates them into a single integer
    let total_gc: usize = reader
        .records()
        .map(|rec_res| {
            // Convert each sequence from bytes to a String
            let record = rec_res.unwrap();
            String::from_utf8_lossy(record.seq()).to_string()
        })
        .filter(|seq| seq.len() >= 50)
        .map(|seq| {
            // Compute GC content for each filtered sequence
            seq.chars()
                .filter(|&c| c == 'G' || c == 'C')
                .count()
        })
        .sum();

    // 3) Print the result
    println!("Total GC content in sequences >= 50 nt: {}", total_gc);

    Ok(())
}