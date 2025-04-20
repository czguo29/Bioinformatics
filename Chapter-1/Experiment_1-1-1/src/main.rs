use bio::io::fasta;
use rayon::prelude::*;
use std::error::Error;

/// Reads a FASTA file and returns a vector of sequences (owned `Vec<String>`).
/// No references to local data are returned, sidestepping E0515 errors.
fn load_sequences() -> Result<Vec<String>, Box<dyn Error>> {
    let reader = fasta::Reader::from_file(r"C:\Users\Christyane Zabdi\Downloads\CZ-UI S2\Sem 2\BGVR\Chapter 1\Experiment_1-1-1\src\example.fasta")?;
    let mut sequences = Vec::new();

    // Read each record and push into a local vector (no Mutex needed if single-thread read)
    for record_result in reader.records() {
        let record = record_result?;
        let seq_str = String::from_utf8_lossy(record.seq()).to_string();
        sequences.push(seq_str);
    }

    Ok(sequences)
}

fn main() -> Result<(), Box<dyn Error>> {
    // 1) Load all sequences from the FASTA file into an owned vector
    let seqs = load_sequences()?;

    // 2) Process them in parallel with Rayon (if needed)
    seqs.par_iter().for_each(|seq| {
        let gc_count = seq.chars().filter(|&c| c == 'G' || c == 'C').count();
        println!("Length: {}, GC: {}", seq.len(), gc_count);
    });

    println!("Successfully processed {} sequences.", seqs.len());
    Ok(())
}