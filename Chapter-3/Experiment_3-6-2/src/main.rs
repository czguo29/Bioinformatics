use serde::{Serialize, Deserialize};
use std::fs;

/// PartialFMIndex is a placeholder for partial data from an FM-index
#[derive(Debug, Serialize, Deserialize)]
struct PartialFMIndex {
    reference_chunk: String,
    bwt_data: Vec<u8>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Suppose we chunk the reference by a user-defined size
    let ref_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| r"C:\Users\Christyane Zabdi\Downloads\CZ-UI S2\Sem 2\BGVR\Chapter 3\Experiment_3-6-2\src\reference.fasta".to_string());
    let chunk_size: usize = std::env::args()
        .nth(2)
        .unwrap_or_else(|| "1000000".to_string())
        .parse()?;

    // Load the file contents, chunk them, and build partial FM-index placeholders.
    let content = fs::read_to_string(&ref_path)?;
    let chunks = content.as_bytes().chunks(chunk_size);

    let mut partial_indexes = Vec::new();
    for (i, chunk) in chunks.enumerate() {
        let partial = PartialFMIndex {
            reference_chunk: format!("chunk_{}", i),
            bwt_data: chunk.to_vec(),
        };
        partial_indexes.push(partial);
    }

    fs::write("partial_fm_indexes.json", serde_json::to_string_pretty(&partial_indexes)?)?;
    println!("Created partial_fm_indexes.json with {} chunks.", partial_indexes.len());
    Ok(())
}