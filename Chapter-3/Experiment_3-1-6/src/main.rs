use needletail::parse_fastx_file;
use serde::{Serialize, Deserialize};
use std::path::PathBuf;
use std::fs;
use std::io::BufWriter;
use std::error::Error;
#[derive(Debug, Serialize, Deserialize)]
struct SeqRecord {
    id: String,
    seq_len: usize,
    gc_content: f64,
}
fn main() -> Result<(), Box<dyn Error>> {
    // Take a FASTA file path from command line or use "example.fasta" by default.
    let fasta_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| r"C:\Users\Christyane Zabdi\Downloads\CZ-UI S2\Sem 2\BGVR\Chapter 3\Experiment_3-1-6\src\example.fasta".to_string());
    // Read all records sequentially from the FASTA/FASTQ file into a vector.
    let mut reader = parse_fastx_file(&PathBuf::from(&fasta_path))?;
    let mut records = Vec::new();
    while let Some(record) = reader.next() {
        let rec = record?;
        let seq = rec.seq(); // This returns a `Cow<'_, [u8]>`.
        // Convert the Cow<'_, [u8]> into a &[u8].
        let gc_val = calc_gc_content(seq.as_ref());
        let seq_id = String::from_utf8_lossy(rec.id()).to_string();
        records.push(SeqRecord {
            id: seq_id,
            seq_len: seq.len(),
            gc_content: gc_val,
        });
    }
    // Perform any transformations on each record here (e.g., in parallel if desired).
    let file = fs::File::create("output.json")?;
    let writer = BufWriter::new(file);
    serde_json::to_writer(writer, &records)?;
    println!("Processed {} sequences", records.len());
    Ok(())
}
fn calc_gc_content(seq: &[u8]) -> f64 {
    let gc_count = seq.iter().filter(|&&c| c == b'G' || c == b'C').count();
    gc_count as f64 / (seq.len().max(1) as f64)
}
