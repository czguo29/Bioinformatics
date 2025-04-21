use anyhow::{Context, Result};
use clap::Parser;
use env_logger;
use log::{error, info};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
#[derive(Parser, Debug)]
#[command(name = "coverage_tool")]
#[command(author, version, about = "A robust coverage analysis tool using rust-htslib")]
struct Cli {
    /// Path to the indexed BAM file
    #[arg(long, value_name = "FILE")]
    bam: PathBuf,
    /// One or more regions specified as chr:start-end
    #[arg(long, value_name = "REGION", required = true, num_args = 1..)]
    region: Vec<String>,
    /// Optional output file to store coverage data
    #[arg(long, value_name = "OUT_FILE")]
    output: Option<PathBuf>,
}
fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    info!("Starting coverage analysis with input: {:?}", cli.bam);
    info!("Processing the following regions: {:?}", cli.region);
    // Use a thread-safe reference-counted pointer to a mutex-protected vector.
    // Each thread will push its local coverage counts into this shared vector.
    let coverage_map = Arc::new(Mutex::new(Vec::new()));
    // Parallel iteration over regions
    cli.region.par_iter().for_each(|region| {
        // Open the BAM file using IndexedReader for random-access fetching
        let mut bam = match bam::IndexedReader::from_path(&cli.bam) {
            Ok(reader) => reader,
            Err(e) => {
                error!("Failed to open BAM file {}: {}", cli.bam.display(), e);
                return; // Skip this region if the file cannot be opened
            }
        };
        // Fetch only the relevant region from the BAM file
        if let Err(e) = bam.fetch(region) {
            error!("Failed to fetch region {}: {}", region, e);
            return; // Skip this region if fetch fails
        }
        // Collect local coverage counts for each record in this region
        let mut local_counts = Vec::new();
        for record_result in bam.records() {
            match record_result {
                Ok(record) => {
                    // Here we measure coverage by counting the length of the sequence.
                    // More sophisticated analyses can be done, such as counting aligned bases.
                    local_counts.push(record.seq().len());
                }
                Err(e) => {
                    error!("Error reading record in region {}: {}", region, e);
                }
            }
        }
        // Lock the shared coverage map for thread-safe access and extend it with local counts
        let mut shared_map = coverage_map.lock().expect("Mutex is poisoned");
        shared_map.extend(local_counts);
    });
    // Final coverage results are stored in coverage_map
    let final_map = coverage_map.lock().expect("Mutex is poisoned");
    let total_reads = final_map.len();
    info!(
        "Completed coverage analysis for {} region(s). Total reads processed: {}",
        cli.region.len(),
        total_reads
    );
    // Optionally write coverage data to an output file if specified
    if let Some(out_path) = &cli.output {
        let mut file = File::create(out_path)
            .with_context(|| format!("Could not create output file {}", out_path.display()))?;
        writeln!(file, "Coverage data for {} reads", total_reads)?;
        info!("Coverage data has been written to {}", out_path.display());
    }
    Ok(())
}
