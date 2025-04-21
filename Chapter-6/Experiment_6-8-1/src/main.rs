use anyhow::{Context, Result};
use clap::Parser;
use log::{info, error};
use rust_htslib::bam::{IndexedReader, Read};
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// A command-line interface for computing coverage in a specified genomic region of a BAM file.
#[derive(Parser, Debug)]
#[command(name = "coverage_tool", version, about = "A robust coverage calculation tool in Rust")]
struct Args {
    /// Path to the indexed BAM file. A matching BAI index is required.
    #[arg(long, value_name = "BAM_FILE")]
    bam: String,

    /// Genomic region to compute coverage over (e.g., "chr1:1-100000").
    #[arg(long, value_name = "REGION")]
    region: String,

    /// Output file to store the coverage result.
    #[arg(long, value_name = "OUTPUT_FILE")]
    out: String,
}

/// Reads the specified region from the BAM file, counting how many reads map there.
fn compute_coverage(bam_path: &str, region: &str) -> Result<u64> {
    let bam_indexed = IndexedReader::from_path(Path::new(bam_path))
        .with_context(|| format!("Failed to open or index BAM file: {bam_path}"))?;
    bam_indexed
        .fetch(region)
        .with_context(|| format!("Failed to fetch region {region}"))?;

    let mut coverage_count = 0u64;
    for record_result in bam_indexed.records() {
        let _record = record_result.with_context(|| {
            format!("Error reading record in region {region} from file {bam_path}")
        })?;
        coverage_count += 1;
    }
    Ok(coverage_count)
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Args::parse();
    info!("Computing coverage for region {} in BAM {}", args.region, args.bam);

    // Perform coverage calculation
    let coverage_result = match compute_coverage(&args.bam, &args.region) {
        Ok(cov) => cov,
        Err(e) => {
            error!("Coverage computation failed: {:?}", e);
            return Err(e);
        }
    };

    // Write the coverage result to the specified output file
    let mut out_file = File::create(&args.out)
        .with_context(|| format!("Could not create output file: {}", args.out))?;
    writeln!(out_file, "{}", coverage_result)?;
    info!("Coverage result: {}", coverage_result);

    Ok(())
}
