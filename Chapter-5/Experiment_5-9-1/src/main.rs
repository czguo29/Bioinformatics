use rayon::prelude::*;
use serde::{Serialize, Deserialize};
use clap::{Arg, Command};
use anyhow::{Context, Result};
use std::fs;
use std::path::PathBuf;
use std::io::{BufReader, BufWriter};
use std::collections::VecDeque;

/// A structure to represent a simple variant in JSON form.
#[derive(Debug, Serialize, Deserialize)]
struct Variant {
    chrom: String,
    pos: u64,
    ref_base: String,
    alt_base: String,
    score: f64,
}

/// A container for partial variant outputs.
#[derive(Debug, Serialize, Deserialize)]
struct PartialVariants {
    variants: Vec<Variant>,
}

/// Demonstrates parallel alignment by distributing read files among available threads.
/// In HPC scenarios, ephemeral tasks might handle these reads in parallel or chunks.
fn parallel_align(read_files: &[String], reference: &str) -> Vec<String> {
    read_files
        .par_iter()
        .map(|reads| {
            // Placeholder logic; actual pipeline might use rust-htslib for alignment.
            format!("{} aligned to {}", reads, reference)
        })
        .collect()
}

/// Demonstrates a parallel variant detection step that merges partial alignment outputs.
/// In real HPC usage, ephemeral tasks could each return partial variant calls for merging.
fn detect_variants(aligned: &[String]) -> Vec<Variant> {
    aligned
        .par_iter()
        .map(|chunk| {
            // Placeholder logic returning a single test variant.
            Variant {
                chrom: "chr1".to_string(),
                pos: 12345,
                ref_base: "A".to_string(),
                alt_base: "G".to_string(),
                score: 42.0,
            }
        })
        .collect()
}

/// Merges two sets of variants by simple concatenation.
fn merge_variants(mut existing: Vec<Variant>, mut new_variants: Vec<Variant>) -> Vec<Variant> {
    existing.append(&mut new_variants);
    existing
}

fn main() -> Result<()> {
    let matches = Command::new("multi_step_pipeline")
        .version("0.1.0")
        .about("Align reads and call variants in HPC concurrency-friendly steps.")
        .arg(Arg::new("reads")
            .short('r')
            .long("reads")
            .required(true)
            .help("Paths to read files for alignment")
            .num_args(1..))
        .arg(Arg::new("ref")
            .short('f')
            .long("ref")
            .required(true)
            .help("Path to the reference file"))
        .arg(Arg::new("out")
            .short('o')
            .long("out")
            .default_value("final_variants.json")
            .help("Name of the final merged variants JSON file"))
        .arg(Arg::new("chunk_size")
            .long("chunk-size")
            .default_value("5")
            .help("Number of read files to process per chunk in HPC usage"))
        .arg(Arg::new("partial_dir")
            .long("partial-dir")
            .default_value("partial_variants")
            .help("Directory in which partial variant results will be stored"))
        .get_matches();

    // Gather command-line arguments.
    let read_paths: Vec<String> = matches.get_many::<String>("reads")
        .unwrap()
        .map(String::from)
        .collect();
    let ref_path = matches.get_one::<String>("ref").unwrap();
    let out_path = matches.get_one::<String>("out").unwrap();
    let chunk_size: usize = matches.get_one::<String>("chunk_size")
        .unwrap()
        .parse()
        .context("Invalid chunk_size")?;
    let partial_dir = matches.get_one::<String>("partial_dir").unwrap();

    // Create the partial directory for ephemeral HPC usage, if it does not exist.
    fs::create_dir_all(partial_dir)
        .with_context(|| format!("Failed to create or open partial directory: {}", partial_dir))?;

    // Break the read file paths into chunks for partial processing.
    let mut queue: VecDeque<String> = read_paths.into_iter().collect();
    let mut chunk_index = 0;

    while !queue.is_empty() {
        let this_chunk: Vec<String> = queue.drain(..chunk_size.min(queue.len())).collect();
        if this_chunk.is_empty() {
            break;
        }

        // In a real HPC pipeline, ephemeral tasks might each receive a chunk, do alignment, and store partial output.
        let alignment_results = parallel_align(&this_chunk, ref_path);
        let variants = detect_variants(&alignment_results);
        let container = PartialVariants { variants };

        let chunk_file = format!("{}/partial_variants_{:04}.json", partial_dir, chunk_index);
        let out_file = fs::File::create(&chunk_file)
            .with_context(|| format!("Failed to create partial variant file {}", chunk_file))?;
        serde_json::to_writer_pretty(BufWriter::new(out_file), &container)
            .with_context(|| format!("Failed to serialize partial variants to {}", chunk_file))?;

        println!(
            "Processed chunk {} with {} read files. Stored partial variants in {}.",
            chunk_index,
            this_chunk.len(),
            chunk_file
        );
        chunk_index += 1;
    }

    // Merge all partial variant outputs into one final JSON.
    let dir_entries = fs::read_dir(partial_dir)
        .with_context(|| format!("Failed to read directory for partial variants: {}", partial_dir))?;
    let mut merged_results = Vec::new();

    for entry in dir_entries {
        let path = entry?.path();
        if path.file_name().map_or(false, |f| f.to_string_lossy().starts_with("partial_variants_")) {
            let file = fs::File::open(&path)
                .with_context(|| format!("Failed to open partial variants file: {:?}", path))?;
            let partial: PartialVariants = serde_json::from_reader(BufReader::new(file))
                .with_context(|| format!("Failed to parse partial variants from: {:?}", path))?;
            merged_results = merge_variants(merged_results, partial.variants);
        }
    }

    // Write out the final merged variants.
    fs::write(out_path, serde_json::to_string_pretty(&merged_results)?)
        .with_context(|| format!("Failed to write final merged variants to {}", out_path))?;

    println!(
        "Merged {} variants into the final output: {}",
        merged_results.len(),
        out_path
    );
    Ok(())
}