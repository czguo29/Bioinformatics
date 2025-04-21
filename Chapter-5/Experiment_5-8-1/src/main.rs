use anyhow::{Context, Result};
use rayon::prelude::*;
use rust_htslib::{bam, bam::Read};
use serde::{Serialize, Deserialize};
use std::collections::HashMap;
use std::fs::{File, create_dir_all};
use std::io::{BufReader, BufWriter};
use std::path::PathBuf;
use clap::Parser;

/// Represents a single correction patch for a contig at a given position.
#[derive(Serialize, Deserialize, Debug)]
struct CorrectionPatch {
    contig: String,
    position: u64,
    ref_base: u8,
    suggested_base: u8,
    coverage: u64,
}

/// Command-line arguments for chunked polishing.
#[derive(Parser, Debug)]
#[command(name = "polish_patches")]
#[command(about = "Collects partial polishing patches from short-read alignments")]
struct Args {
    /// The sorted BAM file containing short-read alignments.
    #[arg(long)]
    bam_input: PathBuf,

    /// Number of records to read from the BAM file in each chunk.
    #[arg(long, default_value_t = 10_000)]
    chunk_size: usize,

    /// Directory for partial correction patch outputs.
    #[arg(long, default_value = "partial_patches")]
    partial_outdir: PathBuf,

    /// Final JSON file with merged patch corrections.
    #[arg(long, default_value = "merged_corrections.json")]
    merged_output: PathBuf,
}

/// A container for serialized partial corrections to be merged later.
#[derive(Serialize, Deserialize, Debug)]
struct PartialPatches {
    patches: Vec<CorrectionPatch>,
}

/// Decodes 2-bit encoded nucleotides from rust-htslib.
fn decode_base(encoded: u8) -> u8 {
    match encoded {
        1 => b'A',
        2 => b'C',
        4 => b'G',
        8 => b'T',
        _ => b'N',
    }
}

/// Reads a chunk of BAM records up to the specified size.
fn read_chunk(bam_reader: &mut bam::Reader, chunk_size: usize) -> Result<Vec<bam::Record>> {
    let mut chunk = Vec::with_capacity(chunk_size);
    for _ in 0..chunk_size {
        if let Some(result) = bam_reader.records().next() {
            chunk.push(result?);
        } else {
            break;
        }
    }
    Ok(chunk)
}

/// Processes a list of BAM records in parallel, grouped by TID, generating correction patches.
fn collect_correction_patches(
    records: &[bam::Record],
    header: &bam::HeaderView,
) -> Vec<CorrectionPatch> {
    // Group records by TID.
    let mut contig_map: HashMap<i32, Vec<bam::Record>> = HashMap::new();
    for record in records {
        contig_map.entry(record.tid()).or_default().push(record.clone());
    }

    // For each TID group, detect mismatches in parallel.
    let patches: Vec<CorrectionPatch> = contig_map
        .par_iter()
        .flat_map(|(tid, recs)| {
            let contig_name = String::from_utf8_lossy(header.tid2name(*tid as u32)).to_string();
            let mut mismatch_map: HashMap<u64, (u8, u8, u64)> = HashMap::new();

            for r in recs {
                let start_pos = r.pos() as u64;
                let seq = r.seq();
                for (i, base_enc) in seq.iter().enumerate() {
                    let global_pos = start_pos + i as u64;
                    let read_base = decode_base(*base_enc);
                    if read_base != b'N' {
                        // For demonstration, assume reference base is 'A'.
                        if read_base != b'A' {
                            // Increment coverage for mismatched positions.
                            let entry = mismatch_map.entry(global_pos).or_insert((b'A', read_base, 0));
                            entry.2 += 1;
                        }
                    }
                }
            }

            let mut local_patches = Vec::new();
            for (pos, (ref_base, suggested, coverage)) in mismatch_map {
                local_patches.push(CorrectionPatch {
                    contig: contig_name.clone(),
                    position: pos,
                    ref_base,
                    suggested_base: suggested,
                    coverage,
                });
            }
            local_patches
        })
        .collect();

    patches
}

/// Merges two vectors of patches by simple concatenation. 
fn merge_patch_vectors(mut acc: Vec<CorrectionPatch>, mut new_patches: Vec<CorrectionPatch>) -> Vec<CorrectionPatch> {
    acc.append(&mut new_patches);
    acc
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Create the output directory for partial patches.
    create_dir_all(&args.partial_outdir)
        .with_context(|| format!("Failed to create partial patches directory at {:?}", args.partial_outdir))?;

    // Open the BAM file using rust-htslib.
    let mut bam_reader = bam::Reader::from_path(&args.bam_input)
        .with_context(|| format!("Failed to open BAM file {:?}", args.bam_input))?;
    let header = bam_reader.header().clone();

    // Read and process records in chunks.
    let mut chunk_index = 0;
    loop {
        let chunk = read_chunk(&mut bam_reader, args.chunk_size)?;
        if chunk.is_empty() {
            break;
        }

        let partial_patches = collect_correction_patches(&chunk, &header);
        let container = PartialPatches { patches: partial_patches };

        // Write the partial patches to disk for HPC ephemeral tasks.
        let outfile_path = args
            .partial_outdir
            .join(format!("partial_patches_chunk_{}.json", chunk_index));
        let out_file = File::create(&outfile_path)
            .with_context(|| format!("Failed to create partial patch file {:?}", outfile_path))?;
        serde_json::to_writer(BufWriter::new(out_file), &container)
            .with_context(|| format!("Failed to serialize partial patches to {:?}", outfile_path))?;

        println!(
            "Processed chunk {} with {} records. Partial patches stored at {:?}.",
            chunk_index,
            chunk.len(),
            outfile_path
        );
        chunk_index += 1;
    }

    // Merge all partial patch files into a single final output.
    let dir_entries = std::fs::read_dir(&args.partial_outdir)
        .with_context(|| format!("Failed to read partial patches directory {:?}", args.partial_outdir))?;

    let mut merged_patches = Vec::new();
    for entry in dir_entries {
        let path = entry?.path();
        if path.file_name().map_or(false, |p| p.to_string_lossy().starts_with("partial_patches_chunk_")) {
            let file = File::open(&path)
                .with_context(|| format!("Failed to open partial patches at {:?}", path))?;
            let partial: PartialPatches =
                serde_json::from_reader(BufReader::new(file))
                    .with_context(|| format!("Failed to parse partial patches from {:?}", path))?;
            merged_patches = merge_patch_vectors(merged_patches, partial.patches);
        }
    }

    // Write final merged corrections to JSON.
    let merged_file = File::create(&args.merged_output)
        .with_context(|| format!("Failed to create merged patch file {:?}", args.merged_output))?;
    serde_json::to_writer(BufWriter::new(merged_file), &merged_patches)
        .with_context(|| format!("Failed to write merged patches to {:?}", args.merged_output))?;

    println!("Merged {} total correction patches into {:?}.", merged_patches.len(), args.merged_output);
    Ok(())
}
