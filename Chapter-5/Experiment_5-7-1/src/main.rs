use anyhow::{Context, Result};
use std::sync::{Arc, Mutex};
use rust_htslib::bam::{self, Read, Reader, Record};
use rayon::prelude::*;
use serde::{Serialize, Deserialize};
use std::collections::HashMap;
use std::{
    fs::{File, create_dir_all},
    io::{BufReader, BufWriter},
    path::PathBuf,
};
use clap::Parser;

#[derive(Serialize, Deserialize, Debug)]
struct TranscriptCount {
    transcript_id: String,
    count: u64,
}

#[derive(Debug)]
struct GtfLookup {
    intervals: Vec<(String, u64, u64, String)>,
}

impl GtfLookup {
    fn new(gtf_path: &str) -> Self {
        println!("Loading annotation from {}", gtf_path);
        // TODO: Implementasi parsing GTF seharusnya ada di sini
        GtfLookup { intervals: vec![] }
    }

    fn transcripts_for_region(&self, chrom: &str, start: u64, end: u64) -> Vec<String> {
        let transcripts: Vec<String> = self.intervals.iter()
            .filter_map(|(c, s, e, tid)| {
                if c == chrom && !(end < *s || start > *e) {
                    Some(tid.clone())
                } else {
                    None
                }
            })
            .collect();
        
        if transcripts.is_empty() {
            println!("No transcripts found for region: {}:{}-{}", chrom, start, end);
        }
        transcripts
    }
}

#[derive(Parser, Debug)]
struct Args {
    #[arg(long)]
    bam_input: PathBuf,
    #[arg(long)]
    annotation: PathBuf,
    #[arg(long, default_value_t = 10_000)]
    chunk_size: usize,
    #[arg(long, default_value = "partial_counts")]
    partial_outdir: PathBuf,
    #[arg(long, default_value = "merged_counts.json")]
    merged_output: PathBuf,
}

#[derive(Serialize, Deserialize, Debug)]
struct PartialTranscriptCounts {
    counts: HashMap<String, u64>,
}

fn parallel_count_transcripts(records: &[bam::Record], tid2name_map: &HashMap<u32, String>, gtf: &GtfLookup) -> HashMap<String, u64> {
    records
        .par_iter()
        .fold(HashMap::new, |mut local_map, record| {
            if let Some(chrom) = tid2name_map.get(&(record.tid() as u32)) {
                let start = record.pos() as u64;
                let end = start + record.seq().len() as u64;
                let overlapping_tids = gtf.transcripts_for_region(chrom, start, end);
                for tid in overlapping_tids {
                    *local_map.entry(tid).or_insert(0) += 1;
                }
            }
            local_map
        })
        .reduce(HashMap::new, |mut global_map, local_map| {
            for (k, v) in local_map {
                *global_map.entry(k).or_insert(0) += v;
            }
            global_map
        })
}

fn main() -> Result<()> {
    let args = Args::parse();

    let bam = Reader::from_path(&args.bam_input).context("Failed to open BAM file")?;
    let header = bam.header().clone();
    let tid2name_map: HashMap<u32, String> = (0..header.target_count())
        .map(|tid| (tid, String::from_utf8_lossy(header.tid2name(tid)).to_string()))
        .collect();

    let gtf = GtfLookup::new(args.annotation.to_str().unwrap());
    let mut bam_reader = Reader::from_path(&args.bam_input).context("Failed to open BAM file")?;
    let mut records = vec![];
    for record in bam_reader.records() {
        records.push(record.context("Failed to read BAM record")?);
    }

    let result: Vec<_> = records
        .par_iter()
        .map(|record| {
            let chrom = tid2name_map.get(&(record.tid() as u32)).cloned().unwrap_or_default();
            let start = record.pos() as u64;
            let end = start + record.seq().len() as u64;
            (chrom, start, end)
        })
        .collect();

    println!("{:?}", result);

    let dir_entries = std::fs::read_dir(&args.partial_outdir)
        .with_context(|| format!("Failed to read partial output directory {:?}", args.partial_outdir))?;

    let mut merged_counts = HashMap::new();
    for entry in dir_entries {
        let path = entry?.path();
        if path.file_name().map_or(false, |p| p.to_string_lossy().starts_with("partial_counts_chunk_")) {
            let file = File::open(&path).context("Failed to open partial counts file")?;
            let partial: PartialTranscriptCounts = serde_json::from_reader(BufReader::new(file)).context("Failed to parse JSON")?;
            for (tid, count) in partial.counts {
                *merged_counts.entry(tid).or_insert(0) += count;
            }
        }
    }

    let final_vec: Vec<TranscriptCount> = merged_counts.into_iter()
        .map(|(tid, count)| TranscriptCount { transcript_id: tid, count })
        .collect();

    let merged_file = File::create(&args.merged_output).context("Failed to create merged output file")?;
    serde_json::to_writer(BufWriter::new(merged_file), &final_vec).context("Failed to write merged transcript counts")?;

    println!("Merged counts for {} transcripts into {:?}", final_vec.len(), args.merged_output);

    Ok(())
}