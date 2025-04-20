use bio::io::fasta;
use rayon::prelude::*;
use rust_htslib::{
    bam,
    bam::Read as BamRead,
    bcf::{self, Read as BcfRead},
};
use std::error::Error;
use std::path::Path;

/// Count occurrences of a given `motif` in `seq`, allowing overlaps.
fn count_occurrences(seq: &str, motif: &str) -> usize {
    let mut count = 0;
    let mut start = 0;
    while let Some(pos) = seq[start..].find(motif) {
        count += 1;
        start += pos + 1;
    }
    count
}

/// Read FASTQ (or FASTA) via bio and tally occurrences of `motif` in parallel.
fn process_fastq(fastq_path: &Path, motif: &str) -> Result<usize, Box<dyn Error>> {
    let reader = fasta::Reader::from_file(fastq_path)?;

    // Collect sequences as Strings
    let sequences: Vec<String> = reader
        .records()
        .filter_map(|r| match r {
            Ok(rec) => Some(String::from_utf8_lossy(rec.seq()).to_string()),
            Err(_) => None,
        })
        .collect();

    // Parallel counting
    let total_occurrences: usize = sequences
        .par_iter()
        .map(|seq| count_occurrences(seq, motif))
        .sum();

    Ok(total_occurrences)
}

/// Read SAM/BAM records using rust-htslib's iterator approach, then parallel motif counting.
fn process_bam(bam_path: &Path, motif: &str) -> Result<usize, Box<dyn Error>> {
    // Declare `bam_reader` as mutable to satisfy the `records()` method requirement.
    let mut bam_reader = bam::Reader::from_path(bam_path)?;
    let mut sequences = Vec::new();

    for record_result in bam_reader.records() {
        let record = record_result?;
        let seq_string = String::from_utf8_lossy(&record.seq().as_bytes()).to_string();
        sequences.push(seq_string);
    }

    // Parallel counting
    let total_occurrences: usize = sequences
        .par_iter()
        .map(|seq| count_occurrences(seq, motif))
        .sum();

    Ok(total_occurrences)
}

/// Read VCF/BCF via rust-htslib, using `record.alleles()` to get REF+ALT. Then parallel motif counting.
fn process_vcf(vcf_path: &Path, motif: &str) -> Result<usize, Box<dyn Error>> {
    let mut vcf_reader = bcf::Reader::from_path(vcf_path)?;
    let mut variants = Vec::new();

    // Each VCF record can have multiple alleles. The first is REF, the rest are ALT.
    for rec_result in vcf_reader.records() {
        let record = rec_result?;
        let alleles = record.alleles();

        // REF allele is alleles[0]
        let ref_allele = String::from_utf8_lossy(alleles[0]).to_string();

        // ALT alleles are the subsequent ones
        let alt_alleles: Vec<String> = alleles[1..]
            .iter()
            .map(|alt| String::from_utf8_lossy(alt).to_string())
            .collect();
        let alt_string = alt_alleles.join(",");

        variants.push(format!("{}{}", ref_allele, alt_string));
    }

    // Parallel counting
    let total_occurrences: usize = variants
        .par_iter()
        .map(|v| count_occurrences(v, motif))
        .sum();

    Ok(total_occurrences)
}

/// Main orchestrator
fn main() -> Result<(), Box<dyn Error>> {
    let motif = "GATTACA";

    let fastq_path = Path::new("example.fastq");
    let bam_path = Path::new("example.bam");
    let vcf_path = Path::new("example.vcf");

    // FASTQ
    let fastq_matches = process_fastq(fastq_path, motif)?;
    println!("FASTQ: Found motif '{}' {} times.", motif, fastq_matches);

    // BAM
    let bam_matches = process_bam(bam_path, motif)?;
    println!("BAM: Found motif '{}' {} times.", motif, bam_matches);

    // VCF
    let vcf_matches = process_vcf(vcf_path, motif)?;
    println!("VCF: Found motif '{}' {} times.", motif, vcf_matches);

    Ok(())
}