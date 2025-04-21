use serde::{Serialize, Deserialize}; // For serializing/deserializing JSON
use rayon::prelude::*; // For parallel iterators
use std::fs::File;
use std::io::{BufReader, BufWriter};

#[derive(Serialize, Deserialize)]
struct EqtlAssoc {
    snp_id: String,
    gene_id: String,
    p_value: f64,
}

#[derive(Serialize, Deserialize)]
struct PeakCall {
    chrom: String,
    start: usize,
    end: usize,
    peak_id: String,
}

#[derive(Serialize, Deserialize)]
struct MotifHit {
    chrom: String,
    position: usize,
    motif_id: String,
}

// Changed to `pub struct` to match visibility requirements of `merge_multiomics_data`.
#[derive(Serialize, Deserialize)]
pub struct IntegratedResult {
    pub snp_id: String,
    pub gene_id: String,
    pub p_value: f64,
    pub peak_id: Option<String>,
    pub motif_id: Option<String>,
}

pub fn merge_multiomics_data(
    eqtl_files: Vec<String>,
    peak_files: Vec<String>,
    motif_files: Vec<String>
) -> Vec<IntegratedResult> {
    // Load partial eQTL associations in parallel
    let eqtls = eqtl_files.par_iter()
        .flat_map(|path| load_eqtl_assoc(path))
        .collect::<Vec<_>>();

    // Load partial peak calls in parallel (currently unused, so we prefix with _)
    let _peaks = peak_files.par_iter()
        .flat_map(|path| load_peak_calls(path))
        .collect::<Vec<_>>();

    // Load partial motif hits in parallel (currently unused, so we prefix with _)
    let _motifs = motif_files.par_iter()
        .flat_map(|path| load_motif_hits(path))
        .collect::<Vec<_>>();

    // Dummy integration logic
    // In real code, we'd unify by genomic coordinates or known variant IDs
    eqtls.par_iter().map(|eqtl| {
        let has_peak = eqtl.snp_id.contains("peak");
        let has_motif = eqtl.snp_id.contains("motif");
        IntegratedResult {
            snp_id: eqtl.snp_id.clone(),
            gene_id: eqtl.gene_id.clone(),
            p_value: eqtl.p_value,
            peak_id: if has_peak { Some("PEAK123".to_string()) } else { None },
            motif_id: if has_motif { Some("MOTIF456".to_string()) } else { None },
        }
    }).collect()
}

fn load_eqtl_assoc(path: &str) -> Vec<EqtlAssoc> {
    let f = File::open(path).expect("Failed to open eQTL file");
    serde_json::from_reader(BufReader::new(f)).expect("Failed to parse eQTL JSON")
}

fn load_peak_calls(path: &str) -> Vec<PeakCall> {
    let f = File::open(path).expect("Failed to open peak file");
    serde_json::from_reader(BufReader::new(f)).expect("Failed to parse Peak JSON")
}

fn load_motif_hits(path: &str) -> Vec<MotifHit> {
    let f = File::open(path).expect("Failed to open motif file");
    serde_json::from_reader(BufReader::new(f)).expect("Failed to parse Motif JSON")
}

fn main() {
    let eqtl_files = vec![r"C:\Users\Christyane Zabdi\Downloads\CZ-UI S2\Sem 2\BGVR\Chapter 4\Experiment_4-8-1\src\eqtl_part1.json".to_string(), r"C:\Users\Christyane Zabdi\Downloads\CZ-UI S2\Sem 2\BGVR\Chapter 4\Experiment_4-8-1\src\eqtl_part2.json".to_string()];
    let peak_files = vec![r"C:\Users\Christyane Zabdi\Downloads\CZ-UI S2\Sem 2\BGVR\Chapter 4\Experiment_4-8-1\src\peak_part1.json".to_string()];
    let motif_files = vec![r"C:\Users\Christyane Zabdi\Downloads\CZ-UI S2\Sem 2\BGVR\Chapter 4\Experiment_4-8-1\src\motif_part1.json".to_string()];

    let integrated_data = merge_multiomics_data(eqtl_files, peak_files, motif_files);
    let out_file = File::create("integrated.json").expect("Cannot create output file");
    serde_json::to_writer_pretty(BufWriter::new(out_file), &integrated_data)
        .expect("Failed to write integrated results");
}