use std::collections::HashMap;
use bio::io::fasta;
use ndarray::Array2;
use nalgebra::DMatrix;
use rayon::prelude::*;

///
/// Builds a De Bruijn graph using a k-mer size of `k`. Each thread processes 
/// a subset of sequences in parallel, accumulating a local HashMap of node->neighbors.
/// The local HashMaps are then merged (reduced) into a global HashMap. 
///
/// This design scales well to millions of reads because it avoids locking a single 
/// shared data structure. Instead, each thread builds partial results and 
/// reduces them at the end, minimizing contention and ensuring memory safety.
///
fn build_de_bruijn(k: usize, sequences: &[String]) -> HashMap<String, Vec<String>> {
    // Parallel map each sequence to a local HashMap<String, Vec<String>>
    // and then reduce (merge) all of these local maps into a single global map.
    sequences
        .par_iter()
        .map(|seq| {
            let mut local_map = HashMap::new();

            for i in 0..seq.len().saturating_sub(k) {
                let node = &seq[i..i + k];
                let edge = &seq[i + 1..i + k + 1];
                local_map
                    .entry(node.to_string())
                    .or_insert_with(Vec::new)
                    .push(edge.to_string());
            }

            local_map
        })
        // The initial empty HashMap (identity) for the reduce step.
        .reduce(
            || HashMap::new(),
            // Merging function: combine two HashMaps by appending edge lists.
            |mut acc, local_map| {
                for (key, mut edges) in local_map {
                    acc.entry(key).or_insert_with(Vec::new).append(&mut edges);
                }
                acc
            },
        )
}

fn main() {
    // Read sequences from a FASTA file in the 'src' directory using rust-bio's fasta::Reader.
    let reader = fasta::Reader::from_file("src/reads.fasta")
        .expect("Cannot open FASTA file in 'src' directory");
    let sequences: Vec<String> = reader
        .records()
        .map(|r| {
            let record = r.expect("Invalid FASTA record");
            String::from_utf8(record.seq().to_vec()).expect("Sequence is not valid UTF-8")
        })
        .collect();

    // Choose a reasonable k-mer size (21 is common for short-read assemblies). 
    // In production, pick k based on read length, coverage, and error profile.
    let k = 21;
    let graph = build_de_bruijn(k, &sequences);

    // Demonstrate HPC-oriented numeric crates: create small example matrices using nalgebra and ndarray.
    let test_matrix = DMatrix::<f32>::from_element(5, 5, 1.0);
    let nd_array = Array2::<f32>::ones((5, 5));

    println!("Constructed De Bruijn graph with {} nodes.", graph.len());
    println!("nalgebra matrix: {} x {}", test_matrix.nrows(), test_matrix.ncols());
    println!("ndarray shape: {} x {}", nd_array.nrows(), nd_array.ncols());
}
