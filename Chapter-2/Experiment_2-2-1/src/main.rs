use std::collections::HashMap;
use bio::io::fasta;
use ndarray::Array2;
use nalgebra::DMatrix;

/// Stores a potential or weight for an edge in the Markov Random Field (MRF).
/// In a real-world scenario, this could be a more complex data structure 
/// or a function that depends on the nucleotides involved.
#[derive(Debug)]
struct MRFEdge {
    potential: f32,
}

/// Represents the Markov Random Field as a graph-like structure.
/// Each node is identified by a (sequence_id, position_in_sequence) tuple.
/// `adjacency` maps each node to a list of neighboring nodes (with associated edges).
#[derive(Debug)]
struct MRF {
    adjacency: HashMap<(usize, usize), Vec<((usize, usize), MRFEdge)>>,
}

/// Constructs a simplified pairwise MRF from a set of sequences.
/// For each position in each sequence, an edge is created to the next position.
/// This can be extended to model more complex interactions (e.g., long-range edges).
fn build_mrf_from_sequences(sequences: &[String]) -> MRF {
    // A HashMap where the key is a node (seq_id, position) 
    // and the value is a vector of (neighbor_node, edge_info).
    let mut adjacency: HashMap<(usize, usize), Vec<((usize, usize), MRFEdge)>> = HashMap::new();

    // Enumerate each sequence so we know which sequence we're in (seq_id),
    // and then iterate over each character in the sequence (except the last one).
    for (seq_id, seq) in sequences.iter().enumerate() {
        for i in 0..seq.len().saturating_sub(1) {
            // The current position (node_a) and the next position (node_b).
            let node_a = (seq_id, i);
            let node_b = (seq_id, i + 1);

            // We'll assign a simple constant potential to each edge for demonstration.
            // This might be calculated based on the two nucleotides in a more advanced application.
            let edge = MRFEdge { potential: 1.0 };

            // Insert the edge into the adjacency structure.
            // node_a connects to node_b with the assigned potential.
            adjacency
                .entry(node_a)
                .or_insert_with(Vec::new)
                .push((node_b, edge));
        }
    }

    // Return the constructed MRF that captures pairwise adjacency between sequence positions.
    MRF { adjacency }
}

fn main() {
    // Create a Reader for the FASTA file named "reads.fasta" located in the "src" folder.
    // If the file is missing or inaccessible, the program will exit with an error message.
    let reader = fasta::Reader::from_file("src/reads.fasta")
        .expect("Error opening FASTA file in the 'src' directory");

    // Read all sequences from the FASTA into a vector of strings.
    // The 'bio' crate handles FASTA parsing, and we convert each record's byte slice to UTF-8 text.
    let seqs: Vec<String> = reader
        .records()
        .map(|rec| {
            let r = rec.expect("Invalid record");
            String::from_utf8(r.seq().to_vec()).expect("Non-UTF8 sequence")
        })
        .collect();

    // Build the MRF from the parsed sequences.
    let mrf = build_mrf_from_sequences(&seqs);

    // Demonstrate using nalgebra and ndarray crates for potential HPC tasks or matrix operations.
    let matrix_nalgebra = DMatrix::<f32>::from_element(10, 10, 1.0);
    let matrix_ndarray = Array2::<f32>::ones((10, 10));

    // Print out information about the constructed MRF and the example matrices.
    println!(
        "Constructed MRF with {} nodes.",
        mrf.adjacency.len()
    );
    println!(
        "nalgebra matrix dimensions: {} x {}",
        matrix_nalgebra.nrows(),
        matrix_nalgebra.ncols()
    );
    println!(
        "ndarray matrix dimensions: {} x {}",
        matrix_ndarray.nrows(),
        matrix_ndarray.ncols()
    );

    // Show details for up to 5 nodes. Each node has edges that link it to subsequent positions.
    // The potential field in MRFEdge is printed to illustrate that it is being used.
    for (node, edges) in mrf.adjacency.iter().take(5) {
        println!("Node {:?} has edges:", node);
        for (neighbor, edge) in edges {
            println!("  -> {:?} with potential = {}", neighbor, edge.potential);
        }
    }
}