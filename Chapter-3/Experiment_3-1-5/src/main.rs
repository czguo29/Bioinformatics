use mpi::traits::*;
use std::collections::HashMap;
use serde::{Serialize, Deserialize};
use bincode;
use std::fs::File;
use std::io::Write;
use std::env;
/// A chunk of a genomic index built by one MPI rank (e.g., a node in an HPC cluster).
/// It’s like a piece of a DNA search engine for part of the genome.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct PartialIndex {
    rank_id: i32,                     // Which rank built this chunk
    index_data: HashMap<String, usize>, // Mini-database of index data
}
/// Builds a partial index for this rank’s slice of the genome.
/// In practice, this might index DNA sequences (e.g., for a suffix array).
fn build_local_index(rank: i32) -> PartialIndex {
    let mut data = HashMap::new();
    data.insert(format!("key_rank_{}", rank), rank as usize); // Simulate index entry
    PartialIndex {
        rank_id: rank,
        index_data: data,
    }
}
/// Merges partial indexes into a global index.
/// Like combining puzzle pieces into a full genome index.
fn merge_indexes(chunks: Vec<PartialIndex>) -> HashMap<String, usize> {
    let mut global_map = HashMap::new();
    for chunk in chunks {
        for (k, v) in chunk.index_data {
            global_map.insert(k, v); // Copy entries to global index
        }
    }
    global_map
}
fn main() {
    // Log MPI environment to confirm communicator setup
    if let Ok(rank) = env::var("OMPI_COMM_WORLD_RANK") {
        println!("OMPI_COMM_WORLD_RANK: {}", rank);
    }
    if let Ok(size) = env::var("OMPI_COMM_WORLD_SIZE") {
        println!("OMPI_COMM_WORLD_SIZE: {}", size);
    }
    // Initialize MPI for parallel processing
    let universe = mpi::initialize().expect("Failed to initialize MPI");
    let world = universe.world();
    let rank = world.rank(); // Rank ID (0, 1, 2, 3 for -np 4)
    let size = world.size(); // Number of ranks (4 for -np 4)
    println!("Rank {} of {} starting", rank, size);
    // Build and serialize this rank’s index
    let local_chunk = build_local_index(rank);
    let serialized_chunk = bincode::serialize(&local_chunk)
        .expect("Failed to serialize partial index");
    // Sync ranks before communication
    world.barrier();
    // Rank 0 collects all indexes; others send
    let mut gathered_chunks = Vec::new();
    if rank == 0 {
        gathered_chunks.push(serialized_chunk); // Rank 0’s own chunk
        for source_rank in 1..size {
            // Receive serialized chunk from each rank
            let (received_chunk, _status) = world
                .process_at_rank(source_rank)
                .receive_vec::<u8>();
            gathered_chunks.push(received_chunk);
        }
    } else {
        // Non-root ranks send to rank 0
        world.process_at_rank(0).send(&serialized_chunk[..]);
    }
    // Sync after communication
    world.barrier();
    // Only rank 0 merges and saves
    if rank == 0 {
        // Deserialize chunks
        let deserialized_chunks: Vec<PartialIndex> = gathered_chunks
            .into_iter()
            .map(|bytes| bincode::deserialize(&bytes).expect("Failed to deserialize"))
            .collect();
        // Merge into global index
        let global_index = merge_indexes(deserialized_chunks);
        // Save to output.txt
        let mut file = File::create("output.txt")
            .expect("Failed to create output.txt");
        writeln!(file, "Final merged index: {:?}", global_index)
            .expect("Failed to write to output.txt");
        // Print for verification
        println!("Final merged index on rank 0: {:?}", global_index);
    }
}
