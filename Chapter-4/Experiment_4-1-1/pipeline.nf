#!/usr/bin/env nextflow

/*
 * nextflow.config can be placed in the same directory as this script
 * to configure containers, executors, etc. For a stand-alone pipeline,
 * you only need this file and the Rust program code (Cargo.toml, main.rs).
 */

params.synthetic_fastq = (params.synthetic_fastq ?: 'synthetic_reads.fastq')

process buildPWMandMRF {
    // If Rust is installed locally, you can compile on the fly.
    // If containerized, your Dockerfile or environment must have Rust installed.
    input:
    file 'synthetic_reads.fastq' from params.synthetic_fastq

    output:
    file 'pwm_results.txt'
    file 'mrf_results.txt'

    """
    # Compile the Rust code (assumes Cargo.toml & src/main.rs are in 'rust_code' directory)
    cd rust_code
    cargo build --release

    # Run the Rust binary with the synthetic FASTQ file
    ./target/release/bioinformatics_tools ../synthetic_reads.fastq pwm_results.txt mrf_results.txt
    """
}

workflow {
    main:
        buildPWMandMRF()
}