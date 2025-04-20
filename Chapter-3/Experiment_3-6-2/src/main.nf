[nextflow/main.nf]
#!/usr/bin/env nextflow

params.ref   = "${params.ref   ?: 'reference.fa'}"
params.reads = "${params.reads ?: 'reads/*.fq'}"

process BUILD_INDEX {
    container 'myregistry.io/bioinformatics/index-builder:latest'
    input:
    path reference_input from params.ref

    output:
    path "partial_fm_indexes.json"

    """
    index-builder ${reference_input} 1000000
    """
}

process ALIGN_READS {
    container 'myregistry.io/bioinformatics/aligner:latest'
    input:
    path read_file from channel.fromPath(params.reads)
    path partial_index from BUILD_INDEX
    output:
    path "alignment_results.json"

    """
    aligner ${read_file}
    """
}

process SUMMARIZE {
    container 'myregistry.io/bioinformatics/summarizer:latest'
    input:
    path results from ALIGN_READS.collect()
    output:
    path "final_summary.json"

    """
    summarizer ${results.join(' ')}
    """
}

workflow {
    main:
        SUMMARIZE(ALIGN_READS(BUILD_INDEX()))
}
