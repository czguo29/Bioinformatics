params.bam_list   = 'bams.txt'
params.region_list = 'regions.txt'

process coverageCalc {
    // Container image that includes the compiled Rust coverage tool and its dependencies
    container 'myrust/coverage:1.0'
    input:
    val bam_file from bamChannel
    val region from regionChannel

    output:
    file("coverage_${bam_file}_${region}.txt")

    """
    # The Rust coverage binary is invoked here, using command-line arguments
    # to specify the BAM file, region, and output file.
    rust_coverage_tool \
      --bam ${bam_file} \
      --region ${region} \
      --out coverage_${bam_file}_${region}.txt
    """
}

process mergeCoverage {
    input:
    file coverage_files from coverageCalc.out.collect()
    output:
    file "merged_coverage.txt"

    """
    # The partial coverage results are concatenated into a single file,
    # simplifying subsequent analysis or archival.
    cat ${coverage_files.join(" ")} > merged_coverage.txt
    """
}

workflow {
    // Read BAM file paths from bams.txt
    bamChannel = Channel.fromPath(params.bam_list).splitText()
    
    // Read genomic regions from regions.txt
    regionChannel = Channel.fromPath(params.region_list).splitText()

    // Calculate coverage for each (BAM, region) pair in parallel
    coverageCalc(bamChannel, regionChannel)

    // Merge all coverage outputs into a single file
    mergeCoverage(coverageCalc.out)
}
