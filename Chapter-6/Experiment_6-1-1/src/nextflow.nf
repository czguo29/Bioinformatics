nextflow.enable.dsl=2
process coverageAnalysis {
    tag "${region}"
    input:
    val region
    output:
    file("coverage_${region.replace(':', '_')}.txt")
    script:
    """
    /mnt/c/Users/trian/BGVR/chapter_06/experiment_61/target/debug/coverage_tool \
        --bam input.bam \
        --region ${region} \
        --output coverage_${region.replace(':', '_')}.txt
    """
}
workflow {
    // Define the channel inside the workflow
    def region_list = Channel.from("1:1-1000000", "1:1000001-2000000")
    // Execute the process with the channel as input
    coverageAnalysis(region_list)
}
