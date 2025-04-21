params.bam_list = 'bams.txt'
params.bcf_file = 'cohort.bcf'
params.gff_file = 'annotations.gff'

process splitBcf {
    input:
    val bcf from bcfChannel

    output:
    file "split_*.bcf"

    """
    bcftools index -t ${bcf}
    for chr in \$(bcftools idxstats ${bcf} | awk '{print \$1}' | grep -v '^*'); do
        bcftools view -r \$chr ${bcf} -O b -o split_\${chr}.bcf
    done
    """
}

process integrateData {
    input:
    file bcf_file from splitBcf.out
    val bam_file from bamChannel
    file gff_file from gffChannel

    output:
    file("integrated_${bam_file.baseName}_${bcf_file.baseName}.json")

    """
    rust_integrate_tool \
      --bam ${bam_file} \
      --bcf ${bcf_file} \
      --gff ${gff_file} \
      --out integrated_${bam_file.baseName}_${bcf_file.baseName}.json
    """
}

process mergeIntegrations {
    input:
    file integrated_files from integrateData.out.collect()
    output:
    file "final_integration.json"

    """
    rust_merge_integration ${integrated_files.join(" ")} final_integration.json
    """
}

workflow {
    bcfChannel = Channel.value(params.bcf_file)
    gffChannel = Channel.value(params.gff_file)
    bamChannel = Channel.fromPath(params.bam_list).splitText()

    splitBcf(bcfChannel)
    integrateData(splitBcf.out, bamChannel, gffChannel)
    mergeIntegrations(integrateData.out)
}
