process UMITOOLS_EXTRACT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi_tools:1.1.6--py39hbcbf7aa_0':
        'biocontainers/umi_tools:1.1.6--py39hbcbf7aa_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*umi_trimmed.fastq.gz"), emit: reads
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        seqkit replace -p " " -r "_" ${reads} | gzip > ${prefix}_space_removed.fastq.gz
        umi_tools extract \\
            ${args} \\
            --stdin=${prefix}_space_removed.fastq.gz \\
            --stdout=${prefix}_umi_trimmed.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            umitools: \$(umitools --version)
        END_VERSIONS
        """
    }
    else {
        """
        seqkit replace -p " " -r "_" ${reads[0]} | gzip > ${prefix}_R1_space_removed.fastq.gz
        seqkit replace -p " " -r "_" ${reads[1]} | gzip > ${prefix}_R2_space_removed.fastq.gz
        umi_tools extract \\
            ${args} \\
            --stdin=${prefix}_R1_space_removed.fastq.gz --read2-in=${prefix}_R2_space_removed.fastq.gz \\
            --stdout=${prefix}_R1_umi_trimmed.fastq.gz --read2-out=${prefix}_R2_umi_trimmed.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            umitools: \$(umitools --version)
        END_VERSIONS
        """
    }

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if (meta.single_end) {
        echo '' | gzip > ${prefix}_umi_trimmed.fastq.gz
    } else {
        echo '' | gzip > ${prefix}_R1_umi_trimmed.fastq.gz
        echo '' | gzip > ${prefix}_R2_umi_trimmed.fastq.gz
    }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umitools: \$(umitools --version)
    END_VERSIONS
    """
}
