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
        [ ! -f  ${prefix}_UMI.fastq.gz ] && ln -s ${reads} ${prefix}_UMI.fastq.gz
        umi_tools extract \\
            ${args} \\
            --stdin=${prefix}_UMI.fastq.gz \\
            --stdout=${prefix}_temp.fastq
        sed -r '1~4 s/^(@.*)-([ACTGN]{3})(.*)/\\1 \\3:\\2/' ${prefix}_temp.fastq | gzip > ${prefix}_umi_trimmed.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            umitools: \$(echo \$(umi_tools --version) | sed 's/^.*: //')
        END_VERSIONS
        """
    }
    else {
        """
        [ ! -f  ${prefix}_UMI_R1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_UMI_R1.fastq.gz
        [ ! -f  ${prefix}_UMI_R2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_UMI_R2.fastq.gz
        umi_tools extract \\
            ${args} \\
            --stdin=${prefix}_UMI_R1.fastq.gz --read2-in=${prefix}_UMI_R2.fastq.gz \\
            --stdout=${prefix}_R1_temp.fastq --read2-out=${prefix}_R2_temp.fastq
        sed -r '1~4 s/^(@.*)-([ACTGN]{6})(.*)/\\1\\3:\\2/' ${prefix}_R1_temp.fastq | gzip > ${prefix}_R1_umi_trimmed.fastq.gz
        sed -r '1~4 s/^(@.*)-([ACTGN]{6})(.*)/\\1\\3:\\2/' ${prefix}_R2_temp.fastq | gzip > ${prefix}_R2_umi_trimmed.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            umitools: \$(echo \$(umi_tools --version) | sed 's/^.*: //')
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
        umitools: \$(echo \$(umi_tools --version) | sed 's/^.*: //')
    END_VERSIONS
    """
}
