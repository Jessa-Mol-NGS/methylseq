process FRAGMENTOMICS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/04/04b086608e967caf920b0af1b20c388933a2e600bb657fc4c4a119aeb15e6b6a/data':
        'community.wave.seqera.io/library/pysam_pandas:baaaf9c9ffb0e943' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*insert.sizes.csv.gz"), emit: insert_sizes
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    template 'extract-insert-sizes.py' 
   
    stub:
    def args = task.ext.args ?: ''
    """   
    touch ${meta.id}.insert.sizes.csv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( echo \$(python --version) | sed 's/^Python //' )
    END_VERSIONS
    """
}
