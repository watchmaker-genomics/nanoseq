process RESTRANDER {
    tag "$meta.id"
    label 'process_medium'

    container "${'912684371407.dkr.ecr.us-west-2.amazonaws.com/restrander:1.2'}"

    input:
    tuple val(meta), path(reads), path(input_config)

    output:
    tuple val(meta), path("*_restrander.fq.gz"), emit: reads
    tuple val(meta), path("*_restrander-unknowns.fq.gz"), emit: unknown_reads
    tuple val(meta), path("*.restrander.json"), emit: metrics
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // _restrander-unknowns.fq.gz

    script:
    """
    prefix=\${task.ext.prefix:-${meta.id}}

    /restrander \\
        ${reads} \\
        \${prefix}_restrander.fq.gz \\
        ${input_config} > \${prefix}.restrander.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        restrander: v1.0.1
    END_VERSIONS
    """
}
