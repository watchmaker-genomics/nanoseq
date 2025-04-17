process RSEQC_GENEBODYCOVERAGE {
    label 'process_high'
    container "quay.io/biocontainers/rseqc:3.0.1--py37h516909a_1"

    input:
    tuple path(bam), path(bai), path(bed12)

    output:
    path("*.pdf")                  , emit: pdf
    path("*.geneBodyCoverage.txt") , emit: rna_txt_ch
    path("versions.yml")           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def name = bam.getName().replaceAll(/\.bam$/, '')

    """
    geneBody_coverage.py \\
        $args \\
        --refgene=$bed12 \\
        --input=$bam  \\
        --minimum_length=100 \\
        --out-prefix=${name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(geneBody_coverage.py --version | sed -e "s/geneBody_coverage.py //g")
    END_VERSIONS
    """
}
