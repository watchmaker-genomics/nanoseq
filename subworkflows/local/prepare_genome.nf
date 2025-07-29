/*
 * Prepare genome/transcriptome before alignment
 */

include { GET_CHROM_SIZES  } from '../../modules/local/get_chrom_sizes'
include { GTF2BED          } from '../../modules/local/gtf2bed'
include { SAMTOOLS_FAIDX   } from '../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_GENOME {
    take:
    ch_fastq

    main:
    // Get unique list of all fasta files - reference FASTA is at position 2
    ch_fastq
        .filter { it[2] }
        .map { it -> [ it[2], it[6].toString() ] }  // [ fasta, annotation_str ]
        .unique()
        .set { ch_fastq_sizes }

    /*
     * Make chromosome sizes file
     */
    GET_CHROM_SIZES ( ch_fastq_sizes )
    ch_chrom_sizes = GET_CHROM_SIZES.out.sizes
    samtools_version = GET_CHROM_SIZES.out.versions

    // Get unique list of all gtf files - GTF is at position 3 in the tuple
    ch_fastq
        .filter { it[3] }
        .map { it -> [ it[3], it[6] ] }  // [ gtf, annotation_str ]
        .unique()
        .set { ch_fastq_gtf }

    /*
     * Convert GTF to BED12
     */
    GTF2BED ( ch_fastq_gtf )
    ch_gtf_bed = GTF2BED.out.gtf_bed
    gtf2bed_version = GTF2BED.out.versions

    ch_chrom_sizes
        .join(ch_gtf_bed, by: 1, remainder:true)
        .map { it -> [ it[1], it[2], it[0] ] }
        .cross(ch_fastq) { it -> it[-1] }
        .flatten()
        .collate(10)
        .map { it -> [ it[5], it[0], it[6], it[1], it[8], it[9] ]} // [ fasta, sizes, gtf, bed, is_transcripts, annotation_str ]
        .unique()
        .set { ch_fasta_index }

    /*
     * Convert GTF to BED12
     */
    ch_fastq
        .filter { it[2] }
        .map { it -> [ it[0], it[2] ] }  // [ meta, fasta ]
        .unique()
        .set { ch_fasta }

    SAMTOOLS_FAIDX ( ch_fasta )
    ch_fai = SAMTOOLS_FAIDX.out.fai

    emit:
    ch_fasta
    ch_fai
    ch_fasta_index
    ch_gtf_bed
    samtools_version
    gtf2bed_version
}
