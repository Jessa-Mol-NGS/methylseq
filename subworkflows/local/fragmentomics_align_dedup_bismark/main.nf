// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join

include { BISMARK_ALIGN                                  } from '../../../modules/nf-core/bismark/align/main'
include { BISMARK_DEDUPLICATE                            } from '../../../modules/nf-core/bismark/deduplicate/main'
include { SAMTOOLS_SORT                                  } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX                                 } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_NAME            } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_PROPERLY_PAIRED } from '../../../modules/nf-core/samtools/view/main'      
include { FRAGMENTOMICS                                  } from '../../../modules/local/fragmentomics/main'
include { BISMARK_REPORT                                 } from '../../../modules/nf-core/bismark/report/main'
include { BISMARK_SUMMARY                                } from '../../../modules/nf-core/bismark/summary/main'

workflow FRAGMENTOMICS_ALIGN_DEDUP_BISMARK {

    take:
    ch_reads             // channel: [ val(meta), [ reads ] ]
    ch_fasta             // channel: [ val(meta), [ fasta ] ]
    ch_bismark_index     // channel: [ val(meta), [ bismark index ] ]
    skip_deduplication   // boolean: whether to deduplicate alignments

    main:
    ch_alignments                 = Channel.empty()
    ch_alignment_reports          = Channel.empty()
    ch_bismark_report             = Channel.empty()
    ch_bismark_summary            = Channel.empty()
    ch_insert_sizes               = Channel.empty()
    ch_multiqc_files              = Channel.empty()
    ch_versions                   = Channel.empty()

    /*
     * Align with bismark
     */
    BISMARK_ALIGN (
        ch_reads,
        ch_fasta,
        ch_bismark_index
    )
    ch_alignments        = BISMARK_ALIGN.out.bam
    ch_alignment_reports = BISMARK_ALIGN.out.report.map{ meta, report -> [ meta, report, [] ] }
    ch_versions = ch_versions.mix(BISMARK_ALIGN.out.versions)

    if (!skip_deduplication) {
        /*
        * Run deduplicate_bismark
        */
        BISMARK_DEDUPLICATE (
            BISMARK_ALIGN.out.bam
        )
        ch_alignments        = BISMARK_DEDUPLICATE.out.bam
        ch_alignment_reports = BISMARK_ALIGN.out.report.join(BISMARK_DEDUPLICATE.out.report)
        ch_versions          = ch_versions.mix(BISMARK_DEDUPLICATE.out.versions)
    }

    /*
    * Generate bismark sample reports
    */
    BISMARK_REPORT (
        ch_alignment_reports.map { meta, report1, report2 -> [ meta, report1, report2, [], [] ] }
    )
    ch_bismark_report = BISMARK_REPORT.out.report
    ch_versions       = ch_versions.mix(BISMARK_REPORT.out.versions)
    
    /*
    * Generate bismark summary report
    */
    BISMARK_SUMMARY (
        BISMARK_ALIGN.out.bam.collect{ meta, bam -> bam.name },
        ch_alignment_reports.collect{ meta, align_report, dedup_report -> align_report },
        ch_alignment_reports.collect{ meta, align_report, dedup_report -> dedup_report }.ifEmpty([]),
        [], []
        )
    ch_bismark_summary = BISMARK_SUMMARY.out.summary
    ch_versions        = ch_versions.mix(BISMARK_SUMMARY.out.versions)

    /*
     * MODULE: Run samtools sort on aligned or deduplicated bam
     */
    SAMTOOLS_SORT (
        ch_alignments,
        [[:],[]] // [ [meta], [fasta]]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    /*
     * MODULE: Run samtools index on aligned or deduplicated bam
     */
    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    /*
    * Remaining fragmentomics workflow is only intended for paired-end reads
    */
    ch_alignments_branched = ch_alignments
        .branch { meta, bam ->
                single_end : meta.single_end
                    return [ meta, bam ]
                paired_end : !meta.single_end
                    return [ meta, bam ]
            }
    /*
    * Sort the BAM file by read name for insert size calculation
    */
    SAMTOOLS_SORT_NAME (
        ch_alignments_branched.paired_end,
        [[:], []] // [ [meta], [fasta]]        
    )

    /*
    * Only keep properly paired reads for insert size calculation
    */
    SAMTOOLS_VIEW_PROPERLY_PAIRED (
        SAMTOOLS_SORT_NAME.out.bam.map { meta, bam -> [ meta, bam, [] ] },
        [[:], []], // [ [meta2], [fasta] ]
        [], // [qname]
        "" // [index_format]        
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_PROPERLY_PAIRED.out.versions)

    /* 
    * Extract insert sizes and sequences 
    */
    FRAGMENTOMICS (
        SAMTOOLS_VIEW_PROPERLY_PAIRED.out.bam
    )
    ch_insert_sizes = FRAGMENTOMICS.out.insert_sizes
    ch_versions = ch_versions.mix(FRAGMENTOMICS.out.versions)    
    
    /*
    * Collect MultiQC inputs
    */
    ch_multiqc_files = ch_bismark_summary
                            .mix(ch_alignment_reports.collect{ meta, align_report, dedup_report -> align_report })
                            .mix(ch_alignment_reports.collect{ meta, align_report, dedup_report -> dedup_report })
    
    emit:
    bam                        = SAMTOOLS_SORT.out.bam         // channel: [ val(meta), [ bam ] ]
    bai                        = SAMTOOLS_INDEX.out.bai        // channel: [ val(meta), [ bai ] ]
    bismark_report             = ch_bismark_report             // channel: [ val(meta), [ report ] ]
    bismark_summary            = ch_bismark_summary            // channel: [ val(meta), [ summary ] ]
    insert_sizes               = ch_insert_sizes               // channel: [ val(meta), [ csv.gz ] ]
    multiqc                    = ch_multiqc_files              // path: *{html,txt}
    versions                   = ch_versions                   // path: *.version.txt

}
