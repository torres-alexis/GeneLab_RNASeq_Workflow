process PUBLISH_STAGED_ANALYSIS {
    publishDir path: { "${ch_outdir}" },
        pattern: '{00-RawData/**,01-TG_Preproc/**,02-Bowtie2_Alignment/**,02-STAR_Alignment/**,03-FeatureCounts/**,03-RSEM_Counts/**,05-DESeq2_DGE/**,Metadata/**}',
        mode: params.publish_dir_mode

    input:
        val(ch_outdir)
        path(runsheet)
        path(staged)
        val(entry_point)

    output:
        path("Metadata/*"), emit: metadata

    script:
        def suf = params.assay_suffix
        def dest
        def sample_cut = ""
        if (entry_point == 'raw_reads') {
            dest = '00-RawData/Fastq'
        } else if (entry_point == 'trimmed_reads') {
            dest = '01-TG_Preproc/Fastq'
        } else if (entry_point == 'bam_files') {
            dest = params.mode == 'microbes' ? '02-Bowtie2_Alignment' : '02-STAR_Alignment'
            sample_cut = params.mode == 'microbes' ? "${suf}.bam" : "${suf}_Aligned.toTranscriptome.out.bam"
        } else if (entry_point == 'genes_results') {
            dest = '03-RSEM_Counts'
            sample_cut = "${suf}.genes.results"
        } else if (entry_point == 'counts_table') {
            dest = params.mode == 'microbes' ? '03-FeatureCounts' : '03-RSEM_Counts'
        } else if (entry_point == 'dge_table') {
            dest = '05-DESeq2_DGE'
        }
        if (sample_cut) {
            """
            mkdir -p Metadata
            mv ${runsheet} Metadata/
            for f in ${staged}; do
                b=\$(basename "\$f")
                s=\${b%${sample_cut}}
                mkdir -p ${dest}/\$s
                mv "\$f" ${dest}/\$s/
            done
            """
        } else {
            """
            mkdir -p Metadata ${dest}
            mv ${runsheet} Metadata/
            mv ${staged} ${dest}/
            """
        }
}
