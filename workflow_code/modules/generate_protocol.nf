process GENERATE_PROTOCOL {
    input:
        val(ch_meta)
        val(strandedness)
        path(software_versions_yaml)
        val(reference_source)
        val(reference_version)
        tuple(path(reference_fasta), path(reference_gtf))
        path(runsheet)
        path(qc_metrics)

    output:
        path("protocol${params.assay_suffix}.txt"), emit: protocol

    script:
        def mode = params.mode == 'microbes' ? '--mode microbes' : ''
        def ref_source = reference_source ? "--reference_source ${reference_source}" : ''
        def ref_version = reference_version ? "--reference_version ${reference_version}" : ''
        def assay_suffix_arg = params.assay_suffix ? "--assay_suffix ${params.assay_suffix}" : ''
        def qc_arg = qc_metrics.name != "NO_FILE" ? "--qc_metrics ${qc_metrics}" : ''
        def drop_on = params.drop_unalignable && params.mode != 'microbes' && params.entry_point in ['raw_reads', 'trimmed_reads', 'bam_files']

        """
        generate_protocol.py \
        ${mode} \
        --outdir . \
        --software_table ${software_versions_yaml} \
        ${assay_suffix_arg} \
        --paired_end ${ch_meta.paired_end} \
        --has_ercc ${ch_meta.has_ercc} \
        --workflow_version ${workflow.manifest.version} \
        --strandedness ${strandedness} \
        --organism "${ch_meta.organism_sci}" \
        ${ref_source} \
        ${ref_version} \
        --reference_fasta ${reference_fasta} \
        --reference_gtf ${reference_gtf} \
        --runsheet ${runsheet} \
        --drop_unalignable ${drop_on} \
        --unalignable_threshold ${params.unalignable_threshold} \
        ${qc_arg}
        """
}
