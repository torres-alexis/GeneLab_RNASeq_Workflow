process PARSE_QC_METRICS {
    input:
        val(osd_accession)
        val(meta)
        path(isa_zip)
        path(all_multiqc_output)
        path(rsem_counts)
        path(runsheet_path)

    output:
        path("qc_metrics${params.assay_suffix}.csv"), emit: file
        path("isa_archive/*"), optional: true, emit: metadata_copy

    script:
        def assay_suffix_arg = params.assay_suffix ? "--assay_suffix ${params.assay_suffix}" : ""
        def mode_param = params.mode ? "--mode ${params.mode}" : ""
        def osd_arg = osd_accession ? "--osd-num ${ osd_accession }" : ""
        """
        mkdir -p isa_archive
        if [ -e '${isa_zip}' ]; then
            cp -L '${isa_zip}' isa_archive/
        fi
        cp -L '${runsheet_path}' isa_archive/
        parse_multiqc.py ${ osd_arg } ${ assay_suffix_arg } ${ meta.paired_end ? '--paired' : '' } ${ mode_param } --runsheet ${ runsheet_path }
        """
}
