process UPDATE_ASSAY_TABLE {
    cache false

    input:
        tuple path(ch_processed_directory), val(accession), val(mode), val(suffix)
    output:
        path("a_*.txt"), emit: assay_table, optional: true

    script:
    def mode_param = mode == "microbes" ? "--mode microbes" : ""
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def assay_suffix_arg = suffix ? "--assay_suffix ${quote.call(suffix)}" : ""
    """
    update_assay_table.py --outdir ${quote.call(ch_processed_directory)} ${mode_param} ${assay_suffix_arg} --glds_accession ${quote.call(accession)}
    """
}
