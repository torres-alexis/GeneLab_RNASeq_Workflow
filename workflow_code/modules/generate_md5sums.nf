process GENERATE_MD5SUMS {
    cache false

    tag { type }

    input:
        tuple path(ch_outdir), val(type), val(suffix)
    output:
        path("*md5sum${suffix}.tsv"), emit: md5sum

    script:
        def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
        def assay_suffix_arg = suffix ? "--assay_suffix ${quote.call(suffix)}" : ""
        """
        generate_md5sums.py --outdir ${quote.call(ch_outdir)} ${assay_suffix_arg} --${type}
        """
}
