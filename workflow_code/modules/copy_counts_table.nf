process COPY_COUNTS_TABLE {

    input:
        path("?.csv")
        path("samples.txt")

    output:
        path("*Unnormalized_Counts*.csv"), emit: counts_table

    script:
        def counts_name = params.mode == "microbes" ? 
            "FeatureCounts_Unnormalized_Counts${params.assay_suffix}.csv" :
            "RSEM_Unnormalized_Counts${params.assay_suffix}.csv"
        """
        order_table_samples.py --input 1.csv --output ${counts_name} --samples samples.txt
        """
}

