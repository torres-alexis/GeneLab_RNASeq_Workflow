process COPY_DGE_TABLE {

    input:
        path("?.csv")
        path("samples.txt")
        path("runsheet.csv")

    output:
        path("differential_expression${params.assay_suffix}.csv"), emit: dge_table

    script:
        def dge_table_name = "differential_expression${params.assay_suffix}.csv"
        """
        order_table_samples.py --input 1.csv --output ${dge_table_name} --samples samples.txt --runsheet runsheet.csv
        """
}

