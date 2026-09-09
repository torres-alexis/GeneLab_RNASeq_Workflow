process COPY_BAMS {
    tag "Sample: ${meta.id}"
    
    publishDir path: {
        def base = params.mode == "microbes" ? "02-Bowtie2_Alignment" : "02-STAR_Alignment"
        return "${publishdir}/${base}/${meta.id}"
    },
        mode: params.publish_dir_mode

    input:
        val(publishdir)
        tuple val(meta), path("?.bam")

    output:
        tuple val(meta), path("${meta.id}*.bam"), emit: bam_files

    script:
        def bam_name = params.mode == "microbes" ? 
            "${meta.id}${params.assay_suffix}.bam" :
            "${meta.id}${params.assay_suffix}_Aligned.toTranscriptome.out.bam"
        """
        cp -P 1.bam ${bam_name}
        """
}

