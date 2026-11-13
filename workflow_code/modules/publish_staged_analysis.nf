/*
    This process publishes the staged raw reads to the output directory.
    This is only used for the stage-only workflow.
*/

process PUBLISH_STAGED_ANALYSIS {
    publishDir path: { "${ch_outdir}" },
    pattern: '00-RawData/**',
    mode: params.publish_dir_mode

    input:
        val(ch_outdir)
        path(raw_reads)
    
    output:
        path("00-RawData/Fastq/*"), emit: raw_reads
    
    script:
        """
        # Create directory structure
        mkdir -p 00-RawData/Fastq
        
        # Move raw reads to 00-RawData/Fastq directory
        mv ${raw_reads} 00-RawData/Fastq/

        """
} 