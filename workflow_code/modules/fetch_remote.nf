process FETCH_REMOTE_READS {
    tag "Sample: ${meta.id}"

    input:
        tuple val(meta), val(uris)

    output:
        tuple val(meta), path("?.gz"), emit: reads

    script:
        def u2 = (uris.size() > 1) ? uris[1] : ""
        def nrec = params.truncate_to ? params.truncate_to.toInteger() : 0
        if ( meta.paired_end && nrec ) {
        """
        fastq_first_n.sh ${nrec} 1.gz '${uris[0]}'
        fastq_first_n.sh ${nrec} 2.gz '${u2}'
        """
        } else if ( meta.paired_end ) {
        """
        fetch_uri.sh '${uris[0]}' 1.gz
        fetch_uri.sh '${u2}' 2.gz
        """
        } else if ( nrec ) {
        """
        fastq_first_n.sh ${nrec} 1.gz '${uris[0]}'
        """
        } else {
        """
        fetch_uri.sh '${uris[0]}' 1.gz
        """
        }
}

process FETCH_REMOTE_BAM {
    tag "Sample: ${meta.id}"

    input:
        tuple val(meta), val(uris)

    output:
        tuple val(meta), path("?.bam"), emit: bam_files

    script:
        """
        fetch_uri.sh '${uris[0]}' 1.bam
        """
}

process FETCH_REMOTE_GENES_RESULTS {
    tag "Sample: ${meta.id}"

    input:
        tuple val(meta), val(uris)

    output:
        tuple val(meta), path("?.genes.results"), emit: genes_results

    script:
        """
        fetch_uri.sh '${uris[0]}' 1.genes.results
        """
}
