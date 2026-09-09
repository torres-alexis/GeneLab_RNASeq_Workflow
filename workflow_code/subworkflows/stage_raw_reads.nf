include { COPY_READS } from '../modules/copy_reads.nf'

// Copy (and optionally truncate) FASTQs. type is "raw" or "trimmed".
workflow STAGE_READS {
    take:
        ch_outdir
        ch_samples
        type

    main:
        truncate_to = params.truncate_to
        if ( truncate_to ) {
            ch_samples | map { sample -> sample[0].paired_end ? [sample[0], sample[1][0], sample[1][1]] : [sample[0], sample[1][0]]}
                 | branch { row ->
                   paired: row.size() == 3
                   single: row.size() == 2
                 }
                 | set{ ch_read_pointers }

            // TO DO: Move splitFastq into a worker process. Driver + S3 hangs on NF < 25.10.6.

            ch_read_pointers.paired | splitFastq(pe: true, decompress: true, compress: true, limit: truncate_to, by: truncate_to, file: true)
                                    | map { sample -> [ sample[0], [ sample[1], sample[2] ] ]}
                                    | set { ch_reads }
            ch_read_pointers.single | splitFastq(decompress: true, compress: true, limit: truncate_to, by: truncate_to, file: true)
                                    | map { sample -> [ sample[0], [ sample[1] ] ]}
                                    | mix( ch_reads )
                                    | set { ch_reads }

            COPY_READS(ch_outdir, ch_reads, type)
        } else {
            ch_samples | map { sample -> sample[0].paired_end ? [sample[0], [ sample[1][0], sample[1][1] ]] : [sample[0], [sample[1][0]]]}
                         | set { ch_reads }

            COPY_READS(ch_outdir, ch_reads, type)
        }

        COPY_READS.out.reads | map { sample -> sample[1] } | collect | set { ch_all_reads }
        COPY_READS.out.reads | map { sample -> sample[0].id }
                        | collectFile(name: "samples.txt", sort: true, newLine: true)
                        | set { samples_txt }

    emit:
        reads = COPY_READS.out.reads
        ch_all_reads = ch_all_reads
        samples_txt = samples_txt
}
