include { COPY_READS } from '../modules/copy_reads.nf'
include { FETCH_REMOTE_READS } from '../modules/fetch_remote.nf'

def is_remote_uri(p) {
    def s = p == null ? "" : p.toString().trim()
    return s.contains("://")
}

workflow STAGE_READS {
    take:
        ch_outdir
        ch_samples
        type

    main:
        ch_split = ch_samples.branch { meta, files ->
            remote: files && files.size() > 0 && is_remote_uri(files[0])
            local: true
        }

        FETCH_REMOTE_READS( ch_split.remote )
        ch_remote = FETCH_REMOTE_READS.out.reads.map { meta, fq ->
            [meta, fq instanceof List ? fq : [fq]]
        }

        ch_local = ch_split.local.map { meta, files ->
            meta.paired_end ? [meta, [files[0], files[1]]] : [meta, [files[0]]]
        }

        COPY_READS( ch_outdir, ch_remote.mix(ch_local), type )

        COPY_READS.out.reads | map { sample -> sample[1] } | collect | set { ch_all_reads }
        COPY_READS.out.reads | map { sample -> sample[0].id }
                        | collectFile(name: "samples.txt", sort: true, newLine: true)
                        | set { samples_txt }

    emit:
        reads = COPY_READS.out.reads
        ch_all_reads = ch_all_reads
        samples_txt = samples_txt
}
