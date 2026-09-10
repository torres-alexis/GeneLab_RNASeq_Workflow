include { RNASEQ } from './workflows/rnaseq.nf'
include { STAGE } from './subworkflows/stage.nf'

include { GENERATE_MD5SUMS } from './modules/generate_md5sums.nf'
include { UPDATE_ASSAY_TABLE } from './modules/update_assay_table.nf'
include { PUBLISH_STAGED_ANALYSIS } from './modules/publish_staged_analysis.nf'

def dp_tools_plugin() {
    return params.dp_tools_plugin
        ? file(params.dp_tools_plugin)
        : file(params.mode == 'microbes'
            ? "${projectDir}/bin/dp_tools__NF_RCP_Bowtie2"
            : "${projectDir}/bin/dp_tools__NF_RCP")
}

def print_banner() {
    if (params.version) {
        println """${workflow.manifest.name}
Workflow Version: ${workflow.manifest.version}"""
        exit 0
    }

    println """
${workflow.manifest.name}
Workflow Version: ${workflow.manifest.version}
""".stripIndent()

    if (params.limit_samples_to || params.truncate_to || params.force_single_end || params.genome_subsample) {
        println("WARNING: Debugging options enabled!")
        println("Sample limit: ${params.limit_samples_to ?: 'Not set'}")
        println("Read truncation: ${params.truncate_to ? "First ${params.truncate_to} records" : 'Not set'}")
        println("Reference genome subsampling: ${params.genome_subsample ? "Region '${params.genome_subsample}'" : 'Not set'}")
        println("Force single-end analysis: ${params.force_single_end ? 'Yes' : 'No'}")
    } else {
        println("No debugging options enabled")
    }
}

def as_list(x) {
    if ( x == null ) return []
    if ( x instanceof Collection && !(x instanceof CharSequence) ) return x as List
    return [x]
}

workflow {
    main:
        if (params.stage_only && params.post_processing) {
            error "Specify only one of --stage_only or --post_processing"
        }
        if (params.stage_only) {
            STAGE_ONLY()
            ch_pub = STAGE_ONLY.out.published
        } else if (params.post_processing) {
            POST_PROCESSING()
            ch_pub = POST_PROCESSING.out.published
        } else {
            print_banner()
            RNASEQ(
                params.outdir ? channel.fromPath(params.outdir, checkIfExists: true) : null,
                channel.value(dp_tools_plugin()),
                channel.value(params.reference_table),
                params.accession ? channel.value(params.accession) : null,
                params.isa_archive_path ? channel.fromPath(params.isa_archive_path) : null,
                params.runsheet_path ? channel.fromPath(params.runsheet_path) : null,
                channel.value(params.api_url),
                channel.value(params.reference_store_path),
                channel.value(params.derived_store_path)
            )
            ch_pub = RNASEQ.out.published
        }
    publish:
        published = ch_pub
}

output {
    published {
        path { dest, f -> f >> dest }
    }
}

workflow STAGE_ONLY {
    main:
        print_banner()
        STAGE(
            params.outdir ? channel.fromPath(params.outdir, checkIfExists: true) : null,
            channel.value(dp_tools_plugin()),
            params.accession ? channel.value(params.accession) : null,
            params.isa_archive_path ? channel.fromPath(params.isa_archive_path) : null,
            params.runsheet_path ? channel.fromPath(params.runsheet_path) : null,
            channel.value(params.api_url)
        )
        def ep = params.entry_point
        if (ep == 'raw_reads') {
            ch_staged = STAGE.out.raw_reads.map { _meta, files -> files }.flatten()
        } else if (ep == 'trimmed_reads') {
            ch_staged = STAGE.out.trimmed_reads.map { _meta, files -> files }.flatten()
        } else if (ep == 'bam_files') {
            ch_staged = STAGE.out.bam_files.map { _meta, files -> files }.flatten()
        } else if (ep == 'genes_results') {
            ch_staged = STAGE.out.genes_results.map { _meta, files -> files }.flatten()
        } else if (ep == 'counts_table') {
            ch_staged = STAGE.out.counts_table
        } else if (ep == 'dge_table') {
            ch_staged = STAGE.out.dge_table
        }
        PUBLISH_STAGED_ANALYSIS(
            STAGE.out.runsheet_path,
            ch_staged.collect(),
            channel.value(ep)
        )
        ch_root = params.accession ? STAGE.out.glds_accession : channel.value('results')
        ch_pub = STAGE.out.published
        if ( params.runsheet_path ) {
            ch_pub = ch_pub.mix(
                PUBLISH_STAGED_ANALYSIS.out.metadata_dir.combine(ch_root).map { d, root ->
                    ["${root}/Metadata".toString(), d]
                }
            )
        }
    emit:
        published = ch_pub
}

workflow POST_PROCESSING {
    main:
        print_banner()
        ch_processed_directory = channel.fromPath("${params.outdir}/${params.accession}", checkIfExists: true)
        UPDATE_ASSAY_TABLE(ch_processed_directory)
        GENERATE_MD5SUMS(ch_processed_directory.combine(channel.of('raw', 'processed')))
        def root = params.accession
        ch_pub = UPDATE_ASSAY_TABLE.out.assay_table.flatMap { f ->
            as_list(f).collect { x -> ["${root}/GeneLab/updated_curation_tables/${x.name}".toString(), x] }
        }.mix(
            GENERATE_MD5SUMS.out.md5sum.map { f -> ["${root}/GeneLab/${f.name}".toString(), f] }
        )
    emit:
        published = ch_pub
}
