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

workflow {
    if (params.stage_only && params.post_processing) {
        error "Specify only one of --stage_only or --post_processing"
    }
    if (params.stage_only) {
        STAGE_ONLY()
    } else if (params.post_processing) {
        POST_PROCESSING()
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
        PUBLISH_STAGED_ANALYSIS(
            STAGE.out.ch_outdir,
            STAGE.out.runsheet_path,
            STAGE.out.raw_reads.map { sample -> sample[1] }.collect()
        )
}

workflow POST_PROCESSING {
    main:
        print_banner()
        ch_processed_directory = channel.fromPath("${params.outdir}/${params.accession}", checkIfExists: true)
        UPDATE_ASSAY_TABLE(ch_processed_directory)
        GENERATE_MD5SUMS(ch_processed_directory)
}
