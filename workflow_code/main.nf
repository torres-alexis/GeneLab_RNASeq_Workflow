include { RNASEQ } from './workflows/rnaseq.nf'
include { STAGE } from './subworkflows/stage.nf'

include { CAPTURE_PROCESSING_RECORDS } from './modules/capture_processing_records.nf'
include { PACKAGE_PROCESSING_INFO } from './modules/package_processing_info.nf'
include { GET_ACCESSIONS as GET_POST_ACCESSIONS } from './modules/get_accessions.nf'

include { GENERATE_MD5SUMS } from './modules/generate_md5sums.nf'
include { UPDATE_ASSAY_TABLE } from './modules/update_assay_table.nf'

def dp_tools_plugin() {
    return params.dp_tools_plugin
        ? file(params.dp_tools_plugin)
        : file(params.mode == 'microbes'
            ? "${projectDir}/bin/dp_tools__NF_RCP_Bowtie2"
            : "${projectDir}/bin/dp_tools__NF_RCP")
}

def enabled(value) {
    return value != null && value.toString().toBoolean()
}

def print_banner() {
    if (enabled(params.version)) {
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
    main:
        if ((params.help && params.help.toString() != 'false') || enabled(params.helpFull)) {
            exit 0
        }
        if (enabled(params.stage_only) && enabled(params.post_processing)) {
            error "Specify only one of --stage_only or --post_processing"
        }
        if (enabled(params.stage_only)) {
            STAGE_ONLY()
            ch_pub = STAGE_ONLY.out.published
        } else if (enabled(params.post_processing)) {
            POST_PROCESSING()
            ch_pub = POST_PROCESSING.out.published
        } else {
            print_banner()
            def run = [run_name: workflow.runName, session_id: workflow.sessionId.toString(),
                       command: workflow.commandLine, entry_point: params.entry_point,
                       accession: params.accession, mode: params.mode, assay_suffix: params.assay_suffix,
                       profile: workflow.profile, nextflow_version: nextflow.version.toString(),
                       launch_dir: workflow.launchDir.toString(), project_dir: workflow.projectDir.toString(),
                       work_dir: workflow.workDir.toString(), output_dir: workflow.outputDir.toString()]
            if (!workflow.preview && !workflow.stubRun) WorkflowProvenance.begin(workflow.outputDir, run)
            RNASEQ(
                channel.value(workflow.outputDir),
                channel.value(dp_tools_plugin()),
                channel.value(params.reference_table),
                params.accession ? channel.value(params.accession) : null,
                params.isa_archive ? channel.fromPath(params.isa_archive) : null,
                params.runsheet ? channel.fromPath(params.runsheet) : null,
                channel.value(params.api_url),
                channel.value(params.reference_store_path),
                channel.value(params.derived_store_path)
            )
            ch_pub = RNASEQ.out.published
            if (!workflow.preview && !workflow.stubRun) {
                CAPTURE_PROCESSING_RECORDS(ch_pub.count(), workflow.outputDir, params.assay_suffix, run,
                    RNASEQ.out.processing_samples, RNASEQ.out.processed_directory)
            }
        }
    publish:
        published = ch_pub
    onComplete:
        if (!workflow.preview && !workflow.stubRun &&
            !enabled(params.post_processing) && !enabled(params.stage_only) &&
            (!params.help || params.help.toString() == 'false') && !enabled(params.helpFull) && !enabled(params.version)) {
            try {
                WorkflowProvenance.finish(workflow.outputDir, workflow.runName, workflow.sessionId.toString(),
                    workflow.success && !(workflow.failOnIgnore && workflow.stats.ignoredCount > 0))
            } catch (Exception e) {
                log.error "Could not finalize processing records: ${e.message}"
                nextflow.Global.session.abort(e)
            }
        }
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
            channel.value(workflow.outputDir),
            channel.value(dp_tools_plugin()),
            params.accession ? channel.value(params.accession) : null,
            params.isa_archive ? channel.fromPath(params.isa_archive) : null,
            params.runsheet ? channel.fromPath(params.runsheet) : null,
            channel.value(params.api_url)
        )
        ch_pub = STAGE.out.published
    emit:
        published = ch_pub
}

workflow POST_PROCESSING {
    main:
        print_banner()
        def requested = params.accession?.toString()?.toUpperCase()
        def current_records = workflow.outputDir.resolve('processing_info')
        def current_record = current_records.resolve('workflow_run.json')
        def record = WorkflowProvenance.postRecord(workflow.outputDir, requested)
        def stored_root = !record && current_record.exists()
            ? new groovy.json.JsonSlurper().parseText(current_record.text).output_root?.toString()?.toUpperCase()
            : null
        def different_root = (requested ==~ /GLDS-\d+/) && (stored_root ==~ /GLDS-\d+/) && requested != stored_root
        if (requested && ((!record && current_record.exists() && !different_root) ||
                          (requested.startsWith('OSD-') && !record?.output_root))) {
            GET_POST_ACCESSIONS(channel.value(requested), channel.value(params.api_url))
            ch_run = GET_POST_ACCESSIONS.out.accessions_txt.map { txt ->
                def accessions = txt.readLines().collect { it.trim() }
                def matched = record ?: WorkflowProvenance.postRecord(workflow.outputDir, requested, accessions)
                [matched?.output_root ?: accessions[1], matched]
            }
        } else {
            ch_run = channel.value([record?.output_root ?: requested ?: 'results', record])
        }
        ch_context = ch_run.map { root, matched ->
            def records = matched?.capture_status == 'complete'
                ? current_records.resolve('runs').resolve(matched.run_name as String)
                : current_records
            def package_records = matched || (!current_record.exists() && records.resolve('samples.txt').exists())
            [directory: file("${workflow.outputDir}/${root}", checkIfExists: true), root: root,
             mode: matched?.containsKey('mode') ? matched.mode : params.mode,
             suffix: matched?.containsKey('assay_suffix') ? matched.assay_suffix : params.assay_suffix,
             records: package_records ? records : null]
        }
        ch_root = ch_context.map { it.root }
        UPDATE_ASSAY_TABLE(ch_context.filter { it.root.startsWith('GLDS-') }.map {
            [it.directory, it.root, it.mode, it.suffix]
        })
        GENERATE_MD5SUMS(ch_context.combine(channel.of('raw', 'processed')).map { context, type ->
            [context.directory, type, context.suffix]
        })
        ch_package = ch_context.filter { it.records }
        PACKAGE_PROCESSING_INFO(ch_package.map { it.records }, ch_package.map { it.suffix })
        ch_pub = UPDATE_ASSAY_TABLE.out.assay_table.flatten().combine(ch_root).map { f, root ->
            ["${root}/GeneLab/updated_curation_tables/${f.name}".toString(), f]
        }.mix(
            GENERATE_MD5SUMS.out.md5sum.combine(ch_root).map { f, root -> ["${root}/GeneLab/${f.name}".toString(), f] },
            PACKAGE_PROCESSING_INFO.out.zip.combine(ch_root).map { f, root ->
                ["${root}/GeneLab/${f.name}".toString(), f]
            }
        )
    emit:
        published = ch_pub
}
