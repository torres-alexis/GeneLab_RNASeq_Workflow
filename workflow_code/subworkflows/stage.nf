include { PARSE_RUNSHEET } from './parse_runsheet.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'
include { GET_ACCESSIONS } from '../modules/get_accessions.nf'
include { STAGE_READS } from './stage_raw_reads.nf'
include { GET_OSDR_FILE_LIST } from '../modules/get_osdr_file_list.nf'
include { DOWNLOAD_OSDR_READS } from '../modules/download_osdr_reads.nf'
include { DOWNLOAD_OSDR_BAM } from '../modules/download_osdr_bam.nf'
include { DOWNLOAD_OSDR_GENES_RESULTS } from '../modules/download_osdr_genes_results.nf'
include { DOWNLOAD_OSDR_COUNTS_TABLE } from '../modules/download_osdr_counts_table.nf'
include { DOWNLOAD_OSDR_DGE_TABLE } from '../modules/download_osdr_dge_table.nf'
include { COPY_BAMS } from '../modules/copy_bams.nf'
include { COPY_GENES_RESULTS } from '../modules/copy_genes_results.nf'
include { COPY_COUNTS_TABLE } from '../modules/copy_counts_table.nf'
include { COPY_DGE_TABLE } from '../modules/copy_dge_table.nf'
include { validateParameters } from 'plugin/nf-schema'

def ENTRY_PARSE_TYPE = [
    raw_reads:     'raw',
    trimmed_reads: 'trimmed',
    bam_files:     'bam',
    genes_results: 'genes',
    counts_table:  'counts',
    dge_table:     'dge'
]

def samples_txt_from(ch_pairs) {
    ch_pairs | map { it[0].id } | collectFile(name: "samples.txt", sort: true, newLine: true)
}

/**
 * One staging path for every entry point.
 *
 * Accession / ISA / outdir is shared. Then:
 *   - OSDR pull if accession-only (no user runsheet) and entry != raw_reads
 *     (raw uses ISA runsheet paths; later entries download the derived artifact)
 *   - else parse the entry's columns and COPY / STAGE_READS
 *   - counts/dge also honor params.counts_table_path / params.dge_table_path
 *
 * file() staging stays as-is. Worker scheme-detect (proteomics STAGE_INPUT) is a later branch.
 */
workflow STAGE {
    take:
        ch_outdir
        dp_tools_plugin
        accession
        isa_archive_path
        runsheet_path
        api_url

    main:
        def ep = params.entry_point
        def parse_type = ENTRY_PARSE_TYPE[ep]
        if (!parse_type) {
            error "Unknown entry_point '${ep}'. Expected one of: ${ENTRY_PARSE_TYPE.keySet()}"
        }

        Channel.empty() | set { osd_accession }
        Channel.empty() | set { glds_accession }

        if ( accession ) {
            GET_ACCESSIONS( accession, api_url )
            osd_accession = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[0].trim() }
            glds_accession = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[1].trim() }
            ch_outdir = ch_outdir.combine(glds_accession).map { outdir, glds -> "$outdir/$glds" }
        } else {
            ch_outdir = ch_outdir.map { it + "/results" }
        }
        ch_outdir = ch_outdir.first()

        Channel.empty() | set { isa_archive }
        if ( runsheet_path == null ) {
            if ( isa_archive_path == null ) {
                FETCH_ISA( ch_outdir, osd_accession, glds_accession )
                isa_archive = FETCH_ISA.out.isa_archive
            } else {
                isa_archive = isa_archive_path
            }
            ISA_TO_RUNSHEET( ch_outdir, osd_accession, glds_accession, isa_archive, dp_tools_plugin )
            runsheet_path = ISA_TO_RUNSHEET.out.runsheet
        } else if ( isa_archive_path != null ) {
            isa_archive = isa_archive_path
        }

        if ( params.validate_params ) {
            validateParameters()
        }

        Channel.empty() | set { samples }
        Channel.empty() | set { samples_txt }
        Channel.empty() | set { raw_reads }
        Channel.empty() | set { trimmed_reads }
        Channel.empty() | set { bam_files }
        Channel.empty() | set { genes_results }
        Channel.empty() | set { counts_table }
        Channel.empty() | set { dge_table }

        def no_runsheet = (params.runsheet_path == null)
        def osdr_derived = no_runsheet && ep != 'raw_reads'
        def counts_override = (ep == 'counts_table' && params.counts_table_path)
        def dge_override = (ep == 'dge_table' && params.dge_table_path)

        // ISA / --counts_table_path / --dge_table_path: metadata only. Don't require raw cols.
        def type = (osdr_derived || counts_override || dge_override) ? 'meta' : parse_type
        PARSE_RUNSHEET( runsheet_path, Channel.value(type) )
        runsheet_path = PARSE_RUNSHEET.out.runsheet

        if ( counts_override ) {
            samples = PARSE_RUNSHEET.out.samples.map { meta, files -> meta }
            COPY_COUNTS_TABLE(ch_outdir, file(params.counts_table_path))
            counts_table = COPY_COUNTS_TABLE.out.counts_table
        } else if ( dge_override ) {
            samples = PARSE_RUNSHEET.out.samples.map { meta, files -> meta }
            COPY_DGE_TABLE(ch_outdir, file(params.dge_table_path))
            dge_table = COPY_DGE_TABLE.out.dge_table
        } else if ( osdr_derived ) {
            GET_OSDR_FILE_LIST( osd_accession )
            if ( ep == 'trimmed_reads' ) {
                samples = PARSE_RUNSHEET.out.samples
                DOWNLOAD_OSDR_READS(
                    ch_outdir,
                    osd_accession,
                    glds_accession,
                    samples.map { meta, files -> meta },
                    "trimmed",
                    GET_OSDR_FILE_LIST.out.file_list
                )
                trimmed_reads = DOWNLOAD_OSDR_READS.out.trimmed_reads
                samples_txt = samples_txt_from(trimmed_reads)
            } else if ( ep == 'bam_files' ) {
                samples = PARSE_RUNSHEET.out.samples
                DOWNLOAD_OSDR_BAM(
                    ch_outdir,
                    osd_accession,
                    glds_accession,
                    samples.map { meta, files -> meta },
                    GET_OSDR_FILE_LIST.out.file_list
                )
                bam_files = DOWNLOAD_OSDR_BAM.out.bam_files
                samples_txt = samples_txt_from(bam_files)
            } else if ( ep == 'genes_results' ) {
                samples = PARSE_RUNSHEET.out.samples
                DOWNLOAD_OSDR_GENES_RESULTS(
                    ch_outdir,
                    osd_accession,
                    glds_accession,
                    samples.map { meta, files -> meta },
                    GET_OSDR_FILE_LIST.out.file_list
                )
                genes_results = DOWNLOAD_OSDR_GENES_RESULTS.out.genes_results
                samples_txt = samples_txt_from(genes_results)
            } else if ( ep == 'counts_table' ) {
                samples = PARSE_RUNSHEET.out.samples.map { meta, files -> meta }
                DOWNLOAD_OSDR_COUNTS_TABLE(
                    ch_outdir,
                    osd_accession,
                    glds_accession,
                    GET_OSDR_FILE_LIST.out.file_list
                )
                counts_table = DOWNLOAD_OSDR_COUNTS_TABLE.out.counts_table
            } else if ( ep == 'dge_table' ) {
                samples = PARSE_RUNSHEET.out.samples.map { meta, files -> meta }
                DOWNLOAD_OSDR_DGE_TABLE(
                    ch_outdir,
                    osd_accession,
                    glds_accession,
                    GET_OSDR_FILE_LIST.out.file_list
                )
                dge_table = DOWNLOAD_OSDR_DGE_TABLE.out.dge_table
            }
        } else if ( ep == 'raw_reads' ) {
            samples = PARSE_RUNSHEET.out.samples
            STAGE_READS( ch_outdir, samples, Channel.value("raw") )
            raw_reads = STAGE_READS.out.reads
            samples_txt = STAGE_READS.out.samples_txt
        } else if ( ep == 'trimmed_reads' ) {
            samples = PARSE_RUNSHEET.out.samples
            STAGE_READS( ch_outdir, samples, Channel.value("trimmed") )
            trimmed_reads = STAGE_READS.out.reads
            samples_txt = STAGE_READS.out.samples_txt
        } else if ( ep == 'bam_files' ) {
            samples = PARSE_RUNSHEET.out.samples
            COPY_BAMS(ch_outdir, samples)
            bam_files = COPY_BAMS.out.bam_files
            samples_txt = samples_txt_from(bam_files)
        } else if ( ep == 'genes_results' ) {
            samples = PARSE_RUNSHEET.out.samples
            COPY_GENES_RESULTS(ch_outdir, samples)
            genes_results = COPY_GENES_RESULTS.out.genes_results
            samples_txt = samples_txt_from(genes_results)
        } else if ( ep == 'counts_table' ) {
            samples = PARSE_RUNSHEET.out.samples.map { meta, files -> meta }
            COPY_COUNTS_TABLE(ch_outdir, PARSE_RUNSHEET.out.table)
            counts_table = COPY_COUNTS_TABLE.out.counts_table
        } else if ( ep == 'dge_table' ) {
            samples = PARSE_RUNSHEET.out.samples.map { meta, files -> meta }
            COPY_DGE_TABLE(ch_outdir, PARSE_RUNSHEET.out.table)
            dge_table = COPY_DGE_TABLE.out.dge_table
        }

    emit:
        ch_outdir       = ch_outdir
        samples         = samples
        samples_txt     = samples_txt
        runsheet_path   = runsheet_path
        isa_archive     = isa_archive
        osd_accession   = osd_accession
        glds_accession  = glds_accession
        raw_reads       = raw_reads
        trimmed_reads   = trimmed_reads
        bam_files       = bam_files
        genes_results   = genes_results
        counts_table    = counts_table
        dge_table       = dge_table
}
