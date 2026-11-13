include { STAGE } from '../subworkflows/stage.nf'

include { PARSE_ANNOTATIONS_TABLE } from '../modules/parse_annotations_table.nf'
include { DOWNLOAD_REFERENCES } from '../modules/download_references.nf'
include { SUBSAMPLE_GENOME } from '../modules/subsample_genome.nf'
include { DOWNLOAD_ERCC } from '../modules/download_ercc.nf'
include { CONCAT_ERCC } from '../modules/concat_ercc.nf'
include { ERCC_ANALYSIS } from '../modules/ercc_analysis.nf'
include { GTF_TO_PRED } from '../modules/gtf_to_pred.nf'
include { PRED_TO_BED } from '../modules/pred_to_bed.nf'

include { FASTQC as RAW_FASTQC } from '../modules/fastqc.nf'
include { FASTQC as TRIMMED_FASTQC } from '../modules/fastqc.nf'
include { GET_MAX_READ_LENGTH } from '../modules/get_max_read_length.nf'
include { TRIMGALORE } from '../modules/trimgalore.nf'

include { BUILD_STAR_INDEX } from '../modules/build_star_index.nf'
include { ALIGN_STAR } from '../modules/align_star.nf'
include { BUILD_RSEM_INDEX } from '../modules/build_rsem_index.nf'
include { QUANTIFY_STAR_GENES } from '../modules/quantify_star_genes.nf'
include { COUNT_ALIGNED } from '../modules/count_aligned.nf'
include { QUANTIFY_RSEM_GENES } from '../modules/quantify_rsem_genes.nf'

include { BUILD_BOWTIE2_INDEX } from '../modules/build_bowtie2_index.nf'
include { ALIGN_BOWTIE2 } from '../modules/align_bowtie2.nf'
include { GET_GTF_FEATURES } from '../modules/get_gtf_features.nf'
include { FEATURECOUNTS } from '../modules/featurecounts.nf'
include { QUANTIFY_FEATURECOUNTS_GENES } from '../modules/quantify_featurecounts_genes.nf'

include { SORT_AND_INDEX_BAM } from '../modules/sort_and_index_bam.nf'

include { INFER_EXPERIMENT; GENEBODY_COVERAGE; INNER_DISTANCE; READ_DISTRIBUTION } from '../modules/rseqc.nf'
include { ASSESS_STRANDEDNESS } from '../modules/assess_strandedness.nf'

include { EXTRACT_RRNA } from '../modules/extract_rrna.nf'
include { REMOVE_RRNA } from '../modules/remove_rrna.nf'
include { REMOVE_RRNA_FEATURECOUNTS } from '../modules/remove_rrna_featurecounts.nf'
include { REMOVE_RRNA_COUNTS_TABLE } from '../modules/remove_rrna_counts_table.nf'
include { DGE_DESEQ2 } from '../modules/dge_deseq2.nf'
include { DGE_DESEQ2 as DGE_DESEQ2_RRNA_RM } from '../modules/dge_deseq2.nf'
include { ANNOTATE_DGE_TABLE } from '../modules/annotate_dge_table.nf'

include { MULTIQC as RAW_READS_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as TRIMMED_READS_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as ALIGN_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as COUNT_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as GENEBODY_COVERAGE_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as INFER_EXPERIMENT_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as INNER_DISTANCE_MULTIQC } from '../modules/multiqc.nf'
include { MULTIQC as READ_DISTRIBUTION_MULTIQC } from '../modules/multiqc.nf'

include { PARSE_QC_METRICS } from '../modules/parse_qc_metrics.nf'
include { VV_RAW_READS;
    VV_TRIMMED_READS;
    VV_STAR_ALIGNMENT;
    VV_BOWTIE2_ALIGNMENT;
    VV_RSEQC;
    VV_RSEM_COUNTS;
    VV_FEATURECOUNTS;
    VV_DGE_DESEQ2;
    VV_CONCAT_FILTER } from '../modules/vv.nf'
include { SOFTWARE_VERSIONS } from '../modules/software_versions.nf'
include { GENERATE_PROTOCOL } from '../modules/generate_protocol.nf'

def strandedness_set() {
    return params.strandedness in ['none', 'forward', 'reverse']
}

def convert_strandedness(p) {
    if ( p == 'forward' ) { return 'sense' }
    if ( p == 'reverse' ) { return 'antisense' }
    return 'unstranded'
}

def rewrite_figshare(url) {
    if ( url && url.toString().contains('figshare.com/ndownloader/files/') ) {
        def file_id = (url.toString() =~ /.*\/files\/([a-zA-Z0-9]+).*/)[0][1]
        return "https://api.figshare.com/v2/file/download/${file_id}"
    }
    return url
}

def stage_gene_annotations(src) {
    if ( !src ) {
        return []
    }
    return file(rewrite_figshare(src.toString()), checkIfExists: true)
}

workflow RNASEQ {
    take:
        ch_outdir
        dp_tools_plugin
        annotations_csv_url_string
        accession
        isa_archive_path
        runsheet_path
        api_url
        reference_store_path
        derived_store_path

    main:
        def ep = params.entry_point
        def microbes = params.mode == 'microbes'
        dge_script = "${projectDir}/bin/dge_deseq2.Rmd"
        ercc_notebook = "${projectDir}/bin/combined_ercc_analysis.ipynb"
        ch_multiqc_config = params.multiqc_config ? channel.fromPath( params.multiqc_config ) : channel.fromPath("NO_FILE")

        println "Entry point: '${ep}' mode: '${params.mode}'"
        if ( microbes && ep == 'genes_results' ) {
            error "entry_point 'genes_results' is eukaryotic-only (RSEM). Use counts_table for microbes."
        }
        if ( ep == 'bam_files' && !strandedness_set() ) {
            error "--strandedness is required for entry_point bam_files (none, forward, or reverse)"
        }

        STAGE(
            ch_outdir,
            dp_tools_plugin,
            accession,
            isa_archive_path,
            runsheet_path,
            api_url
        )
        ch_outdir = STAGE.out.ch_outdir
        samples = STAGE.out.samples
        samples_txt = STAGE.out.samples_txt
        runsheet_path = STAGE.out.runsheet_path
        isa_archive = STAGE.out.isa_archive
        osd_accession = STAGE.out.osd_accession

        if ( ep in ['counts_table', 'dge_table'] ) {
            samples | first | set { ch_meta }
        } else {
            samples | first | map { meta, _files -> meta } | set { ch_meta }
        }
        ch_meta | map { meta -> meta.organism_sci } | set { organism_sci }

        PARSE_ANNOTATIONS_TABLE( annotations_csv_url_string, organism_sci )
        gene_annotations = params.gene_annotations_file
            ? channel.value(stage_gene_annotations(params.gene_annotations_file))
            : PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url.map { url -> stage_gene_annotations(url) }

        if ( params.reference_fasta && params.reference_gtf ) {
            channel.value( params.reference_source ) | set { reference_source }
            channel.value( params.reference_version ) | set { reference_version }
            channel.value( params.reference_fasta ) | set { reference_fasta_url }
            channel.value( params.reference_gtf ) | set { reference_gtf_url }
        } else {
            reference_source = PARSE_ANNOTATIONS_TABLE.out.reference_source
            reference_version = PARSE_ANNOTATIONS_TABLE.out.reference_version
            reference_fasta_url = PARSE_ANNOTATIONS_TABLE.out.reference_fasta_url
            reference_gtf_url = PARSE_ANNOTATIONS_TABLE.out.reference_gtf_url
        }

        channel.empty() | set { genome_references }
        channel.empty() | set { genome_references_pre_ercc }
        channel.empty() | set { genome_bed }
        channel.empty() | set { ch_versions }

        if ( ep != 'dge_table' ) {
            DOWNLOAD_REFERENCES( reference_store_path, organism_sci, reference_source, reference_version, reference_fasta_url, reference_gtf_url )
            genome_references_pre_subsample = DOWNLOAD_REFERENCES.out.reference_files

            if ( params.genome_subsample ) {
                SUBSAMPLE_GENOME( derived_store_path, organism_sci, genome_references_pre_subsample, reference_source, reference_version )
                SUBSAMPLE_GENOME.out.build | flatten | toList | set { genome_references_pre_ercc }
            } else {
                genome_references_pre_subsample | flatten | toList | set { genome_references_pre_ercc }
            }

            DOWNLOAD_ERCC( ch_meta.map { meta -> meta.has_ercc }, reference_store_path ).ifEmpty([file("ERCC92.fa"), file("ERCC92.gtf")]) | set { ch_maybe_ercc_refs }
            CONCAT_ERCC( reference_store_path, organism_sci, reference_source, reference_version, genome_references_pre_ercc, ch_maybe_ercc_refs, ch_meta.map { meta -> meta.has_ercc } )
                .ifEmpty { genome_references_pre_ercc.value } | set { genome_references }

            if ( ep in ['raw_reads', 'trimmed_reads', 'bam_files'] ) {
                GTF_TO_PRED(
                    derived_store_path,
                    organism_sci,
                    reference_source,
                    reference_version,
                    genome_references | map { refs -> refs[1] }
                )
                PRED_TO_BED(
                    derived_store_path,
                    organism_sci,
                    reference_source,
                    reference_version,
                    GTF_TO_PRED.out.genome_pred
                )
                genome_bed = PRED_TO_BED.out.genome_bed
                ch_versions = ch_versions.mix(GTF_TO_PRED.out.versions).mix(PRED_TO_BED.out.versions)
            }
        }

        channel.empty() | set { raw_mqc_data }
        channel.empty() | set { raw_mqc_zip }
        channel.empty() | set { trimmed_mqc_data }
        channel.empty() | set { trimmed_mqc_zip }
        channel.empty() | set { sorted_bam }
        channel.empty() | set { bam_to_transcriptome }
        channel.empty() | set { reads_per_gene }
        channel.empty() | set { bam_only_files }
        channel.empty() | set { align_mqc_data }
        channel.empty() | set { align_mqc_zip }
        channel.empty() | set { genes_results }
        channel.empty() | set { rsem_publishables }
        channel.empty() | set { ercc_counts }
        channel.empty() | set { star_publishables }
        channel.empty() | set { count_mqc_data }
        channel.empty() | set { count_mqc_zip }
        channel.empty() | set { counts }
        channel.empty() | set { infer_mqc_data }
        channel.empty() | set { infer_mqc_zip }
        channel.empty() | set { genebody_mqc_data }
        channel.empty() | set { genebody_mqc_zip }
        channel.empty() | set { inner_mqc_data }
        channel.empty() | set { inner_mqc_zip }
        channel.empty() | set { readdist_mqc_data }
        channel.empty() | set { readdist_mqc_zip }
        channel.empty() | set { dge_table }
        channel.empty() | set { dge_table_rrnarm }
        channel.empty() | set { counts_rrnarm }
        channel.empty() | set { norm_counts }
        channel.empty() | set { qc_counts }

        nf_version = '"NEXTFLOW":\n    nextflow: '.concat("${nextflow.version}\n")
        ch_versions = channel.value(nf_version).mix(ch_versions)

        if ( ep in ['raw_reads', 'trimmed_reads'] ) {
            if ( ep == 'raw_reads' ) {
                RAW_FASTQC( ch_outdir.map { dir -> dir + "/00-RawData/FastQC_Reports" }, STAGE.out.raw_reads )
                RAW_FASTQC.out.fastqc | map { qc -> [ qc[1], qc[2] ] }
                    | flatten
                    | collect
                    | set { raw_fastqc_zip }

                GET_MAX_READ_LENGTH( raw_fastqc_zip )
                max_read_length = GET_MAX_READ_LENGTH.out.length | map { n -> n.toString().toInteger() }

                TRIMGALORE( ch_outdir.map { dir -> dir + "/01-TG_Preproc" }, STAGE.out.raw_reads )
                trimmed_reads = TRIMGALORE.out.reads

                TRIMMED_FASTQC( ch_outdir.map { dir -> dir + "/01-TG_Preproc/FastQC_Reports" }, trimmed_reads )
                TRIMMED_FASTQC.out.fastqc | map { qc -> [ qc[1], qc[2] ] }
                    | flatten
                    | collect
                    | set { trimmed_fastqc_zip }

                RAW_READS_MULTIQC( ch_outdir.map { dir -> dir + "/00-RawData/MultiQC_Reports" }, samples_txt, raw_fastqc_zip, ch_multiqc_config, "raw_")
                TRIMMED_READS_MULTIQC( ch_outdir.map { dir -> dir + "/01-TG_Preproc/MultiQC_Reports" }, samples_txt, trimmed_fastqc_zip | concat( TRIMGALORE.out.reports ) | collect, ch_multiqc_config, "trimmed_")

                raw_mqc_data = RAW_READS_MULTIQC.out.data
                raw_mqc_zip = RAW_READS_MULTIQC.out.zipped_data
                ch_versions = ch_versions.mix(RAW_FASTQC.out.versions).mix(TRIMGALORE.out.versions).mix(RAW_READS_MULTIQC.out.versions)
            } else {
                trimmed_reads = STAGE.out.trimmed_reads

                TRIMMED_FASTQC( ch_outdir.map { dir -> dir + "/01-TG_Preproc/FastQC_Reports" }, trimmed_reads )
                TRIMMED_FASTQC.out.fastqc | map { qc -> [ qc[1], qc[2] ] }
                    | flatten
                    | collect
                    | set { trimmed_fastqc_zip }

                GET_MAX_READ_LENGTH( trimmed_fastqc_zip )
                max_read_length = GET_MAX_READ_LENGTH.out.length | map { n -> n.toString().toInteger() }

                TRIMMED_READS_MULTIQC( ch_outdir.map { dir -> dir + "/01-TG_Preproc/MultiQC_Reports" }, samples_txt, trimmed_fastqc_zip, ch_multiqc_config, "trimmed_")
                ch_versions = ch_versions.mix(TRIMMED_READS_MULTIQC.out.versions)
            }
            trimmed_mqc_data = TRIMMED_READS_MULTIQC.out.data
            trimmed_mqc_zip = TRIMMED_READS_MULTIQC.out.zipped_data

            if ( microbes ) {
                BUILD_BOWTIE2_INDEX( derived_store_path, organism_sci, reference_source, reference_version, genome_references, ch_meta )
                ALIGN_BOWTIE2( ch_outdir.map { dir -> dir + "/02-Bowtie2_Alignment" }, trimmed_reads, BUILD_BOWTIE2_INDEX.out.index_dir )
                SORT_AND_INDEX_BAM( ch_outdir.map { dir -> dir + "/02-Bowtie2_Alignment" }, ALIGN_BOWTIE2.out.bam )
                ALIGN_MULTIQC( ch_outdir.map { dir -> dir + "/02-Bowtie2_Alignment/MultiQC_Reports" }, samples_txt, ALIGN_BOWTIE2.out.alignment_logs | collect, ch_multiqc_config, "align_")

                sorted_bam = SORT_AND_INDEX_BAM.out.sorted_bam
                bam_only_files = SORT_AND_INDEX_BAM.out.sorted_bam.map { row -> row[1] } | toSortedList()
                align_mqc_data = ALIGN_MULTIQC.out.data
                align_mqc_zip = ALIGN_MULTIQC.out.zipped_data
                ch_versions = ch_versions.mix(ALIGN_BOWTIE2.out.versions)
            } else {
                BUILD_STAR_INDEX( derived_store_path, organism_sci, reference_source, reference_version, genome_references, ch_meta, max_read_length )
                ALIGN_STAR( ch_outdir.map { dir -> dir + "/02-STAR_Alignment" }, trimmed_reads, BUILD_STAR_INDEX.out.index_dir )
                SORT_AND_INDEX_BAM( ch_outdir.map { dir -> dir + "/02-STAR_Alignment" }, ALIGN_STAR.out.bam_by_coord )
                ALIGN_MULTIQC( ch_outdir.map { dir -> dir + "/02-STAR_Alignment/MultiQC_Reports" }, samples_txt, ALIGN_STAR.out.alignment_logs | collect, ch_multiqc_config, "align_")

                sorted_bam = SORT_AND_INDEX_BAM.out.sorted_bam
                bam_to_transcriptome = ALIGN_STAR.out.bam_to_transcriptome
                reads_per_gene = ALIGN_STAR.out.reads_per_gene
                bam_only_files = SORT_AND_INDEX_BAM.out.bam_only_files
                align_mqc_data = ALIGN_MULTIQC.out.data
                align_mqc_zip = ALIGN_MULTIQC.out.zipped_data
                ch_versions = ch_versions.mix(ALIGN_STAR.out.versions).mix(SORT_AND_INDEX_BAM.out.versions)
            }

            GENEBODY_COVERAGE( ch_outdir.map { dir -> dir + "/RSeQC_Analyses/02_geneBody_coverage" }, sorted_bam, genome_bed )
            INFER_EXPERIMENT( ch_outdir.map { dir -> dir + "/RSeQC_Analyses/03_infer_experiment" }, sorted_bam, genome_bed )
            INNER_DISTANCE( ch_outdir.map { dir -> dir + "/RSeQC_Analyses/04_inner_distance" }, sorted_bam, genome_bed, max_read_length )
            READ_DISTRIBUTION( ch_outdir.map { dir -> dir + "/RSeQC_Analyses/05_read_distribution" }, sorted_bam, genome_bed )

            infer_expt_out = INFER_EXPERIMENT.out.log | map { row -> row[1] } | collect
            ASSESS_STRANDEDNESS( infer_expt_out )
            if ( strandedness_set() ) {
                strandedness = channel.value(convert_strandedness(params.strandedness))
            } else {
                strandedness = ASSESS_STRANDEDNESS.out | map { f -> f.text.split(":")[0] }
            }

            INFER_EXPERIMENT_MULTIQC( ch_outdir.map { dir -> dir + "/RSeQC_Analyses/MultiQC_Reports" }, samples_txt, INFER_EXPERIMENT.out.log | map { row -> row[1] } | collect, ch_multiqc_config, "infer_exp_")
            GENEBODY_COVERAGE_MULTIQC( ch_outdir.map { dir -> dir + "/RSeQC_Analyses/MultiQC_Reports" }, samples_txt, GENEBODY_COVERAGE.out.log | map { row -> row[1] } | collect, ch_multiqc_config, "geneBody_cov_")
            INNER_DISTANCE_MULTIQC( ch_outdir.map { dir -> dir + "/RSeQC_Analyses/MultiQC_Reports" }, samples_txt, INNER_DISTANCE.out.log | map { row -> row[1] } | collect, ch_multiqc_config, "inner_dist_")
            READ_DISTRIBUTION_MULTIQC( ch_outdir.map { dir -> dir + "/RSeQC_Analyses/MultiQC_Reports" }, samples_txt, READ_DISTRIBUTION.out.log | map { row -> row[1] } | collect, ch_multiqc_config, "read_dist_")

            infer_mqc_data = INFER_EXPERIMENT_MULTIQC.out.data
            infer_mqc_zip = INFER_EXPERIMENT_MULTIQC.out.zipped_data
            genebody_mqc_data = GENEBODY_COVERAGE_MULTIQC.out.data
            genebody_mqc_zip = GENEBODY_COVERAGE_MULTIQC.out.zipped_data
            inner_mqc_data = INNER_DISTANCE_MULTIQC.out.data
            inner_mqc_zip = INNER_DISTANCE_MULTIQC.out.zipped_data
            readdist_mqc_data = READ_DISTRIBUTION_MULTIQC.out.data
            readdist_mqc_zip = READ_DISTRIBUTION_MULTIQC.out.zipped_data
            ch_versions = ch_versions.mix(INFER_EXPERIMENT.out.versions)
                .mix(GENEBODY_COVERAGE.out.versions)
                .mix(INNER_DISTANCE.out.versions)
                .mix(READ_DISTRIBUTION.out.versions)
        } else {
            strandedness = channel.value(convert_strandedness(params.strandedness))
        }

        if ( !microbes && ep in ['raw_reads', 'trimmed_reads', 'bam_files', 'genes_results'] ) {
            if ( ep == 'genes_results' ) {
                genes_results = STAGE.out.genes_results
                rsem_counts = genes_results | map { row -> row[1] } | collect
                QUANTIFY_RSEM_GENES( ch_outdir.map { dir -> dir + "/03-RSEM_Counts" }, samples_txt, rsem_counts )
            } else {
                if ( ep in ['raw_reads', 'trimmed_reads'] ) {
                    QUANTIFY_STAR_GENES( ch_outdir.map { dir -> dir + "/02-STAR_Alignment" }, samples_txt, reads_per_gene | toSortedList, strandedness )
                    star_publishables = QUANTIFY_STAR_GENES.out.publishables
                }

                def tx_bam = ep == 'bam_files' ? STAGE.out.bam_files : bam_to_transcriptome
                BUILD_RSEM_INDEX( derived_store_path, organism_sci, reference_source, reference_version, genome_references, ch_meta )
                COUNT_ALIGNED( ch_outdir.map { dir -> dir + "/03-RSEM_Counts" }, tx_bam, BUILD_RSEM_INDEX.out.index_dir, strandedness )
                genes_results = COUNT_ALIGNED.out.genes_results
                rsem_counts = COUNT_ALIGNED.out.counts | map { row -> row[1] } | collect
                QUANTIFY_RSEM_GENES( ch_outdir.map { dir -> dir + "/03-RSEM_Counts" }, samples_txt, rsem_counts )

                COUNT_MULTIQC( ch_outdir.map { dir -> dir + "/03-RSEM_Counts/MultiQC_Reports" }, samples_txt, rsem_counts, ch_multiqc_config, "RSEM_count_")
                count_mqc_data = COUNT_MULTIQC.out.data
                count_mqc_zip = COUNT_MULTIQC.out.zipped_data
                ch_versions = ch_versions.mix(COUNT_ALIGNED.out.versions).mix(COUNT_MULTIQC.out.versions)
            }
            rsem_publishables = QUANTIFY_RSEM_GENES.out.publishables
            qc_counts = QUANTIFY_RSEM_GENES.out.publishables
            ercc_counts = QUANTIFY_RSEM_GENES.out.publishables.map { unnorm, _nz -> unnorm }
        }

        if ( microbes && ep in ['raw_reads', 'trimmed_reads', 'bam_files'] ) {
            def bam_list = ep == 'bam_files' ? STAGE.out.bam_files.map { row -> row[1] }.collect() : bam_only_files
            GET_GTF_FEATURES( genome_references )
            gtf_features = GET_GTF_FEATURES.out.gtf_features.map { f -> f.text.trim() }
            FEATURECOUNTS( ch_outdir.map { dir -> dir + "/03-FeatureCounts" }, ch_meta, genome_references, gtf_features, strandedness, bam_list )
            QUANTIFY_FEATURECOUNTS_GENES( ch_outdir.map { dir -> dir + "/03-FeatureCounts" }, samples_txt, FEATURECOUNTS.out.counts )
            COUNT_MULTIQC( ch_outdir.map { dir -> dir + "/03-FeatureCounts/MultiQC_Reports" }, samples_txt, FEATURECOUNTS.out.summary | collect, ch_multiqc_config, "FeatureCounts_")

            counts = FEATURECOUNTS.out.counts
            count_mqc_data = COUNT_MULTIQC.out.data
            count_mqc_zip = COUNT_MULTIQC.out.zipped_data
            ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions).mix(COUNT_MULTIQC.out.versions)
        }

        if ( !microbes && ep == 'counts_table' ) {
            ercc_counts = STAGE.out.counts_table
        }

        if ( !microbes && ep != 'dge_table' ) {
            ERCC_ANALYSIS(
                ch_outdir.map { dir -> dir + "/ERCC_Analysis" },
                ch_meta.map { meta -> meta.has_ercc },
                STAGE.out.glds_accession.ifEmpty( params.accession ?: '' ),
                isa_archive.ifEmpty { [] },
                channel.value( params.assay_suffix ),
                ercc_counts,
                ch_meta.map { meta -> meta.organism_sci },
                ercc_notebook
            )
        }

        if ( ep == 'dge_table' ) {
            ANNOTATE_DGE_TABLE( ch_outdir, gene_annotations, ch_meta, STAGE.out.dge_table )
            ch_versions = ch_versions.mix(ANNOTATE_DGE_TABLE.out.versions)
        } else {
            EXTRACT_RRNA( organism_sci, genome_references | map { refs -> refs[1] } )

            if ( ep == 'counts_table' ) {
                def counts_dir = microbes ? "/03-FeatureCounts" : "/03-RSEM_Counts"
                DGE_DESEQ2( ch_outdir, ch_meta, gene_annotations, runsheet_path, STAGE.out.counts_table, dge_script, "" )
                REMOVE_RRNA_COUNTS_TABLE( ch_outdir.map { dir -> dir + counts_dir }, STAGE.out.counts_table, EXTRACT_RRNA.out.rrna_ids )
                counts_rrnarm = REMOVE_RRNA_COUNTS_TABLE.out.counts_rrnarm
                DGE_DESEQ2_RRNA_RM( ch_outdir, ch_meta, gene_annotations, runsheet_path, counts_rrnarm, dge_script, "_rRNArm" )
            } else if ( microbes ) {
                DGE_DESEQ2( ch_outdir, ch_meta, gene_annotations, runsheet_path, counts, dge_script, "" )
                REMOVE_RRNA_FEATURECOUNTS( ch_outdir.map { dir -> dir + "/03-FeatureCounts" }, counts, EXTRACT_RRNA.out.rrna_ids )
                counts_rrnarm = REMOVE_RRNA_FEATURECOUNTS.out.counts_rrnarm
                DGE_DESEQ2_RRNA_RM( ch_outdir, ch_meta, gene_annotations, runsheet_path, counts_rrnarm, dge_script, "_rRNArm" )
            } else {
                DGE_DESEQ2( ch_outdir, ch_meta, gene_annotations, runsheet_path, genes_results.map { row -> row[1] } | collect, dge_script, "" )
                REMOVE_RRNA( ch_outdir.map { dir -> dir + "/03-RSEM_Counts" }, EXTRACT_RRNA.out.rrna_ids, genes_results )
                counts_rrnarm = REMOVE_RRNA.out.genes_results_rrnarm | toSortedList
                DGE_DESEQ2_RRNA_RM( ch_outdir, ch_meta, gene_annotations, runsheet_path, counts_rrnarm, dge_script, "_rRNArm" )
            }

            dge_table = DGE_DESEQ2.out.dge_table
            dge_table_rrnarm = DGE_DESEQ2_RRNA_RM.out.dge_table
            norm_counts = DGE_DESEQ2.out.norm_counts
            ch_versions = ch_versions.mix(DGE_DESEQ2.out.versions)
            if ( microbes ) {
                qc_counts = DGE_DESEQ2.out.norm_counts | map { row -> row[1] }
            }
        }

        all_multiqc = raw_mqc_data
            | concat( trimmed_mqc_data )
            | concat( align_mqc_data )
            | concat( genebody_mqc_data )
            | concat( infer_mqc_data )
            | concat( inner_mqc_data )
            | concat( readdist_mqc_data )
            | concat( count_mqc_data )
            | collect

        if ( ep in ['raw_reads', 'trimmed_reads', 'bam_files', 'genes_results'] ) {
            PARSE_QC_METRICS(
                ch_outdir,
                osd_accession,
                ch_meta,
                isa_archive.ifEmpty { [] },
                all_multiqc,
                qc_counts,
                runsheet_path
            )
        }

        channel.empty() | set { vv_logs }

        if ( ep == 'raw_reads' ) {
            VV_RAW_READS( dp_tools_plugin, ch_outdir, ch_meta, runsheet_path, raw_mqc_zip )
            vv_logs = vv_logs.mix(VV_RAW_READS.out.log)
            ch_versions = ch_versions.mix(VV_RAW_READS.out.versions)
        }
        if ( ep in ['raw_reads', 'trimmed_reads'] ) {
            VV_TRIMMED_READS( dp_tools_plugin, ch_outdir, ch_meta, runsheet_path, trimmed_mqc_zip )
            vv_logs = vv_logs.mix(VV_TRIMMED_READS.out.log)
            ch_versions = ch_versions.mix(VV_TRIMMED_READS.out.versions)

            VV_RSEQC(
                dp_tools_plugin,
                ch_outdir,
                ch_meta,
                runsheet_path,
                genebody_mqc_zip,
                infer_mqc_zip,
                channel.empty() | mix(inner_mqc_zip) | collect | ifEmpty({ file("PLACEHOLDER") }),
                readdist_mqc_zip
            )
            vv_logs = vv_logs.mix(VV_RSEQC.out.log)
        }
        if ( ep in ['raw_reads', 'trimmed_reads'] && !microbes ) {
            VV_STAR_ALIGNMENT(
                dp_tools_plugin,
                ch_outdir,
                ch_meta,
                runsheet_path,
                align_mqc_zip,
                star_publishables,
                bam_only_files | collect
            )
            vv_logs = vv_logs.mix(VV_STAR_ALIGNMENT.out.log)
        }
        if ( ep in ['raw_reads', 'trimmed_reads'] && microbes ) {
            VV_BOWTIE2_ALIGNMENT(
                dp_tools_plugin,
                ch_outdir,
                ch_meta,
                runsheet_path,
                align_mqc_zip,
                bam_only_files | collect
            )
            vv_logs = vv_logs.mix(VV_BOWTIE2_ALIGNMENT.out.log)
        }
        if ( ep in ['raw_reads', 'trimmed_reads', 'bam_files'] && !microbes ) {
            VV_RSEM_COUNTS(
                dp_tools_plugin,
                ch_outdir,
                ch_meta,
                runsheet_path,
                count_mqc_zip,
                counts_rrnarm,
                rsem_publishables
            )
            vv_logs = vv_logs.mix(VV_RSEM_COUNTS.out.log)
            ch_versions = ch_versions.mix(VV_RSEM_COUNTS.out.versions)
        }
        if ( ep in ['raw_reads', 'trimmed_reads', 'bam_files'] && microbes ) {
            VV_FEATURECOUNTS(
                dp_tools_plugin,
                ch_outdir,
                ch_meta,
                runsheet_path,
                count_mqc_zip,
                counts_rrnarm
            )
            vv_logs = vv_logs.mix(VV_FEATURECOUNTS.out.log)
            ch_versions = ch_versions.mix(VV_FEATURECOUNTS.out.versions)
        }
        if ( ep != 'dge_table' ) {
            VV_DGE_DESEQ2(
                dp_tools_plugin,
                ch_outdir,
                ch_meta,
                runsheet_path,
                dge_table,
                dge_table_rrnarm
            )
            vv_logs = vv_logs.mix(VV_DGE_DESEQ2.out.log)
            ch_versions = ch_versions.mix(VV_DGE_DESEQ2.out.versions)

            VV_CONCAT_FILTER( ch_outdir, vv_logs | collect )
        }

        ch_versions
            | unique
            | collectFile( newLine: true )
            | set { ch_final_software_versions }
        SOFTWARE_VERSIONS( ch_outdir, ch_final_software_versions )

        if ( ep != 'dge_table' ) {
            GENERATE_PROTOCOL(
                ch_outdir,
                ch_meta,
                strandedness,
                SOFTWARE_VERSIONS.out.software_versions_yaml,
                reference_source,
                reference_version,
                genome_references_pre_ercc,
                runsheet_path
            )
        }
}
