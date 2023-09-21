#!/usr/bin/env nextflow

/************************************************************************
===========================
singlecell-catch-nf
===========================

Authors
- Sergej Nowoshilow (sergej.nowoshilow@boehringer-ingelheim.com)
************************************************************************/

params.libraries = ""
params.chunkSize = 1_000_000
params.maxDist = 2
params.minReads = 10
params.majorityVote = 90
params.refName = "Day0"
//params.crindex = "/data/gcbds/vie/bioinf/nowoshilow/Data/indices/CellRanger/GRCh38_no_alt/"
params.crindex = "/data/gcbds/users/nowoshil/Projects/CaTCH2.0/KRASi_ON_vs_OFF/cellranger_reference/GRCh38_no_alt_CaTCH/"
params.max_mt_percent = 10
params.min_detected_features = 500
params.hvg_cutoff = 0.1

params.stopOnWarnings = true


def helpMessage() {
    log.info"""
    ======================================================================
      singlecell-catch-nf

      The pipeline performs an analysis of the PCR amplified CaTCH library.


      Version: ${workflow.manifest.version}
      Contact: Sergej Nowoshilow (sergej.nowoshilow@boehringer-ingelheim.com)
    ======================================================================

      Usage:
      nextflow run PIPELINES/singlecell-catch-nf/main.nf --libraries <list of libraries and FASTQ files> 
                                                         -profile ${params.profileList}

      Mandatory arguments:
        --libraries             CSV file with the following columns: 
                                    SampleName      name of the sample (can appear in multiple lines, in case
                                                    the library was sequenced in several runs)
                                    Treatment       sample treatment
                                    Replicate       replicate (even if a single replicate is present, this column
                                                    cannot be missing or be empty)
                                    LibraryType     either GEX or scCaTCH
                                    R1              path to the R1 read
                                    R2              path to the corresponding R2 read  

        --chunkSize             number of reads per chunk (default: ${params.chunkSize})

      Optional arguments:

        --help                  print this help message

        --outputDir             specifies the output directory. Default: ${params.outputDir}
        --reportsDir            specifies the reports directory. Default: ${params.reportsDir}

      ${params.profileDescription}
    """.stripIndent()
}


// Show the help message if --help is specified or if essential arguments are not provided 
if (params.help) {
    helpMessage()
    exit 0
}

if (!params.libraries) {
    log.info("ERROR: --libraries is not specified")
    helpMessage()
    exit 1
}


// Check all required input files before they are fed to a channel because this 
// causes the pipeline to fail immediately on AWS before starting up the machines
if(!file(params.libraries).exists()) {
    log.info("ERROR: the file '${params.libraries}' does not exist")
    exit 2
}


// Print execution parameters to stdout
log.info """
 ======================================================================
 | singlecell-catch-nf
 |
 | Version: ${workflow.manifest.version}
 ----------------------------------------------------------------------
 |
 | Mandatory arguments
 |   libraries           : ${params.libraries}
 |   chunk size          : ${params.chunkSize}
 |
 | Optional arguments
 |   outputDir           : ${params.outputDir}
 |   reportsDir          : ${params.reportsDir}
 |
 ======================================================================
""".stripIndent()


// ---------------------------------------------------------------------
// Pipeline Channels and Processes
// ---------------------------------------------------------------------

// For more information about syntax, please refer to the nextflow documentation at https://www.nextflow.io/docs/latest/index.html
stopOnWarn = (params.stopOnWarnings) ? "yes" : "no"


Channel
    .fromPath(params.libraries)
    .splitCsv(sep: "\t", header: true)
    .into { Ch_csv_GEX; Ch_csv_scCaTCH; Ch_csv_preprocess }


/************************************************************************
                STEP 1: Run CellRanger count
************************************************************************/
Ch_csv_GEX
    .filter { it.LibraryType == "GEX" }
    .branch {
        raw: (new File(it.R1)).isFile()
        precomputed: (new File(it.R1)).isDirectory()
    }
    .set { Ch_csv_GEX_split }

Ch_csv_GEX_split.raw
    .map { row -> tuple(row.SampleName, file(row.R1), file(row.R2)) }
    .groupTuple(by: 0)
    .combine( Channel.fromPath(params.crindex) )
    .set { Ch_cellranger_input }


process runCellrangerCount {

    cache "lenient"

    tag "${sampleName}"

    publishDir "outputs/cellranger/", mode: "copy"

    input:
        tuple val(sampleName), file("inputs/R1_*"), file("inputs/R2_*"), path(index) from Ch_cellranger_input

    output:
        path "${sampleName}"
        tuple val(sampleName), file("${sampleName}/analysis/tsne/gene_expression_2_components/projection.csv") into (Ch_cell_ids_from_raw, 
                                                                                                                     Ch_cell_ids_analytics_from_raw)
        tuple val(sampleName), path("${sampleName}/filtered_feature_bc_matrix") into Ch_cell_data_from_raw

    script:
    """
    # Find all reads, sort them by name to ensure that the paired files are on the consecutive lines,
    # and then create symlinks with proper names (SampleName_S1_R1_xxx.fastq.gz)
    IDX=1
    BKP=\${IFS}
    IFS=\$'\\n'
    for LINE in \$(find inputs/ -name "R[12]_*" -exec readlink -f {} \\; | sort | paste - -);
    do
        SUFFIX=\$(printf "%03d" \${IDX})

        R1=\$(echo \${LINE} | cut -f1)
        NEW_NAME=${sampleName}_S1_R1_\${SUFFIX}.fastq.gz
        ln -s \${R1} inputs/\${NEW_NAME}

        R2=\$(echo \${LINE} | cut -f2)
        NEW_NAME=${sampleName}_S1_R2_\${SUFFIX}.fastq.gz
        ln -s \${R2} inputs/\${NEW_NAME}

        IDX=\$((IDX + 1))
    done

    IFS=\${BKP}

    cellranger count \
        --disable-ui \
        --jobmode local \
        --localcores ${task.cpus} \
        --localmem 96 \
        --transcriptome ${index} \
        --id ${sampleName} \
        --fastqs inputs

    mv ${sampleName} rundir
    mv rundir/outs ${sampleName}
    """
}


Ch_csv_GEX_split.precomputed
    .map { row -> tuple(row.SampleName, file(row.R1)) }
    .set { Ch_cellranger_precomputed }

process useCellrangerData {

    cache "lenient"

    tag "${sampleName}"

    publishDir "outputs/cellranger/", mode: "copy"

    input:
        tuple val(sampleName), path("cr_data") from Ch_cellranger_precomputed

    output:
        path "${sampleName}"
        tuple val(sampleName), file("${sampleName}/analysis/tsne/gene_expression_2_components/projection.csv") into (Ch_cell_ids_from_precomputed, 
                                                                                                                     Ch_cell_ids_analytics_from_precomputed)
        tuple val(sampleName), path("${sampleName}/filtered_feature_bc_matrix") into Ch_cell_data_from_precomputed

    script:
        """
        mv cr_data ${sampleName}
        """
}

Ch_cell_ids_from_raw
    .mix(Ch_cell_ids_from_precomputed)
    .set { Ch_cell_ids }

Ch_cell_ids_analytics_from_raw
    .mix(Ch_cell_ids_analytics_from_precomputed)
    .set { Ch_cell_ids_analytics }

Ch_cell_data_from_raw
    .mix(Ch_cell_data_from_precomputed)
    .set { Ch_cell_data }



/************************************************************************
                STEP 2: Count CaTCH barcodes in chunks separately
************************************************************************/
Ch_csv_scCaTCH
    .filter { it.LibraryType == "scCaTCH" }
    .map { row -> tuple(row.SampleName, file(row.R1), file(row.R2)) }
    .splitFastq(by: params.chunkSize, file: true, compress: true, pe: true)
    .combine(Ch_cell_ids, by: 0)
    .set { Ch_count_input }


process countBarcodesInChunks {

    input:
        tuple val(sampleName), file(r1), file(r2), file(cellIDs) from Ch_count_input

    output:
        tuple val(sampleName), file('counts'), file('reads') into Ch_counts_chunks_out

    script:
    """
    countBarcodesInChunks.py \
        --r1 ${r1} \
        --r2 ${r2} \
        --cellIDs ${cellIDs} \
        --counts counts \
    | tee log \
    | grep -Po "Read [0-9,]+ single cell entries" \
    | cut -d" " -f2 > reads 
    """
}


/************************************************************************
                STEP 3: Merge the chunks data
************************************************************************/
Ch_counts_chunks_out
    .groupTuple(by: 0)
    .set { Ch_counts_chunks_merge }


process mergeBarcodesInChunks {

    publishDir "outputs/libraries/unfiltered", 
        pattern: "*.sclib", 
        saveAs: {filename -> (new File(filename).name)},
        mode: "copy"

    publishDir "outputs/libraries/unfiltered", 
        pattern: "*.stats", 
        saveAs: {filename -> (new File(filename).name)},
        mode: "copy"

    input:
        tuple val(sampleName), file("counts/file*"), file("reads/file*") from Ch_counts_chunks_merge

    output:
        tuple val(sampleName), file('*.sclib') into (Ch_merged_libraries, Ch_merged_analytics)
        tuple val(sampleName), file('*.stats') into Ch_unfiltered_stats

    script:
    """
    find counts -name "file*" > librarieslist
    find reads -name "file*" > readcountslist

    mergeChunkCounts.py \
        --libraries librarieslist \
        --readcounts readcountslist \
        --outlib ${sampleName}.sclib \
    | tee ${sampleName}.stats
    """
}


/************************************************************************
    STEP 4: Collapse similar barcodes and remove the background noise
************************************************************************/
process collapseAndFilterBarcodes {

    publishDir "outputs/libraries/collapsed", 
        pattern: "*.collapsed.sclib", 
        saveAs: {filename -> (new File(filename).name)},
        mode: "copy"

    publishDir "outputs/libraries/collapsed", 
        pattern: "*.collapsed.stats", 
        saveAs: {filename -> (new File(filename).name)},
        mode: "copy"

    input:
        tuple val(sampleName), file(library) from Ch_merged_libraries

    output:
        tuple val(sampleName), file('*.collapsed.sclib') into (Ch_collapsed_libraries, Ch_collapsed_analytics)
        tuple val(sampleName), file('*.collapsed.stats') into Ch_collapsed_stats

    script:
    """
    collapseCaTCHbarcodes.py \
        --library ${library} \
        --maxdist ${params.maxDist} \
        --minsupport ${params.minReads} \
        --outlib ${sampleName}.collapsed.sclib \
    | tee ${sampleName}.collapsed.stats
    """
}


/************************************************************************
                    STEP 5: Resolve multiplets
************************************************************************/
process resolveMultiplets {

    publishDir "outputs/libraries/resolved_multiplets", 
        pattern: "*.resolved_multiplets.sclib", 
        saveAs: {filename -> (new File(filename).name)},
        mode: "copy"

    publishDir "outputs/libraries/resolved_multiplets", 
        pattern: "*.resolved_multiplets.stats", 
        saveAs: {filename -> (new File(filename).name)},
        mode: "copy"

    input:
        tuple val(sampleName), file(library) from Ch_collapsed_libraries

    output:
        tuple val(sampleName), file("*.resolved_multiplets.sclib") into (Ch_resolved_multiplets_libraries, Ch_resolved_multiplets_analytics)
        tuple val(sampleName), file("*.resolved_multiplets.stats") into Ch_resolved_stats

    script:
    """
    resolveMultiplets.py \
        --library ${library} \
        --majority ${params.majorityVote} \
        --outlib ${sampleName}.resolved_multiplets.sclib \
    | tee ${sampleName}.resolved_multiplets.stats
    """
}


/************************************************************************
                    STEP 6: Generate reports
************************************************************************/
process generateReports {
    
    publishDir "outputs/reports", 
        pattern: "*.CaTCHbarcodes", 
        saveAs: {filename -> (new File(filename).name)},
        mode: "copy"

    publishDir "outputs/reports", 
        pattern: "*.cells", 
        saveAs: {filename -> (new File(filename).name)},
        mode: "copy"

    input:
        tuple val(sampleName), file(library) from Ch_resolved_multiplets_libraries

    output:
        file('*.CaTCHbarcodes') into Ch_report_CaTCHbarcodes
        tuple val(sampleName), file('*.cells') into Ch_report_cells

    script:
    """
    generateOutputTables.py \
        --library ${library} \
        --CaTCH ${sampleName}.CaTCHbarcodes \
        --cells ${sampleName}.cells
    """
}


/************************************************************************
                    STEP 7: Analytics report
************************************************************************/
Ch_cell_ids_analytics
    .join(Ch_merged_analytics, by: 0)
    .join(Ch_collapsed_analytics, by: 0)
    .join(Ch_resolved_multiplets_analytics, by: 0)
    .set { Ch_analytics_in }

process generateAnalyticsPlots {

    publishDir "outputs/analytics/plots", 
        pattern: "*.png", 
        saveAs: {filename -> (new File(filename).name)},
        mode: "copy"

    input:
        tuple val(sampleName), file(cell_ids), file(unfiltered), file(collapsed), file(resolved) from Ch_analytics_in

    output:
        tuple val(sampleName), file("*.png") into Ch_analytics_out

    script:
    """
    touch dummy.png
    """
}


/************************************************************************
                    STEP 8: Generate SingleCellExperiment object
************************************************************************/
Channel
    .fromPath("/home/nowoshil/Repositories/nf-pipelines/pipelines-singlecell-catch-nf/docker/scripts/R/preprocessData.R")
    .set { Ch_script }

Ch_cell_data
    .combine(Ch_report_cells, by: 0)
    .combine(Ch_script)
    .set { Ch_preprocess_input }

process preprocessSingleCellData {

    publishDir "outputs/sce/unfiltered", 
        pattern: "*.rda", 
        saveAs: {filename -> (new File(filename).name)},
        mode: "copy"

    input:
        tuple val(sampleName), path(featureMatrix), file(catchBarcodes), file(script) from Ch_preprocess_input

    output:
        file("*.sce.unfiltered.rda") into (Ch_basic_sce, Ch_de_barcodes_sce)

    script:
    """
    #Rscript --vanilla /tools/scripts/R/preprocessData.R 
    Rscript --vanilla ${script} \
       --sample ${sampleName} \
       --data10X ${featureMatrix} \
       --catchBC ${catchBarcodes} \
       --max_mt ${params.max_mt_percent} \
       --min_features ${params.min_detected_features} \
       --hvg_cutoff ${params.hvg_cutoff} \
       --out ${sampleName}.sce.unfiltered.rda
    """
}


/************************************************************************
                    STEP 9: Generate overview plots
************************************************************************
process createOverviewPlots {
    publishDir "outputs/plots", 
        pattern: "*.pdf", 
        saveAs: {filename -> (new File(filename).name)},
        mode: "copy"

    input:
        file(sce) from Ch_basic_sce

    output:
        file("overview.pdf")

    script:
    """
    Rscript --vanilla /tools/scripts/R/create_overview_plots.R \
        --sce ${sce} \
        --out overview.pdf \
        --format pdf \
        --width 25 \
        --height 10
    """
}

process createBarcodeEnrichmentPlots {
    publishDir "outputs/plots", 
        pattern: "*.jpeg", 
        saveAs: {filename -> (new File(filename).name)},
        mode: "copy"

    input:
        file(sce) from Ch_de_barcodes_sce

    output:
        file("*.jpeg")

    script:
    """
    Rscript --vanilla /tools/scripts/R/plot_enriched_and_depleted_BCs.R \
        --sce ${sce} \
        --plots_per_row 5 \
        --format jpeg \
        --width 400 \
        --height 300 \
        --outdir .
    """
}*/