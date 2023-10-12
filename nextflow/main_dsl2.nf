#!/usr/bin/env nextflow

/************************************************************************
===========================
singlecell-catch-nf
===========================

Authors
- Sergej Nowoshilow (sergej.nowoshilow@boehringer-ingelheim.com)
- Joerg Fallmann (joerg.fallmann@imp.ac.at)
************************************************************************/

//Version Check
nextflowVersion = '>=20.01.0.5264'
nextflow.enable.dsl=2

//define unset Params
def get_always(parameter){
    if (!params.containsKey(parameter)){
        params.put(parameter, null)
    }
    return params[parameter]
}

//Params from CL
libraries = get_always('libraries')
chunkSize = get_always('chunkSize') ?: 1_000_000
maxDist = get_always('maxDist') ?: 2
miReads = get_always('minReads') ?: 10
majorityVote = get_always('majorityVote') ?: 90
refName = get_always('refName') ?: "Day0"
crindex = get_always('crindex') ?: "${workflow.workDir}/../GENOMES/Human/INDICES/cellranger_t2t"
max_mt_percent = get_always('max_mt_percent') ?: 10
min_detected_features = get_always('min_detected_features') ?: 500
hvg_cutoff = get_always('hvg_cutoff') ?: 0.1
reportsDir = get_always('reportsDir') ?: "${workflow.workDir}/../REPORTS"
outputDir = get_always('outputDir') ?: "${workflow.workDir}/../scCaTCH_nf_OUTPUT"

stopOnWarnings = get_always('stopOnWarnings') ?: true


def helpMessage() {
    log.info"""
    ======================================================================
      singlecell-catch-nf

      The pipeline performs an analysis of the PCR amplified CaTCH library.


      Version: ${workflow.manifest.version}
      Contact: Sergej Nowoshilow (sergej.nowoshilow@boehringer-ingelheim.com), Joerg Fallmann (joerg.fallmann@imp.ac.at)
    ======================================================================

      Usage:
      nextflow run sccatch/nextflow/main_dsl2.nf --libraries <list of libraries and FASTQ files> 

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

      Optional arguments:
        --chunkSize             number of reads per chunk (default: ${params.chunkSize})
        --index                 Path to cellranger index file (default: ${params.crindex})
        --baseline              Name of reference day/condition (default: ${params.refName})
        --vote                  Number of votes needed for majority voting (default: ${params.majorityVote})
        --outputDir             specifies the output directory (default: ${params.outputDir})
        --reportsDir            specifies the reports directory.(default: ${params.reportsDir})
        --help                  print this help message

    """.stripIndent()
}


// Show the help message if --help is specified or if essential arguments are not provided 
if (params.help) {
    helpMessage()
    exit 0
}

if (!libraries) {
    log.info("ERROR: --libraries is not specified")
    helpMessage()
    exit 1
}


// Check all required input files before they are fed to a channel because this 
// causes the pipeline to fail immediately on AWS before starting up the machines
if(!file(libraries).exists()) {
    log.info("ERROR: the file '${libraries}' does not exist")
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
 |   libraries           : ${libraries}
 |   chunk size          : ${chunkSize}
 |
 | Optional arguments
 |   outputDir           : ${outputDir}
 |   reportsDir          : ${reportsDir}
 |
 ======================================================================
""".stripIndent()


// ---------------------------------------------------------------------
// Pipeline Channels and Processes
// ---------------------------------------------------------------------

// For more information about syntax, please refer to the nextflow documentation at https://www.nextflow.io/docs/latest/index.html
stopOnWarn = (stopOnWarnings) ? "yes" : "no"


/************************************************************************
                STEP 1: Run CellRanger count
************************************************************************/

process runCellrangerCount{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("filtered_feature_bc_matrix") > 0)       "OUTPUT/${sampleName}/CellRanger/filtered_feature_bc_matrix"
        else if (filename.indexOf("projection.csv") >0)          "OUTPUT/${sampleName}/CellRanger/tSNEs/gene_expression_2_components/projection.csv"
        else                                                     "OUTPUT/${sampleName}/CellRanger/${file(filename).getName()}"
    }

    //publishDir "outputs/cellranger/", mode: "copy"

    input:
        tuple val(sampleName), path("inputs/R1_*"), path("inputs/R2_*"), path(index)
    
    
    output:
        path "${sampleName}", emit: name
        tuple val(sampleName), path("${sampleName}/analysis/tsne/gene_expression_2_components/projection.csv"), emit: cell_ids_from_raw
        tuple val(sampleName), path("${sampleName}/filtered_feature_bc_matrix"), emit: cell_data_from_raw
    
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


process useCellrangerData{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${workflow.workDir}/../", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("filtered_feature_bc_matrix") >0)       "OUTPUT/${sampleName}/CellRanger/filtered_feature_bc_matrix"
        else if (filename.indexOf("projection.csv") >0)          "OUTPUT/${sampleName}/CellRanger/tSNEs/gene_expression_2_components/projection.csv"
        else                                                     "OUTPUT/${sampleName}/CellRanger/${file(filename).getName()}"
    }

    input:
        tuple val(sampleName), path("cr_data")

    output:
        path "${sampleName}", emit: name
        tuple val(sampleName), path("${sampleName}/analysis/tsne/gene_expression_2_components/projection.csv"), emit: cell_ids_from_precomputed
        tuple val(sampleName), path("${sampleName}/filtered_feature_bc_matrix"), emit: cell_data_from_precomputed

    script:
        """
        mv cr_data ${sampleName}
        """
}


/************************************************************************
                STEP 2: Count CaTCH barcodes in chunks separately
************************************************************************/

process countBarcodesInChunks{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${workflow.workDir}/../", mode: 'link',
    saveAs: {filename ->
        if (filename == "counts")       "OUTPUT/${sampleName}/CellRanger/counts/counts"
        else if (filename == "reads")          "OUTPUT/${sampleName}/CellRanger/counts/reads"
        else                                                     "OUTPUT/${sampleName}/CellRanger/counts/${file(filename).getName()}"
    }

    input:
        tuple val(sampleName), file(r1), file(r2), file(cellIDs)

    output:
        tuple val(sampleName), file('counts'), file('reads'), emit: counts_chunks_out

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

process mergeBarcodesInChunks{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${workflow.workDir}/../", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".sclib") > 0)       "OUTPUT/${sampleName}/CellRanger/libraries/unfiltered/${file(filename).getName()}"
        else if (filename.indexOf(".stats") > 0)             "OUTPUT/${sampleName}/CellRanger/libraries/unfiltered/${file(filename).getName()}"
        else                                                     "OUTPUT/${sampleName}/CellRanger/libraries/unfiltered/${file(filename).getName()}"
    }

    input:
        tuple val(sampleName), file("counts/file*"), file("reads/file*")

    output:
        tuple val(sampleName), file('*.sclib'), emit: merged_libraries
        tuple val(sampleName), file('*.stats'), emit: unfiltered_stats

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
process collapseAndFilterBarcodes{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${workflow.workDir}/../", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".collapsed.sclib") > 0)       "OUTPUT/${sampleName}/CellRanger/libraries/collapsed/${file(filename).getName()}"
        else if (filename.indexOf(".collapsed.stats") > 0)             "OUTPUT/${sampleName}/CellRanger/libraries/collapsed/${file(filename).getName()}"
        else                                                     "OUTPUT/${sampleName}/CellRanger/libraries/collapsed/${file(filename).getName()}"
    }
    
    input:
        tuple val(sampleName), file(library)

    output:
        tuple val(sampleName), file('*.collapsed.sclib'), emit: collapsed_libraries
        tuple val(sampleName), file('*.collapsed.stats'), emit: collapsed_stats

    script:
    """
    collapseCaTCHbarcodes.py \
        --library ${library} \
        --maxdist ${maxDist} \
        --minsupport ${minReads} \
        --outlib ${sampleName}.collapsed.sclib \
    | tee ${sampleName}.collapsed.stats
    """
}


/************************************************************************
                    STEP 5: Resolve multiplets
************************************************************************/
process resolveMultiplets{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${workflow.workDir}/../", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".resolved_multiplets.sclib") > 0)       "OUTPUT/${sampleName}/CellRanger/libraries/resolved_multiplets/${file(filename).getName()}"
        else if (filename.indexOf(".resolved_multiplets.stats") > 0)             "OUTPUT/${sampleName}/CellRanger/libraries/resolved_multiplets/${file(filename).getName()}"
        else                                                     "OUTPUT/${sampleName}/CellRanger/libraries/resolved_multiplets/${file(filename).getName()}"
    }
    
    input:
        tuple val(sampleName), file(library)

    output:
        tuple val(sampleName), file("*.resolved_multiplets.sclib"), emit: resolved_multiplets_libraries
        tuple val(sampleName), file("*.resolved_multiplets.stats"), emit: resolved_stats

    script:
    """
    resolveMultiplets.py \
        --library ${library} \
        --majority ${majorityVote} \
        --outlib ${sampleName}.resolved_multiplets.sclib \
    | tee ${sampleName}.resolved_multiplets.stats
    """
}


/************************************************************************
                    STEP 6: Generate reports
************************************************************************/

process generateReports{
    
    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${workflow.workDir}/../", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".CaTCHbarcodes") > 0)       "OUTPUT/${sampleName}/CellRanger/reports/${file(filename).getName()}"
        else if (filename.indexOf(".cells") > 0)             "OUTPUT/${sampleName}/CellRanger/reports/${file(filename).getName()}"
        else                                                     "OUTPUT/${sampleName}/CellRanger/libraries/resolved_multiplets/${file(filename).getName()}"
    }

    input:
        tuple val(sampleName), file(library)

    output:
        file('*.CaTCHbarcodes'), emit: report_CaTCHbarcodes
        tuple val(sampleName), file('*.cells'), emit: report_cells

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

process generateAnalyticsPlots{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${workflow.workDir}/../", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".png") > 0)       "OUTPUT/${sampleName}/CellRanger/analytics/plots/${file(filename).getName()}"
        else                                     "OUTPUT/${sampleName}/CellRanger/analytics/plots/${file(filename).getName()}"
    }

    input:
        tuple val(sampleName), file(cell_ids), file(unfiltered), file(collapsed), file(resolved)

    output:
        tuple val(sampleName), file("*.png"), emit: analytics_out

    script:
    """
    touch dummy.png
    """
}


/************************************************************************
                    STEP 8: Generate SingleCellExperiment object
************************************************************************/

process preprocessSingleCellData{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${workflow.workDir}/../", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".rda") > 0)       "OUTPUT/${sampleName}/CellRanger/sce/unfiltered/${file(filename).getName()}"
        else                                     "OUTPUT/${sampleName}/CellRanger/sce/unfiltered/${file(filename).getName()}"
    }

    input:
        tuple val(sampleName), path(featureMatrix), file(catchBarcodes), file(script)

    output:
        file("*.sce.unfiltered.rda"), emit: basic_sce

    script:
    """
    #Rscript --vanilla /tools/scripts/R/preprocessData.R 
    Rscript --vanilla ${script} \
       --sample ${sampleName} \
       --data10X ${featureMatrix} \
       --catchBC ${catchBarcodes} \
       --max_mt ${max_mt_percent} \
       --min_features ${min_detected_features} \
       --hvg_cutoff ${hvg_cutoff} \
       --out ${sampleName}.sce.unfiltered.rda
    """
}


/************************************************************************
                    STEP 9: Generate overview plots
************************************************************************/
process createOverviewPlots{
    
    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${workflow.workDir}/../", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".pdf") > 0)       "OUTPUT/${sampleName}/CellRanger/plots/${file(filename).getName()}"
        else                                     "OUTPUT/${sampleName}/CellRanger/plots/${file(filename).getName()}"
    }

    input:
        file(sce)

    output:
        file("overview.pdf"), emit: pdf

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

process createBarcodeEnrichmentPlots{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${workflow.workDir}/../", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".jpeg") > 0)       "OUTPUT/${sampleName}/CellRanger/plots/${file(filename).getName()}"
        else                                     "OUTPUT/${sampleName}/CellRanger/plots/${file(filename).getName()}"
    }

    input:
        file(sce)

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
}

/***********************************************************************
                        MAIN WORKFLOW
************************************************************************/

workflow{
    main:

        /**********************************************************
                STEP 0: Prepare Input
        ***********************************************************/
        
        (Ch_csv_GEX, Ch_csv_scCaTCH, Ch_csv_preprocess) = Channel.fromPath(libraries).splitCsv(sep: "\t", header: true)

        Ch_csv_GEX_split = Ch_csv_GEX.filter { it.LibraryType == "GEX" }.branch{
                raw: (new File(it.R1)).isFile()
                precomputed: (new File(it.R1)).isDirectory()
            }
        

        /**********************************************************
                STEP 1: Run CellRanger count
        ***********************************************************/

        Ch_cellranger_input = Ch_csv_GEX_split.raw.map { row -> tuple(row.SampleName, file(row.R1), file(row.R2)) }.groupTuple(by: 0).combine( Channel.fromPath(crindex) ).set {  }

        runCellrangerCount(Ch_cellranger_input)

        Ch_cellranger_precomputed = Ch_csv_GEX_split.precomputed.map { row -> tuple(row.SampleName, file(row.R1)) }
        
        useCellrangerData(Ch_cellranger_precomputed)

    
        /**************************************************************
                STEP 2: Count CaTCH barcodes in chunks separately
        ***************************************************************/
        
        Ch_count_input = Ch_csv_scCaTCH.filter { it.LibraryType == "scCaTCH" }.map { row -> tuple(row.SampleName, file(row.R1), file(row.R2)) }.splitFastq(by: chunkSize, file: true, compress: true, pe: true).combine(Ch_cell_ids, by: 0)

        countBarcodesInChunks(Ch_count_input)


        /**************************************************************
                STEP 3: Merge the chunks data
        ***************************************************************/

        Ch_counts_chunks_merge = countBarcodesInChunks.out.counts_chunks_out.groupTuple(by: 0)

        mergeBarcodesInChunks(Ch_counts_chunks_merge)


        /**************************************************************
                STEP 4: Collapse similar barcodes and \
                remove the background noise
        ***************************************************************/

        collapseAndFilterBarcodes(mergeBarcodesInChunks.out.merged_libraries)


        /**************************************************************
                STEP 5: Resolve multiplets
        ***************************************************************/

        resolveMultiplets(collapseAndFilterBarcodes.out.collapsed_libraries)


        /**************************************************************
                STEP 6: Generate reports
        ***************************************************************/

        generateReports(resolveMultiplets.out.resolved_multiplets_libraries)


        /**************************************************************
                STEP 7: Analytics report
        ***************************************************************/
        
        Ch_cell_ids = runCellrangerCount.out.cell_ids_from_raw.mix(useCellrangerData.out.cell_ids_from_precomputed)

        Ch_cell_data = runCellrangerCount.out.cell_data_from_raw.mix(useCellrangerData.out.cell_data_from_precomputed)

        Ch_analytics_in =  Ch_cell_ids.join(mergeBarcodesInChunks.out.merged_libraries, by: 0).join(collapseAndFilterBarcodes.out.collapsed_libraries, by: 0).join(resolveMultiplets.out.resolved_multiplets_libraries, by: 0)

        generateAnalyticsPlots(Ch_analytics_in)


        /**************************************************************
                STEP 8: Generate SingleCellExperiment object
        ***************************************************************/

        //Ch_script = Channel.fromPath("/home/nowoshil/Repositories/nf-pipelines/pipelines-singlecell-catch-nf/docker/scripts/R/preprocessData.R").set { Ch_script }
        Ch_preprocess_input = Ch_cell_data.combine(generateReports.out.report_cells, by: 0).combine('/tools/scripts/R/preprocessData.R')

        preprocessSingleCellData(Ch_preprocess_input)


        /**************************************************************
                STEP 9: Generate overview plots
        ***************************************************************/

        createOverviewPlots(preprocessSingleCellData.out.basic_sce)
        createBarcodeEnrichmentPlots(preprocessSingleCellData.out.basic_sce)

        
    /*
    emit:
    createOverviewPlots.out.pdf
    createBarcodeEnrichmentPlots.out.jpeg
    */
}
