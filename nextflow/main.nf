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
nextflowVersion = '>=23.04.2'
nextflow.enable.dsl=2

//define unset Params
def get_always(parameter){
    if (!params.containsKey(parameter)){
        params.put(parameter, null)
    }
    return params[parameter]
}

//parse chemistry params into defined paramstring
def chemistry_params_to_str(xtra){
    xtra = xtra.trim().replaceAll(" \n", "")
    tmpmap = xtra.tokenize("--").collectEntries{ 
               it.split(" ",2).with{ 
                   [ (it[0].replaceAll("--", "")): (it.size()<2) ? null : it[1] ?: null ] 
                }
            }
    //print(tmpmap)
    map = ["bcStart" : tmpmap["soloCBstart"].toInteger()-1, "bcLength" : tmpmap["soloCBlen"].toInteger(), "umiStart" : tmpmap["soloUMIstart"].toInteger()-1, "umiLength" : tmpmap["soloUMIlen"].toInteger()]
    xtra = map.collect{ k, v -> v ? "--" + k + " " + v : "" }.join(" ")
    return xtra
}

//Params from CL
absDir = workflow.launchDir
scriptDirR = get_always('scriptDirR') ?: params.scriptDirR
scriptDirPy = get_always('scriptDirPy') ?: params.scriptDirPy
binDir = get_always('binDir') ?: params.binDir
libraries = get_always('libraries')
chunkSize = get_always('chunkSize') ?: params.chunkSize
maxDist = get_always('maxDist') ?: params.maxDist
maxDistCaTCH = get_always('maxDistCaTCH') ?: params.maxDistCaTCH
maxDistUMIs = get_always('maxDistUMIs') ?: params.maxDistUMIs
clusterMethodCaTCH = get_always('clusterMethodCaTCH') ?: params.clusterMethodCaTCH
clusterMethodUMIs = get_always('clusterMethodUMIs') ?: params.clusterMethodUMIs
minReads = get_always('minReads') ?: params.minReads
majorityVote = get_always('majorityVote') ?: params.majorityVote
qcparams = get_always('fastqc_params') ?: params.qcParams
mapperbin = get_always('mapper') ?: params.mapperBin
runqc = get_always('withQC') ?: params.runQC
crparams = get_always('cellranger_params') ?: params.crParams
starparams = get_always('star_params') ?: params.starMappingParams
idxparams = get_always('idx_params') ?: params.idxParams
whitelist = get_always('whitelist') ?: params.whitelist
refName = get_always('refName') ?: params.refName
mapindex = get_always('index') ?: params.mapIndex
mapref = get_always('reference') ?: params.mapRef
mapanno = get_always('annotation') ?: params.mapAnno
filtering = get_always('filter') ?: params.filter
organism = get_always('organism') ?: params.organism
uniqueCaTCH = get_always('uniqueCaTCH') ?: params.uniqueCaTCH
minBC = get_always('min_detected_barcodes') ?: params.minDetectedBarcodes
singletCutoff = get_always('singlet_cutoff') ?: params.singletCutoff
bc1Cutoff = get_always('bc1_cutoff') ?: params.bc1Cutoff
bc2Cutoff = get_always('bc2_cutoff') ?: params.bc2Cutoff
maxMtPercent = get_always('max_mt_percent') ?: params.maxMtPercent
minDetectedFeatures = get_always('min_detected_features') ?: params.minDetectedFeatures
hvgCutoff = get_always('hvg_cutoff') ?: params.hvgCutoff
pvalCutoff = get_always('pval_cutoff') ?: params.pvalCutoff
lfcCutoff = get_always('lfc_cutoff') ?: params.lfcCutoff
markerfile = get_always('markers') ?: params.markerFile
reportsDir = get_always('reportsDir') ?: params.reportsDir
outputDir = get_always('outputDir') ?: params.outputDir
stopOnWarnings = get_always('stopOnWarnings') ?: true

// normalize string-like params to avoid Boolean/Null concatenation issues
crparams   = (crparams   instanceof Boolean || crparams   == null) ? '' : crparams.toString()
starparams = (starparams instanceof Boolean || starparams == null) ? '' : starparams.toString()
idxparams  = (idxparams  instanceof Boolean || idxparams  == null) ? '' : idxparams.toString()
qcparams   = (qcparams   instanceof Boolean || qcparams   == null) ? '' : qcparams.toString()


def helpMessage() {
    log.info"""
    ======================================================================
      singlecell-catch-nf

      The pipeline performs an analysis of the PCR amplified CaTCH library.


      Version: ${workflow.manifest.version}
      Contact: Sergej Nowoshilow (sergej.nowoshilow@boehringer-ingelheim.com), Joerg Fallmann (joerg.fallmann@imp.ac.at)
    ======================================================================

      Usage:
      nextflow run CaTCHseq/nextflow/main_dsl2.nf --libraries <list of libraries and FASTQ files> 

      Mandatory arguments:
        --libraries             CSV file with the following columns: 
                                    SampleName      name of the sample (can appear in multiple lines, in case
                                                    the library was sequenced in several runs)
                                    Condition       condition for this sample (e.g. timepoint, KO, treatment)
                                    Replicate       replicate (even if a single replicate is present, this column cannot be missing or be empty)
                                    LibraryType     either GEX or CaTCHseq
                                    R1              path to the R1 read
                                    R2              path to the R2 read (if available)
                                    CellNumber      number of expected cells (NA if not available, number for soft constraint, number! for hard constraint)
                                    Chemistry       chemistry used (e.g. 10X, DropIn, SmartSeq), this may not set adequate parameters for mappers automatically, make sure you check them accordingly

            Optional arguments:
                --outputDir             specifies the output directory (default: ${outputDir})
                --reportsDir            specifies the reports directory.(default: ${reportsDir})
                --scriptDirR            specifies the path to the R scripts directory, do not change if running with docker, otherwise set to path on CaTCHseq git repo.(default: /tools/scripts/R/ which is valid for docker instance; set to \${absDir}/CaTCHseq/docker/scripts/R/ for instances not running docker)
                --scriptDirPy           specifies the path to the Python scripts directory, do not change if running with docker, otherwise set to path on CaTCHseq git repo.(default: /tools/scripts/python/ which is valid for docker instance; set to \${absDir}/CaTCHseq/docker/scripts/python/ for instances not running docker)
                --binDir                specifies the path to the binary directory, do not change if running with docker, otherwise set to path on CaTCHseq git repo.(default: /usr/bin/local which is valid for docker instance; set to '' for instances not running docker)
                --mapper                Which mapper to run (default: CellRanger, optional: STAR)
                --index                 Path to mapper index directory (default: ${mapindex}, NEEDS TO BE SET ALSO TO CREATE NEW INDEX, new index will be stored at given path)
                --reference             Path to reference fasta.gz
                --annotation            Path to annotation gtf.gz 
                --whitelist             Path to barcode whitelist
                --withQC                Boolean, run FastQC and MultiQC (default: ${runqc})
                --fastqc_params         Optional parameters for FASTQC
                --cellranger_params     Optional parameters for CellRanger (default: ${crparams})
                --star_params           Optional parameters for STAR mapping
                --idx_params            Optional parameters for STAR index generation
                --filter                Postprocess filtered counts (default: ${filter})
                --minReads              Minimum reads to keep cell (default: ${minReads})
                --chunkSize             number of reads per chunk (default: ${chunkSize})
                --maxDist               maximum distance for barcode merging (default: ${maxDist})
                --maxDistCaTCH          maximum distance for CaTCH barcode collapsing (default: ${maxDistCaTCH})
                --maxDistUMIs           maximum distance for UMI collapsing (default: ${maxDistUMIs})
                --clusterMethodCaTCH    umi_tools cluster method for CaTCH barcode collapsing (default: ${clusterMethodCaTCH})
                --clusterMethodUMIs     umi_tools cluster method for UMI collapsing (default: ${clusterMethodUMIs})
                --majorityVote          Number of votes needed for majority voting (default: ${majorityVote})
                --uniqueCaTCH           Collapse CaTCH barcodes to unique or keep counts (default: ${uniqueCaTCH})
                --min_detected_barcodes Minimum number of CaTCH barcode reads per cell filter (default: ${minBC})
                --singlet_cutoff        Min ratio for sum of CaTCH barcodes 1 to classify cell as 'Singlet' (default: ${singletCutoff})
                --bc1_cutoff            Min ratio for CaTCH barcode 1 to classify cell as 'Dual_Integration' (default: ${bc1Cutoff})
                --bc2_cutoff            Min ratio for CaTCH barcode 2 to classify cell as 'Dual_Integration' (default: ${bc2Cutoff})
                --max_mt_percent        Maximum mitochondrial read percentage (default: ${maxMtPercent})
                --min_detected_features Minimum detected features cutoff (default: ${minDetectedFeatures})
                --hvg_cutoff            Highly variable genes cutoff (default: ${hvgCutoff})
                --pval_cutoff           P-value cutoff for DE analysis (default: ${pvalCutoff})
                --lfc_cutoff            Log-fold-change cutoff for DE analysis (default: ${lfcCutoff})
                --organism              Identifier for organism (choice: ["Human", "Mouse"], default: ${organism})
                --baseline              Name of reference day/condition (default: ${refName})
                --marker                RDS file of cellcycle markers (default: ${markerfile}, NamedList with gene names for each stage [G1S, S, G2M, M, MG1, G0], S and G2M are needed)
                --stopOnWarnings        Stop pipeline when warnings occur (default: ${stopOnWarnings})
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
 |   libraries               : ${libraries}
 |
 | Optional arguments
 |   outputDir               : ${outputDir}
 |   reportsDir              : ${reportsDir}
 |   scriptDirR              : ${scriptDirR}
 |   scriptDirPy             : ${scriptDirPy}
 |   binDir                  : ${binDir}
 |   mapper                  : ${mapperbin}
 |   index                   : ${mapindex}
 |   reference               : ${mapref}
 |   annotation              : ${mapanno}
 |   whitelist               : ${whitelist}
 |   withQC                  : ${runqc}
 |   fastqc_params           : ${qcparams}
 |   cellranger_params       : ${crparams}
 |   star_params             : ${starparams}
 |   idx_params              : ${idxparams}
 |   filter                  : ${filtering}
 |   minReads                : ${minReads}
 |   chunkSize               : ${chunkSize}
 |   maxDist                 : ${maxDist}
 |   maxDistCaTCH            : ${maxDistCaTCH}
 |   maxDistUMIs             : ${maxDistUMIs}
 |   clusterMethodCaTCH      : ${clusterMethodCaTCH}
 |   clusterMethodUMIs       : ${clusterMethodUMIs}
 |   majorityVote            : ${majorityVote}
 |   uniqueCaTCH             : ${uniqueCaTCH}
 |   min_detected_barcodes   : ${minBC}
 |   singlet_cutoff          : ${singletCutoff}
 |   bc1_cutoff              : ${bc1Cutoff}
 |   bc2_cutoff              : ${bc2Cutoff}
 |   max_mt_percent          : ${maxMtPercent}
 |   min_detected_features   : ${minDetectedFeatures}
 |   hvg_cutoff              : ${hvgCutoff}
 |   pval_cutoff             : ${pvalCutoff}
 |   lfc_cutoff              : ${lfcCutoff}
 |   markerfile              : ${markerfile}
 |   organism                : ${organism}
 |   baseline                : ${refName}
 |   stopOnWarnings          : ${stopOnWarnings}
 |
 ======================================================================
""".stripIndent()

if (runqc){
    log.info("Running QC")
}

if (filtering){
    log.info("Using filtered count output")
}else{
    log.info("Using raw count output")
}

// ---------------------------------------------------------------------
// Pipeline Channels and Processes
// ---------------------------------------------------------------------

// For more information about syntax, please refer to the nextflow documentation at https://www.nextflow.io/docs/latest/index.html
stopOnWarn = (stopOnWarnings) ? "yes" : "no"

/************************************************************************
                SANITY CHECK SAMPLESHEET
************************************************************************/

process check_samplesheet {
    tag "$samplesheet"

    input:
    path samplesheet

    output:
    path 'Valid_*.csv', emit: csv
    path 'sanity_check', emit: check

    script: 
    """
    ${scriptDirPy}checkSampleSheet.py \
        ${samplesheet} \
    && mv ${samplesheet} Valid_${samplesheet}
    """
}


/************************************************************************
                CONCATENATE LANES
************************************************************************/

process concat_lanes {

    tag { sampleName instanceof List ? sampleName.join(',') : sampleName }

    input:
    tuple val(sampleName), path("inputs/R1_?"), path("inputs/R2_?"), val(cells_expected), val(chemistry)

    output:
    tuple val(sampleName), path("${sampleName}_R1.fastq.gz",includeInputs:false), path("${sampleName}_R2.fastq.gz",includeInputs:false), val(cells_expected), val(chemistry)

    script:
    """
    cat inputs/R1_* > "${sampleName}_R1.fastq.gz"
    cat inputs/R2_* > "${sampleName}_R2.fastq.gz"
    """
}


/************************************************************************
                STEP 0: Run read QC
************************************************************************/

process qc_raw{
   
    //conda "$MAPENV"+".yaml"
    //cpus THREADS
	cache 'lenient'
    //label 'big_mem'
    //validExitStatus 0,1
    tag "${sampleName}"

    publishDir "${absDir}/" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "OUTPUT/QC/FASTQC/${file(filename).getName()}"
        else if (filename.indexOf("html") > 0)    "OUTPUT/QC/FASTQC/${file(filename).getName()}"
        else null
    }
    
    input:
    tuple val(sampleName), path(read1), path(read2)

    output:
    path("*.zip"), emit: zip
    path("*.html"), emit: html

    script:
    if (binDir){
        fqc = binDir+"fastqc"
    } else{
        fqc = "fastqc"
    }
    """
    ${fqc} --quiet -t ${task.cpus} $qcparams --noextract -f fastq $read1 $read2 && 
    for fqc in *_fastqc.{zip,html}
    do
        mv "\$fqc" "${sampleName}_\$fqc"
    done
    """
}

process mqc{
    //conda "$MAPENV"+".yaml"
    //cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1
    tag "${sampleName}"

    publishDir "${absDir}/" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "OUTPUT/QC/MULTI/${file(filename).getName()}"
        else if (filename.indexOf("html") > 0)    "OUTPUT/QC/MULTI/${file(filename).getName()}"
        else null
    }

    input:
    path(fastqcs) //, stageAs: "?/*")

    output:
    path "*.zip", includeInputs:false, emit: mqc
    path "*.html", includeInputs:false, emit: html

    script:
    if (binDir != ''){
        mqc = binDir+"multiqc"
    } else{
        fqc = "multiqc"
    }
    """
    touch $fastqcs; export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; ${mqc} -f -k json -z -o \${PWD} .
    """
}

/************************************************************************
                STEP 1: Run CellRanger count
************************************************************************/

process Cellranger_idx{

    //conda "$MAPENV"+".yaml"
    //cpus THREADS
	cache 'lenient'
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${absDir}/" , mode: 'copyNoFollow', overwrite: true,
    saveAs: {filename ->
        if (filename.indexOf("Log.out") > 0)       "OUTPUT/CellRanger/LOGS/${file(filename).getName()}"
        else if (filename.indexOf(".idx") > 0)     "${mapindex}.idx"
        else                                       "${mapindex}"
    }

    input:
    path genome
    path anno

    output:
    path "${file(mapindex).getName()}", emit: idx
    path "${file(mapindex).getName()}*", emit: idx_extra
    path "star.idx", emit: idxlink
    path "*.out", emit: idxlog

    script:
    gen =  file(genome).getName()
    an  = file(anno).getName()
    IDX = file(mapindex).getName()
    taskmem = task.memory.toGiga()

    """
    zcat ${gen} > tmp.fa \
        && zcat ${an} > tmp_anno \
        && cellranger mkref ${idxparams}\
        --genome=${IDX} \
        --nthreads ${task.cpus} \
        --memgb ${taskmem} \
        --fasta tmp.fa \
        --genes tmp_anno \
        && rm -f tmp.fa tmp_anno \
        && ln -s ${IDX}/star star.idx
    """
    //cellranger mkgtf $anno $filt --attribute=gene_biotype:protein_coding &&
}


process runCellrangerCount{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${absDir}/" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("feature_bc_matrix") > 0)       "OUTPUT/CellRanger/${file(filename).getName()}"
        else if (filename.indexOf("projection.csv") >0)      "OUTPUT/CellRanger/${sampleName}/tSNEs/gene_expression_2_components/projection.csv"
        else                                                 "OUTPUT/CellRanger/${file(filename).getName()}"
    }

    //publishDir "outputs/cellranger/", mode: "copy"

    input:
        tuple val(sampleName), path('inputs/'), path('inputs/'), val(cells_expected), val(chemistry), path(index)

    
    output:
        path "${sampleName}", emit: name
        tuple val(sampleName), path("${sampleName}/filtered_feature_bc_matrix/barcodes.tsv.gz"), emit: cell_ids_filtered
        tuple val(sampleName), path("${sampleName}/raw_feature_bc_matrix/barcodes.tsv.gz"), emit: cell_ids_raw
        tuple val(sampleName), path("${sampleName}_filtered_feature_bc_matrix"), emit: cell_data_filtered
        tuple val(sampleName), path("${sampleName}_raw_feature_bc_matrix"), emit: cell_data_raw
    
    script:
    // Check chemistry specific settings
    chemistry = chemistry.unique()[0]
    cells_expected = cells_expected.unique()[0]
    if (!chemistry.contains("10X")){
        log.error("Running CellRanger on chemistry different than 10X is not supported, please check your settings and sample sheet.")

    }
    if (cells_expected != "NA"){
        if (cells_expected.contains("!")){
            crparams = crparams + ' --force-cells ' + cells_expected.replaceAll('!', '')
        }else{
            crparams = crparams + ' --expect-cells ' + cells_expected 
        }
    }
    taskmem = task.memory.toGiga()

    """
    # Find all reads, sort them by name to ensure that the paired files are on the consecutive lines,
    # and then create symlinks with proper names (SampleName_S1_R1_xxx.fastq.gz)
    IDX=1
    BKP=\${IFS}
    IFS=\$'\\n'
    for LINE in \$(find inputs/ -regextype posix-extended -regex ".*R[12][\\._].*fastq.gz" -exec readlink -f {} \\; | sort | paste - -);
    do
        SUFFIX=\$(printf "%03d" \${IDX})

        mkdir -p tomap
        R1=\$(echo \${LINE} | cut -f1|sed 's|*/||g')
        NEW_NAME=${sampleName}_S1_R1_\${SUFFIX}.fastq.gz
        ln -sf \${R1} tomap/\${NEW_NAME}

        R2=\$(echo \${LINE} | cut -f2|sed 's|*/||g')
        NEW_NAME=${sampleName}_S1_R2_\${SUFFIX}.fastq.gz
        ln -sf \${R2} tomap/\${NEW_NAME}

        IDX=\$((IDX + 1))
    done

    IFS=\${BKP}

    cellranger count \
        --disable-ui \
        ${crparams} \
        --jobmode local \
        --localcores ${task.cpus} \
        --localmem ${taskmem} \
        --transcriptome ${index} \
        --id ${sampleName} \
        --fastqs tomap

    mv ${sampleName} rundir
    mv rundir/outs ${sampleName}
    ln -fs ${sampleName}/filtered_feature_bc_matrix ${sampleName}_filtered_feature_bc_matrix
    ln -fs ${sampleName}/raw_feature_bc_matrix ${sampleName}_raw_feature_bc_matrix
    cat ${sampleName}/analysis/tsne/gene_expression_2_components/projection.csv.gz |gzip > ${sampleName}/analysis/tsne/gene_expression_2_components/projection.csv.gz
    """
}


process useCellrangerData{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("feature_bc_matrix") >0)       "OUTPUT/CellRanger/${file(filename).getName()}"
        else if (filename.indexOf("projection.csv") >0)     "OUTPUT/CellRanger/${sampleName}/tSNEs/gene_expression_2_components/projection.csv"
        else                                                "OUTPUT/CellRanger/${file(filename).getName()}"
    }

    input:
        tuple val(sampleName), path("cr_data")

    output:
        path "${sampleName}", emit: name
        tuple val(sampleName), path("${sampleName}/analysis/tsne/gene_expression_2_components/projection.csv"), emit: cell_ids_from_precomputed_filtered
        tuple val(sampleName), path("${sampleName}/raw_feature_bc_matrix/barcodes.tsv.gz"), emit: cell_ids_from_precomputed_raw
        tuple val(sampleName), path("${sampleName}/filtered_feature_bc_matrix"), emit: cell_data_from_precomputed_filtered
        tuple val(sampleName), path("${sampleName}/raw_feature_bc_matrix"), emit: cell_data_from_precomputed_raw        

    script:
        """
        mv cr_data ${sampleName}
        ln -fs ${sampleName}/filtered_feature_bc_matrix ${sampleName}_filtered_feature_bc_matrix
        ln -fs ${sampleName}/raw_feature_bc_matrix ${sampleName}_raw_feature_bc_matrix
        zcat ${sampleName}_raw_feature_bc_matrix/barcodes.tsv.gz > ${sampleName}_raw_feature_bc_matrix/barcodes.tsv 
        """
}

/************************************************************************
                STEP 1 (Optional): Run STARsolo count
************************************************************************/

process star_idx{

    //conda "$MAPENV"+".yaml"
    //cpus THREADS
	cache 'lenient'
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${absDir}/" , mode: 'copyNoFollow', overwrite: true,
    saveAs: {filename ->
        if (filename.indexOf("Log.out") > 0)       "OUTPUT/STAR/LOGS/${file(filename).getName()}"
        else if (filename.indexOf(".idx") > 0)     "${mapindex}.idx"
        else                                       "${mapindex}"
    }

    input:
    path genome
    path anno

    output:
    path "${IDX}", emit: idx
    path "${IDX}.idx", emit: idxlink
    path "*.out", emit: idxlog

    script:
    gen =  genome.getName()
    an  = anno.getName()
    IDX = file(gen).getSimpleName()+'_idx'
    taskmem      = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''   
    """
    zcat ${gen} > tmp.fa && zcat ${an} > tmp_anno && mkdir -p ${IDX} && STAR ${idxparams} --runThreadN ${task.cpus} --runMode genomeGenerate --outTmpDir STARTMP --genomeDir ${IDX} --genomeFastaFiles tmp.fa --sjdbGTFfile tmp_anno ${taskmem} && mv -f ${IDX}/*.out ${IDX}.Log.out && ln -s ${IDX} ${IDX}.idx
    """
}

process create_dummy_whitelist{
    cache 'lenient'
    label 'big_mem'

    output:
    path("dummy.txt"), emit: dummy

    script:
    """
    touch dummy.txt
    """
}

process star_mapping{
    //conda "$MAPENV"+".yaml"
    //cpus THREADS
	cache 'lenient'
    label 'big_mem'
    tag "${sampleName}"
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_unmapped") > 0)          "OUTPUT/STAR/UNMAPPED/"+"${file(filename).getName()}"
        else if (filename.indexOf(".sam.gz") >0)        "OUTPUT/STAR/MAPPED/"+"${file(filename).getName()}"
        else if (filename.indexOf(".bam") >0)           "OUTPUT/STAR/MAPPED/"+"${filename}"
        else if (filename.indexOf(".tab") >0)           "OUTPUT/STAR/MAPPED/"+"${filename}"
        else if (filename.indexOf("Log*.out") >0)       "OUTPUT/STAR/LOGS/${file(filename).getName()}"
        else if (filename.indexOf("*.log") >0)          "OUTPUT/STAR/LOGS/${file(filename).getName()}"
        else if (filename.indexOf("Summary.csv") >0)    "OUTPUT/STAR/SUMMARY/${sampleName}_${file(filename).getName()}"
        else                                            "OUTPUT/STAR/${filename}"
    }

    input:
    tuple val(sampleName), path(r1), path(r2), val(cells_expected), val(chemistry), path(idx), path(whitelist)
    
    output:
    path "${sampleName}", emit: name
    tuple val(sampleName), path("${sampleName}"), emit: out
    tuple val(sampleName), path("${sampleName}_filtered_feature_bc_matrix"), emit: cell_data_filtered
    tuple val(sampleName), path("${sampleName}_raw_feature_bc_matrix"), emit: cell_data_raw
    tuple val(sampleName), path("${sampleName}_filtered_barcodes.tsv.gz"), emit: cell_ids_filtered
    tuple val(sampleName), path("${sampleName}_raw_barcodes.tsv.gz"), emit: cell_ids_raw
    tuple val(sampleName), path("*_mapped.sam.gz"), emit: sam
    tuple val(sampleName), path("*.bam"), emit: bam
    tuple val(sampleName), path("*.bai"), emit: bai
    tuple val(sampleName), path("*.log"), emit: logs
    tuple val(sampleName), path("*Log*.out"), emit: xtralogs
    tuple val(sampleName), path("*.tab"), emit: sjtab
    tuple val(sampleName), path("*_unmapped.fastq.gz", includeInputs:false), emit: unmapped
    //tuple val(sampleName), path("Summary.csv"), emit: qc

    script:
    idxdir = idx.toRealPath()
    extraparams = ''
    
    chemistry = chemistry.unique()[0]
    cells_expected = cells_expected.unique()[0]
    
    // Check whitelist
    if( (whitelist.size() > 0 ) && (chemistry != 'ScaleBio')){
        starparams = starparams + " --soloCBwhitelist ${whitelist}"
    }else{
        starparams = starparams + " --soloCBwhitelist None"
    }
    // Check chemistry specific settings
    if (chemistry.contains("10X")){
        extraparams = params.star_10X
    }else if (chemistry == "Droplet"){
        log.info("Running StarSolo on unspecified Droplet chemistry, please ensure your STARsolo parameters fit the protocol, please check and adapt default settings in mappers.config file in the conf directory of the nextflow subdirectory of this pipeline. There is no automatic sanity check!")
        extraparams = params.star_droplet
    } else if (chemistry == 'ScaleBio'){
        starparams = params.starMappingParamsScale
        extraparams = params.star_scalebio
    } else if (chemistry == "Smart"){
        log.info("Running StarSolo on chemistry different than Droplet based, please ensure your STARsolo parameters fit the protocol, please check and adapt default settings in mappers.config file in the conf directory of the nextflow subdirectory of this pipeline. There is no automatic sanity check!")
        extraparams = params.star_smart
    } else{
        log.error("Unknown chemistry! Please choose between Droplet (10X or ScaleBio) and Smart.")
        exit('Unknown chemistry! Please choose between Droplet (10X or ScaleBio) and Smart.')
    }
    // Check expected cell count
    if (cells_expected != "NA"){
        cellsexpCR = cells_expected.replaceAll('!', '') + ' 0.99 10'
        cellsexpED = cells_expected.replaceAll('!', '') + ' 0.99 10 45000 90000 500 0.01 20000 0.01 10000'
        starparams = starparams.replaceAll('CellRanger2.2', 'CellRanger2.2 '+cellsexpCR).replaceAll('EmptyDrops_CR', 'EmptyDrops_CR '+cellsexpED)
        extraparams = extraparams.replaceAll('CellRanger2.2', 'CellRanger2.2 '+cellsexpCR).replaceAll('EmptyDrops_CR', 'EmptyDrops_CR '+cellsexpED)
    }
    
    // Build params
    starparams = starparams + ' ' + extraparams

    //// Convert extraparams to barcode indices json for postprocessing tools
    //xtra = chemistry_params_to_map(extraparams)
    //def json = new groovy.json.JsonBuilder()
    //json rootKey: xtra
    //xtra =  groovy.json.JsonOutput.prettyPrint(json.toString())

    // Check read order
    if ( starparams.contains('--soloBarcodeMate' )){
        if ( starparams.contains('--soloBarcodeMate 1' )){
            read1 = r2
            read2 = r1
            starparams = starparams.replaceAll('--soloBarcodeMate 1', '' )
        }else if ( starparams.contains('--soloBarcodeMate 2' )){
            read1 = r1
            read2 = r2
            starparams = starparams.replaceAll('--soloBarcodeMate 2', '' )
        }else{
            error('specified --soloBarcodeMate with unknown read')
        }
    }else{
        read1 = r2
        read2 = r1
    }

    //Get specific Feature counts, we use the first feature in --soloFeatures as default
    if ( starparams.contains('--soloFeatures') ){
        sfeature = starparams.split('--soloFeatures')[1].split(' ')[1]
    }else{
        sfeature = 'Gene'
    }

    of = sampleName+'.Aligned.sortedByCoord.out.bam'
    gf = of.replaceAll(/\Q.Aligned.sortedByCoord.out.bam\E/,"_mapped.sam.gz")
    gb = of.replaceAll(/\Q.Aligned.sortedByCoord.out.bam\E/,"_mapped.bam")

    """
    STAR ${starparams} --runThreadN ${task.cpus} --genomeDir ${idxdir} --readFilesCommand zcat --readFilesIn <(cat ${read1}) <(cat ${read2}) --outFileNamePrefix ${sampleName}. --outReadsUnmapped Fastx &&samtools view -h ${of} | gzip > ${gf} && touch ${sampleName}.Unmapped.out.mate1 ${sampleName}.Unmapped.out.mate2 && cat ${sampleName}.Unmapped.out.mate1 | paste - - - - |tr \"\\t\" \"\\n\"| gzip > ${sampleName}_R1_unmapped.fastq.gz && cat ${sampleName}.Unmapped.out.mate2| paste - - - - |tr \"\\t\" \"\\n\"| gzip > ${sampleName}_R2_unmapped.fastq.gz && mv *Log.out ${sampleName}_mapping.log

    mv ${of} ${gb}
    samtools index ${gb}
    mv ${sampleName}.Solo.out ${sampleName}
    gzip -f ${sampleName}/*/filtered/*
    gzip -f ${sampleName}/*/raw/*
    ln -s ${sampleName}/${sfeature}/filtered ${sampleName}_filtered_feature_bc_matrix
    ln -s ${sampleName}/${sfeature}/raw ${sampleName}_raw_feature_bc_matrix
    ln -s ${sampleName}/${sfeature}/filtered/barcodes.tsv.gz ${sampleName}_filtered_barcodes.tsv.gz
    ln -s ${sampleName}/${sfeature}/raw/barcodes.tsv.gz ${sampleName}_raw_barcodes.tsv.gz
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

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename == "Counts")       "OUTPUT/Counts/Chunks/${sampleName}/counts"
        else if (filename == "Reads")   "OUTPUT/Counts/Chunks/${sampleName}/reads"
        else                            "OUTPUT/Counts/Chunks/${sampleName}/${file(filename).getName()}"
    }

    input:
        tuple val(sampleName), path(r1), path(r2), val(chemistry), path(cellIDs)

    output:
        tuple val(sampleName), path('Counts'), path('Reads'), emit: counts_chunks_out

    script:
    // Check chemistry specific settings
    extraparams = ''
    //chemistry = chemistry[0]

    if (chemistry.contains("10X")){
        extraparams = chemistry_params_to_str(params.star_10X)
        log.debug("${sampleName}: Counting 10X chemistry with following parameters: ${extraparams}")
    }else if (chemistry == "Droplet"){
        extraparams = chemistry_params_to_str(params.star_droplet)
        log.debug("${sampleName}: Counting Droplet chemistry with following parameters: ${extraparams}")
    } else if (chemistry == 'ScaleBio'){
        extraparams = chemistry_params_to_str(params.star_scalebio)
        log.debug("${sampleName}: Counting ScaleBio chemistry with following parameters: ${extraparams}")
    } else if (chemistry == "Smart"){
        extraparams = chemistry_params_to_str(params.star_smart)
        log.debug("${sampleName}: Counting Smart chemistry with following parameters: ${extraparams}")
    } else{
        extraparams = ''
        log.debug("${sampleName}: Counting unset chemistry with default 10X parameters: ${extraparams}")
    }

    """
    ${scriptDirPy}countBarcodesInChunks.py \
        --r1 ${r1} \
        --r2 ${r2} \
        --cellIDs ${cellIDs} \
        --counts Counts \
        ${extraparams} \
    | tee log \
    | grep -Po "Read [0-9,]+ single cell entries" \
    | cut -d" " -f2 > Reads 
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

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".sclib") > 0)         "OUTPUT/Counts/libraries/unfiltered/${file(filename).getName()}"
        else if (filename.indexOf(".stats") > 0)    "OUTPUT/Counts/libraries/unfiltered/${file(filename).getName()}"
        else                                        "OUTPUT/Counts/libraries/unfiltered/${file(filename).getName()}"
    }

    input:
        tuple val(sampleName), path("Counts/file*"), path("Reads/file*")

    output:
        tuple val(sampleName), path('*.sclib'), emit: merged_libraries
        tuple val(sampleName), path('*.stats'), emit: unfiltered_stats

    script:
    """
    find Counts -name "file*" > librarieslist
    find Reads -name "file*" > readcountslist

    ${scriptDirPy}mergeChunkCounts.py \
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

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".collapsed.sclib") > 0)         "OUTPUT/Counts/libraries/collapsed/${file(filename).getName()}"
        else if (filename.indexOf(".collapsed.stats") > 0)    "OUTPUT/Counts/libraries/collapsed/${file(filename).getName()}"
        else                                                  "OUTPUT/Counts/libraries/collapsed/${file(filename).getName()}"
    }
    
    input:
        tuple val(sampleName), path(library)

    output:
        tuple val(sampleName), path('*.collapsed.sclib'), emit: collapsed_libraries
        tuple val(sampleName), path('*.collapsed.stats'), emit: collapsed_stats

    script:
    unique = "false"
    if (uniqueCaTCH){
        unique = "true"
    }
    catchMax = maxDistCaTCH ?: maxDist
    umiMax = maxDistUMIs ?: maxDist
    catchMethod = clusterMethodCaTCH ? "--cluster-method-catch ${clusterMethodCaTCH}" : ""
    umiMethod = clusterMethodUMIs ? "--cluster-method-umis ${clusterMethodUMIs}" : ""
    """
    ${scriptDirPy}collapseCaTCHbarcodes.py \
        --library ${library} \
        --maxdist ${maxDist} \
        --maxdist-catch ${catchMax} \
        --maxdist-umis ${umiMax} \
        --minsupport ${minReads} \
        --unique ${unique} \
        ${catchMethod} \
        ${umiMethod} \
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

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".resolved_multiplets.sclib") > 0)         "OUTPUT/Counts/libraries/resolved_multiplets/${file(filename).getName()}"
        else if (filename.indexOf(".resolved_multiplets.stats") > 0)    "OUTPUT/Counts/libraries/resolved_multiplets/${file(filename).getName()}"
        else                                                            "OUTPUT/Counts/libraries/resolved_multiplets/${file(filename).getName()}"
    }
    
    input:
        tuple val(sampleName), path(library)

    output:
        tuple val(sampleName), path("*.resolved_multiplets.sclib"), emit: resolved_multiplets_libraries
        tuple val(sampleName), path("*.resolved_multiplets.stats"), emit: resolved_stats

    script:
    """
    ${scriptDirPy}resolveMultiplets.py \
        --library ${library} \
        --majority ${majorityVote} \
        --outlib ${sampleName}.resolved_multiplets.sclib \
    | tee ${sampleName}.resolved_multiplets.stats
    """
}


/************************************************************************
                    STEP 6: Generate tables
************************************************************************/

process generateTables{
    
    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".CaTCHbarcodes") > 0)       "OUTPUT/Reports/${file(filename).getName()}"
        else if (filename.indexOf(".cells") > 0)          "OUTPUT/Reports/${file(filename).getName()}"
        else                                              "OUTPUT/Reports/${file(filename).getName()}"
    }

    input:
        tuple val(sampleName), path(library)

    output:
        tuple val(sampleName), path('*.CaTCHbarcodes'), emit: report_CaTCHbarcodes
        tuple val(sampleName), path('*.cells'), emit: report_cells

    script:
    """
    ${scriptDirPy}generateOutputTables.py \
        --library ${library} \
        --CaTCH ${sampleName}.CaTCHbarcodes \
        --cells ${sampleName}.cells \
        --unique ${uniqueCaTCH} \
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

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".png") > 0)       "OUTPUT/Reports/plots/${file(filename).getName()}"
        else                                    "OUTPUT/Reports/plots/${file(filename).getName()}"
    }

    input:
        tuple val(sampleName), path(cell_ids), path(unfiltered), path(collapsed), path(resolved)

    output:
        tuple val(sampleName), path("*.png"), emit: analytics_out

    script:
    """
    touch ${sampleName}_dummy.png
    """
}


/************************************************************************
                    STEP 8: Generate SingleCellExperiment/Seurat objects
************************************************************************/

process preprocessSingleCellData{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "Preprocess_SingleCellData_${sampleTag}"

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if(filtering){
            if (filename.indexOf(".rds.gz") > 0)        "OUTPUT/SCE/filtered/${file(filename).getName()}"
            else if (filename.indexOf(".pdf") > 0)      "OUTPUT/Plots/Overview/${file(filename).getName()}"
            else                                        "OUTPUT/SCE/filtered/${file(filename).getName()}"
        } else{
            if (filename.indexOf(".rds.gz") > 0)        "OUTPUT/SCE/raw/${file(filename).getName()}"
            else if (filename.indexOf(".pdf") > 0)      "OUTPUT/Plots/Overview/${file(filename).getName()}"
            else                                        "OUTPUT/SCE/raw/${file(filename).getName()}"
        }
    }

    input:
        tuple val(sampleName), val(sampleTag), path(featureMatrix), path(catchBarcodes), path(gtf)

    output:
        tuple val(sampleTag), path("*_filtered_seurat_sce.rds.gz"), emit: basic_seurat_sce
        tuple val(sampleTag), path("*_filtered_sce.rds.gz"), emit: basic_sce
        tuple val(sampleTag), path("*_unfiltered_seurat_sce.rds.gz"), emit: basic_raw_seurat_sce
        tuple val(sampleTag), path("*_unfiltered_sce.rds.gz"), emit: basic_raw_sce
        tuple val(sampleTag), path("*.pdf"), emit: basic_sce_qc, optional: true

    script:
    if (filtering){
        featurematrix = "${featureMatrix}"
        outname = 'CaTCHseq.prefiltered'
    }else{
        featurematrix = "${featureMatrix}"
        outname = 'CaTCHseq'
    }
    outprefix = "${sampleTag}_${outname}"
    """
    SAMPLES='${sampleName}'
    BCS='${catchBarcodes}'
    FEATURES='${featurematrix}'

    Rscript --vanilla ${scriptDirR}preprocessData.R \
       --sample \$SAMPLES \
       --data10X \$FEATURES \
       --catchBC \$BCS \
       --baseCond ${refName} \
       --annotation ${gtf} \
       --minBC ${minBC} \
       --bc1Cut ${bc1Cutoff} \
       --bc2Cut ${bc2Cutoff} \
       --singletCut ${singletCutoff} \
       --max_mt ${maxMtPercent} \
       --min_features ${minDetectedFeatures} \
       --hvg_cutoff ${hvgCutoff} \
       --out ${outprefix} \
       --marker ${markerfile} \
       --libpath ${scriptDirR}
    """
}


/************************************************************************
                    STEP 9: Generate overview plots
************************************************************************/
process createOverviewPlots{
    
    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleTag}"

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".pdf") > 0)       "OUTPUT/Plots/Overview/${file(filename).getName()}"
        else                                    "OUTPUT/Plots/Overview/${file(filename).getName()}"
    }

    input:
        tuple val(sampleTag), path(sce)

    output:
        tuple val(sampleTag), path("*overview.pdf"), emit: pdf

    script:
    if (filtering){
        outname = "${sampleTag}_CaTCHseq.prefiltered"
    }else{
        outname = "${sampleTag}_CaTCHseq"
    }
    """
    Rscript --vanilla ${scriptDirR}create_overview_plots.R \
        --sce ${sce} \
        --baseCond ${refName} \
        --out ${outname}_overview \
        --format pdf \
        --width 25 \
        --height 10 \
        --libpath ${scriptDirR}
    """
}

/***************************************************************************
                STEP 10 (OPTIONAL): Calculate CaTCH barcode enrichment
****************************************************************************/

process calculateBarcodeEnrichment{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleTag}"

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".pdf") > 0)       "OUTPUT/Plots/${file(filename).getName()}"
        else                                    "OUTPUT/DE/BarCodes/${file(filename).getName()}"
    }

    input:
        tuple val(sampleTag), path(sce)
        val(check)

    output:
        tuple val(sampleTag), path("*.pdf"), emit: pdf, optional: true
        tuple val(sampleTag), path("*.tsv.gz"), emit: tables, optional: true
        tuple val(sampleTag), path("*.rds.gz"), emit: rds

    script:
     if (filtering){
        outname = "${sampleTag}_CaTCHseq.prefiltered"
    }else{
        outname = "${sampleTag}_CaTCHseq"
    }
    """
    Rscript --vanilla ${scriptDirR}identify_de_catch_barcodes.R \
        --sce ${sce} \
        --baseCond ${refName} \
        --plots_per_row 5 \
        --format pdf \
        --width 400 \
        --height 300 \
        --pcut ${pvalCutoff}\
        --fcut ${lfcCutoff} \
        --out ${outname} \
        --libpath ${scriptDirR}
    """
}


/***************************************************************************
                STEP 11 (OPTIONAL): Find DE genes
****************************************************************************/

process identifyDEGenes{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleTag}"

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".pdf") > 0)       "OUTPUT/Plots/${file(filename).getName()}"
        else                                    "OUTPUT/DE/GENES/${file(filename).getName()}"
    }

    input:
        tuple val(sampleTag), path(sce)
        val(check)

    output:
        tuple val(sampleTag), path("*.pdf"), emit: pdf, optional: true
        tuple val(sampleTag), path("*.tsv.gz"), emit: tables, optional: true
        tuple val(sampleTag), path("*.rds.gz"), emit: rds

    script:
     if (filtering){
        outname = "${sampleTag}_CaTCHseq.prefiltered"
    }else{
        outname = "${sampleTag}_CaTCHseq"
    }
    """
    Rscript --vanilla ${scriptDirR}identify_de_genes.R \
        --sce ${sce} \
        --baseCond ${refName} \
        --format pdf \
        --width 400 \
        --height 300 \
        --organism ${organism}\
        --pcut ${pvalCutoff}\
        --fcut ${lfcCutoff} \
        --out ${outname} \
        --libpath ${scriptDirR}
    """
}


/***********************************************************************
                        MAIN WORKFLOW
************************************************************************/

workflow{
    main:

        /**********************************************************
                Validate SampleSheet
        ***********************************************************/

        if (file(libraries).exists()){
            if (file(libraries).isFile()){
                if (libraries.endsWith(".csv")){
                    if (file(libraries).size() == 0){
                        log.error("SampleSheet is empty!")
                    }
                }else{
                    log.error("SampleSheet is without .csv ending!")
                }
            }else{
                log.error("SampleSheet is not a file!")
            }
        }

        check_samplesheet(Channel.fromPath(libraries))

        /**********************************************************
                STEP 0: Prepare Input and Indices and run QC
        ***********************************************************/
        
        Ch_csv = check_samplesheet.out.csv.splitCsv(sep: ",", header: true)
        //Ch_csv.subscribe {  println "CSV: $it"  }

        Ch_csv_GEX_split = Ch_csv.filter { it.LibraryType == "GEX" }.branch{
                raw: (new File(it.R1)).isFile()
                precomputed: (new File(it.R1)).isDirectory()
            }
        //Ch_csv_GEX_split.raw.subscribe {  println "GEX: $it"  }
        
        //Create a list of unique conditions
        Ch_Conditions = Ch_csv.map { row -> row.Condition }.distinct().branch{                 
                multi:  it != refName
                single: it == refName
        }
        //Ch_Conditions.single.subscribe {  println "Single Conditions: $it"  }
        //Ch_Conditions.multi.subscribe {  println "Multi Conditions: $it"  }
        //println(" IDX "+mapindex+" WL "+whitelist)

        if (mapindex){
            if (!file(mapindex).exists()){
                if (mapperbin == 'CellRanger'){
                    Cellranger_idx(Channel.fromPath(mapref), Channel.fromPath(mapanno))
                    Ch_mapping_idx = Cellranger_idx.out.idx
                }else if (mapperbin == 'STAR'){
                    star_idx(Channel.fromPath(mapref), Channel.fromPath(mapanno))
                    Ch_mapping_idx = star_idx.out.idxlink
                }
            }else{
                if (mapperbin == 'CellRanger'){
                    Ch_mapping_idx = Channel.fromPath(mapindex)
                }else if (mapperbin == 'STAR'){
                    Ch_mapping_idx = Channel.fromPath(mapindex+".idx")
                }
            }
        } else {
            Ch_mapping_idx = Channel.empty()
        }
        //Ch_mapping_idx.subscribe {  println "IDX: $it"  }

        if (whitelist){
            if (file(whitelist).exists()){
                Ch_whitelist = Channel.fromPath(whitelist)
            
            } else {         
            create_dummy_whitelist()   
            Ch_whitelist = create_dummy_whitelist.out.dummy
            }
        } else {         
            create_dummy_whitelist()   
            Ch_whitelist = create_dummy_whitelist.out.dummy
        }
        //Ch_whitelist.subscribe {  println "Whitelist: $it"  }

        Ch_map_input = Ch_csv_GEX_split.raw.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), file(row.R1), file(row.R2), row.CellNumber, row.Chemistry ) }.groupTuple(by: 0)
        
        if (mapperbin == 'CellRanger'){
            Ch_map_precomputed = Ch_csv_GEX_split.precomputed.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), file(row.R1)) }.groupTuple(by: 0)            
        }
        //Ch_map_input.subscribe {  println "INPUT: $it"  }

        if(runqc){
            Ch_QC_input = Ch_csv.filter { new File(it.R1).isFile() }.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), file(row.R1), file(row.R2)) }.groupTuple(by: 0)
            //Ch_QC_input.subscribe {  println "QC: $it"  }
            qc_raw(Ch_QC_input)
            mqc(qc_raw.out.zip.collect())//.flatten().filter( ~/.zip/ ))
        }

        /**********************************************************
                STEP 0.1: Concat Lanes
        ***********************************************************/

        concat_lanes(Ch_map_input)

        /**********************************************************
                STEP 1: Count Reads
        ***********************************************************/

        if (mapperbin == 'CellRanger'){
            runCellrangerCount(concat_lanes.out.combine( Ch_mapping_idx ))
            useCellrangerData(Ch_map_precomputed)

            if (filtering){
                Ch_count_input = Ch_csv.filter { it.LibraryType == "CaTCHseq" }.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), file(row.R1), file(row.R2), row.Chemistry ) }.splitFastq(by: chunkSize, file: true, compress: true, pe: true).combine(runCellrangerCount.out.cell_ids_filtered.mix(useCellrangerData.out.cell_ids_from_precomputed_filtered), by: 0)
            }else{
                Ch_count_input = Ch_csv.filter { it.LibraryType == "CaTCHseq" }.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), file(row.R1), file(row.R2), row.Chemistry ) }.splitFastq(by: chunkSize, file: true, compress: true, pe: true).combine(runCellrangerCount.out.cell_ids_raw.mix(useCellrangerData.out.cell_ids_from_precomputed_raw), by: 0)
            }
        }else if (mapperbin == 'STAR'){
            star_mapping(concat_lanes.out.combine( Ch_mapping_idx ).combine( Ch_whitelist ))
            if (filtering){
                Ch_count_input = Ch_csv.filter { it.LibraryType == "CaTCHseq" }.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), file(row.R1), file(row.R2), row.Chemistry ) }.splitFastq(by: chunkSize, file: true, compress: true, pe: true).combine(star_mapping.out.cell_ids_filtered, by: 0)
            }else{
                Ch_count_input = Ch_csv.filter { it.LibraryType == "CaTCHseq" }.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), file(row.R1), file(row.R2), row.Chemistry ) }.splitFastq(by: chunkSize, file: true, compress: true, pe: true).combine(star_mapping.out.cell_ids_raw, by: 0)
            }
            
        }
        //Ch_count_input.subscribe {  println "InputCount: $it"  }
        
        /**************************************************************
                STEP 2: Count CaTCH barcodes in chunks separately
        ***************************************************************/
        
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
                STEP 6: Generate tables
        ***************************************************************/

        generateTables(resolveMultiplets.out.resolved_multiplets_libraries)


        /**************************************************************
                STEP 7: Analytics report
        ***************************************************************/

        if (mapperbin == 'CellRanger'){
            if (filtering){
                Ch_cell_ids = runCellrangerCount.out.cell_ids_filtered.mix(useCellrangerData.out.cell_ids_from_precomputed_filtered)

                Ch_cell_data = runCellrangerCount.out.cell_data_filtered.mix(useCellrangerData.out.cell_data_from_precomputed_filtered)
            } else {
                Ch_cell_ids = runCellrangerCount.out.cell_ids_raw.mix(useCellrangerData.out.cell_ids_from_precomputed_raw)

                Ch_cell_data = runCellrangerCount.out.cell_data_raw.mix(useCellrangerData.out.cell_data_from_precomputed_raw)
            }
        } else if (mapperbin == 'STAR'){            
            if (filtering){
                Ch_cell_ids = star_mapping.out.cell_ids_filtered

                Ch_cell_data = star_mapping.out.cell_data_filtered                
            } else {
                Ch_cell_ids = star_mapping.out.cell_ids_raw

                Ch_cell_data = star_mapping.out.cell_data_raw
            }
        }
    
        Ch_analytics_in =  Ch_cell_ids.join(mergeBarcodesInChunks.out.merged_libraries, by: 0).join(collapseAndFilterBarcodes.out.collapsed_libraries, by: 0).join(resolveMultiplets.out.resolved_multiplets_libraries, by: 0)

        generateAnalyticsPlots(Ch_analytics_in)


        /**************************************************************
                STEP 8: Generate SingleCellExperiment object
        ***************************************************************/

        Ch_preprocess_input = Ch_cell_data.join(generateTables.out.report_cells, by: 0)
                                         .map { sampleName, featureMatrix, catchBarcodes ->
                                             def sampleTag = sampleName.tokenize('_')[0]
                                             tuple(sampleName, sampleTag, featureMatrix, catchBarcodes, file(mapanno))
                                         }

        preprocessSingleCellData(Ch_preprocess_input)

        /**************************************************************
                STEP 9: Generate overview plots
        ***************************************************************/
        createOverviewPlots(preprocessSingleCellData.out.basic_seurat_sce)
        
        /**************************************************************
                STEP 10 & 11 (OPTIONAL): Run DE Analysis for Barcodes and Genes
        ***************************************************************/
        def Ch_multi_conditions = Ch_Conditions.multi.collect()
        
        calculateBarcodeEnrichment(preprocessSingleCellData.out.basic_seurat_sce, Ch_multi_conditions)        
        identifyDEGenes(preprocessSingleCellData.out.basic_seurat_sce, Ch_multi_conditions)        
        
    //emit:
    //createOverviewPlots.out.pdf
    //createBarcodeEnrichmentPlots.out.pdf   
}
