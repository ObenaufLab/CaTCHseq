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

//parse chemistry params into defined paramstring
def chemistry_params_to_str(xtra){
    xtra = xtra.trim().replaceAll(" \n", "&")
    tmpmap = xtra.tokenize("&").collectEntries{ 
               it.split(" ",2).with{ 
                   [ (it[0].replaceAll("--", "")): (it.size()<2) ? null : it[1] ?: null ] 
                }
            }
    map = ["bcStart" : tmpmap["soloCBstart"], "bcLength" : tmpmap["soloCBlen"], "umiStart" : tmpmap["soloUMIstart"], "umiLength" : tmpmap["soloUMIlen"]]
    xtra = map.collect{ k, v -> v ? "--" + k + " " + v : "" }.join(" ")
    return xtra
}

//Params from CL
absDir = workflow.launchDir
scriptDirR = get_always('scriptDirR') ?: '/tools/scripts/R/'
scriptDirPy = get_always('scriptDirPy') ?: '/tools/scripts/python/'
binDir = get_always('binDir') ?: '/usr/local/bin/'
libraries = get_always('libraries')
chunkSize = get_always('chunkSize') ?: 1_000_000
maxDist = get_always('maxDist') ?: 2
minReads = get_always('minReads') ?: 10
majorityVote = get_always('majorityVote') ?: 90
qcparams = get_always('fastqc_params') ?: ''
mapperbin = get_always('mapper') ?: "CellRanger"
runqc = get_always('withQC') ?: null
crparams = get_always('cellranger_params') ?: ""
starparams = get_always('star_params') ?: "--soloFeatures Gene GeneFull SJ Velocyto --soloMultiMappers EM --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outSAMtype BAM SortedByCoordinate --outSAMprimaryFlag AllBestScore"
idxparams = get_always('idx_params') ?: ""
whitelist = get_always('whitelist') ?: null
refName = get_always('refName') ?: "Day0"
mapindex = get_always('index') ?: null
mapref = get_always('reference') ?: null
mapanno = get_always('annotation') ?: null
filtering = get_always('filter') ?: null
organism = get_always('organism') ?: "Human"
max_mt_percent = get_always('max_mt_percent') ?: 10
min_detected_features = get_always('min_detected_features') ?: 500
hvg_cutoff = get_always('hvg_cutoff') ?: 0.1
pval_cutoff = get_always('pcal_cutoff') ?: 0.1
lfc_cutoff = get_always('lfc_cutoff') ?: 1
markerfile = get_always('markers') ?: '/tools/data/R/stagemarkers_xue2020.rds'
reportsDir = get_always('reportsDir') ?: "${absDir}/REPORTS"
outputDir = get_always('outputDir') ?: "${absDir}/scCaTCH_nf_OUTPUT"

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
                                    Condition       condition for this sample (e.g. timepoint, KO, treatment)
                                    Replicate       replicate (even if a single replicate is present, this column
                                                    cannot be missing or be empty)
                                    LibraryType     either GEX or scCaTCH
                                    R1              path to the R1 read
                                    R2              path to the R2 read (if available)
                                    CellNumber      number of expected cells (or NA if not available)
                                    Chemistry       chemistry used (e.g. 10X, DropIn, SmartSeq), this does not set adequate parameters for mappers automatically, make sure you provide them accordingly

      Optional arguments:
        --chunkSize             number of reads per chunk (default: ${chunkSize})
        --index                 Path to mapper index directory (default: ${mapindex}, NEEDS TO BE SET ALSO TO CREATE NEW INDEX, new index will be stored at given path)
        --mapper                Which mapper to run (default: CellRanger, optional: STAR)
        --withQC                Boolean, run FastQC and MultiQC (default: ${runqc})
        --reference             Path to reference fasta.gz for mapper
        --annotation            Path to annotation gtf.gz for mapper
        --whitelist             Path to barcode whitelist
        --fastqc_params         Optional parameters for FASTQC
        --star_params           Optional parameters for STAR mapping
        --idx_params            Optional parameters for STAR index generation
        --filter                Postprocess filtered counts (default: ${filter})
        --organism              Identifier for organism (choice: ["Human", 
        "Mouse"], default: "Human")
        --baseline              Name of reference day/condition (default: ${refName})
        --marker                RDS file of cellcycle markers (default: ${marker}, NamedList with gene names for each stage [G1S, S, G2M, M, MG1, G0], S and G2M are needed)
        --vote                  Number of votes needed for majority voting (default: ${majorityVote})
        --outputDir             specifies the output directory (default: ${outputDir})
        --reportsDir            specifies the reports directory.(default: ${reportsDir})
        --scriptDirR            specifies the path to the R scripts directory, do not change if running with docker, otherwise set to path on sccatch git repo.(default: /tools/scripts/R/ which is valid for docker instance; set to \${absDir}/sccatch/docker/scripts/R/ for instances not running docker)
        --scriptDirPy            specifies the path to the Python scripts directory, do not change if running with docker, otherwise set to path on sccatch git repo.(default: /tools/scripts/python/ which is valid for docker instance; set to \${absDir}/sccatch/docker/scripts/python/ for instances not running docker)
        --binDir            specifies the path to the binary directory, do not change if running with docker, otherwise set to path on sccatch git repo.(default: /usr/bin/local which is valid for docker instance; set to '' for instances not running docker)
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
 |   filter                  : ${filtering}
 |   withQC                  : ${runqc}
 |   whitelist               : ${whitelist}
 |   mapper                  : ${mapperbin}
 |   max_mt_percent          : ${max_mt_percent}
 |   min_detected_features   : ${min_detected_features}
 |   markerfile              : ${markerfile}
 |   chunk size              : ${chunkSize}
 |   organism                : ${organism}
 |   baseline                : ${refName}
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

    publishDir "${absDir}/" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("Log.out") > 0)       "OUTPUT/CellRanger/LOGS/${file(filename).getName()}"
        else if (filename == ${mapindex})          "${mapindex}"
        else                                       "${mapindex}/${file(filename).getName()}"
    }

    input:
    path genome
    path anno

    output:
    path "${file(mapindex).getName()}", emit: idx
    path "${file(mapindex).getName()}*", emit: idx_extra
    path "*.out", emit: idxlog

    script:
    gen =  file(genome).getName()
    an  = file(anno).getName()
    IDX = file(mapindex).getName()
    """
    zcat $gen > tmp.fa && zcat $an > tmp_anno && cellranger mkref --genome=${IDX} --fasta tmp.fa --genes tmp_anno && rm -f tmp.fa tmp_anno
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
        else if (filename.indexOf("projection.csv") >0)          "OUTPUT/CellRanger/${sampleName}/tSNEs/gene_expression_2_components/projection.csv"
        else                                                     "OUTPUT/CellRanger/${file(filename).getName()}"
    }

    //publishDir "outputs/cellranger/", mode: "copy"

    input:
        tuple val(sampleName), path("inputs/R1_*"), path("inputs/R2_*"), val(cells_expected), val(chemistry), path(index)
    
    
    output:
        path "${sampleName}", emit: name
        tuple val(sampleName), path("${sampleName}/analysis/tsne/gene_expression_2_components/projection.csv"), emit: cell_ids_filtered
        tuple val(sampleName), path("${sampleName}/raw_feature_bc_matrix/barcodes.tsv"), emit: cell_ids_raw
        tuple val(sampleName), path("${sampleName}_filtered_feature_bc_matrix"), emit: cell_data_filtered
        tuple val(sampleName), path("${sampleName}_raw_feature_bc_matrix"), emit: cell_data_raw
    
    script:
    // Check chemistry specific settings
    chemistry = chemistry[0]
    cells_expected = cells_expected[0]
    if (chemistry != "10X"){
        log.error("Running CellRanger on chemistry different than 10X is not supported, please check your settings and sample sheet.")

    }
    if (cells_expected != "NA"){
        if (cells_expected.contains("!")){
            crparams = crparams + ' --force-cells ' + cells_expected.replaceAll('!', '')
        }else{
            crparams = crparams + ' --expect-cells ' + cells_expected 
        }
    }
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
        ln -sf \${R1} inputs/\${NEW_NAME}

        R2=\$(echo \${LINE} | cut -f2)
        NEW_NAME=${sampleName}_S1_R2_\${SUFFIX}.fastq.gz
        ln -sf \${R2} inputs/\${NEW_NAME}

        IDX=\$((IDX + 1))
    done

    IFS=\${BKP}

    cellranger count \
        --disable-ui \
        ${crparams} \
        --jobmode local \
        --localcores ${task.cpus} \
        --localmem 96 \
        --transcriptome ${index} \
        --id ${sampleName} \
        --fastqs inputs

    mv ${sampleName} rundir
    mv rundir/outs ${sampleName}
    ln -fs ${sampleName}/filtered_feature_bc_matrix ${sampleName}_filtered_feature_bc_matrix
    ln -fs ${sampleName}/raw_feature_bc_matrix ${sampleName}_raw_feature_bc_matrix
    zcat ${sampleName}_raw_feature_bc_matrix/barcodes.tsv.gz > ${sampleName}_raw_feature_bc_matrix/barcodes.tsv 
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
        else if (filename.indexOf("projection.csv") >0)          "OUTPUT/CellRanger/${sampleName}/tSNEs/gene_expression_2_components/projection.csv"
        else                                                     "OUTPUT/CellRanger/${file(filename).getName()}"
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
        else if (filename.indexOf(".idx") > 0)      "${mapindex}.idx"
        else                                                     "${mapindex}"
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

    """
    zcat $gen > tmp.fa && zcat $an > tmp_anno && mkdir -p $IDX && STAR $idxparams --runThreadN ${task.cpus} --runMode genomeGenerate --outTmpDir STARTMP --genomeDir $IDX --genomeFastaFiles tmp.fa --sjdbGTFfile tmp_anno && mv -f ${IDX}/*.out ${IDX}.Log.out && ln -s ${IDX} ${IDX}.idx
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
        if (filename.indexOf("_unmapped") > 0)       "OUTPUT/STAR/UNMAPPED/"+"${file(filename).getName()}"
        else if (filename.indexOf(".sam.gz") >0)     "OUTPUT/STAR/MAPPED/"+"${filename.replaceAll(/\Q.Aligned.out.sam.gz\E/,"")}_mapped.sam.gz"
        else if (filename.indexOf("Aligned.sortedByCoord.out.bam.bam") >0)     "OUTPUT/STAR/MAPPED/"+"${filename.replaceAll(/\Q.Aligned.sortedByCoord.out.bam\E/,"")}_mapped.bam"
        else if (filename.indexOf(".tab") >0)        "OUTPUT/STAR/MAPPED/"+"${filename}"
        else if (filename.indexOf("Log.out") >0)        "OUTPUT/STAR/LOGS/${file(filename).getName()}"
        else if (filename.indexOf("Summary.csv") >0)        "OUTPUT/STAR/SUMMARY/${sampleName}_${file(filename).getName()}"
        else                                            "OUTPUT/STAR/${filename}"
    }

    input:
    tuple val(sampleName), path(r1), path(r2), val(cells_expected), val(chemistry), path(idx), path(whitelist)
    
    output:
    path "${sampleName}", emit: name
    tuple val(sampleName), path("${sampleName}"), emit: out
    tuple val(sampleName), path("${sampleName}_filtered_feature_bc_matrix"), emit: cell_data_filtered
    tuple val(sampleName), path("${sampleName}_raw_feature_bc_matrix"), emit: cell_data_raw
    tuple val(sampleName), path("${sampleName}/Gene/filtered/barcodes.tsv"), emit: cell_ids_filtered
    tuple val(sampleName), path("${sampleName}/Gene/raw/barcodes.tsv"), emit: cell_ids_raw
    tuple val(sampleName), path("*_mapped.sam.gz"), emit: sam
    tuple val(sampleName), path("*.bam"), emit: bam
    tuple val(sampleName), path("*Log.out"), emit: logs
    tuple val(sampleName), path("*.tab"), emit: sjtab
    tuple val(sampleName), path("*_unmapped.fastq.gz", includeInputs:false), emit: unmapped
    //tuple val(sampleName), path("Summary.csv"), emit: qc

    script:
    idxdir = idx.toRealPath()
    extraparams = ''
    
    chemistry = chemistry[0]
    cells_expected = cells_expected[0]
    
    // Check whitelist
    if( (whitelist.size() >0) && (chemistry != 'ScaleBio')){
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
        starparams = starparams + ' --nExpectedCells ' + cells_expected.replaceAll('!', '')
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

    of = sampleName+'.Aligned.sortedByCoord.out.bam'
    gf = of.replaceAll(/\Q.Aligned.sortedByCoord.out.bam\E/,"_mapped.sam.gz")

    """
    STAR ${starparams} --runThreadN ${task.cpus} --genomeDir ${idxdir} --readFilesCommand zcat --readFilesIn ${read1} ${read2} --outFileNamePrefix ${sampleName}. --outReadsUnmapped Fastx &&samtools view -h ${of} | gzip > ${gf} && touch ${sampleName}.Unmapped.out.mate1 ${sampleName}.Unmapped.out.mate2 && cat ${sampleName}.Unmapped.out.mate1 | paste - - - - |tr \"\\t\" \"\\n\"| gzip > ${sampleName}_R1_unmapped.fastq.gz && cat ${sampleName}.Unmapped.out.mate2| paste - - - - |tr \"\\t\" \"\\n\"| gzip > ${sampleName}_R2_unmapped.fastq.gz && for f in *.Log.*.out; do mv "\$f" "\$(echo "\$f" | sed 's/.Log.*.out/.log/')"; done && \
    mv ${sampleName}.Solo.out ${sampleName} && \
    ln -s ${sampleName}/Gene/filtered ${sampleName}_filtered_feature_bc_matrix && \
    ln -s ${sampleName}/Gene/raw ${sampleName}_raw_feature_bc_matrix
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
        else if (filename == "Reads")          "OUTPUT/Counts/Chunks/${sampleName}/reads"
        else                                                     "OUTPUT/Counts/Chunks/${sampleName}/${file(filename).getName()}"
    }

    input:
        tuple val(sampleName), val(chemistry), path(r1), path(r2), path(cellIDs)

    output:
        tuple val(sampleName), path('Counts'), path('Reads'), emit: counts_chunks_out

    script:
    // Check chemistry specific settings
    extraparams = ''
    chemistry = chemistry[0]

    if (chemistry.contains("10X")){
        extraparams = chemistry_params_to_str(params.star_10X)
        log.info("${sampleName}: Counting 10X chemistry with following parameters: ${extraparams}")
    }else if (chemistry == "Droplet"){
        extraparams = chemistry_params_to_str(params.star_droplet)
        log.info("${sampleName}: Counting Droplet chemistry with following parameters: ${extraparams}")
    } else if (chemistry == 'ScaleBio'){
        extraparams = chemistry_params_to_str(params.star_scalebio)
        log.info("${sampleName}: Counting ScaleBio chemistry with following parameters: ${extraparams}")
    } else if (chemistry == "Smart"){
        extraparams = chemistry_params_to_str(params.star_smart)
        log.info("${sampleName}: Counting Smart chemistry with following parameters: ${extraparams}")
    } else{
        extraparams = ''
        log.info("${sampleName}: Counting unset chemistry with default 10X parameters: ${extraparams}")
    }

    """
    python ${scriptDirPy}countBarcodesInChunks.py \
        --r1 ${r1} \
        --r2 ${r2} \
        --cellIDs ${cellIDs} \
        --counts Counts \
        --bc_features \
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
        if (filename.indexOf(".sclib") > 0)       "OUTPUT/Counts/libraries/unfiltered/${file(filename).getName()}"
        else if (filename.indexOf(".stats") > 0)             "OUTPUT/Counts/libraries/unfiltered/${file(filename).getName()}"
        else                                                     "OUTPUT/Counts/libraries/unfiltered/${file(filename).getName()}"
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

    python ${scriptDirPy}mergeChunkCounts.py \
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
        if (filename.indexOf(".collapsed.sclib") > 0)       "OUTPUT/Counts/libraries/collapsed/${file(filename).getName()}"
        else if (filename.indexOf(".collapsed.stats") > 0)             "OUTPUT/Counts/libraries/collapsed/${file(filename).getName()}"
        else                                                     "OUTPUT/Counts/libraries/collapsed/${file(filename).getName()}"
    }
    
    input:
        tuple val(sampleName), path(library)

    output:
        tuple val(sampleName), path('*.collapsed.sclib'), emit: collapsed_libraries
        tuple val(sampleName), path('*.collapsed.stats'), emit: collapsed_stats

    script:
    """
    python ${scriptDirPy}collapseCaTCHbarcodes.py \
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

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".resolved_multiplets.sclib") > 0)       "OUTPUT/Counts/libraries/resolved_multiplets/${file(filename).getName()}"
        else if (filename.indexOf(".resolved_multiplets.stats") > 0)             "OUTPUT/Counts/libraries/resolved_multiplets/${file(filename).getName()}"
        else                                                     "OUTPUT/Counts/libraries/resolved_multiplets/${file(filename).getName()}"
    }
    
    input:
        tuple val(sampleName), path(library)

    output:
        tuple val(sampleName), path("*.resolved_multiplets.sclib"), emit: resolved_multiplets_libraries
        tuple val(sampleName), path("*.resolved_multiplets.stats"), emit: resolved_stats

    script:
    """
    python ${scriptDirPy}resolveMultiplets.py \
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

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".CaTCHbarcodes") > 0)       "OUTPUT/Reports/${file(filename).getName()}"
        else if (filename.indexOf(".cells") > 0)             "OUTPUT/Reports/${file(filename).getName()}"
        else                                                     "OUTPUT/Reports/${file(filename).getName()}"
    }

    input:
        tuple val(sampleName), path(library)

    output:
        tuple val(sampleName), path('*.CaTCHbarcodes'), emit: report_CaTCHbarcodes
        tuple val(sampleName), path('*.cells'), emit: report_cells

    script:
    """
    python ${scriptDirPy}generateOutputTables.py \
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

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".png") > 0)       "OUTPUT/Reports/plots/${file(filename).getName()}"
        else                                     "OUTPUT/Reports/plots/${file(filename).getName()}"
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
                    STEP 8: Generate SingleCellExperiment object
************************************************************************/

process preprocessSingleCellData{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if(filtering){
            if (filename.indexOf(".rds.gz") > 0)       "OUTPUT/SCE/filtered/${file(filename).getName()}"
            else if (filename.indexOf(".pdf") > 0)       "OUTPUT/Plots/Overview/${file(filename).getName()}"
            else                                     "OUTPUT/SCE/filtered/${file(filename).getName()}"
        } else{
            if (filename.indexOf(".rds.gz") > 0)       "OUTPUT/SCE/raw/${file(filename).getName()}"
            else if (filename.indexOf(".pdf") > 0)       "OUTPUT/Plots/Overview/${file(filename).getName()}"
            else                                     "OUTPUT/SCE/raw/${file(filename).getName()}"
        }
    }

    input:
        path(featureMatrix)
        path(catchBarcodes)
        //path(script)

    output:
        path("scCaTCH*.rds.gz"), emit: basic_sce
        //path("*.sce.prefiltered.tsne.gz"), emit: basic_sce_tsne
        //path("*.sce.prefiltered.metadata.gz"), emit: basic_sce_metadata
        path("*.pdf"), emit: basic_sce_qc, optional: true

    script:
    if (filtering){
        featurematrix = "filtered_feature_bc_matrix"
        outname = 'scCaTCH.prefiltered'
    }else{
        featurematrix = "raw_feature_bc_matrix"
        outname = 'scCaTCH'
    }
    """ 
    SAMPLES=''
    BCS=''
    FEATURES=''
    IDX=1
    for bc in \$(find . -name "*.cells"|cut -d'/' -f2);
    do
        SN=\${bc%*.cells}
        SAMPLES+="\${SN},"
        BCS+="\${bc},"
        FEATURES+="\${SN}_${featurematrix},"
        IDX=\$((IDX + 1))
    done
    
    SAMPLES=\${SAMPLES:0:-1}
    FEATURES=\${FEATURES:0:-1}
    BCS=\${BCS:0:-1}

    Rscript --vanilla ${scriptDirR}preprocessData.R \
       --sample \$SAMPLES \
       --data10X \$FEATURES \
       --catchBC \$BCS \
       --max_mt ${max_mt_percent} \
       --min_features ${min_detected_features} \
       --hvg_cutoff ${hvg_cutoff} \
       --out ${outname} \
       --marker ${markerfile}
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

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".pdf") > 0)       "OUTPUT/Plots/Overview/${file(filename).getName()}"
        else                                     "OUTPUT/Plots/Overview/${file(filename).getName()}"
    }

    input:
        //tuple path(sce), path(script)
        path(sce)

    output:
        path("*overview.pdf"), emit: pdf

    script:
    if (filtering){
        outname = 'scCaTCH.prefiltered'
    }else{
        outname = 'scCaTCH'
    }
    """
    Rscript --vanilla ${scriptDirR}create_overview_plots.R \
        --sce ${sce} \
        --out ${outname}_overview \
        --format pdf \
        --width 25 \
        --height 10
    """
}

process calculateBarcodeEnrichment{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".pdf") > 0)       "OUTPUT/Plots/${file(filename).getName()}"
        else                                     "OUTPUT/DE/BarCodes/${file(filename).getName()}"
    }

    input:
        path(sce)

    output:
        path "*.pdf", emit: pdf
        path "*.tsv.gz", emit: tables
        path "*.rds.gz", emit: rds

    script:
     if (filtering){
        outname = 'scCaTCH.prefiltered'
    }else{
        outname = 'scCaTCH'
    }
    """
    Rscript --vanilla ${scriptDirR}plot_enriched_and_depleted_BCs.R \
        --sce ${sce} \
        --baseCond ${refName} \
        --plots_per_row 5 \
        --format pdf \
        --width 400 \
        --height 300 \
        --pcut ${pval_cutoff}\
        --fcut ${lfc_cutoff} \
        --out ${outname}
    """
}

process identifyDEGenes{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".pdf") > 0)       "OUTPUT/Plots/${file(filename).getName()}"
        else                                     "OUTPUT/DE/GENES/${file(filename).getName()}"
    }

    input:
        path(sce)

    output:
        path "*.pdf", emit: pdf
        path "*.tsv.gz", emit: tables
        path "*.rds.gz", emit: rds

    script:
     if (filtering){
        outname = 'scCaTCH.prefiltered'
    }else{
        outname = 'scCaTCH'
    }
    """
    Rscript --vanilla ${scriptDirR}identify_de_genes.R \
        --sce ${sce} \
        --baseCond ${refName} \
        --format pdf \
        --width 400 \
        --height 300 \
        --organism ${organism}\
        --pcut ${pval_cutoff}\
        --fcut ${lfc_cutoff} \
        --out ${outname}
    """
}


/***********************************************************************
                        MAIN WORKFLOW
************************************************************************/

workflow{
    main:

        /**********************************************************
                STEP 0: Prepare Input and Indices and run QC
        ***********************************************************/
        
        Ch_csv = Channel.fromPath(libraries).splitCsv(sep: "\t", header: true)
        //Ch_csv.subscribe {  println "CSV: $it"  }

        Ch_csv_GEX_split = Ch_csv.filter { it.LibraryType == "GEX" }.branch{
                raw: (new File(it.R1)).isFile()
                precomputed: (new File(it.R1)).isDirectory()
            }
        //Ch_csv_GEX_split.raw.subscribe {  println "GEX: $it"  }
        
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

        if (mapperbin == 'CellRanger'){
            Ch_map_input = Ch_csv_GEX_split.raw.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), file(row.R1), file(row.R2), row.CellNumber, row.Chemistry ) }.groupTuple(by: 0).combine( Ch_mapping_idx )
            
            Ch_map_precomputed = Ch_csv_GEX_split.precomputed.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), file(row.R1)) }
        }else if (mapperbin == 'STAR'){
           Ch_map_input = Ch_csv_GEX_split.raw.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), file(row.R1), file(row.R2), row.CellNumber, row.Chemistry ) }.groupTuple(by: 0).combine( Ch_mapping_idx ).combine( Ch_whitelist )
        }
        //Ch_map_input.subscribe {  println "INPUT: $it"  }

        if(runqc){
            Ch_QC_input = Ch_csv.filter { new File(it.R1).isFile() }.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), file(row.R1), file(row.R2)) }.groupTuple(by: 0)
            //Ch_QC_input.subscribe {  println "QC: $it"  }
            qc_raw(Ch_QC_input)
            mqc(qc_raw.out.zip.collect())//.flatten().filter( ~/.zip/ ))
        }

        /**********************************************************
                STEP 1: Count Reads
        ***********************************************************/

        if (mapperbin == 'CellRanger'){
            runCellrangerCount(Ch_map_input)
            useCellrangerData(Ch_map_precomputed)

            if (filtering){
                Ch_count_input = Ch_csv.filter { it.LibraryType == "scCaTCH" }.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), row.Chemistry, file(row.R1), file(row.R2) ) }.splitFastq(by: chunkSize, file: true, compress: true, pe: true).combine(runCellrangerCount.out.cell_ids_filtered.mix(useCellrangerData.out.cell_ids_from_precomputed_filtered), by: 0).groupTuple(by: 0)
            }else{
                Ch_count_input = Ch_csv.filter { it.LibraryType == "scCaTCH" }.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), row.Chemistry, file(row.R1), file(row.R2) ) }.splitFastq(by: chunkSize, file: true, compress: true, pe: true).combine(runCellrangerCount.out.cell_ids_raw.mix(useCellrangerData.out.cell_ids_from_precomputed_raw), by: 0).groupTuple(by: 0)
            }
        }else if (mapperbin == 'STAR'){
            star_mapping(Ch_map_input)
            if (filtering){
                Ch_count_input = Ch_csv.filter { it.LibraryType == "scCaTCH" }.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), row.Chemistry, file(row.R1), file(row.R2) ) }.splitFastq(by: chunkSize, file: true, compress: true, pe: true).combine(star_mapping.out.cell_ids_filtered, by: 0).groupTuple(by: 0)
            }else{
                Ch_count_input = Ch_csv.filter { it.LibraryType == "scCaTCH" }.map { row -> tuple(row.SampleName.replaceAll("_","-")+'_'+row.Condition.replaceAll("_","-")+'_'+row.Replicate.replaceAll("_","-"), row.Chemistry, file(row.R1), file(row.R2) ) }.splitFastq(by: chunkSize, file: true, compress: true, pe: true).combine(star_mapping.out.cell_ids_raw, by: 0).groupTuple(by: 0)
            }
            
        }
    
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
                STEP 6: Generate reports
        ***************************************************************/

        generateReports(resolveMultiplets.out.resolved_multiplets_libraries)


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

        //Ch_script = Channel.fromPath("${scriptDirR}"+"preprocessData.R")  // This will only work if repo is cloned in absdir, for docker we actually want to use /tools/scripts/R/ NEEDS TO BE OPTION DEPENDENT
        //Ch_preprocess_input = Ch_cell_data.map { sample, data -> data }.combine(generateReports.out.report_cells, by: 0).combine(Ch_script).collect().flatten().collate(4)
        //Ch_preprocess_input.subscribe {  println "SCE: $it"  }
        //preprocessSingleCellData(Ch_preprocess_input)

        preprocessSingleCellData(Ch_cell_data.map { sample, data -> data }.collect(), generateReports.out.report_cells.map { sample, data -> data }.collect())

        /**************************************************************
                STEP 9: Generate overview plots
        ***************************************************************/
        createOverviewPlots(preprocessSingleCellData.out.basic_sce)
        
        /**************************************************************
                STEP 10: Run DE Analysis for Barcodes and Genes
        ***************************************************************/
        
        calculateBarcodeEnrichment(preprocessSingleCellData.out.basic_sce)        
        identifyDEGenes(preprocessSingleCellData.out.basic_sce)
    //emit:
    //createOverviewPlots.out.pdf
    //createBarcodeEnrichmentPlots.out.pdf
    
}
