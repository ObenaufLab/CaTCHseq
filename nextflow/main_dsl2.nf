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
absDir = workflow.launchDir
scriptDirR = get_always('scriptDirR') ?: '/tools/scripts/R'
scriptDirPy = get_always('scriptDirPy') ?: '/tools/scripts/python'
binDir = get_always('BinDir') ?: '/usr/local/bin'
libraries = get_always('libraries')
chunkSize = get_always('chunkSize') ?: 1_000_000
maxDist = get_always('maxDist') ?: 2
minReads = get_always('minReads') ?: 10
majorityVote = get_always('majorityVote') ?: 90
qcparams = get_always('fastqc_params') ?: ''
mapperbin = get_always('mapper') ?: "CellRanger"
starparams = get_always('star_params') ?: "--soloType CB_UMI_Simple --soloStrand Unstranded --soloUMIlen 12 --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloCellFilter EmptyDrops_CR --soloFeatures Gene GeneFull SJ Velocyto --soloMultiMappers EM --soloCBwhitelist None --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outSAMtype BAM SortedByCoordinate --outSAMprimaryFlag AllBestScore"
idxparams = get_always('idx_params') ?: ""
whitelist = get_always('whitelist') ?: null
refName = get_always('refName') ?: "Day0"
mapindex = get_always('index') ?: null
mapref = get_always('reference') ?: null
mapanno = get_always('annotation') ?: null
filter = get_always('filter') ?: true
max_mt_percent = get_always('max_mt_percent') ?: 10
min_detected_features = get_always('min_detected_features') ?: 500
hvg_cutoff = get_always('hvg_cutoff') ?: 0.1
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
                                    Treatment       sample treatment
                                    Replicate       replicate (even if a single replicate is present, this column
                                                    cannot be missing or be empty)
                                    LibraryType     either GEX or scCaTCH
                                    R1              path to the R1 read
                                    R2              path to theOptional parameters for STAR mapping corresponding R2 read  

      Optional arguments:
        --chunkSize             number of reads per chunk (default: ${params.chunkSize})
        --index                 Path to mapper index file (default: ${params.mapindex})
        --mapper                Which mapper to run (default: CellRanger, optional: STAR)
        --reference             Path to reference fasta.gz for STAR
        --annotation            Path to annotation gtf.gz for STAR
        --whitelist             Path to barcode whitelist
        --fastqc_params         Optional parameters for FASTQC
        --star_params           Optional parameters for STAR mapping
        --idx_params            Optional parameters for STAR index generation
        --filter                Postprocess filtered counts (default: ${params.filter})
        --baseline              Name of reference day/condition (default: ${params.refName})
        --vote                  Number of votes needed for majority voting (default: ${params.majorityVote})
        --outputDir             specifies the output directory (default: ${params.outputDir})
        --reportsDir            specifies the reports directory.(default: ${params.reportsDir})
        --scriptDirR            specifies the path to the R scripts directory, do not change if running with docker, otherwise set to path on sccatch git repo.(default: /tools/scripts/R/ which is valid for docker instance; set to \${absDir}/sccatch/docker/scripts/R/ for instances not running docker)
        --scriptDirPy            specifies the path to the Python scripts directory, do not change if running with docker, otherwise set to path on sccatch git repo.(default: /tools/scripts/python/ which is valid for docker instance; set to \${absDir}/sccatch/docker/scripts/python/ for instances not running docker)
        --BinDir            specifies the path to the binary directory, do not change if running with docker, otherwise set to path on sccatch git repo.(default: /usr/bin/local which is valid for docker instance; set to '' for instances not running docker)
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
                STEP 0: Run read QC
************************************************************************/

process qc_raw{
   
    //conda "$MAPENV"+".yaml"
    //cpus THREADS
	cache 'lenient'
    //label 'big_mem'
    //validExitStatus 0,1
    tag "${sampleName}"

    publishDir "${absDir}/" , mode: 'link', overwrite: true,
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "OUTPUT/QC/FASTQC/${file(filename).getName()}"
        else if (filename.indexOf("html") > 0)    "OUTPUT/QC/FASTQC/${file(filename).getName()}"
        else null
    }
    
    input:
    tuple val(sampleName), path(read1), path(read2)

    output:
    path "*.{zip,html}", emit: fastqc_results

    script:
    if (${binDir} != ''){
        fqc = "${binDir}"+"fastqc"
    } else{
        fqc = "fastqc"
    }
    """
    ${fqc} --quiet -t ${task.cpus} $qcparams --noextract -f fastq $read1 $read2
    """
}


process mqc{
    //conda "$MAPENV"+".yaml"
    //cpus THREADS
	cache 'lenient'
    //label 'big_mem'
    //validExitStatus 0,1
    tag "${sampleName}"

    publishDir "${absDir}/" , mode: 'link', overwrite: true,
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "OUTPUT/QC/MULTI/${file(filename).getName()}"
        else if (filename.indexOf("html") > 0)    "OUTPUT/QC/MULTI/${file(filename).getName()}"
        else null
    }

    input:
    path fastqcs

    output:
    path "*.zip", emit: mqc
    path "*.html", emit: html

    script:
    if (${binDir} != ''){
        mqc = "${binDir}"+"multiqc"
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
        else                                                     "OUTPUT/CellRanger/INDICES/${file(filename).getName()}"
    }

    input:
    path genome
    path anno

    output:
    path "file(gen).getSimpleName()"+'_idx', emit: idx
    path "*.out", emit: idxlog

    script:
    gen =  genome.getName()
    an  = anno.getName()
    filt  = file(anno).getSimpleName()+'_filtered.gtf'
    IDX = file(gen).getSimpleName()+'_idx'
    """
    zcat $gen > tmp.fa && zcat $an > tmp_anno && cellranger mkgtf $anno $filt --attribute=gene_biotype:protein_coding && cellranger mkref   --genome=$IDX --fasta tmp.fa --genes=${filt} && rm -f tmp.fa tmp_anno
    """
}


process runCellrangerCount{

    //conda "cellranger.yaml"
    cache 'lenient'
    //label 'big_mem'
    tag "${sampleName}"

    publishDir "${absDir}/" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("filtered_feature_bc_matrix") > 0)       "OUTPUT/CellRanger/${sampleName}/filtered_feature_bc_matrix"
        else if (filename.indexOf("projection.csv") >0)          "OUTPUT/CellRanger/${sampleName}/tSNEs/gene_expression_2_components/projection.csv"
        else                                                     "OUTPUT/CellRanger/${file(filename).getName()}"
    }

    //publishDir "outputs/cellranger/", mode: "copy"

    input:
        tuple val(sampleName), path("inputs/R1_*"), path("inputs/R2_*"), path(index)
    
    
    output:
        path "${sampleName}", emit: name
        tuple val(sampleName), path("${sampleName}/analysis/tsne/gene_expression_2_components/projection.csv"), emit: cell_ids_filtered
        tuple val(sampleName), path("${sampleName}/raw_feature_bc_matrix/barcodes.tsv.gz"), emit: cell_ids_raw
        tuple val(sampleName), path("${sampleName}/filtered_feature_bc_matrix"), emit: cell_data_filtered
        tuple val(sampleName), path("${sampleName}/raw_feature_bc_matrix"), emit: cell_data_raw
    
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
        ln -sf \${R1} inputs/\${NEW_NAME}

        R2=\$(echo \${LINE} | cut -f2)
        NEW_NAME=${sampleName}_S1_R2_\${SUFFIX}.fastq.gz
        ln -sf \${R2} inputs/\${NEW_NAME}

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

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("filtered_feature_bc_matrix") >0)       "OUTPUT/CellRanger/${sampleName}/filtered_feature_bc_matrix"
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
        else                                                     "OUTPUT/STAR/INDICES/${file(filename).getName()}"
    }

    input:
    path genome
    path anno

    output:
    path "${IDX}", emit: idx
    path "*.out", emit: idxlog

    script:
    gen =  genome.getName()
    an  = anno.getName()
    IDX = file(gen).getSimpleName()+'_idx'
    """
    zcat $gen > tmp.fa && zcat $an > tmp_anno && mkdir -p $IDX && STAR $idxparams --runThreadN ${task.cpus} --runMode genomeGenerate --outTmpDir STARTMP --genomeDir $IDX --genomeFastaFiles tmp.fa --sjdbGTFfile tmp_anno && mv -f *.out ${IDX}.Log.out
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
    tuple val(sampleName), path(reads), path(whitelist), path(idx)
    
    output:
    path "${sampleName}", emit: name
    tuple val(sampleName), path("${sampleName}.Solo.out"), emit: out
    tuple val(sampleName), path("${sampleName}.Solo.out/Gene/filtered"), emit: cell_data_filtered
    tuple val(sampleName), path("${sampleName}.Solo.out/Gene/raw"), emit: cell_data_raw
    tuple val(sampleName), path("${sampleName}.Solo.out/Gene/filtered/barcodes.tsv"), emit: cell_ids_filtered
    tuple val(sampleName), path("${sampleName}.Solo.out/Gene/raw/barcodes.tsv"), emit: cell_ids_raw
    tuple val(sampleName), path("*_mapped.sam.gz"), emit: sam
    tuple val(sampleName), path("*.bam"), emit: bam
    tuple val(sampleName), path("*_mapped.sam.gz"), emit: sam
    tuple val(sampleName), path("*Log.out"), emit: logs
    tuple val(sampleName), path("*.tab"), emit: sjtab
    tuple val(sampleName), path("*_unmapped.fastq.gz"), includeInputs:false, emit: unmapped
    tuple val(sampleName), path("Summary.csv"), emit: qc

    script:
    idxdir = idx.toRealPath()
    
    r1 = reads[1]
    fn = file(r1)
    r2 = reads[2]
    if (starparams.contains('--soloBarcodeMate 1')){
        t = r2
        r2 = r1
        r1 = t
    }
    of = fn+'.Aligned.sortedByCoord.out.bam'
    gf = of.replaceAll(/\Q.Aligned.sortedByCoord.out.bam\E/,"_mapped.sam.gz")

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

    STAR ${starparams} --runThreadN ${task.cpus} --genomeDir ${idxdir} --readFilesCommand zcat --readFilesIn ${r1} ${r2} --outFileNamePrefix ${sampleName}. --outReadsUnmapped Fastx && && samtools view -h ${of} | gzip > ${gf} && rm -f ${of} && touch ${fn}.Unmapped.out.mate1 ${fn}.Unmapped.out.mate2 && cat ${fn}.Unmapped.out.mate1 | paste - - - - |tr \"\\t\" \"\\n\"| gzip > ${fn}_R1_unmapped.fastq.gz && at ${fn}.Unmapped.out.mate2| paste - - - - |tr \"\\t\" \"\\n\"| gzip > ${fn}_R2_unmapped.fastq.gz && for f in *.Log.*.out; do mv "\$f" "\$(echo "\$f" | sed 's/.Log.*.out/.log/')"; done
    """

}

workflow MAPPING{
    take: collection

    main:
    checkidx = file(MAPUIDX)
    //collection.filter(~/.fastq.gz/)

    if (checkidx.exists()){
        idxfile = Channel.fromPath(MAPUIDX)
        star_mapping(idxfile.combine(collection))
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        annofile = Channel.fromPath(MAPANNO)
        star_idx(genomefile, annofile)
        star_mapping(star_idx.out.idx.combine(collection))
    }


    emit:
    mapped  = star_mapping.out.maps
    logs = star_mapping.out.logs
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
        tuple val(sampleName), path(r1), path(r2), path(cellIDs)

    output:
        tuple val(sampleName), path('Counts'), path('Reads'), emit: counts_chunks_out

    script:
    """
    countBarcodesInChunks.py \
        --r1 ${r1} \
        --r2 ${r2} \
        --cellIDs ${cellIDs} \
        --counts Counts \
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

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".CaTCHbarcodes") > 0)       "OUTPUT/Reports/${file(filename).getName()}"
        else if (filename.indexOf(".cells") > 0)             "OUTPUT/Reports/${file(filename).getName()}"
        else                                                     "OUTPUT/Reports/${file(filename).getName()}"
    }

    input:
        tuple val(sampleName), path(library)

    output:
        path('*.CaTCHbarcodes'), emit: report_CaTCHbarcodes
        tuple val(sampleName), path('*.cells'), emit: report_cells

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

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".png") > 0)       "OUTPUT/Analytics/plots/${file(filename).getName()}"
        else                                     "OUTPUT/Analytics/plots/${file(filename).getName()}"
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
        if (filename.indexOf(".rda") > 0)       "OUTPUT/SCE/filtered/${file(filename).getName()}"
        else                                     "OUTPUT/SCE/filtered/${file(filename).getName()}"
    }

    input:
        tuple val(sampleName), path(featureMatrix), path(catchBarcodes), path(script)

    output:
        path("*.sce.prefiltered.rda.gz"), emit: basic_sce
        path("*.sce.prefiltered.tsne.gz"), emit: basic_sce_tsne
        path("*.sce.prefiltered.metadata.gz"), emit: basic_sce_metadata
        path("*.pdf"), emit: basic_sce_qc

    script:
    """
    Rscript --vanilla ${script} \
       --sample ${sampleName} \
       --data10X ${featureMatrix} \
       --catchBC ${catchBarcodes} \
       --max_mt ${max_mt_percent} \
       --min_features ${min_detected_features} \
       --hvg_cutoff ${hvg_cutoff} \
       --out ${sampleName}.sce.prefiltered
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
        if (filename.indexOf(".pdf") > 0)       "OUTPUT/Plots/${file(filename).getName()}"
        else                                     "OUTPUT/Plots/${file(filename).getName()}"
    }

    input:
        tuple path(sce), path(script)

    output:
        path("*overview.pdf"), emit: pdf

    script:
    """
    Rscript --vanilla ${script} \
        --sce ${sce} \
        --out ${sce}_overview \
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

    publishDir "${absDir}/", mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".jpeg") > 0)       "OUTPUT/Plots/${file(filename).getName()}"
        else                                     "OUTPUT/DE/DESEQ2/${file(filename).getName()}"
    }

    input:
        tuple path(sce), path(script)

    output:
        path "*.jpeg", emit: jpeg
        path "*.tsv.gz", emit: tables

    script:
    """
    Rscript --vanilla ${script} \
        --sce ${sce} \
        --baseCond ${refName} \
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
                STEP 0: Prepare Input and Indices and run QC
        ***********************************************************/
        
        Ch_csv = Channel.fromPath(libraries).splitCsv(sep: "\t", header: true)

        Ch_csv_GEX_split = Ch_csv.filter { it.LibraryType == "GEX" }.branch{
                raw: (new File(it.R1)).isFile()
                precomputed: (new File(it.R1)).isDirectory()
            }

        if (mapindex != null){
            if (!mapindex.exists()){
                if (mapperbin == 'CellRanger'){
                    Cellranger_idx(Channel.fromPath(mapref), Channel.fromPath(mapanno))
                    Ch_mapping_idx = Cellranger_idx.out.idx
                }else if (mapperbin == 'STAR'){
                    star_idx(Channel.fromPath(mapref), Channel.fromPath(mapanno))
                    Ch_mapping_idx = star_idx.out.idx
                }
            } else{
                Ch_mapping_idx = Channel.fromPath(mapindex)
            }
        } else {
            Ch_mapping_idx = Channel.empty()
        }

        if (whitelist != null){
            if (whitelist.exists()){
                Ch_whitelist = Channel.fromPath(whitelist)
            }
        } else{
            Ch_whitelist = Channel.empty()
        }

        Ch_QC_input = Ch_csv.filter { new File(it.R1).isFile() }.map { row -> tuple(row.SampleName+'_'+row.Replicate, file(row.R1), file(row.R2)) }.groupTuple(by: 0)
        qc_raw(Ch_QC_input)
        mqc(qc_raw.out.fastqc_results.collect())

        /**********************************************************
                STEP 1: Count Reads
        ***********************************************************/

        if (mapperbin == 'CellRanger'){
            Ch_cellranger_input = Ch_csv_GEX_split.raw.map { row -> tuple(row.SampleName+'_'+row.Replicate, file(row.R1), file(row.R2)) }.groupTuple(by: 0).combine( Ch_mapping_idx )

            runCellrangerCount(Ch_cellranger_input)

            Ch_cellranger_precomputed = Ch_csv_GEX_split.precomputed.map { row -> tuple(row.SampleName+'_'+row.Replicate, file(row.R1)) }
        
            useCellrangerData(Ch_cellranger_precomputed)

            if (filter){
                Ch_count_input = Ch_csv.filter { it.LibraryType == "scCaTCH" }.map { row -> tuple(row.SampleName+'_'+row.Replicate, file(row.R1), file(row.R2)) }.splitFastq(by: chunkSize, file: true, compress: true, pe: true).combine(runCellrangerCount.out.cell_ids_filtered.mix(useCellrangerData.out.cell_ids_from_precomputed_filtered), by: 0)
            }else{
                Ch_count_input = Ch_csv.filter { it.LibraryType == "scCaTCH" }.map { row -> tuple(row.SampleName+'_'+row.Replicate, file(row.R1), file(row.R2)) }.splitFastq(by: chunkSize, file: true, compress: true, pe: true).combine(runCellrangerCount.out.cell_ids_raw.mix(useCellrangerData.out.cell_ids_from_precomputed_raw), by: 0)
            }


        }else if (mapperbin == 'STAR'){
            Ch_star_input = Ch_csv_GEX_split.raw.map { row -> tuple(row.SampleName+'_'+row.Replicate, file(row.R1), file(row.R2)) }.groupTuple(by: 0).combine( Ch_whitelist.combine( Ch_mapping_idx ))
            star_mapping(Ch_star_input)
            if (filter){
                Ch_count_input = Ch_csv.filter { it.LibraryType == "scCaTCH" }.map { row -> tuple(row.SampleName+'_'+row.Replicate, file(row.R1), file(row.R2)) }.splitFastq(by: chunkSize, file: true, compress: true, pe: true).combine(star.mapping.out.cell_ids_filtered)
            }else{
                Ch_count_input = Ch_csv.filter { it.LibraryType == "scCaTCH" }.map { row -> tuple(row.SampleName+'_'+row.Replicate, file(row.R1), file(row.R2)) }.splitFastq(by: chunkSize, file: true, compress: true, pe: true).combine(star.mapping.out.cell_ids_raw)
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
            if (filter){
                Ch_cell_ids = runCellrangerCount.out.cell_ids_filtered.mix(useCellrangerData.out.cell_ids_from_precomputed_filtered)

                Ch_cell_data = runCellrangerCount.out.cell_data_filtered.mix(useCellrangerData.out.cell_data_from_precomputed_filtered)
            } else {
                Ch_cell_ids = runCellrangerCount.out.cell_ids_raw.mix(useCellrangerData.out.cell_ids_from_precomputed_raw)

                Ch_cell_data = runCellrangerCount.out.cell_data_raw.mix(useCellrangerData.out.cell_data_from_precomputed_raw)
            }
        } else if (mapperbin == 'STAR'){            
            if (filter){
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

        Ch_script = Channel.fromPath("${absDir}/sccatch/docker/scripts/R/preprocessData.R")  // This will only work if repo is cloned in absdir, for docker we actually want to use /tools/scripts/R/ NEEDS TO BE OPTION DEPENDENT
        Ch_preprocess_input = Ch_cell_data.combine(generateReports.out.report_cells, by: 0).combine(Ch_script)

        preprocessSingleCellData(Ch_preprocess_input)

        /**************************************************************
                STEP 9: Generate overview plots
        ***************************************************************/

        createOverviewPlots(preprocessSingleCellData.out.basic_sce.combine(Channel.fromPath("${absDir}/sccatch/docker/scripts/R/create_overview_plots.R")))
        createBarcodeEnrichmentPlots(preprocessSingleCellData.out.basic_sce.combine(Channel.fromPath("${absDir}/sccatch/docker/scripts/R/plot_enriched_and_depleted_BCs.R")))        
        
    /*
    emit:
    createOverviewPlots.out.pdf
    createBarcodeEnrichmentPlots.out.jpeg
    */
}
