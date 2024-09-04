#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Analysis of ABO using Oxford nanopore minION sequencing data.
 *
 * Authors:
 * - Fredrick Mobegi <fredrickmaati@gmail.com>
 * - Ben Matern <ben.matern@gmail.com>
 * - Mathijs Groeneweg <m.groeneweg@mumc.nl>
 */

/*
* Include resources
*/

include {get_file; is_null; param_unexpected_error} from './modules/cli'
include {check_env; publish_software} from './modules/versions'
/* 
* use function as function_newname for functions used in more than one process: 
* download as download_database;
* download as download_reference;
*/

include {
    run_fastqc;
	process_exon6;
	process_exon7;
	run_multiqc;
    compile_results;
    rename_samples;
} from './modules/processes'

// MultiQC reporting
def multiqc_report = []

// Chanel for multiqc config/yaml files 
ch_multiqc_config = Channel.fromPath(file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true))
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

// // Chanel for sample deobfuscation
// ch_deobfuscation = Channel.fromPath(file(params.renaming_file, checkIfExists: true))

def helpMessage() {
    log.info """
        Analysis of ABO using Oxford nanopore minION sequencing data
        ======================================================================
        Usage:
        nextflow run main.nf \\
			--reads 'input_fastq/*.fastq' \\
			--exonsix 'exon6_reference.fasta' \\
			--exonseven 'exon7_reference.fasta' \\
            --database 'ABO_Database.fasta' \\
            --outdir '/dapath/to/outdir'

        -------------------------------------------------------------------------------------------------------------
        Mandatory/recommended arguments:
        -------------------------------------------------------------------------------------------------------------
            --reads [file]                Path to input fastq data (must be surrounded with quotes)
            --exonsix [str]               Name of exon6 reference sequence
            --exonseven [str]             Name of exon7 reference sequence
            --database [str]              ABO types reference database (ABO_Database.fasta)
            --outdir [path]               Output directory (can be an absolute path)
			
        -------------------------------------------------------------------------------------------------------------
        Other options:
        -------------------------------------------------------------------------------------------------------------
            -w/--work-dir [file]          The temporary directory for intermediate files
            --help                        This usage statement.

    """.stripIndent()
}

// Show help
if (params.help) {
    helpMessage()
    exit 0
}


def licenseMessage() {
    log.info"""
    ## License
    This code is released under the Apache 2.0 license.
    It is free to use and modify for your own perposes.
	
	## Copyright
	Copyright (c) 2023 | Fredrick Mobegi | PathWest Laboratory Medicine WA |.
	This work is for research reference purposes only and it may contain links to legally privileged data. 
	Except as permitted by the copyright law applicable to you, you may not reproduce or 
	communicate any of the content herein, including files downloadable from the linked page(s), 
	without written permission of the copyright owner(s) or senior PathWest management personnel.
    """.stripIndent()
    log.info ""
}


def versionMessage() {
    log.info "${manifest.version}"
}


def contactMessage() {
    log.info"""
    ## Contact us

    Please raise an issue on github and we will address it ASAP.
      https://github.com/fmobegi/abo-analysis
    """.stripIndent()
    log.info ""
}

process publish_it {

    label "cpu_low"
    label "memory_low"
    label "time_short"
    label "posix"

    publishDir "${params.outdir}", saveAs: { name }

    input:
    tuple val(name), path("infile")

    output:
    path "infile", includeInputs: true

    script:
    """
    """
}


workflow check_duplicates {

    take:
    nostrip
    reads

    main:
    if ( ! nostrip ) {
        reads
            .map { f -> [f.baseName, f] }
            .toList()
            .map { li ->
                // Find any duplicated basenames
                duplicates = li.countBy { it[0] }.findResults { it.value > 1 ? it.key : null }

                // Find the filenames with the duplicated basenames.
                duplicated = li.findAll { duplicates.contains(it[0]) }.collect { n, f -> f }

                // Raise an error if there are duplicates.
                if (duplicated.size() > 0) {
                    log.error(
                        "Some filenames are duplicated after stripping the file extensions.\n" +
                        "This names are used as output filenames.\n" +
                        "The offending filenames are: ${duplicated}\n" +
                        "Please either rename these files or use the '--nostrip' option to disable the extension stripping for output filenames."
                    )
                    exit 1
                }
            }
    }

}

workflow validate_input {

    main:
	//-------------- Prepare reads for mapping ------------------------------
    if ( params.reads ) {
        Channel
			.fromPath(params.reads, checkIfExists: true, type: "file")
			.set { reads_ch }
    }

    else {
        log.error "No sequencing reads supplied. Please provide data using the `--reads` parameter"
        exit 1
    }
	
	//-------------- Check for duplicated filenames--------------------------
 	dups = check_duplicates(params.nostrip, reads_ch)

    //------ Set exon6 reference sequence ------------------------
    if ( params.exonsix ) {
        Channel
            .fromPath(params.exonsix, checkIfExists: true, type: "file")
            .set { exonsix_ch }
	}
    
    else {
        log.error "Please provide a reference sequence for exon 6!!"
        exit 1
    }

    //------ Set exon7 reference sequence ------------------------
    if ( params.exonseven ) {
        Channel
            .fromPath(params.exonseven, checkIfExists: true, type: "file")
            .set { exonseven_ch }
	}
    
    else {
        log.error "Please provide a reference sequence for exon 7!!"
        exit 1
    }
    
    //------ Set ABO_Database.fasta sequence ------------------------
    if ( params.database ) {
        Channel
            .fromPath(params.database, checkIfExists: true, type: "file")
            .set { database_ch }
	}
    
    else {
        log.error "Please provide a reference database for ABO types!!"
        exit 1
    }
	
    emit:
        reads_ch
        exonsix_ch
        exonseven_ch
        database_ch
}


workflow.onComplete {
    /*final logs*/
    if (workflow.success) {
        log.info "\n---------------------------\nAnalysis Completed"
        log.info "---------------------------\nPipeline execution summary\n---------------------------"
    }

    File woc = new File("${params.outdir}/workflow.oncomplete.txt")
    Map endSummary = [:]
	endSummary['Version']        = 			workflow.manifest.version
    endSummary['RunName']        = 			custom_runName ?: workflow.runName
    endSummary['Completed on']   = 			workflow.complete
    endSummary['Duration']       = 			workflow.duration
    endSummary['Success']        = 			workflow.success
    endSummary['Exit status']    = 			workflow.exitStatus
    endSummary['Error report']   = 			workflow.errorReport ?: '-'
    endSummary['Error message']  =			(workflow.errorMessage ?: '-')
    endSummary['Commandline']    =          workflow.commandLine
    endSummary['Project Dirrectory'] = workflow.projectDir
    endSummary['summary'] = summary
    endSummary['Date Started'] = workflow.start
    endSummary['Date Completed'] = workflow.complete
    endSummary['Pipeline script file path'] = workflow.scriptFile
    endSummary['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) endSummary['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) endSummary['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) endSummary['Pipeline Git branch/tag'] = workflow.revision
    String endWfSummary = endSummary.collect { k,v -> "${k.padRight(40, '.')}: $v" }.join("\n")
    //println endWfSummary
    String execInfo = "${fullSum}\nExecution summary\n${logSep}\n${endWfSummary}\n${logSep}\n"
    woc.write(execInfo)
}


// Main Workflow
workflow {
	// Welcome message
	println """\
	Analysis of ABO using Oxford nanopore minION sequencing data
	======================================================================

	User parameters provided for this run are:
	----------------------------------------------------
	Mandatory/recommended arguments
	----------------------------------------------------
	Reference fasta files and samples fastq files:
	 exon6          :         ${params.exonsix}
	 exon7          :         ${params.exonseven}
	 reads          :         ${params.reads}
	 database       :         ${params.database}

	----------------------------------------------------
	Other options
	----------------------------------------------------  
	 outdir         :         ${params.outdir}
	""".stripIndent()
	
	main:
    if ( params.help ) {
        helpMessage()
        licenseMessage()
        contactMessage()
        exit 0
    }

    if ( params.license ) {
        licenseMessage()
        contactMessage()
        exit 0
    }

    if ( params.version ) {
        versionMessage()
        exit 0
    }
		
	// Validate user input and format any files if required required etc.
	input = validate_input()

	// Check that all of the software is installed and finds the version info where possible.
    versions = check_env()
    
    // Collect fastq reads for mapping and fastqc 
	if ( params.reads ) {
			reads = Channel
			.fromPath(params.reads, checkIfExists: true, type: "file")
			.map {file -> [file.simpleName, file]}
			.tap { readsForMapping }
	}
    
	// To gather all QC reports for MultiQC
    ch_reports  = Channel.empty()
	
	// Generate QC summary using FastQC
	run_fastqc(reads)
	
	// Collect reports from FastQC
	ch_reports  = ch_reports.mix(run_fastqc.out.zip.collect{it[1]}.ifEmpty([]))

    // Create chanel to serialize SNP reports
    ch_SNP_reports  = Channel.empty()
	
	// Process reads using exon 6 reference 
    exon6_map = process_exon6(reads, params.exonsix, params.database)
    ch_SNP_reports = ch_SNP_reports.mix(exon6_map.txt.collect().ifEmpty([]))
    
    // Process reads using exon 7 reference 
    exon7_map = process_exon7(reads, params.exonseven, params.database)
	ch_SNP_reports = ch_SNP_reports.mix(exon7_map.txt.collect().ifEmpty([]))

	// Create final results
    // ch_SNP_reports.collect().view() // check chanel files
    compile_ch = compile_results(params.outdir, ch_SNP_reports.collect())

    // chanel with final export file
    ch_export_file = Channel.empty()
                        .mix(compile_ch.final_export.collect().ifEmpty([]))
    
    // Rename order # with grid number if samples deobfuscation file in provided
    if (!params.skip_renaming) {
        rename_samples(params.renaming_file, ch_export_file.collect())
    }

  	// Publish software versions to text
	publish_software()
	ch_versions = publish_software.out.txt

    // Generate samples mapping for samtools output files
    // samples_mapping(multiqc_config, params.outdir)
    // ch_multiqc_config = Channel.empty().mix(samples_mapping.out.yaml.collect{it[1]}.ifEmpty([]))

    // MultiQC report
	if (!params.skip_multiqc){
		ch_multiqc_files = Channel.empty()
			.mix(ch_versions.collect(),
				 ch_multiqc_custom_config.collect().ifEmpty([]),
				 ch_reports.collect(),
				 ch_multiqc_config)

        run_multiqc(ch_multiqc_files.collect(), params.logo)
		multiqc_report = run_multiqc.out.report.toList()
    }

}
