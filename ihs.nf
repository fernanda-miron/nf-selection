#!/usr/bin/env nextflow

/*================================================================
The Aguilar Lab presents...

The iHS calculation pipeline

- A tool for the automated computing of iHS

==================================================================
Version: 0.1
Project repository:
==================================================================
Authors:

- Bioinformatics Design
 Israel Aguilar-Ordonez (iaguilaror@gmail.com)
 María Fernanda Mirón Toruño (fernandamiront@gmail.com)

- Bioinformatics Development
 María Fernanda Mirón Toruño (fernandamiront@gmail.com)

- Nextflow Port
 Israel Aguilar-Ordonez (iaguilaror@gmail.com)
 María Fernanda Mirón Toruño (fernandamiront@gmail.com)

=============================
Pipeline Processes In Brief:

Pre-processing:
vcf_phasing
genetic_map_generation

Core-processing:
ihs_computing

Pos-processing:

Analysis:

================================================================*/

/* Define the help message as a function to call when needed *//////////////////////////////
def helpMessage() {
	log.info"""
  ==========================================
  The iHS calculation pipeline
  v${version}
  ==========================================

	Usage:

	nextflow run ${pipeline_name}.nf --input_ihs <path to ihs design file> --input_pbs <path to pbs design file> [--output_dir path to results ]

	  --input_ihs	<- A design file with the path's required for ihs

		--input_pbs	<- A design file with the path's required for pbs

	  --output_dir  <- directory where results, intermediate and log files will bestored;
	      default: same dir where vcf files are

	  -resume	   <- Use cached results if the executed project has been run before;
	      default: not activated
	      This native NF option checks if anything has changed from a previous pipeline execution.
	      Then, it resumes the run from the last successful stage.
	      i.e. If for some reason your previous run got interrupted,
	      running the -resume option will take it from the last successful pipeline stage
	      instead of starting over
	      Read more here: https://www.nextflow.io/docs/latest/getstarted.html#getstart-resume
	  --help           <- Shows Pipeline Information
	  --version        <- Show version
	""".stripIndent()
}

/*//////////////////////////////
  Define pipeline version
  If you bump the number, remember to bump it in the header description at the begining of this script too
*/
version = "0.1"

/*//////////////////////////////
  Define pipeline Name
  This will be used as a name to include in the results and intermediates directory names
*/
pipeline_name = "nf_ihs_computing"

/*
  Initiate default values for parameters
  to avoid "WARN: Access to undefined parameter" messages
*/
params.input_ihs = false // default is false to not trigger help message
params.input_pbs = false // default is false to not trigger help message
params.help = false //default is false to not trigger help message automatically at every run
params.version = false //default is false to not trigger version message automatically at every run

/*//////////////////////////////
  If the user inputs the --help flag
  print the help message and exit pipeline
*/
if (params.help){
	helpMessage()
	exit 0
}

/*//////////////////////////////
  If the user inputs the --version flag
  print the pipeline version
*/
if (params.version){
	println "${pipeline_name} v${version}"
	exit 0
}

/*//////////////////////////////
  Define the Nextflow version under which this pipeline was developed or successfuly tested
  Updated by iaguilar at MAY 2021
*/
nextflow_required_version = '20.01.0'
/*
  Try Catch to verify compatible Nextflow version
  If user Nextflow version is lower than the required version pipeline will continue
  but a message is printed to tell the user maybe it's a good idea to update her/his Nextflow
*/
try {
	if( ! nextflow.version.matches(">= $nextflow_required_version") ){
		throw GroovyException('Your Nextflow version is older than Pipeline required version')
	}
} catch (all) {
	log.error "-----\n" +
			"  This pipeline requires Nextflow version: $nextflow_required_version \n" +
      "  But you are running version: $workflow.nextflow.version \n" +
			"  The pipeline will continue but some things may not work as intended\n" +
			"  You may want to run `nextflow self-update` to update Nextflow\n" +
			"============================================================"
}

/*//////////////////////////////
  INPUT PARAMETER VALIDATION BLOCK
*/

/* Check if the input directory is provided
    if it was not provided, it keeps the 'false' value assigned in the parameter initiation block above
    and this test fails
*/
 if ( !params.input_ihs | !params.input_pbs ) {
  log.error " Please provide the --input_ihs AND --input_pbs \n\n" +
  " For more information, execute: nextflow run nf_ihs_computing --help"
  exit 1
}

/*
Output directory definition
Default value to create directory is the parent dir of --input_dir
*/
params.output_dir = file(params.input_ihs).getParent() //!! maybe creates bug, should check

/*
  Results and Intermediate directory definition
  They are always relative to the base Output Directory
  and they always include the pipeline name in the variable pipeline_name defined by this Script

  This directories will be automatically created by the pipeline to store files during the run
*/
results_dir = "${params.output_dir}/${pipeline_name}-results/"
intermediates_dir = "${params.output_dir}/${pipeline_name}-intermediate/"

/*
Useful functions definition
*/

/*//////////////////////////////
  LOG RUN INFORMATION
*/
log.info"""
==========================================
The iHS computing pipeline
v${version}
==========================================
"""
log.info "--Nextflow metadata--"
/* define function to store nextflow metadata summary info */
def nfsummary = [:]
/* log parameter values beign used into summary */
/* For the following runtime metadata origins, see https://www.nextflow.io/docs/latest/metadata.html */
nfsummary['Resumed run?'] = workflow.resume
nfsummary['Run Name']			= workflow.runName
nfsummary['Current user']		= workflow.userName
/* string transform the time and date of run start; remove : chars and replace spaces by underscores */
nfsummary['Start time']			= workflow.start.toString().replace(":", "").replace(" ", "_")
nfsummary['Script dir']		 = workflow.projectDir
nfsummary['Working dir']		 = workflow.workDir
nfsummary['Current dir']		= workflow.launchDir
nfsummary['Launch command'] = workflow.commandLine
log.info nfsummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "\n\n--Pipeline Parameters--"
/* define function to store nextflow metadata summary info */
def pipelinesummary = [:]
/* log parameter values beign used into summary */
pipelinesummary['iHS: not phased']			= params.notphased
pipelinesummary['iHS: cutoff']			= params.cutoff
pipelinesummary['iHS: maff']			= params.maff
pipelinesummary['iHS: genetic map']			= params.genetic_map
pipelinesummary['PBS: mart_annotation']			= params.mart
pipelinesummary['iHS: mart_annotation']			= params.imart
pipelinesummary['iHS: merge value']			= params.imerged
pipelinesummary['PBS: merge value']			= params.pmerged
pipelinesummary['Input data']			= params.input_ihs
pipelinesummary['Input data']			= params.input_pbs
pipelinesummary['Results Dir']		= results_dir
pipelinesummary['Intermediate Dir']		= intermediates_dir
/* print stored summary info */
log.info pipelinesummary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "==========================================\nPipeline Start"

/*//////////////////////////////
  PIPELINE START
*/
/* Activar modo DSL2*/
nextflow.enable.dsl=2

/*
	READ INPUTS
*/

/* Load files  into channel */
/* Load files from iHS design file for not phased vcfs*/
if (params.notphased) {
	Channel
	    .fromPath(params.input_ihs)
	    .splitCsv(header:true)
			.map{ row -> [ row.chromosome, file(row.path_vcf), file(row.path_hap), file(row.path_legend), file(row.path_sample), file(row.path_genetic_map), file(row.path_strand_exclude)] }
	    .set{ samples_ihs}
}

/* Load files from iHS design file for genetic map generation */
if (!params.genetic_map) {
	Channel
			.fromPath(params.input_ihs)
			.splitCsv(header:true)
			.map{ row -> [ row.chromosome, file(row.path_legend), file(row.path_sample), file(row.path_genetic_map)] }
			.set{ samples_genetic_map}
}

/* Load files from iHS design file for phased vcfs*/
if (!params.notphased)
Channel
		.fromPath(params.input_ihs)
		.splitCsv(header:true)
		.map{ row -> [ row.chromosome, file(row.path_vcf)] }
		.set{ samples_phased}

/* Load files from iHS design file for phased vcfs*/
if (params.genetic_map) {
	Channel
			 .fromPath(params.input_ihs)
			 .splitCsv(header:true)
			 .map{ row -> [ row.chromosome, file(row.path_genetic_map)] }
			 .set{ samples_with_map}
}

/* Load maff value for iHS if applied */
if (params.maff) maff = Channel.value(params.maff)

/* Load cutoff value for iHS plotting if applied */
if (params.cutoff) cutoff = Channel.value(params.cutoff)

/* Load mart file for annotation if applied */
if (params.mart) mart_file = Channel.fromPath(params.mart)

/* Load mart file for annotation if applied */
if (params.imart) imart_file = Channel.fromPath(params.imart)

/* Load mart file for annotation if applied */
if (params.imerged) imerged = Channel.value(params.imerged)

/* Load mart file for annotation if applied */
if (params.pmerged) pmerged = Channel.value(params.pmerged)

/* Load rscript for iHS ihs_treatment */
r_script_ihs = Channel.fromPath("nf_modules/core-rscripts/ihs_treatment.R")

/* Load files from PBS design file */
Channel
    .fromPath(params.input_pbs)
    .splitCsv(header:true)
		.map{ row -> [ file(row.path_vcf), file(row.path_pop1), file(row.path_pop2), file(row.path_popout)] }
    .set{ samples_pbs}

/* Load rscript for PBS calculation and treatment */
r_script_pbs = Channel.fromPath("nf_modules/core-rscripts/pbs_calculator.R")
r_script_format_pbs = Channel.fromPath("nf_modules/core-rscripts/pbs_format.R")
r_script_format_ihs = Channel.fromPath("nf_modules/core-rscripts/ihs_format.R")
r_script_merge = Channel.fromPath("nf_modules/core-rscripts/circus.R")

/* Import modules
*/
 include {phasing_with_ref; vcf_to_hap; generating_map;
	 ihs_computing; add_chromosome; merging_chromosomes;
	 fst_calculation; fst_calculation_2; fst_calculation_3;
	 af_1; af_2; af_3; pbs_by_snp; ggf_format; pbs_annotation;
	 ihs_ggf_format; ihs_annotation; merged_results} from './nf_modules/modules.nf'

/*
* main pipeline logic
*/

 workflow  {
	 if (params.notphased) {
		 p1 = phasing_with_ref(samples_ihs)
	 } else {
		 p1 = samples_phased
	 }
	 p2 = vcf_to_hap(p1)
	 if (params.genetic_map) {
		 p3 = samples_with_map
	 } else {
		 p3 = generating_map(samples_genetic_map)
	 }
	 p4 = p2.combine(p3, by: 0)
	 if (params.maff) {
		 p5 = ihs_computing(p4, maff)
	 } else {
		 p5 = ihs_computing(p4, 0.01)
	 }
	 p6 = add_chromosome(p5)
	 p7 = p6.collect()
	 if (params.cutoff) {
		 p8 = merging_chromosomes(p7, r_script_ihs, cutoff)
	 } else {
		 p8 = merging_chromosomes(p7, r_script_ihs, 2)
	 }
	 if (params.imart) {
		 p8_a = ihs_ggf_format(p8.ihs_tsv, imart_file, r_script_format_ihs)
		 p8_b = ihs_annotation(p8_a.ihs_gff, p8_a.biomart_gff)
	 }
	 p9 = fst_calculation(samples_pbs)
	 p10 = fst_calculation_2(samples_pbs)
	 p11 = fst_calculation_3(samples_pbs)
	 p12 = af_1(samples_pbs)
	 p13 = af_2(samples_pbs)
	 p14 = af_3(samples_pbs)
	 p15 = p9.mix(p10,p11,p12,p13,p14).toList()
	 p16 = pbs_by_snp(p15, r_script_pbs)
	 if (params.mart){
		 p17 = ggf_format(p16.png_tsv, mart_file, r_script_format_pbs)
		 p18 = pbs_annotation(p17.pbs_gff, p17.biomart_gff)
	 }
	 if (params.imerged && params.pmerged) {
		 p19 = merged_results(p16.png_tsv, p8.ihs_tsv, pmerged, imerged ,r_script_merge)
	 } else {
		 p19 = merged_results(p16.png_tsv, p8.ihs_tsv, 0.2, 2, r_script_merge)
	 }
 }
