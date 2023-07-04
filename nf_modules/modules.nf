#!/usr/bin/env nextflow

/*Modulos de prueba para  nextflow */
results_dir = "./test/results"
intermediates_dir = "./test/results/intermediates"

process phasing_with_ref {

publishDir "${results_dir}/phasing_with_ref", mode:"copy"

	input:
	tuple val(chromosome), path(path_vcf), path(path_genetic_map), path(path_reference_vcf), path(path_index_reference), path(path_ancestral), path(manifest)

	output:
	tuple val(chromosome), path("${chromosome}.phased.with.ref.vcf"), path(path_genetic_map), path(path_ancestral), path(manifest)

	"""
	tabix -p vcf ${path_vcf}

	awk '{print \$4, \$1, \$3}' ${path_genetic_map} > file_1
	echo -e "pos\tchr\tcM" > file_2
	cat file_2 file_1 > genetic_map

	shapeit4.2 --input ${path_vcf} \
	--map genetic_map \
	--region chr${chromosome} \
	--reference ${path_reference_vcf} \
	--output ${chromosome}.phased.with.ref.vcf
	"""
}

process ancestral_annotation{

publishDir "${results_dir}/annotated_vcf/", mode:"copy"

	input:
	file java_script
	tuple val(chromosome), path(path_files), path(path_genetic_map), path(path_ancestral), path(manifest)

	output:
	tuple val(chromosome), path("chr${chromosome}_aa.vcf")

	"""
	samtools faidx ${path_ancestral}

	java -jar ${java_script} \
	-m ${manifest} \
	${path_files} |\
	bcftools annotate -x '^INFO/AA' > chr${chromosome}_aa.vcf
	"""
}


process ancestral_vcf{

publishDir "${results_dir}/ancestral_annotated_vcf/", mode:"copy"

	input:
	file java_annotation_script
	file annotation_script
	tuple val(chromosome), path(path_ancestral_files)

	output:
	tuple val(chromosome), path("final_annotated_${chromosome}.vcf")

	"""
	java -jar ${java_annotation_script} vcffilterjdk -f ${annotation_script} ${path_ancestral_files} > final_annotated_${chromosome}.vcf
	"""
}

process ihs_rehh {

publishDir "${results_dir}/ihs_compute/", mode:"copy"

	input:
	file r_script_ihs
	tuple val(chromosome), path(path_ancestral_files)

	output:
	tuple path("ihs_${chromosome}.tsv")

	"""
	Rscript --vanilla ${r_script_ihs} ${path_ancestral_files} ihs_${chromosome}.tsv
	"""
}

process ihs_images {

publishDir "${results_dir}/all_chr_ihs/", mode:"copy"

	input:
	path(path_files)
	file rscript

	output:
	path "*.png", emit: ihs_plots
	path "final_ihs_onepercent.tsv", emit: ihs_tsv_one_percent
	path "final_ihs.tsv", emit: ihs_tsv

	"""
	Rscript --vanilla ihs_final.R . final_ihs.tsv final_ihs_onepercent.tsv final_ihs.histogram.png
	"""
}

process fst_calculation {

publishDir "${results_dir}/fst_results_pop1_pop2/", mode:"copy"

	input:
	tuple path(path_vcf), path(path_pop1), path(path_pop2), path(path_popout)

	output:
	path "*.fst"

	"""
	vcftools --vcf ${path_vcf} \
					 --weir-fst-pop ${path_pop1} \
					 --weir-fst-pop ${path_pop2} \
					 --out pop1pop2
	"""
}

process fst_calculation_2 {

publishDir "${results_dir}/fst_results_pop1_popout/", mode:"copy"

	input:
	tuple path(path_vcf), path(path_pop1), path(path_pop2), path(path_popout)

	output:
	path "*.fst"

	"""
	vcftools --vcf ${path_vcf} \
					 --weir-fst-pop ${path_pop1} \
					 --weir-fst-pop ${path_popout} \
					 --out pop1popout
	"""
}

process fst_calculation_3 {

publishDir "${results_dir}/fst_results_pop2_popout/", mode:"copy"

	input:
	tuple path(path_vcf), path(path_pop1), path(path_pop2), path(path_popout)

	output:
	path "*.fst"

	"""
	vcftools --vcf ${path_vcf} \
					 --weir-fst-pop ${path_pop2} \
					 --weir-fst-pop ${path_popout} \
					 --out pop2popout
	"""
}

process af_1 {

	publishDir "${results_dir}/af_pop1/",mode:"copy"

	input:
	tuple path(path_vcf), path(path_pop1), path(path_pop2), path(path_popout)

	output:
	file "*.frq"

	"""
	vcftools --vcf ${path_vcf} --keep ${path_pop1} --freq --out pop1
	"""
}

process af_2 {

	publishDir "${results_dir}/af_pop2/",mode:"copy"

	input:
	tuple path(path_vcf), path(path_pop1), path(path_pop2), path(path_popout)

	output:
	file "*.frq"

	"""
	vcftools --vcf ${path_vcf} --keep ${path_pop2} --freq --out pop2
	"""
}

process af_3 {

	publishDir "${results_dir}/af_pop3/",mode:"copy"

	input:
	tuple path(path_vcf), path(path_pop1), path(path_pop2), path(path_popout)

	output:
	file "*.frq"

	"""
	vcftools --vcf ${path_vcf} --keep ${path_popout} --freq --out pop3
	"""
}

process pbs_by_snp {

	publishDir "${results_dir}/pbs_by_snp/",mode:"copy"

	input:
	file p15
	file r_script_pbs

	output:
	path "*.png", emit: png_pbs
	path "one_percent_pbs.tsv", emit: one_percent_tsv
	path "pbs.tsv", emit: pbs_tsv

	"""
	Rscript --vanilla pbs_calculator_final.R .
	"""
}

process merged_results {

	publishDir "${results_dir}/pbs_vs_ihs/",mode:"copy"

	input:
	file p16
	file p8
	file r_script_merged

	output:
	path "circus.png", emit: png_file
	path "pbs_vs_ihs.tsv", emit: tsv_file
	path "final_bed", emit: bed_file

	"""
	Rscript --vanilla merging_pbs_ihs.R ${p16} ${p8}
	"""
}

process filter_vcf {

	publishDir "${results_dir}/filtered_vcf/",mode:"copy"

	input:
	tuple path(path_vcf), path(path_pop1), path(path_pop2), path(path_popout)
	file bed_file

	output:
	path "final_results.recode.vcf", emit: vcf_file

	"""
	vcftools --vcf ${path_vcf} --bed ${bed_file} --out final_results --recode
	"""
}

process annotation {

	publishDir "${results_dir}/annotation/",mode:"copy"

	input:
	path annovar
	file vcf

	output:
	path "*annotation*", emit: vcf_file

	"""
	mv annovar/* .
	perl table_annovar.pl -vcfinput ${vcf} . -buildver hg38 -out annotation -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring . -polish -xref example/gene_xref.txt
	"""
}