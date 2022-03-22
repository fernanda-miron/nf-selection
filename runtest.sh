input_ihs="test/data/ihs_files/design_file.csv"
input_pbs="test/data/pbs_files/design_file.csv"
output_directory="test/results"

echo -e "======\n Testing NF execution \n======" \
&& rm -rf $output_directory \
&& nextflow run ihs.nf \
	--input_ihs $input_ihs \
	--input_pbs $input_pbs \
	--output_dir $output_directory \
	--cutoff 2 \
	--maff 0.01 \
	--mart "test/data/pbs_files/mart_export.txt" \
	--imart "test/data/ihs_files/mart_export.txt" \
	--genetic_map \
	-resume \
	-with-report $output_directory/`date +%Y%m%d_%H%M%S`_report.html \
	-with-dag $output_directory/`date +%Y%m%d_%H%M%S`.DAG.html \
&& echo -e "======\n Basic pipeline TEST SUCCESSFUL \n======"
