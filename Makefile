test.run :
	bsub -R"select[mem>10000] rusage[mem=10000]" -M10000 -o farm_test.out -e farm_test.err perl Control_Script.pl -i Test_inputGenes_file -x Test_expression_metaFile -p "8TFs_SpreadC" -r 2 -t 1 -w -R /software/R-3.2.2/bin/R 
