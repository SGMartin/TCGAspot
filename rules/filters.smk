
rule filter_maf_files:
	input:
		get_maf_file
	output:
		filtered_maf         = TABLESDIR + '/MAF/{project}/{project}_filtered.csv',
		filtered_maf_metrics = TABLESDIR + '/MAF/{project}/{project}_metrics.csv'
	threads:
		get_resource('filter_maf_files', 'threads')
	resources:
		mem=get_resource('filter_maf_files', 'mem')
	conda:
		"../envs/tcgaspot.yaml"
	script:
		"../scripts/maf_filter.py"

rule filter_cnv_files:
	input:
		get_cnv_file,
		metadata = INDIR + '/METADATA/cnv_metadata.json',
		affy_db  = 'databases/generated/affy_snp_6.0_translation.csv'
	output:
		filtered_cnv 		 = TABLESDIR + '/CNV/{project}/{project}_cnv_filtered.csv',
		filtered_cnv_metrics = TABLESDIR + '/CNV/{project}/{project}_metrics.csv'
	threads:
		get_resource('filter_cnv_files', 'threads')
	resources:
		mem=get_resource('filter_cnv_files', 'mem')
	conda:
		"../envs/tcgaspot.yaml"
	script:
		"../scripts/cnv_filter.py"

rule filter_mrna_files:
	input:
		get_mrna_file,
		mrna_db  = 'databases/generated/RNAseq_transcripts_translation.csv',
		metadata = INDIR + '/METADATA/cnv_metadata.json'
	output:
		filtered_expression = TABLESDIR + '/MRNA/{project}/{project}_expr_filtered.csv'
	threads:
		get_resource('filter_mrna_files', 'threads')
	resources:
		mem=get_resource('filter_mrna_files', 'mem')
	conda:
		"../envs/tcgaspot.yaml"
	script:
		"../scripts/expression_filter.py"
	