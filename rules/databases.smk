rule rebuild_vulcan_database:
	input:
	output:
		"reference/generated/vulcan_gene_db.csv",
		"reference/generated/vulcan_treatments_db.csv"
	threads: 
		get_resource('rebuild_vulcan_database', 'threads')
	resources:
		mem=get_resource('rebuild_vulcan_database', 'mem')
	script:
		"./scripts/get_vulcan_database.py"

rule rebuild_snp_array_dictionary:
	input:
		"reference/tcga/affy_SNP6.0_ensg.tsv"
	output:
		"reference/generated/affy_snp_6.0_translation.csv"
	threads:
		get_resource('rebuild_snp_array_translation', 'threads')
	resources:
		mem=get_resource('rebuild_snp_array_translation', 'mem')
	script:
		"./scripts/get_affy_translation.py"

rule rebuild_rnaseq_dictionary:
	input:
		"reference/tcga/RNAseq_transcripts.csv"
	output:
		"reference/generated/RNAseq_transcripts_translation.csv"
	threads:
		get_resource('rebuild_rnaseq_dictionary', 'threads')
	resources:
		mem=get_resource('rebuild_rnaseq_dictionary', 'mem')
	script:
		"./scripts/get_rnaseq_translation.py"

rule rebuild_gscore_database:
	input:
		"reference/gscore/hugo_all_genes.tsv",
		"reference/gscore/TumorPortal.csv",
		"reference/CancerGeneCensus.tsv",
		"reference/gscore/srep02650-s3.csv",
		"reference/gscore/gene_essentiality_score.tsv",
		"reference/gscore/oncoscape_all_matrix_highscore.tsv"
	output:
		"reference/generated/genes_gscore.csv"
	threads:
		get_resource('rebuild_gscore_database', 'threads')
	resources:
		mem=get_resource('rebuild_gscore_database', 'mem')
	shell:
		"./scripts/calculate_gscores.py {input} {output}"		
