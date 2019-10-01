import glob
import pandas as pd
from snakemake.utils import min_version

### Global parameters and config ###
min_version('5.4')

configfile: "config.yaml"

INDIR 		   = config['inputdir']
OUTDIR  	   = config['outdir']
LOGDIR   	   = config['logdir']
PROJECTS	   = pd.read_csv(config['projects'], sep='\t')['Project'] #TODO: clean this up

# Global scope functions

#TODO: log the key error
def get_resource(rule,resource):
	'''
	Attempt to parse config.yaml to retrieve resources available for each rule.
	It will revert to default is a key error is found.
	'''
	try:
		return config["rules"][rule]["res"][resource]
	except KeyError:
		return config["rules"]["default"]["res"][resource]

#TODO: refactor this
def get_maf_file(wildcards):
	project_dir = f'{INDIR}/MAF/' + wildcards.project + '/*.maf'

	return glob.glob(project_dir)

def get_cnv_file(wildcards):
	project_dir = f'{INDIR}/CNV/' + wildcards.project + '/*focal_score_by_genes.txt'
	return glob.glob(project_dir)

def get_mrna_file(wildcards):
	project_dir = f'{INDIR}/MRNA/' + wildcards.project + '/*htseq_fpkm.tsv'
	return glob.glob(project_dir)

#### RULES ####

rule all:
	input:
		OUTDIR + '/PLOTS/cases_druggable.svg'
		#expand(OUTDIR + '/MRNA/{project}/{project}_expr_filtered.csv', project=PROJECTS)
		
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
	threads: 1
	resources:
		mem=2048
	shell:
		"./scripts/calculate_gscores.py {input} {output}"		

rule filter_maf_files:
	input:
		get_maf_file
	output:
		filtered_maf         = OUTDIR + '/MAF/{project}/{project}_filtered.csv',
		filtered_maf_metrics = OUTDIR + '/MAF/{project}/{project}_metrics.csv'
	threads:
		get_resource('filter_maf_files', 'threads')
	resources:
		mem=get_resource('filter_maf_files', 'mem')
	shell:
		"./scripts/maf_filter.py {input} {output}"

rule filter_cnv_files:
	input:
		get_cnv_file,
		metadata = INDIR + '/METADATA/cnv_metadata.json',
		affy_db  = 'reference/generated/affy_snp_6.0_translation.csv'
	output:
		filtered_cnv 		 = OUTDIR + '/CNV/{project}/{project}_cnv_filtered.csv',
		filtered_cnv_metrics = OUTDIR + '/CNV/{project}/{project}_metrics.csv'
	threads:
		get_resource('filter_cnv_files', 'threads')
	resources:
		mem=get_resource('filter_cnv_files', 'mem')
	shell:
		"./scripts/cnv_filter.py {input} {output}"

rule filter_mrna_files:
	input:
		get_mrna_file,
		mrna_db  = 'reference/generated/RNAseq_transcripts_translation.csv',
		metadata = INDIR + '/METADATA/cnv_metadata.json'
	output:
		filtered_expression = OUTDIR + '/MRNA/{project}/{project}_expr_filtered.csv'
	threads:
		get_resource('filter_mrna_files', 'threads')
	resources:
		mem=get_resource('filter_mrna_files', 'mem')
	shell:
		"./scripts/expression_filter.py {input} {output}"
	
rule generate_cases_table:
	input:
		rules.filter_maf_files.output.filtered_maf,
		rules.filter_cnv_files.output.filtered_cnv,
		'reference/CancerGeneCensus.tsv',
	output:
		cases_table         = OUTDIR + '/MERGED/{project}/cases_table.csv',
		cases_table_metrics = OUTDIR + '/MERGED/{project}/cases_table_metrics.csv'
	threads:
		get_resource('generate_cases_table', 'threads')
	resources:
		mem=get_resource('generate_cases_table','mem')
	shell:
		"./scripts/generate_patients_table.py {input} {output}"

rule check_gain_of_function_events:
	input:
		rules.generate_cases_table.output.cases_table,
		rules.filter_mrna_files.output.filtered_expression,
	output:
		cases_table_corrected = OUTDIR + '/MERGED/{project}/cases_table_corrected.csv'
	threads:1
	resources:
		mem=2048
	shell:
		"./scripts/gain_of_function_correction.py {input} {output}"

rule vulcanspot_annotation:
	input:
		OUTDIR + '/MERGED/{project}/cases_table_corrected.csv',
		'reference/tcga/tcga-vulcan.tsv',
		'reference/generated/vulcan_treatments_db.csv'
	output:
		OUTDIR + '/MERGED/{project}/cases_table_vulcan_annotated.csv'
	threads:
		get_resource('vulcanspot_annotation', 'threads')
	resources:
		mem=get_resource('vulcanspot_annotation', 'mem')
	shell:
		"./scripts/vulcanspot_annotation.py {input} {output}"

#TODO: config for this?
rule generate_summary:
	input:
		expand(OUTDIR + '/MERGED/{project}/cases_table_vulcan_annotated.csv', project=PROJECTS)
	output:
		summary = OUTDIR + '/summary.csv'
	resources:
		mem=2048
	shell:
		"awk 'NR == 1 || FNR > 1' {input} > {output}"

rule generate_summary_plots:
	input:
		rules.generate_summary.output.summary
	output:
		OUTDIR + '/PLOTS/cases_druggable.svg',
		OUTDIR + '/PLOTS/alterations_classified.svg',
		OUTDIR + '/PLOTS/alterations_count_local.svg',
		OUTDIR + '/PLOTS/alterations_count_pancancer.svg'
	
	threads:
		get_resource('generate_summary_plots', 'threads')
	resources:
		mem=get_resource('generate_summary_plots', 'mem')
	shell:
		"./scripts/generate_summary_plots.py {input} {OUTDIR}/PLOTS"