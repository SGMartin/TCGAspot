import glob
import pandas as pd
from snakemake.utils import min_version

##### Global parameters and config #####
min_version('5.4')

configfile: "config.yaml"

INDIR 		   = config['inputdir']
TABLESDIR      = config['tablesdir']
LOGDIR   	   = config['logdir']
PLOTDIR        = config['plotdir']
PROJECTS	   = pd.read_csv(config['projects'], sep='\t')['Project'] #TODO: clean this up


##### Global scope functions #####

#TODO: log the key error
def get_resource(rule,resource):
	'''
	Attempt to parse config.yaml to retrieve resources available for each rule.
	It will revert to default is a key error is found.
	'''
	try:
		return config["rules"][rule]["res"][resource]
	except KeyError:
		print(f'Failed to resolve resource config for {rule}/{resource}: using default parameters')
		return config["rules"]["default"]["res"][resource]

#TODO: refactor this
def get_maf_file(wildcards):
	project_dir = f'{INDIR}/MAF/' + wildcards.project + '/*.maf'

	return glob.glob(project_dir)

def get_cnv_file(wildcards):
	project_dir = f'{INDIR}/CNV/' + wildcards.project + '/*focal_score_by_genes.txt'
	return glob.glob(project_dir)

def get_mrna_file(wildcards):
	project_dir = f'{INDIR}/MRNA/' + wildcards.project + '/*htseq_counts.tsv'
	return glob.glob(project_dir)


##### RULES #####

## load rules ##
include: 'rules/databases.smk'
include: 'rules/filters.smk'
include: 'rules/plots.smk'

rule all:
	input:
		#expand(TABLESDIR + '/MERGED/{project}/cases_table_corrected.csv', project=PROJECTS)
		TABLESDIR + '/synthetic_lethal.csv',
		PLOTDIR + '/summary/alterations_classified_0.8.svg',
		PLOTDIR + '/vulcanspot/cases_druggable.svg',
		PLOTDIR + '/pandrugs/cases_druggable.svg',
		TABLESDIR + '/summary.csv'

rule generate_cases_table:
	input:
		rules.filter_maf_files.output.filtered_maf,
		rules.filter_cnv_files.output.filtered_cnv,
		'databases/CancerGeneCensus.tsv',
		'databases/generated/genes_gscore.csv'
	output:
		cases_table         = TABLESDIR + '/MERGED/{project}/cases_table.csv',
		cases_table_metrics = TABLESDIR + '/MERGED/{project}/cases_table_metrics.csv'
	threads:
		get_resource('generate_cases_table', 'threads')
	resources:
		mem=get_resource('generate_cases_table','mem_mb')
	conda:
		"./envs/tcgaspot.yaml"
	script:
		"./scripts/generate_patients_table.py"

rule perform_cnv_diff_expr:
	input:
		cases_table      = rules.generate_cases_table.output.cases_table,
		expression_table = rules.filter_mrna_files.output.filtered_expression
	output:
		tested_cnv = TABLESDIR + '/MERGED/{project}/cnv_diff_expr.csv'
	threads:
		get_resource('perform_cnv_diff_expr', 'threads')
	resources:
		mem_mb=get_resource('perform_cnv_diff_expr', 'mem_mb')
	conda:
		"./envs/deseq2.yaml"
	script:
		"./scripts/cnv_diff_expr.R"

rule check_gain_of_function_events:
	input:
		cases_table      = rules.generate_cases_table.output.cases_table,
		tested_cnv	     = rules.perform_cnv_diff_expr.output.tested_cnv
	output:
		cases_corrected = TABLESDIR + '/MERGED/{project}/cases_table_corrected.csv'
	threads:
		get_resource('check_gain_of_function_events', 'threads')
	resources:
		mem_mb=get_resource('check_gain_of_function_events', 'mem_mb')
	conda:
		"./envs/tcgaspot.yaml"
	script:
		"./scripts/gain_of_function_correction.py"

rule vulcanspot_annotation:
	input:
		TABLESDIR + '/MERGED/{project}/cases_table_corrected.csv',
		'databases/tcga/tcga-vulcan.tsv',
		'databases/generated/vulcan_treatments_db.csv',
		'databases/generated/genes_gscore.csv'
	output:
		TABLESDIR + '/MERGED/{project}/cases_table_vulcan_annotated.csv'
	threads:
		get_resource('vulcanspot_annotation', 'threads')
	resources:
		mem=get_resource('vulcanspot_annotation', 'mem_mb')
	conda:
		"./envs/tcgaspot.yaml"
	script:
		"./scripts/vulcanspot_annotation.py"

rule generate_summary:
	input:
		expand(TABLESDIR + '/MERGED/{project}/cases_table_vulcan_annotated.csv', project=PROJECTS)
	output:
		summary = TABLESDIR + '/summary.csv'
	threads:
		get_resource('generate_summary', 'threads')
	resources:
		mem_mb=get_resource('generate_summary', 'mem_mb')
	shell:
		"awk 'NR == 1 || FNR > 1' {input} > {output}"


rule test_synthetic_pairs:
	input:
		'databases/generated/vulcan_treatments_db.csv',
		rules.generate_summary.output.summary
	output:
		TABLESDIR + '/synthetic_lethal.csv'
	threads: 20
	script:
		"./scripts/synthetic_letal_validation.py"