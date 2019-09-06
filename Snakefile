import glob
import pandas as pd
from snakemake.utils import min_version

### Global parameters and config ###
min_version('5.4')

configfile: "config.yaml"

INDIR 		   = config['inputdir']
OUTDIR  	   = config['outdir']
LOGDIR   	   = config['logdir']

PROJECT,MAF,  = glob_wildcards(f'{INDIR}/MAF/' + '{project}/{maf}.maf')
CPROJECT,CNV, = glob_wildcards(f'{INDIR}/CNV/' + '{cproject}/{cnv}.focal_score_by_genes.txt')

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

#### RULES ####

rule all:
	input:
		expand(OUTDIR + '/MERGED/{project}/cases_table_vulcan_annotated.csv', project=PROJECT)
		

rule rebuild_vulcan_database:
	input:
	output:
		"reference/generated/vulcan_db.csv"
	threads: 
		get_resource('rebuild_vulcan_database', 'threads')
	resources:
		mem=get_resource('rebuild_vulcan_database', 'mem')
	script:
		"./scripts/get_vulcan_database.py"

rule rebuild_snp_array_dictionary:
	input:
		"reference/affy_SNP6.0_ensg.tsv"
	output:
		"reference/generated/affy_snp_6.0_translation.csv"
	threads:
		get_resource('rebuild_snp_array_translation', 'threads')
	resources:
		mem=get_resource('rebuild_snp_array_translation', 'mem')
	script:
		"./scripts/get_affy_translation.py"

rule filter_maf_files:
	input:
		get_maf_file
	output:
		filtered_maf = OUTDIR + '/MAF/{project}/{project}_filtered.csv',
		filtered_cnv = OUTDIR + '/MAF/{project}/{project}_metrics.csv'
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

rule generate_cases_table:
	input:
		rules.filter_maf_files.output.filtered_maf,
		rules.filter_cnv_files.output.filtered_cnv,
		'reference/CancerGeneCensus.tsv',
		'reference/tcga-vulcan.tsv'
	output:
		OUTDIR + '/MERGED/{project}/cases_table.csv',
		OUTDIR + '/MERGED/{project}/cases_table_metrics.csv'
	threads:
		get_resource('generate_cases_table', 'threads')
	resources:
		mem=get_resource('generate_cases_table','mem')
	shell:
		"./scripts/generate_patients_table.py {input} {output}"


rule vulcanspot_annotation:
	input:
		OUTDIR + '/MERGED/{project}/cases_table.csv',
		'reference/generated/vulcan_db.csv'
	output:
		OUTDIR + '/MERGED/{project}/cases_table_vulcan_annotated.csv'
	
	threads:
		get_resource('vulcanspot_annotation', 'threads')
	resources:
		mem=get_resource('vulcanspot_annotation', 'mem')
	shell:
		"./scripts/vulcanspot_annotation.py {input} {output}"