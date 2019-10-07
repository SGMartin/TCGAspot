rule generate_summary_plots:
	input:
		TABLESDIR + '/summary.csv'
	output:
		PLOTDIR + '/summary/alterations_classified_0.svg',
		PLOTDIR + '/summary/alterations_classified_0.2.svg',	
		PLOTDIR + '/summary/alterations_classified_0.4.svg',	
		PLOTDIR + '/summary/alterations_classified_0.6.svg',	
		PLOTDIR + '/summary/alterations_classified_0.8.svg'		
	threads:
		get_resource('generate_summary_plots', 'threads')
	resources:
		mem_mb=get_resource('generate_summary_plots', 'mem_mb')
	conda:
		"../envs/tcgaspot.yaml"
	script:
		"../scripts/generate_summary_plots.py"


rule generate_vulcanspot_plots:
	input:
		TABLESDIR + '/summary.csv',
		'databases/generated/vulcan_treatments_db.csv'
	output:
		PLOTDIR + '/vulcanspot/cases_druggable.svg',
		PLOTDIR + '/vulcanspot/alterations_count_local.svg',
		PLOTDIR + '/vulcanspot/alterations_count_pancancer.svg',
		PLOTDIR + '/vulcanspot/drug_sources.svg'
	threads:
		get_resource('generate_vulcanspot_plots', 'threads')
	resources:
		mem_mb=get_resource('generate_vulcanspot_plots', 'mem_mb')
	conda:
		"../envs/tcgaspot.yaml"
	script:
		"../scripts/generate_plots_vulcanspot.py"

rule generate_pandrugs_plots:
	input:
		TABLESDIR + '/summary.csv',
		'databases/pandrugs/Pandrugs_Feb2018.tsv'
	output:
		PLOTDIR + '/pandrugs/cases_druggable.svg'
	threads:
		get_resource('generate_pandrugs_plots', 'threads')
	resources:
		mem_mb=get_resource('generate_pandrugs_plots', 'mem_mb')
	conda:
		"../envs/tcgaspot.yaml"
	script:
		"../scripts/pandrugs_annotation_plots.py"
