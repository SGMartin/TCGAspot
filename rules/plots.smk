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
		mem=get_resource('generate_summary_plots', 'mem')
	shell:
		"./scripts/generate_summary_plots.py {input} {PLOTDIR}/summary"

rule generate_vulcanspot_plots:
	input:
		TABLESDIR + '/summary.csv',
		'databases/generated/vulcan_treatments_db.csv'
	output:
		PLOTDIR + '/vulcanspot/cases_druggable.svg',
		PLOTDIR + '/vulcanspot/alterations_count_local.svg',
		PLOTDIR + '/vulcanspot/alterations_count_pancancer.svg',
		PLOTDIR + '/vulcanspot/drug_sources.svg'
	threads:1
	shell:
		"./scripts/generate_plots_vulcanspot.py {input} {PLOTDIR}/vulcanspot"

rule generate_pandrugs_plots:
	input:
		TABLESDIR + '/summary.csv',
		'databases/pandrugs/Pandrugs_Feb2018.tsv'
	output:
		PLOTDIR + '/pandrugs/cases_druggable.svg'
	threads:1
	shell:
		"./scripts/pandrugs_annotation_plots.py {input} {PLOTDIR}/pandrugs"