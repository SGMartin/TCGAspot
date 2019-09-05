configfile: "config.yaml"



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

rule rebuild_snp_array_translation:
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



'''
rule filter_maf:
	input: 
		raw_maf  = "../../RAW/PAN-TCGA/MAF/COAD/{maf}.maf"
	output:
		filtered = "../../GENERATED/TCGA/MAF/COAD/{maf}_filtered.csv",
		metrics  = "../../GENERATED/TCGA/MAF/COAD/{maf}_filtered_metrics.csv"
	threads: 1
	shell:
		"./scripts/maf_filter.py {input.raw_maf} {output.filtered} {output.metrics}"

rule filter_cnv:
	input:
		cnv 	 = "../../RAW/PAN-TCGA/CNV/COAD/{cnv}.txt",
		metadata = "../../RAW/PAN-TCGA/METADATA/cnv_metadata.json",
		affy_db  = "../../GENERATED/TCGA/affy_snp_6.0_translation.csv"
	
	output:
		filtered_cnv = "../../GENERATED/TCGA/CNV/COAD/{cnv}_filtered.csv",
		metrics 	 = "../../GENERATED/TCGA/CNV/COAD/{cnv}_metrics.csv"
	
	threads:1
	shell:
		"./scripts/cnv_filter.py {input.cnv} {input.metadata} {input.affy_db} {output.filtered_cnv} {output.metrics}"	

rule filter_annotated_cases:
	input:
		maf   = "../../GENERATED/TCGA/MAF/COAD/{maf}_filtered.csv",
		annot = "../../RAW/PAN-TCGA/MAF/COAD/annotations.txt"
	output:
		cases   = "../../GENERATED/TCGA/MAF/COAD/{maf}_filtered_cases.csv",
		metrics = "../../GENERATED/TCGA/MAF/COAD/{maf}_filtered_metrics.csv"
	threads:1
	shell:
		"./scripts/filter_excluded_cases.py {input.maf} {input.annot} {output.cases} {output.metrics}"

rule filter_annotated_cases_2:
	input:
		maf   = "../../GENERATED/TCGA/CNV/COAD/{maf}_filtered.csv",
		annot = "../../RAW/PAN-TCGA/CNV/COAD/annotations.txt"
	output:
		cases   = "../../GENERATED/TCGA/CNV/COAD/{maf}_filtered_cases.csv",
		metrics = "../../GENERATED/TCGA/CNV/COAD/{maf}_filtered_metrics.csv"
	threads:1
	shell:
		"./scripts/filter_excluded_cases.py {input.maf} {input.annot} {output.cases} {output.metrics}"
'''