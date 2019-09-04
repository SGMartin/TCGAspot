
rule get_vulcan_database:
	input:
	output:
		"../../GENERATED/TCGA/vulcan_db.tsv"
	threads: 1
	script:
		"./scripts/get_vulcan_database.py"

rule get_snp_array_translation:
	input:
		"../../RAW/Databases/affy_SNP6.0_ensg.tsv"
	output:
		"../../GENERATED/TCGA/affy_snp_6.0_translation.csv"
	threads:1
	script:
		"./scripts/get_affy_translation.py"

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
