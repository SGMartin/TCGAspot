"""
Section: SL in TCGA data

Purpose: Annotate merged CNV and maf with data from Cancer Gene Census
and SL identified by VulcanSpot.
"""

import numpy as np
from   pandarallel import pandarallel
import pandas as pd

def main():

	# --------------- SNAKE INPUT -------------- #
	cancer_census 		 = snakemake.input[0]
	context_translation  = snakemake.input[1]
	tcga_data            = snakemake.input[2]
	vulcan_data          = snakemake.input[3]

	# --------------  SNAKE OUTPUT ------------- # 
	cores_to_use  = snakemake.threads 
	where_to_save = snakemake.output[0]

	# --------------  METHOD CALLS ------------- #

	tcga_gene_census   		= annotate_cancer_gene_census(cancer_census, tcga_data)
	tcga_translated_tissues = translate_tcga_ccle_tissues(context_translation,
														  tcga_gene_census)
	# Parallel apply operation
	pandarallel.initialize(nb_workers=cores_to_use, shm_size_mb=4096)
	
	tcga_translated_tissues['Consequence'] = tcga_translated_tissues.parallel_apply(func=annotate_gof_lof,axis=1)
	
	tcga_cleared = clear_inconsistent_alterations(tcga_translated_tissues)

	tcga_vulcan_clean = annotate_vulcan_data(vulcan_data, tcga_cleared)
	tcga_vulcan_clean.to_csv(where_to_save, sep='\t', index=False)

def annotate_cancer_gene_census(cancer_census, tcga_data) -> pd.DataFrame:
	"""
	Annotate roles based on cancer gene census.
	"""
	cancer_gene_census = pd.read_csv(cancer_census,
									 sep='\t',
									 usecols=['Gene Symbol', 'Role in Cancer']
								    )

	annotated_mutations = pd.read_csv(tcga_data, sep='\t', low_memory=False)
	
	cancer_gene_census['Role in Cancer'].fillna(value='None', inplace=True)
	cancer_gene_census = cancer_gene_census.set_index('Gene Symbol')['Role in Cancer'].to_dict()

	annotated_mutations['Role'] = annotated_mutations['Hugo_Symbol'].map(cancer_gene_census)
	annotated_mutations['Role'].fillna(value='None', inplace=True)
	
	del cancer_gene_census
	return annotated_mutations


def translate_tcga_ccle_tissues(context_translation: str,
								tcga_data: pd.DataFrame) -> pd.DataFrame:

	# dictionary to translate from TCGA project name to V.spot context

	context_translation = pd.read_csv(context_translation, sep='\t')
	context_translation = context_translation.set_index('TCGA')['Vulcanspot'].to_dict()

	tcga_data['Context']      = tcga_data['Project'].map(context_translation)

	return tcga_data


def annotate_gof_lof(tcga_data: pd.Series) -> str:

	truncated_protein = ['De_novo_Start_OutOfFrame', 'Frame_Shift_Del',
						 'Frame_Shift_Ins', 'In_Frame_Del','In_Frame_Ins',
						 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site',
						 'Start_Codon_Del', 'Start_Codon_Ins',
						 'Stop_Codon_Del', 'Stop_Codon_Ins'
						]
	
	is_truncated = tcga_data['Variant_Classification'] in truncated_protein
	is_missense  = tcga_data['Variant_Classification'] == 'Missense_Mutation'
	has_cnv 	 = tcga_data['copy_number'] != 'None'
	is_oncogene  = 'oncogene' in tcga_data['Role']

	result = 'Unknown'

	if has_cnv:
		if tcga_data['copy_number'] == 'cnv_loss':
			result = 'LoF'

		else: # might be a GoF, check for somatic point mutations
			if (is_truncated == False) & (is_missense == False):
				result = 'GoF'
			else:
				result = 'LoF'
	else:
		if is_truncated  & (tcga_data['VAF'] >= 0.7):
			result = 'LoF'
		
		if is_missense:
			if is_oncogene & (tcga_data['VAF'] >= 0.2):
				result = 'GoF'
			else:
				if tcga_data['VAF'] >= 0.7:
					result = 'LoF'

	return result


def clear_inconsistent_alterations(tcga_data: pd.DataFrame) -> pd.DataFrame:

	is_GoF = tcga_data['Consequence'] == 'GoF'
	is_TSG = tcga_data['Role'].str.contains('TSG')

	tcga_cleared = tcga_data[~(is_GoF & is_TSG)]

	return tcga_cleared


def annotate_vulcan_data(vulcan_table: str,
						 tcga_data: pd.DataFrame) -> pd.DataFrame:

	# Load available genes in vulcanSpot
	vulcan_genes = pd.read_csv(vulcan_table,
							   sep='\t',
							   usecols=['geneA', 'tissue', 'alteration'])

	vulcan_genes.drop_duplicates(keep='first', inplace=True)
	
	# Annotate local vulcan alterations #
	tcga_data = pd.merge(tcga_data, vulcan_genes,
						  how='left',
						  left_on=['Hugo_Symbol', 'Context', 'Consequence'], 
						  right_on=['geneA', 'tissue', 'alteration'], 
						  indicator=True
						)

	tcga_data.drop(['tissue', 'geneA', 'alteration'], inplace=True, axis=1)

	both_to_bool = { 'both': True, 'left_only' : False}	
	tcga_data['Vulcan_Local'] = tcga_data['_merge'].map(both_to_bool)
	tcga_data.drop('_merge', inplace=True, axis=1)
	
	# annotate global vulcan alterations
	#TODO: Refactor this
	vulcan_pancancer_genes = vulcan_genes[vulcan_genes['tissue'] == 'PANCANCER']

	tcga_data = pd.merge(tcga_data, vulcan_pancancer_genes,
						  how='left',
						  left_on=['Hugo_Symbol', 'Consequence'], 
						  right_on=['geneA', 'alteration'], 
						  indicator=True
						)

	tcga_data.drop(['tissue', 'geneA', 'alteration'], inplace=True, axis=1)

	both_to_bool = { 'both': True, 'left_only' : False}	
	tcga_data['Vulcan_Pancancer'] = tcga_data['_merge'].map(both_to_bool)
	tcga_data.drop('_merge', inplace=True, axis=1)

	del vulcan_genes

	return tcga_data


if __name__ == "__main__":
	main()
