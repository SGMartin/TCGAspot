#!/usr/bin/env python3
"""
Section: SL in TCGA data

Purpose: Generate a table of genetic alterations based on CNV and MAF data
from the TCGA. Patients genetic landscape is inferred from these sources.
"""
import sys

import pandas as pd


def main(input_maf:str, input_cnv:str, cancer_census:str, context_translation:str):

	# Loading input dataframes #
	maf_to_merge = pd.read_csv(input_maf, sep=',', low_memory=False)
	cnv_to_merge = pd.read_csv(input_cnv, sep=',', low_memory=False)

	merged_alterations = merge_cnv_to_mutations(maf_to_merge, cnv_to_merge)

	# Annotate with CGC and translate projects to CCLE tissues

	annotated_tcga			 = annotate_cancer_gene_census(cancer_census,
												   annotated_alterations)

	context_translated_tcga  = translate_tcga_ccle_tissues(context_translation,
														   annotated_tcga
														  )

	# Classify alterations in GoF/LoF/Unknown
	context_translated_tcga['Consequence'] = context_translated_tcga.apply(func=annotate_gof_lof,
																		   axis='columns'
																		   )
	
	# clean inconsistent alterations
	tcga_func_annotated = delete_inconsistent_alterations(context_translated_tcga)

	# SAVING #

def merge_cnv_to_mutations(maf: pd.DataFrame, cnv: pd.DataFrame) -> pd.DataFrame:
	'''
	Merges maf and cnv together by outer join and then fills missing rows
	with None instead of NaN.
	'''

	merged_maf_cnv =  maf.merge(
							right=cnv,
					        how='outer',
						    left_on=['case_id', 'Hugo_Symbol'],
							right_on=['case_id', 'Gene Symbol']
							)
	
	merged_maf_cnv.drop('Gene Symbol', axis=1, inplace=True)
	
	# Fill NaN in Project column for cnv
	cases_projects = maf.loc[:,['case_id', 'Project']]
	cases_projects = cases_projects[~cases_projects.duplicated()] # keeps only first case entry

	fill_dict = cases_projects.set_index('case_id')['Project'].to_dict()

	merged_maf_cnv['Project'] = merged_maf_cnv['case_id'].map(fill_dict)

	# Fill nas in copy_number and variant_classification due to merging.
	merged_maf_cnv['Variant_Classification'].fillna(value='None', inplace=True)
	merged_maf_cnv['copy_number'].fillna(value='None', inplace=True)

	return merged_maf_cnv

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
	'''
	Translates from TCGA project names to the most suitable VulcanSpot context
	based on an external equivalence table built manually.
	'''
	context_translation = pd.read_csv(context_translation, sep='\t')
	context_translation = context_translation.set_index('TCGA')['Vulcanspot'].to_dict()

	tcga_data['Context']      = tcga_data['Project'].map(context_translation)

	return tcga_data


def annotate_gof_lof(tcga_data: pd.Series) -> str:
	'''
	Uses Variant_Classification columns, copy number column when available
	and VAF to predict alteration consequence based on VulcanSpot criteria.
	Returns a dataframe where all alterations are classified either in GoF,
	LoF or Unknown.
	'''

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


def delete_inconsistent_alterations(tcga_data: pd.DataFrame) -> pd.DataFrame:
	'''
	Delete rows where a GoF was predicted on a gene annotated as tumor suppressor
	gene. These although rare (<1%) might happen due to how TCGA GISTIC-2 data
	was generated.
	'''

	is_GoF = tcga_data['Consequence'] == 'GoF'
	is_TSG = tcga_data['Role'].str.contains('TSG')

	tcga_cleared = tcga_data[~(is_GoF & is_TSG)]

	return tcga_cleared



if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
