"""
Section: SL in TCGA data

Purpose: Generate a table of genetic alterations based on CNV and MAF data
from the TCGA. Patients genetic landscape is inferred from these sources.
"""

import pandas as pd


def main():

	# ---------- SNAKEMAKE INPUT ----------- #
	input_maf         = snakemake.input[0]
	input_cnv         = snakemake.input[1]
	input_notation    = snakemake.input[2]

	# ---------- SNAKEMAKE OUTPUT ---------- #
	valid_maf_patients = snakemake.output[0]
	valid_cnv_patients = snakemake.output[1]
	combined_maf_cnv   = snakemake.output[2]
	metrics_report 	   = snakemake.output[3]
    # --------------------------------------- #

	mutation_data   = pd.read_csv(input_maf, sep='\t')
	cnv_data        = pd.read_csv(input_cnv, sep='\t')
	annotations     = pd.read_csv(input_notation, sep='\t')

	invalid_cases   = annotations['entity_id']

	# ------ FILTERING ---------# #TODO: Refactor these
	filtered_maf = filter_maf_patients(mutation_data, invalid_cases)
	filtered_cnv = filter_cnv_patients(cnv_data, invalid_cases)

	filtered_maf.to_csv(valid_maf_patients, sep='\t', index=False)
	filtered_cnv.to_csv(valid_cnv_patients, sep='\t', index=False)

	# --------- MERGING --------#
	merged_maf_cnv = merge_cnv_to_mutations(filtered_maf, cnv_data)
	merged_maf_cnv.to_csv(combined_maf_cnv, sep='\t', index=False)

	# ------- Report metrics --------- #
	raw_maf_patients 				= mutation_data['case_id'].nunique()
	filtered_maf_patients			= filtered_maf['case_id'].nunique()
	raw_cnv_patients				= len(cnv_data.columns)
	filtered_cnv_patients			= len(filtered_cnv.columns)

	raw_cnv_patients -=1 		# accounting for gene column
	filtered_cnv_patients -=1
	
	report = {'Item':['maf_cases', 'maf_cases_filtered', 'cnv_cases', 'cnv_filtered'],
			  'Count':[raw_maf_patients, filtered_maf_patients, raw_cnv_patients, filtered_cnv_patients]
			 }
	
	report = pd.DataFrame(report).to_csv(metrics_report, sep='\t')

def filter_maf_patients(mutations: pd.DataFrame, invalid_cases ) -> pd.DataFrame:
	"""
	Remove patients whose id is present in annotation file.
	"""
	cases_to_remove = mutations['case_id'].isin(invalid_cases)

	# keep cases whose id is not in annotation.txt
	mutations = mutations[~cases_to_remove]

	return mutations

def filter_cnv_patients(cnv: pd.DataFrame, invalid_cases: pd.Series) -> pd.DataFrame:
	'''
	Remove patients whose id is present in annotation
	'''	
	# remember that multiple sample cases are duplicated
	filtered_cnv = cnv[cnv.columns.difference(invalid_cases)]

	return filtered_cnv

# TODO: consider spliting this method in two.
def merge_cnv_to_mutations(maf: pd.DataFrame, cnv: pd.DataFrame) -> pd.DataFrame:

	"""
	Merges cnv and mutations together in a single dataframe. Replaces
	GISTIC-2 gene level copy number scores with cnv_gain for +1 and cnv_loss
	for -1"
	"""
   # melting transforms a wide DF to a LONG format
	filtered_cnv = pd.melt(frame=cnv,
	                       id_vars=['Gene Symbol'],
						   var_name='case_id',
						   value_name='copy_number'
						)
	
	# drop neutrals and rename
	filtered_cnv = filtered_cnv[filtered_cnv.loc[:,'copy_number'] != 0]

	gain_and_loss = {1: 'cnv_gain', -1: 'cnv_loss'}
	filtered_cnv = filtered_cnv.replace({'copy_number' : gain_and_loss})

	filtered_cnv.rename(
						{'Gene Symbol':'Hugo_Symbol'},
	 					axis='columns',
						inplace=True
					   )
					
	# merge both dataframes
	merged_maf_cnv =  maf.merge(
							right=filtered_cnv,
					        how='outer',
						    left_on=['case_id', 'Hugo_Symbol'],
							right_on=['case_id', 'Hugo_Symbol']
							)
	
	# Fill NaN in Project column for cnv
	cases_projects = maf.loc[:,['case_id', 'Project']]
	cases_projects = cases_projects[~cases_projects.duplicated()] # keeps only first case entry

	fill_dict = cases_projects.set_index('case_id')['Project'].to_dict()

	merged_maf_cnv['Project'] = merged_maf_cnv['case_id'].map(fill_dict)

	# Fill nas in copy_number and variant_classification due to merging.
	merged_maf_cnv['Variant_Classification'].fillna(value='None', inplace=True)
	merged_maf_cnv['copy_number'].fillna(value='None', inplace=True)

	return merged_maf_cnv


if __name__ == "__main__":
	main()
