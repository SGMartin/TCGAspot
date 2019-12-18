#!/usr/bin/env python3
'''
Synthetic lethality validation: This module attempts to validate the presence
in the TCGA of inferred synthethic lethal gene pairs predicted by VulcanSpot.

To do so it tests whether the p. of both components of the pair being a LoF
is significantly lower than expected if they were not associated.
'''

#TODO: Note to Santi: keep in mind that at the moment only druggable pairs are 
# present in the database. Roughly 13k pairs.

from scipy 		 import stats
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

import pandas as pd
import numpy as np



def main():

	## SNAKEMAKE I/O ##
	vulcan_db 	 = snakemake.input[0]
	tcga_summary = snakemake.input[1]

	# Load datasets
	vulcan_db = pd.read_csv(vulcan_db,
							sep=',',
							usecols=['geneA', 'tissue', 'alteration', 'geneB']
							)
	
	# Drop duplicates due to removing drug info.
	vulcan_db = vulcan_db.drop_duplicates(keep='first')
	
	tcga_summary = pd.read_csv(tcga_summary,
							   sep=',',
							   usecols=['Hugo_Symbol', 'case_id',
							   			'Context', 'Consequence'
									   ]
							   )
	
	# Do not validate pairs with genA = genB as it makes no sense
	candidate_pairs = vulcan_db[vulcan_db['geneA'] != vulcan_db['geneB']]

	candidate_pairs['pvalue'] = candidate_pairs.apply(func=validate_synthetic_pair,
													  axis='columns',
													  args=(tcga_summary,)
													  )

	# Attempt to correct the p-value by FDR: Benjamini and Hochberg
	candidate_pairs['pvalue-c'] = multipletests(candidate_pairs['pvalue'], alpha=0.25, method='fdr_bh')[1]
	candidate_pairs.to_csv(snakemake.output[0], sep=',', index=False)
	

def validate_synthetic_pair(synthetic_pair: pd.Series, tcga):
	'''
	This method checks if a given synthethic lethal pair is less frequent than
	expected in validation_source dataframe, thus hinting the pair is indeed
	lethal. 

	Returns an UNCORRECTED-FDR p-value of the test
	'''
	geneA = synthetic_pair['geneA']
	geneB = synthetic_pair['geneB']

	alteration = synthetic_pair['alteration']
	context    = synthetic_pair['tissue']

	print('Validating pair', geneA, alteration, geneB, 'in', context)

	if context == 'PANCANCER':
		patients = tcga # ignore context for these cases
		
	else:
		# select matched context patients
		patients = tcga[tcga['Context'] == context]
		
	# patients with gen A mutated with the same functional impact
	has_geneA = (patients['Hugo_Symbol'] == geneA) & (patients['Consequence'] == alteration)
	
	# patients with mutated gene B
	has_geneB = (patients['Hugo_Symbol'] == geneB) & (patients['Consequence'] == 'LoF')


	# get counts for each case: 1,2,3,4 in this order

	genA_genB_LoF = patients[has_geneA & has_geneB]['case_id'].nunique()
	genA_gen_B_WT = patients[~has_geneA & ~has_geneB]['case_id'].nunique()

	genA_LoF_genB_WT = patients[has_geneA & ~has_geneB]['case_id'].nunique()
	genA_WT_genB_LoF = patients[~has_geneA & has_geneB]['case_id'].nunique()


	'''						  gene B LoF | gene B other alts + WT
	--------------------------------------------------
	gene A Lof 			   |	1		 | 3
	--------------------------------------------------
	gene A WT + other alts |	4	 	 | 2
	--------------------------------------------------
	'''

	gene_table = [[genA_genB_LoF, genA_LoF_genB_WT],
				  [genA_WT_genB_LoF, genA_gen_B_WT]]

	oddsratio, pvalue = fisher_exact(table=gene_table, alternative='less')

	return pvalue

if __name__ == "__main__":
	main()

