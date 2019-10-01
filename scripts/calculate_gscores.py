#!/usr/bin/env python3
'''
Gscore retrieval module: This scripts generates a new master table to assess
a gene Gscore. It is included for completeness purposes, but it should not be
run on an analysis basis, but to update the Gscore (for instance, if the sources
of information change).
'''

import sys

import pandas as pd
import numpy  as np

def main(gene_list, tumorportal, cancercensus,
		 drivers, essenciality, oncoscape, where_to_save):
	
	# get a list of genes
	all_genes = get_gene_list(gene_list)

	# get all scores 
	tp_gscore  = get_scores_from_tumorportal(tumorportal)
	gcg_gscore = get_scores_from_gene_census(cancercensus)
	dr_gscore  = get_scores_from_predicted_drivers(drivers)
	ess_gscore = get_scores_from_essenciality(essenciality)
	onc_gscore = get_scores_from_oncoscape(oncoscape)

	all_genes = pd.DataFrame(data=all_genes, columns=['Hugo_Symbol'])

	all_genes['tumorportal'] = all_genes['Hugo_Symbol'].map(tp_gscore)
	all_genes['cancercensus'] = all_genes['Hugo_Symbol'].map(gcg_gscore)
	all_genes['driver'] = all_genes['Hugo_Symbol'].map(dr_gscore)
	all_genes['essenciality'] = all_genes['Hugo_Symbol'].map(ess_gscore)
	all_genes['oncoscape'] = all_genes['Hugo_Symbol'].map(onc_gscore)

	all_genes.fillna(0, inplace=True)
	
	# Aggregate all rows but the first
	all_genes['Gscore'] = all_genes.iloc[:,1:].aggregate(func=np.sum, axis='columns')
	all_genes['Gscore'] = all_genes['Gscore'].round(4)
	
	all_genes.to_csv(where_to_save, sep=',', index=False)

def get_gene_list(gene_db: str) -> list:
	'''
	This method loads the database of known gene names and returns a list
	of valid, unique entries.
	'''

	genes = pd.read_csv(gene_db, sep='\t', usecols=['Approved symbol', 'Status'])
	genes = genes[genes['Status'] == 'Approved']

	gene_list = list(genes['Approved symbol'].unique())	

	del genes 
	return gene_list


def get_scores_from_tumorportal(tumor_db: str) -> dict:
	'''
	Retrieves pairs of gene-scores based on TumorPortal data. Returns a dict.
	of key-value pairs where key=gene name, value= gscore based on TumorPortal.
	'''

	tumor_portal = pd.read_csv(filepath_or_buffer=tumor_db,
							   sep='\t',
							   index_col='Gene',
							   usecols=['Gene','Significance']
							   )

	significance_to_score = {'Near significance': 0.025,
							 'Significantly mutated': 0.05,
							 'Highly significantly mutated': 0.1
							}
	
	tumor_portal['Score'] = tumor_portal['Significance'].map(significance_to_score)

	scores = tumor_portal['Score'].to_dict()

	del tumor_portal
	return scores

def get_scores_from_gene_census(cancercensus: str) -> dict:
	'''
	Retrieves pairs of gene-scores based on Cancer Gene Census data. Returns a 
	dict. of key-value pairs where key=gene name, value= gscore based on GCG.
	'''
	cancer_gene_census = pd.read_csv(cancercensus, sep='\t', usecols=['Gene Symbol'])
	cancer_gene_census['Score'] = 0.1
	cancer_gene_census.set_index(keys='Gene Symbol', inplace=True)

	scores = cancer_gene_census['Score'].to_dict()

	del cancer_gene_census
	return scores


def get_scores_from_predicted_drivers(drivers:str) -> dict:
	'''
	Retrieves pairs of gene-scores based on driver prediction data. Returns a 
	dict. of key-value pairs where key=gene name, value= gscore based on
	Tamborero et.al.(2013)
	'''

	drivers = pd.read_csv(drivers, sep=',',
						  index_col='Gene Symbol',
						  usecols=['Gene Symbol', 'Putative Driver Category']
						  )


	driver_to_score = {'High Confidence Driver': 0.1,
					   'Candidate driver': 0.05
					  }
	
	drivers['Score'] = drivers['Putative Driver Category'].map(driver_to_score)
	scores = drivers['Score'].to_dict()

	del drivers
	return scores



def get_scores_from_essenciality(essenciality:str) -> dict:
	'''
	Get scores based on gene essenciality data and returns a dict. of key-value
	pairs where keys are genes and values the gscore according to this source.
	'''

	essencialities = pd.read_csv(essenciality, sep='\t', 
								 index_col='symbol',
								 usecols=['symbol', 'max_min_score']
								 )

	essencialities['Score'] = essencialities['max_min_score'] * 0.4

	scores = essencialities['Score'].to_dict()

	del essencialities
	return scores

def get_scores_from_oncoscape(oncoscape:str) -> dict:
	'''
	Get scores based on oncoscape data and returns a dict. of  key-value
	pairs where keys are genes and values the gscore according to this source.
	'''

	oncoscape = pd.read_csv(oncoscape,
							sep='\t',
							index_col='SYMBOL',
						    usecols=['SYMBOL', 'MAX_ONC_TSG']
							)

	# divide() method is another alias for truediv() in addition to div()
	oncoscape['Score'] = oncoscape['MAX_ONC_TSG'].divide(4) * 0.3
	scores = oncoscape['Score'].to_dict()

	del oncoscape
	return scores


if __name__ == "__main__":
	main(sys.argv[1],sys.argv[2],sys.argv[3],
		 sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7]
		 )

