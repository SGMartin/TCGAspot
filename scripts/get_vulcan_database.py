#!/usr/bin/env python3

"""
Section: SL in TCGA data

Purpose: Get a local copy of vulcan db. 

"""
import pandas as pd
import requests

#TODO: Optimize queries. Recycle gene lists
def main():

	# Snakemake workflow var.
	genes_table_dir  = snakemake.output[0]
	vulcan_table_dir = snakemake.output[1]

	genes_table  = build_vulcan_genes_database()
	genes_table.to_csv(genes_table_dir, sep=',', index=False)
	
	del genes_table

	vulcan_table = build_vulcan_database()
	vulcan_table.to_csv(vulcan_table_dir, sep=',', index=False)


def build_vulcan_genes_database() -> pd.DataFrame:

	rows = list()

	a_genes = get_available_gene_list('A')
	b_genes = get_available_gene_list('B')
	
	all_genes = a_genes + b_genes

	for gene in all_genes:
		print('Getting details for',gene)
		gene_details = get_gene_details(gene)
		rows.append(gene_details)
	
	genes_table = pd.concat([pd.DataFrame([i], columns=['gene', 'gscore', 'druggable']) for i in rows], ignore_index=True)

	genes_table['A'] = genes_table['gene'].isin(a_genes)
	genes_table['B'] = genes_table['gene'].isin(b_genes)

	return genes_table

	
def build_vulcan_database() -> pd.DataFrame:

	VulcanTable = pd.DataFrame(columns=['geneA','tissue','alteration',
										'geneB', 'drug', 'Dscore', 'LINCS']
							  )
	genes_of_interest = get_available_gene_list('A')

	total_gene_count = len(genes_of_interest)
	queried = 1

	for a_gene in genes_of_interest:
		
		print('Querying', queried, ' of', total_gene_count)

		ask_to = 'http://vulcanspot.org/api/genes/' + a_gene + '/treatments'
		g_response = requests.get(ask_to).json()

		# input gene available tissues
		a_gene_contexts = g_response['data'][a_gene].keys()

		for context in a_gene_contexts:
			alteration = list(g_response['data'][a_gene][context]['alterations'].keys())[0]

			genes_B = g_response['data'][a_gene][context]['alterations'][alteration].keys()

			for genetic_dep in genes_B:

				query = g_response['data'][a_gene][context]['alterations'][alteration][genetic_dep]

				# 0 if not found
				rnai_score   = query.get('evidence',{}).get('RNAi', {}).get('score', 0)
				crispr_score = query.get('evidence',{}).get('CRISPR', {}).get('score', 0)
				
				if rnai_score >= 0.00 or crispr_score >= 0.00:
					if 'drugs' in query:
						for drug, drugData in query['drugs'].items():
							if drug != 'validated':
								Dscore = 0
								LINCS  = 0
                                                                                
								if drugData is not False: 
									if 'PANDRUGS' in drugData:
										Dscore = drugData['PANDRUGS'].get('score')
                                            
								if 'LINCS' in drugData:
									LINCS  = drugData['LINCS'].get('score')

								VulcanTable = VulcanTable.append({
                                            'geneA': a_gene,
                                            'tissue': context,
											'alteration': alteration,
                                            'geneB': genetic_dep,
											'drug': drug,
                                            'Dscore': Dscore,
                                            'LINCS' : LINCS},
											ignore_index=True
											 )

		queried +=1

	return VulcanTable
			
def get_available_gene_list(gene_class:str = 'A') ->list():
	'''
	Returns a list of available genes of a certain class in VulcanSpot. gene_class
	can be either A or B. Defaults to A. 
	'''
	list_of_genes = list()
	
	vulcan_request = F'http://vulcanspot.org/api/genes?class={gene_class}'
	
	vulcan_request = requests.get(vulcan_request).json()

	for gene in vulcan_request['data']:
		list_of_genes.append(gene['key'])
	
	return list_of_genes

def get_gene_details(input_gene: str) -> dict():
	'''
	Asks VulcanSpot for a certain gene details and returns a dictionary of 
	key-value pairs with relevant info.
	'''
	gene_details = dict()

	query = F'http://vulcanspot.org/api/genes/{input_gene}'
	response = requests.get(query).json()
	
	gene_details['gene'] 	  = response['data'][0]['symbol']
	gene_details['gscore']	  = response['data'][0]['score']
	gene_details['druggable'] = response['data'][0]['druggable']

	return gene_details


if __name__ == "__main__":
	main()