"""
Section: SL in TCGA data

Purpose: Get a local copy of vulcan db. 

"""
import pandas as pd
import requests

def main():

	# Snakemake workflow var.
	where_to_save = snakemake.output[0]

	v_table = build_vulcan_database()
	v_table.to_csv(where_to_save, sep=',', index=False)


def build_vulcan_database() -> pd.DataFrame:

	VulcanTable = pd.DataFrame(columns=['geneA','tissue','alteration',
										'geneB', 'drug', 'Dscore', 'LINCS']
							  )

	genes_of_interest = list()

	response = requests.get('http://vulcanspot.org/api/genes?class=A').json()

	for gene in response['data']:
		genes_of_interest.append(gene['key'])

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
			


if __name__ == "__main__":
	main()