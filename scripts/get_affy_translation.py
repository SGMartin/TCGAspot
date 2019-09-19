'''
Affy translation module: This scripts generates a new master table to translate
ENSEMBL id to Hugo Symbol. It shouldn't be run on an analysis basis and is 
included for completeness purposes.
'''
import mygene
import pandas as pd

def main():
	"""
	Translates Ensembl IDs to HUGO using mygene queries.
	"""
	# SNAKEMAKE INPUT #
	snp_array_db 		= snakemake.input[0]

	# SNAKEMAKE OUTPUT #
	translated_affy_snp = snakemake.output[0]
	
	# Load data and trim version number from ensembl id

	print("Annotating Ensembl IDs to HUGO...")

	ensembl_snp_array = pd.read_csv(snp_array_db, sep='\t')

	ensembl_snp_array['Gene Symbol'] = ensembl_snp_array['Gene Symbol'].str.split('.').str[0].str.strip()

	MyGene     = mygene.MyGeneInfo()
	gene_query = ensembl_snp_array['Gene Symbol']

	annotated_names = MyGene.getgenes(gene_query,
									  fields='symbol',
									  as_dataframe=True
									 )
	
	ensembl_snp_array['Gene Symbol'] = annotated_names['symbol'].values

	del annotated_names
	
	ensembl_snp_array.to_csv(translated_affy_snp, sep=',', index=False)

if __name__ == "__main__":
	main()

