'''
This script is part of CI unit tests suite. It checks wether the pipeline
is working as intended against a mock dataset.
'''
import sys
import pandas as pd


def main(summary: pd.DataFrame):

	print('Begin testing...')

	summary = pd.read_csv(summary, sep=',')

	'''
	KRAS   should be GoF
	CDKN2A should be LoF     
	FRAS   should be unknown <---- maf, not in CGC
	APC    should be unknown <---- maf
	TP53   should be LoF     <---- cnv
	'''
	assert summary[summary['Hugo_Symbol'] == 'KRAS']['Consequence'].all() == 'GoF',   "KRAS missclassified"
	assert summary[summary['Hugo_Symbol'] == 'CDKN2A']['Consequence'].all() == 'LoF', "CDKN2A missclassified"
	assert summary[summary['Hugo_Symbol'] == 'FRAS']['Consequence'].all() == 'Unknown', "FRAS missclassified"
	assert summary[summary['Hugo_Symbol'] == 'APC']['Consequence'].all() == 'Unknown', "APC missclassified"
	assert summary[summary['Hugo_Symbol'] == 'TP53']['Consequence'].all() == 'GoF', "TP53 missclassified"

	print('Done')


if __name__ == "__main__":
	main(sys.argv[1])