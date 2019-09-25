#!/usr/bin/env python3
'''
Expression filter module: This script
'''

import sys

import pandas as pd
import numpy  as np
from scipy import stats
from scipy.stats import zscore

def main():
	'''
	'''

def get_upper_quartile_expression(melted_data: pd.DataFrame) -> pd.DataFrame:
	'''
	'''
	# calculate intra-aliquot .95 percentile and convert it to a dict of 
	# aliquot-percentile key value pairs
	limit = melted_data.groupby('aliquot')['zscore'].quantile(q=0.95).to_dict()

	melted_data['threshold'] = melted_data['aliquot'].map(limit)

	melted_data['over_threshold'] = melted_data['zscore'] >= melted_data['threshold']


def melt_zscore(raw_data: pd.DataFrame) -> pd.DataFrame:
	'''
	'''
	# melt for easier manip.
	melted = raw_data.melt(id_vars='Ensembl_ID',
						   var_name='aliquot',
						   value_name='expression')
	
	melted['zscore'] = melted.groupby('Ensembl_ID')['expression'].transform(lambda x: zscore(x))

	# Drop rows with NaNs produced by all values being 0,
	# so STD = 0 and zscore = x / 0
	melted.dropna(axis='index', how='any', inplace=True)

	return melted
 

if __name__ == "__main__":
	main(sys.argv[1:])	