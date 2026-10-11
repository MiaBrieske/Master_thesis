# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 13:42:22 2026

@author: miabr
"""

"""
This script (def plot_first_n_tf_tg_pairs) was developed with the help of AI tools (Claude). 
All sections of code were carefully reviewed, tested, and scientifically validated.
"""

__usage__ = """ python3 scatter_plots.py
					--TF_TG_list <LIST WITH TF AND RELATED TG>
					--compressed_counts <DIRECTORY TO COUNTTABLES>
					--out_dir <OUTPUT DIRECTORY>
					--plot <NUMBER OF SCATTER PLOTS TO CREATE>
					--organism <plants OR fungi OR bacteria OR animals>
					"""

### IMPORTS ###
import pandas as pd
import os, sys
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
from scipy.stats import gaussian_kde
import random

def dict_TF_TG(input_file, organism):
	with open(input_file, 'r') as f:
		TF = []
		TG = []
		TF_TG_dict = defaultdict(list)
		TG_dict = defaultdict(list)
		f.readline()						#header
		line = f.readline()
		while line:
			line = line.strip().split('\t')
			if len(line) < 2:
				print(f"Skipping malformed line (fewer than 2 columns): {line}")
				line = f.readline()
				continue
			tf = line[0]
			tg = line[1]
			if organism == 'plants':
				tg = tg.upper()
			tg_name = line[2]
			if tf not in TF:
				TF.append(tf)
			if tg not in TG:
				TG.append(tg)
			if tg_name == '':
				tg_name = tg
			TF_TG_dict[tf].append(tg_name)
			TG_dict[tg].append(tg_name)
			line = f.readline()
	return TF_TG_dict

def plot_first_n_tf_tg_pairs(compressed_tpm_file, TF_TG_dict, plots_dir, n):

	df = pd.read_csv(compressed_tpm_file, sep='\t')

	if 'RefMember' not in df.columns or 'Orthologue' not in df.columns:
		raise ValueError("Compressed TPM file must contain 'RefMember' and 'Orthologue' columns.")
		
	n_before = len(df)
	df = df.drop_duplicates()
	if len(df) < n_before:
		print(f"Warning: {n_before - len(df)} duplicate rows removed from compressed TPM file.")
	# Set index for easy lookup
	df_indexed = df.set_index('Orthologue')
	print(df_indexed)
	
	ref_to_ids = df_indexed.groupby('RefMember').groups
	# Flatten TF-TG pairs
	all_pairs = []
	for tf, tgs in TF_TG_dict.items():
		for tg in tgs:
			all_pairs.append((tf, tg))
	
	if n is not None and n < len(all_pairs):
		pairs = random.sample(all_pairs, n)
	else:
		pairs = all_pairs
	
	missing_tf_refs = []
	missing_tg_refs = []
	matching_tf_refs = []
	matching_tg_refs = []
	os.makedirs(plots_dir, exist_ok=True)
	plotted = 0

	for tf_ref, tg_ref in pairs:
		tf_ids = list(ref_to_ids.get(tf_ref, []))
		tg_ids = list(ref_to_ids.get(tg_ref, []))

		if not tf_ids:
			print(f"Skip {tf_ref}-{tg_ref}: TF not found in compressed TPM file.")
			if tf_ref not in missing_tf_refs:
				missing_tf_refs.append(tf_ref)
			continue
		elif tf_ids:
			matching_tf_refs.append(tf_ref)
			
		if not tg_ids:
			print(f"Skip {tf_ref}-{tg_ref}: TG not found in compressed TPM file.")
			if tg_ref not in missing_tg_refs:
				missing_tg_refs.append(tg_ref)
			continue
		elif tg_ids:
			matching_tg_refs.append(tg_ref)

		for tf_id in tf_ids:
			for tg_id in tg_ids:
				tf_values = df_indexed.loc[tf_id]
				tg_values = df_indexed.loc[tg_id]
	
				# Select only numeric sample columns (exclude RefMember if present as column)
				numeric_cols = tf_values.drop('RefMember', errors='ignore').apply(pd.to_numeric, errors='coerce')
				tf_num = numeric_cols.dropna()
				tg_num = tg_values.loc[tf_num.index].drop('RefMember', errors='ignore').apply(pd.to_numeric, errors='coerce').dropna()
			
				if len(tf_num) == 0 or len(tg_num) == 0:
					continue
				
				x = tf_num.to_numpy(dtype=float)
				y = tg_num.to_numpy(dtype=float)
				
				plt.figure(figsize=(5, 5))
				
				kde_ok = False
				if len(x) >= 3 and np.ptp(x) > 0 and np.ptp(y) > 0:
					try:
						kde = gaussian_kde([x, y])
						xi, yi = np.mgrid[x.min():x.max():100j, y.min():y.max():100j]
						zi = kde(np.vstack([xi.flatten(), yi.flatten()]))
						plt.contourf(xi, yi, zi.reshape(xi.shape), cmap='Blues', alpha=0.6)
						kde_ok = True
					except np.linalg.LinAlgError:
						print(f"WARNING: KDE for {tf_ref} vs. {tg_ref} not possible, only scatter.")
				
				if kde_ok:
					scatter_color = 'black'
					scatter_alpha = 0.3
					s = 8
				else:
					scatter_color = None
					scatter_alpha = 0.7
					s = 10
				plt.scatter(x, y, s=s, alpha=scatter_alpha, color=scatter_color)
				plt.xlabel(f"{tf_ref} Expression")
				plt.ylabel(f"{tg_ref} Expression")
				plt.title(f"{tf_ref} vs. {tg_ref}")
				plt.tight_layout()
				if '/' in tf_ref or ' ' in tf_ref or '|' in tf_ref or ';' in tf_ref or '&' in tf_ref:
					tf_ref = tf_ref.replace('/', '_')
					tf_ref = tf_ref.replace(' ', '_')
					tf_ref = tf_ref.replace('|', '_')
					tf_ref = tf_ref.replace(';', '_')
					tf_ref = tf_ref.replace('&', '_')
				if '/' in tg_ref or ' ' in tg_ref or '|' in tg_ref or ';' in tg_ref or '&' in tg_ref:
					tg_ref = tg_ref.replace('/', '_')
					tg_ref = tg_ref.replace(' ', '_')
					tg_ref = tg_ref.replace('|', '_')
					tg_ref = tg_ref.replace(';', '_')
					tg_ref = tg_ref.replace('&', '_')
				out_path = os.path.join(plots_dir, f"{tf_ref}_vs_{tg_ref}.png")
				plt.savefig(out_path)
				plt.close()
				plotted += 1

	print(f"Finished: processed {len(pairs)} TF-TG pairs, created {plotted} scatter plots in '{plots_dir}'.")
	missing_tf_ref = []
	missing_tg_ref = []
	for element in missing_tf_refs:
		if element not in matching_tf_refs:
			missing_tf_ref.append(element)
	for element in missing_tg_refs:
		if element not in matching_tg_refs:
			missing_tg_ref.append(element)
	print(f"Missing TFs: {missing_tf_ref}")
	print(f"Missing TGs: {missing_tg_ref}")

def main(arguments):
	
	TF_TG_list = arguments[arguments.index('--TF_TG_list') + 1]
	compressed_counts = arguments[arguments.index('--compressed_counts') + 1]
	out_dir = arguments[arguments.index('--out_dir') + 1]
	organism = arguments[arguments.index('--organism') + 1]
	n = int(arguments[arguments.index('--plot') + 1])
	
	TF_TG_dict = dict_TF_TG(TF_TG_list, organism)
	
	os.makedirs(out_dir, exist_ok=True)
	plot_first_n_tf_tg_pairs(compressed_counts, TF_TG_dict, out_dir, n)
		
	print("Plots are done.")
		
if '--TF_TG_list' in sys.argv and '--compressed_counts' in sys.argv and '--out_dir' in sys.argv and '--organism' in sys.argv and '--plot' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )