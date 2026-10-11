# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 11:29:21 2026

@author: miabr
"""

"""
This script (def plot_first_n_tf_tg_pairs) was developed with the help of AI tools (Claude). 
All sections of code were carefully reviewed, tested, and scientifically validated.
"""


__usage__ = """ python3 spearman_correlation.py
					--orthologe_TF <DIRECTORY TO THE TF ORTHOLOGE FILES (z.B. MYB)>
					--orthologe_TF_2 <OPTIONAL: ZWEITES TF ORTHOLOGE DIRECTORY, z.B. bHLH>
					--species <FILE WITH SPECIES FOR ANALYSIS>
					--TF_TG_list <LIST WITH TF AND RELATED TG>
					--orthoFinder_file <SINGLE ORTHOFINDER RESULTS FILE>
					--counts_dir <DIRECTORY TO COUNTTABLES>
					--ref <REFERENCE SPECIES NAME e.g. Arabidopsis_thaliana>
					--out_dir <OUTPUT DIRECTORY>
					--plot <OPTIONAL: erstellt zusätzlich Scatter-Plots für die Referenz-Spezies>
					--organism <plants OR fungi OR bacteria OR animals>
					--p_correction <CORRECTION p-value DEFAULT: fdr_bh; ALTERNATIVE: bonferroni>
					"""

### pip install scipy ###
### IMPORTS ###
import pandas as pd
import os, sys
import matplotlib.pyplot as plt
from collections import defaultdict
import re
from statsmodels.stats.multitest import multipletests
from scipy.stats import spearmanr
import numpy as np

def species_list(species_file):
	with open(species_file, 'r') as f:
		species = []
		line = f.readline()
		while line:
			line = line.strip()
			species.append(line)
			
			line = f.readline()
			
	return species

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
	return TF_TG_dict, TF, TG, TG_dict

def find_TF_file(orthologe_TF_path, name):
	
	for entry in os.listdir(orthologe_TF_path):
		if entry == name:
			results_path = os.path.join(orthologe_TF_path, entry, "RESULTS")
			if os.path.isdir(results_path):
				print(results_path)
				return results_path
		else:
			parts = entry.split('_', 1)
			if re.match(r"^\d{5}$", parts[0]):
				if len(parts) < 2:
					continue
				if parts[1] == name:
					results_path = os.path.join(orthologe_TF_path, entry, "RESULTS")
					if os.path.isdir(results_path):
						print(results_path)
						return results_path
			else:
				parts = entry.split('.')
				if len(parts) < 2:
					continue
				if parts[0] == name:
					results_path = os.path.join(orthologe_TF_path, entry, "RESULTS")
					if os.path.isdir(results_path):
						print (results_path)
						return results_path
	return None

def find_result_file(results_path):

	if results_path is None:
		print("No results_path")
		return None
	for file in os.listdir(results_path):
		parts = file.split('_')
		if parts[1:3] == ['group', 'around']:
			orthologe_file = os.path.join(results_path, file)
			print(orthologe_file)
			return orthologe_file
		elif parts[2:5] == ['group', 'around', 'reference.txt']:
			orthologe_file = os.path.join(results_path, file)
			print(orthologe_file)
			return orthologe_file

def translate_TF(orthologe_file, TF):
	with open(orthologe_file, 'r') as f: 
		new_TFs = defaultdict(list)
		f.readline()
		line = f.readline()
		while line:
			line = line.strip().split('\t')
			ref = line[0].replace('/', '|')
			if len(line) < 3:
				line = f.readline()
				continue
			else:
				new = line[2]
			if ref in TF:
				print(ref)
				if ';' in new:
					new = new.split(';')
					for element in new:
						new_TFs[ref].append(element)
				else:
					new_TFs[ref].append(new)
			line = f.readline()
	return new_TFs

def translate_TG(orthoFinder_file, TG, name):
	with open(orthoFinder_file, 'r') as f: 
		new_TGs = defaultdict(list)
		f.readline()
		line = f.readline()
		while line:
			line = line.strip().split('\t')
			if len(line) < 2:
				print(f"Malformed line: {line}")
			else:
				species = line[1].split('.')[0]
				ref = line[2]
				new = line[3]
				if species == name:
					ref = ref.split(', ')
					for ref_id in ref:
						if ref_id in TG:
							if ', ' in new:
								new = new.split(', ')
								for element in new:
									new_TGs[ref_id].append(element)
							else:
								new_TGs[ref_id].append(new)
			line = f.readline()
	return new_TGs

def new_tpm_file(tpms_file, out, new_TFs, new_TGs, name, reference, TG_dict):
	print(new_TFs)
	with open(tpms_file, 'r') as f:
		header = f.readline()
		header = header.strip().split('\t')
		new = ['RefMember', 'Orthologue'] + header[1:]
		new_header = '\t'.join(new)
		out.write(f"{new_header}\n")
		
		line = f.readline()
		
		while line:
			TF_found = False
			line = line.strip().split('\t')
			gene = line[0]
			if gene.startswith('lcl_'):
				print("ID prefix different (lcl_)!")
				gene = gene.replace('_', '|', 1)
			samples = '\t'.join(line[1:])
			for key, values in new_TFs.items():
				for value in values:
					if gene == value:
						out.write(f"{key}\t{gene}\t{samples}\n")
						TF_found = True
			if not TF_found:
				for key, values in new_TGs.items():
					for value in values:
						if gene == value:
							for element in TG_dict[key]:
								out.write(f"{element}\t{gene}\t{samples}\n")
			line = f.readline()

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
	pairs = []
	for tf, tgs in TF_TG_dict.items():
		for tg in tgs:
			pairs.append((tf, tg))
			if len(pairs) >= n:
				break
		if len(pairs) >= n:
			break
	
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
			
				plt.figure(figsize=(5, 5))
				plt.scatter(tf_num, tg_num, s=10, alpha=0.7)
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

def calculate_spearman_correlations(compressed_tpm_file, TF_TG_dict, out_file, min_samples=5, p_correction='fdr_bh'):
	"""
	Berechnet die Spearman-Korrelation für alle TF-TG-Paare aus einer
	compressed TPM-Datei, korrigiert die p-Werte für multiples Testen 
	(Standard: Benjamin-Hochberg FDR; alternativ 'bonferroni') und 
	schreibt das Ergebnis als Tabelle raus.
	"""
	df = pd.read_csv(compressed_tpm_file, sep='\t')

	if 'RefMember' not in df.columns or 'Orthologue' not in df.columns:
		raise ValueError("Compressed TPM file must contain 'RefMember' and 'Orthologue' columns.")
	
	n_before = len(df)
	df = df.drop_duplicates()
	if len(df) < n_before:
		print(f"Warning: {n_before - len(df)} duplicate rows removed from compressed TPM file.")
	
	df_indexed = df.set_index('Orthologue')
	
	ref_to_ids = df_indexed.groupby('RefMember').groups
	
	results = []

	for tf_ref, tgs in TF_TG_dict.items():
		tf_ids = list(ref_to_ids.get(tf_ref, []))
		if not tf_ids:
			continue

		for tg_ref in tgs:
			tg_ids = list(ref_to_ids.get(tg_ref, []))
			if not tg_ids:
				continue

			for tf_id in tf_ids:
				for tg_id in tg_ids:
					tf_values = df_indexed.loc[tf_id]
					tg_values = df_indexed.loc[tg_id]

					tf_num = tf_values.drop('RefMember', errors='ignore').apply(pd.to_numeric, errors='coerce').dropna()
					tg_num = tg_values.drop('RefMember', errors='ignore').apply(pd.to_numeric, errors='coerce').dropna()

					common_idx = tf_num.index.intersection(tg_num.index)
					if len(common_idx) < min_samples:
						continue

					corr, p_value = spearmanr(tf_num[common_idx], tg_num[common_idx])

					results.append({
						'TF_ref': tf_ref,
						'TF_id': tf_id,
						'TG_ref': tg_ref,
						'TG_id': tg_id,
						'spearman_r': corr,
						'p_value': p_value,
						'n_samples': len(common_idx)
					})

	result_df = pd.DataFrame(results)
	
	if not result_df.empty:
		
		valid_mask = result_df['p_value'].notna()
		n_dropped = (~valid_mask).sum()
		if n_dropped > 0:
			print(f"Warning: {n_dropped} Paare mit NaN p_value (z.B. konstante TPM-Werte) "
					f"werden von der Multiple-Testing-Korrektur ausgeschlossen.")
		
		result_df['p_value_adj'] = np.nan
		result_df['significant'] = False
		
		if valid_mask.any():
			reject, p_adj, _, _ = multipletests(result_df.loc[valid_mask, 'p_value'], method=p_correction)
			result_df.loc[valid_mask, 'p_value_adj'] = p_adj
			result_df.loc[valid_mask, 'significant'] = reject
	else:
		result_df['p_value_adj'] = []
		result_df['significant'] = []
		
	result_df.to_csv(out_file, sep='\t', index=False)
	print(f"Calculating Spearman correlation: {len(result_df)} pairs -> '{out_file}'"
			f"(p-value correction: {p_correction})")
	return result_df

def main(arguments):
	
	orthologe_TF = arguments[arguments.index('--orthologe_TF') + 1]
	species_file = arguments[arguments.index('--species') + 1]
	TF_TG_list = arguments[arguments.index('--TF_TG_list') + 1]
	orthoFinder_file = arguments[arguments.index('--orthoFinder_file') + 1]
	counts_dir = arguments[arguments.index('--counts_dir') + 1]
	reference = arguments[arguments.index('--ref') + 1]
	out_dir = arguments[arguments.index('--out_dir') + 1]
	organism = arguments[arguments.index('--organism') + 1]
	
	if '--p_correction' in arguments:
		p_correction = arguments[arguments.index('--p_correction') + 1]
	else:
		p_correction = 'fdr_bh'
		
	if '--orthologe_TF_2' in arguments:
		orthologe_TF_2 = arguments[arguments.index('--orthologe_TF_2') + 1]
		tf_directories = [orthologe_TF, orthologe_TF_2]
	else:
		tf_directories = [orthologe_TF]
	
	species = species_list(species_file)
	print(species)
	TF_TG_dict, TF, TG, TG_dict = dict_TF_TG(TF_TG_list, organism)
	
	for name in species:
		print(name)
		x = 0
		for directory in tf_directories:
			x += 1
			results_path = find_TF_file(directory, name)
			orthologe_file = find_result_file(results_path)
			if x == 1:
				new_TFs_1 = translate_TF(orthologe_file, TF)
				print("new_TFs_1 created")
			if x == 2:
				new_TFs_2 = translate_TF(orthologe_file, TF)
				print("new_TFs_2 created")
				new_TFs = defaultdict(list)
				for d in [new_TFs_1, new_TFs_2]:
					for key, values in d.items():
						new_TFs[key].extend(values)
		if x == 1:
			new_TFs = new_TFs_1
			print(new_TFs)
		if name != reference:
			new_TGs = translate_TG(orthoFinder_file, TG, name)
		else:
			new_TGs = {}
			for element in TG:
				new_TGs[element] = [element]
				
		os.makedirs(out_dir, exist_ok=True)
		out_file = os.path.join(out_dir, f"{name}_compressed.tpms.txt")
		counts_file = os.path.join(counts_dir, name, f"{name}.tpms.txt")
		if not os.path.isfile(counts_file):
			counts_file = os.path.join(counts_dir, name, f"{name}.tpms.tsv")
		with open(out_file, 'w') as out:
			new_tpm_file(counts_file, out, new_TFs, new_TGs, name, reference, TG_dict)
		
		corr_out_file = os.path.join(out_dir, f"{name}_spearman_corr.txt")
		calculate_spearman_correlations(out_file, TF_TG_dict, corr_out_file, p_correction=p_correction)
	print("New counttables are done.")
	
	if '--plot' in arguments:
		n = int(arguments[arguments.index('--plot') + 1])
		compressed_tpm_file = os.path.join(out_dir, f'{reference}_compressed.tpms.txt')
		plots_dir = os.path.join(out_dir, 'scatter_plots')
		os.makedirs(plots_dir, exist_ok=True)
		plot_first_n_tf_tg_pairs(compressed_tpm_file, TF_TG_dict, plots_dir, n)
		
		print("Plots are done.")
		
if '--orthologe_TF' in sys.argv and '--species' in sys.argv and '--TF_TG_list' in sys.argv and '--orthoFinder_file' in sys.argv and '--counts_dir' in sys.argv and '--ref' in sys.argv and '--out_dir' in sys.argv and '--organism' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )