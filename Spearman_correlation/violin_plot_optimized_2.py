# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 15:15:11 2026

@author: miabr
"""

__usage__ = """ python3 violin_plot_optimized_2.py
					--TF_TG_list <LIST WITH TF AND RELATED TG>
					--spearman_corr_dir <PATH TO SPEARMAN CORRELATION DIR>
					--out_dir <OUTPUT DIRECTORY>
					--organism <plants OR fungi OR bacteria OR animals>
					OPTIONAL:
					--n <NUMBER OF RANDOM SPECIES TO ANALYZE>
					--sig-only <ONLY SIGNIFICANT TF-TG PAIRS, NO VALUE NEEDED>
					--p-col <p_value_adj (DEFAULT) OR p_value, FOR --sig-only>
					--alpha <SIGNIFICANCE THRESHOLD; DEFAULT: 0.05>
					"""

### IMPORTS ###
import math
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys
import pandas as pd
from collections import defaultdict
import random

def dict_TF_TG(input_file, organism):
	with open(input_file, 'r') as f:
		TF = set()
		TG = set()
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
			TF.add(tf)
			TG.add(tg)
			if tg_name == '':
				tg_name = tg
			TF_TG_dict[tf].append(tg_name)
			TG_dict[tg].append(tg_name)
			line = f.readline()
	return TF_TG_dict

def build_pair_dataframe(TF_TG_dict, spearman_corr, n, sig_only=False, p_col='p_value_adj', alpha=0.05):
	
	desired_pairs = set()
	for tf, tgs in TF_TG_dict.items():
		for tg in tgs:
			desired_pairs.add((tf,tg))
	
	files = []
	for f in os.listdir(spearman_corr):
		if '_spearman_corr.txt' in f:
			files.append(f)
	if n:
		files = random.sample(files,n)
	
	frames = []
	for file in files:
		path = os.path.join(spearman_corr, file)
		species = file.replace('_spearman_corr.txt', '')
		
		header_cols = pd.read_csv(path, sep='\t', nrows=0).columns
		required_cols = {'TF_ref', 'TG_ref', 'spearman_r', 'TF_id', 'TG_id'}
		if not required_cols.issubset(header_cols):
			print(f"WARNING: {file} skipped - missing columns {required_cols - set(header_cols)}")
			continue
		
		file_p_col = None
		if sig_only:
			if p_col in header_cols:
				file_p_col = p_col
			elif 'p_value' in header_cols:
				print(f"WARNING: {file}: column '{p_col}' not found -> 'p_value' used instead!")
				file_p_col = 'p_value'
			else:
				print(f"WARNING: {file}: neither '{p_col}' nor 'p_value' -> ignoring --sig-only for this file!")
				
		usecols = ['TF_ref', 'TG_ref', 'spearman_r', 'TF_id', 'TG_id']
		if file_p_col is not None:
			usecols.append(file_p_col)
			
		df_file = pd.read_csv(
			path,
			sep='\t',
			usecols=usecols,
			dtype={'TF_ref': str, 'TG_ref': str, 'TF_id': str, 'TG_id': str}
			)
		
		pair_series = pd.Series(
			list(zip(df_file['TF_ref'], df_file['TG_ref'])),
			index=df_file.index
			)
		mask = pair_series.isin(desired_pairs)
		matched = df_file.loc[mask].dropna(subset=['spearman_r'])
		matched = matched[matched['TF_id'] != matched['TG_id']]
		
		if file_p_col is not None:
			matched = matched[matched[file_p_col] <= alpha]
			
		if not matched.empty:
			frames.append(pd.DataFrame({
				'TF-TG pair' : matched['TF_ref'] + '-' + matched['TG_ref'],
				'spearman correlation' : matched['spearman_r'].astype(float),
				'species' : species
				}))
	if frames:
		return pd.concat(frames, ignore_index=True)
	return pd.DataFrame(columns=['TF-TG pair', 'spearman correlation', 'species'])

def violin_plot(df, output_dir, means_file, organism, group_size=10):
	pairs = df['TF-TG pair'].drop_duplicates().tolist()
	
	n_plots = math.ceil(len(pairs) / group_size)
	
	df['spearman_abs'] = df['spearman correlation'].abs()
	# Calculate means
	means_df = df.groupby('TF-TG pair').agg(
		mean=('spearman correlation', 'mean'), 
		median=('spearman correlation', 'median'), 
		n=('spearman correlation', 'count'),			# Number of TF-TG pairs with orthologues
		m=('species', 'nunique')						# Number of species with the sampe TF-TG pair
	).reindex(pairs)
	
	
	abs_stats = df.groupby('TF-TG pair')['spearman_abs'].agg(
		mean_abs='mean', median_abs='median'
		).reindex(pairs)
	means_df = means_df.join(abs_stats)
	
	if means_file is not None:
		means_df.to_csv(means_file, sep='\t')
		print(f"Means saved as: {means_file}")
	
	for i in range(n_plots):
		start = i * group_size
		end = start + group_size
		
		pair_subset = pairs[start:end]
		
		df_subset = df[df['TF-TG pair'].isin(pair_subset)].copy()
		
		plt.figure(figsize=(14, 7))
		
		color_map = {
			"plants" : "lightgreen",
			"bacteria" : "lightskyblue",
			"fungi" : "plum",
			"animals" : "orange"}
		
		color = color_map[organism]
		
		ax = sns.violinplot(
			data = df_subset,
			x='TF-TG pair',
			y='spearman correlation',
			color=color,
			order=pair_subset,
			inner='box',
			cut=0
			)
		
		plt.ylim(-1, 1)
		
		# n-Werte über jeder Violine einblenden
		for x_pos, pair in enumerate(pair_subset):
			n = means_df.loc[pair, 'n'] 
			m = means_df.loc[pair, 'm']
			ax.text(
				x_pos, 1.1, 	# x = Position der Violine, y = knapp über ylim
				f"n={int(n)}\nspecies={int(m)}",
				ha='center', va='bottom',
				fontsize=13, rotation=0
				)
		
		locs, labels = plt.xticks()						# locs = position; labels
		label = []
		for t in labels:
			element = '-'.join(t.get_text().split('-')[-2:])
			label.append(element)
		plt.xticks(locs, label, fontsize=18, rotation=90, fontstyle='italic')
		plt.yticks(fontsize=18)
		plt.xlabel('TF-TG pair', fontsize=18)
		plt.ylabel('Spearman correlation', fontsize=18)
		plt.title(f'TF-TG pairs {start + 1}-{min(end, len(pairs))}', pad=60, fontsize=18)
		plt.tight_layout()
		
		
		if output_dir is not None:
			output_file = f"{output_dir}/violinplot_{i + 1:03d}.png"
			plt.savefig(output_file, dpi=300, bbox_inches="tight")
			print(f"Violin plot saved as: {output_file}")
			
		plt.show()
		plt.close()
	
def main(arguments):
	
	TF_TG_list = arguments[arguments.index('--TF_TG_list') + 1]
	organism = arguments[arguments.index('--organism') + 1]
	spearman_corr = arguments[arguments.index('--spearman_corr_dir') + 1]
	out_dir = arguments[arguments.index('--out_dir') + 1]
	
	if '--n' in arguments:
		n = int(arguments[arguments.index('--n') + 1])
	else:
		n = None
	
	if '--sig-only' in arguments:
		sig_only = True
	else:
		sig_only = False
	
	if '--p-col' in arguments:
		p_col = arguments[arguments.index('--p-col') + 1]
	else:
		p_col = 'p_value_adj'
	if p_col not in ('p_value_adj', 'p_value'):
		sys.exit(f"--p-col has to be 'p_value_adj' or 'p_value', not '{p_col}'")
	
	if '--alpha' in arguments:
		alpha = float(arguments[arguments.index('--alpha') + 1])
	else:
		alpha = 0.05
	
	TF_TG_dict = dict_TF_TG(TF_TG_list, organism)
	df = build_pair_dataframe(TF_TG_dict, spearman_corr, n, sig_only=sig_only, p_col=p_col, alpha=alpha)
	
	means_file = os.path.join(out_dir, f"{organism}_means.tsv")
	violin_plot(df, out_dir, means_file, organism, group_size=10)
	
	print('Done.')
	
if '--TF_TG_list' in sys.argv and '--spearman_corr_dir' in sys.argv and '--out_dir' in sys.argv and '--organism' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )