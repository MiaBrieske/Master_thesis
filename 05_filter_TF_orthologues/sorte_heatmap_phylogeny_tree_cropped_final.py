# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 14:22:28 2026

@author: miabr
"""

"""
This script was developed with the help of AI tools (Claude).
All sections of code were carefully reviewed, tested, and scientifically validated.
"""

__usage__ = """ python3 sorte_heatmap_phylogeny_tree_cropped_final.py
					--species_tree <SPECIES TREE FILE FROM ORTHOFINDER>
					--csv_heatmap <CSV FROM INITIAL HEATMAP>
					--species_classification_list <SPECIES ORDER OR FAMILY LIST>
					--out_file <OUTPUT FILE>
					--classification_type <family OR order; HEADER FROM species_classification_list APART FROM species>
					"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from Bio import Phylo
import sys, os

def main(arguments):
	
	SPECIES_TREE_FILE = arguments[arguments.index('--species_tree') + 1]
	PRESENCE_ABSENCE_CSV = arguments[arguments.index('--csv_heatmap') + 1]
	SPECIES_FAMILY_TSV = arguments[arguments.index('--species_classification_list') + 1]
	OUTPUT_FILE = arguments[arguments.index('--out_file') + 1]
	CLASSIFICATION_TYPE = arguments[arguments.index('--classification_type') + 1]
	
	out_dir = OUTPUT_FILE.rsplit('/', 1)[0]
	os.makedirs(out_dir, exist_ok=True)
	
	if '--species_filter_list' in sys.argv:
		SPECIES_FILTER_LIST = arguments[arguments.index('--species_filter_list') + 1]   # eine Art pro Zeile
		USE_SPECIES_FILTER = True
	else:
		USE_SPECIES_FILTER = False
		
	if '--TF_filter' in sys.argv:
		REFMEMBER_FILTER_LIST = arguments[arguments.index('--TF_filter') + 1]
		USE_REFMEMBER_FILTER = True
	else:
		USE_REFMEMBER_FILTER = False
	
	###  ---  Baum laden  ---  ###
	tree = Phylo.read(SPECIES_TREE_FILE, "newick")
	tree.ladderize()
	species_order = [leaf.name.removesuffix(".pep") for leaf in tree.get_terminals()]
	
	###  ---  Presence/Absence-Daten laden  ---  ###
	df = pd.read_csv(PRESENCE_ABSENCE_CSV, index_col=0)
	
	missing = set(species_order) - set(df.columns)
	if missing:
	    print("Missing:", missing)
	
	species_order_full = [s for s in species_order if s in df.columns]
	
	###  ---  Familien-Zuordnung laden  ---  ###
	family_df = pd.read_csv(SPECIES_FAMILY_TSV, sep="\t", header=0, names=["species", CLASSIFICATION_TYPE])
	family_map = {k.removesuffix(".pep"): v for k, v in zip(family_df["species"], family_df[CLASSIFICATION_TYPE])}
	
	missing_fam = [s for s in species_order_full if s not in family_map]
	if missing_fam:
	    print("No {CLASSIFICATION_TYPE} found for:", missing_fam)
	
	###  ---  Familien -> Farben, IMMER aus dem KOMPLETTEN Datensatz  ---  ###
	all_families = sorted(set(family_map[s] for s in species_order_full if s in family_map))
	palette = sns.color_palette("tab20", n_colors=len(all_families))
	family_color = dict(zip(all_families, palette))
	family_index = {fam: i for i, fam in enumerate(all_families)}
	
	###  ---  Optional: nur Arten aus der Filter-Liste behalten  ---  ###
	if USE_SPECIES_FILTER:
	    with open(SPECIES_FILTER_LIST, "r") as f:
	        allowed_species = {line.strip().removesuffix(".pep") for line in f if line.strip()}
	
	    not_in_data = allowed_species - set(species_order_full)
	    if not_in_data:
	        print("Aus der Species-Filter-Liste nicht in den Daten gefunden:", not_in_data)
	
	    species_order_present = [s for s in species_order_full if s in allowed_species]
	    print(f"{len(species_order_present)} Arten nach Filterung übrig.")
	else:
	    species_order_present = species_order_full
	
	df_sorted = df[species_order_present]
	
	###  ---  Optional: nur bestimmte RefMember behalten  ---  ###
	if USE_REFMEMBER_FILTER:
	    with open(REFMEMBER_FILTER_LIST, "r") as f:
	        lines = f.readlines()[1:]
	        allowed_refmembers = [line.strip() for line in lines if line.strip()]
	
	    not_in_data_ref = set(allowed_refmembers) - set(df_sorted.index)
	    if not_in_data_ref:
	        print("Aus der RefMember-Filter-Liste nicht in den Daten gefunden:", not_in_data_ref)
	
	    # Reihenfolge der Liste wird beibehalten (statt Original-Reihenfolge der CSV)
	    refmembers_present = [r for r in allowed_refmembers if r in df_sorted.index]
	    df_sorted = df_sorted.loc[refmembers_present]
	    print(f"{len(refmembers_present)} RefMember nach Filterung übrig.")
	
	###  ---  Matrix: Presence -> Familien-Index (aus dem FIXEN Mapping), Absence -> NaN  ---  ###
	matrix = df_sorted.copy().astype(float)
	for col in matrix.columns:
	    fam = family_map.get(col)
	    fam_idx = family_index[fam] if fam is not None else np.nan
	    matrix[col] = matrix[col].apply(lambda v: fam_idx if v == 1 else np.nan)
	
	cmap = ListedColormap([family_color[f] for f in all_families])
	
	###  ---  Cladogramm-Koordinaten berechnen  ---  ###
	leaf_x = {name: i + 0.5 for i, name in enumerate(species_order_present)}
	
	def get_xy(clade):
	    if clade.is_terminal():
	        name = clade.name.removesuffix(".pep")
	        x = leaf_x.get(name, np.nan)
	        return x, 0
	    else:
	        child_coords = [get_xy(c) for c in clade.clades]
	        xs = [c[0] for c in child_coords if not np.isnan(c[0])]
	        ys = [c[1] for c in child_coords]
	        x = np.mean(xs) if xs else np.nan
	        y = max(ys) + 1
	        clade._plot_children = child_coords
	        clade._plot_xy = (x, y)
	        return x, y
	
	def draw_clade(clade, ax):
	    if clade.is_terminal():
	        return
	    x, y = clade._plot_xy
	    child_xs = [c[0] for c in clade._plot_children]
	    child_ys = [c[1] for c in clade._plot_children]
	    ax.plot([min(child_xs), max(child_xs)], [y, y], color="black", lw=1)
	    for cx, cy in zip(child_xs, child_ys):
	        ax.plot([cx, cx], [cy, y], color="black", lw=1)
	    for c in clade.clades:
	        draw_clade(c, ax)
	
	get_xy(tree.root)
	
	### --- Schriftgrößen / Layout --- ###
	FONT_AXIS_LABEL = 28       # "Species" und "RefMember"
	FONT_TICK_X = 20            # Speziesnamen auf x-Achse
	FONT_TICK_Y = 9            # RefMember-Namen auf y-Achse
	FONT_LEGEND = 20           # Familiennamen in der Legende
	FONT_LEGEND_TITLE = 22     # Titel "Family"
	FONT_TITLE = 30            # Gesamttitel
	
	###  ---  Plot: Dendrogramm oben, Heatmap unten  ---  ###
	fig = plt.figure(figsize=(len(matrix.columns) * 0.75 + 6, len(matrix) * 0.35 + 8))
	gs = gridspec.GridSpec(2, 1, height_ratios=[4, len(matrix) * 0.35 + 4], hspace=0.02)
	
	ax_tree = fig.add_subplot(gs[0])
	draw_clade(tree.root, ax_tree)
	ax_tree.set_xlim(0, len(species_order_present))
	ax_tree.axis("off")
	
	ax_heat = fig.add_subplot(gs[1])
	sns.heatmap(
	    matrix, cmap=cmap, vmin=-0.5, vmax=len(all_families) - 0.5,
	    mask=matrix.isna(), cbar=False,
	    linewidths=0.5, linecolor="gray", ax=ax_heat
	)
	ax_heat.set_facecolor("white")
	ax_heat.set_xlabel("Species", fontsize=FONT_AXIS_LABEL, labelpad=15)
	ax_heat.set_ylabel("RefMember", fontsize=FONT_AXIS_LABEL, labelpad=15)
	
	ax_heat.tick_params(axis='x', which='major', labelsize=FONT_TICK_X, pad=3)
	ax_heat.tick_params(axis='y', which='major', labelsize=FONT_TICK_Y, pad=3)
	
	ax_heat.set_yticklabels(ax_heat.get_yticklabels(), rotation=0)
	
	families_in_plot = sorted(set(family_map[s] for s in species_order_present if s in family_map))
	legend_handles = [mpatches.Patch(color=family_color[f], label=f) for f in families_in_plot]
	
	'''
	# Legend vertical
	ax_heat.legend(
	    handles=legend_handles, title="Species family", loc='lower center', ncol=8, columnspacing=1.2, handletextpad=0.5,
	    bbox_to_anchor=(1.02, 1.1), borderaxespad=0, fontsize=FONT_LEGEND, title_fontsize=FONT_LEGEND_TITLE, frameon=True
	)
	
	'''
	# Legend horizontal
	ax_heat.legend(
	    handles=legend_handles, title=f"Species {CLASSIFICATION_TYPE}", loc='upper left',
	    bbox_to_anchor=(1.02, 1), borderaxespad=0, fontsize=FONT_LEGEND, title_fontsize=FONT_LEGEND_TITLE, frameon=True
	)
	
	
	fig.suptitle("Orthologue Presence/Absence plot\n(sorted phylogenetically)", fontsize=FONT_TITLE, fontweight='bold', y=0.995)
	plt.savefig(OUTPUT_FILE, format='svg', bbox_inches='tight', pad_inches=0.2)
	plt.close()
	
	print(f"Heatmap saced as {OUTPUT_FILE}")
	
if '--species_tree' in sys.argv and '--csv_heatmap' in sys.argv and '--species_classification_list' in sys.argv and '--out_file' in sys.argv and '--classification_type' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )