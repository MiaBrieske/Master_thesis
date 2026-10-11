# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 08:31:58 2026

@author: miabr
"""

"""
This script was developed with the help of AI tools (Claude). 
All sections of code were carefully reviewed, tested, and scientifically validated.
"""

"""
Builds one grouped bar plot (TG_id on the x-axis, one bar per TF_id) per
input file, from
files with the columns
	TF_ref  TF_id  TG_ref  TG_id  spearman_r  p_value  p_value_adj  n_samples

Short version of make_heatmaps_spearman_3.py: only the TF/TG correlation
filter (--strong-corr-cutoff / --variable-corr-cutoff, at least one of the
two is required), the pair-wise variability filter (--pair-variable-cutoff,
optional) and the top-N-TFs-per-species feature (--top-n-tfs, optional) can
be configured from the command line. Everything else uses fixed, reasonable
defaults (see the constants below).

TF_id are grouped by their TF_ref (same TF_ref stand next to each other),
TG_id likewise by TG_ref. A separator line is drawn between groups (also
inside the bar plot) and the respective ref name is written next to the axis.
The TG order is additionally re-sorted by hierarchical clustering of the TG
profiles.
"""

__usage__ = """ python3 bar_plots.py
					--files <COMMA-SEPARATED LIST OF TSV/CSV FILES>
					--input-dir <OPTIONAL: FOLDER TO READ ALL MATCHING FILES FROM>
					--pattern <OPTIONAL: FILE PATTERN FOR --input-dir, DEFAULT: *.tsv>
					--outdir <OPTIONAL: OUTPUT FOLDER, DEFAULT: ./barplots>
					--strong-corr-cutoff <ONLY KEEP TF_ref AND TG_ref (EACH A GENE) IF
						mean(|spearman_r|) OVER ALL THEIR TRANSCRIPTS AND CORRELATION
						PARTNERS IS >= CUTOFF>
					--variable-corr-cutoff <ONLY KEEP TF_ref AND TG_ref (EACH A GENE) IF
						THE RANGE (max - min) OF THEIR spearman_r VALUES (OVER ALL
						TRANSCRIPTS) IS >= CUTOFF>
					--corr-filter-combine <OPTIONAL: 'or' (DEFAULT) OR 'and' -- IF BOTH
						CUTOFFS ABOVE ARE SET, A TF_ref/TG_ref MUST FULFILL ONE (or) OR
						BOTH (and) CRITERIA>
					--pair-variable-cutoff <OPTIONAL: ONLY KEEP A SINGLE TF_ref-TG_ref PAIR
						IF THE RANGE (max - min) OF ITS OWN spearman_r VALUES (OVER ALL
						TF_id x TG_id TRANSCRIPT COMBINATIONS OF *THAT PAIR ONLY*) IS
						>= CUTOFF. UNLIKE --variable-corr-cutoff, THIS DOES NOT LOOK AT
						THE WHOLE TF_ref/TG_ref (ACROSS ALL ITS PARTNERS), ONLY AT THIS
						ONE PAIR. APPLIED *AFTER* --strong-corr-cutoff/--variable-corr-cutoff,
						WHICH CAN SERVE AS A CHEAP PRE-FILTER TO SHRINK THE DATA FIRST.
						DEFAULT: FEATURE OFF>
					--top-n-tfs <OPTIONAL: FROM THE FILTERED DATA OF ALL FILES, ONLY KEEP
						THE top-n-tfs TF_ref THAT OCCUR IN THE MOST SPECIES (= FILES);
						FOR THESE, ONE EXTRA BAR PLOT PER TF AND SPECIES IS CREATED IN
						<outdir>/top_tfs. DEFAULT: FEATURE OFF>
					--tg-batch-size <OPTIONAL: MAX AMOUNT OF TG_id (BARS ON THE X-AXIS)
						PER PLOT -- IF EXCEEDED, THE PLOT IS SPLIT INTO SEVERAL FILES
						(..._tgbatch1-of-N.png etc). ALSO APPLIED TO THE --top-n-tfs
						PLOTS. DEFAULT: NO LIMIT>
					--tf-batch-size <OPTIONAL: MAX AMOUNT OF TF_id (BARS PER TG_id GROUP)
						PER PLOT -- IF EXCEEDED, THE PLOT IS SPLIT INTO SEVERAL FILES
						(..._tfbatch1-of-N.png etc). DEFAULT: NO LIMIT>

					Note: at least --files OR --input-dir must be given.
					Note: at least one of --strong-corr-cutoff / --variable-corr-cutoff
					must be given.

					Examples:
					python3 make_heatmaps_spearman_short.py --files file1.tsv --strong-corr-cutoff 0.5
					python3 make_heatmaps_spearman_short.py --input-dir ./data --variable-corr-cutoff 1.2
					python3 make_heatmaps_spearman_short.py --files f1.tsv,f2.tsv,f3.tsv --strong-corr-cutoff 0.5 --top-n-tfs 10
					python3 make_heatmaps_spearman_short.py --files f1.tsv --strong-corr-cutoff 0.1 --pair-variable-cutoff 0.8
					"""

### FIXED SETTINGS (not configurable from the command line in this short version) ###
SEP = '\t'
VALUE_COL = 'spearman_r'
P_COL = 'p_value_adj'          # column used to restrict the correlation filter's
ALPHA = 0.05                   # statistics to significant rows; falls back to 'p_value'
VMIN, VMAX = -1.0, 1.0
DPI = 200
GROUP_FONTSIZE = 9
ANNOT = False

### IMPORTS ###
import colorsys
import glob
import os
import sys
import textwrap
from functools import partial
from multiprocessing import Pool
from types import SimpleNamespace

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # runs without a GUI
import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy
import scipy.cluster.hierarchy as sch


KNOWN_FLAGS = {
	'--files', '--input-dir', '--pattern', '--outdir',
	'--strong-corr-cutoff', '--variable-corr-cutoff', '--corr-filter-combine',
	'--pair-variable-cutoff', '--top-n-tfs', '--tf-batch-size', '--tg-batch-size',
}


def check_for_unknown_flags(arguments):
	"""
	Rejects any '--...' token that isn't a recognized flag, INSTEAD OF
	silently ignoring it (which would otherwise make a typo like
	'--corr_filter_combine' instead of '--corr-filter-combine' fall back to
	that flag's default value without any warning). arguments[0] is the
	script path and is skipped; every other '--...'-looking token must be a
	known flag.
	"""
	unknown = [a for a in arguments[1:] if a.startswith('--') and a not in KNOWN_FLAGS]
	if unknown:
		sys.exit(
			f"Unknown flag(s): {', '.join(unknown)}\n"
			f"Known flags are: {', '.join(sorted(KNOWN_FLAGS))}"
		)


def parse_arguments(arguments):
	"""
	Reads all parameters manually from sys.argv (arguments.index('--flag') + 1).
	Returns a SimpleNamespace object that the rest of the code can use just
	like an argparse namespace, via args.attribute_name.
	"""
	check_for_unknown_flags(arguments)

	if '--files' in arguments:
		files = arguments[arguments.index('--files') + 1].split(',')
	else:
		files = []

	if '--input-dir' in arguments:
		input_dir = arguments[arguments.index('--input-dir') + 1]
	else:
		input_dir = None

	if '--pattern' in arguments:
		pattern = arguments[arguments.index('--pattern') + 1]
	else:
		pattern = '*.tsv'

	if '--outdir' in arguments:
		outdir = arguments[arguments.index('--outdir') + 1]
	else:
		outdir = './barplots'

	if '--strong-corr-cutoff' in arguments:
		strong_corr_cutoff = float(arguments[arguments.index('--strong-corr-cutoff') + 1])
	else:
		strong_corr_cutoff = None

	if '--variable-corr-cutoff' in arguments:
		variable_corr_cutoff = float(arguments[arguments.index('--variable-corr-cutoff') + 1])
	else:
		variable_corr_cutoff = None

	if strong_corr_cutoff is None and variable_corr_cutoff is None:
		sys.exit("At least one of --strong-corr-cutoff or --variable-corr-cutoff has to be given.")
	if strong_corr_cutoff is not None and strong_corr_cutoff < 0:
		sys.exit("--strong-corr-cutoff has to be >= 0")
	if variable_corr_cutoff is not None and variable_corr_cutoff < 0:
		sys.exit("--variable-corr-cutoff has to be >= 0")

	if '--corr-filter-combine' in arguments:
		corr_filter_combine = arguments[arguments.index('--corr-filter-combine') + 1]
	else:
		corr_filter_combine = 'or'
	if corr_filter_combine not in ('or', 'and'):
		sys.exit(f"--corr-filter-combine has to be 'or' or 'and', not '{corr_filter_combine}'")

	if '--pair-variable-cutoff' in arguments:
		pair_variable_cutoff = float(arguments[arguments.index('--pair-variable-cutoff') + 1])
	else:
		pair_variable_cutoff = None
	if pair_variable_cutoff is not None and pair_variable_cutoff < 0:
		sys.exit("--pair-variable-cutoff has to be >= 0")

	if '--top-n-tfs' in arguments:
		top_n_tfs = int(arguments[arguments.index('--top-n-tfs') + 1])
	else:
		top_n_tfs = None
	if top_n_tfs is not None and top_n_tfs < 1:
		sys.exit("--top-n-tfs has to be >= 1")

	if '--tf-batch-size' in arguments:
		tf_batch_size = int(arguments[arguments.index('--tf-batch-size') + 1])
	else:
		tf_batch_size = None
	if tf_batch_size is not None and tf_batch_size < 1:
		sys.exit("--tf-batch-size has to be >= 1")

	if '--tg-batch-size' in arguments:
		tg_batch_size = int(arguments[arguments.index('--tg-batch-size') + 1])
	else:
		tg_batch_size = None
	if tg_batch_size is not None and tg_batch_size < 1:
		sys.exit("--tg-batch-size has to be >= 1")

	return SimpleNamespace(
		files=files,
		input_dir=input_dir,
		pattern=pattern,
		outdir=outdir,
		strong_corr_cutoff=strong_corr_cutoff,
		variable_corr_cutoff=variable_corr_cutoff,
		corr_filter_combine=corr_filter_combine,
		pair_variable_cutoff=pair_variable_cutoff,
		top_n_tfs=top_n_tfs,
		tf_batch_size=tf_batch_size,
		tg_batch_size=tg_batch_size,
	)


def collect_files(args):
	files = list(args.files)
	if args.input_dir:
		files += glob.glob(os.path.join(args.input_dir, args.pattern))
	seen = set()
	uniq = []
	for f in files:
		if f not in seen:
			uniq.append(f)
			seen.add(f)
	return uniq


def build_ordered_groups(df, ref_col, id_col, id_order=None):
	"""
	Returns a list of (ref, ids) tuples -- ids each sorted (or, if id_order is
	given, in that order), refs sorted alphabetically.
	"""
	sub = df[[ref_col, id_col]].drop_duplicates()

	if id_order is None:
		sub = sub.sort_values([ref_col, id_col])
	else:
		order_index = {v: i for i, v in enumerate(id_order)}
		sub = sub.assign(_order=sub[id_col].map(order_index)).sort_values([ref_col, '_order'])

	groups = []
	for ref, ids in sub.groupby(ref_col, sort=True)[id_col]:
		groups.append((ref, ids.tolist()))
	return groups


def build_groups_from_id_order(df, ref_col, id_col, id_order):
	"""
	Builds (ref, ids) groups from a given id_order (e.g. after clustering).
	Equal ref values are only merged into one group if they follow each other
	directly in id_order.
	"""
	sub = df[[ref_col, id_col]].drop_duplicates().set_index(id_col).loc[list(id_order)]

	groups = []
	current_ref = None
	current_ids = []
	for tid, ref in zip(sub.index, sub[ref_col]):
		if ref != current_ref:
			if current_ids:
				groups.append((current_ref, current_ids))
			current_ref = ref
			current_ids = []
		current_ids.append(tid)
	if current_ids:
		groups.append((current_ref, current_ids))
	return groups


def batch_groups(groups, batch_size):
	"""
	Splits a list of (ref, ids) groups into batches so that every batch
	contains at most batch_size ids in total. If a single group is itself
	bigger than batch_size, it is split into its own batch(es) rather than
	being torn apart across a shared batch with other groups.
	"""
	if batch_size is None:
		return [groups]

	batches = []
	current_batch = []
	current_size = 0

	def flush():
		nonlocal current_batch, current_size
		if current_batch:
			batches.append(current_batch)
			current_batch = []
			current_size = 0

	for ref, ids in groups:
		if len(ids) > batch_size:
			flush()
			for i in range(0, len(ids), batch_size):
				batches.append([(ref, ids[i:i + batch_size])])
			continue
		if current_batch and current_size + len(ids) > batch_size:
			flush()
		current_batch.append((ref, ids))
		current_size += len(ids)

	flush()
	return batches if batches else [[]]


def corr_keep_set(df, ref_col, value_col, strong_cutoff, variable_cutoff, combine):
	"""
	Returns the set of ref_col values (TF_ref OR TG_ref) that (depending on
	which cutoffs are set) show a strong and/or variable correlation.
	"""
	stats = df.groupby(ref_col)[value_col].agg(
		mean_abs=lambda s: s.abs().mean(),
		range_=lambda s: s.max() - s.min(),
		n_values=lambda s: s.count(),
	)

	masks = []
	if strong_cutoff is not None:
		masks.append(stats['mean_abs'] >= strong_cutoff)
	if variable_cutoff is not None:
		masks.append(stats['range_'] >= variable_cutoff)

	keep_mask = masks[0]
	for m in masks[1:]:
		if combine == 'or':
			keep_mask = keep_mask | m
		else:
			keep_mask = keep_mask & m
	keep_mask = keep_mask.fillna(False)

	return set(keep_mask[keep_mask].index)


def filter_tf_tg_by_correlation(df, value_col, strong_cutoff, variable_cutoff, combine, filepath, p_col, alpha):
	"""
	Keeps only whole TF_ref AND whole TG_ref (each = a gene) that (depending
	on which cutoffs are set) show a strong and/or variable correlation over
	ALL their transcripts and correlation partners in THIS file combined --
	see corr_keep_set. TF_ref and TG_ref filters are computed independently of
	each other; a row is only removed if NEITHER its TF_ref NOR its TG_ref
	fulfills the criterion.
	"""
	if p_col is not None:
		df_for_stats = df[df[p_col] <= alpha]
	else:
		df_for_stats = df

	keep_tf_refs = corr_keep_set(df_for_stats, 'TF_ref', value_col, strong_cutoff, variable_cutoff, combine)
	keep_tg_refs = corr_keep_set(df_for_stats, 'TG_ref', value_col, strong_cutoff, variable_cutoff, combine)

	before_tf_refs, before_tf_ids = df['TF_ref'].nunique(), df['TF_id'].nunique()
	before_tg_refs, before_tg_ids = df['TG_ref'].nunique(), df['TG_id'].nunique()

	df = df[df['TF_ref'].isin(keep_tf_refs) | df['TG_ref'].isin(keep_tg_refs)]

	after_tf_refs, after_tf_ids = df['TF_ref'].nunique(), df['TF_id'].nunique()
	after_tg_refs, after_tg_ids = df['TG_ref'].nunique(), df['TG_id'].nunique()

	if p_col is not None:
		sig_note = ' (calculated from significant values only)'
	else:
		sig_note = ' (warning: no p-value column -- all values considered)'

	print(f"INFO: {filepath}: TF correlation filter{sig_note} -> {after_tf_refs}/{before_tf_refs} TF_ref "
			f"kept ({after_tf_ids}/{before_tf_ids} TF_id transcripts)")
	print(f"INFO: {filepath}: TG correlation filter{sig_note} -> {after_tg_refs}/{before_tg_refs} TG_ref "
			f"kept ({after_tg_ids}/{before_tg_ids} TG_id transcripts)")
	return df


def filter_pairs_by_variability(df, value_col, cutoff, filepath):
	"""
	Keeps only individual TF_ref-TG_ref PAIRS whose own range (max - min) of
	spearman_r values -- computed ONLY from that pair's own TF_id x TG_id
	transcript combinations, not from the whole TF_ref or TG_ref across all
	their OTHER partners -- is >= cutoff. This is a stricter, pair-local
	variant of --variable-corr-cutoff: a TF_ref/TG_ref can pass the
	whole-gene filter (e.g. because it is highly variable against some OTHER
	partner) and still be dropped here if THIS SPECIFIC pair is flat.
	Intended to run after --strong-corr-cutoff/--variable-corr-cutoff, which
	can pre-shrink the data (cheaper, whole-gene check) before this more
	expensive per-pair check runs on what's left.
	"""
	if cutoff is None:
		return df

	pair_keys = pd.Series(list(zip(df['TF_ref'], df['TG_ref'])), index=df.index)

	stats = df.groupby(['TF_ref', 'TG_ref'])[value_col].agg(
		range_=lambda s: s.max() - s.min(),
		n_values='count',
	)
	keep_pairs = set(stats[stats['range_'] >= cutoff].index)

	before_pairs = df[['TF_ref', 'TG_ref']].drop_duplicates().shape[0]
	df = df[pair_keys.isin(keep_pairs)]
	after_pairs = df[['TF_ref', 'TG_ref']].drop_duplicates().shape[0] if not df.empty else 0

	print(f"INFO: {filepath}: pair-wise variability filter -> {after_pairs}/{before_pairs} "
			f"TF_ref-TG_ref pairs kept (range(r) of that pair alone >= {cutoff})")
	return df


def _safe_linkage_order(data, method='ward', metric='euclidean'):
	"""
	data: DataFrame, rows = items to order, columns = features.
	Returns the row labels re-ordered by hierarchical clustering. Handles
	NaNs (row-wise mean imputation, then 0.0 for still-NaN rows -- otherwise
	scipy's linkage fails with "condensed distance matrix must contain only
	finite values") and degenerate cases with <= 1 row.
	"""
	if data.shape[0] <= 1:
		return list(data.index)

	data = data.fillna(data.mean(axis=1), axis=0)
	data = data.fillna(0.0)

	Y = sch.linkage(data.values, method=method, metric=metric)
	leaves = sch.leaves_list(Y)

	return [data.index[i] for i in leaves]


def cluster_order_within_groups(pivot, groups, method='ward', metric='euclidean'):
	"""
	Orders the TG_ids (pivot columns) via hierarchical clustering while
	keeping every TG_ref group's TG_ids contiguous, so TG_ids sharing a
	TG_ref always cluster together:
	  1. the TG_ref groups themselves are ordered by clustering their mean
	     profiles (mean correlation profile across TF_ids);
	  2. within each TG_ref group, the individual TG_ids are ordered by
	     clustering their own profiles.
	Returns (flat_tg_id_order, ordered_groups) where ordered_groups is the
	list of (TG_ref, [TG_id, ...]) tuples in the resulting order.
	"""
	data = pivot.T  # index = TG_id, columns = TF_id

	if len(groups) > 1:
		group_profiles = pd.DataFrame(
			{ref: data.loc[ids].mean(axis=0) for ref, ids in groups}
		).T
		ordered_refs = _safe_linkage_order(group_profiles, method=method, metric=metric)
		groups_by_ref = dict(groups)
		ordered_groups_in = [(ref, groups_by_ref[ref]) for ref in ordered_refs]
	else:
		ordered_groups_in = groups

	final_order = []
	ordered_groups = []
	for ref, ids in ordered_groups_in:
		ordered_ids = _safe_linkage_order(data.loc[ids], method=method, metric=metric)
		final_order.extend(ordered_ids)
		ordered_groups.append((ref, ordered_ids))

	return final_order, ordered_groups


def wrap_ref_label(text, span, chars_per_unit, min_chars=6, max_lines=3):
	"""
	Wraps a ref label (TF_ref / TG_ref) if it is wider than its group's span
	allows, so neighbouring group labels don't overlap each other. span is
	the number of ids (TF_id / TG_id) in the group.
	"""
	max_chars = max(min_chars, int(round(span * chars_per_unit)))
	wrapped = textwrap.fill(str(text), width=max_chars)
	lines = wrapped.split('\n')
	if len(lines) > max_lines:
		lines = lines[:max_lines]
		lines[-1] = lines[-1].rstrip() + '…'
		wrapped = '\n'.join(lines)
	return wrapped, len(lines)


def add_group_labels(ax, fig, tf_groups, tg_groups, fontsize, group_lw=1.5):
	"""
	Draws separator lines and ref labels (TF_ref / TG_ref) outside the
	plot axes. The gap to the axes is computed DYNAMICALLY from the
	actual rendered size of the TF_id/TG_id tick labels (measured via the
	renderer), so long ids never overlap with the ref group labels. Long ref
	labels are wrapped to fit their group's span, and neighbouring groups are
	staggered (alternating extra offset, sized to fit the wrapped text) so
	the ref labels don't overlap each other either.
	"""
	fig.canvas.draw()  # make sure tick label positions are up to date
	renderer = fig.canvas.get_renderer()

	ax_bbox_px = ax.get_window_extent(renderer=renderer)
	pt_per_px = 72.0 / fig.dpi

	# --- how far do the existing TG_id tick labels extend below the axes? ---
	xtick_labels = [t for t in ax.get_xticklabels() if t.get_text()]
	if xtick_labels:
		min_y0 = min(t.get_window_extent(renderer=renderer).y0 for t in xtick_labels)
		x_extent_pt = max(0.0, ax_bbox_px.y0 - min_y0) * pt_per_px
	else:
		x_extent_pt = 0.0

	# --- how far do the existing TF_id tick labels extend left of the axes? ---
	ytick_labels = [t for t in ax.get_yticklabels() if t.get_text()]
	if ytick_labels:
		min_x0 = min(t.get_window_extent(renderer=renderer).x0 for t in ytick_labels)
		y_extent_pt = max(0.0, ax_bbox_px.x0 - min_x0) * pt_per_px
	else:
		y_extent_pt = 0.0

	line_gap = 6                   # gap between tick labels and separator line (points)
	text_gap = 4                   # gap between separator line and ref text (points)
	line_height = fontsize * 1.25  # approx. line height of the (possibly wrapped) ref text

	# wrap TG_ref / TF_ref labels so they fit inside their group's span
	tg_wrapped = [(ref, ids, *wrap_ref_label(ref, len(ids), chars_per_unit=2.3)) for ref, ids in tg_groups]
	tf_wrapped = [(ref, ids, *wrap_ref_label(ref, len(ids), chars_per_unit=3.0)) for ref, ids in tf_groups]

	max_lines_tg = max((n for *_, n in tg_wrapped), default=1)
	max_lines_tf = max((n for *_, n in tf_wrapped), default=1)

	# the stagger offset between the two alternating label levels must be big
	# enough to fit the tallest (possibly multi-line) label, or the two
	# levels would overlap each other
	stagger_step_tg = text_gap + line_height * max_lines_tg + 8
	stagger_step_tf = text_gap + line_height * max_lines_tf + 8

	tg_label_gap = x_extent_pt + line_gap
	tf_label_gap = y_extent_pt + line_gap

	# ---- TG_ref: line + label BELOW the TG_id ticklabels ----
	x_pos = 0
	for i, (ref, ids, label, n_lines) in enumerate(tg_wrapped):
		start = x_pos - 0.5
		end = x_pos + len(ids) - 0.5
		center = (start + end) / 2

		if x_pos != 0:
			ax.axvline(start, color='black', linewidth=group_lw)

		dy_line = -(tg_label_gap + stagger_step_tg * (i % 2))
		trans_line = offset_copy(ax.get_xaxis_transform(), fig=fig, y=dy_line, units='points')
		ax.plot([start, end], [0, 0], transform=trans_line, color='black', linewidth=1.2, clip_on=False)

		trans_text = offset_copy(ax.get_xaxis_transform(), fig=fig, y=dy_line - text_gap, units='points')
		ax.text(center, 0, label, transform=trans_text, ha='center', va='top',
				fontsize=fontsize, fontweight='bold', fontstyle='italic', linespacing=1.25, clip_on=False)

		x_pos += len(ids)

	# ---- TF_ref: line + label LEFT of the TF_id ticklabels, staggered ----
	y_pos = 0
	for i, (ref, ids, label, n_lines) in enumerate(tf_wrapped):
		start = y_pos - 0.5
		end = y_pos + len(ids) - 0.5
		center = (start + end) / 2

		if y_pos != 0:
			ax.axhline(start, color='black', linewidth=group_lw)

		dx = -(tf_label_gap + stagger_step_tf * (i % 2))
		trans_line = offset_copy(ax.get_yaxis_transform(), fig=fig, x=dx, units='points')
		ax.plot([0, 0], [start, end], transform=trans_line, color='black', linewidth=1.2, clip_on=False)

		trans_text = offset_copy(ax.get_yaxis_transform(), fig=fig, x=dx - text_gap, units='points')
		ax.text(0, center, label, transform=trans_text, ha='right', va='center',
				fontsize=fontsize, fontweight='bold', linespacing=1.25, clip_on=False)

		y_pos += len(ids)


def build_tf_colors(tf_groups, hue_start=0.0):
	"""
	Assigns one color per TF_id such that:
	  - each TF_ref GROUP gets its own clearly distinct base hue (evenly
	    spaced around the color wheel, e.g. group 1 = green family, group 2 =
	    red family, ...), so different genes are never confusable;
	  - each TF_id WITHIN a group is a different lightness shade of that
	    SAME hue (e.g. 4 transcripts of one TF_ref = 4 shades of green), so
	    transcripts of the same gene read as visually related.
	Returns a dict {tf_id: (r, g, b)}.
	"""
	n_groups = len(tf_groups)
	colors = {}
	for gi, (ref, ids) in enumerate(tf_groups):
		hue = (hue_start + gi / max(n_groups, 1)) % 1.0
		n_ids = len(ids)
		for ii, tf_id in enumerate(ids):
			# darker shade for the first transcript, lighter for the last;
			# a single-transcript group gets one medium-toned shade
			lightness = 0.5 if n_ids == 1 else 0.32 + 0.36 * (ii / (n_ids - 1))
			colors[tf_id] = colorsys.hls_to_rgb(hue, lightness, 0.68)
	return colors


def draw_barplot(pivot, tf_groups, tg_groups, title, group_lw=1.5):
	"""
	Grouped bar plot: one x position per TG_id (grouped/labeled by TG_ref,
	same as the heatmap columns), with one bar per TF_id at that position
	(colored by TF_ref group -- see build_tf_colors -- with a legend). Bar
	height = the value in `pivot` (spearman_r). Replaces the former TF_id x
	TG_id heatmap.
	"""
	# (tf_id, tf_ref) pairs, flattened in the same order as tf_groups, so the
	# legend can show both -- e.g. "TFA.1 (TFA)" -- instead of just TF_id
	tf_id_ref_pairs = [(tid, ref) for ref, ids in tf_groups for tid in ids]
	tf_ids = [tid for tid, ref in tf_id_ref_pairs]
	tg_ids = [tid for ref, ids in tg_groups for tid in ids]

	n_tf = max(len(tf_ids), 1)
	n_tg = max(len(tg_ids), 1)

	fig_w = max(7, 0.55 * n_tg + 0.15 * n_tf + 3)
	fig_h = 6

	fig, ax = plt.subplots(figsize=(fig_w, fig_h))

	tf_colors = build_tf_colors(tf_groups)

	bar_width = min(0.8 / n_tf, 0.8)
	x = np.arange(n_tg)

	for j, (tf_id, tf_ref) in enumerate(tf_id_ref_pairs):
		offset = (j - (n_tf - 1) / 2) * bar_width
		heights = pivot.loc[tf_id, tg_ids].values.astype(float)
		# label shows TF_id AND TF_ref, e.g. "TFA.1 (TFA)"; if TF_id and
		# TF_ref are identical (single-transcript case) it's shown just once
		label = tf_id if tf_id == tf_ref else f"{tf_id} ({tf_ref})"
		ax.bar(x + offset, heights, width=bar_width, color=tf_colors[tf_id], label=label)

	ax.axhline(0, color='black', linewidth=0.8)
	ax.set_xticks(x)
	ax.set_xticklabels(tg_ids, rotation=90, fontsize=8)
	ax.set_xlim(-0.5, n_tg - 0.5)
	ax.set_ylim(VMIN, VMAX)
	ax.set_ylabel(VALUE_COL)

	for spine in ('top', 'right'):
		ax.spines[spine].set_visible(False)

	# TG_ref separators + labels below the TG_id ticks (reuses the same
	# dynamic-gap/wrap logic as the heatmap; tf_groups=[] skips the
	# (no longer needed) left-hand TF_ref labels)
	add_group_labels(ax, fig, [], tg_groups, fontsize=GROUP_FONTSIZE, group_lw=group_lw)

	ax.set_title(title, pad=20)

	ncol = 1 if n_tf <= 25 else 2
	ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0), fontsize=7, ncol=ncol, title='TF_id (TF_ref)')

	fig.tight_layout()

	return fig, ax


def read_and_filter_file(filepath, args):
	"""
	Reads one file and applies the TF/TG correlation filter
	(--strong-corr-cutoff / --variable-corr-cutoff). Returns (species, df) or
	None if the file has to be skipped.
	"""
	header_cols = set(pd.read_csv(filepath, sep=SEP, nrows=0).columns)
	required_cols = {'TF_ref', 'TF_id', 'TG_ref', 'TG_id', VALUE_COL}
	missing = required_cols - header_cols
	if missing:
		print(f"WARNING: {filepath} skipped - missing columns {missing}")
		return None

	needed_cols = set(required_cols)
	needed_cols.add(P_COL)
	needed_cols.add('p_value')

	df = pd.read_csv(filepath, sep=SEP, usecols=lambda c: c in needed_cols and c in header_cols)

	p_col_resolved = P_COL
	if p_col_resolved not in df.columns:
		if 'p_value' in df.columns:
			print(f"WARNING: {filepath}: column '{p_col_resolved}' not found -> 'p_value' used instead!")
			p_col_resolved = 'p_value'
		else:
			print(f"WARNING: {filepath}: neither '{P_COL}' nor 'p_value' -> ignoring significance restriction for this file!")
			p_col_resolved = None

	df = filter_tf_tg_by_correlation(
		df, VALUE_COL, args.strong_corr_cutoff, args.variable_corr_cutoff,
		args.corr_filter_combine, filepath, p_col_resolved, ALPHA
	)
	if df.empty:
		print(f"WARNING: {filepath}: no data left after TF/TG correlation filter -> skipped")
		return None

	# pair-wise variability filter runs AFTER the (cheaper) whole-gene filter
	# above, which shrinks the data first
	df = filter_pairs_by_variability(df, VALUE_COL, args.pair_variable_cutoff, filepath)
	if df.empty:
		print(f"WARNING: {filepath}: no data left after pair-wise variability filter -> skipped")
		return None

	species = os.path.basename(filepath).replace('_spearman_corr.txt', '')
	return species, df


def make_barplot(filepath, args):
	loaded = read_and_filter_file(filepath, args)
	if loaded is None:
		return
	species, df = loaded

	tf_groups_all = build_ordered_groups(df, 'TF_ref', 'TF_id')
	tg_groups_all = build_ordered_groups(df, 'TG_ref', 'TG_id')

	# split into batches (alphabetical order) BEFORE clustering, so clustering
	# then only has to re-sort the TGs within each batch
	tf_batches = batch_groups(tf_groups_all, args.tf_batch_size)
	tg_batches = batch_groups(tg_groups_all, args.tg_batch_size)

	n_tf_batches = len(tf_batches)
	n_tg_batches = len(tg_batches)
	multi_batch = n_tf_batches > 1 or n_tg_batches > 1

	base = os.path.splitext(os.path.basename(filepath))[0]
	os.makedirs(args.outdir, exist_ok=True)

	for tf_i, tf_batch in enumerate(tf_batches, start=1):
		tf_ids = [tid for ref, ids in tf_batch for tid in ids]

		for tg_i, tg_batch in enumerate(tg_batches, start=1):
			tg_ids = [tid for ref, ids in tg_batch for tid in ids]

			sub = df[df['TF_id'].isin(tf_ids) & df['TG_id'].isin(tg_ids)]
			pivot = sub.pivot_table(index='TF_id', columns='TG_id', values=VALUE_COL, aggfunc='mean')
			pivot = pivot.reindex(index=tf_ids, columns=tg_ids)

			# re-sort the TGs within this batch by hierarchical clustering, keeping
			# TG_ids that share a TG_ref contiguous (grouped by design, not just by
			# a shared column position)
			tg_clustered_order, tg_batch_clustered = cluster_order_within_groups(pivot, tg_batch)
			pivot = pivot.reindex(index=tf_ids, columns=tg_clustered_order)

			title = f"{species}\nTF x TG ({VALUE_COL})"
			parts = []
			if args.strong_corr_cutoff is not None:
				parts.append(f"mean|r| >= {args.strong_corr_cutoff}")
			if args.variable_corr_cutoff is not None:
				parts.append(f"range(r) >= {args.variable_corr_cutoff}")
			joiner = ' OR ' if args.corr_filter_combine == 'or' else ' AND '
			title += f"\nTF/TG filter: {joiner.join(parts)}"
			if args.pair_variable_cutoff is not None:
				title += f"\npair-wise filter: range(r) of pair >= {args.pair_variable_cutoff}"
			if multi_batch:
				title += f"\n(batch: TF {tf_i}/{n_tf_batches}, TG {tg_i}/{n_tg_batches})"

			suffix = ''
			if n_tf_batches > 1:
				suffix += f'_tfbatch{tf_i}-of-{n_tf_batches}'
			if n_tg_batches > 1:
				suffix += f'_tgbatch{tg_i}-of-{n_tg_batches}'

			outpath = os.path.join(args.outdir, f"{base}_barplot{suffix}.png")

			fig, ax = draw_barplot(pivot, tf_batch, tg_batch_clustered, title)

			fig.savefig(outpath, dpi=DPI, bbox_inches='tight')
			plt.close(fig)
			print(f"Created: {outpath}")


def make_barplot_safe(filepath, args):
	try:
		make_barplot(filepath, args)
	except Exception as e:
		print(f"ERROR: {filepath}: {e}", file=sys.stderr)


def count_species_per_tf(species_dfs):
	"""
	species_dfs: list of (species, df) tuples (already correlation-filtered).
	Returns a pandas Series: TF_ref -> number of distinct species it occurs
	in, sorted descending.
	"""
	rows = []
	for species, df in species_dfs:
		for tf_ref in df['TF_ref'].unique():
			rows.append((tf_ref, species))
	pairs = pd.DataFrame(rows, columns=['TF_ref', 'species'])
	counts = pairs.groupby('TF_ref')['species'].nunique().sort_values(ascending=False)
	return counts


def keep_top_tfs(counts, top_n):
	"""
	Keeps the top_n TF_ref with the highest species count. Ties at the cutoff
	value are all kept, so the result can contain more than top_n entries.
	"""
	if len(counts) > top_n:
		cutoff = counts.iloc[top_n - 1]
		top_tfs = set(counts[counts >= cutoff].index)
	else:
		top_tfs = set(counts.index)

	print(f"Top {top_n} TFs (present in most species):")
	for tf_ref, count in counts.head(top_n).items():
		print(f"  {tf_ref}: {count} species")

	return top_tfs


def draw_tf_species_barplots(species_dfs, top_tfs, args):
	"""
	For every (species, TF_ref) combination where TF_ref is in top_tfs and
	present in that species' data, draws one bar plot: that TF (all its
	TF_id transcripts) against all its TG_id correlation partners in that
	species.
	"""
	out_dir = os.path.join(args.outdir, 'top_tfs')
	os.makedirs(out_dir, exist_ok=True)

	for species, df in species_dfs:
		species_tfs = sorted(t for t in df['TF_ref'].unique() if t in top_tfs)

		for tf_ref in species_tfs:
			tf_sub = df[df['TF_ref'] == tf_ref]
			tf_ids = sorted(tf_sub['TF_id'].unique())
			tf_groups = [(tf_ref, tf_ids)]

			tg_groups_all = build_ordered_groups(tf_sub, 'TG_ref', 'TG_id')
			tg_batches = batch_groups(tg_groups_all, args.tg_batch_size)
			n_tg_batches = len(tg_batches)

			for tg_i, tg_batch in enumerate(tg_batches, start=1):
				tg_ids = [tid for ref, ids in tg_batch for tid in ids]

				pivot = tf_sub.pivot_table(index='TF_id', columns='TG_id', values=VALUE_COL, aggfunc='mean')
				pivot = pivot.reindex(index=tf_ids, columns=tg_ids)

				title = f"{species} - {tf_ref}"
				if n_tg_batches > 1:
					title += f"\n(TG batch {tg_i}/{n_tg_batches})"
				fig, ax = draw_barplot(pivot, tf_groups, tg_batch, title)

				suffix = f'_tgbatch{tg_i}-of-{n_tg_batches}' if n_tg_batches > 1 else ''
				filename = f"top_tf_{species}_{tf_ref}_barplot{suffix}.png"
				outpath = os.path.join(out_dir, filename)

				fig.savefig(outpath, dpi=DPI, bbox_inches='tight')
				plt.close(fig)
				print(f"Created: {outpath}")


def run_top_tf_pipeline(files, args):
	"""
	Collects the correlation-filtered data of ALL input files, keeps only the
	--top-n-tfs TF_ref that occur in the most species (= files), and draws
	one bar plot per TF and species for those.
	"""
	species_dfs = []
	for filepath in files:
		loaded = read_and_filter_file(filepath, args)
		if loaded is None:
			continue
		species_dfs.append(loaded)

	if not species_dfs:
		print("WARNING: --top-n-tfs given, but no data available -> skipped")
		return

	counts = count_species_per_tf(species_dfs)
	top_tfs = keep_top_tfs(counts, args.top_n_tfs)
	draw_tf_species_barplots(species_dfs, top_tfs, args)


def main(arguments):
	args = parse_arguments(arguments)
	files = collect_files(args)
	if not files:
		sys.exit("No input files found. Use --files or --input-dir/--pattern.")

	n_workers = min(len(files), os.cpu_count() or 1)
	with Pool(n_workers) as pool:
		pool.map(partial(make_barplot_safe, args=args), files)

	if args.top_n_tfs is not None:
		run_top_tf_pipeline(files, args)

	print('Done.')


if '--files' in sys.argv or '--input-dir' in sys.argv:
	main(sys.argv)
else:
	sys.exit(__usage__)