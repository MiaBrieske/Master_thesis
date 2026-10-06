# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 09:40:05 2026

@author: miabr
"""

__usage__ = """ python3 find_TGs_matching_TFs_final.py
					--TFLink_file <TFLink file -> TFLink_Mus_musculus_interactions_All_simpleFormat_v1.0_altered.tsv -> AFTER filter_evidence_pairs_final.py>
					--TF_list <TF LIST FROM find_mybs_final.py>
					--out_file <OUT FILE>
					"""

import sys
	
def main(arguments):
	
	TG_file = arguments[arguments.index('--TG_file') + 1]
	TF_list = arguments[arguments.index('--TF_list') + 1]
	out = arguments[arguments.index('--out_file') + 1]
	
	with open(TF_list, 'r') as f:
		TF_list = []
		f.readline()
		line = f.readline()
		while line:
			line = line.strip()
			TF_list.append(line)
			line = f.readline()
	
	with open(TG_file, 'r') as f, open(out, 'w') as o:
		f.readline()
		o.write("Transcription factor\tTranscript ID\tTarget gene\n")
		line = f.readline()
		while line:
			line = line.strip().split('\t')
			if len(line) < 5:
				print(line)
			else:
				TG_ID = line[3]
				TF = line[4]
				TG = line[5]
				if TF in TF_list:
					o.write(f"{TF}\t{TG_ID}\t{TG}\n")
			line = f.readline()

if '--TG_file' in sys.argv and '--TF_list' in sys.argv and '--out_file' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )