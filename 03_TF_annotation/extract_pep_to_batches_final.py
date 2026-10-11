# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 13:38:29 2025

@author: miabr
"""

"""
This script was developed with the help of AI tools (Claude). 
All sections of code were carefully reviewed, tested, and scientifically validated.
"""

import os, sys
import glob
import tarfile
import shutil
import argparse
from datetime import datetime

__usage__ = """ python3 extract_pep_to_batches_final.py 
					--folder_path <PATH TO PLANT FOLDER>
					OPTIONAL:
					--batch_out <BATCH OUTPUT FOLDER>
					--batch_size <BATCH SIZE>
"""

def extract_pep_to_batches(plant_folder, batch_size=50, dest_base="myb_pep_batch"):
    """Extrahiert .pep.* Dateien aus tar.gz Archiven in ./plants/ und verteilt sie in Batches."""
    
    # Alle .tar.gz in ./plants/ finden
    pattern = os.path.join(plant_folder, "*.tar.gz")
    tar_files = glob.glob(pattern)
    
    print(f"Gefunden: {len(tar_files)} .tar.gz Archive in {plant_folder}")
    if len(tar_files) == 0:
        print(f"Kein .tar.gz in {plant_folder} gefunden!")
        print("Prüfe: ls {plant_folder}/*.tar.gz")
        return
    
    # In Batches aufteilen
    for i in range(0, len(tar_files), batch_size):
        batch_tars = tar_files[i:i + batch_size]
        batch_num = i // batch_size + 1
        batch_folder = f"{dest_base}_{batch_num}"
        os.makedirs(batch_folder, exist_ok=True)
        
        extracted = 0
        for tar_path in batch_tars:
            try:
                with tarfile.open(tar_path, "r:gz") as tar:
                    for member in tar.getmembers():
                        # Prüfe ob .pep. als 2. Teil im Dateinamen (z.B. species.pep.fasta)
                        name = member.name
                        parts = member.name.split(".")
                        print(f"  {name} → Teile: {parts}")
                        if "pep" in parts:
                            f = tar.extractfile(member)
                            pep_name = os.path.basename(member.name)
                            ziel_datei = os.path.join(batch_folder, pep_name)
                            
                            if f:
                                with open(ziel_datei, 'wb') as out_file:
                                    shutil.copyfileobj(f, out_file)
                                extracted += 1
                                print(f"  → Extrahiert: {pep_name}")
                                break  # Nur 1 .pep pro tar.gz
            except Exception as e:
                print(f"  Fehler bei {os.path.basename(tar_path)}: {e}")
                continue
        
        print(f"Batch {batch_num}: {len(batch_tars)} tar.gz → {extracted} .pep → {batch_folder}/")
    
    total_batches = len(tar_files) // batch_size + (1 if len(tar_files) % batch_size else 0)
    print(f"\nFertig! {total_batches} Batch-Ordner mit je {batch_size} .pep Dateien.")

def main(arguments):
    
    folder_path = arguments[arguments.index("--folder_path") + 1]
    if "--batch_out" in sys.argv:
        batch_out = arguments[arguments.index("--batch_out") + 1]
    else:
        batch_out = "pep_batch"
    
    if "--batch_size" in sys.argv:
       batch_size = int(arguments[arguments.index("--batch_size") + 1])
    else:
       batch_size = 50
    
    start = datetime.now()
    extract_pep_to_batches(folder_path, batch_size, batch_out)
    end = datetime.now()
    print(f"Runtime: {end - start}")

if '--folder_path' in sys.argv:
    main( sys.argv )
else:
    sys.exit( __usage__ )

