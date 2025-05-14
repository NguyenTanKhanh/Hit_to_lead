#!/usr/bin/env python3
import os
import subprocess
import csv
import re
import sys
import argparse

ascii_art = r"""
____   ____.__                _________                                   
\   \ /   /|__| ____ _____   /   _____/ ___________   ____   ____   ____  
 \   Y   / |  |/    \\__  \  \_____  \_/ ___\_  __ \_/ __ \_/ __ \ /    \ 
  \     /  |  |   |  \/ __ \_/        \  \___|  | \/\  ___/\  ___/|   |  \
   \___/   |__|___|  (____  /_______  /\___  >__|    \___  >\___  >___|  /
                   \/     \/        \/     \/            \/     \/     \/ 
"""
print(ascii_art)

parser = argparse.ArgumentParser(description="Multi-ligand docking using AutoDock Vina")
parser.add_argument('--receptor', required=True, help='Path to receptor PDBQT file')
parser.add_argument('--ligand_dir', required=True, help='Directory containing ligand PDBQT files')
parser.add_argument('--output_dir', required=True, help='Directory to save docking results')
parser.add_argument('--center', nargs=3, type=float, required=True, help='Grid center (x y z)')
parser.add_argument('--size', nargs=3, type=float, required=True, help='Grid size (x y z)')
parser.add_argument('--exhaustiveness', type=int, default=8, help='Docking exhaustiveness')
args = parser.parse_args()

RECEPTOR = args.receptor
DOCKING_DIR = args.ligand_dir
OUTPUT_DIR = args.output_dir
CENTER_X, CENTER_Y, CENTER_Z = args.center
SIZE_X, SIZE_Y, SIZE_Z = args.size
EXHAUSTIVENESS = str(args.exhaustiveness)

CONSOLE_OUTPUT = os.path.join(OUTPUT_DIR, "docking_console_output.txt")
OUTPUT_CSV = os.path.join(OUTPUT_DIR, "docking_scores.csv")
DEBUG_LOG = os.path.join(OUTPUT_DIR, "debug_log.txt")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Clear previous console output
with open(CONSOLE_OUTPUT, 'w') as console_file:
    pass

def write_debug(message):
    with open(DEBUG_LOG, 'a') as debug_file:
        debug_file.write(message + "\n")
    print(message)

write_debug(f"Grid box: center=({CENTER_X}, {CENTER_Y}, {CENTER_Z}), size=({SIZE_X}, {SIZE_Y}, {SIZE_Z})")

try:
    with open(OUTPUT_CSV, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Ligand", "Best Affinity"])
        csvfile.flush()
        os.fsync(csvfile.fileno())

        if not os.path.isdir(DOCKING_DIR):
            write_debug("Error: Docking directory does not exist: " + DOCKING_DIR)
            sys.exit(1)

        ligand_count = 0
        for ligand in os.listdir(DOCKING_DIR):
            if ligand.endswith(".pdbqt") and '_docked' not in ligand:
                ligand_count += 1
                ligand_path = os.path.join(DOCKING_DIR, ligand)
                ligand_name = os.path.splitext(ligand)[0]
                docked_file = os.path.join(OUTPUT_DIR, f"{ligand_name}_docked.pdbqt")

                write_debug("Processing ligand: " + ligand_name)

                command = [
                    "./vina_1.2.7_linux_x86_64",
                    "--receptor", RECEPTOR,
                    "--ligand", ligand_path,
                    "--out", docked_file,
                    "--center_x", str(CENTER_X),
                    "--center_y", str(CENTER_Y),
                    "--center_z", str(CENTER_Z),
                    "--size_x", str(SIZE_X),
                    "--size_y", str(SIZE_Y),
                    "--size_z", str(SIZE_Z),
                    "--exhaustiveness", EXHAUSTIVENESS
                ]

                try:
                    vina_output = subprocess.check_output(command, stderr=subprocess.STDOUT).decode('utf-8')
                    write_debug("Vina executed for " + ligand_name)

                    with open(CONSOLE_OUTPUT, 'a') as console_file:
                        console_file.write(f"Docking results for {ligand_name}:\n")
                        console_file.write(vina_output + "\n\n")

                    match = re.search(r'^\s*1\s+([-\d.]+)', vina_output, re.MULTILINE)
                    if match:
                        best_score = float(match.group(1))
                        csv_writer.writerow([ligand_name, best_score])
                        write_debug(f"Written: {ligand_name}, {best_score}")
                    else:
                        write_debug("No score for " + ligand_name)
                        csv_writer.writerow([ligand_name, "N/A"])

                    csvfile.flush()
                    os.fsync(csvfile.fileno())

                except subprocess.CalledProcessError as e:
                    write_debug(f"Error with {ligand_name}: {e}")
                    with open(CONSOLE_OUTPUT, 'a') as console_file:
                        console_file.write(f"Error docking {ligand_name}: {e}\n\n")
                    csv_writer.writerow([ligand_name, "Error"])
                    csvfile.flush()
                    os.fsync(csvfile.fileno())

        write_debug(f"Total ligands processed: {ligand_count}")

except Exception as e:
    write_debug("Unexpected error: " + str(e))

write_debug(f"Completed. Results: {OUTPUT_CSV}")
