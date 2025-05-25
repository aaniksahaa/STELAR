#!/usr/bin/env python3

import os
import subprocess
import time
import csv
from pathlib import Path
import re

# Define the Java command prefix (from Eclipse)
JAVA_CMD = [
    "/home/aaniksahaa/.p2/pool/plugins/org.eclipse.justj.openjdk.hotspot.jre.full.linux.x86_64_21.0.7.v20250502-0916/jre/bin/java",
    "-Dfile.encoding=UTF-8",
    "-Dstdout.encoding=UTF-8",
    "-Dstderr.encoding=UTF-8",
    "-classpath",
    "/mnt/H/Research/STELAR-extension/STELAR/main:/mnt/H/Research/STELAR-extension/STELAR/Data-set/11-taxon/astral.5.6.1.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/11-taxon:/mnt/H/Research/STELAR-extension/STELAR/Data-set/11-taxon/SuperTriplets_v1.1.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/15-taxon/lib/main.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/15-taxon/lib/JSAP-2.1.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/15-taxon/lib/colt.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/15-taxon/astral.5.6.1.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/15-taxon/STELAR.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/15-taxon:/mnt/H/Research/STELAR-extension/STELAR/Data-set/15-taxon/SuperTriplets_v0.31.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/15-taxon/SuperTriplets_v1.1.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/37-taxon/astral.5.6.1.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/37-taxon/STELAR.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/37-taxon:/mnt/H/Research/STELAR-extension/STELAR/Data-set/37-taxon/SuperTriplets_v1.1.jar:/mnt/H/Research/STELAR-extension/STELAR/DynaDup/main.jar:/mnt/H/Research/STELAR-extension/STELAR/main.jar:/mnt/H/Research/STELAR-extension/STELAR/mgd.jar:/mnt/H/Research/STELAR-extension/STELAR/miscellaneous-scripts/others/phylonet_v2_4.jar:/mnt/H/Research/STELAR-extension/STELAR/STELAR.jar",
    "-XX:+ShowCodeDetailsInExceptionMessages",
    "phylonet.coalescent.MGDInference_DP"
]

# Define outer folders (taxa numbers)
folders = ["37-taxon"]  # Adjust as needed: 11-taxon, 37-taxon, 48-taxon, etc.
fresh = 1  # Set to 0 to skip existing output files, 1 to overwrite

# Define methods to run
methods = ["base", "weighted_2_terminal", "weighted_3_terminal"]

# Define inner folder names for each taxa number
inner_folder_names = {
    "11-taxon": [
        "estimated_Xgenes_strongILS/estimated_5genes_strongILS",
        "estimated_Xgenes_strongILS/estimated_15genes_strongILS",
        "estimated_Xgenes_strongILS/estimated_25genes_strongILS",
        "estimated_Xgenes_strongILS/estimated_50genes_strongILS",
        "estimated_Xgenes_strongILS/estimated_100genes_strongILS"
    ],
    "15-taxon": [
        "100gene-100bp/estimated-genetrees",
        "100gene-1000bp/estimated-genetrees",
        "1000gene-100bp/estimated-genetrees",
        "1000gene-1000bp/estimated-genetrees"
    ],
    "37-taxon": [
        "estimated-genetrees/1X-200-500",
        "estimated-genetrees/1X-200-1000",
        "estimated-genetrees/1X-400-500",
        "estimated-genetrees/1X-400-1000",
        "estimated-genetrees/1X-800-500",
        "estimated-genetrees/1X-800-1000",
        "estimated-genetrees/2X-200-500"
    ],
    "48-taxon-latest": ["estimated_genetrees/1X-1000-500"],
    "100-taxon": ["inner100"],
    "200-taxon": ["inner200"],
    "500-taxon": ["model.500.2000000.0.000001"],
    "biological": ["nuclear"],
    "mammalian": ["424genes"],
    "amniota": ["aa", "nt"]
}

# Define number of replicates for each folder
replicates = {
    "11-taxon": 20,
    "15-taxon": 10,
    "37-taxon": 20,
    "48-taxon-latest": 20,
    "100-taxon": 8,
    "200-taxon": 10,
    "500-taxon": 1,
    "biological": 1,
    "mammalian": 1,
    "amniota": 1
}

# Create global CSV files for summary_new
summary_dir = Path("./datasets/summary_new")
summary_dir.mkdir(parents=True, exist_ok=True)
timediffcsv_main = summary_dir / "average_diffs.csv"
rf_csv_main = summary_dir / "average_rfs.csv"

with open(timediffcsv_main, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["summary_new_method", "folder", "inner_folder", "average_diff"])

with open(rf_csv_main, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["summary_new_method", "folder", "inner_folder", "average_rf"])

# Loop through each taxa folder
for folder in folders:
    print(f"Processing folder: {folder}")

    # Create local CSV files for each taxa number
    local_summary_dir = Path(f"./datasets/{folder}/summary_new")
    local_summary_dir.mkdir(parents=True, exist_ok=True)
    timediffcsv = local_summary_dir / "average_diffs.csv"
    rf_csv = local_summary_dir / "average_rfs.csv"

    with open(timediffcsv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["summary_new_method", "folder", "inner_folder", "average_diff"])

    with open(rf_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["summary_new_method", "folder", "inner_folder", "average_rf"])

    # Get inner folders and number of replicates
    inner_folders = inner_folder_names[folder]
    R = replicates[folder]

    # Loop through each inner folder
    for inner_folder in inner_folders:
        print(f"Processing inner folder: {inner_folder}")

        # Loop through each method
        for method in methods:
            print(f"Processing method: {method}")

            # Initialize lists to store diffs and RFs for averaging
            diffs = []
            rfs = []

            # Loop through replicates
            for j in range(1, R + 1):
                print(f"Processing replicate R{j} for method {method}")

                # Define paths
                gt_folder = f"{inner_folder}/R{j}"
                input_file = f"./datasets/{folder}/{gt_folder}/all_gt.tre.rooted"
                output_file = f"./datasets/{folder}/{gt_folder}/stelar_outputs/stelar_output_{method}.tre"
                true_tree = f"./datasets/{folder}/true_tree_trimmed"

                # Special handling for certain folders
                if folder in ["100-taxon", "200-taxon", "500-taxon"]:
                    true_tree = f"./datasets/{folder}/true-species-trees/R{j}/sp-cleaned"
                elif folder == "11-taxon-new":
                    true_tree = f"./datasets/{folder}/higher-ILS/true-speciestrees/R{j}.true.tre"
                elif folder == "1000-taxon":
                    gt_folder = f"{inner_folder}/{'0' + str(j) if j < 10 else j}"
                    true_tree = f"./datasets/{folder}/true-species-trees/{'0' + str(j) if j < 10 else j}/s_tree.trees"
                    input_file = f"./datasets/{folder}/{gt_folder}/stelar_inputs/stelar_input.tre"
                    output_file = f"./datasets/{folder}/{gt_folder}/stelar_outputs/stelar_output_{method}.tre"

                # Create output directory
                output_dir = Path(f"./datasets/{folder}/{gt_folder}/stelar_outputs")
                output_dir.mkdir(parents=True, exist_ok=True)

                # Skip if output exists and fresh is 0
                if os.path.exists(output_file) and fresh == 0:
                    print(f"Output exists, skipping: {output_file}")
                    continue

                # Ensure input file exists
                if not os.path.exists(input_file):
                    print(f"Error: Input file does not exist: {input_file}")
                    continue

                # Run Java codebase with method and measure time
                start_time = time.time()
                try:
                    subprocess.run(
                        JAVA_CMD + ["-m", method, "-i", input_file, "-o", output_file],
                        check=True,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True
                    )
                except subprocess.CalledProcessError as e:
                    print(f"Error running Java command: {e.stderr}")
                    continue
                end_time = time.time()
                diff = end_time - start_time
                print(f"Time taken for method {method}: {diff:.6f} seconds")

                # Calculate RF distance
                if os.path.exists(output_file) and os.path.exists(true_tree):
                    try:
                        result = subprocess.run(
                            ["python3", "./RF/getFpFn.py", "-e", output_file, "-t", true_tree],
                            check=True,
                            capture_output=True,
                            text=True
                        )
                        # Extract RF distance from tuple output (e.g., "(0,2)")
                        match = re.match(r"\((\d+),(\d+)\)", result.stdout.strip())
                        if match:
                            rf_distance = int(match.group(2))
                            print(f"RF Distance for method {method}: {rf_distance}")

                            # Store diff and RF
                            diffs.append(diff)
                            rfs.append(rf_distance)
                        else:
                            print(f"Error: Invalid RF output format: {result.stdout}")
                    except (subprocess.CalledProcessError, ValueError) as e:
                        print(f"Error calculating RF distance: {e}")
                else:
                    print("Error: Output or true tree file missing for RF calculation")

            # Calculate and store averages
            if diffs:
                average_diff = sum(diffs) / len(diffs)
                average_rf = sum(rfs) / len(rfs) if rfs else 0
                print(f"Average DIFF for {inner_folder} ({method}): {average_diff:.6f}")
                print(f"Average RF for {inner_folder} ({method}): {average_rf:.6f}")

                # Write to local CSV
                with open(timediffcsv, "a", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow([method, folder, inner_folder, f"{average_diff:.6f}"])

                with open(rf_csv, "a", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow([method, folder, inner_folder, f"{average_rf:.6f}"])

                # Write to global CSV
                with open(timediffcsv_main, "a", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow([method, folder, inner_folder, f"{average_diff:.6f}"])

                with open(rf_csv_main, "a", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow([method, folder, inner_folder, f"{average_rf:.6f}"])

            # Add newline to CSVs
            for csv_file in [timediffcsv, rf_csv, timediffcsv_main, rf_csv_main]:
                with open(csv_file, "a", newline="") as f:
                    f.write("\n")

print("Processing complete. Results saved in CSV files.")