#!/bin/bash
# Script to run Java codebase for inferring species trees from gene trees and compute RF scores.
# Supports multiple methods (base, weighted_2_terminal, weighted_3_terminal) and handles complex Java command.

# Define the Java command prefix (from Eclipse)
JAVA_CMD="/home/aaniksahaa/.p2/pool/plugins/org.eclipse.justj.openjdk.hotspot.jre.full.linux.x86_64_21.0.7.v20250502-0916/jre/bin/java \
-Dfile.encoding=UTF-8 \
-Dstdout.encoding=UTF-8 \
-Dstderr.encoding=UTF-8 \
-classpath /mnt/H/Research/STELAR-extension/STELAR/main:/mnt/H/Research/STELAR-extension/STELAR/Data-set/11-taxon/astral.5.6.1.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/11-taxon:/mnt/H/Research/STELAR-extension/STELAR/Data-set/11-taxon/SuperTriplets_v1.1.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/15-taxon/lib/main.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/15-taxon/lib/JSAP-2.1.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/15-taxon/lib/colt.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/15-taxon/astral.5.6.1.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/15-taxon/STELAR.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/15-taxon:/mnt/H/Research/STELAR-extension/STELAR/Data-set/15-taxon/SuperTriplets_v0.31.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/15-taxon/SuperTriplets_v1.1.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/37-taxon/astral.5.6.1.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/37-taxon/STELAR.jar:/mnt/H/Research/STELAR-extension/STELAR/Data-set/37-taxon:/mnt/H/Research/STELAR-extension/STELAR/Data-set/37-taxon/SuperTriplets_v1.1.jar:/mnt/H/Research/STELAR-extension/STELAR/DynaDup/main.jar:/mnt/H/Research/STELAR-extension/STELAR/main.jar:/mnt/H/Research/STELAR-extension/STELAR/mgd.jar:/mnt/H/Research/STELAR-extension/STELAR/miscellaneous-scripts/others/phylonet_v2_4.jar:/mnt/H/Research/STELAR-extension/STELAR/STELAR.jar \
-XX:+ShowCodeDetailsInExceptionMessages phylonet.coalescent.MGDInference_DP"

# Define outer folders (taxa numbers)
folders=("37-taxon")  # Adjust as needed: 11-taxon, 37-taxon, 48-taxon, etc.
fresh=1  # Set to 0 to skip existing output files, 1 to overwrite

# Define methods to run
methods=("base" "weighted_2_terminal" "weighted_3_terminal")

# Define inner folder names for each taxa number
# Define inner folder names for each taxa number as space-separated strings
declare -A innerFolderNames
innerFolderNames["11-taxon"]="estimated_Xgenes_strongILS/estimated_5genes_strongILS estimated_Xgenes_strongILS/estimated_15genes_strongILS estimated_Xgenes_strongILS/estimated_25genes_strongILS estimated_Xgenes_strongILS/estimated_50genes_strongILS estimated_Xgenes_strongILS/estimated_100genes_strongILS"
innerFolderNames["15-taxon"]="100gene-100bp/estimated-genetrees 100gene-1000bp/estimated-genetrees 1000gene-100bp/estimated-genetrees 1000gene-1000bp/estimated-genetrees"
innerFolderNames["37-taxon"]="estimated-genetrees/1X-200-500 estimated-genetrees/1X-200-1000 estimated-genetrees/1X-400-500 estimated-genetrees/1X-400-1000 estimated-genetrees/1X-800-500 estimated-genetrees/1X-800-1000 estimated-genetrees/2X-200-500"
innerFolderNames["48-taxon-latest"]="estimated_genetrees/1X-1000-500"
innerFolderNames["100-taxon"]="inner100"
innerFolderNames["200-taxon"]="inner200"
innerFolderNames["500-taxon"]="model.500.2000000.0.000001"
innerFolderNames["biological"]="nuclear"
innerFolderNames["mammalian"]="424genes"
innerFolderNames["amniota"]="aa nt"





# Define number of replicates for each folder
declare -A replicates
replicates["11-taxon"]=20
replicates["15-taxon"]=10
replicates["37-taxon"]=20
replicates["48-taxon-latest"]=20
replicates["100-taxon"]=8
replicates["200-taxon"]=10
replicates["500-taxon"]=1
replicates["biological"]=1
replicates["mammalian"]=1
replicates["amniota"]=1

# Create global CSV files for summary_new
mkdir -p ./datasets/summary_new
timediffcsv_main=./datasets/summary_new/average_diffs.csv
rf_csv_main=./datasets/summary_new/average_rfs.csv
printf "%s,%s,%s,%s\n" "summary_new_method" "folder" "inner_folder" "average_diff" > "$timediffcsv_main"
printf "%s,%s,%s,%s\n" "summary_new_method" "folder" "inner_folder" "average_rf" > "$rf_csv_main"

# Loop through each taxa folder
for folder in "${folders[@]}"; do
    echo "Processing folder: $folder"

    # Create local CSV files for each taxa number
    mkdir -p "./datasets/$folder/summary_new"
    timediffcsv="./datasets/$folder/summary_new/average_diffs.csv"
    rf_csv="./datasets/$folder/summary_new/average_rfs.csv"
    printf "%s,%s,%s,%s\n" "summary_new_method" "folder" "inner_folder" "average_diff" > "$timediffcsv"
    printf "%s,%s,%s,%s\n" "summary_new_method" "folder" "inner_folder" "average_rf" > "$rf_csv"

    # Get inner folders and number of replicates
    # inner_folders=("${innerFolderNames[$folder]}")
    IFS=' ' read -r -a inner_folders <<< "${innerFolderNames[$folder]}"

    R=${replicates[$folder]}

    # Loop through each inner folder
    for inner_folder in "${inner_folders[@]}"; do
        echo "Processing inner folder: $inner_folder"

        # Loop through each method
        for method in "${methods[@]}"; do
            echo "Processing method: $method"

            # Initialize arrays to store sums and counts for averaging
            sum_diffs=0
            sum_rfs=0
            count_diffs=0

            # Loop through replicates
            for ((j=1; j<=R; j++)); do
                echo "Processing replicate R$j for method $method"

                # Define paths
                gt_folder="$inner_folder/R$j"
                input="./datasets/$folder/$gt_folder/all_gt.tre.rooted"
                output="./datasets/$folder/$gt_folder/stelar_outputs/stelar_output_${method}.tre"
                true_tree="./datasets/$folder/true_tree_trimmed"

                # Special handling for certain folders
                if [[ "$folder" == "100-taxon" || "$folder" == "200-taxon" || "$folder" == "500-taxon" ]]; then
                    true_tree="./datasets/$folder/true-species-trees/R$j/sp-cleaned"
                elif [[ "$folder" == "11-taxon-new" ]]; then
                    true_tree="./datasets/$folder/higher-ILS/true-speciestrees/R$j.true.tre"
                elif [[ "$folder" == "1000-taxon" ]]; then
                    if [[ $j -lt 10 ]]; then
                        gt_folder="$inner_folder/0$j"
                        true_tree="./datasets/$folder/true-species-trees/0$j/s_tree.trees"
                    else
                        gt_folder="$inner_folder/$j"
                        true_tree="./datasets/$folder/true-species-trees/$j/s_tree.trees"
                    fi
                    input="./datasets/$folder/$gt_folder/stelar_inputs/stelar_input.tre"
                    output="./datasets/$folder/$gt_folder/stelar_outputs/stelar_output_${method}.tre"
                fi

                # Create output directory
                mkdir -p "./datasets/$folder/$gt_folder/stelar_outputs"

                # Skip if output exists and fresh is 0
                if [[ -f "$output" && $fresh -eq 0 ]]; then
                    echo "Output exists, skipping: $output"
                    continue
                fi

                # Ensure input file exists
                if [[ ! -f "$input" ]]; then
                    echo "Error: Input file does not exist: $input"
                    continue
                fi

                # Run Java codebase with method and measure time
                START=$(date +%s.%N)
                $JAVA_CMD -m "$method" -i "$input" -o "$output"
                END=$(date +%s.%N)
                DIFF=$(echo "$END - $START" | bc)
                echo "Time taken for method $method: $DIFF seconds"

                # Calculate RF distance
                if [[ -f "$output" && -f "$true_tree" ]]; then
                    index_to_extract=2
                    tuple=$(python3 ./RF/getFpFn.py -e "$output" -t "$true_tree")
                    RFdistance=$(echo "$tuple" | cut -d',' -f"$index_to_extract" | tr -d '()')
                    echo "RF Distance for method $method: $RFdistance"

                    # Update sums and count
                    sum_diffs=$(echo "$sum_diffs + $DIFF" | bc)
                    sum_rfs=$(echo "$sum_rfs + $RFdistance" | bc)
                    ((count_diffs++))
                else
                    echo "Error: Output or true tree file missing for RF calculation"
                fi
            done

            # Calculate and store averages
            if [[ $count_diffs -gt 0 ]]; then
                average_diff=$(echo "$sum_diffs / $count_diffs" | bc -l)
                average_rf=$(echo "$sum_rfs / $count_diffs" | bc -l)
                printf "Average DIFF for %s (%s): %.6f\n" "$inner_folder" "$method" "$average_diff"
                printf "Average RF for %s (%s): %.6f\n" "$inner_folder" "$method" "$average_rf"

                # Write to local CSV
                printf "%s,%s,%s,%.6f\n" "$method" "$folder" "$inner_folder" "$average_diff" >> "$timediffcsv"
                printf "%s,%s,%s,%.6f\n" "$method" "$folder" "$inner_folder" "$average_rf" >> "$rf_csv"

                # Write to global CSV
                printf "%s,%s,%s,%.6f\n" "$method" "$folder" "$inner_folder" "$average_diff" >> "$timediffcsv_main"
                printf "%s,%s,%s,%.6f\n" "$method" "$folder" "$inner_folder" "$average_rf" >> "$rf_csv_main"
            fi
        done

        # Add newline to CSVs
        printf "\n" >> "$timediffcsv"
        printf "\n" >> "$rf_csv"
        printf "\n" >> "$timediffcsv_main"
        printf "\n" >> "$rf_csv_main"
    done
done

echo "Processing complete. Results saved in CSV files."
