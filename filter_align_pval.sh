
#!/bin/bash
####################################################################
# Purpose: Filter Genebass beta results at ~7 million variants 
#       Description: Remove ChrX ChrY variants, and variants that are not in our UKB 200K WES RV dataset.
# Date: Feb 28, 2024 (UPDATED)
# by: Keona Pang 
####################################################################   
# Execute:
# trait="Lipoprotein_A" 
# bash /UKBIOBANK/pangk/Genebass_results/SCRIPT_2_filter_align_pval.sh $trait
####################################################################

trait=$1

# trait="Glycated_haemoglobin" #Calcium, Creatinine, Direct_bilirubin, Gamma_glutamyltransferase, Lipoprotein_A, albumin, Aspartate_aminotransferase, ApoB, ApoA
echo "================================================================================"
start_date=$(date)
echo "Script: Filter raw $trait pval files from GENEBASS."
echo "Script started at: $start_date"
echo ""

# zcat $file | cut -f3,4 --complement | gzip > ${trait}_temp.txt.gz && mv ${trait}_temp.txt.gz $file 
###########################################################################################
# PART 1 - Processing raw p-val file 
# Objective: convert to space-delimited, remove brackets/commas, creates PLINK_SNP_NAME,  RENAMES headers
# Output: space delimited, edited "500k_${trait}_pval.txt.gz"  
###########################################################################################

dir="/UKBIOBANK/pangk/Genebass_results"
cd ${dir}/1_pval_raw
echo "[PART 1] Directory: ${dir}/1_pval_raw"

# Input pval.txt.gz file 
file="500k_${trait}_pval.txt.gz"
num_lines=$(zcat "$file" | wc -l)
num_columns=$(zcat "$file" | awk -F'\t' '{print NF; exit}')
echo "Dimension of raw $file is $num_lines x $num_columns"
echo ""
zcat $file | head -n 5
echo ""

# Process: convert to space-delimited, remove brackets, replace commas with colons, remove 2nd column, add header
echo "Removing brackets"
zcat $file | tr -d '[]"' | gzip > ${trait}_temp.txt.gz && mv ${trait}_temp.txt.gz $file 
zcat $file | head -n 2
echo ""
num_lines=$(zcat "$file" | wc -l)
num_columns=$(zcat "$file" | awk -F'\t' '{print NF; exit}')
echo "Dimension of raw $file is $num_lines x $num_columns"
echo ""

echo "Replaced colon with commas"
zcat $file | awk -F' ' 'BEGIN{OFS=FS} {gsub(",", ":", $2); print}' | gzip > ${trait}_temp.txt.gz && mv ${trait}_temp.txt.gz $file # substitute comma with ":""
zcat $file | head -n 2
echo ""
num_lines=$(zcat "$file" | wc -l)
num_columns=$(zcat "$file" | awk -F' ' '{print NF; exit}')
echo "Dimension of raw $file is $num_lines x $num_columns"
echo ""

echo "Joined locus and alleles to form PLINK SNP NAME in col 1"
zcat $file | awk -F' ' 'BEGIN{OFS=FS} {$1=$1":"$2; $2=""; print}' | gzip > ${trait}_temp.txt.gz && mv ${trait}_temp.txt.gz $file # joins col 1 locus and 2 allele
zcat $file | head -n 2
zcat $file | tail -n 2
echo ""
num_lines=$(zcat "$file" | wc -l)
num_columns=$(zcat "$file" | awk -F' ' '{print NF; exit}')
echo "Dimension of raw $file is $num_lines x $num_columns"
echo ""

echo "Removed prefix "chr" from column 1 PLINK IDs"
zcat $file | awk -F' ' 'BEGIN{OFS=FS} NR>1{gsub(/^chr/, "", $1)}1' | gzip > ${trait}_temp.txt.gz && mv ${trait}_temp.txt.gz $file # remove prefix "chr"
zcat $file | head -n 2
zcat $file | tail -n 2
num_lines=$(zcat "$file" | wc -l)
num_columns=$(zcat "$file" | awk -F' ' '{print NF; exit}')
echo "Dimension of raw $file is $num_lines x $num_columns"
echo ""
echo "------------- done reformatting ---------------------"
echo ""
# renames column header to PLINK_SNP_NAME {trait}_pval
echo "Renaming headers to: PLINK_SNP_NAME and trait_pval"
echo ""
trait_name=${file#500k_} # Remove prefix
trait_name=${trait%.txt.gz} # Remove suffix
zcat $file | awk -F' ' -v trait_header=$trait_name 'BEGIN{OFS=FS} NR==1{$0="PLINK_SNP_NAME " trait_header}1' | gzip > ${trait}_temp.txt.gz && mv ${trait}_temp.txt.gz $file 
zcat $file | head -n 2
echo ""

# Count 'NA' in p-value column of the file 
num_NA=$(zcat "$file" | awk -F' ' 'NR>1 && $2=="NA"{count++} END{print count+0}')
echo "Number of NA in pvalues of $file before removal: $num_NA"

# Check if there is NA
if [ "$num_NA" -gt 0 ]; then
    echo "Number of NA among pvalues of $file before removal: $num_NA"

    # Remove all lines with NA in the second column
    output_file="${filter_dir}/$(basename $file)" # i.e. 500k_Alkaline_phosphatase_pval.txt.gz
    zcat $file | awk -F' ' 'BEGIN{OFS=FS} NR>1 && !($1 ~ /^Y:/ || $1 ~ /^X:/) && $2!="NA"' | gzip > $output_file

    # Count NA after removal
    num_NA_after=$(zcat "$output_file" | awk -F' ' 'NR>1 && $2=="NA"{count++} END{print count+0}')
    echo "Number of NA in pvalues of $file after removal: $num_NA_after"
else
    echo ""
    echo "No 'NA' (missing values) among the pvalues of $file"
fi
echo ""

num_lines=$(zcat "$file" | wc -l)
num_columns=$(zcat "$file" | awk -F' ' '{print NF; exit}')
echo "$(basename $file) dimension after re-formatting: $num_lines x $num_columns" 
echo ""
echo "END of Part 1: process and overwrite raw ${trait} pval files in /1_pval_raw."
echo "========================================================================"
echo ""

###########################################################################################
# PART 2 FILTERING: Remove Chr X/Y from p-value files and count number of NA (missing) values
#     - Output: Column vector of p-vals, but number of lines are reduced 
###########################################################################################

filter_dir="${dir}/2_pval_filtered"
cd ${filter_dir}
echo "[PART 2] Directory: ${filter_dir}"
echo ""

# Remove Chr X and Y 
output_file="${filter_dir}/$(basename $file)" # i.e. 500k_Alkaline_phosphatase_pval.txt.gz
zcat ${dir}/1_pval_raw/$file | awk -F' ' 'BEGIN{OFS=FS} !($1 ~ /^Y:/ || $1 ~ /^X:/)' | gzip > $output_file
echo "Removed ChrX and ChrY" 
echo ""
zcat $output_file | head -n 2
echo ""

# post-filtering dimensions
num_lines=$(zcat "$output_file" | wc -l)
num_columns=$(zcat "$output_file" | awk -F' ' '{print NF; exit}')
echo "$(basename $file) dimension after ChrY/X removed: $num_lines x $num_columns" 
echo ""
echo "END of Part 2: Removed chrX/Y from raw ${trait} pval files. Results in /2_pval_filtered"
echo "===========================[ Part 2: Filtering done ]====================================="
echo ""

###########################################################################################
# PART 3: [Method B] p-val to FDR conversion using the ENTIRE p-val file of ~7 million 
# Script: /UKBIOBANK/pangk/Genebass_results/2_pval_to_FDR.r
# Input (argument): one 500k_trait_pval.txt.gz file containing PLINK_SNP_NAME and trait_pval
# Output: space-delimited .txt to /UKBIOBANK/pangk/Genebass_results/3_FDR; renames second column to trait_FDR
        # PLINK_SNP_NAME height_FDR
        # 1:925920:A:C 0.771294151937272
# tool: qvalue package in R
###########################################################################################

filter_dir="${dir}/2_pval_filtered"
cd ${filter_dir}
echo "[PART 3: FDR Conversion]: ${filter_dir}"
echo ""
echo "Performing FDR conversion on ${trait} pval file"

# Check the number of columns in the file
num_columns=$(zcat ${filter_dir}/500k_${trait}_pval.txt.gz | awk -F' ' '{print NF; exit}')

# If the number of columns is 2, run run script:
if [ "$num_columns" -eq 2 ]; then
  Rscript ${dir}/2_pval_to_FDR.r ${filter_dir}/500k_${trait}_pval.txt.gz
else
  echo "The file does not contain two columns."
fi
echo ""
echo "Output:"
head -n 3 /UKBIOBANK/pangk/Genebass_results/3_FDR/${trait}_FDR_7m.txt
echo ""
echo "===========================[ Part 3: P-val to FDR conversion DONE ]====================================="

###########################################################################################
#   PART 4 - FILTER to align our PLINK FDR file to the exact order of the hg38 PLINK (my list) 
###########################################################################################
FDR_DIR="${dir}/3_FDR"
cd $FDR_DIR
echo "[PART 4] Directory: ${FDR_DIR}"
echo ""
echo "Execute function align_2col_by_SNP() to align FDR file to UKB 200k variants"
# sed -i 's/\r$//' /UKBIOBANK/pangk/Other_scripts/Functions/align_2col_by_SNP.sh

align_2col_by_SNP() {
    file1=$1 # vector of PLINK_SNP_NAME 
    file2=$2 # larger 2-col annotation table to filter 
    file3=$3 # output a space-delimited annotation table
    # Process the files and write the output to file3
    awk -F' ' 'NR==FNR{a[$1]=$2;next} {print $1, ($1 in a ? a[$1] : "NA")}' "$file2" "$file1" > "$file3"

    # Check if each line in "$file3" has the same number of columns
    min_columns=$(awk '{print NF}' "$file3" | sort -nu | head -n 1)
    max_columns=$(awk '{print NF}' "$file3" | sort -nu | tail -n 1)
    num_rows=$(wc -l < "$file3")
    if [ "$min_columns" -eq "$max_columns" ]; then
        echo "All lines have the same number of columns :))"
        echo "Dimension: $num_rows x $min_columns"
    else
        echo "WARNING: Not all lines have the same number of columns. Number of rows: $num_rows. Min cols: "$min_columns", max cols: "$max_columns""
    fi
}

hg38_genebass="/UKBIOBANK/pangk/20231101_annotation/hg38_noChr23_genebass.txt" # 3,301,826 x 1
trait_FDR="${trait}_FDR_7m.txt"
anno_aligned_PLINK="${FDR_DIR}/${trait}_FDR_noChr23_genebass.txt" #  output file 

align_2col_by_SNP "$hg38_genebass" "$trait_FDR" "$anno_aligned_PLINK"
echo "output: $anno_aligned_PLINK"
echo ""
head -n 3 $anno_aligned_PLINK
echo ""
echo "===========================[ Part 4: Align GENEBASS FDR to UKB dataset ]====================================="
echo ""

###########################################################################################
#   PART 5 - FILTER to align FDR to MAF threshold of simulation dataset 
###########################################################################################

cd $FDR_DIR
trait_FDR="${trait}_FDR_7m.txt"

echo "Filtering FDR values to keep only variants with 0.0001 < MAF < 0.01"
hg38_AF0001="/UKBIOBANK/pangk/20231101_annotation/hg38_noChr23_genebass_AF0.0001.txt" # 0.0001 < MAF < 0.01
FDR_0001="${FDR_DIR}/${trait}_FDR_noChr23_genebass_AF0.0001.txt" #  output file, 157,974 x 2
align_2col_by_SNP $hg38_AF0001 $trait_FDR $FDR_0001
echo ""
head -n 3 $FDR_0001
echo ""
echo "Filtering FDR values to keep only variants with 0.00001 < MAF < 0.01"
hg38_AF00001="/UKBIOBANK/pangk/20231101_annotation/hg38_noChr23_genebass_AF0.00001.txt" # 0.00001 < MAF < 0.01
FDR_00001="${FDR_DIR}/${trait}_FDR_noChr23_genebass_AF0.00001.txt" #  output file 
align_2col_by_SNP $hg38_AF00001 $trait_FDR $FDR_00001
echo ""
head -n 3 $FDR_00001
echo ""
end_date=$(date)
runtime_seconds=$(( $(date -d"$end_date" +%s) - $(date -d"$start_date" +%s) ))
runtime_minutes=$(printf "%.0f\n" $(echo "$runtime_seconds / 60" | bc -l))

echo ""
echo "Start date: $start_date"
echo "End date: $end_date"
echo "Runtime duration: $runtime_minutes minutes"
echo ""

echo "END of script: Three filtered ${trait} FDR files. Results in /3_FDR"
echo "They are ready to be merged to form $trait training dataset."
echo "======== [$trait]: Done SCRIPT_2_filter_align_pval.sh ========="
echo ""
echo "Next, run: /20231101_annotation/SCRIPT_9_add_FDR_tot_3m.sh"


#============================== End of script to process one trait from Genebass ==========================================


###################### to process multiple traits from Genebass at a time ############################
# echo "================================================================================"
# start_date=$(date)
# echo "Purpose: Filter raw GENEBASS beta files from "beta_output_20240115" directory, to remove ChrX/Y, and variants that are not in our UKB 200K WES rare variant dataset."
# echo "Output: Filtered/aligned beta PLINK files are found in beta_filtered_20240115 directory, all dimensions are 3315361 x 2."
# echo "Script started at: $start_date. The script for this code: ${0}."
# echo ""

# ###########################################################################################
# # PART 1 - Processing and overwrites raw p-val file
# # Detail: convert to space-delimited, remove brackets/commas, creates PLINK_SNP_NAME column and RENAMES header
# ###########################################################################################
# cd ${dir}/1_pval_raw

# for file in *_pval.txt.gz
# do
#     num_lines=$(zcat "$file" | wc -l)
#     num_columns=$(zcat "$file" | awk -F'\t' '{print NF; exit}')
#     echo "Dimensions of raw $file is $num_lines lines and $num_columns columns"

#     # Process: convert to space-delimited, remove brackets, replace commas with colons, remove 2nd column, add header
#     zcat $file | tr -d '[]"' | gzip > temp.txt.gz && mv temp.txt.gz $file 
#     zcat $file | awk -F' ' 'BEGIN{OFS=FS} {gsub(",", ":", $2); print}' | gzip > temp.txt.gz && mv temp.txt.gz $file # substitute comma with ":""
#     zcat $file | awk -F' ' 'BEGIN{OFS=FS} {$1=$1":"$2; $2=""; print}' | gzip > temp.txt.gz && mv temp.txt.gz $file # joins col 1 locus and 2 allele
#     zcat $file | awk -F' ' 'BEGIN{OFS=FS} NR>1{gsub(/^chr/, "", $1)}1' | gzip > temp.txt.gz && mv temp.txt.gz $file # remove prefix "chr"

#     trait=${file#500k_} # Remove the prefix
#     trait=${trait%.txt.gz} # Remove the suffix
#     zcat $file | awk -F' ' -v trait=$trait 'BEGIN{OFS=FS} NR==1{$0="PLINK_SNP_NAME " trait}1' | gzip > temp.txt.gz && mv temp.txt.gz $file # renames column headers

# done
# echo "END of Part 1: processing raw beta files. Beta files were overwritten."
# echo "========================================================================"

# ###########################################################################################
# # PART 2 FILTERING: Remove Chr X/Y from p-value files 
# #     - Output: Column vector of p-vals, but number of lines are reduced 
# ###########################################################################################

# for file in *_pval.txt.gz
# do
#     num_lines=$(zcat "$file" | wc -l)
#     num_columns=$(zcat "$file" | awk -F' ' '{print NF; exit}')
#     echo "Dimensions of $file pre-filtering $num_lines lines and $num_columns columns"

#     output_file="${filter_dir}/$(basename $file)"
#     zcat $file | awk -F' ' 'BEGIN{OFS=FS} !($1 ~ /^Y:/ || $1 ~ /^X:/)' | gzip > $output_file

#     num_lines=$(zcat "$output_file" | wc -l)
#     num_columns=$(zcat "$output_file" | awk -F' ' '{print NF; exit}')
#     echo "$output_file dimension ChrY/X removed: $num_lines x $num_columns" 

# done
# echo "END of Part 2: Filtering beta files to remove ChrY and ChrX. All outputs to "pval_filtered_20240212" directory."
# echo "=================================================="
