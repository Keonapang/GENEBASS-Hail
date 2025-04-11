
#!/bin/bash
####################################################################
# Purpose: Filter Genebass beta results at ~7 million variants 
#       Description: Remove ChrX ChrY variants, and variants that are not in our UKB 200K WES RV dataset.
# Date: Feb 28, 2024 (UPDATED)
# by: Keona Pang 
####################################################################   
# Execute:
# nohup /UKBIOBANK/pangk/Genebass_results/SCRIPT_2_filter_align_pval.sh $trait >> /UKBIOBANK/pangk/Genebass_results/logfiles/SCRIPT_2_filter_align_pval.sh.log &

# Only argument to modify before running this script -- CANNOT RUN AUTOMATICALLY. CHECK PART 1's script: 
trait="T2D" 
bash /UKBIOBANK/pangk/Genebass_results/filter_align_pval.sh $trait
