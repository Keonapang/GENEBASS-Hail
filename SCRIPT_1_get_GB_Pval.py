
cd /your/file/path/Genebass

#################### 1. Run Hail on Python ##################################
PYSPARK_SUBMIT_ARGS="--driver-memory 12g --executor-memory 12g pyspark-shell" /usr/bin/python-3.11.6/bin/python3.11
import hail as hl
hl.init()
var_mt = hl.read_matrix_table('variant_results.mt') # (8074878, 4529)

################### option A: Export PVal for continous phenotype ######################
# currently using as of Feb 9, 2024

var_mt = var_mt.select_entries(var_mt.Pvalue)
var_mt = var_mt.key_cols_by()
var_mt = var_mt.key_rows_by()
fields_to_drop = ['n_cases', 'n_controls', 'heritability', 'saige_version', 'inv_normalized', 'n_cases_both_sexes', 'n_cases_defined', 'n_cases_males', 'description_more', 'coding_description', 'category', 'n_cases_females', 'annotation', 'call_stats', 'markerID', 'gene', 'coding', 'modifier', 'trait_type', 'pheno_sex']
var_mt = var_mt.drop(*fields_to_drop)

# trait='Glucose'
# phenocode = '30740'
trait='height'
phenocode = '50'
AA = var_mt.filter_cols(var_mt.phenocode == phenocode) #matrixtable
AA_entries = AA.entries() # converts to hail.table 
AA_entries.show(2)

# Drop 
AA_entries = AA_entries.drop(*['phenocode', 'description'])

# continous trait
continous_pval = AA_entries.head(8074878) 
continous_pval.export(f'/your/file/path/1_pval_raw/500k_{trait}_pval.txt.gz', delimiter='\t') 

############################################################
# Extract the UKB allele frequency 
############################################################
#   used on Feb 14, 2024
# PYSPARK_SUBMIT_ARGS="--driver-memory 12g --executor-memory 12g pyspark-shell" /usr/bin/python-3.11.6/bin/python3.11
# import hail as hl
# hl.init()
# var_mt = hl.read_matrix_table('variant_results.mt') # (8074878, 4529)

# var_mt = var_mt.select_entries(var_mt.AF)

# var_mt = var_mt.key_cols_by()
# var_mt = var_mt.key_rows_by()
# fields_to_drop = ['n_cases', 'n_controls', 'heritability', 'saige_version', 'inv_normalized', 'n_cases_both_sexes', 'n_cases_defined', 'n_cases_males', 'description', 'description_more', 'coding_description', 'category', 'n_cases_females', 'annotation', 'call_stats', 'markerID', 'gene', 'coding', 'modifier', 'trait_type', 'pheno_sex']
# var_mt = var_mt.drop(*fields_to_drop)

# trait="LDL_direct"
# phenocode="30780"

# AA = var_mt.filter_cols(var_mt.phenocode == phenocode) #matrixtable
# AA_entries = AA.entries() # converts to hail.table 
# AA_entries = AA_entries.drop(*['phenocode'])
# AA_entries.describe()
# AA_entries.show(2)

# continous_pval = AA_entries.head(8074878) 
# continous_pval.export(f'500k_AF_all.txt.gz', delimiter='\t') 

#################### Ooption B (old): Export BETA and SE of continous phenotype #######################
# old option 
#2. Filter/drop MatrixTable 
# var_mt = var_mt.select_entries(var_mt.BETA, var_mt.SE)

# var_mt = var_mt.key_cols_by()
# var_mt = var_mt.key_rows_by()
# fields_to_drop = ['n_cases', 'n_controls', 'heritability', 'saige_version', 'inv_normalized', 'n_cases_both_sexes', 'n_cases_defined', 'n_cases_males', 'description', 'description_more', 'coding_description', 'category', 'n_cases_females', 'annotation', 'call_stats', 'markerID', 'gene', 'coding', 'modifier', 'trait_type', 'pheno_sex']
# var_mt = var_mt.drop(*fields_to_drop)

# trait='LDL_direct'
# phenocode = '30780'

# AA = var_mt.filter_cols(var_mt.phenocode == phenocode) #matrixtable
# AA_entries = AA.entries() # converts to hail.table 
# AA_entries = AA_entries.drop(*['phenocode', 'locus', 'alleles'])

# continous_beta = AA_entries.head(8074878) # 4000000
# continous_beta.export(f'500k_{trait}_beta.txt.gz', delimiter='\t') 
