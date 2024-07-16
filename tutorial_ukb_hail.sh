

###########################################################
# sample data from tutorial (https://hail.is/docs/0.2/tutorials/01-genome-wide-association-study.html)
###########################################################
# Date: Dec 2023, 2023

cd /UKBIOBANK/pangk/

# Run Hail on Python
cd /home/pangk/.local/lib/python3.11/ # this dir works 
/usr/bin/python-3.11.6/bin/python3.11 # run python3.11.6
import hail as hl
import hail.expr.aggregators as agg
hl.init()


hl.utils.get_1kg('data/') #Download subset of the 1000 Genomes dataset and sample annotations.
hl.import_vcf('data/1kg.vcf.bgz').write('data/1kg.mt', overwrite=True)
mt = hl.read_matrix_table('data/1kg.mt')
mt.rows().select().show(5)
mt.entry.take(5)
mt.count_cols()


