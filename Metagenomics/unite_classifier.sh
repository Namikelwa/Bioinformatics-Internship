
# taxonomic classification
qiime feature-classifier classify-sklearn \
--i-classifier ../feature-classifier/qime2/unite-ver8-99-classifier-04.02.2020.qza\
--i-reads ./rep-seqs-dada2.qza \
--o-classification ./taxonomy.qza

#Generating taxonomy file
qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

#creating taxonomy barplots
qiime taxa barplot \
--i-table dada2-table.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file zanzibar1.tsv \
--o-visualization taxa-barplot.qzv

# Filtering the unassigned to species level
qiime taxa filter-table \
  --i-table dada2-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-include s_ \
  --o-filtered-table table-species-level.qza

qiime metadata tabulate \
--m-input-file table-species-level.qza \
--o-visualization taxonomy_species.qzv

qiime taxa barplot \
--i-table table-species-level.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file zanzibar1.tsv \
--o-visualization taxa-barplot-species.qzv
