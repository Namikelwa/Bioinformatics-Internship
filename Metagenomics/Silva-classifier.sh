#download the silva database

wget https://data.qiime2.org/2021.4/common/silva-138-99-seqs.qza
wget https://data.qiime2.org/2021.4/common/silva-138-99-tax.qza

#train the feature classifier

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-99-seqs.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier classifier.qza
  
  
  #taxonomic classification
qiime feature-classifier classify-sklearn \
--i-classifier ./classifier.qza \
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
--o-visualization taxa-bar-plots.qzv

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

