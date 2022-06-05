
#Importing data

mkdir qiime
cd qiime

pip install send2trash
module load qiime2

# Running automated python script 
python3 ./manifest.py --input_dir ../trimmed-reads
sort Manifest.csv >> sorted.csv
sed -i '1s/^/sample-id,absolute-filepath,direction\n/' sorted.csv
sed '$d' sorted.csv >> ready.csv
rm Manifest.csv && rm sorted.csv
mv ready.csv Manifest.csv

#Import your data into Qiime2 Artifacts
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path Manifest.csv \
--output-path demux.qza \
--input-format PairedEndFastqManifestPhred33



##Denoising using dada2

#Using dada2
qiime dada2 denoise-paired \
--i-demultiplexed-seqs demux.qza \
--p-trim-left-f 17 \
--p-trim-left-r 17 \
--p-trunc-len-f 250 \
--p-trunc-len-r 200 \
--o-table dada2-table.qza \
--o-representative-sequences rep-seqs-dada2.qza \
--o-denoising-stats dada2-denoise-stats.qza

#Adding metadata and examining count tables
qiime feature-table summarize \
--i-table dada2-table.qza \
--o-visualization dada2-table.qzv \
--m-sample-metadata-file ./zanzibar1.tsv

qiime feature-table tabulate-seqs \
--i-data rep-seqs-dada2.qza \
--o-visualization rep-seqs-dada2.qzv

