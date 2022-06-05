
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


