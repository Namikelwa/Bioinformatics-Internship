
# Fungi-identification project
## Data
The data consists of 26 samples from Zanzibar
## Dependencies

* fastqc/0.11.9
* multiqc/1.4
* trimmomatic/0.39
* qiime2/2020.6

## Workflow

### Quality check
```
#runfastqc
module load fastqc
fastqc -o ../fastqc-results ./*.fastq.gz
cd ../fastqc-results

#run multiqc
multiqc ./
```
From the multiqc [report](file:///home/icipe/Desktop/Fungi_Identification/results/multiqc_report.html).

 I found to have sequences of the length 301bp. With relatively good quality

![postmultiqc_qualityscores](https://user-images.githubusercontent.com/78470415/118993543-4475d680-b98e-11eb-9ed6-c12f343766a7.png)

 However, from the adapter content, we had three troublesome reads or samples.

They showed high adapter content, which are of the type Nextera PE-PE.fa

![fastqc_adapter_content_plot](https://user-images.githubusercontent.com/78470415/118993369-21e3bd80-b98e-11eb-998a-7f3bbefc2288.png)


## Quality control
Trimmed using trimmomatic
```
for i in *_R1_001.fastq.gz;
do
    name=$(basename -s _R1_001.fastq.gz ${i})
    trimmomatic PE -threads 20 ${i} ${name}_R2_001.fastq.gz \
                ${name}.R1.fastq.gz ${name}.R1.untrim.fastq.gz \
                ${name}.R2.fastq.gz ${name}.R2.untrim.fastq.gz \
                LEADING:6 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:100 ILLUMINACLIP:NexteraPE-PE.fa:2:30:10
done
rm -r *untrim*
mv *R*.fastq.gz ../trimmed-reads

```
### Post fastqc
```
fastqc -o postfastqc *.fastq.gz
multiqc ./
```
[post_trimmommatic_fastqc_report](file:///home/icipe/Desktop/Fungi_Identification/scratchpad/Zanzibar2_magroov_L001_R1_pair_001_fastqc.html#M4)   


From the results, the adapters were cut out, and overall quality improved.
Here is the adapter content results.

![Zanzibar 2 adapter content](assets/markdown-img-paste-20210520125748110.png)

However, the lengths are still not similar, ranging from 10 to 301bp. This will, however, be dealt with denoising steps.
## Importing and denoising

To import data into qiime, you will need to have a manifest file; this can be challenging, especially when you have many samples, so doing it manually is out of the question. To understand the use of a manifest file, here is a link to the  qiime2 documentation page on

[manifest files](https://docs.qiime2.org/2021.4/tutorials/importing/)   


## Importing data into qiime
Why qiime? The microbiome data contain sequences, and transitioning from raw sequence data to useful results often requires significant work, but computational tools such as QIIME2 have made it much easier.    
QIIME2 is a microbiome analysis package that could translate raw sequence data into useful statistical results

```
mkdir qiime
cd qiime

pip install send2trash
# Using send2trash, we can send files to the Trash or Recycle Bin instead of permanently deleting them
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

```
## Denoising using Dada2
17 bp are trimmed from each end to remove primers if any present and the reads are truncated at 250bp for forward reads and 200 reads as informed by the post fastqc report.
```
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
```

```
#Adding metadata and examining count tables
qiime feature-table summarize \
--i-table dada2-table.qza \
--o-visualization dada2-table.qzv \
--m-sample-metadata-file ./zanzibar1.tsv

qiime feature-table tabulate-seqs \
--i-data rep-seqs-dada2.qza \
--o-visualization rep-seqs-dada2.qzv
```
Sample ID |	Feature Count
--------- | ---------
Pemba | 158768
T | 147947
UZI02 | 130183
A | 122742
Q | 118121
S | 117305
O | 116026
N | 108091
Zanzibar1 | 106866
H | 104477
G | 103349
F | 102540
B | 101460
HO2 | 97772
V | 93520
M | 93107
R | 92397
P | 92339
TZC | 90426
HO1 | 88255
Zanzibar2 | 80679
C | 77444
D | 73628
U | 72359
HYPOSP1 | 37016
HYPOSP3 | 7381

![image](https://drive.google.com/uc?export=view&id=1D1_zo8LMqxqX4p__zpIjTqZCZ2qtoV4y)    

## Training a Classifier
### SILVA Database
```
#download the silva database

wget https://data.qiime2.org/2021.4/common/silva-138-99-seqs.qza
wget https://data.qiime2.org/2021.4/common/silva-138-99-tax.qza

#train the feature classifier

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-99-seqs.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier classifier.qza
  ```
## Taxonomic classification
```
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

```

### UNITE Database
#Using a trained classifier
```
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

```
#### Results
* Phyla present

**SILVA**                                                                                        
![Genus](https://user-images.githubusercontent.com/85280529/184815477-0a8ffe52-4617-43c9-b4bb-95359733b041.PNG)



**UNITE**

![silva](https://user-images.githubusercontent.com/85280529/184886886-083073d3-f469-401d-a557-c95e2e1858ba.PNG)

![silva species](https://user-images.githubusercontent.com/85280529/184886876-41f8200b-33a7-4ae7-8de8-06245d6faaa8.PNG)


## Discussion

The most abundant fungal genera present are ;
Silva | Unite
----- | -----
Yueomyces | Zygosaccharomyces
Talaromyces | Blumeria
Cladosporium | Issatchenkia
Phoma | Pyronema
Myrothecium | Paraophiobolus

The harmful species were not visible from classification using SILVA database but were present in UNITE classification.
The genera present are;
* Fusarium ;one specie(3.457% relative frequency)
* Penicillium ;one specie(2.432% relative frequency)
* Aspergillus ; one specie(5.263% relative frequency)

Zygosaccharomyces is seen to be abundant genera.Stingless bees rely on  Zygosaccharomyces for survival.The larvae feeds on the filaments in their initial stages of development.
Stingless bees live in symbiotic relationship with most fungi species which include; Zygosaccharomyces and Candida species found in this study.

The huge difference in the taxonomic representation between the two databases would be because;

* UNITE is a web-based database and sequence management environment for the molecular identification of fungi.It targets formal fungal barcode nuclear Internal Transcribed Spacer region.

* SILVA is a comprehensive web resource for quality controlled database of aligned rRNA gene sequences from three domains of life; bacteria,archaea and eukaryota


## Computing diversity


![shannon](https://user-images.githubusercontent.com/85280529/184886738-0db1fa12-2b49-4df0-8675-c3c5c708422e.PNG)

![diversity](https://user-images.githubusercontent.com/85280529/184886804-abc3e5c2-e7ae-40f0-9d2d-1d8a359bdd6a.PNG)

![PCA](https://user-images.githubusercontent.com/85280529/184886824-450f2294-5abd-4cbb-bd75-3b640c70e406.PNG)

 
 Challenge | Solution
 --------- | --------
 basename command | *basename -s _R1_001.fastq.gz file-name*
 Trimmomatic jar file | provided path  */opt/apps/trimmomatic/0.39/trimmomatic-0.39.jar*
 Adapter path | */opt/apps/trimmomatic/0.39/adapters/NexteraPE-PE.fa*
 Illumina clip | provide the argues NexteraPE-PE.fa:2:30:10
 Pip command | unload all modules ,then *pip install send2trash*
 Qiime feature-classifier| Running the command on a single line
 


