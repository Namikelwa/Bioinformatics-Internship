# RNA-seq basic analysis steps

RNA-seq is used for differential gene expression data analysis.It offers transcriptional profiling of data  
Steps involved to do the analysis:  

## Obtaining data
Raw data can be obtained from online repositories/databases i.e NCBI        
If for instance you have many datasets you can access the accession lists using esearch together with your experiment id on the command line or you can go to the database i.e NCBI select the datasets and download the acc numbers into a txt file.           
Then the data can be downloaded with the help of commandline arguments *wget,fastq-dump,fasterq-dump*

Fetching raw data from NCBI       
The first step is to download the accession list into a txt file from the command line using **esearch**         
```
esearch -db sra -query <input exp id> | efetch --format runinfo | cut -d "," -f 1 > SraAccList.txt # provide the txt output file

# esearch command uses the ESearch utility to search on the NCBI database for query and finds the unique identifiers for all query that match the search query.
#-db database
#efetch downloads selected records in a style designated by -format
# -f flag is used to select the field you want i.e field with the accession numbers
```
Once you have the accession list,the next step is to download the data.          
**fastq-dump** is a tool in SRAtoolkit used for downloading sequenced reads from NCBI Sequence Read Archive(SRA).The data is dowloaded in fastq format.Here we are using the two options *--gzip* for compressing the sequence reads and *--split-files* to give both forward and reverse reads since the reads are from Illumina platform hence they are interleaved.

Getting the data one data-set at a time


```
fastq-dump --gzip --split-files <Accession number>
```
Compiling the data
```

for i in $(<input>);    # path to the SraAcclist.txt that contains a list of the accesion numbers
do
    echo $i #prints the acc numbers to the commandline
    fastq-dump --gzip --split-files $i  #fastq-dump gets data in fastq format 
done
```
## Pre-Quality analysis
quality check
After having the raw reads with you quality analysis is done so as to know the quality of your data for better downstream analysis.

**fastqc** checks the quality of the raw reads from high throughput sequencing platforms like Illumina.It provides one output file for each sample,a html report summary of the quality of reads in graphs and tables.The html report can be viewed on the browser and gives information about your sequences per sample.
**multiqc** makes fastqc output more manageable by compiling them and generating one report.
it has 2 outputs a  single HTML report with plots to visualize and compare various QC metrics between the samples and a data file
It can be done on the forward as  well as reverse reads   
After fastqc the outputs the 
```
fastqc -o <path to output dir>  <path to raw data with fastq.gz>
           or      
for i in `<path to the rawreads>`;
do
  	fastqc $i
done
multiqc ./
# forward and reverse multiqc
multiqc ./R1_fastqc.gz
multiqc ./R2_fastqc.gz
```

## Trimming
Filtering/trimming
This is a step done when the quality of your reads is not good     
Tools utilized here include trimmomatic-a tool for trimming Illumina sequenced data reads,Cutadapt-a tool for quality control of high-throughput sequencing reads..the functions include an adapter, primer, and poly-A tail removal.bowtie e.t.c       
Low quality reads together with adapters are removed after quality assessment.**Trimmomatic** is a  Java executable software used to trim and crop reads.     
```
#It's crusial to provide the path to your trimmomatic java jar to check if you have your java jar use this command:
 compgen -ac | grep trimmomatic
java -jar <path to trimmomatic.jar> PE <input 1> <input 2>] <paired output 1> <unpaired output 1> <paired output 2> <unpaired output 2> <options>
```
```
#basename is a command that strips trailing suffix in a file name
for i in *_1.fastq.gz;
do
    name=$(basename ${i} _1.fastq.gz) #flag -s for suffix can be used after basename followed by the suffix then ${i}
    trimmomatic PE -threads 20
                   ${i} ${name}_2.fastq.gz\#input files
                   ${name}_1.trim.fastq.gz ${name}_1un.trim.fastq.gz \ #output files for forward reads
                    ${name}_2.trim.fastq.gz ${name}_2un.trim.fastq.gz \ #output files for reverse reads
                    HEADCROP:11
done
PE - paired end
HEADCROP -removes the first 11 bases of the reads

```
## Post-Quality analysis  
The trimmed reads are then checked for quality   
fastqc  
multiqc   

## Mapping
Requires a reference genome in fa format (obtained from database) and the trimmed reads/raw reads if trimming wasn't done        
tools used include hisat2 which requires building an index for the reference genome         
**hisat2** is a fast and sensitive splice-aware aligner that compresses the genome using an indexing scheme to reduce the amount of space needed to store the genome. This also makes the genome quick to search, using a whole-genome index.          
We use samtools to convert the output file from sam to bam format(downsream anlysis requiresbam file) and to index the bam files.Indexing creates a searchable index of sorted bam files required in some programs.     

```
wget <input link to the reference fasta file>

#Building a reference genome index
hisat2-build -p25 <input path to fasta file>  <output path with an.idx extension>

#Run hisat2 using indexed reference

for i in $(<input file>);#provide path to the trimmed reads
do
   echo ${i}#just to find out from the cmd line whether the script is working
   hisat2 -p25 -x <path to the indexed ref genome>
               -1 ${i}_1.fastqc.gz -2 ${i}_2.fastqc.gz\ #forward and reverse reads and adding extensions
               -S <path to aligned reads>/${i}_hisat.sam # -S to tell that the file is in sam format aand adding the sam extension
   samtools view -Sb <path to aligned reads>/${i}_hisat.sam  | samtools sort  > <output path>/${i}_hisat_sorted.bam
   samtools index <output path>/${i}_hisat_sorted.bam
   rm <output path>/${i}_hisat.sam # we are removing sam files because they are huge hence take alot of space
done
```


**STAR Aligner**(Spliced Transcripts Alignment to a Reference)
STAR is a splice aware aligner designed to specifically address many of the challenges of RNA-seq data.It shows high accuracy and mapping speed.To generate an idx file it requires
a ref genome file in fasta format as well as aref gff/gtf file found together in the database.  
Alignemnt in STAR involves two steps;

1. Creating genome index
```
STAR --runThreadN 6  # number of threads\
    --runMode genomeGenerate \
    --genomeDir <input>r #path to store genome indices\
    --genomeFastaFiles <input ref genome> \
    --sjdbGTFfile <input corresponding gff file>\
    --sjdbOverhang 99 #readlength-1 --sjdbGTFtagExonParentTranscript gene
    
    # --sjdbGTFtagExonParentTranscript gene -allows star to use a gff file coz appparently it requires a gtf file
    
```

2. Mapping reads to the genome
```
for i in $(<input file/path to the file>);
do
    STAR --genomeDir  \ #path to indexed ref genome
    --readFilesIn  ${i}_1.fastq.gz ${i}_2.fastq.gz\#forward and reverse reads and provide the extensions
    --readFilesCommand zcat  \
    --outSAMtype BAM SortedByCoordinate \# bam file
    --quantMode GeneCounts \
    --outFileNamePrefix <output dir>/${i}
done
#zcat decompress the files

```

## Abundance estimation

Once you have your aligned reads,**htseq** is used to give counts of reads mapped to each feature.A feature is an interval on a chromosome.
```
htseq-count -t exon -i gene_id -f bam <input bam file/s>  <input gff file>  > <output file with a txt extension>
#-t featuretype
#-i attribute
#-f format
#a wild card can be used for multiple bam files
```

## Differential analysis

Differential analysis involves using read counts to perform statistical analysis to discover quantitative changes in gene expression levels between experimental groups;exposed and non-exposed.**DESeq2** is used for differential analysis. 
This is usually done in R

