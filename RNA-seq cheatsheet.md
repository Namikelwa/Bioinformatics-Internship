RNA-seq is used for differential gene expression data analysis.it offers transcriptional profiling of data
The first step to take in starting RNAseq analysisis is
obtaining data
raw data can be obtained from online repositories/databases i.e NCBI
the data can be accessesed by commandline arguments wget,fastq-dump...
If for instance you have many datasets you can access the accesion lists bu using esearch and your experiment id on the command line or you can go to the database i.e NCBI select the datasets and download the acc numbers into a txt file.

quality check
After having the raw reads with you quality analysis is done so as to know the quality of your data 
 fastqc tool
fastqc generates one output file, a html report which can be viewed on the browser and gives information about your sequences per sample
multiqc
this is a tool which compiles fastqc reports into one report
it makes fastqc outputs more manageable
it has 2 outputs a  single HTML report with plots to visualize and compare various QC metrics between the samples and a data file


filtering/trimming
this is a step done when the quality of your breads is not good 
tools utilized here include trimmomatic, a tool for trimming Illumina sequenced data reads,bowtie,Cutadapt  a tool for quality control of high-throughput sequencing reads. The functions include an adapter, primer, and poly-A tail removal
post fastqc
This is doneto visualize the trimming it is done  on the trimmed data followed by multiqc

mapping
requires a reference genome(obtained from database) and the trimmed reads
tools used include hisat2 which requires building an index for the reference genome 
abundance estimation

differential analysis
