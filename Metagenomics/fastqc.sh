#runfastqc
module load fastqc
fastqc -o ../fastqc-results ./*.fastq.gz
cd ../fastqc-results

#run multiqc
multiqc ./
