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
