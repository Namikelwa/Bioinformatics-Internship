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
