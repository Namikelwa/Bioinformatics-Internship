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


#Postfastq
fastqc -o postfastqc *.fastq.gz
multiqc ./
