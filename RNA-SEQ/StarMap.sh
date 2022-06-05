##Creating genome indexes

STAR --runThreadN 6  # number of threads\
    --runMode genomeGenerate \
    --genomeDir <input>r #path to store genome indices\
    --genomeFastaFiles <input ref genome> \
    --sjdbGTFfile <input corresponding gff file>\
    --sjdbOverhang 99 #readlength-1 --sjdbGTFtagExonParentTranscript gene
    
    # --sjdbGTFtagExonParentTranscript gene -allows star to use a gff file coz appparently it requires a gtf file
    
   
   
   #Mapping
   
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
