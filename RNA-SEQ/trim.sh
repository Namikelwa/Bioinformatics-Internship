#It's crusial to provide the path to your trimmomatic java jar to check if you have your java jar use this command:
 compgen -ac | grep trimmomatic
java -jar <path to trimmomatic.jar> PE <input 1> <input 2>] <paired output 1> <unpaired output 1> <paired output 2> <unpaired output 2> <options>

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
