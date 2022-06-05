esearch -db sra -query <input exp id> | efetch --format runinfo | cut -d "," -f 1 > SraAccList.txt # provide the txt output file

# esearch command uses the ESearch utility to search on the NCBI database for query and finds the unique identifiers for all query that match the search query.
#-db database
#efetch downloads selected records in a style designated by -format
# -f flag is used to select the field you want i.e field with the accession numbers

#downloading data one at a time
fastq-dump --gzip --split-files <Accession number>

#Getting the data downloaded at once
for i in $(<input>);    # path to the SraAcclist.txt that contains a list of the accesion numbers
do
    echo $i #prints the acc numbers to the commandline
    fastq-dump --gzip --split-files $i  #fastq-dump gets data in fastq format 
done
