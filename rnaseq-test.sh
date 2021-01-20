###############
# DOWNLOAD DATA
###############
# download homd genomes for reference
wget http://www.homd.org/ftp/HOMD_prokka_genomes/gff/ALL_genomes.gff
wget http://www.homd.org/ftp/HOMD_prokka_genomes/fna/ALL_genomes.fna
# test dataset
# CHANGE THIS TO YOUR DATASET
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR340/004/SRR3404944/SRR3404944_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR340/004/SRR3404944/SRR3404944_2.fastq.gz
# get 16S rRNA database
wget https://bioinfo.lifl.fr/sortmerna/material/set1-database.fasta.zip
unzip set1-database.fasta.zip

############################
# ACTIVATE CONDA ENVIRONMENT
############################
conda env create -f environment.yml
conda activate rnaseq-test
# to deactivate:
# conda deactivate 

#########################
# QUALTIY FILTER AND TRIM
#########################
# first generate quality score metrics
ls *fastq.gz | parallel 'fastqc {}'
# download html files to local machine, open in web browser: 
# scp allie@stella.clemson.edu:/home/allie/project/rnaseq_test/*html .
# trim low quality sequences
ls *1*fastq.gz | sed 's/_.*.fastq.gz//' | parallel 'cutadapt -o {}.1.trim.fastq -p {}.2.trim.fastq --trim-n --minimum-length 50 --max-n 0 -q 30,30 {}_1.fastq.gz {}_2.fastq.gz 1>{}.trim.out'
# merge reads and further quality filter
ls *1.trim.fastq | sed 's/\..*trim.fastq//' | parallel 'pear -f {}.1.trim.fastq -r {}.2.trim.fastq -o {}.merge.fastq  -q 30 -j 4'
# All trimmed files were merged using cat and saved as project_merged_trimmed.fastq

#############
# rRNA FILTER
#############
# identify rRNA sequences
ls *.fastq.gz_trim.fastq | sed 's/.fastq.gz_trim.fastq//' | parallel 'sortmerna --ref set1-database.fasta --reads {} --fastx --workdir /home/prakrit/merged_fastq_files/sortmerna_{}' 
# filter rRNA sequences from dataset
cd sortmerna/out/
seqtk seq -a aligned.fastq > aligned.fasta 

# how many 16S sequences?
grep "^>" aligned.fasta -c 
#ERR348407: 6104949 
#ERR348408: 7901918 
#ERR348409: 5440810 
#ERR348410: 6576728 
#ERR348411: 5476873 
#ERR348412: 7595046 
#ERR348413: 7765705 
#ERR348414: 5318011 
#ERR348415: 6163219 
#ERR348416: 5354819 
#ERR348417: 5278476 
#ERR348418: 6819399
#ERR348419: 4990145 
#ERR348420: 4638782 
#ERR348421: 5078364 
#ERR348422: 5435179 
# All the aligned files were eventually merged using cat and saved as project_aligned.16s.ids
grep "^>" aligned.fasta | sed 's/>//' | awk '{print $1}' > aligned.16S.ids
cd ../..
seqtk subseq SRR1646851.merge.fastq.assembled.fastq sortmerna/out/aligned.16S.ids > filtered.16S.fastq
seqtk subseq ERR3484xx.fastq.gz_trim.fastq sortmerna_ERR3484xx/out/aligned.16S.ids > project_filtered.16S.fastq # I did this for each read individually, hence, "xx"
# remove 16S seqs from full dataset
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' SRR1646851.merge.fastq.assembled.fastq | grep -Ff sortmerna/out/aligned.16S.ids - | tr "\t" "\n" > filtered.fastq
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' ERR3484xx.fastq.gz_trim.fastq | grep -Ff sortmerna_ERR3484xx/out/aligned.16S.ids - | tr "\t" "\n" > ERR348412_filtered.fastq
# I did this for each read individually, hence, "xx"
###########
# ALIGNMENT
###########
# convert GFF to GTF format
gffread ALL_genomes.gff -T -o ALL_genomes.gtf
# generate reference index 
# NOTE:
# be careful to check no one is running something else when running this command or the computer will crash
# NOTE:
# this step takes a long time but only has to be run once
STAR --runMode genomeGenerate \
	--genomeFastaFiles ALL_genomes.fna  \
	--runThreadN 8 \
	--limitGenomeGenerateRAM 66959267424 \
	--sjdbGTFfile ALL_genomes.gtf \
	--genomeChrBinNbits 15
# map to reference genomes
STAR --runThreadN 4 \
	--genomeDir GenomeDir \
	--readFilesIn filtered.fastq \
	--outFileNamePrefix star-results \
	--outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped unmapped.bam \
	--quantMode TranscriptomeSAM GeneCounts \
	--alignIntronMax 1 \
	--chimOutType SeparateSAMold 

#################
# GET GENE COUNTS
#################
# convert bam to sam formatted file
samtools view -h -o \
	star-resultsAligned.toTranscriptome.out.sam \
	star-resultsAligned.toTranscriptome.out.bam
# what genes were identified?
wget http://www.homd.org/ftp/HOMD_prokka_genomes/tsv/ALL_genomes.tsv
# pull from locus id
grep -v "^@" star-resultsAligned.toTranscriptome.out.sam | awk -F"\t" '{print $3}' | while read line; do grep -w -m 1 $line ALL_genomes.tsv ; done > test.txt








## TO DO: Add alignment to human genome step to identify/filter human transcripts?
## TO DO: Taxonomic analysis using filtered 16S rRNA reads
## TO DO: Add in fungal/viral genomes for analysis


