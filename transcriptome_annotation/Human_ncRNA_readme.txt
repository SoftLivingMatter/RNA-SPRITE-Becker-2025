Grabbing snRNA and snoRNA sequences from ensembl

FastaToTbl /projects/BRANGWYNNE/tools/genomics/transcriptomes/bowtie/Human_cDNA_and_ncRNA/Homo_sapiens.GRCh38.ncrna.fa | grep -E 'snoRNA|snRNA' | TblToFasta>Homo_sapiens.GRCh38.snRNAsnoRNA.fa

rRNA from Sofi’s curated fasta
/projects/BRANGWYNNE/tools/genomics/transcriptomes/bowtie/human_custom_transcripts/Hs_repetitive_sequences_20200109_SQcurate.fa

Renamed sequences so I can grab them more easily:
ENST_5ETS_45SrRNA
ENST_ITS1_45SrRNA
ENST_ITS2_45SrRNA
ENST_3ETS_45SrRNA
ENST_5_8SrRNA
ENST_18SrRNA
ENST_28SrRNA
ENST_5SrRNA

tRNA
From database http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi38/Hsapi38-seq.html, high confidence mature sequences

There are too many for me to edit the names, just going to do this later in R

New Bowtie library with abundant, repetitive RNAs for first round mapping
Conda activate myenv
bowtie2-build Homo_sapiens.GRCh38.snRNAsnoRNArRNAtRNA.fa Homo_sapiens.GRCh38.snRNAsnoRNArRNAtRNA

