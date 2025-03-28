require(SparseM)
require(Matrix)
require(tibble)
require(ggplot2)
require(tidyverse)
require(dplyr)

#if sequences have strong homology, it's hard to tell if reads with the same barcode were from the same RNA molecule or from distinct molecules. 
#we aggregated observations into a single gene ID if the overlap score was > 50
#you only need to check homology once for each sequence pair. We kept the results from each alignment test, and filtered sig1x to only test for pairs that had not previously been tested. 

dir = "working_directory/prefix "

edges1x = read_csv(paste0(dir,"20clustsizemax no_readcount_z dist_strand gxgpairwise.csv"))

head(edges1x)

#these are the pairs you are testing for homology, we check everything with zscore > 3
sig1x = edges1x %>% filter(zscore > 4, pearson_resid > 0)

nrow(sig1x)
min(sig1x$interactions)
min(sig1x$pearson_resid)

#matches each gene ID with a transcript ID
transcripts = read_csv(paste0(dir,"PConly 1trans per GeneID with coverage info.csv"))
head(transcripts)

to_join = transcripts %>% select(GeneID_curate, Transcript_stable_ID, Gene_name_curate)

colnames(to_join) = paste0(colnames(to_join), "_A")

sig1x = left_join(sig1x, to_join, by = c("geneA" = "GeneID_curate_A"))

to_join = transcripts %>% select(GeneID_curate, Transcript_stable_ID, Gene_name_curate)

colnames(to_join) = paste0(colnames(to_join), "_B")

sig1x = left_join(sig1x, to_join, by = c("geneB" = "GeneID_curate_B"))

head(sig1x)

max(sig1x$interactions)
max(sig1x$pearson_resid)
min(sig1x$pval_adj)

trans <- unique(c(sig1x$Transcript_stable_ID_A, sig1x$Transcript_stable_ID_B))

#human, this file is too big for me to provide. Download the sequences from ensembl
#http://useast.ensembl.org/info/data/ftp/index.html
sequences = read_csv("Homo_sapiens.GRCh38.cdnaANDncrna sequences by transcriptID.csv.gz")

#rat, this file is too big for me to provide. Download the sequences from ensembl
#sequences = read_csv("Rattus_norvegicus.mRatBN7.2cdnaANDncrna sequences by transcriptID.csv.gz")

head(sequences)

sequences = sequences %>% filter(Transcript_stable_ID %in% trans) %>% mutate(transcript_length = nchar(sequence)) 

min(sequences$transcript_length)
max(sequences$transcript_length)

sequences = sequences %>% filter(transcript_length < 30000, transcript_length > 10)

seqcheck = sequences$Transcript_stable_ID

sig1x = sig1x %>% filter(Transcript_stable_ID_A %in% seqcheck, Transcript_stable_ID_B %in% seqcheck) %>% arrange(desc(zscore))
nrow(sig1x) 

results <- data.frame()

require("Biostrings")
for(i in 1:nrow(sig1x)){
  
  selected_row <- sig1x[i, ]
  
  Transcript_stable_ID_A = selected_row$Transcript_stable_ID_A
  Transcript_stable_ID_B = selected_row$Transcript_stable_ID_B
  
  temp = sequences %>% filter(Transcript_stable_ID == Transcript_stable_ID_A)
  sequenceA = temp$sequence
  
  temp = sequences %>% filter(Transcript_stable_ID == Transcript_stable_ID_B)
  sequenceB = temp$sequence
  
  alignment <- pairwiseAlignment(pattern = sequenceA, subject = sequenceB, type = "overlap")
  selected_row$align_score = score(alignment)
  
  results <- bind_rows(results, selected_row)
}

head(results)

write_csv(results, "overlap homology results.csv")