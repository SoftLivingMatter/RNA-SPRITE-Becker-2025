require(tidyverse)
require(dplyr)
require(ggplot2)
require(SparseM)
require(Matrix)
require(tibble)

#metadata for human
meta = read_csv("transcriptome_annotation/240321_human_AllGene_types ensembl meta GeneID_curate.csv.gz")

#metadata for rat
#meta = read_csv("transcriptome_annotation/240321_rat_AllGene_types ensembl meta GeneID_curate.csv.gz")

setwd("working_directory")

#output folder for the pre-processing snakemake
snakemake_out = "snakemake_output/"

#minimum number of observsation per gene ID, set this based on your dataset
min_reads = 500

prefix = paste0("file_prefix", " minreads", min_reads)

pdf(paste0(prefix," R plots.pdf"))


#Step 1: count raw reads, collapse by GeneID_curate

#these are raw count outputs from the abundant ncRNA custom transcriptome
readcountsBT1 = read_tsv(paste0(snakemake_out,"alignments/raw_reads/OVERALL.aligned.Hsrepeat2.tsv"), col_names = c("treatment", "Transcript_stable_ID", "rawcountBT"))
head(readcountsBT1)

#these are raw count outputs from STAR whole transcriptome
readcountsBT2 = read_tsv(paste0(snakemake_out,"alignments/raw_reads/OVERALL.Aligned.toTranscriptome.out.tsv"), col_names = c("treatment", "Transcript_stable_ID", "rawcountBT"))
head(readcountsBT2)

rawcounts = rbind(readcountsBT1, readcountsBT2)

#taking the version numbers off of transcript ID
rawcounts = rawcounts %>% mutate(Transcript_stable_ID = sub("\\..*","", Transcript_stable_ID))

rawcounts = rawcounts %>% group_by(Transcript_stable_ID) %>% summarise(rawcountBT = sum(rawcountBT))

rawcounts = rawcounts %>% filter(Transcript_stable_ID != "*")

head(rawcounts)
print("all unique transcripts")
nrow(rawcounts)

to_join = meta %>% select(Transcript_stable_ID, GeneID_curate, Gene_type_curate, Gene_name_curate) %>% distinct()

rawcounts = left_join(rawcounts, to_join, by = "Transcript_stable_ID")
head(rawcounts)

#filling in missing info for the hand annotated ncRNA transcripts (human version)
rawcounts = rawcounts %>% mutate(GeneID_curate = ifelse(grepl("RNA",Transcript_stable_ID), sub("ENST", "ENSG", Transcript_stable_ID), GeneID_curate),
                                 Gene_name_curate = ifelse(grepl("RNA",Transcript_stable_ID), sub("ENST_", "", Transcript_stable_ID), Gene_name_curate),
                                 Gene_type_curate = ifelse(grepl("rRNA",Transcript_stable_ID), "rRNA", ifelse(grepl("tRNA", Transcript_stable_ID), "tRNA", Gene_type_curate)))


#filling in missing info for the hand annotated ncRNA transcripts (rat version)
# rawcounts = rawcounts %>% mutate(GeneID_curate = ifelse(grepl("RNA",Transcript_stable_ID), sub("ENSRNOT", "ENSRNOG", Transcript_stable_ID), GeneID_curate),
#                                  Gene_name_curate = ifelse(grepl("RNA",Transcript_stable_ID), sub("ENSRNOT", "", Transcript_stable_ID), Gene_name_curate),
#                                  Gene_type_curate = ifelse(grepl("rRNA",Transcript_stable_ID), "rRNA", ifelse(grepl("tRNA", Transcript_stable_ID), "tRNA",Gene_type_curate)))

print("any NA GeneID_curate?")
any(is.na(rawcounts$GeneID_curate))

#collapsing by GeneID_curate
rawcountsBG_DF = rawcounts %>% group_by(GeneID_curate) %>% summarise(rawcountBG = sum(rawcountBT), Gene_type_curate = first(Gene_type_curate), Gene_name_curate = first(Gene_name_curate))

head(rawcountsBG_DF)
print("all unique genes")
nrow(rawcountsBG_DF)

sum(rawcountsBG_DF$rawcountBG)

write_csv(rawcountsBG_DF, paste0(prefix, " raw read counts by GeneID.csv"))

temp = rawcountsBG_DF %>% group_by(Gene_type_curate) %>% summarise(GeneIDs = n(), total_count = sum(rawcountBG))
temp

quantile(rawcountsBG_DF$rawcountBG)
hist(log10(rawcountsBG_DF$rawcountBG), 200)

#removing mitochondrial-encoded genes. These are transcribed into 1 RNA molecule, if you want to leave them in, collapse into a single GeneID_curate
MTgenes = grep("^MT_.*", meta$Gene_name_curate, value = TRUE)

#rat version
#MTgenes = grep("^Mt_.*", meta$Gene_name_curate, value = TRUE)

print("Mitochondrial genes")
MTgenes

#filter to look at nuclear encoded mRNAs only. 
#If you want to look at other types of RNAs, then you can change the filtering here to do that
PC = rawcountsBG_DF %>% filter(Gene_type_curate == "protein_coding", !(Gene_name_curate %in% MTgenes))
quantile(PC$rawcountBG)
hist(log10(PC$rawcountBG), 200)

PCnomito_geneIDs = PC$GeneID_curate

to_join = rawcountsBG_DF %>% select(GeneID_curate, rawcountBG) %>% distinct()

head(to_join)

rawcounts = left_join(rawcounts, to_join, by = "GeneID_curate")

head(rawcounts)
print("these should be the same")
length(unique(rawcounts$Transcript_stable_ID))
nrow(rawcounts)

temp =  rawcounts %>% select(GeneID_curate, Gene_type_curate, Gene_name_curate, rawcountBG) %>% distinct()
print("these should be the same")
length(unique(temp$GeneID_curate))
nrow(temp)

write_csv(rawcounts, paste0(prefix, " raw read counts by GeneID and TranscriptID.csv"))

rm(readcountsBT1, readcountsBT2, PC, temp, to_join)


#Step 2: filter alignments and make a cluster by gene matrix

#recommend starting with a small cluster file to test code
ti = read_csv(paste0(snakemake_out,"/clusters/name/full_barcode_no_singles.csv.gz"))
head(ti)

#cleaning up transcript IDs human
ti = ti %>% mutate(Transcript_stable_ID = sub(":.*","",sub(".*_ENST","ENST",sub("\\..*","", alignment))))

#cleaning up transcript IDs rat
#ti = ti %>% mutate(Transcript_stable_ID = sub(":.*","",sub(".*_ENSRNOT","ENSRNOT",sub("\\..*","", alignment))))

to_join = rawcounts %>% select(Transcript_stable_ID, GeneID_curate) %>% distinct()

ti = left_join(ti, to_join, by = "Transcript_stable_ID")

head(ti)
print("any NA GeneID_curate?")
any(is.na(ti$GeneID_curate))

ti = ti %>% select(barcode, GeneID_curate) %>% distinct()

ti = ti %>% group_by(barcode) %>% mutate(cluster_size = n()) %>% ungroup() %>% filter(cluster_size > 1)

#grouping by cluster to look at distribution of cluster sizes
bycluster = ti %>% group_by(barcode) %>% summarise(cluster_size = cluster_size[1])

bycluster = bycluster %>% mutate(full_barcode = barcode) %>% separate(barcode, into = paste0("tag",1:5), sep = "[.]")

bycluster = bycluster %>% mutate(sonication = sub("ODD_","", sub("_rep.*", "", tag4)), rep = sub("_.*","", sub(".*Cov_","", tag4)))

head(bycluster)

write_csv(bycluster, paste0(prefix," ByCluster.csv"))

rm(bycluster)

to_join = rawcounts %>% select(GeneID_curate, Gene_type_curate, rawcountBG) %>% distinct()

print("these should be the same")
length(unique(to_join$GeneID_curate))
nrow(to_join)

print("Any NA rawcountBG?")
any(is.na(to_join$rawcountBG))

ti = left_join(ti, to_join, by = c("GeneID_curate"))
head(ti)

unique(ti$Gene_type_curate)

#cluster size max cutoff calculated based on all possible alignments (including low abundance alignments, excluding alignments sharing a GeneID_curate)
#doing a cluster size filter to limit the influence of any single cluster on the dataset. We were worried large clusters might be debris. 
ti = ti %>% filter(cluster_size <= 20, rawcountBG >= min_reads)

#filtered gene types
unique(ti$Gene_type_curate)

nrow(ti)
#getting rid of clusters that are likely from the nucleus (have nuclear alignments) human
nuc = ti %>% filter(Gene_type_curate %in% c("snRNA", "snoRNA", "scaRNA") | grepl("ITS", GeneID_curate) | grepl("ETS", GeneID_curate) | GeneID_curate %in% c("ENSG00000229807", "ENSG00000277027", "ENSG00000245532", "ENSG00000241743", "ENSG00000251562")) %>% select(barcode) %>% distinct()

#getting rid of clusters that are likely from the nucleus (have nuclear alignments) rat
#nuc = ti %>% filter(Gene_type_curate %in% c("snRNA", "snoRNA", "scaRNA") | grepl("ITS", GeneID_curate) | GeneID_curate == "ENSRNOG00000057527") %>% select(barcode) %>% distinct()

nrow(nuc)
nuc_clusters = nuc$barcode

#if you want to keep clusters with nuclear alignments then don't do this filter
ti = ti %>% filter(!(barcode %in% nuc_clusters))

rm(nuc_clusters, nuc)
head(ti)
nrow(ti)

#protein coding interactions only, change the filter above in "PC =" if you want to look at more than just mRNA
ti = ti %>% filter(GeneID_curate %in% PCnomito_geneIDs) %>% group_by(barcode) %>% mutate(cluster_size = n()) %>% ungroup() %>% filter(cluster_size > 1)

PCbarcodes = ti %>% select(barcode) %>% distinct()

head(PCbarcodes)

write_csv(PCbarcodes, paste0(prefix," PC barcode.csv"))

#making the cluster by gene matrix
i_factor = factor(ti$barcode)
j_factor = factor(ti$GeneID_curate)

clu_gene_curate = sparseMatrix(
  i = as.numeric(i_factor),
  j = as.numeric(j_factor),
  x = 1, #no normalization for cluster size
  dimnames = list(levels(i_factor), levels(j_factor))
)

saveRDS(clu_gene_curate, file = paste0(prefix, " 20clustsizemax clustersxgeneID.rds"))

dim(clu_gene_curate)

rm(ti, to_join, i_factor, j_factor)



# #Step 3 map alignments in mRNA-mRNA clusters to transcript position bins to find a single transcript sequence per GeneID_curate that has the best coverage across the entire sequence
ti = read_csv(paste0(snakemake_out, "clusters/position/full_barcode_no_singles.csv.gz"))
head(ti)

# #this is from the gene level filtering already done in step 1
head(PCbarcodes)
goodbarcodes = PCbarcodes$barcode

ti = ti %>% filter(barcode %in% goodbarcodes)
rm(PCbarcodes)
rm(goodbarcodes)

# #cleaning up transcript ids human
ti = ti %>% mutate(Transcript_stable_ID = sub(":.*","",sub(".*_ENST","ENST",sub("\\..*","", alignment))))

# #cleaning up transcript ids rat
#ti = ti %>% mutate(Transcript_stable_ID = sub(":.*","",sub(".*_ENSRNOT","ENSRNOT",sub("\\..*","", alignment))))

head(ti)

head(rawcounts)
to_join = rawcounts %>% select(Transcript_stable_ID, GeneID_curate, Gene_type_curate, rawcountBG) %>% distinct()

ti = left_join(ti, to_join, by = "Transcript_stable_ID")
head(ti)

# #cluster size max cutoff calculated based on all possible alignments (including low abundance alignments)
ti = ti %>% filter(GeneID_curate %in% PCnomito_geneIDs, rawcountBG >= min_reads) %>% select(-Gene_type_curate, -rawcountBG)

to_join = meta %>% select(Transcript_stable_ID, Transcript_length, Gene_type) %>% distinct()

ti = left_join(ti, to_join, by = "Transcript_stable_ID")
head(ti)

ti = ti %>% filter(Gene_type == "protein_coding")

ti = ti %>% mutate(position_start = sub("-.*", "", sub(".*:","", alignment)), position_end = gsub(".*-","", alignment))
ti$position_start = as.numeric(ti$position_start)
ti$position_end = as.numeric(ti$position_end)

# #break the transcript length into 10 bins, round up to the nearest integer
ti = ti %>% mutate(position_mid = (position_start+position_end)/2, positionbin = ceiling(position_mid/(Transcript_length/10)), position_norm = (position_mid/Transcript_length)*100)

head(ti)

#looking at coverage by calculating the standard deviation of the position of alignments across the length of each transcript sequence
bins = ti %>% group_by(Transcript_stable_ID) %>% 
  mutate(transcount = n(), 
         norm_position_sd = sd(position_norm, na.rm =T), 
         norm_position_mean = mean(position_norm, na.rm =T)) %>% 
  ungroup() %>% 
  group_by(Transcript_stable_ID, positionbin) %>% 
  summarise(Transcript_length = first(Transcript_length), 
            norm_position_sd = first(norm_position_sd), 
            norm_position_mean = first(norm_position_mean), 
            transcount = first(transcount), 
            bincount = n(), 
            count_prop = bincount/transcount) %>% 
  ungroup()

head(bins)

bins$positionbin = as.factor(bins$positionbin)

#You can use this dataframe to look at coverage across the length of transcripts
write_csv(bins, paste0(prefix," PConly alignment by position bin.csv"))

bins = bins %>% group_by(Transcript_stable_ID) %>% mutate(max_count_prop = max(count_prop, na.rm = T), max_bin = ifelse(count_prop == max_count_prop, T, F)) %>% ungroup()

BT = bins %>% filter(max_bin == T) %>% group_by(Transcript_stable_ID) %>% slice_sample(n = 1) %>% ungroup()

rm(bins)

to_join = rawcounts %>% filter(rawcountBG >= min_reads, GeneID_curate %in% PCnomito_geneIDs)

BT = left_join(to_join, BT, by = "Transcript_stable_ID")

BT$norm_position_sd[is.na(BT$norm_position_sd)] <- 0
BT$transcount[is.na(BT$transcount)] <- 0
BT$max_count_prop[is.na(BT$max_count_prop)] <- 0

BT = BT %>% group_by(GeneID_curate) %>% mutate(genecount = sum(transcount), max_transcount = max(transcount), BT_over_BG = transcount/genecount, maxBT_over_BG = max(BT_over_BG), max_norm_position_sd = max(norm_position_sd, na.rm = T), minBG_maxBT_count_prop = min(max_count_prop)) %>% ungroup()

head(BT)
print("these should be the same")
nrow(BT)
length(unique(BT$Transcript_stable_ID))

ggplot(BT, aes(positionbin, norm_position_sd)) + geom_boxplot()
ggplot(BT, aes(positionbin, count_prop)) + geom_boxplot()

TscriptsBG = BT %>% group_by(GeneID_curate) %>% mutate(max_norm_position_sd = max(norm_position_sd), keep = ifelse(norm_position_sd < 15, ifelse(norm_position_sd == max_norm_position_sd, T, F), T)) %>% ungroup() %>% filter(keep == T) %>% group_by(GeneID_curate) %>% slice_max(transcount) %>% slice_max(rawcountBT) %>% slice_sample(n = 1) %>% ungroup()

head(TscriptsBG)
print("these should be the same")
length(unique(TscriptsBG$GeneID_curate))
nrow(TscriptsBG)
length(unique(BT$GeneID_curate)) - length(unique(TscriptsBG$GeneID_curate))

hist(TscriptsBG$norm_position_sd)
hist(TscriptsBG$max_count_prop)
hist(log10(TscriptsBG$transcount))
hist(log10(TscriptsBG$genecount))
ggplot(TscriptsBG, aes(max_count_prop, norm_position_sd)) + geom_point(alpha = 0.05)

#This file will have the resulting 1 transcript sequence per gene ID
write_csv(TscriptsBG, paste0(prefix, " PConly 1trans per GeneID with coverage info.csv"))




#Step4: Negative binomial regression with pairwise results table

#creating a gene by gene interaction table
interaction_counts = t(clu_gene_curate) %*% clu_gene_curate

#converting interaction counts to a tibble
interaction_counts_tib = as_tibble(as(interaction_counts,"matrix"))
interaction_counts_tib$geneA = rownames(interaction_counts)

head(interaction_counts_tib)

print("clusters x genes")
dim(interaction_counts_tib)

rm(clu_gene_curate)

#making a table with pairwise interaction counts
#setting self interaction to zero
edges2x = interaction_counts_tib %>% pivot_longer(-geneA,names_to = "geneB", values_to = "interactions")
edges2x = edges2x %>% mutate(interactions = ifelse(geneA == geneB, 0, interactions))

head(edges2x)
length(unique(edges2x$geneA))

length(unique(edges2x$geneB))

rm(interaction_counts_tib, interaction_counts)

#getting total interactions per gene
ByGeneA = edges2x %>% group_by(geneA) %>% summarise(total_interactionsA = sum(interactions))

head(ByGeneA)
print("total genes filtered")
nrow(ByGeneA)

ByGeneA = left_join(ByGeneA, rawcountsBG_DF, by = c("geneA" = "GeneID_curate")) %>% rename(total_interactions = total_interactionsA)
head(ByGeneA)

nrow(ByGeneA)

write_csv(ByGeneA, paste0(prefix, " interactions and reads byGeneID PConly.csv"))

rm(rawcountsBG_DF)

#each gene pair seen once, no self interaction
edges1x = edges2x %>% mutate(geneA_factor = factor(geneA))
edges1x = edges1x %>% mutate(geneB_factor = factor(geneB, levels(edges1x$geneA_factor))) %>%
  filter(as.numeric(geneA_factor) > as.numeric(geneB_factor)) %>%
  filter(geneA != geneB) %>%
  select(-geneA_factor, -geneB_factor)

nrow(edges1x)
#n*(n-1)/2
nrow(edges1x)/nrow(edges2x)
#~0.5
length(unique(edges1x$geneA))
#first geneID is only geneB, last is only geneA
length(unique(edges1x$geneB))

setdiff(unique(edges1x$geneA), unique(edges1x$geneB))
setdiff(unique(edges1x$geneB), unique(edges1x$geneA))
head(edges1x)

rm(edges2x)

#Adding in count and other info for geneA
edges1x = left_join(edges1x, ByGeneA, by = c("geneA")) %>% rename(total_interactionsA = total_interactions, read_countA = rawcountBG)

head(edges1x)

#for geneB
edges1x = left_join(edges1x, ByGeneA, by = c("geneB" = "geneA")) %>% rename(total_interactionsB = total_interactions, read_countB = rawcountBG)

head(edges1x)

test = edges1x %>% filter(geneA == geneB)
sum(test$interactions)

#interactions in the whole dataset
All_interactions = sum(ByGeneA$total_interactions)/2
All_interactions
sum(edges1x$interactions)

rm(ByGeneA)

#z transform of log read counts so it's easier for the fitting algorithm
zscore <- function(x) { (x - mean(x, na.rm=T) ) / sd(x, na.rm=T) }

edges1x = edges1x %>% mutate(expected_prelim = (total_interactionsA*total_interactionsB)/All_interactions, read_count_z = zscore(log(read_countA)+log(read_countB)))

head(edges1x)

to_join = meta %>% group_by(GeneID_curate) %>% summarise(Gene_name_curate = first(Gene_name_curate), Chromosome_scaffold_name = first(Chromosome_scaffold_name), Gene_start_bp = first(Gene_start_bp), Gene_end_bp = first(Gene_end_bp), Strand = first(Strand))
nrow(to_join)
length(unique(to_join$GeneID_curate))

head(to_join)

edges1x = left_join(edges1x, to_join, by = c("geneA" = "GeneID_curate")) %>% rename(Gene_nameA = Gene_name_curate, ChromoA = Chromosome_scaffold_name, Gene_start_bpA = Gene_start_bp, Gene_end_bpA = Gene_end_bp, StrandA = Strand)

edges1x = left_join(edges1x, to_join, by = c("geneB" = "GeneID_curate")) %>% rename(Gene_nameB = Gene_name_curate, ChromoB = Chromosome_scaffold_name, Gene_start_bpB = Gene_start_bp, Gene_end_bpB = Gene_end_bp, StrandB = Strand)

head(edges1x)

rm(to_join)

#calculating distance from center point of genes
edges1x = edges1x %>% mutate(distance = ifelse(ChromoA == ChromoB, abs(((Gene_start_bpA+Gene_end_bpA)/2) - ((Gene_start_bpB+Gene_end_bpB)/2)), Inf))

edges1x = edges1x %>% mutate(distance = ifelse(is.na(distance), Inf, distance))

min(edges1x$distance)
max(edges1x$distance)

head(edges1x)

test = edges1x %>% filter(is.na(distance))
nrow(test)

edges1x = edges1x %>% mutate(dist_bin = ifelse(distance < 1e4,"kb_0to10", ifelse(distance < 5e4,"kb_10to50", ifelse(distance < 1e5, "kb_50to100", ifelse(distance < 5e5, "kb_100to500", "over_500kb")))))

unique(edges1x$dist_bin)

head(edges1x)

edges1x$dist_bin = factor(edges1x$dist_bin, level = c("over_500kb", "kb_0to10", "kb_10to50", "kb_50to100", "kb_100to500"))

bin_info = edges1x %>% group_by(dist_bin) %>% summarise(count = n())

bin_info

sum(bin_info$count)

rm(bin_info)

#seeing if strand matters
edges1x = edges1x %>% mutate(strand_cat = ifelse(dist_bin != "over_500kb", ifelse(StrandA == StrandB, "same", "diff"), ""),
                             dist_strand = paste0(dist_bin, strand_cat))

#flag
edges1x$dist_strand = factor(edges1x$dist_strand, level = c("over_500kb", "kb_0to10same", "kb_0to10diff", "kb_10to50same", "kb_10to50diff", "kb_50to100same", "kb_50to100diff", "kb_100to500same", "kb_100to500diff"))


unique(edges1x$dist_strand)
print("dist_strand NAs")
any(is.na(edges1x$dist_strand))

test = edges1x %>% filter(!(dist_bin %in% c("diff_chromos", "over_1MB", "over_500kb")))

sum(test$interactions)
sum(edges1x$interactions)

print("proportion interactions within 500kb")
sum(test$interactions)/sum(edges1x$interactions)


require(adaglm) #devtools::install_github("davidaknowles/adaglm")

#this is the base model, it takes into account total interactions for both mRNAs and their genomic distance
nb_glm_fit = adaglm(
  interactions ~ log(expected_prelim) + dist_strand,
  data = edges1x,
  verbosity = 2,
  batch_size = 100000,
  learning_rate = 0.01,
  loglik_tol = 0.01,
  epochs = 100)

summary(nb_glm_fit)

head(edges1x)

#pearson residuals, difference between observed and expected
PR = residuals(nb_glm_fit)
hist(asinh(PR), 200)
length(PR)

#Expected values from the full model
expected = fitted(nb_glm_fit)
hist(log(expected))
length(expected)

#observed values 
length(edges1x$interactions)
hist(log10(edges1x$interactions))

#Adding more info to dataframe, removing uneeded
pvals = pointwise_pvalues(nb_glm_fit, edges1x$interactions)
pvals = as.tibble(pvals)

to_save = edges1x
to_save$pval = pvals$p_two_sided
hist(-log10(to_save$pval))

to_save$pearson_resid = residuals(nb_glm_fit)[,1]
to_save$expected = fitted(nb_glm_fit)[,1]

to_save = to_save %>% mutate(ORint = interactions/expected_prelim, 
                             OR_NB = interactions/expected, 
                             pval_adj = p.adjust(pval, method = "BH"),
                             zscore = ifelse(pearson_resid > 0, qnorm(1 - (pval/2)), qnorm(pval/2))) %>% mutate(
                             zscore = ifelse(is.infinite(zscore), max(zscore[!(is.infinite(zscore))]), zscore))
hist(log10(to_save$OR_NB))

sample = sample_n(to_save, 100000)

to_save = to_save %>% select(geneA, geneB, interactions, pearson_resid, pval, zscore, pval_adj, OR_NB, ORint, distance)

head(to_save)

write_csv(to_save, paste0(prefix, " 20clustsizemax no_readcount_z dist_strand gxgpairwise.csv"))

#this is an alternative model that also takes into account total observations per gene.

# rm(nb_glm_fit)
# 
# nb_glm_fit = adaglm(
#   interactions ~ log(expected_prelim) + dist_strand + read_count_z,
#   data = edges1x,
#   verbosity = 2,
#   batch_size = 100000,
#   learning_rate = 0.01,
#   loglik_tol = 0.01,
#   epochs = 100)
# #
# summary(nb_glm_fit)
# 
# head(edges1x)
# 
# #pearson residuals, difference between observed and expected
# PR = residuals(nb_glm_fit)
# hist(asinh(PR), 200)
# length(PR)
# 
# #Expected values from the full model
# expected = fitted(nb_glm_fit)
# hist(log(expected))
# length(expected)
# 
# #observed values 
# length(edges1x$interactions)
# hist(log10(edges1x$interactions))
# 
# #Adding more info to dataframe, removing uneeded
# pvals = pointwise_pvalues(nb_glm_fit, edges1x$interactions)
# pvals = as.tibble(pvals)
# 
# to_save = edges1x
# to_save$pval = pvals$p_two_sided
# hist(-log10(to_save$pval))
# 
# to_save$pearson_resid = residuals(nb_glm_fit)[,1]
# to_save$expected = fitted(nb_glm_fit)[,1]
# 
# to_save = to_save %>% mutate(ORint = interactions/expected_prelim, 
#                              OR_NB = interactions/expected, 
#                              pval_adj = p.adjust(pval, method = "BH"),
#                              zscore = ifelse(pearson_resid > 0, qnorm(1 - (pval/2)), qnorm(pval/2)), 
#                              zscore = ifelse(is.infinite(zscore), max(zscore[!(is.infinite(zscore))]), zscore))
# hist(log10(to_save$OR_NB))
# 
# #sample = sample_n(to_save, 100000)
# 
# to_save = to_save %>% select(geneA, geneB, interactions, pearson_resid, pval, zscore, pval_adj, OR_NB, ORint, distance)
# 
# head(to_save)
# 
# write_csv(to_save, paste0(prefix, " 20clustsizemax w_readcount_z dist_strand gxgpairwise.csv"))

rm(nb_glm_fit, edges1x)

ggplot(sample, aes(asinh(expected), asinh(interactions))) + geom_point(alpha = 0.2) + theme_classic()

ggplot(sample, aes(asinh(expected), asinh(expected_prelim))) + geom_point(alpha = 0.2) + theme_classic()

ggplot(sample, aes(log10(read_countA), log10(total_interactionsA))) + geom_point(alpha = 0.01) + theme_classic()

ggplot(sample, aes((log(total_interactionsA)+log(total_interactionsB)), asinh(pearson_resid))) + geom_point(alpha = 0.05) + theme_classic()

ggplot(sample, aes(asinh(OR_NB), asinh(pearson_resid))) + geom_point(alpha = 0.05) + theme_classic()

ggplot(sample, aes(log10(read_countA), asinh(pearson_resid))) + geom_point(alpha = 0.05) + theme_classic()

ggplot(sample, aes(log10(total_interactionsA), asinh(pearson_resid))) + geom_point(alpha = 0.05) + theme_classic()

ggplot(sample, aes(asinh(pearson_resid), zscore)) + geom_point(alpha = 0.05) + theme_classic()

ggplot(sample, aes(asinh(interactions), zscore)) + geom_point(alpha = 0.05) + theme_classic()

ggplot(sample, aes(asinh(read_countA), zscore)) + geom_point(alpha = 0.05) + theme_classic()

rm(sample)

# #close the pdf device
dev.off()
