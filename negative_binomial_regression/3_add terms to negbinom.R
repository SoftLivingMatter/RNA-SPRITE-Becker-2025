require(tidyverse)
require(dplyr)
require(ggplot2)
require(SparseM)
require(Matrix)
require(tibble)
require(foreach)
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::first)
conflicted::conflicts_prefer(dplyr::rename)
conflicted::conflicts_prefer(base::intersect)
conflicted::conflicts_prefer(base::setdiff)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::slice)

dir = "working_directory/prefix "

setwd("working_directory")

meta = read_csv(paste0(dir, "1 trans per geneID with meta.csv"))

edges1x = read_csv(paste0(dir, "20clustsizemax no_readcount_z dist_strand gxgpairwise.csv"))

#metadata for human
meta = read_csv("transcriptome_annotation/240321_human_AllGene_types ensembl meta GeneID_curate.csv.gz")

#metadata for rat
#meta = read_csv("transcriptome_annotation/240321_rat_AllGene_types ensembl meta GeneID_curate.csv.gz")

#Human
pairs2x_string = read_csv("transcriptome_annotation/24.07.29 Human stringdb phys scores over700 GeneIDcurate 2rows_per_pair only.csv")

pairs2x_DB = read_csv("transcriptome_annotation/24.08.08 CORUM and complex_portal GeneIDcurate pairs Human studies only 2rows_per_pair only.csv")

pairs2xH = read_csv("transcriptome_annotation/24.09.09 RNA_RNA_duplexes only PConly Hu_only pairs2x.csv")

pairs2xHM = read_csv("transcriptome_annotation/24.09.09 RNA_RNA_duplexes only PConly Hu and Ms to Hu pairs2x.csv")

#rat
#pairs2x_string = read_csv("transcriptome_annotation/24.07.28 Rat stringdb phys scores over700 GeneIDcurate 2rows_per_pair bymax_scores.csv")

#pairs2x_DB = read_csv("transcriptome_annotation/24.08.08 Rt CORUM and complex_portal GeneIDcurate pairs 2rows_per_pair only.csv")

#pairs2xM = read_csv("transcriptome_annotation/24.09.09 RNA_RNA_duplexes only PC_only Ms_only to Rt pairs2x.csv")

#pairs2xHM = read_csv("transcriptome_annotation/24.09.09 RNA_RNA_duplexes only PC_only Hu and Ms to Rt pairs2x.csv")



#adding meta to pairwise data
head(edges1x)

edges1x = edges1x %>% select(geneA, geneB, interactions, distance)

print("any(is.na(edges1x))")
any(is.na(edges1x))

meta = meta %>% mutate_if(is.numeric, ~replace(., is.na(.), 0))

to_join = meta %>% select(geneA, total_interactions)

# #Adding in count and other info for geneA
edges1x = left_join(edges1x, to_join, by = c("geneA")) %>% rename(total_interactionsA = total_interactions)

# #for geneB
edges1x = left_join(edges1x, to_join, by = c("geneB" = "geneA")) %>% rename(total_interactionsB = total_interactions)
head(edges1x)

edges1x = edges1x %>% select(geneA, geneB, interactions, distance, total_interactionsA, total_interactionsB)

All_interactions = sum(edges1x$interactions)
sum(edges1x$interactions)
sum(meta$total_interactions)/2

edges1x = edges1x %>% mutate(expected_prelim = (total_interactionsA*total_interactionsB)/All_interactions)

head(edges1x)

# add in genomic distance bins, this should be the only thing besides the files to change between datasets
edges1x = edges1x %>% mutate(dist_bin = ifelse(distance < 1e4,"kb_0to10", ifelse(distance < 5e4,"kb_10to50", ifelse(distance < 1e5, "kb_50to100", ifelse(distance < 5e5, "kb_100to500", "over_500kb")))))

unique(edges1x$dist_bin)

head(edges1x)

edges1x$dist_bin = factor(edges1x$dist_bin, level = c("over_500kb", "kb_0to10", "kb_10to50", "kb_50to100", "kb_100to500"))

bin_info = edges1x %>% group_by(dist_bin) %>% summarise(count = n())

bin_info

sum(bin_info$count)

rm(bin_info)

to_join = meta2 %>% group_by(Gene_stable_ID) %>% summarise(Strand = first(Strand))

edges1x = left_join(edges1x, to_join, by = c("geneA" = "Gene_stable_ID")) %>% rename(StrandA = Strand)

# #for geneB
edges1x = left_join(edges1x, to_join, by = c("geneB" = "Gene_stable_ID")) %>% rename(StrandB = Strand)

head(edges1x)

#seeing if strand matters
edges1x = edges1x %>% mutate(strand_cat = ifelse(dist_bin != "over_500kb", ifelse(StrandA == StrandB, "same", "diff"), ""),
                             dist_strand = paste0(dist_bin, strand_cat))
#flag
edges1x$dist_strand = factor(edges1x$dist_strand, level = c("over_500kb", "kb_0to10same", "kb_0to10diff", "kb_10to50same", "kb_10to50diff", "kb_50to100same", "kb_50to100diff", "kb_100to500same", "kb_100to500diff"))

unique(edges1x$dist_strand)
unique(edges1x$strand_cat)
unique(edges1x$dist_bin)

print("any(is.na(edges1x$expected_prelim))")
any(is.na(edges1x$expected_prelim))

print("any(is.na(edges1x$interactions))")
any(is.na(edges1x$interactions))

print("any(is.na(edges1x))")
any(is.na(edges1x))

edges1x = edges1x %>% select(geneA, geneB, interactions, expected_prelim, dist_strand) %>% arrange(geneA, geneB)


#adding complex info
#double checked this is every pair twice
colnames(pairs2x_string)
nrow(pairs2x_string)

temp = pairs2x_string %>% select(GeneID_curateA, GeneID_curateB) %>% distinct()
nrow(temp)

colnames(pairs2x_DB)
nrow(pairs2x_DB)

temp = pairs2x_DB %>% select(geneA, geneB) %>% distinct()
nrow(temp)

pairs2x = full_join(pairs2x_DB, pairs2x_string, by = c("geneA" = "GeneID_curateA", "geneB" = "GeneID_curateB"))

pairs2x = pairs2x %>% mutate(same_complex_CP = ifelse(is.na(Complex_ac), 0, 1),
                             same_complex_CORUM = ifelse(is.na(ComplexID), 0, 1),
                             phys_score = ifelse(is.na(phys_score), 0, phys_score),
                             phys_string = ifelse(phys_score > 700, 1, 0),
                             phys2 = ifelse((same_complex_CP + same_complex_CORUM + phys_string) > 1, 1, 0)) 

nrow(pairs2x)
pairs2x = pairs2x %>% select(geneA, geneB, phys2) %>% distinct()
nrow(pairs2x)

edges1x = left_join(edges1x, pairs2x, by = c("geneA", "geneB")) 

edges1x = edges1x %>% mutate_if(is.numeric, ~replace(., is.na(.), 0))

phys2_genes <- unique(c(edges1x %>% filter(phys2 == 1) %>% pull(geneA), edges1x %>% filter(phys2 == 1) %>% pull(geneB)))

edges1x = edges1x %>% mutate(any_complex_phys2 = ifelse(geneA %in% phys2_genes, ifelse(geneB %in% phys2_genes, 2, 1),  ifelse(geneB %in% phys2_genes, 1, 0)))

head(edges1x)

unique(edges1x$phys2)
unique(edges1x$any_complex_phys2)

print("any(is.na(edges1x))")
any(is.na(edges1x))

print("pairs with encoded proteins interacting")
sum(edges1x$phys2 > 0)

head(edges1x)




#RNA-RNA hybridization
pairs2xH = pairs2xH %>% rename(H_assays = assays)

pairs2xHM = full_join(pairs2xHM, pairs2xH, by = c("Human_gene_stable_IDA", "Human_gene_stable_IDB")) %>% 
  mutate(H_assays = ifelse(is.na(H_assays), 0, H_assays),
         M_assays = assays - H_assays) 

head(pairs2xHM)

edges1x = left_join(edges1x, pairs2xHM, by = c("geneA" = "Human_gene_stable_IDA", "geneB" = "Human_gene_stable_IDB")) 

edges1x = edges1x %>% mutate_if(is.numeric, ~replace(., is.na(.), 0))

head(edges1x)


HMgenes <- unique(c(edges1x %>% filter(assays > 0) %>% pull(geneA), edges1x %>% filter(assays > 0) %>% pull(geneB)))

edges1x = edges1x %>% mutate(HuMs_int = as.numeric(assays > 0),
                             any_HuMs_int = ifelse(geneA %in% HMgenes, ifelse(geneB %in% HMgenes, 2, 1),  ifelse(geneB %in% HMgenes, 1, 0)))


head(edges1x)

print("pairs with mRNA-mRNA duplex annotations")
sum(edges1x$HuMs_int > 0)

edges1x = edges1x %>% select(-H_assays, -M_assays, -assays)

unique(edges1x$HuMs_int)
unique(edges1x$any_HuMs_int)

print("any(is.na(edges1x))")
any(is.na(edges1x))
head(edges1x)


require(adaglm) # # devtools::install_github("davidaknowles/adaglm")

nb_glm_fit = adaglm(
  interactions ~ log(expected_prelim) + dist_strand + phys2 + any_complex_phys2, 
  data = edges1x,
  verbosity = 2, 
  batch_size = 100000,
  learning_rate = 0.01,
  loglik_tol = 0.01,
  epochs = 100)

print(summary(nb_glm_fit))

nb_glm_fit = adaglm(
  interactions ~ log(expected_prelim) + dist_strand + HuMs_int + any_HuMs_int, 
  data = edges1x,
  verbosity = 2, 
  batch_size = 100000,
  learning_rate = 0.01,
  loglik_tol = 0.01,
  epochs = 100)

print(summary(nb_glm_fit))