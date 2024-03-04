#!/usr/local/bin/Rscript

defaultW <- getOption("warn")
# message(defaultW)
options(warn = -1)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(epitools))
suppressPackageStartupMessages(library(operators))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(writexl))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library(memes))
suppressPackageStartupMessages(library(bioseq))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(gdata))

############################################REGIONS DEFINITION##############################################################################
p <- arg_parser("Somatic mutation rate analysis")
p <- add_argument(p, "--reg", help="file of regions to read and to edit",short="-r")
p <- add_argument(p, "--root1", help="output file root (default: root1 name of 'file' argument)",short="-b1", default = NULL)
args <- parse_args(p)

if(is.na(args$root1)){args$root1 <- gsub("_1kb_filtered.bed", "", basename(args$reg))}

#Import regions to manipulate in R
setwd("/Volumes/Enzo_lab/materiale_lab_Capasso/Lezioni_Bonfiglio/progetto_CRC/inhouse_public_EGA_TARGET_inhouse_integration/Mutation_rate/results/A3_new")

All_reg <- fread(args$reg, sep = "\t")

Encode_bed <- fread("/Volumes/Enzo_lab/materiale_lab_Capasso/Lezioni_Bonfiglio/progetto_CRC/inhouse_public_EGA_TARGET_inhouse_integration/Mutation_rate/regions/encodetfbsV3.bed", sep = "\t")

Problematic_reg <- fread("/Volumes/Enzo_lab/materiale_lab_Capasso/Lezioni_Bonfiglio/progetto_CRC/inhouse_public_EGA_TARGET_inhouse_integration/Mutation_rate/regions/wgEncodeCrgMapabilityAlign36mer_red.bed.gz", sep = '\t', header = FALSE)

#Set column names
setnames(All_reg , c("chr", "start", "end", "TFBS_ID"))
setnames(Encode_bed, c("chr", "start", "end", "width", "strand", "ENCODE_ID"))
#setnames(CRC_TFBS, c("chr", "start", "end", "TFBS_ID"))
setnames(Problematic_reg, c("chr", "start", "end"))

#Generate pos_ID for All_reg
All_reg$pos_ID <- paste(All_reg$chr, All_reg$start, All_reg$end, sep = "_")

#Genomic range object conversion
#All_reg
All_reg$start <- All_reg$start + 1
All_reg$chr <- paste("chr", All_reg$chr, sep = "")
All_reg.range <- makeGRangesFromDataFrame(All_reg, keep.extra.columns=TRUE)
#seqinfo(All_reg.range) <- Seqinfo(genome="hg19")
All_reg.range <- sort(All_reg.range)
All_reg.table <- as.data.table(All_reg.range)
table(All_reg.table$width)

#Check lenghts (2001)
check_lenght <- All_reg.range@ranges@width == 2001
table(check_lenght)

#Properly modify
All_reg.table <- as.data.table(All_reg.range)
table(All_reg.table$width)
All_reg.table$end <- ifelse(All_reg.table$width == 2000, All_reg.table$end + 1, All_reg.table$end)
All_reg.table$start <- ifelse(All_reg.table$width == 2002, All_reg.table$start + 1, All_reg.table$start)
All_reg.range <- makeGRangesFromDataFrame(All_reg.table, keep.extra.columns=TRUE)
#seqinfo(All_reg.range) <- Seqinfo(genome="hg19")
All_reg.range <- sort(All_reg.range)
All_reg.table <- as.data.table(All_reg.range)
setnames(All_reg.table, "seqnames", "chr")
table(All_reg.table$width)

#Encode_bed
Encode_bed$start <- Encode_bed$start + 1
Encode_bed$chr <- paste("chr", Encode_bed$chr, sep = "")
Encode.range <- makeGRangesFromDataFrame(Encode_bed, keep.extra.columns=TRUE)
seqinfo(Encode.range) <- Seqinfo(genome="hg19")
Encode.range <- sort(Encode.range)

#CRG low mappability regions
Problematic.range <- makeGRangesFromDataFrame(Problematic_reg, keep.extra.columns=TRUE)
seqinfo(Problematic.range) <- Seqinfo(genome="hg19")
Problematic.range <- sort(Problematic.range)

#Get core regions (size = 200, anchor = center) of All_reg.range
All_reg.CORE.range <- resize(All_reg.range, width = 200, fix = "center")

#Check lenghts (200)
check_lenght.CORE <- All_reg.CORE.range@ranges@width == 200
table(check_lenght.CORE)

#Get flanking regions (size = 900, START = TRUE) of All_reg.range
All_reg.FLANKING.start.range <- flank(All_reg.CORE.range, width = 900, TRUE)

#Get flanking regions (size = 900, START = FALSE) of All_reg.range
All_reg.FLANKING.end.range <- flank(All_reg.CORE.range, width = 900, FALSE)

#Check lenghts (900)
check_lenght.FLANK.start <- All_reg.FLANKING.start.range@ranges@width == 900
table(check_lenght.FLANK.start)

#Check lenghts (900)
check_lenght.FLANK.end <- All_reg.FLANKING.end.range@ranges@width == 900
table(check_lenght.FLANK.end)

#################################################################################################################################
#FILTERING REGIONS WITH ENCODE TFBSs IN FLANKING REGIONS
#OPERATIONS ON FLANKING START
#Find Overlap between All_reg.FLANKING.start.range and Encode.range
All_reg.FL.start.over.Encode.range <- mergeByOverlaps(All_reg.FLANKING.start.range, Encode.range, type = "end", maxgap = 900)
#Convert All_reg.FL.start.over.Encode.range to data.table
All_reg.FL.start.over.Encode.table <- as.data.table(All_reg.FL.start.over.Encode.range)
#Retain important columns
All_reg.FL.start.over.Encode.table <- All_reg.FL.start.over.Encode.table[, c("All_reg.FLANKING.start.range.seqnames", "All_reg.FLANKING.start.range.start",
                                                                             "All_reg.FLANKING.start.range.end", "TFBS_ID", "pos_ID", "Encode.range.seqnames", "Encode.range.start",
                                                                             "Encode.range.end", "ENCODE_ID")]
#Add diff.end column
All_reg.FL.start.over.Encode.table$diff.end <- All_reg.FL.start.over.Encode.table$All_reg.FLANKING.start.range.end - All_reg.FL.start.over.Encode.table$Encode.range.end
#Remove rows with negative diff.end
All_reg.FL.start.over.Encode.table <- All_reg.FL.start.over.Encode.table[All_reg.FL.start.over.Encode.table$diff.end >= 0,]
#Check on max diff.end (diff.end <= 900)
which(All_reg.FL.start.over.Encode.table$diff.end > 900)
which(All_reg.FL.start.over.Encode.table$diff.end == 0)

#OPERATIONS ON FLANKING END
#Find Overlap between All_reg.FLANKING.end.range and Encode.range
All_reg.FL.end.over.Encode.range <- mergeByOverlaps(All_reg.FLANKING.end.range , Encode.range, type = "start", maxgap = 900)
#Convert All_reg.FL.end.over.Encode.range to data.table
All_reg.FL.end.over.Encode.table <- as.data.table(All_reg.FL.end.over.Encode.range)
#Retain important columns
All_reg.FL.end.over.Encode.table <- All_reg.FL.end.over.Encode.table[, c("All_reg.FLANKING.end.range.seqnames", "All_reg.FLANKING.end.range.start",
                                                                         "All_reg.FLANKING.end.range.end", "TFBS_ID", "pos_ID", "Encode.range.seqnames", "Encode.range.start",
                                                                         "Encode.range.end", "ENCODE_ID")]
#Add diff.start column
All_reg.FL.end.over.Encode.table$diff.start <- All_reg.FL.end.over.Encode.table$All_reg.FLANKING.end.range.start - All_reg.FL.end.over.Encode.table$Encode.range.start
#Remove rows with positive diff.start
All_reg.FL.end.over.Encode.table <- All_reg.FL.end.over.Encode.table[All_reg.FL.end.over.Encode.table$diff.start <= 0,]
#Check on max diff.start (diff.start >= -900)
which(All_reg.FL.end.over.Encode.table $diff.start < -900)
which(All_reg.FL.end.over.Encode.table $diff.start == 0)

#Subsetting of concatenated regions non-overlapping with ENCODE TFBS in flanking regions.start
All_reg_filtered.bed <- All_reg.table[All_reg.table$pos_ID %!in% All_reg.FL.start.over.Encode.table$pos_ID]

#Subsetting of concatenated regions non-overlapping with other ENCODE TFBS in flanking regions.start and flanking regions.end
All_reg_filtered.bed  <- All_reg_filtered.bed[All_reg_filtered.bed$pos_ID %!in% All_reg.FL.end.over.Encode.table$pos_ID]
which(duplicated(All_reg_filtered.bed$TFBS_ID))

#Bed filtering
All_reg_filtered.bed <- All_reg_filtered.bed[All_reg_filtered.bed$chr != "chrY"]
All_reg_filtered.bed <- All_reg_filtered.bed[All_reg_filtered.bed$chr != "chrX"]
All_reg_filtered.bed <- All_reg_filtered.bed[All_reg_filtered.bed$chr != "chrMT"]
All_reg_filtered.range <- makeGRangesFromDataFrame(All_reg_filtered.bed , keep.extra.columns=TRUE)
#seqinfo(All_reg_filtered.range) <- Seqinfo(genome="hg19")
All_reg_filtered.range <- sort(All_reg_filtered.range)
which(duplicated(All_reg_filtered.range$pos_ID))
which(duplicated(All_reg_filtered.range$TFBS_ID))

#Filtering with CRG bed
All_reg_filtered.range <- subsetByOverlaps(All_reg_filtered.range, Problematic.range, invert = TRUE, minoverlap = 1)
All_reg_filtered.bed <- as.data.table(All_reg_filtered.range)


#Get human genome (hg19)
human_genome <- BSgenome.Hsapiens.UCSC.hg19

#Add +- 1bp for each sequence considered
All_reg_filtered_exp.bed <- All_reg_filtered.bed
All_reg_filtered_exp.bed$start <- All_reg_filtered_exp.bed$start - 1
All_reg_filtered_exp.bed$end <- All_reg_filtered_exp.bed$end + 1

#GENOMICRANGE OF THE SET OF SEQUENCES
All_reg_filtered_exp.range <- makeGRangesFromDataFrame(All_reg_filtered_exp.bed, keep.extra.columns=TRUE)
#seqinfo(All_reg_filtered_exp.range) <- Seqinfo(genome="hg19")
All_reg_filtered_exp.range <- sort(All_reg_filtered_exp.range)
All_reg_filtered_exp.bed.sorted <- as.data.table(All_reg_filtered_exp.range)

#GET SEQUENCES 
All_reg_filtered_exp.bed.sorted$sequences <- as.data.table(get_sequence(All_reg_filtered_exp.range, human_genome))
TFBS_ID_to_rm <- All_reg_filtered_exp.bed.sorted[which(All_reg_filtered_exp.bed.sorted$sequences %like% "N")]$TFBS_ID
All_reg_filtered_exp.bed.sorted <- All_reg_filtered_exp.bed.sorted[All_reg_filtered_exp.bed.sorted$TFBS_ID %!in% TFBS_ID_to_rm]
All_reg_filtered.bed <- All_reg_filtered.bed[All_reg_filtered.bed$TFBS_ID %!in% TFBS_ID_to_rm]

#GENOMICRANGE OF THE SET OF SEQUENCES
All_reg_filtered.range <- makeGRangesFromDataFrame(All_reg_filtered.bed, keep.extra.columns=TRUE)
#seqinfo(All_reg_filtered_exp.range) <- Seqinfo(genome="hg19")
All_reg_filtered.range <- sort(All_reg_filtered.range)


#GENOMICRANGE OF THE SET OF SEQUENCES
All_reg_filtered_exp.range <- makeGRangesFromDataFrame(All_reg_filtered_exp.bed.sorted, keep.extra.columns=TRUE)
#seqinfo(All_reg_filtered_exp.range) <- Seqinfo(genome="hg19")
All_reg_filtered_exp.range <- sort(All_reg_filtered_exp.range)
All_reg_filtered_exp.bed.sorted <- as.data.table(All_reg_filtered_exp.range)

######################################################################################################################################
#GENOMICRANGE OF THE FINAL MUTATION SET --> REMEMBER TO PROPERLY CHANGE "snvs.pass.bed" in paretheses !!!
#setwd("/Volumes/Enzo_lab/materiale_lab_Capasso/Lezioni_Bonfiglio/progetto_CRC/inhouse_public_EGA_TARGET_inhouse_integration/Mutation_rate/mutations")
snvs.pass.bed <- fread("/Volumes/Enzo_lab/materiale_lab_Capasso/Lezioni_Bonfiglio/progetto_CRC/inhouse_public_EGA_TARGET_inhouse_integration/Mutation_rate/mutations/ALL_SNVs.bed", sep = "\t", header = TRUE)
snvs.pass.range <- makeGRangesFromDataFrame(snvs.pass.bed, keep.extra.columns=TRUE)
seqinfo(snvs.pass.range) <- Seqinfo(genome="hg19")
snvs.pass.range <- sort(snvs.pass.range)

##########################################################OBSERVED MUTATION RATE##############################################################################

#Find Overlaps between regions of interest (filtered or not with ENCODE/CRC TFBS) and snvs.pass.range --> REMEMBER TO PROPERLY CHANGE THE INVESTIGATED SET OF REGIONS in paretheses !!!
All_reg_filtered.over.snvs.range <- mergeByOverlaps(All_reg_filtered.range, snvs.pass.range, ignore.strand = TRUE, minoverlap = 1)
All_reg_filtered.over.snvs.table <- as.data.table(All_reg_filtered.over.snvs.range)
All_reg_filtered.over.snvs.table <- All_reg_filtered.over.snvs.table[, c("All_reg_filtered.range.seqnames", "All_reg_filtered.range.start", "All_reg_filtered.range.end",
                                                                         "All_reg_filtered.range.width", "TFBS_ID","pos_ID", "snvs.pass.range.seqnames", "snvs.pass.range.start", "snvs.pass.range.end",
                                                                         "snvs.pass.range.width", "REF","ALT", "ID")]

length(which(duplicated(All_reg_filtered.over.snvs.table$TFBS_ID)))
length(which(duplicated(All_reg_filtered.over.snvs.table$ID)))

#Generate a copy of the table for exp.mut rate calculus and calculate the total number of observed mutation for each region in the set
All_reg_filtered.over.snvs.exp.table <- All_reg_filtered.over.snvs.table
mut_per_reg <- as.data.table(table(All_reg_filtered.over.snvs.exp.table$TFBS_ID), keep.rownames = TRUE)
setnames(mut_per_reg, c("TFBS_ID", "n_mut"))
All_reg_ID <- as.data.table(All_reg_filtered.bed$TFBS_ID)
setnames(All_reg_ID, c("TFBS_ID"))
All_reg_ID_w_mut <- merge.data.table(All_reg_ID, mut_per_reg, by = "TFBS_ID", all.x = TRUE)
All_reg_ID_w_mut[is.na(All_reg_ID_w_mut)] <- 0

#Create mut_pos ID
All_reg_filtered.over.snvs.table$mutation_pos <- (All_reg_filtered.over.snvs.table$snvs.pass.range.start - All_reg_filtered.over.snvs.table$All_reg_filtered.range.start) + 1

#Create chrom_mut_pos_ID
All_reg_filtered.over.snvs.table$chrom_mut_pos_ID <- paste(All_reg_filtered.over.snvs.table$ID, All_reg_filtered.over.snvs.table$mutation_pos, sep = "_")

#Subset of All_reg_filtered.over.snvs.table with duplicated chrom_mut_pos_ID
duplicated_chrom_mut_pos_ID <- All_reg_filtered.over.snvs.table$chrom_mut_pos_ID[which(duplicated(All_reg_filtered.over.snvs.table$chrom_mut_pos_ID))]
All_reg_filtered.over.snvs.dup_chrom_mut_pos_ID.table <- All_reg_filtered.over.snvs.table[All_reg_filtered.over.snvs.table$chrom_mut_pos_ID %in% duplicated_chrom_mut_pos_ID]

#Remove duplicate in chrom_mut_pos_ID
All_reg_filtered.over.snvs.table <- All_reg_filtered.over.snvs.table[!duplicated(All_reg_filtered.over.snvs.table$chrom_mut_pos_ID), ]

#Add class_reg to All_reg_w_OCR_filtered.over.snvs.table
All_reg_filtered.over.snvs.table$class_reg <- ifelse(All_reg_filtered.over.snvs.table$mutation_pos >= 901 & All_reg_filtered.over.snvs.table$mutation_pos <= 1101, "CORE", "FLANKING")

#Check on class_reg
length(which(All_reg_filtered.over.snvs.table[All_reg_filtered.over.snvs.table$class_reg == "CORE"]$mutation_pos < 901))
length(which(All_reg_filtered.over.snvs.table[All_reg_filtered.over.snvs.table$class_reg == "CORE"]$mutation_pos > 1101))
length(which(All_reg_filtered.over.snvs.table[All_reg_filtered.over.snvs.table$class_reg == "FLANKING"]$mutation_pos >= 901 & All_reg_filtered.over.snvs.table[All_reg_filtered.over.snvs.table$class_reg == "FLANKING"]$mutation_pos <= 1101))

#Duplication check
which(duplicated(All_reg_filtered.over.snvs.table$ID))
which(All_reg_filtered.over.snvs.table$mutation_pos < 0)
which(All_reg_filtered.over.snvs.table$mutation_pos > 2001)
length(which(duplicated(All_reg_filtered.over.snvs.table$mutation_pos)))
length(which(duplicated(All_reg_filtered.over.snvs.table$chrom_mut_pos_ID)))

Mut_count <- as.data.table(table(All_reg_filtered.over.snvs.table$mutation_pos))
setnames(Mut_count, c("SNV_pos", "SNV_count"))
Mut_count$SNV_pos <- as.numeric(Mut_count$SNV_pos)
Mut_per_pos <- data.table(Pos = c(1:2001))
Mut_per_pos <- merge.data.table(Mut_count,Mut_per_pos, by.x = "SNV_pos", by.y = "Pos", all.y = TRUE )
Mut_per_pos <- replace(Mut_per_pos, is.na(Mut_per_pos), 0)
Mut_per_pos$norm_mut_count <- round(Mut_per_pos$SNV_count/nrow(All_reg_filtered.bed), digits = 5)
setnames(Mut_per_pos, c("Mut_pos", "SNV_count", "norm_mut_count"))
Mut_per_pos$`distance from summit(bp)` <- c(-1000:1000)
tot_mut_mapped <- sum(Mut_per_pos$SNV_count)
Mut_per_pos$class_reg <- ifelse(Mut_per_pos$`distance from summit(bp)` %in% c(-100:100), "CORE", "FLANKING")
table(Mut_per_pos$class_reg)

Obs_mut_count_Core <- Mut_per_pos[, sum(SNV_count), by = class_reg]

#Import trinuc_prob_aggregate.txt
trinuc_prob_aggregate <- fread("/Volumes/Enzo_lab/materiale_lab_Capasso/Lezioni_Bonfiglio/progetto_CRC/inhouse_public_EGA_TARGET_inhouse_integration/Mutation_rate/mutations/trinuc_prob_aggregate.txt", sep = "\t", header = TRUE)

################################################Fl CALCULUS FOR EACH TRINUCLEOTIDE IN SEQUENCES
All_reg_seq <- dna(All_reg_filtered_exp.bed.sorted$sequences)

#Define cores
n.cores <- 5

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)

#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#Parallelize k-mer splitting
start_time <- Sys.time()
All_reg_splitted_seq <- foreach(
  i = seq_along(All_reg_seq), .combine = 'c') %dopar% {
    seq_split_kmer(All_reg_seq[i], k = 3)
  }
end_time <- Sys.time()

end_time - start_time
parallel::stopCluster(cl = my.cluster)

names(All_reg_splitted_seq) <- All_reg_filtered_exp.bed.sorted$TFBS_ID

#Convert list of splitted sequences in tables of splitted sequences
All_reg_splitted_seq <- lapply(All_reg_splitted_seq, function(x) {as.data.table(x)})
All_reg_splitted_seq <- lapply(All_reg_splitted_seq, function(x) {setnames(x, "trinuc")})
All_reg_splitted_seq_annot1 <- lapply(All_reg_splitted_seq, function(x) {x$pos <- c(1:2001);x})

start_time <- Sys.time()
All_reg_splitted_seq_annot2 <- lapply(All_reg_splitted_seq_annot1, function(x) {merge.data.table(x, trinuc_prob_aggregate, by = "trinuc", all.x = TRUE)})
end_time <- Sys.time()
end_time - start_time

#names(All_reg_splitted_seq_annot2) <- names(All_reg_splitted_seq_annot_sub)
#Order positions
All_reg_splitted_seq_annot3 <- lapply(All_reg_splitted_seq_annot2, function(x) {x[order(pos),]})

#Normalize P for each seqeunce
All_reg_splitted_seq_annot3 <- lapply(All_reg_splitted_seq_annot3, function(x) {x$P_norm <- x$P/sum(x$P);x})

#Extract vector of names
seq_names <- names(All_reg_splitted_seq_annot3)

#Add column with peak names
All_reg_splitted_seq_annot3 <- mapply(cbind, All_reg_splitted_seq_annot3, "TFBS_ID"=seq_names, SIMPLIFY=F)

#Add column with total mutation per sequence
All_reg_splitted_seq_annot4 <- lapply(All_reg_splitted_seq_annot3, function(x) {merge.data.table(x, All_reg_ID_w_mut, by = "TFBS_ID", all.x = TRUE)})

#Obtain list subset with regions containing at least one mutation
All_reg_splitted_seq_annot4_mut <- lapply(All_reg_splitted_seq_annot4, function(x) {x[x$n_mut > 0,]})
All_reg_splitted_seq_annot4_mut <- All_reg_splitted_seq_annot4_mut[sapply(All_reg_splitted_seq_annot4_mut, function(x) dim(x)[1]) > 0]

#Count regions with no mutations
seq_wo_mut <- length(All_reg_splitted_seq_annot4) - length(All_reg_splitted_seq_annot4_mut)

#Get list with only useful columns
All_reg_splitted_seq_annot5_mut <- lapply(All_reg_splitted_seq_annot4_mut, function(x) {x[, c(3,5,6)]})
keep(All_reg_splitted_seq_annot5_mut, All_reg_filtered.bed, All_reg_splitted_seq_annot4, Mut_per_pos, All_reg_filtered_exp.bed.sorted, Obs_mut_count_Core, n.cores, tot_mut_mapped, p, args, sure = TRUE)

#Permutation test
#Create an empty vector
output <- c()

#Defining an empty dataframe
permutations <- data.frame()

#For loop to append sample() results in output vector
start_time <- Sys.time()
for (j in 1:1000) {
for(i in 1:length(All_reg_splitted_seq_annot5_mut)) {
  new_value <- sample(All_reg_splitted_seq_annot5_mut[[i]]$pos, size = All_reg_splitted_seq_annot5_mut[[i]]$n_mut[1], prob = All_reg_splitted_seq_annot5_mut[[i]]$P_norm)
  output <- c(output, new_value)
}
permutations <- rbind(permutations, output)
output <- c()
message(i)}
end_time <- Sys.time()
end_time - start_time

#check perm_1 length
dim(permutations)
tot_mut_mapped

#rename columns
names(permutations)[1:ncol(permutations)]<-paste0('mut_', 1:ncol(permutations))

#annotate with the number of mutations in core regions
permutations$mut_in_core <- rowSums(permutations >= 901 & permutations <= 1101)

#annotate indicating permutations with more expected mutations than observedin core regions
permutations$more_exp_core <- ifelse(permutations$mut_in_core >= Obs_mut_count_Core[[2,2]], 1, 0)
outFileemp <- paste0(args$root1, "_emp_p.txt")
#Calculate empirival p-value
count <- which(permutations$more_exp_core == 1)
empirical_p <- length(count)/1000
class(empirical_p)
cat(empirical_p)
write.table(empirical_p, outFileemp, sep="\t", row.names=F, quote=F,col.names=T)

#Generate contingency table and compute fisher-square 
outFileFE <- paste0(args$root1, "_FE.txt")
outFileFISH <- paste0(args$root1, "_fisher_p.txt")
contingengy.data <- c(2*900*nrow(All_reg_filtered.bed) - sum(Mut_per_pos[class_reg == "FLANKING",]$SNV_count),
                      sum(Mut_per_pos[class_reg == "FLANKING",]$SNV_count),
                      201*nrow(All_reg_filtered.bed) - sum(Mut_per_pos[class_reg == "CORE",]$SNV_count),
                      sum(Mut_per_pos[class_reg == "CORE",]$SNV_count))
rnames <- c("FLANKING", "CORE")
cnames <- c("NOT-MUTATED", "MUTATED")
contingency_matrix <- matrix(contingengy.data,nrow=2,byrow=TRUE,dimnames=list(rnames,cnames))
odds <- oddsratio.fisher(contingency_matrix)
FE <- odds$measure
FISHER <- odds$p.value
write.table(FE, outFileFE, sep="\t",row.names=T, quote=F,col.names=T)
write.table(FISHER, outFileFISH, sep="\t",row.names=T, quote=F,col.names=T)

############################################PLOTTING MUTATION RATES##############################################################################
#Calculus of exp_SNV_count
All_reg_splitted_seq_annot4 <- lapply(All_reg_splitted_seq_annot4, function(x) {x$exp_SNV_count <- x$P_norm*x$n_mut;x})

#Extract list of tables with useful columns
All_reg_splitted_seq_annot5 <- lapply(All_reg_splitted_seq_annot4, function(x) {x[, c("pos", "exp_SNV_count", "TFBS_ID")]})

#Create a master table with all annotated trinuc of all sequences reporting names
All_trinuc_master <- bind_rows(All_reg_splitted_seq_annot5)

#Compute Exp_mut_per_pos
Exp_mut_per_pos <- All_trinuc_master[, sum(exp_SNV_count), by = pos]
setnames(Exp_mut_per_pos, c("Mut_pos", "exp_SNV_count"))

#Add column of Exp_SNV_count to Mut_per_pos table
Mut_per_pos <- merge.data.table(Mut_per_pos, Exp_mut_per_pos, by = "Mut_pos", all.x = TRUE)

#Add column of norm_Exp_SNV_count to Mut_per_pos table
Mut_per_pos$norm_Exp_SNV_count <- Mut_per_pos$exp_SNV_count/nrow(All_reg_filtered.bed)

#Smoothing curves
mod.smsp.obs <- smooth.spline(Mut_per_pos$`distance from summit(bp)`, Mut_per_pos$norm_mut_count)
Mut_per_pos$fit_obs <-  mod.smsp.obs$y
mod.smsp.exp <- smooth.spline(Mut_per_pos$`distance from summit(bp)`, Mut_per_pos$norm_Exp_SNV_count)
Mut_per_pos$fit_exp <-  mod.smsp.exp$y

#plotting mutation counts per nucleotide
outFileplot <- paste0(args$root1, ".png")
p1 <- ggplot(Mut_per_pos, aes(x=`distance from summit(bp)`)) + geom_line(aes(y= norm_mut_count), color = "red", size = 1, alpha = 0.5) + 
  geom_line(aes(y= fit_obs ), color = "red", size = 1) +
  geom_line(aes(y= norm_Exp_SNV_count), color = "black", size = 1, alpha = 0.5) +
  geom_line(aes(y= fit_exp), color = "black", size = 1) +
  theme_classic() + 
  theme(legend.title = element_text(size =20), axis.line = element_line(linetype = "solid", size = 1), 
        axis.text.x = element_text(size=20, angle = 45, hjust=1),axis.text.y = element_text(size=20), axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20), plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) + 
        scale_x_continuous(breaks = seq(min(Mut_per_pos$`distance from summit(bp)`), max(Mut_per_pos$`distance from summit(bp)`), by = 100))

  ggsave(outFileplot, dpi=1000, device = "png", path = "/Volumes/Enzo_lab/materiale_lab_Capasso/Lezioni_Bonfiglio/progetto_CRC/inhouse_public_EGA_TARGET_inhouse_integration/Mutation_rate/results/A3_new/", height = 6, width = 10, units = "in")

outFilereg <- paste0(args$root1, "_final_reg.bed")
#Export useful files for testing and plotting
#fwrite(Mut_per_pos, file = "/srv/ngs/analysis/aievola/Obs_Exp_Mut_rate_enrich_inhouse_R2/results/A2/SKNBE2C_PHOX2B_w_OCR_mut_table.txt", sep = "\t")
#fwrite(permutations, file = "/srv/ngs/analysis/aievola/Obs_Exp_Mut_rate_enrich_inhouse_R2/results/A2/SKNBE2C_PHOX2B_w_OCR_perm.txt", sep = "\t")
write.table(All_reg_filtered.bed, outFilereg, sep="\t", row.names=F, quote=F,col.names=T)
