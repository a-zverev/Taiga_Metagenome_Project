# In this step, we perform Illumina 16S amplicone processing, and export and 
# visualization of 16s sequences, extracted from ONP assembles.


# Set libraries, path to references and working directory
library(dada2)
library(Biostrings)
library(DECIPHER)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(plyr)

setwd('/home/alexey/BI/Praktikum/Taiga/16s_rdna/')

path_train_set = "/home/alexey/tax_n_refs/silva_nr_v132_train_set.fa.gz"
path_train_set_species = "/home/alexey/tax_n_refs/silva_species_assignment_v132.fa.gz"
path <- 'raw/'


# dada2 raw pipeline for single-end reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, trimLeft=19, maxN=0, maxEE=2, rm.phix=TRUE, compress=TRUE, multithread=8)
errF <- learnErrors(filtFs, multithread=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
names(derepFs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
table(nchar(getSequences(seqtab)))
getN <- function(x) sum(getUniques(x))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=8, verbose=TRUE)

# assign taxonomy
taxa.dada2 <- assignTaxonomy(seqtab.nochim, path_train_set , multithread=8)
taxa.dada2 <- addSpecies(taxa.dada2, path_train_set_species)

# determine metadata for amplicone reads
mdat <- data.frame(Filename = c("S26", "S27", "S28", "S29", "S30", "S6",  "S7",  "S8",  "S9", "S10"),
                   Sample = c(rep("N2", 5), rep("N1", 5)),
                   Source = rep("Amplicone", 10))
rownames(mdat) <- mdat$Filename

# sort metadata as seqtable and rename seqtable according new names
mdat <- mdat %>% arrange(factor(Filename, levels = rownames(seqtab.nochim)))
if (all(mdat$Filename == rownames(seqtab.nochim))) {
  rownames(seqtab.nochim) <- rownames(mdat)
  print("Ok")
}

# make a phyloseq object
ps.amp <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(mdat), 
               tax_table(taxa.dada2))

# extract dna seqs in taxa_names end export them
dna <- Biostrings::DNAStringSet(taxa_names(ps.amp))
names(dna) <- taxa_names(ps.amp)
ps.amp <- merge_phyloseq(ps.amp, dna)
taxa_names(ps.amp) <- paste0("ASV", seq(ntaxa(ps.amp)))


#Read 16s from metagenomes, add metadata
ps.mg <- import_biom(BIOMfilename = 'qiime/otus/otu_table.biom', taxa_are_rows = TRUE)
colnames(ps.mg@tax_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

map <- data.frame(Filename = c("N1", "N2"),
                  Sample = c("N1", "N2"),
                  Source = rep("Metagenome", 2))
rownames(map) <- map$Filename

# create a phyloseq object
ps.mg <- phyloseq(otu_table(ps.mg),
                  sample_data(map),
                  tax_table(ps.mg))

# Bargraph function for amplicons
bargraph <- function(ps, rank, threshold){
  #  physeq2 = filter_taxa(ps, function(x) mean(x) > 0.1, TRUE)
  
  ps2 <- tax_glom(ps, taxrank = rank)
  ps3 = transform_sample_counts(ps2, function(x) x / sum(x) )
  data <- psmelt(ps3) # create dataframe from phyloseq object
  data$Plot <- as.character(data[,rank]) #convert to character
  data$Plot[data$Abundance < threshold] <- paste0("<", threshold, " abund.")
  medians <- ddply(data, ~Plot, function(x) c(median=median(x$Abundance)))
  remainder <- medians[medians$median <= threshold,]$Plot
  p <- ggplot(data=data, aes(x=Sample, y=Abundance, fill=Plot))
  p + geom_bar(aes(), stat="identity", position="stack") + theme_light() +
    scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
                                 "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
    theme(legend.position="bottom") + guides() +
    theme(axis.text.x = element_text(angle = 90))
}

# Draw bargraphs
bargraph(ps.amp, "Phylum", 0.05) + facet_grid(. ~ sample_Sample, scales = 'free_x')
plot_bar(ps.mg, fill = "Phylum")

