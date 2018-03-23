# source("http://bioconductor.org/biocLite.R")
# biocLite("GenomicRanges")
# install.packages(tidyverse)
# install.packages("optparse")

library(GenomicRanges)
library(tidyverse)
library(optparse)

option_list = list(
  make_option(c("-G", "--genome_name"), type="character", default=NULL,
              help="Name of genomes (including .fasta)", metavar="character"),
  make_option(c("-S", "--species_name"), type="character", default=NULL,
              help="name of species", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


if (any(is.null(opt))){
  print_help(opt_parser)
  stop("all arguments must be supplied", call.=FALSE)
}

s_name=opt$species_name
g_name=opt$genome_name


# Need check for if CARP folder and fais + igor exist
if (!file.exists(paste('~/fastdir/InsectCarp/', s_name, sep=""))){
  stop("folder containg CARP data doesn't exist", call.=FALSE)
}

setwd(paste('~/fastdir/InsectCarp/', s_name, sep=""))

if (!file.exists(paste(s_name, "_Denovo_TE_Library.fasta.fai", sep = "")) | 
    !file.exists(paste(s_name, ".igor.gff", sep = "")) | 
    !file.exists(paste('~/fastdir/InsectGenomes/', s_name, '/', g_name, ".fai", sep = ""))
    ){
  stop("gff and/or fais doesn't exist", call.=FALSE)
}

# Read in fai of consensus library, igor.gff and fai of genome
fai = read_tsv(paste(s_name, "_Denovo_TE_Library.fasta.fai", sep = ""),
  col_names = c("name", "consensus_length", 1:3)) %>%
  select(-c(3:5))  %>%
  separate(col = 1, into = c("family", "type"), sep = "[: #]") %>%
  mutate(family = gsub(pattern = "family", replacement = "", x = family)) %>%
  mutate(family = as.integer(gsub(pattern = "_consensus", replacement = "", x = family)))

genome_fai = read_tsv( paste('~/fastdir/InsectGenomes/', s_name, '/', g_name, ".fai", sep = ""),
                        col_names = c("contig", "length", 1:3)) %>%
  select(-c(3:5))

igor = read_tsv(paste(s_name, ".igor.gff", sep = ""), 
                 col_names = c("contig", 2, 3, "start", "end", 6, "strand", 8, "family")) %>%
  select(c(1,4,5,7,9)) %>%
  mutate(family = as.integer(gsub(pattern = "Family", replacement = "", x = family)))

# Select families annotated as repeats and >2500bp long
repeats = fai %>%
  filter(!grepl("Retrovirus_like|Unclassified", type), consensus_length > 2500)

if(nrow(repeats)<1){
  stop("No repeats", call.=FALSE)
}

# Select sequences from families annotated as repeats and >2500bp long
igor_repeats = igor %>%
  mutate(repeat_length = end - start + 1L) %>%
  filter(family %in% repeats$family, repeat_length > 2500)

# Collapse overlapping sequences
igor_repeats.gr = GenomicRanges::reduce(with(igor_repeats, GRanges(contig, IRanges(start=start, end=end))))
igor_repeats_bed = as_tibble(igor_repeats.gr) %>%
  select(c(1:3)) %>%
  mutate(start = start - 2000L, end = end + 2000L)

# Correct extended ends and starts if outside of contig and give each family a name
igor_repeats_bed = inner_join(igor_repeats_bed, genome_fai, by = c("seqnames" = "contig")) %>%
  mutate(label = paste("repeat", 1:base::nrow(igor_repeats_bed), sep = "_"))
if(nrow(igor_repeats_bed[igor_repeats_bed$start < 1,])> 0 ){
  igor_repeats_bed[igor_repeats_bed$start < 1,]$start = 1L
}
if(nrow(igor_repeats_bed[igor_repeats_bed$end > igor_repeats_bed$length,]) > 0){
  igor_repeats_bed[igor_repeats_bed$end > igor_repeats_bed$length,]$end = igor_repeats_bed[igor_repeats_bed$end > igor_repeats_bed$length,]$length
}

# Output families as a bed
write_tsv(x = igor_repeats_bed[,c(1:3,5)], path = paste(s_name, "_repeats.bed", sep = ""), col_names = F)