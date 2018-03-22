# source("http://bioconductor.org/biocLite.R")
# biocLite("GenomicRanges")
# install.packages(tidyverse)
# install_github('arendsee/rhmmer')
# install.packages("micropan")
library(micropan)
library(GenomicRanges)
library(tidyverse)
library(rhmmer)

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
  stop("all arguments must be supplied.n", call.=FALSE)
}

s_name=opt$species_name
g_name=opt$genome_name

setwd(paste("~/fastdir/InsectCarp/",s_name,sep=""))

# Read in hmmer data
hmmscan_out = read_domtblout(file = paste(s_name , "_ORFs.fasta.out", sep="")) %>%
  filter(domain_ievalue < 0.5) %>%
  mutate(hmm_lem = hmm_len = hmm_to - hmm_from + 1) %>%
  filter(hmm_len >= domain_len * 0.9)

# Read in types of each query
repeat_types = read_tsv("~/fastdir/krishna_databases/repeat_types.tsv", col_names = c("repeat_name", "repeat_type"))

# Make a list of ranges with hits
hmmscan_list = hmmscan_out %>%
  select(query_name, ali_from, ali_to)
hmmscan_list = GenomicRanges::reduce(with(hmmscan_list, GRanges(seqnames = query_name, IRanges(start=ali_from, end=ali_to)))) %>%
  as_tibble() %>%
  select(1,2,3)

# Make a empty tibble for filtered data
hmmscan_filter_1 = hmmscan_out[0,]
for(i in (1:nrow(hmmscan_list))){
  # Select hits in range, sort by hit length, keep top hit
  hmmscan_filter_1 =  hmmscan_out %>%
    filter(query_name == hmmscan_list$seqnames[i], ali_from >= hmmscan_list$start[i], ali_to <= hmmscan_list$end[i]) %>%
    arrange(!(ali_to - ali_from)) %>%
    dplyr::slice(1) %>%
    bind_rows(hmmscan_filter_1)
}


# Sort filtered list by query name and alignment position, add type of repeat
hmmscan_filter = hmmscan_filter_1 %>%
  arrange(query_name, ali_from, ali_to) %>%
  mutate(repeat_name = gsub("\\|.*", "", query_name))%>%
  inner_join(repeat_types, by = "repeat_name")

hmmscan_filter = hmmscan_filter %>%
  mutate(strand = gsub(".*\\|", "", query_name), start  = gsub(".*\\:", "", strand), end = gsub(".*-", "", start), repeat_length = gsub(".*\\(", "", end)) %>%
  mutate(strand = gsub(".\\:.*", "", strand), start  = gsub("-.*", "", start), end = gsub("\\(.*", "", end), repeat_length = gsub(")", "", repeat_length)) %>%
  mutate(start = as.integer(start), end = as.integer(end), hmm_len = hmm_to - hmm_from + 1) %>%
  arrange(repeat_type, repeat_name, start, ali_from)


CR1s = hmmscan_filter %>%
  filter(repeat_type == "CR1")

CR1_counts = table(CR1s$query_name) %>%
  as_tibble() %>%
  mutate(repeat_name = gsub("\\|.*", "",Var1))
table(CR1_counts$repeat_name) %>%
  as_tibble() %>%
  arrange(-n)
  
intact_CR1s = CR1s %>%
  



# start = start of ORF in repeat; ali_from = start protein in ORF, start of ORF in repeat = start + ali_from - 1, 
# end   = end of ORF in repeat; ali_to = end of protein in ORF, end of ORF in repeat = start + ali_to - 1

             # # Make list of ORFs containing RT
# RT_hits = hmmscan_out %>%
#   filter(domain_name == "RVT_1" | domain_name == "RVT_3") %>%
#   select(4) %>%
#   distinct()
# 
# # Extract ORFs containing RT
# hmmscan_contains_RT = hmmscan_out %>%
#   filter(query_name %in% RT_hits$query_name)
# 
# # Make list of ORFs containing RT and EN
# RT_EN_hits =hmmscan_contains_RT %>%
#   filter(domain_name == "Exo_endo_phos_2" | domain_name == "Exo_endo_phos") %>%
#   select(4) %>%
#   distinct()
# 
# # Extract ORFs containing RT and EN
# hmmscan_contains_RT_EN = hmmscan_out %>%
#   filter(query_name %in% RT_EN_hits$query_name)
# 
# # Extract ORFs containing RT and but not EN
# hmmscan_contains_RT_noEN = hmmscan_out %>%
#   filter(!query_name %in% hmmscan_contains_RT_EN$query_name)
# 
# # Extract ORFs not containing RT
# hmmscan_contains_noRT = hmmscan_out %>%
#   filter(!query_name %in% RT_hits$query_name)
# 
# hmmscan_out_filter = hmmscan_out %>%
#   arrange(query_name, ali_from, ali_to) %>%
#   filter((query_name == lag(query_name) & ali_from <= lag(ali_to)) | (query_name == lead(query_name) & ali_to >= lead(ali_from)))
# 
#   # mutate(query_name_hit = paste(query_name, "_", domain_name, sep=""))
# 
# 
# 
# 
# 
# 
# list = hmmscan_out_filter %>%
#   select(hmm_from, hmm_to, ali_from, ali_to, domain_ievalue, query_name, domain_name)