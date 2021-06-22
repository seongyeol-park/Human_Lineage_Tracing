#probability of mutation overlap by chance  #reviewer2 number2

#args1: out_fd_path
#args2: no of thread


#load library
library(BSgenome)
#find available genomes
#available.genomes()
#I'm interested in hg19
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
genome = BSgenome.Hsapiens.UCSC.hg19
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
out_fd <- args[1]
system(paste0('mkdir -p ', out_fd))
n_thread <- as.numeric(args[2])


make_converted_context <- function(str){
  change_v=c('A','T','G','C')
  names(change_v)<- c('T','A','C','G')
  converted_str <- paste0(map_chr(rev(unlist(strsplit(str,''))), function(x) change_v[x]), collapse='')
  return(converted_str)
}

make_random_segment <- function(n, readlength, genome){  #give random segment table with the number of n*size_factor
  my_chr <- c(1:22,'X','Y')
  my_chr <- gsub(pattern="^", replacement='chr', my_chr)
  size_factor <- seqlengths(genome)[my_chr]/10000000
  merged_tbl <- tibble()
  for (chrom in my_chr){
    seqlength <- seqlengths(genome)[chrom]
    start <- sample(seqlength - readlength + 1, size=round(size_factor[chrom]*n,0), replace=F )
    ranges <- IRanges(start=unlist(start), width=readlength)
    strand='+'
    gr <- GRanges(seqnames=chrom, ranges=ranges, strand=strand)
    seq <- getSeq(genome, gr)
    gr_tbl <- gr %>% as.tibble
    seq_tbl <- seq %>% as.data.frame() %>% as.tibble()
    tbl <- bind_cols(gr_tbl, seq_tbl)
    merged_tbl <- bind_rows(merged_tbl, tbl)
  }
  return(merged_tbl)
}

random_mutation_by_96context <- function(genome, context_tbl){  #context_tbl: substitution, context, count
  context32 <- c("ACA","ACC","ACG","ACT","ATA","ATC","ATG","ATT","CCA","CCC","CCG","CCT","CTA","CTC","CTG","CTT",
                 "GCA","GCC","GCG","GCT","GTA","GTC","GTG","GTT","TCA","TCC","TCG","TCT","TTA","TTC","TTG","TTT")
  n=100 # unit for sampling at one time
  merged_tbl <- tibble()
  add_data <- function(merged_tbl){
    tmp_tbl <- make_random_segment(n, 3, genome)
    merged_tbl <- bind_rows(merged_tbl, tmp_tbl)
    merged_tbl <- merged_tbl %>% mutate(context = ifelse(substr(x, 2, 2) %in% c('C','T'), x, map_chr(x, make_converted_context)))
    return(merged_tbl)
  }
  context32_tbl <- context_tbl %>% group_by(context) %>% summarise(count=sum(count))
  ct_tbl = tibble(context=context32, n=0)
  while (nrow(ct_tbl) < 32 | left_join(ct_tbl, context32_tbl) %>% mutate(gap = n-count) %>% pull(gap) %>% min() < 0){
    merged_tbl <- add_data(merged_tbl)
    ct_tbl <- merged_tbl %>% group_by(context) %>% dplyr::count() %>% filter(context %in% context32)
    ct_tbl <- ct_tbl %>% filter(context %in% context32)
    left_join(ct_tbl, context32_tbl) %>% mutate(gap = n-count) %>% pull(gap) %>% min() %>% print()
  }
  final_tbl = tibble()
  for (i in 1:nrow(context_tbl)){
    t_context = context_tbl$context[i]
    t_count = context_tbl$count[i]
    tmp_tbl <- merged_tbl %>% filter(context == t_context) %>% sample_n(t_count)
    final_tbl <- bind_rows(final_tbl, tmp_tbl)
  }
  final_tbl <- final_tbl %>% mutate(pos = start+1, ref=substr(x, 2,2)) %>% dplyr::rename(chrom = seqnames, original_context = x, context32 = context) %>% 
    select(chrom, pos, ref, original_context, context32) %>% arrange(context32)
  context32_order <- final_tbl$context32 %>% unique()
  alt_res=c()
  for (n in 1:length(context32_order)){
    t_context=context32_order[n]
    t_count = context32_tbl$count[context32_tbl$context == t_context]
    tmp_tbl <- context_tbl %>% filter(context == t_context) %>% separate(substitution, c('ref','alt'), sep='>')
    alt_res = c(alt_res, sample(tmp_tbl$alt, t_count, prob=tmp_tbl$count, replace=T))
  }
  final_tbl$context_alt <- alt_res
  conv_v <- c('A','G','C','T')
  names(conv_v) <- c('T','C','G','A')
  final_tbl <- final_tbl %>% mutate(alt = ifelse(original_context == context32, context_alt, conv_v[context_alt])) %>% 
    mutate(substitution_type = paste0(substr(context32, 2,2), '>', context_alt)) %>%
    select(chrom, pos, ref, alt, substitution_type, original_context, context32) %>% 
    arrange(chrom, pos)
  return(final_tbl)
}


spectra_list<- list.files('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link_prv/sig_SBS15187abcd_a1', 
                          pattern = 'signature_spectra.tsv$', full.names = T )

for (j in 1:1000){
  out_fd2 = paste0(out_fd, '/re_set',j)
  system(paste0('mkdir -p ',out_fd2))
  library(doMC)
  registerDoMC(n_thread)
  foreach(i = 1:length(spectra_list)) %dopar% {
    print(i)
    sample_id <- unlist(strsplit(basename(spectra_list[i]),'\\.'))[1]
    print(sample_id)
    context_tbl <- read_tsv(spectra_list[i], col_types= cols(Original='n', Reconstructed='n', Residual='n'))
    context_tbl <- context_tbl %>% dplyr::rename(substitution = Substitution, context = Trinucleotide, count=Original) %>% select(substitution, context, count)
    final_tbl <- random_mutation_by_96context(BSgenome.Hsapiens.UCSC.hg19,context_tbl)
    final_tbl %>% write_tsv(paste0(out_fd2, '/', sample_id, '.random_subs.tsv'))
  }
}

