library(tidyverse)


#calculate bait coverage

file_list <- list.files('/home/users/team_projects/Lineage_tracing/Targeted_seq/', recursive = T, full.names = T, pattern = '.bait_cov$') 

file_list_sorted <- file_list[grepl('sorted',file_list)]
file_list_1st <- file_list_sorted[!grepl('2nd', file_list_sorted)]
res <- lapply(file_list_1st, function(x) 
  {dt <- read_tsv(x, col_names = c('chrom','start_pos','end_pos','index','cov'), col_types = cols('chrom'='c'))
  sampleid <- str_replace(str_replace(x,'.*/',''),'\\.bait_cov','')
  value <- round(sum(dt$cov)/nrow(dt),2)
  list(sampleid, value)}
  )
res_tbl1 <- tibble(sampleid = as.character(lapply(res, function(x) x[[1]])), bait_cov1 =as.numeric(lapply(res, function(x) x[[2]])))
nrow(res_tbl1)

file_list_2nd <- file_list[grepl('2nd',file_list)]
res <- lapply(file_list_2nd, function(x) {dt <- read_tsv(x, col_names = c('chrom','start_pos','end_pos','index','cov'), col_types = cols('chrom'='c'))
  sampleid <- str_replace(str_replace(x,'.*/',''),'\\.s\\.bait_cov','')
  sampleid <- gsub('\\.bait_cov','',sampleid)
  value <- round(sum(dt$cov)/nrow(dt),2)
  list(sampleid, value)})
res_tbl2 <- tibble(sampleid = as.character(lapply(res, function(x) x[[1]])), bait_cov2 =as.numeric(lapply(res, function(x) x[[2]])))
nrow(res_tbl2)

file_list_merge <- file_list[grepl('merged',file_list)]
res <- lapply(file_list_merge, function(x) 
  {dt <- read_tsv(x, col_names = c('chrom','start_pos','end_pos','index','cov'), col_types = cols('chrom'='c'))
  sampleid <- str_replace(str_replace(x,'.*/',''),'\\.merged\\.bait_cov','')
  value <- round(sum(dt$cov)/nrow(dt),2)
  list(sampleid, value)})
res_tbl3 <- tibble(sampleid = as.character(lapply(res, function(x) x[[1]])), bait_cov_merged =as.numeric(lapply(res, function(x) x[[2]])))
nrow(res_tbl3)
res_tbl3

bait_dt <- left_join(left_join(res_tbl1, res_tbl2), res_tbl3)
nrow(bait_dt)
bait_dt %>% filter(is.na(bait_cov2))
write_tsv(bait_dt, '/home/users/team_projects/Lineage_tracing/Targeted_seq/bait_coverage.txt')

tmp_dt <- bait_dt %>% gather(-sampleid, key=batch, value=bait_coverage)
tmp_dt$batch <- factor(tmp_dt$batch, levels = c('bait_cov1','bait_cov2','bait_cov_merged'))
ggplot(tmp_dt, aes(x=batch, y=log10(bait_coverage)))+
  geom_point()+
  geom_line(aes(group=sampleid), alpha=0.3)+
  scale_y_continuous(limits=c(0,6))
