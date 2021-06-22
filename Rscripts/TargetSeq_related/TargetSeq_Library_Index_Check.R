#target seq library index check

library(tidyverse)
library(readxl)

index_dt <- read_excel('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Targetseq_Library_info_SYPedit.xlsx', sheet ='index_check')
index_dt <- index_dt %>% filter(memo == 'included')
nrow(index_dt) #319
index_dt %>% group_by(index_sequence) %>% count() %>% group_by(n) %>% count()

index_dt %>% filter(weight == 'high') %>% pull(index_sequence) %>% unique() %>% length()


319/8

f_dt <- index_dt %>% filter(macrogen_id %in% c('2004KNM-0042','2006KNM-0026'))
f_dt %>% group_by(index_sequence) %>% count() %>% group_by(n) %>% count()


ct_dt <- index_dt %>% group_by(index_sequence) %>% count() %>% arrange(desc(n))
cluster_dt <- tibble(sample = as.character(), cluster_num = as.integer(), index_sequence = as.character())
n=0
for (seq in ct_dt$index_sequence){
  print(seq)
  sample_list <- index_dt %>% filter(index_sequence == seq) %>% pull(sample)
  for (t_sam in sample_list){
    n= n+1
    m=n%%8
    cluster_dt[n,] <- c(t_sam, m, seq)
  }
}
cluster_dt$cluster_num[cluster_dt$cluster_num == 0] <- '8'
cluster_dt %>% group_by(cluster_num) %>% count()
cluster_dt %>% filter(cluster_num == 1)  %>% group_by(index_sequence) %>% count() %>% group_by(n) %>% count()
cluster_dt %>% filter(cluster_num == 2)  %>% group_by(index_sequence) %>% count() %>% group_by(n) %>% count()
cluster_dt %>% filter(cluster_num == 3)  %>% group_by(index_sequence) %>% count() %>% group_by(n) %>% count()
cluster_dt %>% filter(cluster_num == 4)  %>% group_by(index_sequence) %>% count() %>% group_by(n) %>% count()
cluster_dt %>% filter(cluster_num == 5)  %>% group_by(index_sequence) %>% count() %>% group_by(n) %>% count()
cluster_dt %>% filter(cluster_num == 6)  %>% group_by(index_sequence) %>% count() %>% group_by(n) %>% count()
cluster_dt %>% filter(cluster_num == 7)  %>% group_by(index_sequence) %>% count() %>% group_by(n) %>% count()
cluster_dt %>% filter(cluster_num == 8)  %>% group_by(index_sequence) %>% count() %>% group_by(n) %>% count()

index_dt <- left_join(index_dt,cluster_dt)
index_dt %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Targetseq_index_clustering.tsv')
