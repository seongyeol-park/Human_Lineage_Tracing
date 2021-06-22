#change lineage id L1 <- L2 in DB3
stop('abc')

library(tidyverse)

#summary per sample
meta_dt1 <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_200908.txt')
orig_order <- colnames(meta_dt1)
meta_dt1 <- meta_dt1 %>% mutate(new_lineage_id_snv = ifelse(deadbody == 'DB8' & substr(lineage_id_snv, 1,2) == 'L1', gsub('L1','L2', lineage_id_snv), 
                                                ifelse(deadbody == 'DB8' & substr(lineage_id_snv,1,2) == 'L2', gsub('L2','L1', lineage_id_snv), lineage_id_snv))) %>%
  mutate(new_lineage_id_cindel = ifelse(deadbody == 'DB8' & substr(lineage_id_cindel, 1,2) == 'L1', gsub('L1','L2', lineage_id_cindel), 
                                     ifelse(deadbody == 'DB8' & substr(lineage_id_cindel,1,2) == 'L2', gsub('L2','L1', lineage_id_cindel), lineage_id_cindel))) 

m_meta_dt1 <- meta_dt1 %>% select(-lineage_id_snv, -lineage_id_cindel) %>% dplyr::rename(lineage_id_snv_only = new_lineage_id_snv, lineage_id = new_lineage_id_cindel)
colnames(m_meta_dt1)
orig_order <- gsub('lineage_id_cindel','lineage_id',gsub('lineage_id_snv', 'lineage_id_snv_only',orig_order))
m_meta_dt1 <- m_meta_dt1[,orig_order]
m_meta_dt1 %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_210227.txt')


#summary per lineage
meta_dt2 <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_200908.txt')
orig_order <- colnames(meta_dt2)
m_meta_dt2 <- meta_dt2 %>% mutate(new_lineage_id = ifelse(deadbody == 'DB8' & substr(lineage_id, 1,2) == 'L1', gsub('L1','L2', lineage_id), 
                                                ifelse(deadbody == 'DB8' & substr(lineage_id,1,2) == 'L2', gsub('L2','L1', lineage_id), lineage_id))) %>%
  mutate(new_lineage_id2 = ifelse(deadbody == 'DB8' & substr(lineage_id2, 5,6) == 'L1', gsub('L1','L2', lineage_id2), 
                                        ifelse(deadbody == 'DB8' & substr(lineage_id2,5,6) == 'L2', gsub('L2','L1', lineage_id2), lineage_id2)))
m_meta_dt2 <- m_meta_dt2 %>% select(-lineage_id, -lineage_id2) %>% dplyr::rename(lineage_id = new_lineage_id, lineage_id2 =new_lineage_id2)
colnames(m_meta_dt2)
m_meta_dt2 <- m_meta_dt2[,orig_order]
m_meta_dt2 %>% arrange(deadbody, lineage_id) %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_210227.txt')


#