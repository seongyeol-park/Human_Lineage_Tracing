## Draw indel heatmap and filter using all sample VAF data
library(tidyverse)
library(ComplexHeatmap)
#input paths
DB2_path <- '/home/users/sypark/00_Project/06_LineageTracing/db2/08_indel/DB2_filtered_merged.txt.cut.mvaf'
DB3_path <- '/home/users/sypark/00_Project/06_LineageTracing/db3/04-4_indel/reprocessing_with_newbatch/DB3_filtered_merged.txt.cut.mvaf'
DB5_path <- '/home/users/sypark/00_Project/06_LineageTracing/db5/07_indel/DB5_filtered_merged.txt.cut.mvaf'
DB6_path <- '/home/users/sypark/00_Project/06_LineageTracing/db6/11_indel/DB6_filtered_merged.txt.cut.mvaf'
DB8_path <- '/home/users/sypark/00_Project/06_LineageTracing/db8/07_indel/DB8_filtered_merged.txt.cut.mvaf'
DB9_path <- '/home/users/sypark/00_Project/06_LineageTracing/db9/07_indel/DB9_filtered_merged.txt.cut.mvaf'
DB10_path <- '/home/users/sypark/00_Project/06_LineageTracing/db10/07_indel/DB10_filtered_merged.txt.cut.mvaf'
#load inputs
dt <- read_tsv(DB2_path, col_types = cols(`#CHROM` = 'c'))
dt <- read_tsv(DB3_path, col_types = cols(`#CHROM` = 'c'))
dt <- read_tsv(DB5_path, col_types = cols(`#CHROM` = 'c'))
dt <- read_tsv(DB6_path, col_types = cols(`#CHROM` = 'c'))
dt <- read_tsv(DB8_path, col_types = cols(`#CHROM` = 'c'))
dt <- read_tsv(DB9_path, col_types = cols(`#CHROM` = 'c'))
dt <- read_tsv(DB10_path, col_types = cols(`#CHROM` = 'c'))
#
dt <-  dt %>% unite(var_id, c(`#CHROM`, POS, REF, ALT), sep='_') %>% select(-c(ID, QUAL, FILTER, case_list, called_sample)) 
colnames(dt)
this_ncol <- ncol(dt)
dt$pos_median <- apply(dt[2:this_ncol],1,function(x) median(x[is.na(x) == F & x>0]))
dt[is.na(dt)] <-2
dt$n_more0.2 <- rowSums(dt[2:this_ncol] > 0.2 & dt[2:this_ncol]<2)
dt$n_more0 <- rowSums(dt[2:this_ncol] > 0 & dt[2:this_ncol]<2)
dt$n_NA <- rowSums(dt[2:this_ncol] ==2)
dt$n_eva <- rowSums(dt[2:this_ncol]<2)
dt[2:this_ncol][dt[2:this_ncol]==2] <- NA
#DB2 filter and save
f_dt <- dt %>% filter(pos_median > 0.25 & n_NA < 10 & n_more0.2 < 5  & (n_eva - n_more0) > 0)
f_dt$ids_more0.2 <- apply(f_dt[,2:this_ncol], 1, function(x) paste(na.omit(names(x[x >0.2])), collapse=';'))
f_dt %>% select(var_id, ids_more0.2) %>% separate(var_id, c('#CHROM', 'POS', 'REF', 'ALT'), sep = '_') %>% 
  mutate(ID = '.') %>% select(`#CHROM`, POS, ID, REF, ALT, ids_more0.2) %>% write_tsv(paste0(DB2_path, '.final'))
#DB3 filter and save
DB3_false_list <- c("16_26714335_CA_C","17_17963714_CA_C","X_2200622_A_AGAGAGAGAAGAAAGAGGG")
f_dt <- dt %>% filter(pos_median > 0.25 & n_NA < 7 & n_more0.2 > 0  & (n_eva - n_more0) > 5 & !var_id %in% DB3_false_list)
f_dt$ids_more0.2 <- apply(f_dt[,2:this_ncol], 1, function(x) paste(na.omit(names(x[x >0.2])), collapse=';'))
f_dt %>% select(var_id, ids_more0.2) %>% separate(var_id, c('#CHROM', 'POS', 'REF', 'ALT'), sep = '_') %>% 
  mutate(ID = '.') %>% select(`#CHROM`, POS, ID, REF, ALT, ids_more0.2) %>% write_tsv(paste0(DB3_path, '.final'))
#DB5 filter and save
f_dt <- dt %>% filter(pos_median > 0.25 & n_NA < 7 & n_more0.2 < 10  & (n_eva - n_more0) > 1)
f_dt$ids_more0.2 <- apply(f_dt[,2:this_ncol], 1, function(x) paste(na.omit(names(x[x >0.2])), collapse=';'))
f_dt %>% select(var_id, ids_more0.2) %>% separate(var_id, c('#CHROM', 'POS', 'REF', 'ALT'), sep = '_') %>% 
  mutate(ID = '.') %>% select(`#CHROM`, POS, ID, REF, ALT, ids_more0.2) %>% write_tsv(paste0(DB5_path, '.final'))
#DB6 filter and save
f_dt <- dt %>% filter(pos_median > 0.25 & n_NA < 10 & n_more0.2 > 0 & (n_eva - n_more0) > 12)
f_dt$ids_more0.2 <- apply(f_dt[,2:this_ncol], 1, function(x) paste(na.omit(names(x[x >0.2])), collapse=';'))
f_dt %>% select(var_id, ids_more0.2) %>% separate(var_id, c('#CHROM', 'POS', 'REF', 'ALT'), sep = '_') %>% 
  mutate(ID = '.') %>% select(`#CHROM`, POS, ID, REF, ALT, ids_more0.2) %>% write_tsv(paste0(DB6_path, '.final'))
#DB8 filter and save
f_dt <- dt %>% filter(pos_median > 0.25 & n_NA < 10 & n_more0.2 > 0 & (n_eva - n_more0) > 8)
f_dt$ids_more0.2 <- apply(f_dt[,2:this_ncol], 1, function(x) paste(na.omit(names(x[x >0.2])), collapse=';'))
f_dt %>% select(var_id, ids_more0.2) %>% separate(var_id, c('#CHROM', 'POS', 'REF', 'ALT'), sep = '_') %>% 
  mutate(ID = '.') %>% select(`#CHROM`, POS, ID, REF, ALT, ids_more0.2) %>% write_tsv(paste0(DB8_path, '.final'))
#DB9 filter and save
f_dt <- dt %>% filter(pos_median > 0.25 & n_NA < 10 & n_more0.2 > 2 & (n_eva - n_more0) > 3)
f_dt$ids_more0.2 <- apply(f_dt[,2:this_ncol], 1, function(x) paste(na.omit(names(x[x >0.2])), collapse=';'))
f_dt %>% select(var_id, ids_more0.2) %>% separate(var_id, c('#CHROM', 'POS', 'REF', 'ALT'), sep = '_') %>% 
  mutate(ID = '.') %>% select(`#CHROM`, POS, ID, REF, ALT, ids_more0.2) %>% write_tsv(paste0(DB9_path, '.final'))
#DB10 filter and save
f_dt <- dt %>% filter(pos_median > 0.25 & n_NA < 10 & n_more0.2 > 0 & (n_eva - n_more0) > 10)
f_dt$ids_more0.2 <- apply(f_dt[,2:this_ncol], 1, function(x) paste(na.omit(names(x[x >0.2])), collapse=';'))
f_dt %>% select(var_id, ids_more0.2) %>% separate(var_id, c('#CHROM', 'POS', 'REF', 'ALT'), sep = '_') %>% 
  mutate(ID = '.') %>% select(`#CHROM`, POS, ID, REF, ALT, ids_more0.2) %>% write_tsv(paste0(DB10_path, '.final'))
###tmp
f_dt <- dt %>% filter(var_id == '20_53660398_C_CTA')
#Draw heatmap
dim(f_dt)
mx <- f_dt[1:this_ncol] %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
Heatmap(mx, column_names_gp = gpar(fontsize=12),row_names_gp = gpar(fontsize=12), cluster_columns=F)

#indel nsubs correlation
ggplot(meta_dt, aes(x=n_subs, y=n_indels, color = deadbody))+
  geom_point(alpha=0.5)+
  coord_cartesian(xlim=c(0,10000), ylim=c(0,300))+
  theme_bw(base_size=12)
