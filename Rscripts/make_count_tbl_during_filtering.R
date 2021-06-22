#make filtering count table  #Reviewer2 number1

library(tidyverse)
library(ggsci)
library(cowplot)

meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_210330.txt') %>% filter(current_final_for_lineage == 'Y')

theme_syp = theme_bw() + theme(axis.text = element_text(size=15), axis.title = element_text(size=18), strip.text = element_text(size=15), panel.border = element_blank(), axis.line = element_line(), title = element_text(size=18))

source('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Lineage_Tracing_Color_Palette.R')

fd_tbl = tribble(~deadbody, ~snv_fd, ~indel_fd, ~SV_fd,
                 'DB2','/home/users/sypark/00_Project/06_LineageTracing/db2/07_substitution','/home/users/sypark/00_Project/06_LineageTracing/db2/08_indel','/home/users/sypark/00_Project/06_LineageTracing/db2/09_delly/02_merged_annotated',
                 'DB3','/home/users/sypark/00_Project/06_LineageTracing/db3/04-3_substitution/re_processing_with_new_batch','/home/users/sypark/00_Project/06_LineageTracing/db3/04-4_indel/reprocessing_with_newbatch','/home/users/sypark/00_Project/06_LineageTracing/db3/10_Delly/02_somatic_call/02_merged_call',
                 'DB5','/home/users/sypark/00_Project/06_LineageTracing/db5/06_substitution','/home/users/sypark/00_Project/06_LineageTracing/db5/07_indel','/home/users/sypark/00_Project/06_LineageTracing/db5/10_delly/02_merged_annotated',
                 'DB6','/home/users/sypark/00_Project/06_LineageTracing/db6/07_substitution/filter3','/home/users/sypark/00_Project/06_LineageTracing/db6/11_indel','/home/users/sypark/00_Project/06_LineageTracing/db6/10_Delly/02_merged',
                 'DB8','/home/users/sypark/00_Project/06_LineageTracing/db8/05_substitution','/home/users/sypark/00_Project/06_LineageTracing/db8/07_indel','/home/users/sypark/00_Project/06_LineageTracing/db8/08_delly/02_merged_annotated',
                 'DB9','/home/users/sypark/00_Project/06_LineageTracing/db9/05_substitution','/home/users/sypark/00_Project/06_LineageTracing/db9/07_indel','/home/users/sypark/00_Project/06_LineageTracing/db9/08_delly/02_merged_annotated',
                 'DB10','/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution','/home/users/sypark/00_Project/06_LineageTracing/db10/07_indel','/home/users/sypark/00_Project/06_LineageTracing/db10/08_delly/02_merged_annotated'
                 )


#SNV performance
#2021-04-02 all done: DB2, DB3, DB5, DB6, DB8, DB9, DB10
for (this_db in c('DB3','DB5','DB6')){
  print(this_db)
  #snv count matrix (DB2 done)
  cmd1=paste0('wc -l ',fd_tbl$snv_fd[fd_tbl$deadbody == this_db], '/*common' )
  res1 <- system(cmd1, intern =T)
  res1 <- res1[grepl('total', res1)==F & is.na(res1) == F]
  n_row1 <- lapply(res1, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[1]) %>% as.numeric()
  id1 <- basename(paste0('/home',lapply(res1, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[2]) %>% as.character()))
  id1 <- lapply(id1, function(x) unlist(strsplit(x, '\\.'))[1]) %>% as.character()
  tmp_dt1 <- tibble(sample_id = id1, n_row1 = n_row1)
  
  cmd2=paste0('wc -l ',fd_tbl$snv_fd[fd_tbl$deadbody == this_db], '/*.cmn_fi' )
  res2 <- system(cmd2, intern =T)
  res2 <- res2[grepl('total', res2)==F & is.na(res2) == F]
  n_row2 <- lapply(res2, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[1]) %>% as.numeric()
  id2 <- basename(paste0('/home',lapply(res2, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[2]) %>% as.character()))
  id2 <- lapply(id2, function(x) unlist(strsplit(x, '\\.'))[1]) %>% as.character()
  tmp_dt2 <- tibble(sample_id = id2, n_row2 = n_row2)
  
  if (this_db %in% c('DB3','DB5','DB6')){
    cmd3=paste0('wc -l ',fd_tbl$snv_fd[fd_tbl$deadbody == this_db], '/*.anv.ri_fi.cmn_fi' )
  } else{
    cmd3=paste0('wc -l ',fd_tbl$snv_fd[fd_tbl$deadbody == this_db], '/*.anv.ri_fi' )
  }
  res3 <- system(cmd3, intern =T)
  res3 <- res3[grepl('total', res3)==F & is.na(res3) ==F]
  n_row3 <- lapply(res3, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[1]) %>% as.numeric()
  id3 <- basename(paste0('/home',lapply(res3, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[2]) %>% as.character()))
  id3 <- lapply(id3, function(x) unlist(strsplit(x, '\\.'))[1]) %>% as.character()
  tmp_dt3 <- tibble(sample_id = id3, n_row3 = n_row3)
  
  cmd4=paste0('wc -l ',fd_tbl$snv_fd[fd_tbl$deadbody == this_db], '/*.anv.fi' )
  res4 <- system(cmd4, intern =T)
  res4 <- res4[grepl('total', res4)==F & is.na(res4) ==F]
  n_row4 <- lapply(res4, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[1]) %>% as.numeric()
  id4 <- basename(paste0('/home',lapply(res4, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[2]) %>% as.character()))
  id4 <- lapply(id4, function(x) unlist(strsplit(x, '\\.'))[1]) %>% as.character()
  tmp_dt4 <- tibble(sample_id = id4, n_row4 = n_row4)
  
  res5 <- read_tsv(list.files(fd_tbl$snv_fd[fd_tbl$deadbody == this_db], pattern = 'g_fi$', full.names=T), col_types = cols(`#CHROM`='c'))
  t_id_list <- meta_dt %>% filter(deadbody == this_db) %>% pull(sample_id)
  tmp_dt5 <- tibble(sample_id = as.character(), n_row5 = as.numeric())
  n=0
  for (t_id in t_id_list){
    n=n+1
    print(t_id)
    this_n <- res5 %>% filter(grepl(t_id, case_list)) %>% nrow()
    tmp_dt5[n,] <- list(t_id, this_n)
  }

  snv_count_tbl <- left_join(left_join(left_join(left_join(tmp_dt1, tmp_dt2), tmp_dt3), tmp_dt4), tmp_dt5) %>% filter(is.na(n_row5) == F)
  snv_count_tbl %>% write_tsv(paste0("/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/count_tbl_during_filtering/", this_db, '_snv_count_tbl.tsv'))
}

#indel performance
#2021-04-02 all done:DB2, DB3, DB5, DB6, DB8, DB9, DB10
for (this_db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
  print(this_db)
  #indel count matrix
  cmd1=paste0('wc -l ',fd_tbl$indel_fd[fd_tbl$deadbody == this_db], '/*common' )
  res1 <- system(cmd1, intern =T)
  res1 <- res1[grepl('total', res1)==F]
  n_row1 <- lapply(res1, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[1]) %>% as.numeric()
  id1 <- basename(paste0('/home',lapply(res1, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[2]) %>% as.character()))
  id1 <- lapply(id1, function(x) unlist(strsplit(x, '\\.'))[1]) %>% as.character()
  tmp_dt1 <- tibble(sample_id = id1, n_row1 = n_row1)
  
  cmd1=paste0('wc -l ',fd_tbl$indel_fd[fd_tbl$deadbody == this_db], '/*cmn_fi' )
  res1 <- system(cmd1, intern =T)
  res1 <- res1[grepl('total', res1)==F]
  n_row1 <- lapply(res1, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[1]) %>% as.numeric()
  id1 <- basename(paste0('/home',lapply(res1, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[2]) %>% as.character()))
  id1 <- lapply(id1, function(x) unlist(strsplit(x, '\\.'))[1]) %>% as.character()
  tmp_dt12 <- tibble(sample_id = id1, n_row12 = n_row1)
  
  if (this_db %in% c('DB6','DB8','DB9','DB10')){
    cmd2=paste0('wc -l ',fd_tbl$indel_fd[fd_tbl$deadbody == this_db], '/*.rasm2.ri_fi' )
    res2 <- system(cmd2, intern =T)
    res2 <- res2[grepl('total', res2)==F & is.na(res2)==F]
    n_row2 <- lapply(res2, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[1]) %>% as.numeric()
    id2 <- basename(paste0('/home',lapply(res2, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[2]) %>% as.character()))
    id2 <- lapply(id2, function(x) unlist(strsplit(x, '\\.'))[1]) %>% as.character()
    tmp_dt2 <- tibble(sample_id = id2, n_row2 = n_row2)
    
    cmd3=paste0('wc -l ',fd_tbl$indel_fd[fd_tbl$deadbody == this_db], '/*.rasm2.rinp_fi' )
    res3 <- system(cmd3, intern =T)
    res3 <- res3[grepl('total', res3)==F & is.na(res3)==F]
    n_row3 <- lapply(res3, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[1]) %>% as.numeric()
    id3 <- basename(paste0('/home',lapply(res3, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[2]) %>% as.character()))
    id3 <- lapply(id3, function(x) unlist(strsplit(x, '\\.'))[1]) %>% as.character()
    tmp_dt3 <- tibble(sample_id = id3, n_row3 = n_row3)
    
    cmd4=paste0('wc -l ',fd_tbl$indel_fd[fd_tbl$deadbody == this_db], '/*.rasm2.fi' )
    res4 <- system(cmd4, intern =T)
    res4 <- res4[grepl('total', res4)==F & is.na(res4)==F]
    n_row4 <- lapply(res4, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[1]) %>% as.numeric()
    id4 <- basename(paste0('/home',lapply(res4, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[2]) %>% as.character()))
    id4 <- lapply(id4, function(x) unlist(strsplit(x, '\\.'))[1]) %>% as.character()
    tmp_dt4 <- tibble(sample_id = id4, n_row4 = n_row4)
    
    
  } else{
    cmd2=paste0('wc -l ',fd_tbl$indel_fd[fd_tbl$deadbody == this_db], '/*.rasm.ri_fi' )
    res2 <- system(cmd2, intern =T)
    res2 <- res2[grepl('total', res2)==F & is.na(res2)==F]
    n_row2 <- lapply(res2, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[1]) %>% as.numeric()
    id2 <- basename(paste0('/home',lapply(res2, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[2]) %>% as.character()))
    id2 <- lapply(id2, function(x) unlist(strsplit(x, '\\.'))[1]) %>% as.character()
    tmp_dt2 <- tibble(sample_id = id2, n_row2 = n_row2)
    
    cmd3=paste0('wc -l ',fd_tbl$indel_fd[fd_tbl$deadbody == this_db], '/*.rasm.rinp_fi' )
    res3 <- system(cmd3, intern =T)
    res3 <- res3[grepl('total', res3)==F & is.na(res3)==F]
    n_row3 <- lapply(res3, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[1]) %>% as.numeric()
    id3 <- basename(paste0('/home',lapply(res3, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[2]) %>% as.character()))
    id3 <- lapply(id3, function(x) unlist(strsplit(x, '\\.'))[1]) %>% as.character()
    tmp_dt3 <- tibble(sample_id = id3, n_row3 = n_row3)
    
    cmd4=paste0('wc -l ',fd_tbl$indel_fd[fd_tbl$deadbody == this_db], '/*.rasm.fi' )
    res4 <- system(cmd4, intern =T)
    res4 <- res4[grepl('total', res4)==F & is.na(res4)==F]
    n_row4 <- lapply(res4, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[1]) %>% as.numeric()
    id4 <- basename(paste0('/home',lapply(res4, function(x) unlist(strsplit(gsub(' ','', x), '/home'))[2]) %>% as.character()))
    id4 <- lapply(id4, function(x) unlist(strsplit(x, '\\.'))[1]) %>% as.character()
    tmp_dt4 <- tibble(sample_id = id4, n_row4 = n_row4)
  }
  indel_count_tbl <- left_join(left_join(left_join(left_join(tmp_dt1,tmp_dt12,), tmp_dt2), tmp_dt3), tmp_dt4)
  indel_count_tbl %>% write_tsv(paste0("/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/count_tbl_during_filtering/", this_db, '_indel_count_tbl.tsv'))
}  


#SV performance

#SV folder path write and run
#all done: DB2, DB3, DB5, DB6,DB8,DB9,DB10 
for (this_db in c('DB5')){
  print(this_db)
  #SV count matrix
  print(this_db)
  t_id_list <- meta_dt %>% filter(deadbody == this_db) %>% pull(sample_id)
  
  path_list1 = list.files(fd_tbl$SV_fd[fd_tbl$deadbody == this_db], pattern = 'delly.vcf$', full.names=T)
  tmp_dt1 <- tibble(sample_id = as.character(), n_row1=as.numeric())
  n=0
  for (t_id in t_id_list){
    n = n+1
    t_path <- path_list1[grepl(paste0(t_id,'\\.'), path_list1)]
    res1 <- read_tsv(t_path, comment = '##', col_types = cols(`#CHROM`='c'))
    n_row1 = nrow(res1)  
    tmp_dt1[n,] <- list(t_id, n_row1)
  }
  
  path_list2 = list.files(fd_tbl$SV_fd[fd_tbl$deadbody == this_db], pattern = 'pre_rionly_fi$', full.names=T)
  tmp_dt2 <- tibble(sample_id = as.character(), n_row2=as.numeric())
  n=0
  for (t_id in t_id_list){
    n = n+1
    t_path <- path_list2[grepl(paste0(t_id,'\\.'), path_list2)]
    res2 <- read_tsv(t_path, comment = '##', col_types = cols(`#CHROM`='c'))
    n_row2 = nrow(res2)  
    tmp_dt2[n,] <- list(t_id, n_row2)
  }
  
  path_list3 = list.files(fd_tbl$SV_fd[fd_tbl$deadbody == this_db], pattern = 'pre_fi$', full.names=T)
  tmp_dt3 <- tibble(sample_id = as.character(), n_row3=as.numeric())
  n=0
  for (t_id in t_id_list){
    n = n+1
    t_path <- path_list3[grepl(paste0(t_id,'\\.'), path_list3)]
    res3 <- read_tsv(t_path, comment = '##', col_types = cols(`#CHROM`='c'))
    n_row3 = nrow(res3)  
    tmp_dt3[n,] <- list(t_id, n_row3)
  }
  if(this_db %in% c('DB2','DB5')){
    path_list4 = list.files(fd_tbl$SV_fd[fd_tbl$deadbody == this_db], pattern = 'fi_nOL$', full.names=T)
    tmp_dt4 <- tibble(sample_id = as.character(), n_row4=as.numeric())
    n=0
    for (t_id in t_id_list){
      n = n+1
      t_path <- path_list4[grepl(paste0(t_id,'\\.'), path_list4)]
      res4 <- read_tsv(t_path, comment = '##', col_types = cols(`#CHROM`='c'))
      n_row4 = nrow(res4)  
      tmp_dt4[n,] <- list(t_id, n_row4)
    }
  } else {
    path_list4 = list.files(fd_tbl$SV_fd[fd_tbl$deadbody == this_db], pattern = 'SVvaf.ri_fi$', full.names=T)
    tmp_dt4 <- tibble(sample_id = as.character(), n_row4=as.numeric())
    n=0
    for (t_id in t_id_list){
      n = n+1
      t_path <- path_list4[grepl(paste0(t_id,'\\.'), path_list4)]
      res4 <- read_tsv(t_path, comment = '##', col_types = cols(`#CHROM`='c'))
      n_row4 = nrow(res4)  
      tmp_dt4[n,] <- list(t_id, n_row4)
    }
  }
 
  if (this_db %in% c('DB2','DB5')){
    path_list5 = list.files(fd_tbl$SV_fd[fd_tbl$deadbody == this_db], pattern = 'g_fi$', full.names=T)
    tmp_dt5 <- tibble(sample_id = as.character(), n_row5=as.numeric())
    n=0
    for (t_id in t_id_list){
      n = n+1
      t_path <- path_list5[grepl(paste0(t_id,'\\.'), path_list5)]
      res5 <- read_tsv(t_path, comment = '##', col_types = cols(`#CHROM`='c'))
      n_row5 = nrow(res5)  
      tmp_dt5[n,] <- list(t_id, n_row5)
    }
  } else{
    path_list5 = list.files(fd_tbl$SV_fd[fd_tbl$deadbody == this_db], pattern = 'SVvaf.fi_OL$', full.names=T)
    tmp_dt5 <- tibble(sample_id = as.character(), n_row5=as.numeric())
    n=0
    for (t_id in t_id_list){
      n = n+1
      t_path <- path_list5[grepl(paste0(t_id,'\\.'), path_list5)]
      res5 <- read_tsv(t_path, comment = '##', col_types = cols(`#CHROM`='c'))
      n_row5 = nrow(res5)  
      tmp_dt5[n,] <- list(t_id, n_row5)
    }
  }
  
  
  path_list6 = list.files(fd_tbl$SV_fd[fd_tbl$deadbody == this_db], pattern = 'igv_check$', full.names=T)
  tmp_dt6 <- tibble(sample_id = as.character(), n_row6=as.numeric())
  n=0
  for (t_id in t_id_list){
    n = n+1
    t_path <- path_list6[grepl(paste0(t_id,'\\.'), path_list6)]
    res6 <- read_tsv(t_path, comment = '##', col_types = cols(`#CHROM`='c'))
    n_row6 = nrow(res6)  
    tmp_dt6[n,] <- list(t_id, n_row6)
  }
  
  SV_count_tbl <- left_join(left_join(left_join(left_join(left_join(tmp_dt1, tmp_dt2), tmp_dt3), tmp_dt4), tmp_dt5), tmp_dt6)
  SV_count_tbl %>% write_tsv(paste0("/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/count_tbl_during_filtering/", this_db, '_SV_count_tbl.tsv'))
}
  


#merge results
fd_path='/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/count_tbl_during_filtering/'
snv_list <- list.files(fd_path, pattern = 'snv', full.names = T)
snv_dt <- do.call(rbind, lapply(snv_list, function(x) read_tsv(x)))
snv_dt
colnames(snv_dt) <- c('sample_id', 'snv_initial','snv_step1','snv_step2','snv_step3','snv_step4')

indel_list <- list.files(fd_path, pattern='indel', full.names=T)
indel_dt <- do.call(rbind, lapply(indel_list, function(x) read_tsv(x)))
indel_dt
colnames(indel_dt) <- c('sample_id', 'indel_initial','indel_step1','indel_step2','indel_step3','indel_step4')
indel_dt %>% group_by(sample_id) %>% count() %>% arrange(desc(n))

sv_list <- list.files(fd_path, pattern='SV', full.names=T)
sv_dt <- do.call(rbind, lapply(sv_list, function(x) read_tsv(x)))
sv_dt
colnames(sv_dt) <- c('sample_id', 'SV_initial','SV_step1','SV_step2','SV_step3','SV_step4', 'SV_step5')
sv_dt %>% group_by(sample_id) %>% count() %>% arrange(desc(n))

merged_ct_tbl <- left_join(left_join(left_join(meta_dt %>% select(deadbody, sample_id, n_snv, n_indel), snv_dt), indel_dt), sv_dt)
merged_ct_tbl <- merged_ct_tbl %>% rename(snv_step5 = n_snv, indel_step5 = n_indel)
merged_ct_tbl %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/count_tbl_during_filtering/final_merged_tbl.tsv')

colSums(is.na(merged_ct_tbl))


#plot performance

snv_tbl <- merged_ct_tbl %>% select(deadbody, sample_id, starts_with('snv'))
snv_tbl <- snv_tbl %>% gather(-deadbody, -sample_id, key="filtering_steps", value="count")
snv_tbl$filtering_steps %>% unique()
snv_tbl$filtering_steps <- factor(snv_tbl$filtering_steps, levels = c('snv_initial','snv_step1','snv_step2','snv_step3','snv_step4','snv_step5'))
snv_tbl$deadbody <- factor(snv_tbl$deadbody, levels = c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))
snv_tbl$count %>% summary()

ggplot(snv_tbl, aes(x=filtering_steps, y=log10(count)))+
  geom_boxplot(aes(fill = deadbody))+
  scale_y_continuous(limits=log10(c(1,5000000)), breaks = log10(c(1, 10, 100, 10^3, 10^4, 10^5, 10^6, 4*10^6, 10^7)), labels = c(1, 10, 100, 10^3, 10^4, 10^5, 10^6, 4*10^6,10^7))+
  ylab("No. of variants")+
  scale_x_discrete(labels = c("Initial","Step1","Step2","Step3","Step4","Final"))+
  scale_fill_manual(values =db_pal, name = '')+
  theme_syp+
  theme(axis.title.x = element_blank())
  #ggtitle("Substitutions")


indel_tbl <- merged_ct_tbl %>% select(deadbody, sample_id, starts_with('indel'))
indel_tbl <- indel_tbl %>% gather(-deadbody, -sample_id, key="filtering_steps", value="count")
indel_tbl$filtering_steps %>% unique()
indel_tbl$filtering_steps <- factor(indel_tbl$filtering_steps, levels = c('indel_initial','indel_step1','indel_step2','indel_step3','indel_step4','indel_step5'))
indel_tbl$deadbody <- factor(indel_tbl$deadbody, levels = c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))
indel_tbl$count %>% summary()

ggplot(indel_tbl, aes(x=filtering_steps, y=log10(count)))+
  geom_boxplot(aes(fill = deadbody))+
  scale_y_continuous(limits=log10(c(1,1000000)), breaks = log10(c(1, 10, 100, 10^3, 10^4, 10^5, 7*10^5)), labels = c(1, 10, 100, 10^3, 10^4, 10^5, 7*10^5))+
  ylab("No. of variants")+
  scale_x_discrete(labels = c("Initial","Step1","Step2","Step3","Step4","Final"))+
  scale_fill_manual(values =db_pal, name='')+
  theme_syp+
  theme(axis.title.x = element_blank())
  #ggtitle("Indels")

sv_tbl <- merged_ct_tbl %>% select(deadbody, sample_id, starts_with('SV'))
sv_tbl <- sv_tbl %>% gather(-deadbody, -sample_id, key="filtering_steps", value="count")
sv_tbl$filtering_steps %>% unique()
sv_tbl$filtering_steps <- factor(sv_tbl$filtering_steps, levels = c('SV_initial','SV_step1','SV_step2','SV_step3','SV_step4','SV_step5'))
sv_tbl$deadbody <- factor(sv_tbl$deadbody, levels = c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))
sv_tbl$count %>% summary()

ggplot(sv_tbl, aes(x=filtering_steps, y=log10(count)))+
  geom_boxplot(aes(fill = deadbody))+
  scale_y_continuous(limits=log10(c(1,5000000)), breaks = log10(c(1, 10, 100, 10^3, 10^4, 10^5, 1*10^6)), labels = c(1, 10, 100, 10^3, 10^4, 10^5, 10^6))+
  ylab("No. of variants")+
  scale_x_discrete(labels = c("Initial","Step1","Step2","Step3","Step4","Final"))+
  scale_fill_manual(values =db_pal, name='')+
  theme_syp+
  theme(axis.title.x = element_blank())
#ggtitle("Indels")
