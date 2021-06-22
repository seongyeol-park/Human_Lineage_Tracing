stop('abcd')

#2021-02-27 DB8 L1 <-> L2 changed and applied, lineage_id_cindel -> lineage_id


###########################<Common processing>##############################################
#load library
library(tidyverse)
library(ggsci)
library(ggtree)
library(ggrepel)
library(cowplot)
library(scales)
library(grid)
library(gridExtra)
library(png)

#load dataset
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_210322.txt') %>% filter(current_final_for_lineage == 'Y')
raw_meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_210322.txt')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_210322.txt')
DB_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_DB_210322.txt')
rate_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Mutation_rate_simulation/df_forest_plot.txt')
m_rate_dt <- rate_dt %>% select(-ll, -ul) %>%  spread(key=type, value=mean) %>% dplyr::rename(deadbody = ID)
pmeta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_TargetSeq_210325.txt')



# path assign
nwk_tbl <- tribble(
  ~deadbody, ~nwk_path,
  "DB2", '/home/users/sypark/00_Project/06_LineageTracing/db2/10_point_mutations/DB2_Lineage_count_table.txt.nwk', 
  "DB3", '/home/users/sypark/00_Project/06_LineageTracing/db3/14_point_mutations/DB3_Lineage_count_table.txt.nwk',
  "DB5",'/home/users/sypark/00_Project/06_LineageTracing/db5/11_point_mutations/DB5_Lineage_count_table.txt.nwk',
  "DB6",'/home/users/sypark/00_Project/06_LineageTracing/db6/12_point_mutations/DB6_Lineage_count_table.txt.nwk',
  "DB8",'/home/users/sypark/00_Project/06_LineageTracing/db8/09_point_mutations/DB8_Lineage_count_table.txt.nwk',
  "DB9",'/home/users/sypark/00_Project/06_LineageTracing/db9/09_point_mutations/DB9_Lineage_count_table.txt.nwk',
  "DB10",'/home/users/sypark/00_Project/06_LineageTracing/db10/09_point_mutations/DB10_Lineage_count_table.txt.nwk'
)
pointmt_path_tbl <- tribble(
  ~deadbody, ~pointmt_path,
  "DB2","/home/users/sypark/00_Project/06_LineageTracing/db2/10_point_mutations/",
  "DB3","/home/users/sypark/00_Project/06_LineageTracing/db3/14_point_mutations/",
  "DB5","/home/users/sypark/00_Project/06_LineageTracing/db5/11_point_mutations/",
  "DB6","/home/users/sypark/00_Project/06_LineageTracing/db6/12_point_mutations/",
  "DB8","/home/users/sypark/00_Project/06_LineageTracing/db8/09_point_mutations/",
  "DB9","/home/users/sypark/00_Project/06_LineageTracing/db9/09_point_mutations/",
  "DB10","/home/users/sypark/00_Project/06_LineageTracing/db10/09_point_mutations/"
)

#edit meta_dt
meta_dt <- meta_dt %>% mutate(Source3 = ifelse(is.na(Source_side)==T, Source2, paste0(Source_side,'_',Source2)))
meta_dt$deadbody <- factor(meta_dt$deadbody, levels =c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))

#edit lng_dt
lng_dt$deadbody <- factor(lng_dt$deadbody, levels =c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))

#edit DB_dt
DB_dt$deadbody <- factor(DB_dt$deadbody, levels =c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))

#color setting

#show_col(pal_npg('nrc')(10))
db_pal <- pal_npg('nrc')(10)[c(1:5,9,7)]
names(db_pal) <- unique(meta_dt$deadbody)
cell_pal <- pal_npg('nrc')(6)
names(cell_pal) <- unique(meta_dt$Cell_type)
source2_pal <- pal_d3('category10')(6)
names(source2_pal) <- c('HN','lt_LE','lt_UE','rt_LE','trunk','rt_UE')
timing_pal <- pal_npg('nrc')(10)[c(1,3,4)]
names(timing_pal) <- c('early','early_to_late','late')
pal_combi <- c(pal_npg("nrc")(10), pal_d3('category10')(10), pal_aaas('default')(10))
sig_pal <- pal_combi[c(9,11,1,12,5,17,13)]
names(sig_pal) <- c('v3_1','v3_5','v3_7a','v3_7b','v3_7c','v3_7d','v3_18')


#Fx source
source('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/merged_analysis_functions.R')

#theme setting
theme_syp = theme_bw() + theme(axis.text = element_text(size=15), axis.title = element_text(size=18), strip.text = element_text(size=15), panel.border = element_blank(), axis.line = element_line(), title = element_text(size=18))

########Generate some tables--------------------------------------------------------

#make early point mutation list per DB
if(F){
  lng_dt <- lng_dt %>%  mutate(SBS7total_prop = SBS7a_prop + SBS7b_prop + SBS7c_prop + SBS7d_prop)
  for (this_db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
    print(this_db)
    early_lids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early') %>% .$lineage_id
    lateUV_lids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage != 'early' & SBS7total_prop >=5) %>% .$lineage_id
    lateNonUV_lids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage != 'early' & SBS7total_prop < 5) %>% .$lineage_id
    #merge early lids
    merged_tbl <- tibble()
    for (lid in early_lids){
      this_dt <- read_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',lid, '.pointmts.txt'), col_types = cols(`#CHROM`='c', REF='c',ALT='c'))
      this_dt <- this_dt %>% mutate(lineage_id2 = paste0(this_db, '_', lid))
      merged_tbl <- bind_rows(merged_tbl,this_dt)
    }
    print(nrow(merged_tbl))
    write_tsv(merged_tbl,paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',this_db, '_early_pointmts.txt') )
    merged_tbl %>% filter(nchar(REF) == 1 & nchar(ALT) == 1) %>% write_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',this_db, '_early_SNVs.txt'))
    merged_tbl %>% filter(nchar(REF) != 1 | nchar(ALT) != 1) %>% write_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',this_db, '_early_indels.txt'))
    #merge late UV lids
    merged_tbl <- tibble()
    for (lid in lateUV_lids){
      this_dt <- read_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',lid, '.pointmts.txt'), col_types = cols(`#CHROM`='c', REF='c',ALT='c'))
      this_dt <- this_dt %>% mutate(lineage_id2 = paste0(this_db, '_', lid))
      merged_tbl <- bind_rows(merged_tbl,this_dt)
    }
    print(nrow(merged_tbl))
    write_tsv(merged_tbl,paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',this_db, '_lateUV_pointmts.txt') )
    merged_tbl %>% filter(nchar(REF) == 1 & nchar(ALT) == 1) %>% write_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',this_db, '_lateUV_SNVs.txt'))
    merged_tbl %>% filter(nchar(REF) != 1 | nchar(ALT) != 1) %>% write_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',this_db, '_lateUV_indels.txt'))
    #merge late NonUV lids
    merged_tbl <- tibble()
    for (lid in lateNonUV_lids){
      this_dt <- read_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',lid, '.pointmts.txt'), col_types = cols(`#CHROM`='c', REF='c',ALT='c'))
      this_dt <- this_dt %>% mutate(lineage_id2 = paste0(this_db, '_', lid))
      merged_tbl <- bind_rows(merged_tbl,this_dt)
    }
    print(nrow(merged_tbl))
    write_tsv(merged_tbl,paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',this_db, '_lateNonUV_pointmts.txt') )
    merged_tbl %>% filter(nchar(REF) == 1 & nchar(ALT) == 1) %>% write_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',this_db, '_lateNonUV_SNVs.txt'))
    merged_tbl %>% filter(nchar(REF) != 1 | nchar(ALT) != 1) %>% write_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',this_db, '_lateNonUV_indels.txt'))
  }
}

#Make_Merged_Early_Late_Mut_List
if(F){
    ###DB total early
    DBto_merged_tbl <- tibble()
    for (this_db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
      print(this_db)
      this_dt <- read_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',this_db, '_early_pointmts.txt'), col_types = cols(`#CHROM`='c', REF='c',ALT='c'))
      DBto_merged_tbl <- bind_rows(DBto_merged_tbl, this_dt)
    }
    #save datas
    write_tsv(DBto_merged_tbl, '/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB_total_early_pointmts.txt')
    DBto_merged_tbl %>% filter(nchar(REF)==1 & nchar(ALT)==1) %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB_total_early_SNVs.txt')
    DBto_merged_tbl %>% filter(nchar(REF)!=1 | nchar(ALT)!=1) %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB_total_early_indels.txt')
    ###DB total late UV
    DBto_merged_tbl <- tibble()
    for (this_db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
      print(this_db)
      this_dt <- read_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',this_db, '_lateUV_pointmts.txt'), col_types = cols(`#CHROM`='c', REF='c',ALT='c'))
      DBto_merged_tbl <- bind_rows(DBto_merged_tbl, this_dt)
    }
    #save datas
    write_tsv(DBto_merged_tbl, '/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB_total_lateUV_pointmts.txt')
    DBto_merged_tbl %>% filter(nchar(REF)==1 & nchar(ALT)==1) %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB_total_lateUV_SNVs.txt')
    DBto_merged_tbl %>% filter(nchar(REF)!=1 | nchar(ALT)!=1) %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB_total_lateUV_indels.txt')
    ###DB total late NonUV
    DBto_merged_tbl <- tibble()
    for (this_db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
      print(this_db)
      this_dt <- read_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',this_db, '_lateNonUV_pointmts.txt'), col_types = cols(`#CHROM`='c', REF='c',ALT='c'))
      DBto_merged_tbl <- bind_rows(DBto_merged_tbl, this_dt)
    }
    #save datas
    write_tsv(DBto_merged_tbl, '/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB_total_lateNonUV_pointmts.txt')
    DBto_merged_tbl %>% filter(nchar(REF)==1 & nchar(ALT)==1) %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB_total_lateNonUV_SNVs.txt')
    DBto_merged_tbl %>% filter(nchar(REF)!=1 | nchar(ALT)!=1) %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB_total_lateNonUV_indels.txt')
}

#Make_Merged_Early_Late_Mut_List (Five deadbodies)
if(F){
  ###DB total early
  DBto_merged_tbl <- tibble()
  for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
    print(this_db)
    this_dt <- read_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',this_db, '_early_pointmts.txt'), col_types = cols(`#CHROM`='c', REF='c',ALT='c'))
    DBto_merged_tbl <- bind_rows(DBto_merged_tbl, this_dt)
  }
  #save datas
  write_tsv(DBto_merged_tbl, '/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB368910_total_early_pointmts.txt')
  DBto_merged_tbl %>% filter(nchar(REF)==1 & nchar(ALT)==1) %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB368910_total_early_SNVs.txt')
  DBto_merged_tbl %>% filter(nchar(REF)!=1 | nchar(ALT)!=1) %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB368910_total_early_indels.txt')
  ###DB total late UV
  DBto_merged_tbl <- tibble()
  for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
    print(this_db)
    this_dt <- read_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',this_db, '_lateUV_pointmts.txt'), col_types = cols(`#CHROM`='c', REF='c',ALT='c'))
    DBto_merged_tbl <- bind_rows(DBto_merged_tbl, this_dt)
  }
  #save datas
  write_tsv(DBto_merged_tbl, '/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB368910_total_lateUV_pointmts.txt')
  DBto_merged_tbl %>% filter(nchar(REF)==1 & nchar(ALT)==1) %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB368910_total_lateUV_SNVs.txt')
  DBto_merged_tbl %>% filter(nchar(REF)!=1 | nchar(ALT)!=1) %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB368910_total_lateUV_indels.txt')
  ###DB total late NonUV
  DBto_merged_tbl <- tibble()
  for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
    print(this_db)
    this_dt <- read_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/',this_db, '_lateNonUV_pointmts.txt'), col_types = cols(`#CHROM`='c', REF='c',ALT='c'))
    DBto_merged_tbl <- bind_rows(DBto_merged_tbl, this_dt)
  }
  #save datas
  write_tsv(DBto_merged_tbl, '/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB368910_total_lateNonUV_pointmts.txt')
  DBto_merged_tbl %>% filter(nchar(REF)==1 & nchar(ALT)==1) %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB368910_total_lateNonUV_SNVs.txt')
  DBto_merged_tbl %>% filter(nchar(REF)!=1 | nchar(ALT)!=1) %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB368910_total_lateNonUV_indels.txt')
}

#make link for early_to_late linege
if (F){
  el_lids <- lng_dt %>% filter(time_of_lineage == 'early_to_late') %>% pull(lineage_id2)
  dt <- tibble(lid = el_lids)
  dt <- dt %>% separate(lid, c('deadbody','lineage_id'), sep='_') %>% mutate(path = apply(.['deadbody'], 1, function(x) pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == x])) %>%
    mutate(path2 = paste0(path, 'perLineage/',lineage_id,'.pointmts.txt'), lineage_id2 = paste(deadbody, lineage_id, sep="_"))
  dt
  for (i in 1:nrow(dt)){
    lid2=dt[i,'lineage_id2']
    path2=dt[i,'path2']
    cmd = paste0('ln -s ', path2, ' /home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perLineage_EarlyToLate_link/', lid2, '.pointmts.txt')
    print(cmd)
    system(cmd)
  }
}

#make count table of point mutations per sample
if(F){
  file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/',pattern = 'vcf.anv$', full.names = T)
  dt <- do.call(rbind, lapply(file_list, function(x) {
    t_id = unlist(strsplit(unlist(strsplit(x,'//'))[2],'\\.'))[1]
    read_tsv(x, comment='##',col_types = cols(`#CHROM`='c')) %>% mutate(sample_id=t_id) %>% rename(rd_ad_vaf=t_id)
  }))
  tmp_tbl1 <- dt %>% group_by(sample_id) %>% summarise(n_pointmt = n())
  dt <- dt %>% mutate(mt_type = ifelse(nchar(REF) ==1 & nchar(ALT)==1, 'SNV',ifelse(nchar(REF) > nchar(ALT), 'del','ins')))
  tmp_tbl2 <- dt %>%group_by(sample_id, mt_type) %>% count() %>% spread(key=mt_type, value=n)
  tmp_tbl3 <- dt %>% filter(Func_refGene == 'exonic') %>% group_by(sample_id) %>% summarise(n_exonic = n())
  m_dt <- left_join(tmp_tbl1, tmp_tbl2) %>% left_join(tmp_tbl3) %>%  mutate(n_indel = del+ins) %>% rename(n_del = del, n_ins = ins, n_snv = SNV)
  write_tsv(m_dt, '/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/pointmt_number_perSample.txt' )
}

#########Method related -------------------------------------

#quality control of samples
raw_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_200908.txt')
f_dt <- raw_dt %>% filter(deadbody %in% c('DB2','DB3','DB5','DB6','DB8','DB9','DB10
                                          '))
nrow(f_dt) #394
f_dt$Memo %>% unique()
f_dt <- f_dt %>% filter(!Memo %in% c('serial_daughter','bulk','other_person', 'half_splitted_b', 'duplicated_sequencing','diverge_during_culture', 'multiclonal,poor_quality')) 
nrow(f_dt) #374

f_dt %>% arrange(seqz_cov) %>% select(sample_id, Memo, seqz_cov, peakVAF) %>% View()
f_dt %>% group_by(Memo) %>% dplyr:: count()
#sample QC plot
f_dt <- f_dt %>% mutate(new_group = ifelse(Memo %in% c('multiclonal'),'Multilineages','Single lineage'))
ggplot(f_dt, aes(x=seqz_cov, y=peakVAF, color=current_final_for_lineage, shape=new_group))+
  geom_point(alpha=0.5)+
  scale_shape_manual(values = c("Multilineages" = 4, "Single lineage" = 19))+
  scale_color_manual(values = c("Y" = "green","N"="red"))+
  scale_y_continuous(limits = c(0,65))+
  theme_syp

#merged VAF distribution of all somatic mutations
file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link',pattern = 'vcf$', full.names = T)
vcfs <- do.call(rbind, lapply(file_list, function(x) {
  r <- regexec('.*/(.*).edited.vcf',x)
  m <- regmatches(x,r)
  t_id <- map_chr(m,2)
  dt <-read_tsv(x, comment = '##', col_types = cols(`#CHROM`='c')) %>% mutate(sample_id = t_id)
  colnames(dt)[10] <- 'info'
  dt
  }))
nrow(vcfs) #1827861
vcfs <- vcfs %>% separate(info, c('RD','AD','VAF'), sep=':', convert = T)
ggplot(vcfs, aes(x=VAF))+
  geom_histogram(binwidth=2)+
  xlab("Variant-allele frequency")+ylab("No. of point mutations")+
  theme_syp



#peakVAF distribution per deadbody
ggplot(meta_dt, aes(x=peakVAF, fill=deadbody))+
  geom_histogram(binwidth=2)+
  scale_x_continuous(limits = c(20,80))+
  scale_color_manual(values = db_pal)+
  facet_wrap(~deadbody, scales = "free_y")+
  xlab("Variant-allele frequency (%)")+ylab("Count")+
  theme_syp+theme(legend.position = 'none')


#compare SBS18 proportion in DB10
dt <- bind_rows(meta_dt %>% filter(deadbody == 'DB10') %>% mutate(group = 'Single-cell clones')  %>% select(group, sample_id, SBS18),
  tribble(~group, ~sample_id, ~SBS18,
          'Culture-associated','10_ARL10-4',27.76,
          'Culture-associated','10_LA14-2',25.39,
          'Culture-associated','10_PRL3-1',23.54)
  )
ggplot(dt, aes(x=group, y=SBS18))+
  geom_boxplot(fill='gray')+
  ylab("Proportion of SBS18 (%)")+
  theme_syp+theme(axis.title.x = element_blank())


#####check VAF of point mutations in X and Y in male
if(F){
  male_ids <- meta_dt$sample_id[meta_dt$gender == 'M']
  dt <- do.call(rbind,lapply(male_ids, function(x) read_tsv(paste0('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/',x,'.edited.vcf'), comment='##', col_types = cols(`#CHROM`='c')) %>%
                               separate(get(x), c('RD','AD','VAF'), sep=':', convert = T) %>% filter(`#CHROM` %in% c('X','Y')) %>% mutate(sample_id = x, type = ifelse(nchar(REF)>1 | nchar(ALT) >1, 'indel','snv'),                                                                                                                                  deadbody = paste0('DB',unlist(strsplit(x,'_'))[1]))
  ))
  nrow(dt)
  #check samples with Xgain 
  CNV_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/merged_data/CNVs/final_segments_more10Mbp.rds')
  CNV_tbl %>% filter(chrom %in% c('X','Y') & type == 'amp' & deadbody %in% male_ids) %>% pull(sample_id) #no sample with gain
  #remove PAR
  f_dt <- dt %>% filter((`#CHROM` == 'X' & (POS < 60001 | (POS > 2699520 & POS < 154931044) | POS > 155260560)) |
                          (`#CHROM` == 'Y' & (POS < 10001 | (POS > 2649520 & POS < 59034050) | POS > 59363566)))
  nrow(f_dt)
  f_dt %>% filter(VAF < 100 & RD > 1) %>% nrow()
  f_dt %>% filter(VAF < 100 & RD > 1) %>% View()
  tmp_dt <- left_join(f_dt %>% group_by(sample_id) %>% summarise(n_Total = n()),
                      f_dt %>% filter(VAF < 70 & RD > 1) %>% group_by(sample_id) %>% summarise(n_NotHomo = n())
  ) %>% mutate(prop = n_NotHomo/n_Total) %>% arrange(desc(prop)) 
  tmp_dt[is.na(tmp_dt)] <- 0
  tmp_dt$prop %>% summary()
}


#### check VAF of point mutations in one copy region
if(F){
  CNV_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/merged_data/CNVs/final_segments_more10Mbp.rds')
  del_dt <- CNV_tbl %>% filter(totCN == 1)
  fd_path = '/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/'
  n_total_v = c()
  n_NotHomo_v = c()
  for (i in 1:nrow(del_dt)){
    t_id <- del_dt$sample_id[i]
    t_chrom <- del_dt$chrom[i]
    t_start <- del_dt$start_pos[i] + 1000000
    t_end <- del_dt$end_pos[i] - 1000000
    dt <- read_tsv(paste0(fd_path,t_id,'.edited.vcf'), comment = '##', col_types = cols(`#CHROM` = 'c')) %>% separate(get(t_id), c('RD','AD','AF'), sep=':', convert= T) %>% filter(`#CHROM` == t_chrom & POS >= t_start & POS < t_end) 
    n_total <- dt %>% nrow()
    n_NotHomo <- dt %>% filter(RD >1 & AF < 70) %>% nrow()
    n_total_v = c(n_total_v, n_total)
    n_NotHomo_v = c(n_NotHomo_v, n_NotHomo)
  }
  del_dt$n_total = n_total_v
  del_dt$n_NotHomo = n_NotHomo_v
  sum(del_dt$n_NotHomo)/sum(del_dt$n_total)
}

#Draw DB3 heatmap
if(F){
  vaf_dt <- read_tsv("~sypark/00_Project/06_LineageTracing/db3/04-3_substitution/re_processing_with_new_batch/DB3_45s_merged.txt.sampleN.fi.rmFP.Li_fi.sub_fi2.46sCall")
  sample_list <- meta_dt %>% filter(deadbody == 'DB3') %>% pull(sample_id) %>% paste0(.,'_vafpct')
  colnames(vaf_dt)
  s_vaf_dt <- vaf_dt %>% rename(CHR = `#CHROM`) %>% select(CHR, POS, REF, ALT, sample_list, '3_blood_4s_merged_vafpct')
  colnames(s_vaf_dt)
  s_vaf_mx <- s_vaf_dt %>% mutate(snv_id=paste(CHR,POS,REF,ALT, sep='_')) %>% select(-CHR, -POS, -REF, -ALT, -`3_blood_4s_merged_vafpct`) %>% as.data.frame() %>% column_to_rownames('snv_id') %>% as.matrix()
  library(ComplexHeatmap)
  rt_annot = rowAnnotation(blood = s_vaf_dt$`3_blood_4s_merged_vafpct`)
  Heatmap(s_vaf_mx, col = colorRamp2(c(0,50,100), c("blue","yellow","red")),
          right_annotation = rt_annot, row_names_gp=gpar(fontsize=5), column_names_gp=gpar(fontsize=7),
          cell_fun = function(j,i,x,y,width,height,fill){
            grid.rect(x=x, y=y, width=width, height=height,gp = gpar(col="gray",fill=NA))
          })
}

#draw heatmap with early embryonic mutations
if(F){
  library(pheatmap)
  library(grid)
  data_path=tribble(~deadbody, ~path,
                    'DB2','/home/users/sypark/00_Project/06_LineageTracing/db2/10_point_mutations/perLineage/DB2_early_pointmts.txt.27sCall',
                    'DB3','/home/users/sypark/00_Project/06_LineageTracing/db3/14_point_mutations/perLineage/DB3_early_pointmts.txt.43sCall',
                    'DB5','/home/users/sypark/00_Project/06_LineageTracing/db5/11_point_mutations/perLineage/DB5_early_pointmts.txt.28sCall',
                    'DB6','/home/users/sypark/00_Project/06_LineageTracing/db6/12_point_mutations/perLineage/DB6_early_pointmts.txt.115sCall',
                    'DB8','/home/users/sypark/00_Project/06_LineageTracing/db8/09_point_mutations/perLineage/DB8_early_pointmts.txt.47sCall',
                    'DB9','/home/users/sypark/00_Project/06_LineageTracing/db9/09_point_mutations/perLineage/DB9_early_pointmts.txt.35sCall',
                    'DB10','/home/users/sypark/00_Project/06_LineageTracing/db10/09_point_mutations/perLineage/DB10_early_pointmts.txt.39sCall'
  )
  
  make_downstream_order <- function(t_lid, this_db){
    f_lng_dt <- lng_dt %>% filter(deadbody == this_db)
    t_step <- length(unlist(strsplit(t_lid, '-')))
    f_lng_dt$time_of_lineage <- factor(f_lng_dt$time_of_lineage, levels = c('early','early_to_late','late'))
    res <- f_lng_dt %>% filter(substr(lineage_id, 1, nchar(t_lid)) == t_lid & step_n == t_step +1)%>% arrange(-n_samples, time_of_lineage, lineage_id) %>% pull(lineage_id)
    if(length(res) > 0){
      return(res)
    } else{
      return(t_lid)
    }
  }
  
  plot_list=list()
  for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
    print(this_db)
    dt <- read_tsv(data_path$path[data_path$deadbody == this_db])
    dt <- dt %>% mutate(var_id = paste0(`#CHROM`,':',POS, ' ',REF,'>',ALT)) %>% select(var_id, lineage_id2, ends_with('vafpct'), -starts_with('blood'))
    dt <- left_join(dt,lng_dt %>% select(lineage_id2, lineage_id, n_samples))

    
    t_lid_list <- lng_dt %>% filter(deadbody == this_db & step_n == 1) %>% arrange(-n_samples, lineage_id) %>% pull(lineage_id)
    while (T){
      new_lid_list = c()
      for (t_lid in t_lid_list){
        new_lid_list = c(new_lid_list, make_downstream_order(t_lid, this_db))
      }
      if (length(t_lid_list) == length(new_lid_list)){
        break
      }
      t_lid_list <- new_lid_list
    }
    l_order <- t_lid_list
    s_order <- meta_dt %>% filter(deadbody == this_db) %>% select(lineage_id, sample_id) %>% deframe()
    s_order <- s_order[l_order]
    row_order_list=c()
    for(lid in l_order){
      lid_list <- unlist(strsplit(lid, '-'))
      new_lid_list = c()
      for (i in 1: length(lid_list)){
        new_lid_list <- c(new_lid_list, paste(lid_list[1:i], collapse='-'))
      }
      added_list <- setdiff(new_lid_list, row_order_list)
      row_order_list <- c(row_order_list, added_list)
    }
    f_dt <- dt
    colnames(f_dt) <- gsub('_vafpct','',colnames(f_dt))
    row_order_list <- intersect(row_order_list, f_dt$lineage_id)
    f_dt$lineage_id <- factor(f_dt$lineage_id, levels = row_order_list)
    f_dt <- f_dt %>% arrange(lineage_id) %>% select(var_id, unname(s_order))
    mx <- f_dt %>% as.data.frame() %>%column_to_rownames('var_id') %>% as.matrix()
    n_L1 = meta_dt %>% filter(deadbody == this_db & substr(lineage_id, 1,2)== 'L1') %>% nrow()
    n_L2 = meta_dt %>% filter(deadbody == this_db & substr(lineage_id, 1,2)== 'L2') %>% nrow()
    n_L1
    n_L2
    plot_list[[this_db]] <- pheatmap(mx, cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames = F,main = this_db)
  }
  
  pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/EDFig3b_DB3_annotated.pdf')
  pheatmap(mx, cluster_cols = F, cluster_rows = F, show_rownames = T, show_colnames = T,main = this_db)
  dev.off()
  
  pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/EDFig3b.pdf')
  grid.newpage()
  pushViewport(viewport(x=0, y=0.5, width = 0.3, height = 0.5, just = c('left', 'bottom')))
  print(plot_list[['DB3']], newpage = F)
  popViewport(1)
  pushViewport(viewport(x=0.3, y=0.5, width = 0.3, height = 0.5, just = c('left', 'bottom')))
  print(plot_list[['DB8']], newpage = F)
  popViewport(1)
  pushViewport(viewport(x=0.6, y=0, width = 0.4, height = 1, just = c('left', 'bottom')))
  print(plot_list[['DB6']], newpage = F)
  popViewport(1)
  pushViewport(viewport(x=0, y=0, width = 0.3, height = 0.5, just = c('left', 'bottom')))
  print(plot_list[['DB9']], newpage = F)
  popViewport(1)
  pushViewport(viewport(x=0.3, y=0, width = 0.3, height = 0.5, just = c('left', 'bottom')))
  print(plot_list[['DB10']], newpage = F)
  popViewport(1)
  dev.off()
}

#number of samples ↑, detected early mutations ↑ 
if(F){
  
  ggplot(DB_dt, aes(x=sample_n , y=early_pointmt_n))+
    geom_point(aes(color = deadbody), size=5, alpha=0.7)+
    #geom_text_repel(aes(label = deadbody))+
    geom_smooth(method="lm")+
    scale_color_manual(values = db_pal)+
    scale_x_continuous(limits = c(0,130))+
    scale_y_continuous(limits = c(0, 200))+
    xlab("No. of samples")+ylab("No. of early point mutations")+
    theme_syp
  
  tmp_dt <- DB_dt %>% select(deadbody, sample_n, early_snv_n, early_indel_n) %>% gather(-deadbody, -sample_n, key=type, value=count)
  tmp_dt
  ggplot(tmp_dt , aes(x=sample_n, y=count, shape=type, color=deadbody))+
    geom_point(size=3)+
    geom_smooth(method = "lm", aes(group = type))+
    scale_color_manual(values = db_pal)+
    scale_x_continuous(limits = c(0,130))+
    scale_y_continuous(limits = c(0, 200))+
    xlab("No. of samples")+ylab("No. of early mutations")+
    theme_syp
  }

#plotting included samples
if(F){
  r_meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_200908.txt')
  tot_seq_dt <- r_meta_dt %>% filter(deadbody != 'DB1' & !(Memo %in% c('serial_daughter','half_splitted_b','bulk','other_person','duplicated_sequencing','diverge_during_culture','multiclonal,poor_quality'))) 
  tot_seq_dt %>% rename(Included)
  nrow(tot_seq_dt)
  colnames(tot_seq_dt)
  ggplot(tot_seq_dt, aes(x=seqz_cov, y=peakVAF, color=deadbody, shape=current_final_for_lineage))+
    geom_point(size=3, alpha=0.5)+
    scale_color_manual(values = db_pal)+
    theme_syp
}


#Window of our data (using end point of early branch)
colnames(lng_dt)
tmp <- lng_dt %>% filter(time_of_lineage == 'early_to_late') %>% group_by(deadbody) %>% summarise(median_start = median(start_mtn),
                                                                                                  mean_start = mean(start_mtn))

tmp <- left_join(tmp, rate_dt %>% filter(type == 'M0') %>% rename(deadbody = ID, M0 = mean) %>% select(deadbody, M0))
tmp <- left_join(tmp, rate_dt %>% filter(type == 'M1') %>% rename(deadbody = ID, M1 = mean) %>% select(deadbody, M1))
tmp <- tmp %>% mutate(median_stage = (median_start - M0*2)/M1+2,
               mean_stage = (mean_start - M0*2)/M1+2)
tmp <- tmp %>% mutate(median_cell_number = 2^median_stage, 
               mean_cell_number = 2^mean_stage)





#####################################generate various statistical values-------------
if(F){
  
  #number of samples
  colnames(pmeta_dt)
  pmeta_dt %>% filter(final_included == 'Y') %>% nrow()
  
  
  #number of mutations
  colnames(meta_dt)
  meta_dt %>% filter(is.na(seqz_cov)) %>% pull(sample_id)
  meta_dt$seqz_cov %>% summary()
  nrow(meta_dt)
  meta_dt %>% filter(deadbody %in% c('DB3','DB6','DB8','DB9','DB10')) %>% nrow()
  colnames(lng_dt)
  lng_dt %>% filter(is.na(n_pointmt) == F) %>% pull(n_pointmt) %>% sum()
  lng_dt %>% filter(is.na(n_pointmt) == F) %>% pull(n_subs) %>% sum()
  lng_dt %>% filter(is.na(n_pointmt) == F) %>% pull(n_indels) %>% sum()
  lng_dt %>% filter(is.na(n_pointmt) == F & time_of_lineage != 'early') %>% pull(n_pointmt) %>% sum()
  lng_dt %>% filter(is.na(n_pointmt) == F & deadbody %in% c('DB3','DB6','DB8','DB9','DB10')) %>% pull(n_pointmt) %>% sum()
  lng_dt %>% filter(is.na(n_pointmt) == F & deadbody %in% c('DB3','DB6','DB8','DB9','DB10')) %>% pull(n_subs) %>% sum()
  lng_dt %>% filter(is.na(n_pointmt) == F & deadbody %in% c('DB3','DB6','DB8','DB9','DB10')) %>% pull(n_indels) %>% sum()
  lng_dt %>% filter(is.na(n_pointmt) == F & deadbody %in% c('DB3','DB6','DB8','DB9','DB10') &
                      time_of_lineage != 'late') %>% pull(n_subs) %>% sum()
  
  lng_dt %>% filter(is.na(n_pointmt) == F & deadbody %in% c('DB3','DB6','DB8','DB9','DB10') & 
                      branch_n >0) %>% pull(n_subs) %>% sum()
  lng_dt %>% filter(is.na(n_pointmt) == F & deadbody %in% c('DB3','DB6','DB8','DB9','DB10') & 
                      branch_n >0) %>% pull(n_indels) %>% sum()
  lng_dt %>% filter(is.na(n_pointmt) == F & deadbody %in% c('DB3','DB6','DB8','DB9','DB10') & 
                      time_of_lineage == 'early') %>% pull(n_subs) %>% sum()
  lng_dt %>% filter(is.na(n_pointmt) == F & deadbody %in% c('DB3','DB6','DB8','DB9','DB10') & 
                      time_of_lineage == 'early') %>% pull(n_indels) %>% sum()

  
  #coverage
  pmeta_dt %>% pull(target_mutation_mean_depth) %>% summary()
  raw_meta_dt %>% filter(Cell_type == 'blood' & deadbody != 'DB1') %>% pull(seqz_cov) %>% summary()
  
  #signatures
  early_subs <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/sig_subs_SBS1518/DB368910_total_early_SNVs.txt.signature_exposures.tsv')
  early_subs %>% mutate(prop = Exposure*100/sum(Exposure))
  early_indels <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/sig_indel_final_pdf/DB368910_total_early_indels.txt.signature_exposures.ID12410.tsv')
  early_indels %>% mutate(prop = Exposure*100/sum(Exposure))
  nonUV_subs <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/sig_subs_SBS1518/DB368910_total_lateNonUV_SNVs.txt.signature_exposures.tsv')
  nonUV_subs  %>% mutate(prop = Exposure*100/sum(Exposure))               
  nonUV_indels <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/sig_indel_final_pdf/DB368910_total_lateNonUV_indels.txt.signature_exposures.ID589.tsv')
  nonUV_indels %>% mutate(prop = Exposure*100/sum(Exposure))  
  UV_subs <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/sig_subs_SBS15718pdf/DB368910_total_lateUV_SNVs.txt.signature_exposures.tsv')
  UV_subs  %>% mutate(prop = Exposure*100/sum(Exposure)) 
  UV_indels <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/sig_indel_final_pdf/DB368910_total_lateUV_indels.txt.signature_exposures.ID589.tsv')
  UV_indels %>%  mutate(prop = Exposure*100/sum(Exposure)) 
  
  a = early_subs$Exposure[early_subs$`Signature #` == 'SBS1']
  b = early_subs$Exposure[early_subs$`Signature #` == 'SBS5']
  c = nonUV_subs$Exposure[nonUV_subs$`Signature #` == 'SBS1']
  d = nonUV_subs$Exposure[nonUV_subs$`Signature #` == 'SBS5']
  fisher.test(matrix(c(a,b,c,d), nrow=2, byrow=T))
  
  #statistics per deadbody
  meta_dt %>% group_by(deadbody, tissue_id) %>% count() %>% select(-n) %>% group_by(deadbody) %>% count()
  meta_dt %>% group_by(deadbody) %>% count()
  meta_dt %>% group_by(deadbody, Source_class) %>% count() %>% spread(key=Source_class, value=n)
  pmeta_dt$deadbody <- factor(pmeta_dt$deadbody, levels = c('DB3','DB6','DB8','DB9','DB10'))
  pmeta_dt %>% group_by(deadbody) %>% count()
  pmeta_dt %>% group_by(deadbody, anatomy_class1) %>% count() %>% spread(key=anatomy_class1, value=n)
  pmeta_dt %>% group_by(deadbody, dominant_layer2) %>% count()
  
  
  #L1 L2 ratio
  res_tbl<- tibble(deadbody = as.character(), L1 = as.integer(), L2 = as.integer())
  n=0
  for (t_db in c('DB3','DB6','DB8','DB9','DB10')){
    n = n+1
    L1 <-lng_dt %>% filter(deadbody == t_db & time_of_lineage == 'early_to_late' & grepl(paste0(t_db,'_L1'), lineage_id2)) %>% nrow()
    L2 <- lng_dt %>% filter(deadbody == t_db & time_of_lineage == 'early_to_late' & grepl(paste0(t_db,'_L2'), lineage_id2)) %>% nrow()
    res_tbl[n,] <- list(t_db, L1, L2)
  }
  res_tbl %>% mutate(ratio = apply(.[c("L1","L2")], 1, function(x) max(x)/min(x)))
  
  #random sampling ratio
  x <-c(rep("A", 1000000), rep("B", 1000000))
  res_list=list()
  for(i in 1:1000){
    res_list[[i]] = sample(x, 112)
  }
  dt <- tibble(A = unlist(lapply(res_list, function(x) length(x[x=='A']))), B= unlist(lapply(res_list, function(x) length(x[x=='B']))))
  dt %>% mutate(ratio = apply(., 1, function(x) max(x)/min(x))) %>% filter(ratio > 6.5)
  lng_dt %>% filter(deadbody == 'DB9' & step_n == 1) %>% select(lineage_id2, mean_blood_VAFpct, median_blood_VAFpct)
  
  #DB9_L1-1-2-1-2
  dt <- read_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == 'DB9'],'perLineage/L1-1-2-1-2.pointmts.txt'), col_types = cols(`#CHROM`='c'))
  dt %>% separate(`blood_dp;blood_var;blood_vafpct`, c('blood_dp','blood_var','blood_vafpct'), sep=';', convert = T) %>% filter(blood_var >1) %>% pull(blood_vafpct) %>% mean()
  colnames(meta_dt)
  meta_dt %>% filter(Cell_type == 'skin_fb') %>% pull(n_pointmt) %>% summary()
  sim_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Mutation_rate_simulation/df_forest_plot.txt')
  sim_dt %>% filter(type == 'M1') %>% pull(mean) %>% summary()
  r_meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_200908.txt')
  r_meta_dt$Memo %>% unique()
  r_meta_dt %>% filter(deadbody != 'DB1' & !(Memo %in% c('serial_daughter','half_splitted_b'))) %>% pull(tissue_id) %>% unique() %>% length()
  r_meta_dt %>% filter(deadbody != 'DB1' & !(Memo %in% c('serial_daughter','half_splitted_b'))) %>% group_by(tissue_id) %>% count() %>% pull(n) %>% summary()
  r_meta_dt %>% filter(deadbody != 'DB1' & !(Memo %in% c('serial_daughter','half_splitted_b'))) %>% pull(Memo) %>% unique()
  r_meta_dt %>% filter(deadbody != 'DB1' & !(Memo %in% c('serial_daughter','half_splitted_b','bulk','other_person','dupicated_sequencing','diverge_during_culture'))) %>% nrow()
  r_meta_dt %>% filter(Memo == 'bulk') %>% pull(seqz_cov) %>% summary()
  tot_seq_dt <- r_meta_dt %>% filter(deadbody != 'DB1' & !(Memo %in% c('serial_daughter','half_splitted_b','bulk','other_person','duplicated_sequencing','diverge_during_culture','multiclonal,poor_quality'))) 
  nrow(tot_seq_dt)
  tot_seq_dt %>% pull(tissue_id) %>% unique() %>% length()
  tot_seq_dt %>% group_by(current_final_for_lineage) %>% count()
  tot_seq_dt %>% filter(current_final_for_lineage == 'N') %>% pull(Memo) %>% unique()
  tot_seq_dt %>% filter(current_final_for_lineage == 'N' & Memo == 'low_coverage') %>% nrow()
  tot_seq_dt %>% filter(current_final_for_lineage == 'N' & Memo %in% c('very_lowVAF','lowVAF')) %>% nrow()

}


###################Plotting phylogenetic tree ####################
#Draw individual phylogenetic trees
if(F){
  this_db='DB8'
  Draw_Phylo_Full_DB(this_db)
  }

##Total Tree overview
if(F){
  Draw_totalTree_Overview(meta_dt, lng_dt, nwk_tbl)
  }
#Draw legend
if(F){
  plot.new()
  plot(x=1)
  legend(x=0.8,y=1.2,legend=c("Early","Early to late","Late"),lty=c(1,1,1),col=timing_pal, title=expression(bold("Mutation classification")), cex=1.5, lwd=2)
}

#Total Early Tree overview(coloring at multifurcation)
if(F){
  Draw_totalEarly_Tree_BranchN(meta_dt, lng_dt, nwk_tbl)
}

#Total Early Tree overview (marking bifurcation and multifurcation)
if(F){
  plist <- list()
  n=0
  for (db in c('DB3','DB6','DB8','DB9','DB10')){
    print(db)
    n = n+1
    tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
    m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
    m_lng_dt <- lng_dt %>% filter(deadbody == db) %>% mutate(l_id = lineage_id) %>% rename(taxa = lineage_id) %>% select(-deadbody)
    m_lng_dt <- m_lng_dt %>% mutate(mother_id = apply(.["lineage_id2"], 1, function(x) {
      t_list = unlist(strsplit(x, '-'))
      t_db = unlist(strsplit(x,'_'))[1]
      if (length(t_list) == 1){
        paste0(t_db, '_L0')
      } else{
        paste(t_list[1:length(t_list)-1], collapse='-')
      }} ))
    m_lng_dt <- m_lng_dt %>% mutate(multifurcating = ifelse(is.na(early_down_node_n) == T | early_down_node_n < 2, NA, ifelse(branch_n ==2, 'D','P')))
    m_lng_dt <- m_lng_dt %>% mutate(sister_group = apply(.[c("mother_id","multifurcating", "sister_n")], 1, function(x){
      if (is.na(.$multifurcating[.$lineage_id2 == x[1]])){
        'Last'
      } else if (x[3] == 2){
        "D"
      } else{
        "P"
      }
    } )) 
    m_lng_dt <- m_lng_dt %>% mutate(node_group = ifelse(is.na(multifurcating), NA, ifelse(branch_n == 2, 'D','P')))
    m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
    bin <- ifelse(max(m_lng_dt$end_mtn) > 30000, 10000, ifelse(max(m_lng_dt$end_mtn) > 15000, 5000, 2000))
    plist[[n]] <- ggtree(tree, aes(color = sister_group)) %<+% m_lng_dt +theme_tree2() +
      geom_point(aes(color = node_group))+
      coord_cartesian(xlim = c(-5,35))+
      scale_color_manual(values=c('D'=pal_combi[12],'P'=pal_combi[24], "Last" = "gray"))+
      ggtitle(db)+
      theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  }
  rplot <- plot_grid(plotlist = list(plist[[1]], plist[[3]], plist[[4]], plist[[5]]), nrow=2)
  plot_grid(plotlist = list(plist[[2]], rplot), nrow=1, rel_widths = c(1,1.5))
}



# C>TG proportion in phylogenetic tree
if(F){
  Draw_totalTree_CtoTpGprop(meta_dt, lng_dt, nwk_tbl)
  }
# legend
if(F){
  colfunc <- colorRampPalette(c("yellow", "blue"))
  legend_image <- as.raster(matrix(colfunc(20), ncol=1))
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
  text(x=0.5, y = seq(0,1,0.25), labels = paste0(seq(0,100,25),'%'), cex=2.5)
  rasterImage(legend_image, 0, 0, 0.3,1)
}

# CNV presence in phylogenetic tree
if(F){
  Draw_totalTree_CNV(meta_dt, lng_dt, nwk_tbl)
  }

# SV number in phylogenetic tree
if(F){Draw_totalTree_SV(meta_dt, lng_dt, nwk_tbl)}




#####################Features of early embryonic mutations

#draw plot to determine range of trunk
if(F){
  m_lng_dt <- lng_dt %>% mutate(end_branch = ifelse(branch_n == 0, 'Y', 'N'))
  ggplot(m_lng_dt, aes(x=deadbody, y=end_mtn))+
    geom_jitter(aes(color=end_branch), height=0, width=0.3, alpha=0.5, size=3)+
    coord_cartesian(ylim=c(0,10000))+
    ylab("Position of the end of branch")+
    ggtitle('High values are skipped')+
    theme_bw()
}

#Definition of early mutations
if(F){
  f_lng_dt <- lng_dt %>% filter(end_mtn < 1000)
  ggplot(lng_dt, aes(x=end_mtn))+
    geom_histogram(aes(fill=deadbody), binwidth=10)+
    coord_cartesian(xlim=c(0,10000))+
    scale_fill_manual(values = db_pal)+
    theme_bw()+
    xlab("No. of mutations in phylogenetic tree")+
    ylab("No. of nodes")
}

#correlation between n_earlySNV and n_earlyIndels
if(F){
  ggplot(DB_dt, aes(x=early_snv_n , y=early_indel_n))+
    geom_point(aes(color = deadbody), size=5, alpha=0.7)+
    geom_smooth(method="lm")+
    geom_text_repel(aes(label = deadbody))+
    scale_color_manual(values = db_pal)+
    scale_x_continuous(limits = c(0,100))+
    scale_y_continuous(limits = c(0, 20))+
    xlab("No. of early SNVs")+ylab("No. of early indels")+
    theme_bw()
  
  m_DB_dt <- DB_dt %>% select(deadbody, sample_n, early_snv_n, early_indel_n) %>% gather(-deadbody,-sample_n, key="mutation_type", value="number")
  ggplot(m_DB_dt, aes(x=sample_n, y=number))+
    geom_point(aes(shape=mutation_type, color=deadbody),size=5, alpha=0.7)+
    geom_smooth(aes(group = mutation_type), method="lm", se = T)+
    scale_color_manual(values=db_pal)+
    scale_y_continuous(breaks=seq(0,220,20))+
    scale_shape_manual(values = c(early_indel_n = 17, early_snv_n = 19))+
    xlab("No. of samples")+ylab("No. of detected early mutations")+
    theme_bw()+theme(axis.text= element_text(size=15), axis.title = element_text(size=18))
}

#early vs late signature composition, merged
if(F){
  #substitutions
  file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/sig_subs_SBS1518/', pattern = 'DB368910.*exposures.tsv$', full.names = T)
  e_dt <- read_tsv(file_list[1])
  e_dt <- e_dt %>% rename(Signature = `Signature #`) %>% mutate(Prop = Exposure/sum(Exposure), group = 'early')
  l_dt <- read_tsv(file_list[2])
  l_dt <- l_dt %>% rename(Signature = `Signature #`) %>% mutate(Prop = Exposure/sum(Exposure), group = 'late')
  snv_dt <- bind_rows(e_dt,l_dt)
  snv_dt$Signature = factor(snv_dt$Signature, levels = c('SBS1','SBS5','SBS18'))
  ggplot(snv_dt, aes(x=group, y=Prop, fill=group))+
    geom_bar(stat= "identity")+
    scale_fill_manual(values = c('early' = pal_combi[8], 'late' = pal_combi[4]))+
    theme_syp+theme(axis.title.x = element_blank())+
    facet_grid(cols=vars(Signature))+
    ylab("Proportion")+
    ggtitle("Substitutions")
  #indel
  file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/sig_indel_124589/', pattern = 'DB368910.*exposures.tsv$', full.names = T)
  e_dt <- read_tsv(file_list[1])
  e_dt <- e_dt %>% rename(Signature = `Signature #`) %>% mutate(Prop = Exposure/sum(Exposure), group = 'early')
  l_dt <- read_tsv(file_list[2])
  l_dt <- l_dt %>% rename(Signature = `Signature #`) %>% mutate(Prop = Exposure/sum(Exposure), group = 'late')
  indel_dt <- bind_rows(e_dt,l_dt)
  indel_dt$Signature = factor(indel_dt$Signature, levels = c('ID1','ID2','ID4','ID5','ID8','ID9'))
  ggplot(indel_dt, aes(x=group, y=Prop, fill=group))+
    geom_bar(stat= "identity")+
    scale_fill_manual(values = c('early' = pal_combi[8], 'late' = pal_combi[4]))+
    theme_syp+theme(axis.title.x = element_blank())+
    facet_grid(cols=vars(Signature))+
    ylab("Proportion")+ ggtitle("Indels")
}


#early vs late signature comparison per deadbody
if(F){
  merged_dt <- tibble()
  for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
    e_dt <- read_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/sig_subs_SBS1518_v401/',this_db, '_early_SNVs.txt.signature_exposures.tsv')) 
    e_dt <- e_dt %>% rename(Signature = `Signature #`) %>% mutate(Prop = Exposure/sum(Exposure), deadbody = this_db, group='early', type='SNV') %>% 
      select(deadbody, Signature, Prop, group, type)
    l_dt <- read_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/sig_subs_SBS1518_v401/',this_db, '_lateNonUV_SNVs.txt.signature_exposures.tsv'))
    l_dt <- l_dt %>% rename(Signature = `Signature #`) %>% mutate(Prop = Exposure/sum(Exposure), deadbody = this_db, group = 'late', type='SNV') %>%
      select(deadbody, Signature, Prop, group, type)
    e2_dt <- read_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/sig_indel_124589/',this_db, '_early_indels.txt.signature_exposures.tsv')) 
    e2_dt <- e2_dt %>% rename(Signature = `Signature #`) %>% mutate(Prop = Exposure/sum(Exposure), deadbody = this_db, group='early', type='Indel') %>% 
      select(deadbody, Signature, Prop, group, type)
    l2_dt <- read_tsv(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/sig_indel_124589/',this_db, '_lateNonUV_indels.txt.signature_exposures.tsv'))
    l2_dt <- l2_dt %>% rename(Signature = `Signature #`) %>% mutate(Prop = Exposure/sum(Exposure), deadbody = this_db, group = 'late', type='Indel') %>%
      select(deadbody, Signature, Prop, group, type)
    merged_dt <- bind_rows(merged_dt, e_dt) %>% bind_rows(l_dt) %>% bind_rows(e2_dt) %>% bind_rows(l2_dt)
  }
  ggplot(merged_dt, aes(x=group, y=Prop, color = deadbody))+
    geom_point()+
    geom_line(aes(group = deadbody))+
    facet_grid(cols = vars(Signature), rows=vars(type), scales = "free_y")
}

#CIRCOS rainfall plot
if(F){
  pmt_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/DB_total_early_pointmts.txt', col_types = cols(`#CHROM`='c'))
  pmt_dt <- pmt_dt %>% rename(CHR = `#CHROM`)
  pmt_dt$DB_id <- apply(pmt_dt["lineage_id2"], 1, function(x) unlist(strsplit(x,'_'))[1])
  pmt_dt$color <- db_pal[pmt_dt$DB_id]
  pmt_dt <- pmt_dt %>% mutate(shape = ifelse(nchar(REF) == 1 & nchar(ALT) == 1, 19, 17)) 
  pmt_dt %>% arrange(CHR, POS) %>% filter(CHR=='22') %>% View()
  source('/home/users/sypark/01_Python_files/circos_R/CIRCOS_R_SYPark.R')
  circos_VCS(SNV=pmt_dt, SNV_rainfall = T, SNV_height=0.3)
}

#C>T proportion comparison 
if(F){
  ct_tbl <- lng_dt %>% mutate(SBS7_prop = SBS7a_prop + SBS7b_prop + SBS7c_prop + SBS7d_prop,
                              time_of_lineage2 = ifelse(time_of_lineage == 'early', 'early', ifelse(SBS7_prop < 5, 'late_without_UVdamage','late_withUVdamage'))) %>% filter(is.na(n_CtoTG)==F) %>% 
    group_by(deadbody, time_of_lineage2) %>% summarise(sum_CtoTG = sum(n_CtoTG), sum_snv = sum(n_subs)) %>% mutate(CTGprop = sum_CtoTG/sum_snv)
  library(ggsignif)
  ggplot(ct_tbl, aes(x=time_of_lineage2, y=CTGprop))+
    geom_boxplot(outlier.size=-1)+
    geom_signif(comparisons = list(c("early","late_without_UVdamage"), c("early","late_withUVdamage")), map_signif_level = T, y_position = c(0.375,0.4))+
    geom_jitter(aes(color=deadbody), height=0, width=0.1, size=5, alpha=0.7)+
    scale_color_manual(values = db_pal)+
    ylab("CtoTpG substituion/total subtitution")+
    theme_bw() + theme(axis.title.x = element_blank())
}


###############################mutation rate at early embryogenesis############################
#Frequency of multifurcation
if(F){
  f_dt <- lng_dt %>% filter(deadbody %in% c('DB3','DB6','DB8','DB9','DB10') & time_of_lineage %in% c('early', 'germline') & step_n < 4)
  ct_tbl <- left_join(f_dt %>% group_by(deadbody, step_n, branch_n) %>% count(),  f_dt %>% group_by(deadbody, step_n) %>% summarise(total_n = n()))
  ggplot(ct_tbl, aes(x=step_n, y=n, fill=as.factor(branch_n)))+
    geom_bar(stat="identity")+
    facet_wrap(~deadbody)
  ggplot(ct_tbl, aes(x=step_n, y=n/total_n, fill=as.factor(branch_n)))+
    geom_bar(stat="identity")+
    facet_wrap(~deadbody)
  
  ct_tbl <- left_join(f_dt %>% group_by(step_n, branch_n) %>% count(),  f_dt %>% group_by(step_n) %>% summarise(total_n = n()))
  ct_tbl <- ct_tbl %>% mutate(stage_in_tree = step_n + 1, group = ifelse(branch_n ==2, 'dichotomy', 'polytomy'))
  ggplot(ct_tbl, aes(x=stage_in_tree, y=n/total_n, fill=group))+
    geom_bar(stat="identity")+
    scale_fill_manual(values = c('dichotomy' = pal_combi[19], 'polytomy' = pal_combi[24]))+
    ylab("Proportion") + xlab("Stage in phylogenetic tree")+
    theme_syp
    
  print(this_db)
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  m_lng_dt <- lng_dt %>% filter(deadbody == this_db) %>% rename(taxa = lineage_id) %>% select(-deadbody) %>% mutate(lid = taxa)
  m_lng_dt <- m_lng_dt %>% mutate(is_multifurc = ifelse(sister_n > 2, 'Y', 'N'))
  ggtree(tree, aes(color=is_multifurc)) %<+% m_lng_dt +
    coord_cartesian(xlim = c(-1,35))+
    scale_color_manual(values=c('Y' = 'red', 'N' = 'black'))
  }

# No. of mutations per cell per cell doubling
if(F){
  f_dt <- lng_dt %>% filter(time_of_lineage == 'early') %>% mutate(MPCPCD=n_pointmt/(branch_n -1)) %>% mutate(is_step1 = ifelse(step_n == 1, 'Y', 'N'))
  
  rpois_v <- rpois(10000, lambda=1)
  rpois_dt <- as.tibble(table(rpois_v))
  rpois_dt <- rpois_dt %>% mutate(count=as.numeric(rpois_v)) %>% select(-rpois_v)
  rpois_dt$prop <- rpois_dt$n/sum(rpois_dt$n)
  rpois_dt
  
  ggplot(f_dt, aes(x=MPCPCD))+
    geom_bar(data=rpois_dt, aes(x=count,y=prop),stat="identity")+
    geom_density()
  
  ggplot(f_dt, aes(x=MPCPCD))+
    geom_bar(data=rpois_dt, aes(x=count,y=prop),stat="identity")+
    geom_density()+
    facet_wrap(~deadbody)
  quantile(f_dt$MPCPCD)
  f_dt %>% group_by(deadbody) %>% summarise(med_v = median(MPCPCD), q1 = quantile(MPCPCD)[2], q3 = quantile(MPCPCD)[4])
}
  
#plotting nodal stage and n_nodes filled with n_branch
if(F){
  p <- list()
  n=0
  for (db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
    n = n+1
    m_res_tbl <- res_tbl %>% filter(deadbody == db & branch_n >0) %>% group_by(step_n, branch_n) %>% count() 
    p[[n]] <- ggplot(m_res_tbl, aes(x=step_n+1, y=n))+
      geom_bar(aes(fill = as.factor(branch_n)), stat="identity") +
      scale_x_continuous(breaks=seq(0,10,1))+
      scale_fill_manual(values=num_pal)+
      xlab("Nodal stage")+ylab("No. of nodes")+ggtitle(db)+
      theme_bw(base_size=12)+theme(legend.position = 'none')
  }
  lplot <- plot_grid(plotlist=p[c(1,2,3,5,6,7)], nrow=2, align="hv")
  plot_grid(plotlist=list(lplot, p[[4]]), nrow=1, rel_widths = c(1.5,1))
  plot.new()
  plot(x=1)
  legend("topleft", legend=names(num_pal), fill=num_pal,
         bty="n", border="white", cex=2, ncol=length(num_pal)) 
}

#plotting nodal stage and proportion of n_branch
if(F){
  p <- list()
  n=0
  for (db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
    n = n+1
    m_res_tbl <- res_tbl %>% filter(deadbody == db & branch_n >0) %>% group_by(step_n, branch_n) %>% summarise(count = n()) %>% mutate(perc = count/sum(count))
    p[[n]] <- ggplot(m_res_tbl, aes(x=step_n+1, y=perc*100))+
      geom_bar(aes(fill = as.factor(branch_n)), stat="identity") +
      scale_x_continuous(breaks=seq(1,8,1))+
      scale_fill_manual(values=num_pal)+
      xlab("Nodal stage")+ylab("No. of nodes")+ggtitle(db)+
      theme_bw(base_size=12)+theme(legend.position = 'none')
  }
  lplot <- plot_grid(plotlist=p[c(1,2,3,5,6,7)], nrow=2, align="hv")
  plot_grid(plotlist=list(lplot, p[[4]]), nrow=1, rel_widths = c(1.5,1))
}

#plotting nodal stage and n_pointmt
if(F){
  res_tbl <- lng_dt %>% filter(end_mtn < 300 & deadbody %in% c('DB3','DB6','DB8','DB9','DB10'))
  stat_tbl <- res_tbl %>% group_by(step_n) %>% summarise(mean_pointmt = mean(n_pointmt), median_pointmt = median(n_pointmt))
  stat_tbl
  ggplot(res_tbl, aes(x=as.factor(step_n), y=n_pointmt))+
    geom_violin()+
    geom_jitter(aes(color = deadbody), size=3, alpha=0.7, width=0.2, height=0.2)+
    #geom_point(data = stat_tbl, aes(x=as.factor(step_n), y=median_pointmt), size=5, color = "red")+
    coord_cartesian(xlim=c(2,7))+
    scale_color_manual(values = db_pal)+
    theme_bw()+
    xlab("Stage of branch")+ylab("No. of point mutations")+
    ggtitle("Red dot = median value")+
    facet_wrap(~deadbody, scales = "free_x")
  ggplot(subset(lng_dt,time_of_lineage == 'early'), aes(x=start_mtn, y=n_pointmt))+
    #geom_violin()+
    geom_point(aes(color = deadbody), size=3, alpha=0.7, width=0.2, height=0.2)+
    #geom_point(data = stat_tbl, aes(x=as.factor(step_n), y=median_pointmt), size=5, color = "red")+
    #coord_cartesian(xlim=c(2,7))+
    scale_color_manual(values = db_pal)+
    theme_bw()+
    xlab("Stage of branch")+ylab("No. of point mutations")+
    ggtitle("Red dot = median value")+
    facet_wrap(~deadbody, scales = "free_x")
}

#plot mutation skipping rate
if(F){
  p <- list()
  n =0
  for (db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
    n=n +1
    tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
    m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% mutate(l_id = taxa) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
    m_lng_dt <- lng_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id) %>% select(-deadbody)
    m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
    m_lng_dt$early_mt_skipping_rate[m_lng_dt$early_desc_n < 3] <- NA
    #m_lng_dt <- left_join(m_lng_dt, res_tbl %>% filter(deadbody == db) %>% rename(taxa=lineage_id), by='taxa')
    #total_skip_rate <- round(res_tbl$early_mt_skipping_rate[res_tbl$deadbody == db & res_tbl$lineage_id == 'L0'],2)
    #fill data of nodes if information of it's branches is unique
    node_dt <- m_lng_dt
    p[[n]] <- ggtree(tree) %<+% node_dt + geom_tree()+theme_tree2() + 
      geom_label(aes(label=round(early_mt_skipping_rate,2), fill=early_mt_skipping_rate),hjust=1.1, alpha=0.7, size=2)+
      scale_fill_gradient2(low= "blue", mid = "yellow", high="red", limits=c(0,1), midpoint=0.5)+
      coord_cartesian(xlim = c(-2,30))+
      scale_x_continuous(breaks=seq(0,30,5))+
      ggtitle(db)+
      theme(axis.text = element_text(size=10), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0.2),"cm"))
  }
  lplot <- plot_grid(plotlist=p[c(1,2,3,5,6,7)], nrow=2, align="hv")
  plot_grid(plotlist=list(lplot, p[[4]]), nrow=1, rel_widths = c(1.5,1))
  #legend
  colfunc <- colorRampPalette(c("red", "yellow", "blue"))
  plot(1:2, 1:2, pch = 19, cex=2, col = colfunc(20))
  legend_image <- as.raster(matrix(colfunc(20), ncol=1))
  plot(c(0,2),c(0,2),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
  text(x=0.15, y = seq(0,1,0.5), labels = seq(0,1,0.5), cex=1.5)
  rasterImage(legend_image, 0, 0, 0.1,1)
}

#branching number in phylogenetic tree
if(F){
  p <- list()
  n =0
  for (db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
    n=n +1
    tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
    m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% mutate(l_id = taxa) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
    m_lng_dt <- lng_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id) %>% select(-deadbody)
    m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
    p[[n]] <- ggtree(tree) %<+% m_lng_dt + geom_tree()+theme_tree2() + 
      geom_label(aes(label=branch_n, fill=branch_n),hjust=1.1, alpha=0.7, size=3)+
      scale_fill_gradient2(low= "blue", mid = "yellow", high="red", limits=c(0,7), midpoint=3)+
      coord_cartesian(xlim = c(-2,30))+
      scale_x_continuous(breaks=seq(0,30,5))+
      ggtitle(db)+
      theme(axis.text = element_text(size=10), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0.2),"cm"))
  }
  lplot <- plot_grid(plotlist=p[c(1,2,3,5,6,7)], nrow=2, align="hv")
  plot_grid(plotlist=list(lplot, p[[4]]), nrow=1, rel_widths = c(1.5,1))v
}


#Simulation result summary
if(F){
  mr_tbl <- tribble(~deadbody,~low_MR,~MR,~high_MR,
                    'DB2',NA,NA,NA,
                    'DB3',0.53,0.79,1.12,
                    'DB5',NA,NA,NA,
                    'DB6',0.96,1.11,1.35,
                    'DB8',0.95,1.27,1.65,
                    'DB9',1.42,1.89,2.49,
                    'DB10',NA,NA,NA)
  mr_tbl$deadbody <- factor(mr_tbl$deadbody, levels = c('DB10','DB9','DB8','DB6','DB5','DB3','DB2'))
  ggplot(mr_tbl, aes(x=deadbody, y=MR))+
    geom_point(size=5)+
    geom_errorbar(aes(ymin=low_MR, ymax=high_MR), width=0.2)+
    scale_y_continuous(limits=c(0,3), breaks =seq(0,3,0.5))+
    theme_bw()+theme(axis.title.y = element_blank(), axis.text = element_text(size=15), axis.title.x = element_text(size=18))+
    ylab("Mutation rate (No. of point mutations/cell/division)")+
    coord_flip()
}


# mutatition per cell generation per L1 or L2
if(F){
  m_lng_dt <- lng_dt
  m_lng_dt$high_group <- apply(lng_dt["lineage_id"], 1, function(x) unlist(strsplit(x,'-'))[1])
  m_lng_dt <- m_lng_dt %>% filter(time_of_lineage == 'early' & deadbody %in% c('DB3','DB6','DB8','DB9','DB10')) %>% mutate(n_pmt_pcg = n_pointmt/branch_n)
  m_lng_dt$deadbody = factor(m_lng_dt$deadbody, levels = c('DB3','DB6','DB8','DB9','DB10'))
  ggplot(m_lng_dt, aes(x=deadbody, y=n_pmt_pcg))+
    geom_boxplot(aes(fill = high_group))+
    theme_bw()+labs(fill = "Lineage group")+
    ylab("No. of mutation per cell generation")+
    theme(axis.title.x = element_blank(), axis.text = element_text(size=15), axis.title = element_text(size=18))
  for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
    print(this_db)
    m_lng_dt %>% filter(deadbody == this_db & high_group == 'L1') %>% .$n_pmt_pcg -> a
    m_lng_dt %>% filter(deadbody == this_db & high_group == 'L2') %>% .$n_pmt_pcg -> b
    print(wilcox.test(a,b)$p.value)
  }
}


##################################Asymmetry###################################
#Asymmetic distribution of early lineages
#plot asymmetry rate on tree
if(F){
  Draw_totalTree_AsymRatio(meta_dt, lng_dt, nwk_tbl)
  Make_Bifurc_Asymmetry_PNG(lng_dt, save_path = "/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Asymmetry_png/")
  Draw_indiTree_AsymPNG(meta_dt, lng_dt, nwk_tbl,  "/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Asymmetry_png/", 'DB10', width=2, height=0.75)
  Draw_gradient_legend3(color_v=c("gray50","yellow","red"), start_n=1, change_n=2, end_n=7)
  Draw_Dist_Asym_Obs_Sim(lng_dt, sim_res_path ='/home/users/sypark/00_Project/06_LineageTracing/Rscripts/asymmetry_simulation_1e5_objects.rds')
  }

#plot asymmetry ratio with scatter plots, branch length comparison of dichotomies
if(F){
  subdt <- lng_dt %>% filter(deadbody %in% c('DB3','DB6','DB8','DB9','DB10') & branch_n == 2 & early_desc_n >=6) 
  res_tbl = tibble(lineage_id2 = as.character(), n1 = as.numeric(), n2= as.numeric(), ratio = as.numeric(), 
                   length1 = as.numeric(), length2 = as.numeric(), group1 = as.character(), group2 = as.character())
  i=0
  for (t_lid in subdt$lineage_id2){
    i = i+1
    print(t_lid)
    did1 = gsub('0-','',paste0(t_lid,'-1'))
    did2 = gsub('0-','',paste0(t_lid,'-2'))
    n1 = lng_dt %>% filter(grepl(did1, lineage_id2) & time_of_lineage == 'early_to_late') %>% nrow()
    n2 = lng_dt %>% filter(grepl(did2, lineage_id2) & time_of_lineage == 'early_to_late') %>% nrow()
    t_ratio = max(n1, n2)/min(n1,n2)
    branch1 = lng_dt %>% filter(lineage_id2 == did1) %>% pull(n_pointmt)
    branch2 = lng_dt %>% filter(lineage_id2 == did2) %>% pull(n_pointmt)
    group1 = lng_dt %>% filter(lineage_id2 == did1) %>% pull(time_of_lineage)
    group2 = lng_dt %>% filter(lineage_id2 == did2) %>% pull(time_of_lineage)
    res_tbl[i,] <- list(t_lid, n1, n2, t_ratio, branch1, branch2, group1, group2)
  }
  res_tbl <- res_tbl %>% mutate(deadbody = apply(.["lineage_id2"], 1, function(x) unlist(strsplit(x,'_'))[1]))
  res_tbl$deadbody <- factor(res_tbl$deadbody, levels = c('DB3','DB6','DB8','DB9','DB10'))
  ggplot(res_tbl, aes(x=deadbody, y=ratio, fill = deadbody))+
    geom_violin()+
    geom_jitter(size=3, alpha=0.5,height=0, width=0.05)+
    geom_hline(yintercept = 1, linetype = 'dashed')+
    scale_fill_manual(values = db_pal)+
    scale_y_continuous(breaks = c(1,seq(5,25,5)))+
    ylab("Asymmetric ratio at dichotomy")+
    theme_syp + theme(legend.position = 'none', axis.title.x = element_blank())
  
  res_tbl %>% filter(ratio >=2 & group1 == 'early' & group2 == 'early') %>% nrow()
  f_res_tbl <- res_tbl %>% filter(ratio >=1.5 & group1 == 'early' & group2 == 'early')
  #f_res_tbl <- res_tbl %>% filter(ratio >=2 & group1 == 'early' & group2 == 'early')
  f_res_tbl %>% filter((n1 < n2 & length1 > length2) | (n1 > n2 & length1 < length2))
  f_res_tbl %>% filter(!((n1 < n2 & length1 > length2) | (n1 > n2 & length1 < length2)))
  
}


#Total Early Tree overview (marking asymmetric bifurcation)
if(F){
  plist <- list()
  n=0
  for (db in c('DB3','DB6','DB8','DB9','DB10')){
    print(db)
    n = n+1
    tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
    
    tmp_tbl1 <- f_res_tbl %>% mutate(small_portion = ifelse(n1 < n2, 'first','second')) %>% select(lineage_id2, small_portion) %>%
      mutate(small_id = ifelse(small_portion == 'first', gsub('0-','',paste0(lineage_id2, '-1')), gsub('0-','',paste0(lineage_id2,'-2')))) %>%
      mutate(large_id = ifelse(small_portion == 'first', gsub('0-','',paste0(lineage_id2, '-2')), gsub('0-','',paste0(lineage_id2,'-1')))) %>% select(small_id, large_id) %>%
      gather(key='group', value='lineage_id2')
    
    tmp_tbl2 <- f_res_tbl %>% mutate(node_group = ifelse(length1 == length2, 'same', 
                                             ifelse((n1 < n2 & length1 > length2) | (n1 > n2 & length1 < length2), 'small_long', 'large_long'))) %>%
      select(lineage_id2, node_group)
    tmp_tbl <- full_join(tmp_tbl1, tmp_tbl2) %>% separate(lineage_id2, c('deadbody', 'lineage_id'), sep='_', remove = F)
    m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
    m_lng_dt <- left_join(lng_dt %>% select(deadbody, lineage_id, lineage_id2),tmp_tbl)
    m_lng_dt$group[is.na(m_lng_dt$group)] <- 'other'
    m_lng_dt <- m_lng_dt %>% filter(deadbody == db) %>% mutate(l_id = lineage_id) %>% rename(taxa = lineage_id) %>% select(-deadbody)
    m_lng_dt <- m_lng_dt %>% select(taxa, l_id, lineage_id2, group, node_group)
    bin <- ifelse(max(m_lng_dt$end_mtn) > 30000, 10000, ifelse(max(m_lng_dt$end_mtn) > 15000, 5000, 2000))
    plist[[n]] <-    ggtree(tree, aes(color = group)) %<+% m_lng_dt +theme_tree2() +
      geom_point(aes(color = node_group))+
      coord_cartesian(xlim = c(-5,35))+
      scale_color_manual(values=c('small_id'='red', "large_id" = "blue", 'mother'='darkgray', 'other'='darkgray', 'same'='black', 'small_long'='red','large_long'='blue'))+
      ggtitle(db)+
      theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  }
  rplot <- plot_grid(plotlist = list(plist[[1]], plist[[3]], plist[[4]], plist[[5]]), nrow=2)
  plot_grid(plotlist = list(plist[[2]], rplot), nrow=1, rel_widths = c(1,1.5))
}


######## culture associated mutations

#define culture-associated mutations
file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/merged_data/culture_associated_mutations/point_mutations', pattern = '.scF.scF$', full.names = T)
for (n in 1:length(file_list)){
  dt <- read_tsv(file_list[n], col_types = cols(`#CHROM`='c'))
  t_id <- unlist(strsplit(unlist(strsplit(file_list[n],'point_mutations/'))[2],'\\.mother'))[1]
  print(t_id)
  colnames(dt)
  g <- ggplot(dt, aes(x=mt_CCF, y=dt_CCF))+
    geom_point(size=3, alpha=0.1)+
    xlab("Mother CCF")+ylab("Daughter CCF")+
    ggtitle(t_id) + theme_syp
  print(g)
  CA_dt <- dt %>% filter(mt_CCF < 0.25 & dt_CCF >= 0.25 & dt_CCF < 1.75) 
  nrow(CA_dt) %>% print() #point_mt 491 / 865
  CA_dt %>% filter(nchar(REF) == 1 & nchar(ALT)==1) %>% nrow() %>% print() #SNV 476
  CA_dt %>% filter(nchar(REF) != 1 | nchar(ALT)!=1) %>% nrow() %>% print() #indel 15
  CA_dt %>% write_tsv(paste0('/home/users/sypark/00_Project/06_LineageTracing/merged_data/culture_associated_mutations/point_mutations/',t_id,'.CA_mutations.txt'))
}

#signature of culture-associated mutations
sig <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/culture_associated_mutations/point_mutations/sig_SNV_SBS1518_v200/Three_samples_merged.CA_mutations.txt.signature_exposures.tsv')
sig %>% mutate(prop = Exposure*100/sum(Exposure))



