stop("a")

#load libraries
library(tidyverse)
library(ggtree)
library(ggsci)
library(scales)
library(cowplot)
library(ggrepel)
library(ggpubr)

#file paths
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
if(F){ #previous files
  targetseq_path_tbl <- tribble(
    ~deadbody, ~folder_path, ~Call_name, 
    "DB3",'/home/users/sypark/00_Project/06_LineageTracing/db3/04-5_Target_seq/05_re_processing/','DB3_TargetSeq_Target_list.GSNPadded.60sCall',
    "DB6",'/home/users/sypark/00_Project/06_LineageTracing/db6/14_TargetSeq/','DB6_TargetSeq_Target_list.GSNPadded.143sCall',
    "DB8",'/home/users/sypark/00_Project/06_LineageTracing/db8/11_TargetSeq/','DB8_TargetSeq_Target_list.GSNPadded.44sCall',
    "DB9",'/home/users/sypark/00_Project/06_LineageTracing/db9/11_TargetSeq/','DB9_TargetSeq_Target_list.GSNPadded.43sCall',
    "DB10",'/home/users/sypark/00_Project/06_LineageTracing/db10/11_TargetSeq/','DB10_TargetSeq_Target_list.GSNPadded.89sCall'
  )
}


targetseq_path_tbl <- tribble(
  ~deadbody, ~folder_path, ~Call_name, 
  "DB3",'/home/users/sypark/00_Project/06_LineageTracing/db3/04-5_Target_seq/05_re_processing/','DB3_TargetSeq_Target_list.GSNPadded.60sCall',
  "DB6",'/home/users/sypark/00_Project/06_LineageTracing/db6/14_TargetSeq/','DB6_TargetSeq_Target_list.merged.GSNPadded.143sCall',
  "DB8",'/home/users/sypark/00_Project/06_LineageTracing/db8/11_TargetSeq/','DB8_TargetSeq_Target_list.merged.GSNPadded.44sCall',
  "DB9",'/home/users/sypark/00_Project/06_LineageTracing/db9/11_TargetSeq/','DB9_TargetSeq_Target_list.merged.GSNPadded.43sCall',
  "DB10",'/home/users/sypark/00_Project/06_LineageTracing/db10/11_TargetSeq/','DB10_TargetSeq_Target_list.merged.GSNPadded.89sCall'
)


#load data
pmeta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_TargetSeq_201023.txt')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_200908.txt')
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_200908.txt') %>% filter(current_final_for_lineage == 'Y')
meta_dt <- meta_dt %>% mutate(Source3 = ifelse(is.na(Source_side)==T, Source2, paste0(Source_side,'_',Source2)))
meta_dt$deadbody <- factor(meta_dt$deadbody, levels =c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))

#color setting
db_pal <- pal_npg('nrc')(10)[c(1:5,9,7)]
names(db_pal) <- unique(meta_dt$deadbody)
pmeta_dt$anatomy_class3 %>% unique()
#show_col(c(pal_npg("nrc")(10), pal_d3('category10')(10)))
pal_combi <- c(pal_npg("nrc")(10), pal_d3('category10')(10), pal_aaas('default')(10))
#show_col(pal_combi)
#show_col(pal_npg("nrc")(10)[c(2,4,6,1,5,8,3,6,9)])
#show_col(pal_d3('category10')(10))
#show_col(pal_d3('category10')(10)[c(1,10,5,2,3,9,6,8)])
#anatomy_cl3_pal <-c(pal_d3('category10')(10)[c(1,10,3,9,2,4,5,6)],'black')
#names(anatomy_cl3_pal) <- c('lt_UE','lt_LE','rt_UE','rt_LE','trunk','internal_organ','HN')
#anatomy_cl3_pal <- c('#a6cee3','#1f78b4','#33a02c','#fb9a99','#e31a1c','#ff7f00','#cab2d6','#6a3d9a','#ffff99')
anatomy_cl3_pal <- pal_combi[c(24,11,2,8,12,17,3,13,18)]
names(anatomy_cl3_pal) <- c('lt_LE','lt_UE','lt_trunk','rt_LE','rt_UE','rt_trunk','center_trunk','HN','internal_organ')
#show_col(anatomy_cl3_pal)
source2_pal <- pal_d3('category10')(6)
names(source2_pal) <- c('HN','lt_LE','lt_UE','rt_LE','trunk','rt_UE')

#theme setting
theme_syp = theme_bw() + theme(axis.text = element_text(size=15), axis.title = element_text(size=18), panel.border = element_blank(), axis.line = element_line(), title = element_text(size=18))

#Sourcing functions
source('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/TargetSeq_functions.R')

#load data
if (F){
  merged_m_vaf_tbl <- tibble()
  for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
    print(this_db)
    ts_dt <- read_tsv(paste0(targetseq_path_tbl$folder_path[targetseq_path_tbl$deadbody == this_db], targetseq_path_tbl$Call_name[targetseq_path_tbl$deadbody == this_db]))
    ts_dt <- ts_dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_'))
    early_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early') %>% .$lineage_id2 
    #assign lineage_id
    merged_tbl <- Assign_Lineage_ID(ts_dt, this_db)
    nrow(ts_dt) #DB3:141, DB6:215, DB8:80
    nrow(merged_tbl) #DB3:145, DB6:215, DB8:82
    #make column list
    starting_col=c(7,7,8,8,8)
    names(starting_col) <- c('DB3','DB6','DB8','DB9','DB10')
    col_list <- colnames(merged_tbl)[starting_col[this_db]:(ncol(merged_tbl)-3)] 
    sampleN <- length(col_list) 
    sampleN #DB3:60, DB6:143, DB8:44, DB9:43, DB10:89
    
    #MAKE VAF tbl
    vaf_tbl1 <- Make_VAF_tbl(merged_tbl, this_db)
    #add expected VAF value
    vaf_tbl1 <- Add_Expected_VAF(vaf_tbl1, this_db)
    #make gathered VAF table with VAF and DP value
    m_vaf_tbl1 <- vaf_tbl1 %>% select(var_id, lineage_id2, WGS_bloodVAF, ends_with('_vaf'), expected_VAF, start_mtn, end_mtn, time_of_lineage, step_n) %>% gather(-var_id, -lineage_id2, -WGS_bloodVAF, -start_mtn, -end_mtn, -time_of_lineage, -step_n,-expected_VAF, key='tissue', value='VAF') %>%
      mutate(sample_id = gsub('_vaf','',tissue)) %>% select(-tissue)
    m_vaf_tbl2 <- vaf_tbl1 %>% select(var_id, lineage_id2, ends_with('_dp')) %>% gather(-var_id, -lineage_id2, key='tissue', value='DP') %>%
      mutate(sample_id = gsub('_dp','',tissue)) %>% select(-tissue)
    m_vaf_tbl3 <- vaf_tbl1 %>% select(var_id, lineage_id2, ends_with('_var')) %>% gather(-var_id, -lineage_id2, key='tissue', value='var') %>%
      mutate(sample_id = gsub('_var','',tissue)) %>% select(-tissue)
    m_vaf_tbl <- left_join(left_join(left_join(m_vaf_tbl1, m_vaf_tbl2), m_vaf_tbl3), pmeta_dt)
    merged_m_vaf_tbl <- bind_rows(merged_m_vaf_tbl , m_vaf_tbl %>% mutate(deadbody = this_db))
  }
  #assign expected number mutation
  rank_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Target_seq_related/early_varid_DP100n_rank.txt') %>% filter(is.na(rank) == F)
  rank_dt <- left_join(rank_dt, lng_dt %>% select(lineage_id2, start_mtn, end_mtn, n_pointmt))
  rank_dt <- left_join(rank_dt, rank_dt %>% group_by(lineage_id2) %>% summarise(targeted_n = n()))
  rank_dt <- left_join(rank_dt, rank_dt %>% group_by(lineage_id2, rank) %>% summarise(n_vars_same_rank = n()))
  rank_dt <- rank_dt %>% group_by(lineage_id2) %>% mutate(mean_rank = rank(rank), rank_ctie_first = rank(rank, ties.method = "first")) %>%
    mutate(adjusted_mean_rank = mean_rank * n_pointmt/targeted_n) %>% mutate(assigned_mean_mtn = start_mtn + adjusted_mean_rank) %>%
    mutate(adjusted_ctie_rank = rank_ctie_first * n_pointmt/targeted_n) %>% mutate(assigned_ctie_mtn = start_mtn + adjusted_ctie_rank)
  merged_m_vaf_tbl <- left_join(merged_m_vaf_tbl, rank_dt %>% 
                                  select(var_id, lineage_id2, rank, n_vars_same_rank, mean_rank, rank_ctie_first, adjusted_mean_rank, assigned_mean_mtn, adjusted_ctie_rank, assigned_ctie_mtn))
  saveRDS(merged_m_vaf_tbl, '/home/users/sypark/00_Project/06_LineageTracing/merged_data/Target_seq_related/merged_m_vaf_tbl_new.rds')
}
merged_m_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Target_seq_related/merged_m_vaf_tbl_new.rds')


#first_occuring_variants
first_variants <- tribble(~deadbody, ~lineage_id2, ~var_id,
                          'DB3','DB3_L1','1_74024235_G_T',
                          'DB3','DB3_L1','2_4923793_T_C',
                          'DB3','DB3_L1','7_32925447_G_C',
                          'DB3','DB3_L1','7_129376347_C_T',
                          'DB3','DB3_L1','4_34076338_C_A',
                          'DB3','DB3_L1','16_85372995_G_A',
                          'DB3','DB3_L2','19_19279396_T_G',
                          'DB3','DB3_L2','18_6114473_T_C',
                          'DB6','DB6_L1','13_61660346_G_A',
                          'DB6','DB6_L2','2_82379165_T_A',
                          'DB6','DB6_L2','16_16266502_G_A',
                          'DB8','DB8_L1','5_85367961_G_T',
                          'DB8','DB8_L1','8_123054768_C_T',
                          'DB8','DB8_L1','13_102309939_G_A',
                          'DB8','DB8_L1','22_17072858_C_G', 
                          'DB8','DB8_L2','1_162486606_C_A',
                          'DB9','DB9_L1','11_50359839_A_C',
                          'DB9','DB9_L2','8_58127053_G_A',
                          'DB10','DB10_L1','1_166645826_C_T',
                          'DB10','DB10_L1','4_1293246_C_T', 
                          'DB10','DB10_L1','8_93863618_C_A', 
                          'DB10','DB10_L1','15_94394073_C_A',
                          'DB10','DB10_L2','1_56148328_C_A',
                          'DB10','DB10_L2','6_44692877_C_A', 
                          'DB10','DB10_L2','8_139121031_C_T',
                          'DB10','DB10_L2','9_94617243_T_C', 
                          'DB10','DB10_L2','11_96493593_G_T',
                          'DB10','DB10_L2','15_62973855_A_C'
                          )


this_db = 'DB3'
m_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db)
#---------------------------------------------------------------script-------------------------------------------------------------------

#Tissue type
pmeta_dt %>% View()
pmeta_dt %>% filter(deadbody == "DB10") %>% group_by(deadbody, tissue_class2) %>% count() %>% View()
pmeta_dt %>% filter(deadbody == "DB10") %>% group_by(deadbody, dominant_layer2) %>% count() %>% View()


#Early_mutation detection rate
if(F){
  fmm_vaf_tbl <- merged_m_vaf_tbl %>% filter(DP >= 30 & time_of_lineage == 'early') 
  axis_labels =c(0,0.05,0.1,0.2,0.5,1.0)
  fmm_vaf_tbl %>% group_by(step_n, dominant_layer) %>% summarise(total_n=n(), detected_n=sum(var >= 3, na.rm = T)) %>% filter(dominant_layer != 'mixed') %>%
    ggplot(aes(x=step_n, y=log2(detected_n/total_n+0.1), color=dominant_layer))+
    geom_point()+
    geom_line(aes(group=dominant_layer))+
    scale_x_continuous(breaks=seq(1,10,1))+
    scale_y_continuous(labels=axis_labels, breaks=log2(axis_labels+0.1))+
    ylab("Propotion of mutation_detected_bulk_tissues")+
    theme_syp
}


#write list of varid with DP>=100 n of tissues, preparation for writing VAF drop
if(F){
  merged_m_vaf_tbl %>%filter(DP > 100 & time_of_lineage=='early' & tissue_class != 'cancer_tissue') %>% group_by(deadbody, var_id, lineage_id2) %>% count() %>%
    rename(DP100n = n) %>% arrange(lineage_id2) %>% 
    write_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Target_seq_related/early_varid_DP100n_new.tsv')
}


#Intralineage VAF drop (pairedT and bootstrap)
if(F){  #total available lids
  for(this_db in c('DB3','DB6','DB8','DB9','DB10')){
    m_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db)
    incls <- m_vaf_tbl %>% group_by(lineage_id2, var_id) %>% count() %>% filter(n >= 10) %>% pull(var_id)
    fm_vaf_tbl <- m_vaf_tbl %>% filter(DP > 100 & time_of_lineage=='early' & tissue_class != 'cancer_tissue' & var_id %in% incls)
    lid_list <- fm_vaf_tbl %>% select(lineage_id2, var_id) %>% unique() %>% group_by(lineage_id2) %>% count() %>% filter(n>1) %>% pull(lineage_id2)
    pdf(paste0(targetseq_path_tbl$folder_path[targetseq_path_tbl$deadbody == this_db], '/plot_per_lineageid_grouping_pairedT_incEcto_new.pdf'), width=8, heigh=10)
    for(lid2 in lid_list){
      print(lid2)
      print(VariantGroupingByTissueVAF(fm_vaf_tbl, lid2, 10000))
    }
    dev.off()
  }
}

#Intralineage VAF drop (pairedT and bootstrap total)
if(F){  #total available lids
  #incls <- merged_m_vaf_tbl %>% group_by(lineage_id2, var_id) %>% count() %>% filter(n >= 10) %>% pull(var_id)
  #fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(DP > 100 & time_of_lineage=='early' & tissue_class != 'cancer_tissue' & var_id %in% incls)
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(DP > 100 & time_of_lineage=='early' & tissue_class != 'cancer_tissue')
  lid_list <- fm_vaf_tbl %>% select(lineage_id2, var_id) %>% unique() %>% group_by(lineage_id2) %>% count() %>% filter(n>1) %>% arrange(lineage_id2) %>% pull(lineage_id2) 
  pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/pairedT_bootstrap_new.pdf', width=8, heigh=10)
  for(lid2 in lid_list){
    print(lid2)
    print(VariantGroupingByTissueVAF(fm_vaf_tbl, lid2, 10000))
  }
  dev.off()
}

#Intralineage VAF drop (L1 and L2 only)
if (F){
  pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/varid_medianBootStrap_step1.pdf')
  for(this_db in c('DB3','DB6','DB8','DB9','DB10')){
    m_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db)
    fm_vaf_tbl <- m_vaf_tbl %>% filter(DP > 100 & time_of_lineage=='early')
    lid_list <- fm_vaf_tbl %>% filter(step_n==1 )%>% select(lineage_id2, var_id) %>% unique() %>% group_by(lineage_id2) %>% count() %>% filter(n>1) %>% pull(lineage_id2)
    for(lid2 in lid_list){
      print(lid2)
      print(VariantGroupingByTissueVAF(fm_vaf_tbl, lid2, 10000))
    }
  }
  dev.off()
}

#Intralineage VAF paired t test
if(F){ 
  lids <- merged_m_vaf_tbl %>%  filter(DP> 100 & time_of_lineage=='early' & tissue_class != 'cancer_tissue' & dominant_layer2 !='ectoderm' & dominant_layer2 != 'ecto_mesoderm') %>% select(lineage_id2, var_id) %>% unique() %>% group_by(lineage_id2) %>% count() %>% filter(n>1) %>% pull(lineage_id2)
  pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/pairedT_DP100_noEcto_noCancer.pdf')
  for(lid in lids){
    print(lid)
    f_tbl <- merged_m_vaf_tbl %>%  filter(lineage_id2 == lid & tissue_class != 'cancer_tissue' & dominant_layer2 !='ectoderm' & dominant_layer2 != 'ecto_mesoderm') %>% select(var_id, VAF, sample_id)
    expand.grid(unique(f_tbl$var_id),unique(f_tbl$var_id))
    tmp_tbl <- merged_m_vaf_tbl %>% filter(lineage_id2 == lid) %>% select(var_id, VAF, sample_id) %>% spread(key=var_id, value=VAF)
    library(ggpubr)
    x_order <- f_tbl %>% group_by(var_id) %>% summarise(meanVAF = mean(VAF)) %>% arrange(desc(meanVAF)) %>% pull(var_id)
    mycomparison <- combn(unique(f_tbl$var_id),2, simplify = F)
    g <- ggplot(f_tbl, aes(x=var_id, y=VAF))+
      geom_boxplot()+
      geom_line(aes(group=sample_id), alpha=0.1)+
      stat_compare_means(aes(label = ..p.signif..), comparisons = mycomparison,paired=T, method='t.test')+
      scale_x_discrete(limits = x_order)+
      theme_syp+theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1))+
      ggtitle(lid)
    print(g)
  }
  dev.off()
}

if(F){ #individual plot
  tmp_dt <- merged_m_vaf_tbl %>% filter(DP>=500 & lineage_id2 == 'DB9_L1-2-1' & dominant_layer2 != 'ectoderm' &  dominant_layer2 != 'ecto_mesoderm' & tissue_class !='cancer_tissue')
  x_order <- tmp_dt %>% group_by(var_id) %>% summarise(meanVAF = mean(VAF))  %>% arrange(desc(meanVAF)) %>% pull(var_id)
  ggplot(tmp_dt, aes(x=var_id, y=VAF))+
    geom_line(aes(group=sample_id),alpha=0.1)+
    geom_boxplot()+
    geom_point()+
    scale_x_discrete(limits=x_order)+
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  merged_m_vaf_tbl %>% filter(lineage_id2 == 'DB9_L1-2-1' & dominant_layer2 != 'ectoderm' &  dominant_layer2 != 'ecto_mesoderm' & tissue_class !='cancer_tissue') %>% 
    group_by(var_id) %>% summarise(meanDP=mean(DP))
}

# VAF of L1+L2 barplot absolute
if(F){Draw_L1plusL2_VAF(merged_m_vaf_tbl)}
# VAF of L1+L2 barplot relative
if(F){Draw_L1L2_VAF_relative(merged_m_vaf_tbl)}
# VAF of L1+L2 density plot
if(F){
  m_dt <- left_join(first_variants, merged_m_vaf_tbl %>% filter(DP> 100)) %>% mutate(lineage_id = apply(.["lineage_id2"], 1, function(x) unlist(strsplit(x,'_'))[2])) %>% group_by(deadbody, lineage_id, sample_id) %>%
    summarise(medVAF = median(VAF)) %>% spread(key=lineage_id, value=medVAF) %>% mutate(L12 = L1 + L2) %>% filter(is.na(L12) == F)
  GSNP_mVAF <- merged_m_vaf_tbl %>% filter(DP > 100 & grepl('GSNP', lineage_id2)) %>% pull(VAF) %>% mean(., na.rm=T)
  ggplot(m_dt, aes(x=L12, color=deadbody))+
    geom_density()+
    geom_vline(xintercept = GSNP_mVAF, linetype = "dashed")+
    scale_color_manual(values=db_pal)+
    scale_y_continuous(expand = c(0,0))+
    theme_syp
  }

#Draw No. of dected bulk tissues
if (F){
  ggplot(vaf_tbl1, aes(x=step_n, y=nPositive))+
    geom_jitter(size=5, alpha=0.2, width=0.05, height=0)+
    xlab("Stage in phylogenetic tree")+ylab("No. of bulk tissues with variant reads")+
    ggtitle(this_db)+
    #coord_cartesian(ylim=c(130,145), xlim=c(0,3))+ #DB6
    #geom_text_repel(data = subset(vaf_tbl1, nPositive > 140 & nPositive < 143 & step_n < 4),aes(label = var_id))+ #DB6
    theme_syp
  }

#Draw correlation between expected VAF and actual VAF
if(F){
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(DP >= 100 & time_of_lineage == 'early')
  ggplot(fm_vaf_tbl, aes(x=expected_VAF, y=VAF))+
    geom_point(aes(color=deadbody), alpha=0.3)+
    geom_abline(slope=1, intercept=0)+
    ggtitle("DP >= 100")+
    scale_color_manual(values = db_pal)+
    theme_syp
}

#Draw variance comparison between dominant layers
if(F){
  ms_tbl <- merged_m_vaf_tbl %>% filter(time_of_lineage == 'early' & DP >= 100 & tissue_class != 'cancer_tissue' & dominant_layer != 'mixed') %>% group_by(deadbody, var_id, dominant_layer) %>% summarise(meanVAF = mean(VAF), sdVAF = sd(VAF))
  my_comparison <- list(c("ectoderm", "endoderm"), c("ectoderm","mesoderm"))
  ggplot(ms_tbl, aes(x=dominant_layer, y=sdVAF))+
    geom_boxplot()+
    stat_compare_means(aes(label = ..p.signif..), comparisons = my_comparison)+
    ylab("SD of VAF")+
    theme_syp
  }

#Draw phylo plot per sample
if(F){
  fm_vaf_tbl <- m_vaf_tbl %>% filter(DP > 30)
  pdf(paste0(targetseq_path_tbl$folder_path[targetseq_path_tbl$deadbody == this_db],'plot_per_sample.pdf'), width=14, height=6)
  for (this_id in unique(fm_vaf_tbl$sample_id)){
    DB_Sample_Plot(fm_vaf_tbl, this_id, this_db)
  }
  dev.off()
}

#Draw phylo plot per tissue_class
fm_vaf_tbl <- m_vaf_tbl %>% filter(DP > 30)
Draw_phylo_tissue(fm_vaf_tbl,lng_dt, this_db)


# comparison between skin_epi and skin_derm
#comparison1: histogram
fm_vaf_tbl <- m_vaf_tbl %>% filter(DP > 30)
ggplot()+
  geom_histogram(data = subset(fm_vaf_tbl, tissue_class2 == 'skin_epi'), aes(x=log10(VAF*100+1)), fill="red", alpha=0.2)+
  geom_histogram(data = subset(fm_vaf_tbl, tissue_class2 == 'skin_derm'), aes(x=log10(VAF*100+1)), fill="blue", alpha=0.2)+
  scale_x_continuous(labels = axis_labels, breaks = log10(axis_labels + 0.1), limits = log10(c(0.5,70)+0.1))+
  xlab("VAF")+
  ggtitle(paste0(this_db," red=skin_epi, blue=skin_derm"))+
  theme_syp
#comparison2: dot and line plot
pmeta_dt$tissue_class %>% unique()
f_pmeta_dt <- pmeta_dt %>% filter(tissue_class == 'skin' & tissue_class2 != 'skin') %>% mutate(sample_id2 = gsub('-epi','',gsub('-derm','', sample_id)))
pair_id_list <- f_pmeta_dt  %>% group_by(sample_id2) %>% count() %>% filter(n==2) %>% .$sample_id2
ff_pmeta_dt <- f_pmeta_dt %>% filter(sample_id2 %in% pair_id_list ) 
library(mgsub)
m_vaf_tbl$sample_id2 <- mgsub(m_vaf_tbl$sample_id, c('-epi','-derm','-Epi','Epi','derm'),c('','','','','')) 
fm_vaf_tbl <- m_vaf_tbl %>% filter(tissue_class2 %in% c('skin_epi', 'skin_derm') & DP > 30 & start_mtn > 5)
fm_vaf_tbl <- fm_vaf_tbl %>% mutate(var_id2 = paste(var_id, sample_id2, sep=':'))
axis_labels = c(0,1,2,3,5,10)
ggplot(fm_vaf_tbl, aes(x=tissue_class2, y=log10(VAF*100+1)))+
  geom_point()+
  geom_line(aes(group=var_id2), alpha=0.2)+
  scale_y_continuous(labels = axis_labels, breaks = log10(axis_labels + 0.1), limits = log10(c(0.5,10)+0.1))+
  ylab("VAF")+
  ggtitle(paste0(this_db," start_mtn > 10 & DP > 30"))+
  theme_syp
#comparison3: phylogenetic tree
fm_vaf_tbl <- m_vaf_tbl %>% filter(DP > 30) ## filter VAF TABLE
Draw_Phylo_Epi_Derm(fm_vaf_tbl, lng_dt, this_db)

#Asymmetry analysis
#L1 vs L2
lid1=paste0(this_db,'_L1'); lid2=paste0(this_db,"_L2")
Draw_maxVAF_1v1(m_vaf_tbl, this_db, lid1,lid2)

#plot per lineage id
if(F){
  pdf(paste0(targetseq_path_tbl$folder_path[targetseq_path_tbl$deadbody == this_db], '/plot_per_lineageid.pdf'), width=14, heigh=10)
  for (lid in unique(m_vaf_tbl$lineage_id2)){
    if (m_vaf_tbl %>% filter(lineage_id2 == lid) %>% .$var_id %>% unique() %>% length() < 2) next
    print(lid)
    Draw_Lineage_Plot(m_vaf_tbl, this_db, lid)
  }
  dev.off()
}

# WGS blood VAF vs mesodermal VAF
this_db = 'DB3'

plot_list=list()
n=0
for(this_db in c('DB3','DB6','DB8','DB9','DB10')){
  n=n+1
  f_vaf_tbl <- merged_m_vaf_tbl %>% filter(DP >= 100 & deadbody == this_db & time_of_lineage == 'early')
  x_order <- f_vaf_tbl %>% group_by(var_id) %>% summarise(medVAF = median(VAF)) %>% arrange(desc(medVAF)) %>% pull(var_id)
  f_vaf_tbl$sample_id %>% unique()
  if(this_db == 'DB3'){
    plot_list[[n]] <- ggplot(data = subset(f_vaf_tbl, dominant_layer2 == 'mesoderm'), aes(x=var_id, y=log10(VAF+0.01)))+
      geom_boxplot(outlier.size=-1)+
      geom_point(data=subset(f_vaf_tbl, sample_id == '3_WholeBlood-3'),color='blue')+
      geom_point(aes(y=log10(WGS_bloodVAF+0.01)), color="red")+
      scale_x_discrete(limits= x_order)+
      theme(axis.text.x = element_text(angle=45, hjust=1))+
      ggtitle(this_db)
  }else {
    plot_list[[n]] <- ggplot(data = subset(f_vaf_tbl, dominant_layer2 == 'mesoderm'), aes(x=var_id, y=log10(VAF+0.01)))+
      geom_boxplot(outlier.size=-1)+
      geom_point(aes(y=log10(WGS_bloodVAF+0.01)), color="red")+
      scale_x_discrete(limits= x_order)+
      theme(axis.text.x = element_text(angle=45, hjust=1))+
      ggtitle(this_db)
  }
}
plot_grid(plotlist = plot_list, nrow=2)



#plot VAFratio drop and assigned mtn of each variant
if(F){
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(DP > 100 & time_of_lineage == 'early' & dominant_layer %in% c('ectoderm','mesoderm') & deadbody !='DB3')
  fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group = ifelse(dominant_layer == 'mesoderm', 'mesoderm', paste0(tissue_type, '_epidermis'))) %>% filter(tissue_type != 'formalin_fixed')
  med_tbl <- fm_vaf_tbl %>% group_by(deadbody,group, assigned_mtn) %>% summarise(medVAF = median(VAF)) %>% spread(key=group, value=medVAF) %>% 
    filter(mesoderm >0) %>% mutate(ratio = fresh_frozen_epidermis/mesoderm)
  med_tbl %>% filter(ratio > 10)
  ggplot(med_tbl, aes(x=assigned_mtn, y=ratio))+
    geom_point()+
    geom_smooth()+
    coord_cartesian(ylim=c(0,1.2))+
    facet_wrap(~deadbody)
}

#plot differentiation timing
this_db='DB3'
fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(DP > 100 & tissue_class != 'cancer_tissue' & time_of_lineage=='early' & deadbody == this_db)


cutoff1=5;cutoff2=3
Draw_phylo_timing_TotalAlign_3group(fm_vaf_tbl,lng_dt, this_db, cutoff1, cutoff2)
diff_tbl_3group <- fm_vaf_tbl  %>% group_by(lineage_id2, var_id) %>% summarise(medVAF = median(VAF)) %>% mutate(diff_group = ifelse(medVAF >= cutoff1*0.01, 'early',
                                                                                                                                                                            ifelse(medVAF < cutoff2*0.01, 'late', 'mid')))
Draw_Heatmap_Diff_3group(fm_vaf_tbl, diff_tbl_3group)
Draw_tSNE_Diff_3group(fm_vaf_tbl, diff_tbl_3group)
target_class <- fm_vaf_tbl %>% select(anatomy_class3, sample_id) %>% unique() %>% group_by(anatomy_class3) %>% count() %>% filter(n >= 3) %>% .$anatomy_class3
Draw_Silhouette_Coeff_Anatomy_3group(subset(fm_vaf_tbl, anatomy_class3 %in% target_class), diff_tbl_3group)

tmp1 <- left_join(diff_tbl_3group,merged_m_vaf_tbl %>% select(var_id, start_mtn, rank, n_vars_same_rank) %>% unique())
tmp1 <- tmp1 %>% mutate(estimated_mtn = start_mtn + rank -1 + n_vars_same_rank/2)
ggplot(tmp1, aes(x=var_id, y=estimated_mtn, color=diff_group))+
  geom_point()

tmp1 %>% filter(is.na(estimated_mtn)== T)
mt1 <- tmp1$estimated_mtn[tmp1$diff_group == 'mid' & is.na(tmp1$estimated_mtn) == F]
mt2 <- tmp1$estimated_mtn[tmp1$diff_group == 'late' & is.na(tmp1$estimated_mtn) == F]
summary(mt1)
summary(mt2)
rate_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Mutation_rate_simulation/df_forest_plot.txt')
cut1 <- quantile(mt1,probs=0.25) #6
cut2 <- quantile(mt1,probs=0.75) #7.5
DB6_M0 <- rate_dt %>% filter(ID == 'DB6' & type == 'M0') %>% pull(mean)
DB6_M1 <- rate_dt %>% filter(ID == 'DB6' & type == 'M1') %>% pull(mean)

2+ (cut1 - 2*DB6_M0)/DB6_M1  #4.94
2+ (cut2 - 2*DB6_M0)/DB6_M1  #6.30


Draw_tSNE_Diff_3group_edit(fm_vaf_tbl, diff_tbl_3group)
#Draw_phylo_timing_skin(fm_vaf_tbl,lng_dt, this_db, 0.5)
Draw_phylo_timing_TotalAlign_2group(fm_vaf_tbl, lng_dt, this_db, cutoff1)
diff_tbl_2group <- fm_vaf_tbl %>% group_by(lineage_id2, var_id) %>% summarise(medVAF = median(VAF)) %>% mutate(diff_group = ifelse(medVAF >= cutoff1*0.01, 'early','late'))
Draw_Heatmap_Diff_2group(fm_vaf_tbl, diff_tbl_2group)
#Draw_PCA_Diff(fm_vaf_tbl, diff_tbl)
Draw_tSNE_Diff_2group(fm_vaf_tbl, diff_tbl_2group)
#Draw_tSNE_Diff(subset(fm_vaf_tbl, tissue_class %in% c('muscle','skin')), diff_tbl)
#Draw_Silhouette_Coeff_Anatomy_2group(subset(fm_vaf_tbl, anatomy_class3 %in% target_class), diff_tbl_2group)






####=====================================================Draw DB merged plot======================================
#merged_fm_vaf_tbl <- fm_vaf_tbl
#merged_fm_vaf_tbl <- bind_rows(merged_fm_vaf_tbl, fm_vaf_tbl)
merged_fm_vaf_tbl$deadbody %>% unique()
merged_fm_vaf_tbl$tissue_class2 %>% unique()
merged_fm_vaf_tbl <- merged_fm_vaf_tbl %>% mutate(clone_type = ifelse(tissue_class2 =='skin_epi', 'epidermis', 
                                                                      ifelse(deadbody == 'DB10' & (tissue_class2 == 'cancer'| sample_id == '10_LowerRtLung'), 'cancer', 'other')))

ggplot(merged_fm_vaf_tbl, aes(x=expected_VAF, y=VAF))+
  geom_point(aes(color=deadbody, shape=clone_type), alpha=0.7, size=3)+
  geom_abline(slope=1, intercept=0)+
  scale_color_manual(values = db_pal)+
  scale_shape_manual(values = c('epidermis'=4, 'cancer' = 5, 'other'=20))+
  theme_syp
