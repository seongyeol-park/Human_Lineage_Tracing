stop("a")

#load libraries
library(tidyverse)
library(ggtree)
library(ggsci)
library(scales)
library(cowplot)
library(ggrepel)
library(ggpubr)
library(grid)

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

targetseq_path_tbl <- tribble(
  ~deadbody, ~folder_path, ~Call_name, 
  "DB3",'/home/users/sypark/00_Project/06_LineageTracing/db3/04-5_Target_seq/05_re_processing/','DB3_TargetSeq_Target_list.GSNPadded.60sCall',
  "DB6",'/home/users/sypark/00_Project/06_LineageTracing/db6/14_TargetSeq/','DB6_TargetSeq_Target_list.merged.GSNPadded.143sCall',
  "DB8",'/home/users/sypark/00_Project/06_LineageTracing/db8/11_TargetSeq/','DB8_TargetSeq_Target_list.merged.GSNPadded.44sCall',
  "DB9",'/home/users/sypark/00_Project/06_LineageTracing/db9/11_TargetSeq/','DB9_TargetSeq_Target_list.merged.GSNPadded.43sCall',
  "DB10",'/home/users/sypark/00_Project/06_LineageTracing/db10/11_TargetSeq/','DB10_TargetSeq_Target_list.merged.GSNPadded.89sCall'
)


#load data
#pmeta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_TargetSeq_210325.txt')
pmeta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_TargetSeq_210330.txt')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_210322.txt')
m_rate_dt <- read_tsv("/home/users/sypark/00_Project/06_LineageTracing/meta_data/estimated_mutation_rate_pcpcd.tsv")
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_210330.txt') %>% filter(current_final_for_lineage == 'Y')
meta_dt <- meta_dt %>% mutate(Source3 = ifelse(is.na(Source_side)==T, Source2, paste0(Source_side,'_',Source2)))
meta_dt$deadbody <- factor(meta_dt$deadbody, levels =c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))

#Sourcing color palette
source('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Lineage_Tracing_Color_Palette.R')


#theme setting
theme_syp = theme_bw() + theme(axis.text = element_text(size=15), axis.title = element_text(size=18), panel.border = element_blank(), axis.line = element_line(), title = element_text(size=18),
                               legend.text = element_text(size=15))

#Sourcing functions
source('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/TargetSeq_functions_afterTopUp.R')


#process of making f_vaf_tbl (you don't need to do this again)
if(F){
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
      col_list <- merged_tbl %>% select(starts_with(paste0(gsub('DB','',this_db),'_'))) %>% colnames()
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
    saveRDS(merged_m_vaf_tbl, '/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/merged_m_vaf_tbl.rds')
  }
  
  if(F){
    merged_m_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/merged_m_vaf_tbl.rds')
    #DP per early mutations
    tmp_dt <- merged_m_vaf_tbl %>% filter(time_of_lineage == 'early')
    x_order <- tmp_dt %>% group_by(var_id) %>% summarise(medDP=median(DP)) %>% arrange(desc(medDP)) %>% pull(var_id)
    ggplot(tmp_dt, aes(x=var_id, y=DP, fill=deadbody))+
      geom_boxplot()+
      geom_hline(yintercept=1000, color='red')+
      geom_hline(yintercept=500, color='blue')+
      geom_hline(yintercept=100, color='green')+
      scale_x_discrete(limits=x_order)+
      scale_fill_manual(values = db_pal)+
      theme_syp+
      theme(axis.text.x=element_blank())
    
    #DPper samples
    tmp_dt <- merged_m_vaf_tbl %>% filter(time_of_lineage == 'early')
    x_order <- tmp_dt %>% group_by(sample_id) %>% summarise(medDP=median(DP)) %>% arrange(desc(medDP)) %>% pull(sample_id)
    ggplot(tmp_dt, aes(x=sample_id, y=DP, fill=deadbody))+
      geom_boxplot()+
      scale_x_discrete(limits=x_order)+
      geom_hline(yintercept=1000, color='red')+
      geom_hline(yintercept=500, color='blue')+
      geom_hline(yintercept=100, color='green')+
      scale_fill_manual(values = db_pal)+
      theme_syp+
      theme(axis.text.x=element_blank())
    
    # define low targeted early mutations
    merged_m_vaf_tbl %>% filter(time_of_lineage == 'early') %>% pull(var_id) %>% unique() %>% length() #441
    nonDB3_low_target_ids <- merged_m_vaf_tbl %>% filter(deadbody != 'DB3' & time_of_lineage == 'early') %>% group_by(deadbody, var_id) %>% summarise(medDP=median(DP)) %>%
      filter(medDP < 500) %>% pull(var_id)
    length(nonDB3_low_target_ids)  #32
    DB3_low_target_ids <- merged_m_vaf_tbl %>% filter(deadbody == 'DB3' & time_of_lineage == 'early') %>% group_by(deadbody, var_id) %>% summarise(medDP=median(DP)) %>%
      filter(medDP < 100) %>% pull(var_id)
    length(DB3_low_target_ids)  #5
    # exclude low targeted early mutations
    f_vaf_tbl <- merged_m_vaf_tbl %>% filter(!var_id %in% c(nonDB3_low_target_ids, DB3_low_target_ids))
    f_vaf_tbl %>% filter(time_of_lineage == 'early') %>% pull(var_id) %>% unique() %>% length() #404
    saveRDS(f_vaf_tbl, '/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/f_vaf_tbl.rds')
  }
  if(F){
    f_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/f_vaf_tbl.rds')
    
    if(F){
      #VAF distribution
      ff_vaf_tbl <- f_vaf_tbl %>% filter(time_of_lineage=='early')
      g1 <- ggplot(ff_vaf_tbl, aes(x=VAF))+
        geom_histogram()
      g2 <- ggplot(ff_vaf_tbl, aes(x=log10(VAF+0.00001)))+
        geom_histogram()
      plot_grid(g1, g2, nrow=1, align="h", axis="tb")
    }
    
    #add log10VAF value
    f_vaf_tbl$logVAF = log10(f_vaf_tbl$VAF+0.00001)
    #add lineage_id
    f_vaf_tbl <- f_vaf_tbl %>% separate(lineage_id2, c('deadbody','lineage_id'), remove = F, sep='_') 
    
    #write list of varid and n of tissues, preparation for writing VAF drop
    if(F){
      tmp_tbl1 <- f_vaf_tbl %>%filter(time_of_lineage=='early' & tissue_class != 'cancer_tissue' & dominant_layer2 !='ectoderm' & dominant_layer2 != 'ecto_mesoderm') %>% group_by(deadbody, var_id, lineage_id2) %>% count() %>%
        rename(n_sample_noCancer_noEcto = n) 
      tmp_tbl2 <- f_vaf_tbl %>%filter(time_of_lineage=='early' & tissue_class != 'cancer_tissue' & dominant_layer2 !='ectoderm' & dominant_layer2 != 'ecto_mesoderm') %>% group_by(deadbody, var_id, lineage_id2)  %>%
        summarise(medDP = median(DP))
      tmp_tbl3 <- f_vaf_tbl %>%filter(time_of_lineage=='early' & tissue_class != 'cancer_tissue' & dominant_layer2 !='ectoderm' & dominant_layer2 != 'ecto_mesoderm') %>% group_by(deadbody, var_id, lineage_id2)  %>%
        summarise(meanLogVAF = mean(logVAF))
      left_join(left_join(tmp_tbl1, tmp_tbl2) , tmp_tbl3) %>% arrange(lineage_id2, desc(meanLogVAF)) %>%
        write_tsv('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/table_for_rank.tsv')
    }
    
    #Intralineage VAF paired t test (cancer and ecto excluded)
    if(F){ 
      lids <-f_vaf_tbl %>%  filter(time_of_lineage=='early' & tissue_class != 'cancer_tissue' & dominant_layer2 !='ectoderm' & dominant_layer2 != 'ecto_mesoderm') %>% select(lineage_id2, var_id) %>% unique() %>% group_by(lineage_id2) %>% count() %>% filter(n>1) %>% pull(lineage_id2)
      pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/pairedT_logVAF.pdf', width=15, height=20)
      for(lid in lids){
        print(lid)
        f_tbl <- f_vaf_tbl %>%  filter(lineage_id2 == lid & tissue_class != 'cancer_tissue' & dominant_layer2 !='ectoderm' & dominant_layer2 != 'ecto_mesoderm') %>% select(var_id, logVAF, DP, sample_id)
        expand.grid(unique(f_tbl$var_id),unique(f_tbl$var_id))
        dp_tbl <- f_tbl %>% group_by(var_id) %>% summarise(medDP = median(DP))
        lowest_v <- f_tbl$logVAF %>% min()
        library(ggpubr)
        x_order <- f_tbl %>% group_by(var_id) %>% summarise(meanLogVAF = mean(logVAF)) %>% arrange(desc(meanLogVAF)) %>% pull(var_id)
        mycomparison <- combn(unique(f_tbl$var_id),2, simplify = F)
        g1 <- ggplot(f_tbl, aes(x=var_id, y=logVAF))+
          geom_boxplot()+
          geom_text(data=dp_tbl, aes(y=lowest_v, label=medDP), color="blue")+
          geom_line(aes(group=sample_id), alpha=0.1)+
          stat_compare_means(aes(label = ..p.signif..), comparisons = mycomparison,paired=T, method='t.test')+
          scale_x_discrete(limits = x_order)+
          theme_syp+theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
          ggtitle(lid)
        g2 <- ggplot(f_tbl, aes(x=var_id, y=logVAF))+
          geom_boxplot(outlier.size=-1)+
          geom_jitter(height=0, width=0.1, alpha=0.2, aes(size=DP))+
          geom_line(aes(group=sample_id), alpha=0.1)+
          scale_x_discrete(limits = x_order)+
          theme_syp+theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=45, hjust=1))
        print(plot_grid(g1,g2,ncol=1, align="v", axis="lr"))
      }
      dev.off()
    }
    
    if(T){
      #Using the above plots, assigning rank was performed manually!!!
      #assign rank 
      rank_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/table_for_rank_assigned.tsv')
      #check data integrity
      rank_dt %>% group_by(lineage_id2, rank) %>% dplyr::count() %>% arrange(lineage_id2, rank) %>% group_by(lineage_id2) %>%
        mutate(cumsum=cumsum(n)) %>% mutate(expected = lag(cumsum)+1) %>% mutate(test = ifelse(rank == expected, 'correct','error')) %>%
        filter(test == 'error')  # 0 row, no error
      rank_dt <- left_join(rank_dt, lng_dt %>% select(lineage_id2, start_mtn, end_mtn, n_pointmt))
      rank_dt <- left_join(rank_dt, rank_dt %>% group_by(lineage_id2) %>% summarise(targeted_n = n()))
      rank_dt <- left_join(rank_dt, rank_dt %>% group_by(lineage_id2, rank) %>% summarise(n_vars_same_rank = n()))
      rank_dt <- rank_dt %>% group_by(lineage_id2) %>% mutate(mean_rank = rank(rank), rank_ctie_first = rank(rank, ties.method = "first")) %>%
        mutate(adjusted_mean_rank = mean_rank * n_pointmt/targeted_n) %>% mutate(assigned_mean_mtn = start_mtn + adjusted_mean_rank) %>%
        mutate(adjusted_ctie_rank = rank_ctie_first * n_pointmt/targeted_n) %>% mutate(assigned_ctie_mtn = start_mtn + adjusted_ctie_rank)
      f_vaf_tbl <- left_join(f_vaf_tbl, rank_dt %>% 
                               select(var_id, lineage_id2, rank, n_vars_same_rank, mean_rank, rank_ctie_first, adjusted_mean_rank, assigned_mean_mtn, adjusted_ctie_rank, assigned_ctie_mtn))
    }
    saveRDS(f_vaf_tbl, '/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/f_vaf_tbl_rank.rds')  
  }
}

#load data!!!!!!!!!!!!!!
f_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/f_vaf_tbl_rank.rds')



# VAF of L1+L2 (first figure 1e)
GSNP_mVAF <- f_vaf_tbl %>% filter(grepl('GSNP', lineage_id2)) %>% pull(VAF) %>% mean(., na.rm=T)
ff_vaf_tbl <- f_vaf_tbl %>% filter(step_n==1 & rank ==1) 
ff_vaf_tbl <- ff_vaf_tbl %>% separate(lineage_id2, c('deadbody','lineage_id'), sep='_')
tmp_dt1 <- ff_vaf_tbl %>% filter(lineage_id == 'L1') %>% group_by(deadbody, sample_id) %>% summarise(meanL1 = mean(VAF))
tmp_dt2 <- ff_vaf_tbl %>% filter(lineage_id == 'L2') %>% group_by(deadbody, sample_id) %>% summarise(meanL2 = mean(VAF))
tmp_dt <- left_join(tmp_dt1, tmp_dt2) %>% mutate(L1plusL2 = meanL1 + meanL2)
tmp_dt$deadbody <- factor(tmp_dt$deadbody, levels = c('DB3','DB6','DB8','DB9','DB10'))
pdf("/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/merged_L1_L2_VAF.pdf", width=6, height=4)
ggplot(tmp_dt, aes(x=deadbody, y=L1plusL2*100, color=deadbody))+
  geom_jitter(height=0, width=0.2, alpha=0.5)+
  geom_boxplot(outlier.size=-1, color="red", alpha=0)+
  scale_y_continuous(limits=c(0,100))+
  scale_color_manual(values=db_pal)+
  theme_syp+
  ylab("VAF in bulk tissues(%)")+
  theme(axis.title.x = element_blank())
dev.off()
tmp_dt %>% group_by(deadbody) %>% count()


# VAF of L1 and L2 with expected VAF (first figure 2a)
tmp_dt1 <- f_vaf_tbl %>% filter(lineage_id =='L1') #when drawing L1 mutations
tmp_dt1 <- f_vaf_tbl %>% filter(lineage_id =='L2') #when drawing L2 mutations
p_list=list()
n=0
for (this_db in c("DB3","DB6",'DB8','DB9','DB10')){
  n=n+1
  print(this_db)
  tmp_dt2 <- tmp_dt1 %>% filter(deadbody == this_db)
  tmp_dt2 <- tmp_dt2 %>% group_by(deadbody, expected_VAF, tissue_class2, dominant_layer2, sample_id) %>% summarise(meanVAF = mean(VAF))
  x_order <- tmp_dt2 %>% arrange(desc(meanVAF)) %>% pull(sample_id)
  exp_v <- tmp_dt2$expected_VAF %>% unique()
  p_list[[n]] <- ggplot(tmp_dt2, aes(x=sample_id, y=meanVAF, fill=dominant_layer2))+
    geom_bar(stat="identity")+
    geom_hline(yintercept=exp_v, color='red')+
    scale_x_discrete(limits=x_order)+
    scale_y_continuous(limits=c(0,1))+
    #scale_fill_manual(values = db_pal)+
    ggtitle(this_db)+
    theme_syp+theme(legend.position='right', axis.text.x =element_blank())
}
plot_grid(plotlist=p_list)


#Expected VAF, Target-seq VAF correlation (first figure 2b)
ggplot(f_vaf_tbl, aes(x=expected_VAF, y=VAF))+
  geom_point(aes(color=deadbody), alpha=0.5)+
  geom_abline(slope=1, intercept=0)+
  scale_color_manual(values = db_pal)+
  theme_syp
  


##############germ layer divergence######## 
if(T){
  this_db = 'DB6'
  plot_list=list()
  for (this_db in c('DB3','DB8','DB9','DB10')){
    print(this_db)
    fm_vaf_tbl <- f_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
    fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group2 = dominant_layer2) %>% filter(group2 %in% c('ectoderm','endoderm','mesoderm'))
    ct_tbl <- fm_vaf_tbl %>% group_by(deadbody,group2, var_id) %>% count() %>% group_by(group2) %>% summarise(mean_n = mean(n))
    incl_ids <- ct_tbl %>% filter(mean_n >=2) %>% pull(group2)
    fm_vaf_tbl <- fm_vaf_tbl %>% filter(group2 %in% incl_ids)
    med_tbl <- fm_vaf_tbl %>% group_by(deadbody,group2, var_id) %>% summarise(medVAF = median(VAF), q0 = quantile(VAF)[2], q1=quantile(VAF)[4])
    med_tbl <- left_join(med_tbl, fm_vaf_tbl %>% select(deadbody, var_id, lineage_id2, step_n) %>% unique())
    med_tbl <- med_tbl %>% ungroup() %>%  mutate(first_group = apply(.["lineage_id2"], 1,function(x) unlist(strsplit(x, '-'))[1]))
    var_order <- fm_vaf_tbl %>% group_by(deadbody,var_id) %>% summarise(medVAF = median(VAF)) %>% arrange(desc(medVAF)) %>% pull(var_id)
    #graph with only 3layer ###main figure3a
    f_med_tbl <- med_tbl %>% filter(group2 %in% c('ectoderm','mesoderm','endoderm','meso_endoderm')) #select some data
    plot_list[[this_db]] <- ggplot(f_med_tbl, aes(x=var_id, y=log10(medVAF*100+0.01), color = group2))+  
      geom_point(size=3, alpha=0.3)+
      #geom_errorbar(aes(ymin=log10(q0*100+0.01), ymax=log10(q1*100+0.01)), alpha=0.5)+
      stat_smooth(aes(group = group2))+
      coord_cartesian(ylim=log10(c(0,50)+0.01))+
      scale_color_manual(values=total_pal, name='')+
      ggtitle(this_db)+
      ylab("Median VAF (%)")+xlab("Early mutations")+
      scale_x_discrete(limits = var_order)+
      scale_y_continuous(expand = c(0,0), labels = c(0.01, 0.05, 0.1, 0.5,1,2,3,4,5,10,20,30), breaks= log10(c(0.01, 0.05, 0.1, 0.5,1,2,3,4,5,10,20,30)+0.01))+
      theme_syp+theme(axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks.x = element_blank(), legend.position = 'none')
  }
  pdf("/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/DB38910_germlayer_divergence.pdf", height=4, width=16)
  plot_grid(plotlist = plot_list,nrow=1, align="h")
  dev.off()
  
  pdf("/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/DB6_germlayer_divergence.pdf", height=6, width=8)
  plot_list["DB6"]
  dev.off()
  
  plot_grid(plotlist=list(plot_list[c("DB3","DB8","DB9","DB10")]), nrow=2, align="hv")
  
  
  #more extensive figure in DB6 (Fig.SD3)
  this_db = 'DB6'
  fm_vaf_tbl <- f_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
  fm_vaf_tbl <- fm_vaf_tbl %>% mutate(dominant_layer3 = ifelse(tissue_class2 %in% c('lung','liver','colon'), tissue_class2, dominant_layer2))
  fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group2 = dominant_layer3)
  f_med_tbl <- med_tbl %>% filter(group2 %in% c('ectoderm','mesoderm','liver','lung','colon'))
  my_pal <- c('ectoderm' = '#EBB50A',
              'mesoderm' = '#245DB5',
              "colon" = "#d10a67",
              "liver" = "#047004",
              "lung"= "#670b8e"
  )
  ggplot(f_med_tbl, aes(x=var_id, y=log10(medVAF+0.0001), color = group2))+  
    geom_point(size=3, alpha=0.3)+
    #geom_errorbar(aes(ymin=log10(q0+0.0001), ymax=log10(q1+0.0001)), alpha=0.5)+
    stat_smooth(aes(group = group2))+
    coord_cartesian(ylim=log10(c(0,0.5)+0.0001))+
    scale_color_manual(values=my_pal, name='')+
    #ggtitle(this_db)+
    ylab("Median VAF")+xlab("Early mutations")+
    scale_x_discrete(limits = var_order)+
    scale_y_continuous(expand = c(0,0.05), labels = c(0.01, 0.05, 0.1, 0.5,1,2,3,4,5,10,20,30)*0.01, breaks= log10(c(0.01, 0.05, 0.1, 0.5,1,2,3,4,5,10,20,30)*0.01+0.0001))+
    theme_syp+theme(axis.text.x = element_blank(), panel.grid = element_blank())
  
  #graph with all samples
  ggplot(med_tbl, aes(x=var_id, y=log10(medVAF+0.0001), color = group2))+
    geom_point(size=3, alpha=0.5)+
    geom_errorbar(aes(ymin=log10(q0+0.0001), ymax=log10(q1+0.0001)), alpha=0.5)+
    stat_smooth(aes(group = group2))+
    coord_cartesian(ylim=log10(c(0,0.5)+0.0001))+
    scale_color_manual(values=dl2_pal)+
    ggtitle(this_db)+
    scale_x_discrete(limits = var_order)+
    scale_y_continuous(expand = c(0,0.05), labels = c(0.01, 0.05, 0.1, 0.5,1,2,3,4,5,10,20,30)*0.01, breaks= log10(c(0.01, 0.05, 0.1, 0.5,1,2,3,4,5,10,20,30)*0.01+0.0001))+
    theme_syp+theme(axis.text.x = element_blank(), panel.grid = element_blank())
 
  
  #graph with all datapoint in 3layer
  tmp_tbl <- fm_vaf_tbl %>% filter(dominant_layer2 %in% c('ectoderm','mesoderm','endoderm'))
  ggplot(tmp_tbl, aes(x=var_id, y=log10(VAF+0.0001), color=dominant_layer2))+
    geom_point(alpha=0.3)+
    scale_x_discrete(limits=var_order)+
    scale_y_continuous(expand = c(0,0.05), labels = c(0.01, 0.05, 0.1, 0.5,1,2,3,4,5,10,20,30)*0.01, breaks= log10(c(0.01, 0.05,  0.1, 0.5,1,2,3,4,5,10,20,30)*0.01+0.0001))+
    scale_color_manual(values=dl2_pal)+
    theme_syp+theme(axis.text.x = element_blank(), panel.grid = element_blank())+
    ggtitle(this_db)
  

  
  # graph using mean VAF
  mean_tbl <- fm_vaf_tbl %>% group_by(deadbody,group2, var_id) %>% summarise(meanVAF = mean(VAF), sd = sd(VAF))
  mean_tbl <- mean_tbl %>% mutate(q1= meanVAF + sd, q0 = ifelse(meanVAF -2*sd <0, 0, meanVAF -sd))
  mean_tbl <- left_join(mean_tbl, fm_vaf_tbl %>% select(deadbody, var_id, lineage_id2, step_n) %>% unique())
  mean_tbl <- mean_tbl %>% ungroup() %>%  mutate(first_group = apply(.["lineage_id2"], 1,function(x) unlist(strsplit(x, '-'))[1]))
  var_order <- fm_vaf_tbl %>% group_by(deadbody,var_id) %>% summarise(meanVAF = mean(VAF)) %>% arrange(desc(meanVAF)) %>% pull(var_id)
  ggplot(mean_tbl, aes(x=var_id, y=log10(meanVAF+0.0001), color = group2))+
    geom_point(size=3, alpha=0.5)+
    geom_errorbar(aes(ymin=log10(q0+0.0001), ymax=log10(q1+0.0001)), alpha=0.5)+
    stat_smooth(aes(group = group2))+
    coord_cartesian(ylim=log10(c(0,0.5)+0.0001))+
    scale_color_manual(values=dl2_pal)+
    ggtitle(this_db)+
    scale_x_discrete(limits = var_order)+
    scale_y_continuous(expand = c(0,0.05), labels = c(0.01, 0.05, 0.1, 0.5,1,2,3,4,5,10,20,30)*0.01, breaks= log10(c(0.01, 0.05, 0.1, 0.5,1,2,3,4,5,10,20,30)*0.01+0.0001))+
    theme_syp+theme(axis.text.x = element_blank(), panel.grid = element_blank())
}

#get the timing of  Ecto mesodermal divergence (deprecated)
f_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/f_vaf_tbl_rank.rds')
this_db='DB6'
p10_tbl <- tibble(deadbody = as.character(), min = as.numeric(), min_step = as.numeric(), p5 = as.numeric(), p10= as.numeric(), p25= as.numeric())
m=0
for (this_db in c('DB6','DB8','DB9','DB10')){
  m=m+1
  ff_vaf_tbl <- f_vaf_tbl %>% filter(time_of_lineage == 'early' & tissue_class !='cancer_tissue' & deadbody == this_db)
  ff_vaf_tbl <- ff_vaf_tbl %>% mutate(var_g_id = paste0(lineage_id, '_', rank))
  this_id = ff_vaf_tbl$var_g_id[1]
  res_tbl <- tibble(deadbody = as.character(), var_g_id = as.character(), mean_ecto_vaf = as.numeric(), mean_mesoendo_vaf = as.numeric(), wilcox_p = as.numeric())
  n=0
  for (this_id in ff_vaf_tbl$var_g_id %>% unique()){
    n=n+1
    print(this_id)
    ecto_vaf <- ff_vaf_tbl %>% filter(var_g_id == this_id) %>% filter(dominant_layer2 == 'ectoderm') %>% pull(VAF)
    mesoendo_vaf <- ff_vaf_tbl %>% filter(var_g_id == this_id) %>% filter(dominant_layer2 %in% c('mesoderm','endoderm')) %>% pull(VAF)
    wilcox_p <- wilcox.test(ecto_vaf, mesoendo_vaf)$p.value
    res_tbl[n,] <- list(this_db, this_id, mean(ecto_vaf), mean(mesoendo_vaf), wilcox_p)
  }
  ff_vaf_tbl <- left_join(ff_vaf_tbl, res_tbl)
  tmp <- ff_vaf_tbl %>% filter(mean_ecto_vaf < mean_mesoendo_vaf & wilcox_p < 0.05) %>% select(var_id, assigned_ctie_mtn, step_n) %>% unique()
  min <- tmp$assigned_ctie_mtn %>% min()
  min_step <- tmp$step_n[tmp$assigned_ctie_mtn == min]
  p5 <- tmp$assigned_ctie_mtn %>% quantile(0.05)
  p10 <- tmp$assigned_ctie_mtn %>% quantile(0.1)
  p25 <- tmp$assigned_ctie_mtn %>% quantile(0.25)
  g <- ggplot(tmp, aes(x=assigned_ctie_mtn))+
    geom_density()+
    geom_vline(xintercept=p10, color="red")+
    geom_vline(xintercept=p25, color="blue")+
    ggtitle(paste0(this_db, ' 10%=', p10, ' 25%=', p25))
  #print(g)
  p10_tbl[m,] <- list(this_db, min, min_step, p5, p10, p25)
}
p_tbl <- left_join(p10_tbl, m_rate_dt)
p_tbl <- p_tbl %>% mutate(min_stage = ifelse(min < M0, min/M0 +1, (min-M0)/M1+2),
                          p25_stage = ifelse(p25 < M0, p25/M0+1, (p25 - M0)/M1+2),
                          p10_stage = ifelse(p10 < M0, p10/M0+1, (p10 - M0)/M1+2),
                          p5_stage = ifelse(p5 < M0, p5/M0+1, (p5 - M0)/M1+2))
p_tbl <- p_tbl %>% mutate(final_min_stage = ifelse(min_step > min_stage, min_step, floor(min_stage)))
p_tbl



#get the timing of  Ecto mesodermal divergence2
f_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/f_vaf_tbl_rank.rds')
this_db='DB10'
for (this_db in c('DB6','DB8','DB9','DB10')){
  ff_vaf_tbl <- f_vaf_tbl %>% filter(time_of_lineage == 'early' & tissue_class !='cancer_tissue' & deadbody == this_db)
  ff_vaf_tbl <- ff_vaf_tbl %>% mutate(var_g_id = paste0(lineage_id, '_', rank))
  this_id = ff_vaf_tbl$var_g_id[1]
  res_tbl <- tibble(deadbody = as.character(), var_g_id = as.character(), mean_ecto_vaf = as.numeric(), mean_mesoendo_vaf = as.numeric(), wilcox_p = as.numeric())
  n=0
  for (this_id in ff_vaf_tbl$var_g_id %>% unique()){
    n=n+1
    print(this_id)
    ecto_vaf <- ff_vaf_tbl %>% filter(var_g_id == this_id) %>% filter(dominant_layer2 == 'ectoderm') %>% pull(VAF)
    mesoendo_vaf <- ff_vaf_tbl %>% filter(var_g_id == this_id) %>% filter(dominant_layer2 %in% c('mesoderm','endoderm')) %>% pull(VAF)
    wilcox_p <- wilcox.test(ecto_vaf, mesoendo_vaf)$p.value
    res_tbl[n,] <- list(this_db, this_id, mean(ecto_vaf), mean(mesoendo_vaf), wilcox_p)
  }
  ff_vaf_tbl <- left_join(ff_vaf_tbl, res_tbl)
  tmp <- ff_vaf_tbl %>% select(lineage_id, var_id, rank, start_mtn, rank_ctie_first, mean_ecto_vaf, mean_mesoendo_vaf, wilcox_p) %>% unique()
  tmp <- tmp %>% mutate(assigned_location = start_mtn + rank_ctie_first -1)
  pval_cutoff=0.1
  tmp <- tmp %>% mutate(p_group = ifelse(wilcox_p > pval_cutoff, 'non_sig',ifelse(mean_ecto_vaf > mean_mesoendo_vaf, 'ecto_high','mesoendo_high')))
  g1 <- ggplot(tmp)+
    geom_bar(aes(x=assigned_location, fill=p_group), stat="count")+
    ggtitle(this_db)
  tmp2 <- tmp %>% group_by(assigned_location) %>% summarise(n_nonsig = sum(wilcox_p > pval_cutoff), n_ecto = sum(wilcox_p <= pval_cutoff & mean_ecto_vaf > mean_mesoendo_vaf), n_mesoendo = sum(wilcox_p <=pval_cutoff & mean_ecto_vaf < mean_mesoendo_vaf))
  g2 <- ggplot(tmp2, aes(x=assigned_location, y=n_mesoendo/(n_ecto+n_mesoendo)))+
    geom_point()+
    geom_line(aes(group=1))+
    scale_y_continuous(limits=c(0,1))+
    scale_x_continuous(limits=c(0,35), breaks=c(1:10, seq(15,35,3)))+
    ggtitle(this_db)
  g3 <- ggplot(tmp2, aes(x=assigned_location, y=n_mesoendo/(n_ecto+n_mesoendo+n_nonsig)))+
    geom_point()+
    geom_line(aes(group=1))+
    scale_y_continuous(limits=c(0,1))+
    scale_x_continuous(limits=c(0,35), breaks=c(1:10, seq(15,35,3)))+
    geom_hline(yintercept=c(0.3,0.5), color="red")+
    ggtitle(this_db)
  
  plot_grid(g1, g2,g3, nrow=1) %>% print()
}

res <- tribble(~deadbody, ~mt_m50,
               'DB6',6,
               'DB8',6,
               'DB9',12,
               'DB10',10)
res <- left_join(res, m_rate_dt)
res <- res %>% mutate(stage = ifelse(mt_m50 < M0, mt_m50/M0 +1, (mt_m50-M0)/M1+2) %>% floor())


#########DB6 liver only and lung only comparison
this_db = 'DB6'
fm_vaf_tbl <- f_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
fm_vaf_tbl <- fm_vaf_tbl %>% filter(tissue_class2 != 'lung')
fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group2 = dominant_layer2)
#fm_vaf_tbl <- fm_vaf_tbl %>% filter(group2 %in% c('mesoderm','endoderm','mixed','fresh_frozen_epidermis'))
ct_tbl <- fm_vaf_tbl %>% group_by(deadbody,group2, var_id) %>% count() %>% group_by(group2) %>% summarise(mean_n = mean(n))
print(ct_tbl)
incl_ids <- ct_tbl %>% filter(mean_n >=2) %>% pull(group2)
print(incl_ids)
fm_vaf_tbl <- fm_vaf_tbl %>% filter(group2 %in% incl_ids)

med_tbl <- fm_vaf_tbl %>% group_by(deadbody,group2, var_id) %>% summarise(medVAF = median(VAF), q0 = quantile(VAF)[2], q1=quantile(VAF)[4])
med_tbl <- left_join(med_tbl, fm_vaf_tbl %>% select(deadbody, var_id, lineage_id2, step_n) %>% unique())
med_tbl <- med_tbl %>% ungroup() %>%  mutate(first_group = apply(.["lineage_id2"], 1,function(x) unlist(strsplit(x, '-'))[1]))
var_order <- fm_vaf_tbl %>% group_by(deadbody,var_id) %>% summarise(medVAF = median(VAF)) %>% arrange(desc(medVAF)) %>% pull(var_id)
g1 <- ggplot(med_tbl, aes(x=var_id, y=log10(medVAF+0.0001), color = group2))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymin=log10(q0+0.0001), ymax=log10(q1+0.0001)), alpha=0.5)+
  stat_smooth(aes(group = group2))+
  coord_cartesian(ylim=log10(c(0,0.5)+0.0001))+
  scale_color_manual(values=dl2_pal)+
  ggtitle(paste0(this_db, ' endoderm=liver'))+
  scale_x_discrete(limits = var_order)+
  scale_y_continuous(expand = c(0,0.05), labels = c(0.01, 0.05, 0.1, 0.5,1,2,3,4,5,10,20,30)*0.01, breaks= log10(c(0.01, 0.05, 0.1, 0.5,1,2,3,4,5,10,20,30)*0.01+0.0001))+
  theme_syp+theme(axis.text.x = element_blank(), panel.grid = element_blank())

fm_vaf_tbl <- f_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
fm_vaf_tbl <- fm_vaf_tbl %>% filter(tissue_class2 != 'liver')
fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group2 = dominant_layer2)
#fm_vaf_tbl <- fm_vaf_tbl %>% filter(group2 %in% c('mesoderm','endoderm','mixed','fresh_frozen_epidermis'))
ct_tbl <- fm_vaf_tbl %>% group_by(deadbody,group2, var_id) %>% count() %>% group_by(group2) %>% summarise(mean_n = mean(n))
print(ct_tbl)
incl_ids <- ct_tbl %>% filter(mean_n >=2) %>% pull(group2)
print(incl_ids)
fm_vaf_tbl <- fm_vaf_tbl %>% filter(group2 %in% incl_ids)

med_tbl <- fm_vaf_tbl %>% group_by(deadbody,group2, var_id) %>% summarise(medVAF = median(VAF), q0 = quantile(VAF)[2], q1=quantile(VAF)[4])
med_tbl <- left_join(med_tbl, fm_vaf_tbl %>% select(deadbody, var_id, lineage_id2, step_n) %>% unique())
med_tbl <- med_tbl %>% ungroup() %>%  mutate(first_group = apply(.["lineage_id2"], 1,function(x) unlist(strsplit(x, '-'))[1]))
var_order <- fm_vaf_tbl %>% group_by(deadbody,var_id) %>% summarise(medVAF = median(VAF)) %>% arrange(desc(medVAF)) %>% pull(var_id)
g2 <- ggplot(med_tbl, aes(x=var_id, y=log10(medVAF+0.0001), color = group2))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymin=log10(q0+0.0001), ymax=log10(q1+0.0001)), alpha=0.5)+
  stat_smooth(aes(group = group2))+
  coord_cartesian(ylim=log10(c(0,0.5)+0.0001))+
  scale_color_manual(values=dl2_pal)+
  ggtitle(paste0(this_db, ' endoderm=lung'))+
  scale_x_discrete(limits = var_order)+
  scale_y_continuous(expand = c(0,0.05), labels = c(0.01, 0.05, 0.1, 0.5,1,2,3,4,5,10,20,30)*0.01, breaks= log10(c(0.01, 0.05, 0.1, 0.5,1,2,3,4,5,10,20,30)*0.01+0.0001))+
  theme_syp+theme(axis.text.x = element_blank(), panel.grid = element_blank())

plot_grid(g1,g2, nrow=1, align="h", axis="tb")



#######Ectoderm-mesoendoderm divergence on phylogenetic tree #(first figure 3b)
f_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/f_vaf_tbl_rank.rds')
ff_vaf_tbl <- f_vaf_tbl %>% filter(time_of_lineage == 'early' & tissue_class !='cancer_tissue')
pval_co=0.05
this_db='DB6'
plot_list <- list()
for (this_db in c('DB6','DB8','DB9','DB10')){
  plot_list[[this_db]] <-Wilcox_Test_on_Phylo_Ecto_MesoEndo(ff_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co) # Main Figure 3b
}
pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/Fig3d.pdf', width=3, height=5)
plot_list[["DB6"]]
dev.off()

pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/ecto_mesoendo_diff_on_phylo.pdf', width=16, height=6)
plot_grid(plotlist=plot_list[1:4], nrow=1, align="h")
dev.off()

#legend
pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/ecto_mesoendo_diff_on_phylo_legend.pdf', width=12, height=6)
library(colorRamps)
total_break=50
col_breaks = seq(0, 0.5, length.out=total_break)
a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
rampCol1 <- colorRampPalette(c("#E65996", "#ED8BB6"))(n=a)
rampCol2 <- colorRampPalette(c("#ED8BB6","gray"))(n=b)
rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
my_palette1<- c(rampCol1, rampCol2, rampCol3)
 
plot.new()
plot(c(0,2),c(0,10),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Meso-Endoderm > Ectoderm')
colfunc1 <- rev(rampCol1)
legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
rasterImage(legend_image1, 0, 0, 0.3,5)
colfunc2 <- rev(rampCol2)
legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
rasterImage(legend_image2, 0, 5, 0.3,10)
#colfunc3 <- rev(rampCol3)
#legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
#rasterImage(legend_image3, 0, 10, 0.3,20)
text(x=0.5, y = c(0,5,10,20), labels =c(0,0.05, 0.1, 0.2), cex=2)

total_break=50
col_breaks = seq(0, 0.5, length.out=total_break)
a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1

b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
rampCol1 <- colorRampPalette(c("#245DB5", "#668ECC"))(n=a)
rampCol2 <- colorRampPalette(c("#668ECC","gray"))(n=b)
rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
my_palette2<- c(rampCol1, rampCol2, rampCol3)

plot.new()
plot(c(0,2),c(0,10),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Meso-Endoderm < Ectoderm')
colfunc1 <- rev(rampCol1)
legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
rasterImage(legend_image1, 0, 0, 0.3,5)
colfunc2 <- rev(rampCol2)
legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
rasterImage(legend_image2, 0, 5, 0.3,10)
#colfunc3 <- rev(rampCol3)
#legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
#rasterImage(legend_image3, 0, 10, 0.3,20)
text(x=0.5, y = c(0,5,10,50), labels =c(0,0.05, 0.1, 0.2), cex=2)

dev.off()


if(F){
  Wilcox_Test_on_Phylo_Ecto_Meso(ff_vaf_tbl, 'DB6', lng_dt, nwk_tbl, pval_co)
  Wilcox_Test_on_Phylo_Endo_Meso(ff_vaf_tbl, 'DB6', lng_dt, nwk_tbl, pval_co)
  Wilcox_Test_on_Phylo_Endo_Meso(ff_vaf_tbl, 'DB3', lng_dt, nwk_tbl, pval_co)
}


#######Ectoderm-mesoderm divergence on phylogenetic tree )
ff_vaf_tbl <- f_vaf_tbl %>% filter(time_of_lineage == 'early' & tissue_class !='cancer_tissue')
pval_co=0.05
for (this_db in c('DB6','DB8','DB9','DB10')){ 
  Wilcox_Test_on_Phylo_Ecto_Meso(ff_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co)
}


#######Left-right divergence on phylogenetic tree#############
pval_co=0.05
s_range = "EandT"
ff_vaf_tbl <- f_vaf_tbl %>% filter(time_of_lineage == 'early' & tissue_class !='cancer_tissue')
plot_list <- list()
for (this_db in c('DB3','DB6','DB8','DB9','DB10')){ 
  plot_list[[this_db]] <- Wilcox_Test_on_Phylo_LR(ff_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co, s_range)
}

pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/Lt_Rt_divergence_on_phylo_DB368910.pdf', width=12, height=6)
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

pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/Lt_Rt_divergence_on_phylo_legend.pdf')
library(colorRamps)
total_break=50
col_breaks = seq(0, 0.5, length.out=total_break)
a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1

b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
rampCol1 <- colorRampPalette(c('#023858', "#3690c0"))(n=a)
rampCol2 <- colorRampPalette(c("#3690c0","gray"))(n=b)
#rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
my_palette1<- c(rampCol1, rampCol2)
plot.new()
plot(c(0,2),c(0,15),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Rt. VAF > Lt. VAF')
colfunc1 <- rev(rampCol1)
legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
rasterImage(legend_image1, 0, 0, 0.3,5)
colfunc2 <- rev(rampCol2)
legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
rasterImage(legend_image2, 0, 5, 0.3,10)
#colfunc3 <- rev(rampCol3)
#legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
#rasterImage(legend_image3, 0, 10, 0.3,20)
text(x=0.5, y = c(0,5,10), labels =c(0,0.05, 0.1), cex=2)

total_break=50
col_breaks = seq(0, 0.5, length.out=total_break)
a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
rampCol1 <- colorRampPalette(c('#df1919',"#fcd57b"))(n=a)
rampCol2 <- colorRampPalette(c("#fcd57b","gray"))(n=b)
#rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
my_palette2<- c(rampCol1, rampCol2)
plot.new()
plot(c(0,2),c(0,15),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Rt. VAF < Lt. VAF')
colfunc1 <- rev(rampCol1)
legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
rasterImage(legend_image1, 0, 0, 0.3,5)
colfunc2 <- rev(rampCol2)
legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
rasterImage(legend_image2, 0, 5, 0.3,10)
#colfunc3 <- rev(rampCol3)
#legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
#rasterImage(legend_image3, 0, 10, 0.3,20)
text(x=0.5, y = c(0,5,10), labels =c(0,0.05, 0.1), cex=2)
dev.off()



#Get the timing of  Lt-rt divergence timing
f_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/f_vaf_tbl_rank.rds')
this_db='DB6'
p10_tbl <- tibble(deadbody = as.character(), min = as.numeric(), min_step = as.numeric(), p5 = as.numeric(), p10= as.numeric(), p25= as.numeric())
m=0
for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
  m=m+1
  ff_vaf_tbl <- f_vaf_tbl %>% filter(time_of_lineage == 'early' & tissue_class !='cancer_tissue' & deadbody == this_db)
  ff_vaf_tbl <- ff_vaf_tbl %>% mutate(var_g_id = paste0(lineage_id, '_', rank))
  this_id = ff_vaf_tbl$var_g_id[1]
  res_tbl <- tibble(deadbody = as.character(), var_g_id = as.character(), mean_lt_vaf = as.numeric(), mean_rt_vaf = as.numeric(), wilcox_p = as.numeric())
  n=0
  for (this_id in ff_vaf_tbl$var_g_id %>% unique()){
    n=n+1
    print(this_id)
    lt_vaf <- ff_vaf_tbl %>% filter(var_g_id == this_id) %>% filter(Source_side == 'lt') %>% pull(VAF)
    rt_vaf <- ff_vaf_tbl %>% filter(var_g_id == this_id) %>% filter(Source_side == 'rt') %>% pull(VAF)
    wilcox_p <- wilcox.test(lt_vaf, rt_vaf)$p.value
    res_tbl[n,] <- list(this_db, this_id, mean(lt_vaf), mean(rt_vaf), wilcox_p)
  }
  ff_vaf_tbl <- left_join(ff_vaf_tbl, res_tbl)
  tmp <- ff_vaf_tbl %>% filter(wilcox_p < 0.05) %>% select(var_id, assigned_ctie_mtn, step_n) %>% unique()
  min <- tmp$assigned_ctie_mtn %>% min()
  min_step <- tmp$step_n[tmp$assigned_ctie_mtn == min]
  p5 <- tmp$assigned_ctie_mtn %>% quantile(0.05) %>% round(.,2)
  p10 <- tmp$assigned_ctie_mtn %>% quantile(0.1) %>% round(.,2)
  p25 <- tmp$assigned_ctie_mtn %>% quantile(0.25) %>% round(.,2)
  g <- ggplot(tmp, aes(x=assigned_ctie_mtn))+
    geom_density()+
    geom_vline(xintercept=p10, color="red")+
    geom_vline(xintercept=p25, color="blue")+
    ggtitle(paste0(this_db, ' 10%=', p10, ' 25%=', p25))
  #print(g)
  p10_tbl[m,] <- list(this_db, min, min_step, p5, p10, p25)
}
p_tbl <- left_join(p10_tbl, m_rate_dt)
p_tbl <- p_tbl %>% mutate(min_stage = ifelse(min < M0, min/M0 +1, (min-M0)/M1+2),
                          p25_stage = ifelse(p25 < M0, p25/M0+1, (p25 - M0)/M1+2),
                          p10_stage = ifelse(p10 < M0, p10/M0+1, (p10 - M0)/M1+2),
                          p5_stage = ifelse(p5 < M0, p5/M0+1, (p5 - M0)/M1+2))
p_tbl <- p_tbl %>% mutate(final_min_stage = ifelse(min_step > min_stage, min_step, floor(min_stage)))



######## box plot of L1 and L2 grouped by left and right axis (first figure 3c)
p_list=list()
n=0
for (t_db in c('DB3','DB6','DB8','DB9','DB10')){
  n=n+1
  f_tbl <- f_vaf_tbl %>% filter(deadbody == t_db & step_n == 1 & rank == 1 &
                                  Source_side %in% c('lt','rt') & tissue_class != 'cancer_tissue')
  f_tbl <- f_tbl %>% group_by(sample_id, lineage_id, Source_side) %>% summarise(meanVAF=mean(VAF))

  ##change figure to meanVAF!
  p_list[[n]] <- ggplot(f_tbl, aes(x=lineage_id, y=meanVAF*100, color=Source_side))+
    geom_boxplot(outlier.colour = NA, 
                 position = position_dodge(width=0.9))+
    geom_point(position=position_jitterdodge(dodge.width=0.9, jitter.width=0.2), alpha=0.2, size=5)+
    stat_compare_means(aes(label = ..p.signif..), method="wilcox.test")+
    #stat_compare_means(aes(label = ..p.format..), method="wilcox.test")+
    scale_color_manual(values = lr_pal)+
    scale_y_continuous(limits=c(0,55),breaks=seq(0,50,10))+
    theme_syp+
    theme(axis.title.x = element_blank(), legend.position="None")+
    ggtitle(t_db)+ylab("VAF in bulk tissues (%)")
}
pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/DB368910_first_branch_Lt_Rt_difference.pdf', width=20, height=6)
plot_grid(plotlist=p_list, nrow=1, align="h")
dev.off()


############## Cranio caudal comparison in L1 and L2
p_list=list()
n=0
for (t_db in c('DB3','DB6','DB8','DB9','DB10')){
  n=n+1
  f_tbl <- f_vaf_tbl %>% filter(tissue_class != 'cancer_tissue' & deadbody == t_db & step_n == 1 & rank == 1 &
                                  Source_side_CC %in% c('cranio','caudal'))
  f_tbl <- f_tbl %>% group_by(lineage_id, sample_id, Source_side_CC) %>% summarise(meanVAF = mean(VAF))
  f_tbl$Source_side_CC <- factor(f_tbl$Source_side_CC, levels = c('cranio','caudal'))
  if (t_db != 'DB10'){
    p_list[[n]] <- ggplot(f_tbl, aes(x=lineage_id, y=meanVAF*100, color=Source_side_CC))+
      geom_boxplot(outlier.colour = NA, 
                   position = position_dodge(width=0.9))+
      geom_point(position=position_jitterdodge(dodge.width=0.9), alpha=0.5)+
      #stat_compare_means(aes(label = ..p.signif..), method="wilcox.test")+
      stat_compare_means(aes(label = ..p.format..), method="wilcox.test")+
      scale_color_manual(values=cc_pal)+
      scale_y_continuous(limits=c(0,55),breaks=seq(0,50,10))+
      theme_syp+
      theme(axis.title.x = element_blank(), legend.position = 'none')+
      ggtitle(t_db)+ylab("VAF in bulk tissues (%)")
  } else {
    p_list[[n]] <- ggplot(f_tbl, aes(x=lineage_id, y=meanVAF*100, color=Source_side_CC))+
      geom_boxplot(outlier.colour = NA, 
                   position = position_dodge(width=0.9))+
      geom_point(position=position_jitterdodge(dodge.width=0.9), alpha=0.5)+
      #stat_compare_means(aes(label = ..p.signif..), method="wilcox.test")+
      stat_compare_means(aes(label = ..p.format..), method="wilcox.test")+
      scale_color_manual(values=cc_pal)+
      scale_y_continuous(limits=c(0,55),breaks=seq(0,50,10))+
      theme_syp+
      theme(axis.title.x = element_blank())+
      ggtitle(t_db)+ylab("VAF in bulk tissues (%)")
  }

}
pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/DB368910_first_branch_Cranial_Caudal_difference.pdf', width=20, height=6)
plot_grid(plotlist=p_list, nrow=1, rel_widths = c(1,1,1,1,1.6))
dev.off()
############




#Window of our data
colnames(lng_dt)
tmp <- lng_dt %>% filter(time_of_lineage == 'early_to_late') %>% group_by(deadbody) %>% summarise(median_start = median(start_mtn),
                                                                                                  mean_start = mean(start_mtn))

tmp <- left_join(tmp, rate_dt %>% filter(type == 'M0') %>% rename(deadbody = ID, M0 = mean) %>% select(deadbody, M0))
tmp <- left_join(tmp, rate_dt %>% filter(type == 'M1') %>% rename(deadbody = ID, M1 = mean) %>% select(deadbody, M1))
tmp <- tmp %>% mutate(median_stage = (median_start - M0*2)/M1+2,
                      mean_stage = (mean_start - M0*2)/M1+2)


#########Group by stage (deprecated, we will use group by VAF)
if(F){
  this_db='DB9'
  t_M0 <- rate_dt %>% filter(ID == this_db & type == 'M0') %>% pull(mean)
  t_M1 <- rate_dt %>% filter(ID == this_db & type == 'M1') %>% pull(mean)
  ff_vaf_tbl <- f_vaf_tbl %>% filter(deadbody ==  this_db & time_of_lineage == 'early') 
  diff_tbl <- ff_vaf_tbl %>% select(lineage_id2, var_id, start_mtn, end_mtn, assigned_mean_mtn) %>% unique()
  if(this_db %in% c('DB3')){
    diff_tbl <- diff_tbl %>% mutate(diff_group = 'early')
  } else if (this_db %in% c('DB8', 'DB10')){
    diff_tbl <- diff_tbl %>% mutate(diff_group = ifelse(assigned_mean_mtn < t_M0*2 + t_M1*0, 'early','late'))
  } else if (this_db %in% c('DB9')){
    diff_tbl <- diff_tbl %>% mutate(diff_group = ifelse(assigned_mean_mtn < t_M0*2 + t_M1*1, 'early','late'))
  }else if(this_db == 'DB6') {
    diff_tbl <- diff_tbl %>% mutate(diff_group = ifelse(assigned_mean_mtn < t_M0*2 + t_M1*3, 'early',
                                                        ifelse(assigned_mean_mtn < t_M0*2 + t_M1*5, 'mid','late')))
  }
  diff_tbl %>% group_by(diff_group) %>% dplyr::count()
  fm_vaf_tbl <- f_vaf_tbl %>% filter(tissue_class != 'cancer_tissue' & time_of_lineage=='early' & deadbody == this_db)
  Draw_phylo_GroupByManual(fm_vaf_tbl,lng_dt,this_db, diff_tbl)
  Draw_Heatmap_Diff(fm_vaf_tbl, diff_tbl)
  Draw_tSNE_Diff(fm_vaf_tbl,
                 diff_tbl, 
                 color_type="LR",
                 seed=5)  
  
  Draw_tSNE_Diff(fm_vaf_tbl,
                 diff_tbl, 
                 color_type="AC1",
                 seed=1)  
  #seed/ DB6:3
  
  Draw_tSNE_total(fm_vaf_tbl,   #for DB3
                  color_type='LR',
                  seed=3)
  
  #L1 only and group by stage
  this_db='DB6'
  #t_M0 <- rate_dt %>% filter(ID == this_db & type == 'M0') %>% pull(mean)
  #t_M1 <- rate_dt %>% filter(ID == this_db & type == 'M1') %>% pull(mean)
  ff_vaf_tbl <- f_vaf_tbl %>% filter(deadbody ==  this_db & time_of_lineage == 'early') 
  lid_list <- lng_dt %>% filter(deadbody == this_db & time_of_lineage=='early_to_late') %>% pull(lineage_id)
  i=3
  t_lid = lid_list[i]
  print(t_lid)
  t_lid_v <- unlist(strsplit(t_lid,'-'))
  t_lid_line <- c()
  for (n in 1:length(t_lid_v)){
    t_lid_line = c(t_lid_line, paste(t_lid_v[1:n], collapse='-'))
  }
  t_lid_line
  fff_vaf_tbl <- ff_vaf_tbl %>% filter(lineage_id %in% t_lid_line)
  fff_vaf_tbl
  fff_vaf_tbl$var_id %>% unique() %>% length()
  Draw_tSNE_total(fff_vaf_tbl, 
                  color_type='AC')
  Draw_tSNE_total(fff_vaf_tbl, 
                  color_type='CC')
  
  diff_tbl <- ff_vaf_tbl %>% select(lineage_id2, var_id, start_mtn, end_mtn, assigned_mean_mtn) %>% unique()
  if(this_db %in% c('DB3')){
    diff_tbl <- diff_tbl %>% mutate(diff_group = 'early')
  } else if (this_db %in% c('DB8', 'DB10')){
    diff_tbl <- diff_tbl %>% mutate(diff_group = ifelse(assigned_mean_mtn < t_M0*2 + t_M1*0, 'early','late'))
  } else if (this_db %in% c('DB9')){
    diff_tbl <- diff_tbl %>% mutate(diff_group = ifelse(assigned_mean_mtn < t_M0*2 + t_M1*1, 'early','late'))
  }else if(this_db == 'DB6') {
    diff_tbl <- diff_tbl %>% mutate(diff_group = ifelse(assigned_mean_mtn < t_M0*2 + t_M1*3, 'early',
                                                        ifelse(assigned_mean_mtn < t_M0*2 + t_M1*5, 'mid','late')))
  }
  diff_tbl %>% group_by(diff_group) %>% dplyr::count()
  fm_vaf_tbl <- f_vaf_tbl %>% filter(tissue_class != 'cancer_tissue' & time_of_lineage=='early' & deadbody == this_db)
  
  Draw_tSNE_Diff(fm_vaf_tbl,
                 diff_tbl, 
                 color_type="CC",
                 seed=3)  #DB10:5 
}

# VAF heatmap for L1-2-2-1-1-1 (Figure for R3-A05)
mx <- f_vaf_tbl %>% filter(deadbody == 'DB6' & lineage_id == 'L1-2-2-1-1-1') %>% select(var_id, sample_id, VAF) %>% spread(key=sample_id, value=VAF) %>% as.data.frame %>% column_to_rownames('var_id') %>% as.matrix()
f_vaf_tbl$VAF %>% summary()
annot_cols <- f_vaf_tbl$anatomy_class3 %>% unique() 
annot_cols <- f_vaf_tbl %>% filter(deadbody == 'DB6' & lineage_id == 'L1-2-2-1-1-1') %>% select(sample_id, anatomy_class3) %>% unique() %>% as.data.frame() %>% column_to_rownames('sample_id')


mx <- mx*100
pdf("/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/Fig_R3A05.pdf", width=6, height=4)
pheatmap(mx, annotation_col = annot_cols, show_colnames = F, clustering_method = "ward.D2")
dev.off()

##############tSNE plot (DB6 for Main Fig 3d, others for EDFig 7b)
this_db = 'DB6'
plot_list = list()
for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
  print(this_db)
  f_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/f_vaf_tbl_rank.rds')
  f_vaf_tbl <- f_vaf_tbl %>% mutate(a1t2_class = ifelse(anatomy_class1 == 'internal_organ',tissue_class2, anatomy_class1))
  f_vaf_tbl <- f_vaf_tbl %>% mutate(tissue_class3 = ifelse(tissue_class2 %in% c('skin','fat','muscle','skin_epi','skin_derm','fascia','cancer_tissue','cornea','peripheral_vessel'), 'non_internal_organ',tissue_class2))
  f_vaf_tbl$tissue_class3[f_vaf_tbl$tissue_class3 == 'colon'] <- 'bowel'
  f_vaf_tbl$tissue_class3[f_vaf_tbl$tissue_class3 == 'intestine'] <- 'bowel'
  fm_vaf_tbl <- f_vaf_tbl %>% filter(tissue_class != 'cancer_tissue' & time_of_lineage=='early' & deadbody == this_db)
  n_total <- fm_vaf_tbl$sample_id %>% unique() %>% length()
  n_total
  #grouping
  cutoff1=5;cutoff2=0.5
  if (this_db == 'DB6'){
    diff_tbl <- fm_vaf_tbl  %>% group_by(lineage_id2, var_id) %>% summarise(medVAF = median(VAF)) %>% mutate(diff_group = ifelse(medVAF >= cutoff1*0.01, 'early',
                                                                                                                                 ifelse(medVAF < cutoff2*0.01, 'late', 'mid')))
  }else if (this_db %in% c('DB3','DB10')){
    diff_tbl <- fm_vaf_tbl  %>% group_by(lineage_id2, var_id) %>% summarise(medVAF = median(VAF)) %>% mutate(diff_group= 'early')
  } else {
    diff_tbl <- fm_vaf_tbl  %>% group_by(lineage_id2, var_id) %>% summarise(medVAF = median(VAF)) %>% mutate(diff_group = ifelse(medVAF >= cutoff1*0.01, 'early','late'))
  }
  #Draw_phylo_GroupByManual(fm_vaf_tbl,lng_dt, this_db,diff_tbl)
  nrow(diff_tbl)
  diff_tbl %>% group_by(diff_group) %>% count()                                                                                                                  
  
  
  if(T){
    #make mx in each group
    mx <- fm_vaf_tbl %>% filter(time_of_lineage == 'early') %>% select(var_id, sample_id, VAF) %>% spread(key=sample_id, value=VAF) %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
    mx <- mx[rowSums(mx, na.rm = T)!=0,]
    library(missForest) #NA value imputation using missForest
    Early <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'early' & is.na(diff_tbl$diff_group) == F ],]
    Early <- Early[rowSums(is.na(Early)) < ncol(Early)/2,]
    Early <- missForest(t(Early))$ximp
    Early <- Early[rowSums(Early) >0,]
    Early <- t(scale(t(Early)))
    if('late' %in% diff_tbl$diff_group %>% unique()){
      Late <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'late' & is.na(diff_tbl$diff_group) == F ],]
      Late <- Late[rowSums(is.na(Late)) < ncol(Late)/2,]
      Late <- missForest(t(Late))$ximp
      Late <- Late[rowSums(Late) >0,]
      Late <- t(scale(t(Late)))
    }
    if("mid" %in% diff_tbl$diff_group %>% unique){
      Mid <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'mid' & is.na(diff_tbl$diff_group) == F ],]
      Mid <- Mid[rowSums(is.na(Mid)) < ncol(Mid)/2,]
      Mid <- missForest(t(Mid))$ximp
      Mid <- Mid[rowSums(Mid) >0,]
      Mid <- t(scale(t(Mid)))
    }
    #run tsne
    library(Rtsne)
    if (diff_tbl$diff_group %>% unique()%>% length() == 1){
      timing_group = c('Early')
    }else if("mid" %in% diff_tbl$diff_group %>% unique){
      timing_group = c("Early","Mid","Late")
    } else{
      timing_group = c("Early","Late")
    }
    print(timing_group)
    out_fd =  "/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/"
    #draw left-right coloring

    set.seed(3)
    p_list=list()
    n=0
    for (timing in timing_group){
      print(timing)
      n = n+1
      df_tsne <- Rtsne(get(timing), perplexity=10)
      df_out <- as.data.frame(df_tsne$Y)
      rownames(df_out) <- rownames(get(timing))
      colnames(df_out) <- c('tSNE1','tSNE2')
      dt <- left_join(df_out %>% as.data.frame() %>% rownames_to_column('sample_id') %>% as.tibble(), 
                      fm_vaf_tbl %>% select(sample_id, starts_with('anatomy'), starts_with('tissue'), starts_with("Source")) %>% unique())
      left_join(dt %>% select(sample_id, tSNE1, tSNE2), pmeta_dt %>% select(sample_id, sample_new_id)) %>% select(sample_new_id, tSNE1, tSNE2, sample_id) %>% write_tsv(paste0 (out_fd, this_db, '_tSNE_',timing,'_LR.tsv'))
      if(n < length(timing_group)){
        p_list[[n]] <- ggplot(dt, aes(x=tSNE1, y=tSNE2, color=Source_side))+
          geom_point(size=5, alpha=0.9)+
          scale_color_manual(values=lr_pal)+
          theme_syp+
          theme(legend.position = 'none')
          #ggtitle(timing)  
      } else{
        p_list[[n]] <- ggplot(dt, aes(x=tSNE1, y=tSNE2, color=Source_side))+
          geom_point(size=5, alpha=0.9)+
          scale_color_manual(values=lr_pal)+
          theme_syp+
          theme(legend.position = 'none')
          #ggtitle(timing)  
      }
    }

    
    tc3_pal <-c("grey20","#8388D6","#DB9E51","#D9CDE1","#DDE265","#D5938F","#CCD4AD","#B250E1","#83DA94","#9EDDFF","#73B2CF","#E3478B","#81E4D4","#DC89D2")
    names(tc3_pal)<- c("non_internal_organ","lung","bowel","heart","large_vessel","kidney","intestine","pancreas","liver","spleen","blood","stomach","ovary","uterus")
    set.seed(3)
    for (timing in timing_group){
      print(timing)
      n = n+1
      df_tsne <- Rtsne(get(timing), perplexity=10)
      df_out <- as.data.frame(df_tsne$Y)
      rownames(df_out) <- rownames(get(timing))
      colnames(df_out) <- c('tSNE1','tSNE2')
      dt <- left_join(df_out %>% as.data.frame() %>% rownames_to_column('sample_id') %>% as.tibble(), 
                      fm_vaf_tbl %>% select(sample_id, starts_with('anatomy'), starts_with('tissue'), starts_with("Source")) %>% unique())
      if(n%%length(timing_group) >0){
        p_list[[n]] <- ggplot(dt, aes(x=tSNE1, y=tSNE2, color=tissue_class3, shape=anatomy_class1))+
          geom_point(size=5, alpha=0.9)+
          scale_color_manual(values=tc3_pal)+
          scale_shape_manual(values = c("internal_organ" = 19, "UE" = 17,"LE"=15, "trunk"= 3,"HN"=7))+
          theme_syp+
          theme(legend.position = 'none')
          #ggtitle(timing)  
      } else{
        p_list[[n]] <- ggplot(dt, aes(x=tSNE1, y=tSNE2, color=tissue_class3, shape=anatomy_class1))+
          geom_point(size=5, alpha=0.9)+
          scale_color_manual(values=tc3_pal)+
          scale_shape_manual(values = c("internal_organ" = 19, "UE" = 17,"LE"=15, "trunk"= 3,"HN"=7))+
          theme_syp+
          theme(legend.position = 'none')
          #ggtitle(timing)  
      }
    }
    plot_list[[this_db]] <- plot_grid(plotlist = p_list, nrow=2, align="hv")
    }
}

pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/EDFig7b_DB3.pdf', width=4, height=8)
plot_list[['DB3']]
dev.off()

pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/EDFig7b_DB8.pdf', width=8, height=8)
plot_list[['DB8']]
dev.off()

pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/EDFig7b_DB9.pdf', width=8, height=8)
plot_list[['DB9']]
dev.off()

pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/EDFig7b_DB10.pdf', width=4, height=8)
plot_list[['DB10']]
dev.off()

#draw legend
pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/EDFig7b_legend.pdf', width=8, height=8)
new_pal <- tc3_pal[names(tc3_pal) != "non_internal_organ"]
lgd1 <- ComplexHeatmap::Legend(at = names(new_pal), type = "points", background = NULL,
                    legend_gp = gpar(col = new_pal), size=unit(5,'mm'), grid_height=unit(8,'mm'), ncol=1
                    )
lgd2 <- ComplexHeatmap::Legend(at = c("internal_organ", "UE","LE", "trunk","HN"), type = "points", background = NULL,
                               legend_gp = gpar(col = "grey20" ), size=unit(5,'mm'), pch = c(19,17,15,3,7), grid_height=unit(8,'mm')
)
                              
lgd3 <- ComplexHeatmap::Legend(at = names(lr_pal), type = "points", background = NULL,
                               legend_gp = gpar(col = lr_pal ),  size=unit(5,'mm'), grid_height=unit(8,'mm')
)
grid.newpage()
pushViewport(viewport(x=0, y=0.3, width = 0.5, height = 0.7, just = c('left', 'bottom')))
grid.draw(lgd1)
popViewport(1)
pushViewport(viewport(x=0.5, y=0.5, width = 0.5, height = 0.5, just = c('left', 'bottom')))
grid.draw(lgd2)
popViewport(1)
pushViewport(viewport(x=0, y=0, width = 0.5, height = 0.3, just = c('left', 'bottom')))
grid.draw(lgd3)
popViewport(1)
dev.off()


#calculate the cell stage
diff_tbl <- left_join(diff_tbl, f_vaf_tbl %>% select(assigned_mean_mtn, start_mtn, var_id) %>% unique())
n_group <- diff_tbl$diff_group %>% unique() %>% length()

if (n_group == 3){
  lower.limit <- min(diff_tbl$assigned_mean_mtn[diff_tbl$diff_group %in% c('early','mid')])
  upper.limit <- max(diff_tbl$assigned_mean_mtn[diff_tbl$diff_group %in% c('early','mid')])
  early.density <- density(subset(diff_tbl, diff_group == "early")$assigned_mean_mtn, from = lower.limit, to = upper.limit, n = 2^10)
  mid.density <- density(subset(diff_tbl, diff_group == "mid")$assigned_mean_mtn, from = lower.limit, to = upper.limit, n = 2^10)
  density.difference <- early.density$y - mid.density$y
  em.intersection.point <- early.density$x[which(diff(density.difference > 0) != 0) + 1]
  
  lower.limit <- min(diff_tbl$assigned_mean_mtn[diff_tbl$diff_group %in% c('mid','late')])
  upper.limit <- max(diff_tbl$assigned_mean_mtn[diff_tbl$diff_group %in% c('mid','late')])
  mid.density <- density(subset(diff_tbl, diff_group == "mid")$assigned_mean_mtn, from = lower.limit, to = upper.limit, n = 2^10)
  late.density <- density(subset(diff_tbl, diff_group == "late")$assigned_mean_mtn, from = lower.limit, to = upper.limit, n = 2^10)
  density.difference <- mid.density$y - late.density$y
  ml.intersection.point <- mid.density$x[which(diff(density.difference > 0) != 0) + 1]
  
  g <- ggplot(diff_tbl,aes(assigned_mean_mtn, fill = diff_group)) + geom_density(alpha = 0.2) + 
    geom_vline(xintercept = em.intersection.point, color = "red")+
    geom_vline(xintercept = ml.intersection.point, color="red")+
    ggtitle(this_db)
  print(g)
  
  t_m0 <- rate_dt %>% filter(type== 'M0' & ID == this_db) %>% pull(mean)
  t_m1 <- rate_dt %>% filter(type== 'M1' & ID == this_db) %>% pull(mean)
  print(em.intersection.point)
  print(ml.intersection.point)
  em_stage <- (em.intersection.point - t_m0)/t_m1 +2
  ml_stage <- (ml.intersection.point - t_m0)/t_m1 +2
  print(em_stage)
  print(ml_stage)
} else if (n_group ==2){
  lower.limit <- min(diff_tbl$assigned_mean_mtn[diff_tbl$diff_group %in% c('early','late')])
  upper.limit <- max(diff_tbl$assigned_mean_mtn[diff_tbl$diff_group %in% c('early','late')])
  early.density <- density(subset(diff_tbl, diff_group == "early")$assigned_mean_mtn, from = lower.limit, to = upper.limit, n = 2^10)
  late.density <- density(subset(diff_tbl, diff_group == "late")$assigned_mean_mtn, from = lower.limit, to = upper.limit, n = 2^10)
  density.difference <- early.density$y - late.density$y
  el.intersection.point <- early.density$x[which(diff(density.difference > 0) != 0) + 1]
  g <- ggplot(diff_tbl,aes(assigned_mean_mtn, fill = diff_group)) + geom_density(alpha = 0.2) + 
    geom_vline(xintercept = el.intersection.point, color = "red")+
    ggtitle(this_db)
  print(g)
  t_m0 <- rate_dt %>% filter(type== 'M0' & ID == this_db) %>% pull(mean)
  t_m1 <- rate_dt %>% filter(type== 'M1' & ID == this_db) %>% pull(mean)
  el_stage <- (el.intersection.point - t_m0)/t_m1 +2
  print(el.intersection.point)
  print(el_stage)
  
} else{
  print('NA')
}



# caclulate our window using targeted seq
f_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/f_vaf_tbl_rank.rds')
tmp <- f_vaf_tbl %>% filter(time_of_lineage == 'early') %>% select(deadbody, var_id, lineage_id, assigned_ctie_mtn) %>% unique() %>%
  group_by(deadbody) %>% summarise(p90 = quantile(assigned_ctie_mtn, 0.9),
                                   p75 = quantile(assigned_ctie_mtn, 0.75),
                                   p95 = quantile(assigned_ctie_mtn, 0.95))
tmp <- left_join(tmp, m_rate_dt)
tmp <- tmp %>% mutate(p90_stage = (p90 - M0)/M1+2,
                      p75_stage = (p75 - M0)/M1+2,
                      p95_stage = (p95 - M0)/M1+2)
tmp



