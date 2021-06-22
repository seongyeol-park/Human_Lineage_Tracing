#intended error
stop('abc')

#load library
library(tidyverse)
library(ggsci)
library(scales)
library(ggtree)
library(ggpubr)
library(cowplot)

#load data
merged_m_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Target_seq_related/merged_m_vaf_tbl.rds')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_200908.txt')
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_200908.txt') %>% filter(current_final_for_lineage == 'Y')

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


#color setting
pal_combi <- c(pal_npg("nrc")(10), pal_d3('category10')(10), pal_aaas('default')(10))
show_col(pal_combi[c(20,22,23,24)])
group_pal <- pal_combi[c(20,22,12, 23,24)]
names(group_pal) <- c('mixed','fresh_frozen_epidermis','formalin_fixed_epidermis','endoderm','mesoderm')
dl2_pal <- pal_combi[c(8, 12, 4, 3, 13)]
#show_col(dl2_pal)
names(dl2_pal) <- c('ectoderm','ecto_mesoderm','mesoderm','meso_endoderm','endoderm')
db_pal <- pal_npg('nrc')(10)[c(1:5,9,7)]
names(db_pal) <- unique(meta_dt$deadbody)
lr_pal = c(pal_combi[8], pal_combi[4])
names(lr_pal) = c('rt','lt')
cc_pal = c(pal_combi[24], pal_combi[23])
names(cc_pal) = c('cranio','caudal')
#theme setting
theme_syp = theme_bw() + theme(axis.text = element_text(size=15), axis.title = element_text(size=18), strip.text = element_text(size=15), panel.border = element_blank(), axis.line = element_line(), title = element_text(size=18))

source('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Fate_priming_functions.R')




##############germ layer divergence########
#Main Fig 4A
if(T){
  this_db = 'DB10'
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
  fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group2 = dominant_layer2)
  #fm_vaf_tbl <- fm_vaf_tbl %>% filter(group2 %in% c('mesoderm','endoderm','mixed','fresh_frozen_epidermis'))
  ct_tbl <- fm_vaf_tbl %>% group_by(deadbody,group2, var_id) %>% count() %>% group_by(group2) %>% summarise(mean_n = mean(n))
  print(ct_tbl)
  incl_ids <- ct_tbl %>% filter(mean_n >=2) %>% pull(group2)
  fm_vaf_tbl <- fm_vaf_tbl %>% filter(group2 %in% incl_ids)
  med_tbl <- fm_vaf_tbl %>% group_by(deadbody,group2, var_id) %>% summarise(medVAF = median(VAF), q0 = quantile(VAF)[2], q1=quantile(VAF)[4])
  med_tbl <- left_join(med_tbl, fm_vaf_tbl %>% select(deadbody, var_id, lineage_id2, step_n) %>% unique())
  med_tbl <- med_tbl %>% ungroup() %>%  mutate(first_group = apply(.["lineage_id2"], 1,function(x) unlist(strsplit(x, '-'))[1]))
  var_order <- fm_vaf_tbl %>% group_by(deadbody,var_id) %>% summarise(medVAF = median(VAF)) %>% arrange(desc(medVAF)) %>% pull(var_id)
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
  ggplot(med_tbl, aes(x=var_id, y=medVAF, color = group2))+
    geom_point(size=3, alpha=0.5)+
    stat_smooth(aes(group = group2))+
    scale_color_manual(values=dl2_pal)+
    ggtitle(this_db)+
    scale_x_discrete(limits = var_order)+
    scale_y_continuous(limits=c(0,0.05))+
    theme_syp+theme(axis.text.x = element_blank(), panel.grid = element_blank())
}

#######Ectoderm-mesoendoderm summary plot
#Main Fig 4B
m_ch_dt <- tibble()
pval_co=0.05
for (this_db in c('DB6','DB8','DB9','DB10')){
  changing_lines <- Wilcox_Test_on_Phylo_Ecto(merged_m_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co)
  m_ch_dt <- bind_rows(m_ch_dt, changing_lines)
}
tmp_tbl <- left_join(
  m_ch_dt,
  lng_dt %>% select(deadbody, lineage_id2, early_desc_n)
)
tmp_tbl <- tmp_tbl %>% mutate(new_mtn = start_mtn + rank) %>% group_by(deadbody, new_mtn) %>% summarise(added_lineage_n = sum(early_desc_n))
tmp_tbl <- tmp_tbl %>% mutate(cum_sum = cumsum(added_lineage_n))
tmp_tbl <- bind_rows(tmp_tbl, tibble(deadbody = tmp_tbl$deadbody %>% unique(), new_mtn = 0,added_lineage_n = 0, cum_sum=0))
tmp_tbl <- left_join(tmp_tbl, lng_dt %>% filter(time_of_lineage== 'early_to_late') %>% group_by(deadbody) %>% summarise(total_lineage_n=n()))
ggplot(tmp_tbl, aes(x=new_mtn, y=cum_sum/total_lineage_n, color=deadbody))+
  geom_point(size=3) +
  geom_line(aes(group = deadbody), linetype='dashed')+
  xlab("No. of mutations")+ylab("Proportion of biased lineages")+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(limits=c(0,35))+
  scale_color_manual(values=db_pal)+
  ggtitle('Ectoderm - Mesoendoderm bias timing')+
  theme_syp


########DB3 and DB6. box plot of L1 and L2 grouped by left and right axis
t_db = 'DB3'
f_tbl <- merged_m_vaf_tbl %>% filter(deadbody == t_db & step_n == 1 & DP >= 100 &
                                       Source_side %in% c('lt','rt'))
ggplot(f_tbl, aes(x=lineage_id2, y=VAF, color=Source_side))+
  geom_boxplot(outlier.colour = NA, 
               position = position_dodge(width=0.9))+
  geom_point(position=position_jitterdodge(dodge.width=0.9))+
  theme(axis.title.x = element_blank())+
  scale_x_discrete(labels = c('L1','L2'))+
  ggtitle(t_db)

t_db = 'DB6'
f_tbl <- merged_m_vaf_tbl %>% filter(deadbody == t_db & step_n == 1 & DP >= 100 &
                                       Source_side %in% c('lt','rt'))
ggplot(f_tbl, aes(x=lineage_id2, y=VAF, color=Source_side))+
  geom_boxplot(outlier.colour = NA, 
               position = position_dodge(width=0.9))+
  geom_point(position=position_jitterdodge(dodge.width=0.9))+
  theme(axis.title.x = element_blank())+
  scale_x_discrete(labels = c('L1','L2'))+
  ggtitle(t_db)
  



#########Representative phylogenetic tree of DB6 (EDFig 4A)
#1. Ectoderm-MesoEndoderm
this_db = 'DB6'
pval_co=0.05
Wilcox_Test_on_Phylo_Ecto(merged_m_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co)
#2. Left-Right
this_db='DB6'
pval_co=0.05
s_range='EandT' #c('EandT','Eonly') #Extremities and trunk or Extremities only
Wilcox_Test_on_Phylo_LR(merged_m_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co, s_range)
#3. Craino-Caudal
pval_co=0.05
s_range = 'HTE' #c('HTE','HandLE') # head trunk extremities, head and LEx
this_db='DB6'
Wilcox_Test_on_Phylo_CC(merged_m_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co, s_range)

########Left-righ and Cranio-cauda axis priming per deadbody (EDFig 4B)
m_ch_dt <- tibble()
for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
  pval_co=0.05
  s_range='EandT'
  changing_lines <- Wilcox_Test_on_Phylo_LR(merged_m_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co,s_range)
  changing_lines <- changing_lines %>% mutate(group = 'LRaxis')
  m_ch_dt <- bind_rows(m_ch_dt, changing_lines)
  s_range = 'HTE'
  changing_lines <- Wilcox_Test_on_Phylo_CC(merged_m_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co,s_range)
  changing_lines <- changing_lines %>% mutate(group = 'CCaxis')
  m_ch_dt <- bind_rows(m_ch_dt, changing_lines)
}
tmp_tbl <- left_join(
  m_ch_dt,
  lng_dt %>% select(deadbody, lineage_id2, early_desc_n)
)
tmp_tbl <- tmp_tbl %>% mutate(new_mtn = start_mtn + rank) %>% group_by(deadbody, group, new_mtn) %>% summarise(added_lineage_n = sum(early_desc_n))
tmp_tbl <- tmp_tbl %>% mutate(cum_sum = cumsum(added_lineage_n))
tmp_tbl <- bind_rows(tmp_tbl, tibble(deadbody = rep(tmp_tbl$deadbody %>% unique(),2), group = c(rep('LRaxis',5),rep('CCaxis',5)), new_mtn = 0,added_lineage_n = 0, cum_sum=0))
tmp_tbl <- left_join(tmp_tbl, lng_dt %>% filter(time_of_lineage== 'early_to_late') %>% group_by(deadbody) %>% summarise(total_lineage_n=n()))
plot_list = list()
n=0
for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
  n=n+1
  plot_list[[n]] <- ggplot(subset(tmp_tbl, deadbody == this_db), aes(x=new_mtn, y=cum_sum/total_lineage_n, color=group))+
    geom_point(size=3)+
    geom_line(aes(group = group), linetype='dashed')+
    xlab("No. of early mutations")+ylab("Proportion of primed lineages")+
    scale_y_continuous(limits=c(0,1))+
    scale_x_continuous(limits=c(0,35))+
    ggtitle(this_db)+
    theme_syp
}
plot_grid(plotlist = plot_list, nrow=2)


