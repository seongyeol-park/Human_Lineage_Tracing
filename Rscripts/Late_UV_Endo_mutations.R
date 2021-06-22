# Late-period mutations

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
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_200908.txt') %>% filter(current_final_for_lineage == 'Y')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_200908.txt')
DB_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_DB_200918.txt')

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

#############recurrent mutations 

mt_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perLineage_EarlyToLate_link/DBto_EtoL_pointmts_merged.txt', col_types = cols(`#CHROM`='c'))
if(F){
  mt_dt %>% filter(caseNo >1) %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perLineage_EarlyToLate_link/DBto_EtoL_pointmts_merged.txt.2more')
  mt_dt %>% filter(caseNo >1 & (nchar(REF) >1 | nchar(ALT) >1)) %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perLineage_EarlyToLate_link/DBto_EtoL_pointmts_merged.txt.2more.indel')
}


########################late mutation recurrency
library(dndscv)
mt_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perLineage_EarlyToLate_link/DBto_EtoL_pointmts_merged.txt', col_types = cols(`#CHROM`='c'))
mt_dt <- mt_dt %>% select(case_list, `#CHROM`, POS, REF, ALT) 
colnames(mt_dt) <- c('sampleID','chr','pos','ref','mut')
dndsout <- dndscv(mt_dt)
sel_cv = dndsout$sel_cv
print(head(sel_cv, n=30))
sel_cv$qglobal_cv %>% summary()

annot_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perLineage_EarlyToLate_link/DBto_306lineage_merged.txt.n_cos')
annot_dt %>% arrange(desc(n_COSMIC)) %>% select(`#CHROM`, POS, REF, ALT, Gene_refGene, ExonicFunc_refGene, AAChange_refGene,case_list, n_COSMIC ) %>% slice(1:20) %>% View()
annot_dt %>% filter(Gene_refGene == 'TP53' & Func_refGene == 'exonic') %>% View()


###########################signature related################################
#mutational burden per sample per signature
if(F){
  m_meta_dt <- meta_dt %>% mutate_at(vars(starts_with("v3")), list(~ .*n_snv*0.01))
  ct_tbl <- m_meta_dt %>% select(deadbody,sample_id, starts_with("v3")) %>% gather(-deadbody, -sample_id, key='signature',value='count')
  ct_tbl$sample_id <- factor(ct_tbl$sample_id, levels = m_meta_dt %>% arrange(desc(n_snv)) %>% pull(sample_id))
  
  
  ggplot(ct_tbl, aes(x=sample_id, y=count, fill = signature))+
    geom_bar(stat='identity')+
    scale_fill_manual(values=sig_pal)+
    theme_syp+theme(axis.text.x = element_blank())+
    facet_wrap(~deadbody, scales = "free")
}

#signature count multiple correlation
if(F){
  m_meta_dt <- meta_dt %>% mutate_at(vars(starts_with("v3")), list(~ .*n_snv*0.01))
  mx <- m_meta_dt  %>% select(sample_id, starts_with('v3')) %>% as.data.frame() %>% column_to_rownames('sample_id') %>% as.matrix()
  library(PerformanceAnalytics)
  chart.Correlation(mx)
}

#signtaure 7 and other correlation in skin_fb
if(F){
  m_meta_dt <- meta_dt %>% mutate_at(vars(starts_with("v3")), list(~ .*n_snv*0.01)) %>% mutate(v3_7total = v3_7a + v3_7b + v3_7c + v3_7d)
  f_meta_dt <- m_meta_dt %>% filter(Cell_type == 'skin_fb')
  
  #5 vs 7
  ggplot(f_meta_dt, aes(x=v3_7total, y=v3_5, color = deadbody))+
    geom_point()+
    geom_smooth(method="lm")+
    scale_color_manual(values = db_pal)+
    theme_syp
  lm(f_meta_dt$v3_5 ~ f_meta_dt$v3_7total + f_meta_dt$age) %>% summary()
  
  #1 vs 7
  ggplot(f_meta_dt, aes(x=v3_7total, y=v3_1, color = deadbody))+
    geom_point()+
    geom_smooth(method="lm")+
    scale_color_manual(values = db_pal)+
    theme_syp
  lm(f_meta_dt$v3_1 ~ f_meta_dt$v3_7total + f_meta_dt$age) %>% summary()
  ff_meta_dt <- f_meta_dt %>% filter(v3_1 > 0)
  ggplot(ff_meta_dt, aes(x=v3_7total, y=v3_1, color = deadbody))+
    geom_point()+
    geom_smooth(method="lm")+
    scale_color_manual(values = db_pal)+
    theme_syp
  lm(ff_meta_dt$v3_1 ~ ff_meta_dt$v3_7total + ff_meta_dt$age) %>% summary()
  #18 vs 7
  ggplot(f_meta_dt, aes(x=v3_7total, y=v3_18, color = deadbody))+
    geom_point()+
    geom_smooth(method="lm")+
    scale_color_manual(values = db_pal)+
    theme_syp
  lm(f_meta_dt$v3_18 ~ f_meta_dt$v3_7total + f_meta_dt$age) %>% summary() 
}

#age mutation correlation
if(F){
  m_meta_dt <- meta_dt %>% mutate_at(vars(starts_with("v3")), list(~ .*n_snv*0.01)) %>% mutate(v3_7total = v3_7a + v3_7b + v3_7c + v3_7d)
  m_meta_dt$Cell_type %>% table()
  f_meta_dt <- m_meta_dt %>% filter(Cell_type == 'skin_fb')
  #total SNV
  ggplot(f_meta_dt, aes(x=age, y=n_snv))+
    geom_boxplot(aes(fill = deadbody))+
    geom_smooth(method="lm")+
    scale_x_continuous(limits=c(0,100))+
    scale_fill_manual(values=db_pal)+
    theme_syp+
    ggtitle('total SNV')
  lm(f_meta_dt$n_snv ~ f_meta_dt$age) %>% summary
  #v3_1,5,18
  ggplot(f_meta_dt, aes(x=age, y=v3_1 + v3_5 + v3_18))+
    geom_boxplot(aes(fill = deadbody))+
    geom_smooth(method="lm")+
    scale_x_continuous(limits=c(0,100))+
    scale_y_continuous(limits=c(0,4500))+
    scale_fill_manual(values=db_pal)+
    theme_syp+
    ggtitle('sig1+5+18')
  f_meta_dt <- f_meta_dt %>% mutate(v3_1518 = v3_1 + v3_5 + v3_18)
  lm(f_meta_dt$v3_1518 ~ f_meta_dt$age) %>% summary()
  #signature count correction
  f_meta_dt <- f_meta_dt %>% mutate(v3_5cor = v3_5 - (v3_7total*0.1087)) %>% 
    mutate(v3_18cor = v3_18 + (v3_7total*0.0106))
  ggplot(f_meta_dt, aes(x=age, y=v3_1 + v3_5cor + v3_18cor))+
    geom_boxplot(aes(fill = deadbody))+
    geom_smooth(method="lm")+
    scale_x_continuous(limits=c(0,100))+
    scale_y_continuous(limits=c(0,4500))+
    scale_fill_manual(values=db_pal)+
    theme_syp+
    ggtitle('sig1+5cor+18cor')
  f_meta_dt <- f_meta_dt %>% mutate(v3_15c18c = v3_1 + v3_5cor + v3_18cor)
  lm(f_meta_dt$v3_15c18c ~ f_meta_dt$age) %>% summary()
}

#tissue variation of mutations
if(F){
  m_meta_dt <- meta_dt %>% mutate_at(vars(starts_with("v3")), list(~ .*n_snv*0.01)) %>% mutate(v3_7total = v3_7a + v3_7b + v3_7c + v3_7d) %>%
    mutate(v3_5cor = v3_5 - (v3_7total*0.1087), v3_18cor = v3_18 + (v3_7total*0.0106)) %>%
    mutate(n_snv_endo = v3_1 + v3_5cor + v3_18cor)
  f_meta_dt <- m_meta_dt %>% filter(Cell_type == 'skin_fb')
  multi_ids <- f_meta_dt %>% group_by(tissue_id) %>% count() %>% filter(n>1) %>% pull(tissue_id)
  f_meta_dt <- f_meta_dt %>% filter(tissue_id %in% multi_ids)
  fm_meta_dt <- f_meta_dt %>% select(tissue_id, v3_7total, n_snv_endo) %>% gather(-tissue_id, key='group', value='n_muts')
  x_order <- f_meta_dt %>% group_by(deadbody, tissue_id) %>% summarise(med7 = median(v3_7total)) %>% arrange(deadbody, desc(med7)) %>% pull(tissue_id)
  ggplot(fm_meta_dt, aes(x=tissue_id, y=n_muts, fill=group, color=group))+
    geom_boxplot()+
    scale_x_discrete(limits= x_order)+
    theme_syp+theme(axis.text.x = element_text(angle=45, hjust=1))
}

