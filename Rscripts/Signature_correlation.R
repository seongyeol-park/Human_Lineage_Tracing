#signature correlation
library(tidyverse)
library(PerformanceAnalytics)
library(ggsci)
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_200706.txt') %>% filter(current_final_for_lineage == 'Y')
dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/sig_SBS15187abcd_a4/merged_sig_result.txt')

db_pal <- pal_npg('nrc')(10)[c(1:5,9,7)]
names(db_pal) <- unique(meta_dt$deadbody)

theme_syp = theme(panel.border = element_blank(), axis.line = element_line(), axis.text = element_text(size=15), axis.title = element_text(size = 18))


# proportion table
mx <- dt %>% mutate(sample_id= gsub('.snv','',`#sample_id`), v3_7total = v3_7a + v3_7b + v3_7c + v3_7d) %>% select(sample_id, v3_1, v3_5, v3_7total, v3_18) %>% as.data.frame() %>% column_to_rownames('sample_id') %>% as.matrix()
chart.Correlation(mx)
# count table
mx2 <- dt %>% mutate_at(vars(starts_with('v3')), list(~ .*n_subs*0.01)) %>% mutate(sample_id= gsub('.snv','',`#sample_id`), v3_7total = v3_7a + v3_7b + v3_7c + v3_7d) %>% select(sample_id, v3_1, v3_5, v3_7total, v3_18) %>% as.data.frame() %>% column_to_rownames('sample_id') %>% as.matrix()
chart.Correlation(mx2)

# Count table linear correlation
m_meta_dt <- left_join(meta_dt, dt %>% mutate(sample_id= gsub('.snv','',`#sample_id`), CS_v3 = CS) %>% select(sample_id, starts_with('v3_'), CS_v3))
m_meta_dt <- m_meta_dt %>% mutate_at(vars(starts_with('v3_')), list(~ .*n_snv*0.01))
m_meta_dt <- m_meta_dt %>% mutate(v3_7total = v3_7a + v3_7b + v3_7c + v3_7d)              
#Correlation between v3_1 and v3_7
ggplot(m_meta_dt, aes(x=v3_7total, y =v3_1, color = deadbody))+
  geom_point(alpha=0.8)+  
  scale_color_manual(values = db_pal)+
  ggtitle("Total samples")+
  theme_bw()+
  theme_syp
ggplot(subset(m_meta_dt, v3_7total > 100 & v3_1 > 0 & Cell_type == 'skin_fb' & deadbody %in% c('DB6','DB8','DB9','DB10')), aes(x=v3_7total, y =v3_1, color = deadbody))+
  geom_point()+  
  geom_smooth(method="lm")+
  geom_abline(slope=-0.01, intercept=110, linetype='dashed')+
  scale_color_manual(values=db_pal)+
  ggtitle("Subset")+
  theme_bw() + theme_syp
#Correlation between v3_5 and v3_7
ggplot(m_meta_dt, aes(x=v3_7total, y =v3_5, color = deadbody))+
  geom_point()+
  scale_color_manual(values = db_pal)+
  ggtitle("Total samples")+
  theme_bw()+
  theme_syp
ggplot(subset(m_meta_dt, v3_7total > 100 & Cell_type == 'skin_fb' & deadbody %in% c('DB6','DB8','DB9','DB10')), aes(x=v3_7total, y =v3_5, color = deadbody))+
  geom_point()+
  geom_smooth(method="lm")+
  geom_abline(slope=0.1, intercept=1500, linetype='dashed')+
  scale_color_manual(values=db_pal)+
  ggtitle("Subset")+
  theme_bw() + theme_syp
#Correlation between v3_18 and v3_7
ggplot(m_meta_dt, aes(x=v3_7total, y =v3_18, color = deadbody))+
  geom_point()+
  scale_color_manual(values = db_pal)+
  ggtitle("Total samples")+
  theme_bw()+
  theme_syp
ggplot(subset(m_meta_dt, v3_7total > 100 & Cell_type == 'skin_fb' & deadbody %in% c('DB6','DB8','DB9','DB10')), aes(x=v3_7total, y =v3_18, color = deadbody))+
  geom_point()+
  geom_smooth(method="lm")+
  geom_abline(slope=-0.01, intercept=250, linetype='dashed')+
  scale_color_manual(values=db_pal)+
  ggtitle("Subset")+
  theme_bw() + theme_syp


#add correlated count
m_meta_dt <- m_meta_dt %>% mutate(v3_1cor = v3_1 + v3_7total*0.01, 
                                  v3_5cor = v3_5 - v3_7total*0.1, 
                                  v3_18cor = v3_18 + v3_7total*0.01,
                                  v3_endo = v3_1cor + v3_5cor + v3_18cor)

#age and totalSNV
ggplot(subset(m_meta_dt, Cell_type == 'skin_fb'), aes(x=age, y=n_snv))+
  geom_boxplot(aes(fill = deadbody), outlier.size=-1, alpha=0.5, width=2)+
  geom_jitter(aes(color = deadbody), alpha=0.3, size=1, height=0, width=1)+
  geom_smooth()+
  scale_y_continuous(limits=c(0,30000))+
  scale_x_continuous(limits=c(0,100))+
  scale_fill_manual(values = db_pal)+
  scale_color_manual(values = db_pal)+
  ggtitle("total SNV")+
  theme_bw() + theme_syp

m_meta_dt %>% filter(Cell_type == 'skin_fb') %>% .$age -> age
m_meta_dt %>% filter(Cell_type == 'skin_fb') %>% .$n_snv -> mutation
summary(lm(mutation ~ age))

#age and v3_1 + v3_5 + v3_18
ggplot(subset(m_meta_dt, Cell_type == 'skin_fb'), aes(x=age, y=v3_5+v3_1+v3_18))+
  geom_boxplot(aes(fill = deadbody), outlier.size=-1, alpha=0.5, width=2)+
  geom_jitter(aes(color = deadbody), alpha=0.3, size=1, height=0, width=1)+
  geom_smooth()+
  scale_y_continuous(limits=c(0,4500))+
  scale_x_continuous(limits=c(0,100))+
  scale_fill_manual(values = db_pal)+
  scale_color_manual(values = db_pal)+
  ggtitle("Raw SBS1+5+18")+
  theme_bw() + theme_syp

m_meta_dt %>% filter(Cell_type == 'skin_fb') %>% .$age -> age
m_meta_dt %>% filter(Cell_type == 'skin_fb') %>% mutate(endo_mt = v3_1 + v3_5 + v3_18) %>% .$endo_mt -> mutation
summary(lm(mutation ~ age))

#age and corrected(v3_1 + v3_5 + v3_18)
ggplot(subset(m_meta_dt, Cell_type == 'skin_fb'), aes(x=age, y=v3_5cor+v3_1cor+v3_18cor))+
  geom_boxplot(aes(fill = deadbody), outlier.size=-1, alpha=0.5, width=2)+
  geom_jitter(aes(color = deadbody), alpha=0.3, size=1, height=0, width=1)+
  geom_smooth()+
  scale_y_continuous(limits=c(0,4500))+
  scale_x_continuous(limits=c(0,100))+
  scale_fill_manual(values = db_pal)+
  scale_color_manual(values = db_pal)+
  ggtitle("Corrected SBS1+5+18")+
  theme_bw() + theme_syp

m_meta_dt %>% filter(Cell_type == 'skin_fb') %>% .$age -> age
m_meta_dt %>% filter(Cell_type == 'skin_fb') %>% mutate(cor_endo_mt = v3_1cor + v3_5cor + v3_18cor) %>% .$cor_endo_mt -> mutation
fit <- lm(mutation ~ age)
summary(fit)

#high VAF point mutations
file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/', pattern = '.seqzcn$', full.names = T)
res_tbl <- tibble(sample_id = as.character(), n_pointmt_VAF1sLOH = as.integer())
n=0
for (fn in file_list){
  n = n+1
  sample_id <- unlist(strsplit(unlist(strsplit(fn,'//'))[2],'\\.'))[1]
  print(sample_id)
  dt <- read_tsv(fn, comment='##', col_types = cols(`#CHROM`='c'))
  colnames(dt)[10] <- 'info'
  this_n <- dt %>% separate(info, c('RD','AD','AF'), sep=':', convert = T) %>% filter(minCN > 0 & AF > 95) %>% nrow()
  res_tbl[n,] <- list(sample_id, this_n)
}
m_meta_dt <- left_join(m_meta_dt, res_tbl)
ggplot(m_meta_dt, aes(x=v3_7total, y=n_pointmt_VAF1sLOH))+
  geom_point(aes(color=deadbody),size=3, alpha=0.7)+
  scale_color_manual(values = db_pal)+
  theme_bw()+theme_syp

#IGV inspection
dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/3_RK_001B4.edited.vcf.anv.seqzcn', comment = '##', col_types = cols(`#CHROM`='c'))
colnames(dt)[10] <- 'info'
dt %>% separate(info, c('RD','AD','AF'), sep=':', convert = T) %>% filter(minCN > 0 & AF > 95) %>%  View()



#nSV and UV
DB2_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/db2/09_delly/02_merged_annotated/',pattern ='igv_check$', full.names = T)
DB3_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/db3/10_Delly/02_somatic_call/02_merged_call/', pattern = 'igv_check$', full.names = T)
DB5_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/db5/10_delly/02_merged_annotated/',pattern = 'igv_check$', full.names = T)
DB6_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/db6/10_Delly/02_merged/', pattern = 'igv_check$', full.names = T)
DB8_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/db8/08_delly/02_merged_annotated/', pattern = 'igv_check$', full.names = T)
DB9_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/db9/08_delly/02_merged_annotated/', pattern = 'igv_check$', full.names = T)
DB10_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/db10/08_delly/02_merged_annotated/', pattern = 'igv_check$', full.names = T)

file_list <- c(DB2_list, DB3_list, DB5_list, DB6_list, DB8_list, DB9_list, DB10_list)
res_tbl <- tibble(sample_id = as.character(), n_SV = as.integer(), n_DEL = as.integer())
NA_rc_v = vector()
for (n in 1:length(file_list)){
  sample_id <- unlist(strsplit(unlist(strsplit(file_list[n],'//'))[2],'\\.'))[1]
  print(sample_id)
  dt <- read_tsv(file_list[n], col_types= cols(`#CHR1`='c', CHR2='c'))
  n_SV <- nrow(dt)
  n_DEL <- dt %>% filter(svtype == 'DEL') %>% nrow()
  res_tbl[n,] <- list(sample_id, n_SV, n_DEL)
  
}
res_tbl
m_meta_dt <- left_join(m_meta_dt, res_tbl)
ggplot(m_meta_dt, aes(x=v3_7total, y=n_DEL))+
  geom_point(aes(color = deadbody), size=3, alpha=0.5)+
  theme_bw()+theme_syp

#CNV and UV
m_meta_dt <- m_meta_dt %>% mutate(CNV_presence = ifelse(CNV == 'none', 'N', 'Y'))
ggplot(m_meta_dt, aes(x=CNV_presence, y=v3_7total))+
  geom_boxplot(outlier.size =-1)+
  geom_jitter(height=0, width=0.1, alpha=0.3)+
  theme_bw() + theme_syp

ggplot(m_meta_dt, aes(x=v3_7total, y=v3_5))+
  geom_point(aes(color = CNV_presence), alpha=0.5)+
  theme_bw() + theme_syp

#telomer legnth and UV
ggplot(m_meta_dt, aes(x=v3_7total, y=log10(telomere_length), color=deadbody))+
  geom_point( size=3, alpha=0.5)+
  geom_smooth(method = "lm")+
  scale_color_manual(values = db_pal)+
  theme_bw() + theme_syp
  
