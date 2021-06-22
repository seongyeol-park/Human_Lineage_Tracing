#check_bulk_sample_number_per_deadbody

library(tidyverse)
library(cowplot)

pmeta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_TargetSeq_201023.txt')
pmeta_dt %>% group_by(deadbody, Source_side) %>% dplyr::count()

merged_m_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Target_seq_related/merged_m_vaf_tbl.rds')
merged_m_vaf_tbl <- merged_m_vaf_tbl %>% mutate(Source_LR_CC = paste0(Source_side, '_', Source_side_CC)) 
merged_m_vaf_tbl %>% select(sample_id, deadbody, bait_region_mean_coverage) %>% unique() %>% group_by(deadbody) %>% dplyr::count()

colnames(merged_m_vaf_tbl)
merged_m_vaf_tbl$Source_side %>% unique()
merged_m_vaf_tbl$Source_side_CC %>% unique()
merged_m_vaf_tbl$anatomy_class1 %>% unique()
merged_m_vaf_tbl$dominant_layer2 %>% unique()

lr_pal = c("red","blue","gray","darkgray")
names(lr_pal) <- c("lt","rt","center","unknown")
cc_pal = c("pink","purple","darkgray")
names(cc_pal) <- ("cranio","caudal","unknown")
dl2_pal <- c("blue","yellow","red","lightgreen","orange")
names(dl2_pal) <- c("ectoderm","mesoderm","endoderm","ecto_mesoderm","meso_endoderm")

this_db='DB10'
tmp_tbl1 <- merged_m_vaf_tbl %>% filter(deadbody == this_db)%>% select(sample_id, bait_region_mean_coverage, Source_side, Source_side_CC, dominant_layer2) %>% unique()
g1 <- ggplot(tmp_tbl1, aes(x=reorder(sample_id,-bait_region_mean_coverage), y=bait_region_mean_coverage, fill = Source_side))+
  geom_bar(stat="identity")+
  geom_hline(yintercept=c(100,500), linetype='dashed')+
  scale_fill_manual(values = lr_pal)+
  theme(axis.text.x = element_blank(), axis.title = element_blank(), plot.title = element_text(size=18))+
  ggtitle(this_db)
g2 <- ggplot(tmp_tbl1, aes(x=reorder(sample_id,-bait_region_mean_coverage), y=bait_region_mean_coverage, fill = Source_side_CC))+
  geom_bar(stat="identity")+
  geom_hline(yintercept=c(100,500), linetype='dashed')+
  scale_fill_manual(values=cc_pal)+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
g3 <- ggplot(tmp_tbl1, aes(x=reorder(sample_id,-bait_region_mean_coverage), y=bait_region_mean_coverage, fill = dominant_layer2))+
  geom_bar(stat="identity")+
  geom_hline(yintercept=c(100,500), linetype='dashed')+
  scale_fill_manual(values=dl2_pal)+
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title= element_blank())
plot_grid(g1,g2,g3, ncol=1, align='v', rel_heights = c(1.2,1,1.7))

tmp_tbl1 %>% group_by(Source_side) %>% dplyr::count()
tmp_tbl1 %>% group_by(Source_side_CC) %>% dplyr::count()
tmp_tbl1 %>% group_by(dominant_layer2) %>% dplyr::count()




