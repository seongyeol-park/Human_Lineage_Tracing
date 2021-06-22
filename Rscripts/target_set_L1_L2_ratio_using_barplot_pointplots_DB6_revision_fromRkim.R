#load libraries
library(tidyverse)
library(ggtree)
library(ggsci)
library(scales)
library(cowplot)
library(ggrepel)


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
pmeta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_TargetSeq_210330.txt')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_210322.txt')
rate_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Mutation_rate_simulation/df_forest_plot.txt')
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_210330.txt') %>% filter(current_final_for_lineage == 'Y')
meta_dt <- meta_dt %>% mutate(Source3 = ifelse(is.na(Source_side)==T, Source2, paste0(Source_side,'_',Source2)))
meta_dt$deadbody <- factor(meta_dt$deadbody, levels =c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))

#Sourcing functions
source('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/TargetSeq_functions_afterTopUp.R')


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


#---------------------------------------------------------------script-------------------------------------------------------------------
#load data
f_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/f_vaf_tbl_rank.rds')

GSNP_mVAF <- f_vaf_tbl %>% filter(grepl('GSNP', lineage_id2)) %>% pull(VAF) %>% mean(., na.rm=T)


p_list <- list()
n=0

this_db <- "DB6"

n=n+1
m_vaf_tbl <- f_vaf_tbl %>% filter(deadbody == this_db & tissue_class != 'cancer_tissue')
fm_vaf_tbl <- m_vaf_tbl



# L1-2 vs L1-1
tmp_tbl <- fm_vaf_tbl %>% group_by(sample_id, lineage_id2) %>% summarise(meanVAF = mean(VAF)) 

tmp_tbl <- left_join(tmp_tbl,pmeta_dt %>% dplyr::select(sample_id, dominant_layer2, anatomy_class1)) %>% separate(lineage_id2, c('deadbody','lineage_id'), sep="_") %>% 
  spread(key=lineage_id, value=meanVAF) %>% dplyr::select(sample_id:anatomy_class1, `L1-1`) %>% 
  mutate(`L1-1`=`L1-1`*100, `L1-2`=0, `L1ratio`=0)

tmp_tbl %<>% mutate(dominant_layer2=replace(dominant_layer2, dominant_layer2!="ectoderm" & dominant_layer2!="endoderm" & dominant_layer2!="mesoderm", "mixed"))

g1 <- tmp_tbl %>% 
  ggplot(aes(x=sample_id, y=`L1-2`)) + 
  geom_bar(stat="identity", width=1.0, fill=rgb(249/255, 205/255, 45/255)) + 
  #geom_hline(yintercept=3.6/2, color="red", linetype="dashed")+
  coord_flip() + 
  scale_y_continuous(limits=c(0,100), expand=c(0,0), breaks=c(0, 50,100)) + 
  #scale_y_reverse(limits=c(100,0), expand=c(0,0), breaks=c(50,100)) + 
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(0,0.2,0,0), "cm"), axis.line.x=element_line(color="black"),
        text=element_text(size=25))

g2 <- tmp_tbl %>% 
  ggplot(aes(x=sample_id, y=`L1-1`)) + 
  geom_bar(stat="identity", width=1.0, fill=rgb(94/255, 133/255, 218/255)) + 
  geom_hline(yintercept=3.6/2, color="red", linetype="dashed")+
  coord_flip() + 
  scale_y_reverse(limits=c(100,0), expand=c(0,0), breaks=c(50,100)) + 
  #scale_y_continuous(limits=c(0,100), expand=c(0,0), breaks=c(0, 50,100)) + 
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(0,0,0,0.4), "cm"), axis.line.x=element_line(color="black"),
        text=element_text(size=25))

g3 <- tmp_tbl %>% 
  ggplot(aes(x=sample_id, y=L1ratio)) + 
  geom_point(color="black", size=1) + 
  #geom_hline(yintercept=exp_L1ratio, color="red", linetype="dashed")+
  coord_flip() + 
  #scale_y_reverse() + 
  scale_y_continuous(limits=c(0,100), expand=c(0,0)) + 
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(0,0.3,0,0.2), "cm"), panel.border = element_rect(color="black", fill=NA, size=1))

g4 <- tmp_tbl %>% 
  ggplot(aes(x=sample_id, y=factor(1), fill=anatomy_class1)) + geom_tile() + theme_void() + 
  scale_fill_manual(values=c("HN"=rgb(42/255,61/255,68/255),
                             "internal_organ"=rgb(214/255,124/255,52/255),
                             "LE"=rgb(40/255,142/255,203/255),
                             "trunk"=rgb(160/255,51/255,52/255),
                             "UE"=rgb(104/255,161/255,133/255))) + 
  scale_x_discrete(expand=c(0.01,0.01)) + 
  scale_y_discrete(expand=c(0,0)) + 
  coord_flip() +
  theme(legend.position="none", plot.margin=unit(c(0,0,0,0.1), unit="cm"))

g5 <- tmp_tbl %>% 
  ggplot(aes(x=sample_id, y=factor(1), fill=dominant_layer2)) + geom_tile() + theme_void() + 
  #scale_fill_igv() + 
  scale_fill_manual(values=c("ectoderm"=rgb(62/255,48/255,255/255),
                             "mesoderm"=rgb(193/255,40/255,38/255),
                             "endoderm"=rgb(98/255,140/255,70/255),
                             "mixed"=rgb(234/255, 225/255, 113/255))) +
  scale_x_discrete(expand=c(0.01,0.01)) + 
  scale_y_discrete(expand=c(0,0)) + 
  coord_flip() +
  theme(legend.position="none", plot.margin=unit(c(0,0,0,0.1), unit="cm"))

#plot_grid(g2, g1, g3, g4, g5, nrow=1, rel_widths = c(0.2, 0.2, 0.5, 0.05, 0.05), align = "h")
plot_grid(g2, g1, g4, g5, nrow=1, rel_widths = c(0.45, 0.45, 0.05, 0.05), align = "h")




# L1-2-2 vs L1-2-1

tmp_tbl <- fm_vaf_tbl %>% group_by(sample_id, lineage_id2) %>% summarise(meanVAF = mean(VAF)) 
tmp_tbl <- left_join(tmp_tbl,pmeta_dt %>% dplyr::select(sample_id, dominant_layer2, anatomy_class1)) %>% separate(lineage_id2, c('deadbody','lineage_id'), sep="_") %>% 
  spread(key=lineage_id, value=meanVAF) %>% 
  mutate(L1ratio=`L1-2-1`/(`L1-2-2`+`L1-2-1`), L2ratio=`L1-2-2`/(`L1-2-1`+`L1-2-2`)) %>% filter(is.na(L1ratio)==F & is.na(L2ratio) == F) %>%
  dplyr::select(sample_id:anatomy_class1, `L1-2-1`, `L1-2-2`, L1ratio, L2ratio) #%>%   gather(-c(sample_id,deadbody,dominant_layer, anatomy_class1), key=lineage_id, value=ratio)
tmp_tbl %<>% mutate(`L1-2-1`=`L1-2-1`*100, `L1-2-2`=`L1-2-2`*100, L1ratio=L1ratio*100, L2ratio=L2ratio*100)

tmp_tbl %<>% ungroup() %>% arrange(`L1ratio`) %>% mutate(sample_id=factor(sample_id, levels=sample_id))

tmp_tbl %<>% mutate(dominant_layer2=replace(dominant_layer2, dominant_layer2!="ectoderm" & dominant_layer2!="endoderm" & dominant_layer2!="mesoderm", "mixed"))


g1 <- tmp_tbl %>% 
  ggplot(aes(x=sample_id, y=`L1-2-1`)) + 
  geom_bar(stat="identity", width=1.0, fill=rgb(42/255,61/255,67/255)) + 
  geom_hline(yintercept=43.8/2, color="red", linetype="dashed")+
  coord_flip() + 
  scale_y_continuous(limits=c(0,100), expand=c(0,0), breaks=c(0, 50,100)) + 
  #scale_y_reverse(limits=c(100,0), expand=c(0,0), breaks=c(50,100)) + 
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(0,0.2,0,0), "cm"), axis.line.x=element_line(color="black"),
        text=element_text(size=25))

g2 <- tmp_tbl %>% 
  ggplot(aes(x=sample_id, y=`L1-2-2`)) + 
  geom_bar(stat="identity", width=1.0, fill=rgb(214/255,124/255,52/255)) +  
  geom_hline(yintercept=39.2/2, color="red", linetype="dashed")+
  coord_flip() + 
  scale_y_reverse(limits=c(100,0), expand=c(0,0), breaks=c(50,100)) + 
  #scale_y_continuous(limits=c(0,100), expand=c(0,0), breaks=c(0, 50,100)) + 
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(0,0,0,0.4), "cm"), axis.line.x=element_line(color="black"),
        text=element_text(size=25))

g3 <- tmp_tbl %>% 
  ggplot(aes(x=sample_id, y=L1ratio)) + 
  geom_point(color="black", size=1) + 
  geom_hline(yintercept=43.8/(43.8+39.2)*100, color="red", linetype="dashed")+
  coord_flip() + 
  #scale_y_reverse() + 
  scale_y_continuous(limits=c(0,100), expand=c(0,0)) + 
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(0,0.3,0,0.2), "cm"), panel.border = element_rect(color="black", fill=NA, size=1))

g4 <- tmp_tbl %>% 
  ggplot(aes(x=sample_id, y=factor(1), fill=anatomy_class1)) + geom_tile() + theme_void() + 
  scale_fill_manual(values=c("HN"=rgb(42/255,61/255,68/255),
                             "internal_organ"=rgb(214/255,124/255,52/255),
                             "LE"=rgb(40/255,142/255,203/255),
                             "trunk"=rgb(160/255,51/255,52/255),
                             "UE"=rgb(104/255,161/255,133/255))) + 
  scale_x_discrete(expand=c(0.01,0.01)) + 
  scale_y_discrete(expand=c(0,0)) + 
  coord_flip() +
  theme(legend.position="none", plot.margin=unit(c(0,0,0,0.1), unit="cm"))

g5 <- tmp_tbl %>% 
  ggplot(aes(x=sample_id, y=factor(1), fill=dominant_layer2)) + geom_tile() + theme_void() + 
  #scale_fill_igv() + 
  scale_fill_manual(values=c("ectoderm"=rgb(62/255,48/255,255/255),
                             "mesoderm"=rgb(193/255,40/255,38/255),
                             "endoderm"=rgb(98/255,140/255,70/255),
                             "mixed"=rgb(234/255, 225/255, 113/255))) +
  scale_x_discrete(expand=c(0.01,0.01)) + 
  scale_y_discrete(expand=c(0,0)) + 
  coord_flip() +
  theme(legend.position="none", plot.margin=unit(c(0,0,0,0.1), unit="cm"))

#plot_grid(g2, g1, g3, g4, g5, nrow=1, rel_widths = c(0.2, 0.2, 0.5, 0.05, 0.05), align = "h")
plot_grid(g2, g1, g4, g5, nrow=1, rel_widths = c(0.45, 0.45, 0.05, 0.05), align = "h")


# L1-2-1-2 vs L1-2-1-1
tmp_tbl <- fm_vaf_tbl %>% group_by(sample_id, lineage_id2) %>% summarise(meanVAF = mean(VAF)) 

tmp_tbl <- left_join(tmp_tbl,pmeta_dt %>% dplyr::select(sample_id, dominant_layer2, anatomy_class1)) %>% separate(lineage_id2, c('deadbody','lineage_id'), sep="_") %>% 
  spread(key=lineage_id, value=meanVAF) %>% dplyr::select(sample_id:anatomy_class1, `L1-2-1-1`) %>% 
  mutate(`L1-2-1-1`=`L1-2-1-1`*100, `L1-2-1-2`=0, `L1ratio`=0)

tmp_tbl %<>% ungroup() %>% arrange(desc(`L1ratio`)) %>% mutate(sample_id=factor(sample_id, levels=sample_id))

tmp_tbl %<>% mutate(dominant_layer2=replace(dominant_layer2, dominant_layer2!="ectoderm" & dominant_layer2!="endoderm" & dominant_layer2!="mesoderm", "mixed"))

g1 <- tmp_tbl %>% 
  ggplot(aes(x=sample_id, y=`L1-2-1-1`)) + 
  geom_bar(stat="identity", width=1.0, fill=rgb(40/255, 142/255, 203/255)) + 
  geom_hline(yintercept=34.0/2, color="red", linetype="dashed")+
  coord_flip() + 
  scale_y_continuous(limits=c(0,100), expand=c(0,0), breaks=c(0, 50,100)) + 
  #scale_y_reverse(limits=c(100,0), expand=c(0,0), breaks=c(50,100)) + 
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(0,0.2,0,0), "cm"), axis.line.x=element_line(color="black"),
        text=element_text(size=25))

g2 <- tmp_tbl %>% 
  ggplot(aes(x=sample_id, y=`L1-2-1-2`)) + 
  geom_bar(stat="identity", width=1.0, fill=rgb(40/255, 142/255, 203/255)) + 
  #geom_hline(yintercept=34.0/2, color="red", linetype="dashed")+
  coord_flip() + 
  scale_y_reverse(limits=c(100,0), expand=c(0,0), breaks=c(50,100)) + 
  #scale_y_continuous(limits=c(0,100), expand=c(0,0), breaks=c(0, 50,100)) + 
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(0,0,0,0.4), "cm"), axis.line.x=element_line(color="black"),
        text=element_text(size=25))

g3 <- tmp_tbl %>% 
  ggplot(aes(x=sample_id, y=L1ratio)) + 
  geom_point(color="black", size=1) + 
  #geom_hline(yintercept=exp_L1ratio, color="red", linetype="dashed")+
  coord_flip() + 
  #scale_y_reverse() + 
  scale_y_continuous(limits=c(0,100), expand=c(0,0)) + 
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid = element_blank(),
        plot.margin = unit(c(0,0.3,0,0.2), "cm"), panel.border = element_rect(color="black", fill=NA, size=1))

g4 <- tmp_tbl %>% 
  ggplot(aes(x=sample_id, y=factor(1), fill=anatomy_class1)) + geom_tile() + theme_void() + 
  scale_fill_manual(values=c("HN"=rgb(42/255,61/255,68/255),
                             "internal_organ"=rgb(214/255,124/255,52/255),
                             "LE"=rgb(40/255,142/255,203/255),
                             "trunk"=rgb(160/255,51/255,52/255),
                             "UE"=rgb(104/255,161/255,133/255))) + 
  scale_x_discrete(expand=c(0.01,0.01)) + 
  scale_y_discrete(expand=c(0,0)) + 
  coord_flip() +
  theme(legend.position="none", plot.margin=unit(c(0,0,0,0.1), unit="cm"))

g5 <- tmp_tbl %>% 
  ggplot(aes(x=sample_id, y=factor(1), fill=dominant_layer2)) + geom_tile() + theme_void() + 
  #scale_fill_igv() + 
  scale_fill_manual(values=c("ectoderm"=rgb(62/255,48/255,255/255),
                             "mesoderm"=rgb(193/255,40/255,38/255),
                             "endoderm"=rgb(98/255,140/255,70/255),
                             "mixed"=rgb(234/255, 225/255, 113/255))) +
  scale_x_discrete(expand=c(0.01,0.01)) + 
  scale_y_discrete(expand=c(0,0)) + 
  coord_flip() +
  theme(legend.position="none", plot.margin=unit(c(0,0,0,0.1), unit="cm"))

#plot_grid(g2, g1, g3, g4, g5, nrow=1, rel_widths = c(0.2, 0.2, 0.5, 0.05, 0.05), align = "h")
plot_grid(g2, g1, g4, g5, nrow=1, rel_widths = c(0.45, 0.45, 0.05, 0.05), align = "h")
