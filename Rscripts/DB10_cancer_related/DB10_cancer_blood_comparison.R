
library(tidyverse)


ch_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/clonal_hematopoiesis/CH_final_list.txt.anv.readc')
ch_dt <- ch_dt %>% filter(case_list == "10_Blood_3s_merged")
ch_dt <- ch_dt %>% mutate(vaf = var_readN/(var_readN + ref_readN))
ch_dt <- ch_dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_'))
ch_dt$mt_loca='blood'
colnames(ch_dt)[(ncol(ch_dt)-10):ncol(ch_dt)]
ch_dt <- ch_dt %>% separate("10_Breastcancer_WGS_ref;10_Breastcancer_WGS_var;10_Breastcancer_WGS_ukn;10_Breastcancer_WGS_vafpct", c("tumor_ref","tumor_var","tumor_ukn","tumor_vafpct"), sep=';', convert=T)
nrow(ch_dt)
ggplot(ch_dt)+
  geom_histogram(aes(x=tumor_vafpct))+
  ggtitle('Tumor VAFs of blood mutations (n=255)')



#cancer vcf
dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10_cancer_WGS/somatic_call_with_blood/04_mutect2_strelka2_union/tmp.txt.anv.readc.seqzcn.ccf',
               col_types=cols(`#CHROM`='c'))
colnames(dt)
dt <- dt %>% mutate(vaf= (var_readN/(var_readN + ref_readN))) %>%
  separate("10_Blood_3s_merged_ref;10_Blood_3s_merged_var;10_Blood_3s_merged_ukn;10_Blood_3s_merged_vafpct", c("blood_ref","blood_var","blood_ukn","blood_vaf"), sep=';', convert=T) 
dt$blood_vaf %>% summary()
dt <- dt %>% mutate(var_id = paste(`#CHROM`,POS,REF,ALT,sep='_'))
dt <- left_join(dt, ch_dt %>% select(var_id, mt_loca))
dt$mt_loca[is.na(dt$mt_loca)==T] <- 'tumor_only'
dt$mt_loca %>% table()

tmp1 <- dt %>% select(`#CHROM`, POS, vaf, ref_readN, var_readN, blood_vaf, mt_loca)
tmp2 <- ch_dt %>% mutate(blood_vaf=vaf*100, vaf=tumor_vafpct*0.01, ref_readN = tumor_ref, var_readN = tumor_var) %>% select(`#CHROM`, POS, ref_readN, var_readN,vaf, blood_vaf, mt_loca)
tmp <- bind_rows(tmp1, tmp2) %>% unique()
tmp$vaf %>% summary()
tmp$blood_vaf %>% summary()
tmp %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10_cancer_WGS/somatic_call_with_blood/04_mutect2_strelka2_union/intersection_plus_blooddetected.txt')

#cancer + blood vcf
merge_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10_cancer_WGS/somatic_call_with_blood/04_mutect2_strelka2_union/intersection_plus_blooddetected.txt.seqzcn.ccf',
                     col_types =cols(`#CHROM`='c', CCF='n'))

merge_dt <- merge_dt %>% mutate(group = ifelse(blood_vaf >=1, 'pos','neg'))

t_chr='1'
plot_list=list()
n=0
for(t_chr in c(1:22,'X')){
  n=n+1
plot_list[[n]] <- ggplot(subset(merge_dt, `#CHROM`==t_chr), aes(x=POS, y=mutCN))+
  geom_point(alpha=0.6, aes(color=group))+
  geom_segment(data = subset(seg_dt, chromosome == t_chr), aes(x=start.pos, xend=end.pos, y=A, yend=A), color="red")+
  geom_segment(data = subset(seg_dt, chromosome == t_chr), aes(x=start.pos, xend=end.pos, y=B, yend=B), color="blue")+
  scale_y_continuous(limits=c(0,10), breaks=seq(0,10,1))+
  scale_color_manual(values=c('pos'='red', 'neg'='gray'))+
  #scale_color_manual(values=c('blood'='red','tumor_only'='gray'))+
  #scale_color_gradient2(low='gray', mid='gray',high='red', midpoint=3)+
  theme(legend.position = 'none')+
  xlab(t_chr)
}
library(cowplot)
plot_grid(plotlist = plot_list[1:6], nrow=1, align="h", axis="tb")


merge_dt %>% filter(`#CHROM`=='4' & mutCN >5 & blood_vaf < 5)


ggplot(merge_dt, aes())


##########################

ggplot(dt, aes(x=vaf, y=blood_vaf))+
  geom_point()



#
seg_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10_cancer_WGS/somatic_call_with_blood/01_sequenza/10_Breastcancer_WGS_sequenza.seqz.v3/10_Breastcancer_WGS_sequenza_Sol1_c0.48_p3_gF/10_Breastcancer_WGS_sequenza_Sol1_c0.48_p3_gF_segments.txt',
                   col_types=cols(chromosome='c'))

t_chr='1'
ggplot(subset(dt, `#CHROM`==t_chr), aes(x=POS, y=vaf))+
  geom_point(alpha=0.8)



  geom_point(data=subset(dt, `#CHROM`==t_chr), aes(x=POS, y=vaf*10+10, color=mt_loca), alpha=0.3)+
  geom_segment(data = subset(seg_dt, chromosome == t_chr), aes(x=start.pos, xend=end.pos, y=CNt, yend=CNt), color="red")+
  scale_y_continuous(limits=c(0,20), breaks=seq(0,20,2))+
  xlab(t_chr)+theme(axis.title.y=element_blank())


plot_list=list()
n=0
for(t_chr in c(1:22,'X')){
  n=n+1
  plot_list[[n]] <- ggplot(subset(ch_dt, `#CHROM`==t_chr), aes(x=POS, y=vaf*10))+
    geom_point(alpha=0.8)+
    geom_point(data=subset(dt, `#CHROM`==t_chr), aes(x=POS, y=vaf*10+10, color=mt_loca), alpha=0.3)+
    geom_segment(data = subset(seg_dt, chromosome == t_chr), aes(x=start.pos, xend=end.pos, y=CNt, yend=CNt), color="red")+
    scale_y_continuous(limits=c(0,20), breaks=seq(0,20,2))+
    xlab(t_chr)+theme(axis.title.y=element_blank())+theme(legend.position='none')
  
}
library(cowplot)
plot_grid(plotlist = plot_list[1:6], nrow=1, align="h", axis="tb")







plot_list=list()
n=0
for(t_chr in c(1:22,'X')){
  n=n+1
  plot_list[[n]] <- ggplot(subset(ch_dt, `#CHROM`==t_chr), aes(x=POS, y=vaf*10))+
    geom_point(alpha=0.8)+
    geom_point(data=subset(dt, `#CHROM`==t_chr), aes(x=POS, y=vaf*10+10, color=mt_loca), alpha=0.3)+
    geom_segment(data = subset(seg_dt, chromosome == t_chr), aes(x=start.pos, xend=end.pos, y=CNt, yend=CNt), color="red")+
    scale_y_continuous(limits=c(0,20), breaks=seq(0,20,2))+
    xlab(t_chr)+theme(axis.title.y=element_blank())
  
}
library(cowplot)
plot_grid(plotlist = plot_list[21:23], nrow=1, align="h", axis="tb")

