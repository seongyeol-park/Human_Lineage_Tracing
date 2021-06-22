#Deprecated scripts












####<UV signature >>>###########
Draw_totalTree_Sig7prop(meta_dt, lng_dt, nwk_tbl)
colfunc <- colorRampPalette(c("yellow", "blue"))
legend_image <- as.raster(matrix(colfunc(20), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
text(x=0.5, y = seq(0,1,0.25), labels = paste0(seq(0,100,25),'%'), cex=2.5)
rasterImage(legend_image, 0, 0, 0.3,1)


####Signature analysis
#signature multicorrelation
sig_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/sig_SBS15187abcd_a1/merged_sig_result_v2.txt')
#sig_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/sig_SBS15187abcd_a4_v401/merged_sig_result_v2.txt')
sig_dt <- sig_dt %>% rename(sample_id = `#sample_id`) %>% mutate_at(vars(starts_with('v3')), list(~ .*n_subs*0.01)) %>% mutate(deadbody = apply(.["sample_id"], 1, function(x) paste0('DB',unlist(strsplit(x,'_'))[1])))
mx <- sig_dt  %>% select(-n_subs, -CS, -deadbody) %>% as.data.frame() %>% column_to_rownames('sample_id') %>% as.matrix()
library(PerformanceAnalytics)
chart.Correlation(mx)
sig_dt
ggplot(sig_dt, aes(x=v3_7a, y=v3_5, color=deadbody))+
  geom_point()
m_sig_dt <- sig_dt %>% mutate_at(vars(starts_with('v3')), list(~ .*n_subs*0.01)) %>% gather(-sample_id, -n_subs, -CS, key='signature', value='count')
x_order <- sig_dt %>% arrange(desc(n_subs)) %>% pull(sample_id)
m_sig_dt
ggplot(m_sig_dt, aes(x=sample_id, y=count))+
  geom_bar(aes(fill=signature), stat="identity")+
  scale_x_discrete(limits = x_order)


#all deadbody Sig5 and sig7 correlation 
if(F){
  meta_dt %>% select(deadbody, age) %>% unique()
  meta_dt <- meta_dt %>% mutate(Sig7_ct = n_snv * Sig7_prop, Sig5_ct = n_snv * Sig5_prop)
  ggplot(subset(meta_dt, Cell_type == 'skin_fb'), aes(x=Sig7_ct, y=Sig5_ct))+
    geom_point(aes(color=deadbody))+
    scale_color_manual(values = db_pal)+
    theme_bw()+xlab("No. of Sig7 substitution")+ ylab("No. of Sig5 substitution")+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=15))
  
}

#Calculate linear correlation of each deadbody
meta_dt %>% filter(deadbody == 'DB10') %>% .$Sig7_ct -> a
meta_dt %>% filter(deadbody == 'DB10') %>% .$Sig5_ct -> b
lm(b ~ a)
ggplot(subset(meta_dt, deadbody == 'DB10'), aes(x=Sig7_ct, y=Sig5_ct))+
  geom_point(aes(color= Source2))+
  geom_abline(slope = 0.29, intercept = 1500, )
ggplot(subset(meta_dt, deadbody == 'DB8'), aes(x=Sig7_ct, y=Sig5_ct))+
  geom_point(aes(color= Source2), size=5, alpha=0.7)+
  geom_abline(slope = 0.29, intercept = 3100, linetype ='dashed' )+
  geom_abline(slope = 0.29, intercept = 2100, linetype ='dashed' )
ggplot(subset(m_meta_dt, deadbody == 'DB9'), aes(x=Sig7_ct, y=Sig5_ct))+
  geom_point(aes(color= Source3), size=5, alpha=0.7)+
  geom_abline(slope = 0.29, intercept = 3300, linetype ='dashed' )+
  geom_abline(slope = 0.29, intercept = 2100, linetype ='dashed' )+
  ggtitle("DB9")
ggplot(subset(meta_dt, deadbody == 'DB6'), aes(x=Sig7_ct, y=Sig5_ct))+
  geom_point(aes(color= Cell_type), size=5, alpha=0.7)+
  geom_abline(slope=0.29, intercept=2000, linetype="dashed")+
  geom_abline(slope=0.7, intercept=2000, linetype="dotted")

#exceptional two samples ("6_MS-11_001B3", "6_LMCS-2_001F6") in sig5-7 correlation
m_meta_dt <- meta_dt %>% mutate(Sig1_ct = n_snv * Sig1_prop,
                                Sig5_ct = n_snv * Sig5_prop,
                                Sig7_ct = n_snv * Sig7_prop,
                                corSig5_ct = Sig5_ct - 0.29*(Sig7_ct))

ggplot(subset(m_meta_dt, Cell_type != 'skin_fb'), aes(x=Sig7_ct, y=Sig5_ct))+
  geom_point(aes(shape=Cell_type, color=deadbody))+
  geom_text_repel(data=subset(m_meta_dt, Cell_type != 'skin_fb' & Sig5_ct > 5800), aes(label=sample_id))


#UV signature per anatomical site
if(F){
  m_meta_dt <- meta_dt %>% mutate(Sig1_ct = n_snv * Sig1_prop,
                                  Sig5_ct = n_snv * Sig5_prop,
                                  Sig7_ct = n_snv * Sig7_prop,
                                  corSig5_ct = Sig5_ct - 0.29*(Sig7_ct))
  
  s_meta_dt <- m_meta_dt %>% filter(Cell_type %in% c('skin_fb'))
  s_meta_dt$deadbody <- factor(s_meta_dt$deadbody, levels = c('DB3','DB5','DB10','DB6','DB9','DB2','DB8'))
  s_meta_dt %>% group_by(deadbody, tissue_id) %>% summarise(numb=n(), median=median(Sig7_ct)) %>% arrange(deadbody, desc(median)) %>% .$tissue_id -> x_order
  
  ggplot(s_meta_dt, aes(x=tissue_id, y=Sig7_ct))+
    geom_boxplot(outlier.size=-1)+
    geom_point(aes(color=deadbody))+
    scale_x_discrete(limits=x_order)+
    scale_color_manual(values = db_pal)+
    #scale_y_continuous(limits=c(0,4700))+
    xlab("Skin tissues")+ylab("No. of UV-associated mutations")+
    theme_bw()+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())
  
}


####anatomical poisition and endogenous mutation rate variation
if(F){
  m_meta_dt <- meta_dt %>% mutate(Sig1_ct = n_snv * Sig1_prop,
                                  Sig5_ct = n_snv * Sig5_prop,
                                  Sig7_ct = n_snv * Sig7_prop,
                                  corSig5_ct = Sig5_ct - 0.29*(Sig7_ct))
  
  s_meta_dt <- m_meta_dt %>% filter(Cell_type %in% c('skin_fb'))
  s_meta_dt$deadbody <- factor(s_meta_dt$deadbody, levels = c('DB3','DB5','DB10','DB6','DB9','DB2','DB8'))
  s_meta_dt %>% group_by(deadbody, tissue_id) %>% summarise(numb=n(), median=median(Sig1_ct + corSig5_ct)) %>% arrange(deadbody, desc(median)) %>% .$tissue_id -> x_order
  
  ggplot(s_meta_dt, aes(x=tissue_id, y=Sig1_ct + corSig5_ct))+
    geom_boxplot(outlier.size=-1)+
    geom_point(aes(color=deadbody))+
    scale_x_discrete(limits=x_order)+
    scale_color_manual(values = db_pal)+
    scale_y_continuous(limits=c(0,4700))+
    xlab("Skin tissues")+ylab("No. of endogenous mutations")+
    theme_bw()+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())
}


#anatomical position (proximal to distal) and endogenous mutation rate (Deprecated)
meta_dt %>% group_by(deadbody) %>% summarise(source2_n = length(unique(Source2)))
meta_dt <- meta_dt %>% mutate(Sig1_ct = n_snv * Sig1_prop,
                              Sig5_ct = n_snv * Sig5_prop,
                              Sig7_ct = n_snv * Sig7_prop,
                              corSig5_ct = Sig5_ct - 0.29*(Sig7_ct))
s_meta_dt <- meta_dt %>% filter(deadbody == 'DB6' & Cell_type=='skin_fb')
gsub("lower_arm","arm",gsub("upper_arm","arm",s_meta_dt$Source2))
ggplot(s_meta_dt, aes(x=Source2, y=Sig1_ct + corSig5_ct))+
  geom_boxplot()+
  geom_jitter(height=0, width=0.1)+
  theme(axis.text.x= element_text(angle=45, hjust=1))
target_sources= c('thigh','leg','foot','upper_arm','lower_arm','arm','hand','deltoid','elbow')
p <- list()
this_db = 'DB6'
s_meta_dt <- meta_dt %>% filter(deadbody == this_db & Cell_type=='skin_fb' & Source2 %in% target_sources)
s_meta_dt$Source3 <- gsub('deltoid','shoulder',gsub('elbow','arm',gsub("lower_arm","arm",gsub("upper_arm","arm",s_meta_dt$Source2))))
s_meta_dt$Source3 = factor(s_meta_dt$Source3, levels = c('thigh','leg','foot','shoulder','arm','hand'))
p[[1]] <- ggplot(s_meta_dt, aes(x=Source3, y=Sig1_ct + corSig5_ct))+
  geom_boxplot(aes(fill=Source_class), outlier.size=-1)+
  geom_jitter(height=0, width=0.1)+
  theme_bw()+
  ylab("No. of endogenous mutations")+
  theme(axis.text.x= element_blank(), legend.position = "none", axis.title.x=element_blank(),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))+
  ggtitle(this_db)
p[[2]] <- ggplot(s_meta_dt, aes(x=Source3, y=telomere_length))+
  geom_boxplot(aes(fill=Source_class), outlier.size=-1)+
  geom_jitter(height=0, width=0.1)+
  theme_bw()+
  ylab("Score of telomere length")+
  theme(axis.text.x= element_text(angle=45, hjust=1), legend.position = "none", axis.title.x=element_blank(), plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))
this_db= 'DB8'
s_meta_dt <- meta_dt %>% filter(deadbody == this_db & Cell_type=='skin_fb' & Source2 %in% target_sources)
s_meta_dt$Source3 <- gsub('deltoid','shoulder',gsub('elbow','arm',gsub("lower_arm","arm",gsub("upper_arm","arm",s_meta_dt$Source2))))
s_meta_dt$Source3 = factor(s_meta_dt$Source3, levels = c('thigh','leg','foot','shoulder','arm','hand'))
p[[3]] <- ggplot(s_meta_dt, aes(x=Source3, y=Sig1_ct + corSig5_ct))+
  geom_boxplot(aes(fill=Source_class), outlier.size=-1)+
  geom_jitter(height=0, width=0.1)+
  theme_bw()+
  theme(axis.text.x= element_blank(), legend.position = "none", axis.title=element_blank(),plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))+
  ggtitle(this_db)
p[[4]] <- ggplot(s_meta_dt, aes(x=Source3, y=telomere_length))+
  geom_boxplot(aes(fill=Source_class), outlier.size=-1)+
  geom_jitter(height=0, width=0.1)+
  theme_bw()+
  theme(axis.text.x= element_text(angle=45, hjust=1), legend.position = "none", axis.title=element_blank(), plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"))
plot_grid(plotlist=p[c(1,3,2,4)], nrow=2, align="hv", axis="lrbt")



#CNV variation
db='DB3'
tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% mutate(l_id = taxa) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
m_lng_dt <- lng_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id) %>% select(-deadbody)
m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
max_tbl <- lng_dt %>% group_by(deadbody) %>% summarise(max_value=max(end_mtn)+500)
max_tbl$max_value[max_tbl$deadbody=='DB3'] <- 7500
node_list <- m_lng_dt %>% filter(is.na(Cell_type)== T) %>% .$taxa
for (node_id in node_list){
  print(node_id)
  for(col_id in c('CNV')){
    print(col_id)
    ex_values <- m_lng_dt %>% filter(grepl(node_id, taxa)==T & is.na(get(col_id)) == F) %>% .[[col_id]] %>% unique()
    if (length(ex_values)==1){
      m_lng_dt[m_lng_dt$taxa == node_id,][[col_id]] <- ex_values
      m_lng_dt[grepl(paste0(node_id,'-'),m_lng_dt$taxa)==T,][[col_id]] <- 'none'
    }
  }
}
m_lng_dt$CNV[m_lng_dt$CNV=='none'] <- NA
m_lng_dt %>% filter(start_mtn < 300) %>% nrow() -> a
m_lng_dt %>% filter(start_mtn < 300 & is.na(CNV) ==F ) %>% nrow() -> b
b*100/a
#fill data of nodes if information of it's branches is unique
ggtree(tree) %<+% m_lng_dt + geom_tree()+theme_tree2() +
  geom_label(aes(label=CNV),hjust=1.1, alpha=0.9, size=4)+
  #scale_fill_gradient2(low= "gray50", mid = "yellow", high="red", limits=c(0,60), midpoint=30)+
  coord_cartesian(xlim = c(-2,30))+
  scale_x_continuous(breaks=seq(0,30,5))+
  ggtitle(db)+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20), plot.title = element_text(size=30), plot.margin= unit(c(0.2,0.2,0.2,0.2),"cm"))
ggtree(tree) %<+% m_lng_dt + geom_tree()+theme_tree2() +
  geom_label(aes(label=CNV),hjust=0, alpha=0.8, size=5, fill ="yellow")+
  geom_tiplab(aes(label=Source3, color = Cell_type), hjust=0)+
  #scale_color_manual(values=cell_pal)+
  #scale_fill_gradient2(low= "gray50", mid = "yellow", high="red", limits=c(0,60), midpoint=30)+
  coord_cartesian(xlim = c(1,max_tbl$max_value[max_tbl$deadbody == db]))+
  #coord_cartesian(xlim=c(1,15000))+
  #scale_x_continuous(breaks=seq(15000,25000,5000))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20), plot.title = element_text(size=30), plot.margin= unit(c(0.2,0.2,0.2,0.2),"cm"))

#####UV damage 
db='DB6'
tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% mutate(l_id = taxa) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
m_lng_dt <-lng_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id) %>% select(-deadbody) %>% mutate(CCTTpct = CCTTcount*200/n_subs)
m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
max_tbl <- lng_dt %>% group_by(deadbody) %>% summarise(max_value=max(end_mtn)+500)
max_tbl$max_value[max_tbl$deadbody=='DB3'] <- 7500
#fill data of nodes if information of it's branches is unique
ggtree(tree) %<+% m_lng_dt + geom_tree()+theme_tree2() +
  geom_label(aes(label=round(CCTTpct,2), fill=CCTTpct),hjust=1.1, alpha=0.8, size=5)+
  scale_fill_gradient2(low= "gray50", mid = "yellow", high="red", limits=c(0,60), midpoint=30)+
  coord_cartesian(xlim = c(-2,30))+
  scale_x_continuous(breaks=seq(0,30,5))+
  ggtitle(db)+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20), plot.title = element_text(size=30), plot.margin= unit(c(0.2,0.2,0.2,0.2),"cm"))
ggtree(tree) %<+% m_lng_dt + geom_tree()+theme_tree2() +
  geom_label(aes(label=round(CCTTpct,2), fill=CCTTpct),hjust=1, alpha=0.8, size=5)+
  geom_tiplab(aes(label=Source3, color = Cell_type), hjust=0)+
  scale_color_manual(values=cell_pal)+
  scale_fill_gradient2(low= "gray50", mid = "yellow", high="red", limits=c(0,60), midpoint=30)+
  coord_cartesian(xlim = c(0,max_tbl$max_value[max_tbl$deadbody == db]))+
  #coord_cartesian(xlim=c(11000,28000))+
  #scale_x_continuous(breaks=seq(15000,25000,5000))+
  theme(axis.title = element_text(size=20),axis.text = element_text(size=20), plot.title = element_text(size=30), plot.margin= unit(c(0.2,0.2,0.2,0.2),"cm"))
colfunc <- colorRampPalette(c("yellow", "gray50"))
plot(1:10, 1:10, pch = 19, cex=2, col = colfunc(15))
legend_image <- as.raster(matrix(colfunc(15), ncol=1))
plot(c(0,3),c(0,10),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
text(x=0.15, y = seq(0,3,1), labels = seq(0,30,10), cex=1.5)
rasterImage(legend_image, 0, 0, 0.1,3)
colfunc <- colorRampPalette(c("red", "yellow"))
legend_image <- as.raster(matrix(colfunc(15), ncol=1))
text(x=0.15, y = seq(3,6,1), labels = seq(30,60,10), cex=1.5)
rasterImage(legend_image, 0, 3, 0.1,6)
plot(x=1)
legend("topleft", legend=names(cell_pal), fill=cell_pal,
       bty="n", border="white", cex=1, ncol=6, text.width = c(0.05,0.05, 0.05, 0.05, 0.06, 0.07), x.intersp = 0.2 )

#UV per anatomical position
sub_dt <- meta_dt %>% filter(deadbody == 'DB6' & current_final_for_lineage == 'Y')
sample_order <- sub_dt %>% group_by(Source3) %>% summarise(median_sig7_ct = median(Sig7_ct)) %>% arrange(desc(median_sig7_ct)) %>% .$Source3
ggplot(sub_dt, aes(x=Source3, y=Sig7_ct))+
  geom_jitter(aes(color=Source_class),alpha=0.6, size=3, height=0, width=0.2)+
  scale_x_discrete(limits = sample_order)+
  xlab('Samples sorted by median value')+ylab("No. of Sig7 subs")+
  ggtitle("DB6")+
  theme_bw(base_size=15)+ theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))


#####Age correlation without correction
ggplot(subset(meta_dt, Cell_type == 'skin_fb'), aes(x=age, y=Sig1_ct+Sig5_ct))+
  geom_point(aes(color=deadbody), size=3, alpha=0.5)+
  geom_smooth(method="lm")+
  geom_boxplot(aes(fill=deadbody),alpha=0.5)+
  scale_color_manual(values=db_pal)+
  scale_fill_manual(values=db_pal)+
  scale_x_continuous(limits = c(0,100))+
  scale_y_continuous(limits =c(0,5000))+
  xlab("Age")+ylab("No. of Signature1 + Signature5")+
  theme_bw()

####Age correlation with correction
if(F){
  m_meta_dt <- meta_dt %>% mutate(Sig1_ct = n_snv * Sig1_prop,
                                  Sig5_ct = n_snv * Sig5_prop,
                                  Sig7_ct = n_snv * Sig7_prop,
                                  corSig5_ct = Sig5_ct - 0.29*(Sig7_ct))
  ggplot(subset(m_meta_dt, Cell_type == 'skin_fb'), aes(x=age, y=Sig1_ct+corSig5_ct))+
    geom_point(aes(color=deadbody), size=3, alpha=0.5)+
    geom_smooth(method="lm")+
    geom_boxplot(aes(fill=deadbody),alpha=0.5)+
    scale_color_manual(values=db_pal)+
    scale_fill_manual(values=db_pal)+
    scale_x_continuous(limits = c(0,100))+
    scale_y_continuous(limits = c(0,5000))+
    xlab("Age")+ylab("No. of Signature1 + corrected Signature5")+
    theme_bw()
  summary(lm(m_meta_dt$Sig1_ct+m_meta_dt$Sig5_ct ~ m_meta_dt$age))
  summary(lm(m_meta_dt$Sig1_ct+m_meta_dt$corSig5_ct ~ m_meta_dt$age))
}


#plot cov and peakVAF
library(tidyverse)
library(ggrepel)
library(ggsci)
colnames(meta_dt)
m_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Lineage_tracing_summary_200331.txt')
ggplot(subset(meta_dt,is.na(peakVAF)==F), aes(x=seqz_cov, y=peakVAF))+
  geom_point(aes(color=current_final_for_lineage), alpha=0.5, size=3)+
  geom_text_repel(data=subset(meta_dt, is.na(peakVAF)==F & (peakVAF >= 45 & peakVAF < 55 & current_final_for_lineage == 'N')), aes(label = Memo))+
  scale_color_npg()+xlab("mean coverage") + ylab("peak VAF")+
  theme_bw(base_size = 18)


#plot cov and n_subs
library(ggrepel)
ggplot(subset(meta_dt,is.na(n_subs)==F), aes(x=seqz_cov, y=n_subs))+
  geom_point(aes(color=deadbody), alpha=0.5, size=3)+
  coord_cartesian(ylim=c(0,3000))+
  geom_text_repel(data=subset(meta_dt, is.na(n_subs)==F & n_subs <1100), aes(label = sample_id))

#plot sampleN per DB
library('scale')
db_pal <- pal_npg('nrc')(7)
names(db_pal) <- unique(meta_dt$deadbody)
f_meta_dt <- meta_dt %>% group_by(deadbody) %>% count()
ggplot(f_meta_dt, aes(x=reorder(deadbody,-n), y=n))+
  geom_bar(aes(fill = deadbody), stat="identity")+
  scale_fill_manual(values=db_pal)+xlab('')+ylab("No. of clonal samples")+
  theme_bw(base_size=18)+theme(legend.position = 'none')

#plot n_subs per DB
f_meta_dt <- meta_dt %>% filter(current_final_for_lineage =='Y')
f_meta_dt$deadbody <- factor(f_meta_dt$deadbody, levels = c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))
ggplot(f_meta_dt, aes(x=deadbody, y=n_subs))+
  geom_boxplot(aes(fill=deadbody))+
  geom_point()+
  scale_fill_manual(values = db_pal)+xlab('')+ylab('No. of somatic substitutions')+
  theme_bw(base_size=18)

#plot n_subs filled signature per sample
if(F){
  library(ggsci)
  library(scales)
  library(cowplot)
  sig_pal <- pal_npg('nrc')(3)
  names(sig_pal) <- c('Signature1', 'Signature5','Signature7')
  f_meta_dt <- meta_dt %>% filter(current_final_for_lineage =='Y')
  p <- list()
  n=0
  for (db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
    n <- n+1
    mf_meta_dt <- f_meta_dt %>% filter(deadbody == db) %>% 
      select(sample_id, n_snv, Sig1_ct,Sig5_ct,Sig7_ct) %>% 
      rename(Signature1 = Sig1_ct, Signature5 = Sig5_ct, Signature7 = Sig7_ct) %>% 
      gather(-sample_id, -n_snv, key=signature, value=count)
    mf_meta_dt$signature <- factor(mf_meta_dt$signature, levels = c('Signature1','Signature5','Signature7'))
    p[[n]] <-  ggplot(mf_meta_dt, aes(x=reorder(sample_id,-n_snv), y=count))+
      geom_bar(aes(fill = signature), stat="identity")+
      scale_fill_manual(values=sig_pal)+
      ylab("No. of substitutions")+xlab(db)+
      theme_bw(base_size=12)+theme(axis.text.x=element_text(size=8,angle=45, vjust=1, hjust=1), legend.position = 'none')
  }
  plot_grid(plotlist=p, ncol=3, align="hv")
  plot.new()
  plot(x=1)
  legend("topleft", legend=c("Sig1", "Sig5", "Sig7"), fill=sig_pal,
         bty="n", border="white") 
}


#simulation
n=10000
large_peak=32
A <- rpois(n*(large_peak*0.02), large_peak)
B <- rpois(n*(1-large_peak*0.02), 50-large_peak)
dt <- tibble(value = c(A,B))
ggplot(dt, aes(value))+
  geom_density()+
  scale_x_continuous(limits=c(0,100), breaks=seq(0,100,25))


#merge_and_sort_early_mutations_for_bait
library(tidyverse)
plng_tbl <- tribble(
  ~deadbody, ~plng_path,
  "DB2", '/home/users/sypark/00_Project/06_LineageTracing/db2/07_substitution/perLineage/', 
  "DB3", '/home/users/sypark/00_Project/06_LineageTracing/db3/04-3_substitution/re_processing_with_new_batch/perLineage/',
  "DB5",'/home/users/sypark/00_Project/06_LineageTracing/db5/06_substitution/perLineage/',
  "DB6",'/home/users/sypark/00_Project/06_LineageTracing/db6/07_substitution/filter3/perLineage/',
  "DB8",'/home/users/sypark/00_Project/06_LineageTracing/db8/05_substitution/perLineage/',
  "DB9",'/home/users/sypark/00_Project/06_LineageTracing/db9/05_substitution/perLineage/',
  "DB10",'/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution/perLineage/'
)

for (this_db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
  print(this_db)
  file_list <- list.files(plng_tbl$plng_path[plng_tbl$deadbody==this_db], pattern = '\\.subs.txt$', full.names = T)
  early_ids <- lng_dt %>% filter(deadbody == this_db & end_mtn < 300) %>% .$lineage_id
  res_tbl <- tibble()
  for (i in 1:length(file_list)){
    l_id <- gsub('.subs.txt','',unlist(strsplit(file_list[i],'//'))[2])
    dt <- read_tsv(file_list[i], col_types = cols(`#CHROM` = 'c', ALT = 'c', REF='c',`1000g2015aug_all` = 'c', `1000g2015aug_eas` = 'c'))
    if (l_id %in% early_ids){
      print(l_id)
      dt$lineage_id <- l_id
      dt$deadbody <- this_db
      res_tbl <- bind_rows(res_tbl, dt)
    }
  }
  res_tbl$`#CHROM` <- factor(res_tbl$`#CHROM`, levels = c(seq(1,22,1), 'X','Y'))
  res_tbl <- res_tbl %>% arrange(`#CHROM`,POS)
  dim(res_tbl)
  res_tbl %>% write_tsv(paste0(plng_tbl$plng_path[plng_tbl$deadbody==this_db],this_db,'_early_subs.txt'))
}      

# merge and edit bait design table
library(tidyverse)
db8_tbl <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db8/05_substitution/perLineage/DB8_early_subs_GsnpAdd.txt', col_types = cols(`#CHROM`='c'))
db9_tbl <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db9/05_substitution/perLineage/DB9_early_subs_GsnpAdd.txt', col_types = cols(`#CHROM`='c'))
db10_tbl <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution/perLineage/DB10_early_subs_GsnpAdd.txt', col_types = cols(`#CHROM`='c'))

merged_tbl <- bind_rows(db8_tbl,db9_tbl) %>% bind_rows(db10_tbl)
merged_tbl %>% mutate(distance = abs(POS - GermHeteroSNP)) %>%
  mutate(target_start_pos = ifelse(is.na(GermHeteroSNP)== T, POS-250,
                                   ifelse(distance < 400, pmin(POS,GermHeteroSNP) - (500 - distance)%/%2, POS-250))) %>%
  mutate(target_end_pos = ifelse(is.na(GermHeteroSNP)== T, POS+250,
                                 ifelse(distance < 400, pmax(POS,GermHeteroSNP) + (500 - distance)%/%2, POS+250))) %>%
  mutate(target_size= target_end_pos - target_start_pos +1) %>%
  mutate(target_id = paste0(deadbody, '_', `#CHROM`,':',POS,'_',REF,'_',ALT)) %>%
  select(target_id, target_chrom = `#CHROM`,target_start_pos, target_end_pos, target_size, POS,REF,ALT,lineage_id,deadbody,GermHeteroSNP) %>%
  write_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/bait_design/DB8_DB9_DB10_target_region.txt')





#Draw CIRCOS
library(tidyverse)
context_converter = c('C>A','C>A','C>G','C>G','C>T','C>T','T>A','T>A','T>C','T>C','T>G','T>G')
names(context_converter) = c('CA','GT','CG','GC','CT','GA','TA','AT','TC','AG','TG','AC')

SNV_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db6/07_substitution/filter3/perSample/6_LFAFb-4_001H12.snv.readc', col_types = cols(`#CHROM`='c'))
SNV_dt <- SNV_dt %>% separate(`self_ref_readN;var_readN;vaf%`, c( 'ref_readN','var_readN','vafpct'), sep=';') %>%
  mutate(chr = paste0('chr', `#CHROM`), start = POS, end = POS+1, vaf=as.numeric(vafpct)*0.01, type = context_converter[paste0(REF,ALT)]) %>%
  select(chr, start, end, vaf, type)
SNV_dt

load('/home/users/sypark/00_Project/06_LineageTracing/db6/09_sequenza/6_LFAFb-4_001H12.XX.100_20/6_LFAFb-4_001H12.XX.100_20_sequenza_extract.RData')
seqz_dt <- `6_LFAFb-4_001H12.XX.100_20_sequenza_extract`
rm(`6_LFAFb-4_001H12.XX.100_20_sequenza_extract`)

load('/home/users/sypark/00_Project/06_LineageTracing/db6/09_sequenza/input/6_LFAFb-4_001H12.XX.100_20.v2/6_LFAFb-4_001H12.XX.100_20.v2_sequenza_extract.RData')
seqz_dt <- `6_LFAFb-4_001H12.XX.100_20.v2_sequenza_extract`
rm(`6_LFAFb-4_001H12.XX.100_20.v2_sequenza_extract`)

CNV_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db6/09_sequenza/6_LFAFb-4_001H12.XX.100_20/6_LFAFb-4_001H12.XX.100_20_segments.txt')
CNV_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db6/09_sequenza/input/6_LFAFb-4_001H12.XX.100_20.v2/6_LFAFb-4_001H12.XX.100_20.v2_segments.txt', col_types = cols(chromosome = "c"))
CNV_dt$chr <- paste0("chr",CNV_dt$chromosome)
CNV_dt$start <- CNV_dt$start.pos
CNV_dt$end <- CNV_dt$end.pos
CNV_dt <- select(CNV_dt,chr,start,end,CNt,A,B)

SV_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db6/10_Delly/02_merged/6_LFAFb-4_001H12.delly.vcf.somatic.sort.pon.pre_fi.annotated.DB6_merged.SVvaf.fi.simple_edit.igv_check')
SV_dt <- SV_dt %>% mutate(chr1 = paste0('chr',`#CHR1`), chr2 = paste0('chr',CHR2), pos1 = POS1, pos2=POS2, ori=terinfo) %>% select(chr1, pos1, chr2, pos2, ori)
SV_dt$ori <- c("3to5"="3'-5'",
               "3to3"="3'-3'",
               "5to5"="5'-5'",
               "5to3"="5'-3'")[SV_dt$ori]

source('/home/users/sypark/01_Python_files/circos_R/CIRCOS_R_SYPark.R')
circos_VCS(SNV_dt, CNV_dt, SV_dt, seqz_dt, 1, 2, 'female')

#Blood VAF screening for removal of fasle positive
library(ggrepel)
path_dt <- tribble(~deadbody, ~si_path,
                   'DB2', '/home/users/sypark/00_Project/06_LineageTracing/db2/10_point_mutations/DB2_merged_subs_indels.txt.edit',
                   'DB3','/home/users/sypark/00_Project/06_LineageTracing/db3/14_point_mutations/DB3_merged_subs_indels.txt.edit',
                   'DB5', '/home/users/sypark/00_Project/06_LineageTracing/db5/11_point_mutations/DB5_merged_subs_indels.txt.edit',
                   'DB6','/home/users/sypark/00_Project/06_LineageTracing/db6/12_point_mutations/DB6_merged_subs_indels.txt.edit',
                   'DB8','/home/users/sypark/00_Project/06_LineageTracing/db8/09_point_mutations/DB8_merged_subs_indels.txt.edit',
                   'DB9','/home/users/sypark/00_Project/06_LineageTracing/db8/09_point_mutations/DB8_merged_subs_indels.txt.edit',
                   'DB10','/home/users/sypark/00_Project/06_LineageTracing/db10/09_point_mutations/DB10_merged_subs_indels.txt.edit')

i=6
dt <- read_tsv(as.character(path_dt[i,2]), col_types = cols(`#CHROM`='c'))
dt <- dt %>% separate(`blood_dp;blood_var;blood_vafpct`, c('blood_dp','blood_var','blood_vafpct'), sep = ';', convert = T) %>% 
  mutate(n_samples = apply(dt['sample_ids'], 1, function(x) length(unlist(strsplit(x, ';'))))) %>%
  mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep=':'))
ggplot(dt, aes(x=n_samples, y=blood_vafpct))+
  geom_jitter(alpha=0.5, width=0.2, height=0)+
  geom_text_repel(data = subset(dt, n_samples < 5 & blood_vafpct >30), aes(label = var_id))+
  theme_bw()+
  ggtitle(as.character(path_dt[i,1]))
dt %>% filter(grepl('6:10958045', var_id))


