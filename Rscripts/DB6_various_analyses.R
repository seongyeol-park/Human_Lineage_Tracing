read_tsv('abcd')

#load library
library(tidyverse)
library(ggsci)
theme_set(theme_bw(base_size = 18)) 
theme_set(theme_gray(base_size = 18)) 

#load dataset
library(ggtree)
tree <- read.tree('/home/users/sypark/00_Project/06_LineageTracing/db6/07_substitution/filter3/DB6_Lineage_count_table.txt.nwk')
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Lineage_tracing_summary_200209.txt') %>% filter(deadbody == 'DB6' & current_final_for_lineage == 'Y')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_200209.txt') %>% filter(deadbody == 'DB6')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_200209m.txt') %>% filter(deadbody == 'DB6')



# color setting
library(ggsci)
library("scales")
my_pal = pal_aaas("default")(10)
my_pal = pal_jco("default")(10)
show_col(my_pal)



# DB6 subs count
m_sig_dt <- meta_dt %>% mutate_each(funs(.*n_subs), starts_with('Sig')) %>% select(sample_id, n_subs, CS_signature, starts_with('Sig'))
colnames(m_sig_dt) <- gsub('prop','ct',colnames(m_sig_dt))
gm_sig_dt <- m_sig_dt %>% gather(-sample_id, -n_subs, -CS_signature, key=signature, value=num)
ggplot(gm_sig_dt, aes(x=reorder(sample_id, -n_subs), y=num))+
  geom_bar(aes(fill=signature), stat = "identity")+
  scale_fill_jco()+
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(), axis.text.x = element_text(angle=90),
                     axis.title.x = element_blank())

m_meta_dt <- meta_dt %>% mutate(source_class2 = paste(Cell_type, Source_side, Source_direction, Source_1, sep=':'),
                                corSig5_ct = Sig5_ct - 0.26*Sig7_ct)
multiple_values <- m_meta_dt %>% select(source_class2) %>% table() %>% .[.>1] %>% names() 
ggplot(subset(m_meta_dt, source_class2 %in% multiple_values), aes(x=reorder(source_class2, -Sig7_ct), y=Sig7_ct))+
  geom_jitter(size=3, alpha=0.5, height=0, width=0.2)+
  theme(axis.text.x =element_text(angle=90))

ggplot(subset(m_meta_dt, source_class2 %in% multiple_values), aes(x=reorder(source_class2, -(Sig1_ct + corSig5_ct)), y=Sig1_ct+corSig5_ct))+
  geom_jitter(size=3, alpha=0.5, height=0, width=0.2)+
  theme(axis.text.x =element_text(angle=90))



#DB6 sig7 sig5 correlation
library(ggrepel)
ggplot(m_sig_dt, aes(x=Sig7_ct, y=Sig5_ct))+
  geom_point()+
  geom_segment(aes(x=0, xend=20000, y=2100, yend=7300), linetype="dashed")+
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line())+
  geom_text_repel(data=subset(m_sig_dt, Sig5_ct > 4700 & Sig7_ct < 7000), aes(label = sample_id))

ggplot(m_sig_dt, aes(x=Sig7_ct, y=Sig1_ct))+
  geom_point()+
  scale_x_continuous(breaks=seq(0,20000,1000))+
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line())
  geom_text_repel(data=subset(m_sig_dt, Sig7_ct >))
#y=0.26x+2100
# sig5 correction
m_sig_dt <- m_sig_dt %>% mutate(corSignature5 = Sig5_ct - (Sig7_ct*0.26))
gm_sig_dt <- m_sig_dt %>% gather(-c(sample_id, n_subs, CS_signature, Sig7_ct, Sig5_ct),key=signature, value=num)
ggplot(gm_sig_dt, aes(x=reorder(sample_id, -num), y=num))+
  geom_bar(aes(fill = signature), stat = "identity")+
  scale_fill_jco()+
  scale_y_continuous(limits = c(0,5000))+
  theme_bw()+theme(axis.text.x = element_text(angle=90),axis.title.x = element_blank())
outliers <- c('6_RMCS-34_001F2','6_MS-11_001B3','6_RLCS-21_001B12','6_LtHand-5_001H3','6_RLCS-23_001H9')

fgm_sig_dt <- gm_sig_dt %>% filter(!sample_id %in% outliers)
ggplot(fgm_sig_dt, aes(x=reorder(sample_id, -num), y=num))+
  geom_bar(aes(fill = signature), stat = "identity")+
  scale_fill_jco()+
  scale_y_continuous(limits = c(0,5000))+
  theme_bw()+theme(axis.text.x = element_text(angle=90),axis.title.x = element_blank())

mm_sig_dt <- left_join(m_sig_dt, meta_dt, by='sample_id')
mm_sig_dt
ggplot(mm_sig_dt, aes(x=Source_1, y=Sig1_ct + corSignature5))+
  geom_violin()+
  geom_jitter(aes(color = Source_direction),size=3, alpha=0.5, width=0.1, height=0)+
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(), axis.text.x= element_text(angle=90))
fmm_sig_dt <- mm_sig_dt %>% filter(Source_1 %in% c('arm','deltoid', 'elbow', 'hand','rectus_femoris','leg','foot')) %>% mutate(Source_2 = ifelse(Source_1 == 'arm', paste0(Source_direction, Source_1), Source_1))
fmm_sig_dt$Source_1 <- factor(fmm_sig_dt$Source_1, levels = c('deltoid','elbow','arm','hand', 'rectus_femoris','leg','foot'))
fmm_sig_dt <- fmm_sig_dt %>% filter(!sample_id %in% outliers)
ggplot(fmm_sig_dt, aes(x=Source_1, y=Sig1_ct + corSignature5))+
  geom_violin()+
  geom_jitter(aes(color = Source_1),size=3, alpha=0.5, width=0.1, height=0)+
  scale_color_jco()+
  scale_y_continuous(limits =c (0,4000))+
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(), 
                     axis.text.x= element_text(angle=90, size= 15), axis.text.y = element_text(size=15),
                     axis.title.x= element_blank(), axis.title.y = element_text(size=18))

fmm_sig_dt <- mm_sig_dt %>% filter(Source_1 %in% c('rectus_femoris','leg', 'foot')) 
fmm_sig_dt$Source_1 <- factor(fmm_sig_dt$Source_1, levels = c('rectus_femoris','leg', 'foot'))
ggplot(fmm_sig_dt, aes(x=Source_1, y=Signature1 + corSignature5))+
  geom_violin()+
  geom_jitter(aes(color = Source_1),size=3, alpha=0.5, width=0.1, height=0)+
  scale_color_jco()+
  scale_y_continuous(limits =c (0,4000))+
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(), axis.text.x= element_text(angle=90))


#DB6 phylogenetic tree
m_meta_dt <- meta_dt %>% rename(taxa = lineage_id) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
m_meta_dt <- m_meta_dt %>% mutate(side_source2 = paste(Source_side, Source2, sep='_'))
unique(m_meta_dt$side_source2)

m_lng_dt <- lng_dt %>% filter(deadbody == 'DB6') %>% rename(taxa = lineage_id) %>% select(-deadbody)
m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
#fill data of nodes if information of it's branches is unique
node_list <- m_lng_dt %>% filter(is.na(Cell_type)== T) %>% .$taxa
for (node_id in node_list){
  print(node_id)
  for(col_id in c('Source','Cell_type','Source_side','Source_direction','Source_class','Source2','side_source2')){
    print(col_id)
    ex_values <- m_lng_dt %>% filter(grepl(node_id, taxa)==T & is.na(get(col_id)) == F) %>% .[[col_id]] %>% unique()
    if (length(ex_values)==1){
      m_lng_dt[m_lng_dt$taxa == node_id,][[col_id]] <- ex_values
    }
  }
}
node_dt <- m_lng_dt %>% filter(n_samples >1)

#default trees
ggtree(tree) %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id), size=3) + coord_cartesian(xlim = c(0, 15000))
ggtree(tree) %<+% m_meta_dt + geom_tree()+theme_tree2() + coord_cartesian(xlim = c(0, 20000))
ggtree(tree) %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id), hjust=0) + coord_cartesian(xlim = c(0,30))+
  scale_x_continuous(breaks=seq(0,30,5))+
  theme(axis.text = element_text(size=20))





library(ggtree)
library(tidyverse)

m_meta_dt <- meta_dt %>% filter(deadbody == 'DB6') %>% rename(taxa = lineage_id) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
m_meta_dt <- m_meta_dt %>% mutate(source_side_name = paste(Source_side, Source_class, sep='_'))
m_meta_dt$source_side_name[m_meta_dt$source_side_name=='rt_HN'] <- 'HN'
m_meta_dt$source_side_name[m_meta_dt$source_side_name=='lt_HN'] <- 'HN'
m_meta_dt$source_side_name[m_meta_dt$source_side_name=='NA_HN'] <- 'HN'
m_meta_dt$source_side_name[m_meta_dt$source_side_name=='rt_trunk'] <- 'trunk'
m_meta_dt$source_side_name[m_meta_dt$source_side_name=='lt_trunk'] <- 'trunk'
unique(m_meta_dt$source_side_name)

m_lng_dt <- lng_dt %>% filter(deadbody == 'DB6') %>% rename(taxa = lineage_id) %>% select(-deadbody)
m_meta_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')

library(RColorBrewer)
library(colorRamps)
source_pal = colorRampPalette(brewer.pal(8,"Set1"))(8)
names(source_pal) = unique(m_meta_dt$source_side_name)
show_col(source_pal)
#source_pal2 = colorRampPalette(brewer.pal(8,"Set1"))(8)[c(1,1,3,1,5,2,2,2)]  # combine UE/HN and LE
#names(source_pal2) = unique(m_meta_dt$source_side_name)
#show_col(source_pal2)
#source_pal3 = colorRampPalette(brewer.pal(8,"Set1"))(8)[c(1,2,3,1,5,1,2,4)]  # combine UE/HN and LE
#names(source_pal3) = unique(m_meta_dt$source_side_name)
#show_col(source_pal3)

#default trees

p %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id), size=3) + scale_x_continuous(limits = c(0, 30000))
p %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id), hjust=0) + coord_cartesian(xlim = c(0,30))+
  scale_x_continuous(breaks=seq(0,30,2))+
  theme(axis.text = element_text(size=20))

#Trees with anatomical information
#pdf('/home/users/sypark/00_Project/06_LineageTracing/db6/07_substitution/filter3/DB6_phylo_proximal_40_191204.pdf', height=30, width=20)
p %<+% m_meta_dt + geom_tree(aes(color = source_side_name),alpha=0.7, size=3,show.legend = TRUE)+theme_tree2() +
  geom_tiplab(aes(label = sample_id)) + 
  coord_cartesian(xlim = c(0,30))+ 
  scale_x_continuous(breaks = seq(0,50,2))+
  scale_color_manual(values = source_pal)+
  theme(axis.text = element_text(size=15))
#dev.off()
plot.new()
plot(x=1)
legend(1,1.4,legend=unique(m_meta_dt$source_side_name),
       col=source_pal2, lty=1, cex=2, lwd=4)

#Draw relationship among anatomical sectors
link_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db6/07_substitution/filter3/linked_different_sources_191212.txt')
link_dt
source_list <- sort(unique(c(link_dt$mate1, link_dt$mate2)))
library(circlize)
circos.par("track.height" = 0.1)
circos.initialize(factors = source_list, xlim=c(0,1))
circos.track(factors =source_list, ylim=c(0,1),
             panel.fun = function(x, y) {
             })

for (source_name in source_list){
  circos.text(0.5, 0.5, source_name, sector.index = source_name, track.index = 1)
}
for (n in 1:nrow(link_dt)){
  circos.link(link_dt$mate1[n], sample(30:70, 1, replace = TRUE)/100 , link_dt$mate2[n], sample(30:70, 1, replace = TRUE)/100)
}

# Tree with CNV presence
m_meta_dt <- m_meta_dt %>% mutate(CNVpresence = ifelse(CNV != "none" & is.na(CNV)== F, 'Y',NA))
p %<+% m_meta_dt + geom_tree(aes(color = CNVpresence),alpha=0.6, size=3,show.legend = TRUE)+theme_tree2() +
  geom_tiplab(aes(label = sample_id)) + 
  coord_cartesian(xlim = c(0,30))+ 
  scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.text = element_text(size=15))

m_meta_dt <- m_meta_dt %>% mutate(Xloss = ifelse(CNV =='chrX_loss', 'Y',NA))
p %<+% m_meta_dt + geom_tree(aes(color = Xloss),alpha=0.6, size=3,show.legend = TRUE)+theme_tree2() +
  geom_tiplab(aes(label = sample_id)) + 
  coord_cartesian(xlim = c(0,30))+ 
  scale_x_continuous(breaks = seq(0,30,5))+
  theme(axis.text = element_text(size=15))

#Trees with blood VAF
p %<+% m_meta_dt + geom_tree(aes(color = mean_blood_VAFpct),alpha=0.6, size=3,show.legend = TRUE)+theme_tree2() +
  geom_tiplab(aes(label = sample_id)) + 
  coord_cartesian(xlim = c(0,30))+ 
  scale_x_continuous(breaks = seq(0,50,2))+
  scale_color_gradient(low ="gray",high = "red", breaks = c(0,30))+
  theme(axis.text = element_text(size=15))

colfunc <- colorRampPalette(c("red", "gray"))
legend_image <- as.raster(matrix(colfunc(20), ncol=1))
plot(c(0,2),c(0,30),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
text(x=0.5, y = seq(0,30,10), labels = paste0(seq(0,30,10),'%'), cex=2.5)
rasterImage(legend_image, 0, 0, 0.3,30)


#Trees with CCTT prop
m_meta_dt <- m_meta_dt %>% mutate(CCTTprop = CCTTcount*2/n_subs)
summary(m_meta_dt$CCTTprop)

p %<+% m_meta_dt + geom_tree(aes(color = CCTTprop),alpha=0.6, size=3,show.legend = TRUE)+theme_tree2() +
  geom_tiplab(aes(label = sample_id)) + 
  coord_cartesian(xlim = c(0,30))+ 
  scale_x_continuous(breaks = seq(0,30,2))+
  scale_color_gradient(low ="gray",high = "red", breaks = c(0,30))+
  theme(axis.text = element_text(size=15))

colfunc <- colorRampPalette(c("red", "gray"))
legend_image <- as.raster(matrix(colfunc(20), ncol=1))
plot(c(0,2),c(0,30),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
text(x=0.5, y = seq(0,30,10), labels = paste0(seq(0,30,10),'%'), cex=2.5)
rasterImage(legend_image, 0, 0, 0.3,30)


#chrX subs count and signature
m_meta_dt <- meta_dt %>% mutate(Xloss = ifelse(CNV =='chrX_loss', 'Y',NA)) %>%
  mutate_at(vars(starts_with('chrX_sig')), ~(.*chrX_n_subs) ) %>% rename_at(vars(starts_with('chrX_sig')), ~(paste0(.,'_ct')))
ggplot(m_meta_dt, aes(x=Sig7_ct, y=chrX_sig7_ct))+
  geom_point(aes(color=Xloss))+
  scale_y_continuous(limits=c(0,500))+scale_x_continuous(limits=c(0,10000))+
  geom_abline(slope=165/2500)+geom_abline(slope=82.5/2500, linetype="dashed")
ggplot(m_meta_dt, aes(x=Sig5_ct, y=chrX_sig5_ct))+
  geom_point(aes(color=Xloss))+
  coord_cartesian(xlim=c(0,8000),ylim=c(0,500))+
  geom_abline(slope=165/2500)+geom_abline(slope=82.5/2500, linetype="dashed")
ggplot(m_meta_dt, aes(x=n_subs, y=chrX_n_subs))+
  geom_point(aes(color=Xloss))+
  coord_cartesian(xlim=c(0,30000),ylim=c(0,2000))+
  geom_abline(slope=165/2500)+geom_abline(slope=82.5/2500, linetype="dashed")
ggplot(m_meta_dt, aes(x=n_subs, y=chrX_CtTG))+
  geom_point(aes(color = Xloss))+
  coord_cartesian(xlim=c(0,30000),ylim=c(0,100))+
  geom_abline(slope=30/10000)+geom_abline(slope=15/10000,linetype="dashed")
ggplot(m_meta_dt, aes(x=Sig7_ct, y=chrX_TCtTC))+
  geom_point(aes(color= Xloss))+
  geom_abline(slope=96/5000)+geom_abline(slope=48/5000, linetype="dashed")
ggplot(m_meta_dt, aes(x=Sig7_ct, y=chrX_TCtTC_intragen))+
  geom_point(aes(color= Xloss))+
  geom_abline(slope=28/5000)+geom_abline(slope=14/5000, linetype="dashed")
ggplot(m_meta_dt, aes(x=n_CCtTT, y=chrX_CCtTT))+
  geom_point(aes(color= Xloss))+
  geom_abline(slope=62/1000)+geom_abline(slope=31/1000, linetype="dashed")
ggplot(m_meta_dt, aes(x=n_subs, y=chrX_n_subs_intragen))+
  geom_point(aes(color= Xloss))+
  geom_abline(slope=200/10000)+geom_abline(slope=100/10000, linetype="dashed")


# VAF plot
library(tidyverse)
library(cowplot)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/db6/07_substitution/filter3/perSample/',pattern = 'readc$', full.names = T)
length(file_list)
res_tbl <- tibble(id=character(),peakVAF=numeric())
p<- list()
for (i in 1:length(file_list)){
  print(file_list[i])
  dt <- read_tsv(file_list[i])
  sample_id <- unlist(strsplit(unlist(strsplit(file_list[i],'//'))[2],'\\.'))[1]
  dt <- dt %>% separate(`self_ref_readN;var_readN;vaf%`, c('refn','varn','vaf'), sep=';', convert = T)
  max_i <- which.max(density(dt$vaf)$y)
  p[[i]] <- ggplot(dt)+
    geom_density(aes(vaf), fill="gray")+
    geom_vline(xintercept=density(dt$vaf)$x[max_i], color="red")+
    geom_text(aes(x=80, y=0.05), label = paste0('n = ', nrow(dt)))+
    geom_text(aes(x=80, y=0.045), label = paste0('peakVAF = ',round(density(dt$vaf)$x[max_i],2)))+
    scale_x_continuous(limits=c(0,100), breaks = seq(0,100,25))+
    ggtitle(sample_id)+
    theme_bw()
  res_tbl[i,] <- list(sample_id, round(density(dt$vaf)$x[max_i],2))
}
res_tbl
res_tbl %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/db6/07_substitution/filter3/perSample/peakVAF.txt')
length(p)
png('/home/users/sypark/00_Project/06_LineageTracing/db6/07_substitution/filter3/perSample/peakVAF_plot10.png', width=14, height=7, res=150, units="in")
plot_grid(plotlist=p[109:124], ncol=4, align="hv")
dev.off()


#VAF distribution of 6_MS-8_001G3
dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db6/07_substitution/filter3/perSample/6_MS-8_001G3.snv.readc', col_types = cols(`#CHROM`='c'))
dt <- dt %>% separate(`self_ref_readN;var_readN;vaf%`, c('refn','varn','vaf'), sep=';', convert = T)
dt %>% group_by(`#CHROM`) %>% summarise(chrom_mean_vaf = mean(vaf), n=n()) %>% View()
ggplot(dt, aes(x=`#CHROM`, y=vaf))+
  geom_violin()+
  geom_hline(yintercept=50, color="red")+
  scale_y_continuous(limits=c(0,100))


#Extract coverage from sequenza
file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/db6/09_sequenza/', pattern = 'extract.RData$', recursive = T, full.names = T)
dt <- tibble(sample_id = str_split_fixed(str_split_fixed(str_split_fixed(file_list,pattern ='//', n=2)[,2], pattern = '/', n=2)[,1],'\\.',2)[,1] , avg_t_dp = unlist(lapply(file_list, function(x) get(load(x))$avg.depth.tumor)), avg_n_dp = unlist(lapply(file_list, function(x) get(load(x))$avg.depth.normal)))
mean(dt$avg_n_dp)
dt %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/db6/09_sequenza/cov_from_seqz.txt')
