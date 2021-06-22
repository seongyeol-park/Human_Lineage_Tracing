read_tsv('abcd')

library(tidyverse)
library(ggsci)
library(ggrepel)
library(ComplexHeatmap)

#load dataset
setwd("~/00_Project/06_LineageTracing/db3/04-3_substitution/previous_processing/splitted_incXY/03_perSample/mutalisk_output")
sig_dt <- read_tsv('DB3.perSample.SigsResult')

library(ggtree)
tree <- read.tree('/home/users/sypark/00_Project/06_LineageTracing/db3/04-3_substitution/re_processing_with_new_batch/DB3_Lineage_count_table.txt.nwk')
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Lineage_tracing_summary_200209.txt') %>% filter(deadbody == 'DB3' & current_final_for_lineage == 'Y')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_200209.txt') %>% filter(deadbody == 'DB3')


# DB3 SV count
setwd("~/00_Project/06_LineageTracing/db3/10_Delly/02_somatic_call/02_merged_call")
sv_dt <- read_tsv('DB3_36s.SVcount.edit.txt', comment = '##')
sv_dt <- sv_dt %>% rename(sample = `#sample`) %>% gather(-sample, -total_sv, key=sv_type, value=num)
ggplot(sv_dt, aes(x=reorder(sample, -num), y=num))+
  geom_bar(aes(fill = sv_type), stat = "identity")+
  scale_fill_jco()+
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(), axis.text.x = element_text(angle=90),
                     axis.title.x = element_blank())

# DB3 sub signature
sig_dt
m_sig_dt <- sig_dt %>% rename(sampleid = `#SampleID`) %>% mutate_at(vars(starts_with('Sig')), funs(.*Num_Mts))
gm_sig_dt <- m_sig_dt %>% gather(-sampleid, -Num_Mts, -Num_Sigs, -GOF, -CosSim, key=signature, value=num)
ggplot(gm_sig_dt, aes(x=reorder(sampleid, -Num_Mts), y=num))+
  geom_bar(aes(fill=signature), stat = "identity")+
  scale_fill_jco()+
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=15),
                     axis.title.x = element_blank(), axis.title = element_text(size=18))+
  ylab("No. of substitutions")

m_sig_dt$Sig5 %>% summary()
ggplot(m_sig_dt, aes(x=Sig7, y=Sig5))+
  geom_point(size=3, alpha=0.7)+
  geom_text_repel(data=subset(m_sig_dt, Sig7 > 1000), aes(label = sampleid))+
  scale_y_continuous(breaks=seq(0,10000,1000), limits = c(0,10000))+
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(), axis.text.x = element_text(angle=90), axis.text = element_text(size=15),
                     axis.title = element_text(size=18))+
  xlab("No. of signature7")+ylab("No. of signature5")
m_sig_dt %>% filter(Sig7 > 3000)


# DB3 subs heatmap reprocessed

#vaf_dt <- read_tsv('DB3_45s_merged.txt.sampleN.sub_fi.46sCall.rmFP')
setwd("~/00_Project/06_LineageTracing/db3/04-3_substitution/re_processing_with_new_batch")
vaf_dt <- read_tsv('DB3_45s_merged.txt.sampleN.fi.rmFP.Li_fi.sub_fi2.46sCall')
sample_list <- meta_dt %>% filter(deadbody == 'DB3') %>% pull(sample_id) %>% paste0(.,'_vafpct')
colnames(vaf_dt)
s_vaf_dt <- vaf_dt %>% rename(CHR = `#CHROM`) %>% select(CHR, POS, REF, ALT, sample_list, '3_blood_4s_merged_vafpct')
colnames(s_vaf_dt)
s_vaf_mx <- s_vaf_dt %>% mutate(snv_id=paste(CHR,POS,REF,ALT, sep='_')) %>% select(-CHR, -POS, -REF, -ALT, -`3_blood_4s_merged_vafpct`) %>% as.data.frame() %>% column_to_rownames('snv_id') %>% as.matrix()
library(ComplexHeatmap)
rt_annot = rowAnnotation(blood = s_vaf_dt$`3_blood_4s_merged_vafpct`)
pdf('candidate_heatmap.pdf', height=10, width=7)
pdf('error_detection_heatmap.pdf', height=10, width=7)
Heatmap(s_vaf_mx, col = colorRamp2(c(0,50,100), c("blue","yellow","red")),
        right_annotation = rt_annot, row_names_gp=gpar(fontsize=5), column_names_gp=gpar(fontsize=7),
        cell_fun = function(j,i,x,y,width,height,fill){
          grid.rect(x=x, y=y, width=width, height=height,gp = gpar(col="gray",fill=NA))
        })
dev.off()

# DB3 phylogenetic tree
library(RColorBrewer)
library(colorRamps)
library("scales")

m_meta_dt <- meta_dt %>% rename(taxa = lineage_id) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
m_meta_dt <- m_meta_dt %>% mutate(side_source2 = paste(Source_side, Source2, sep='_'))
m_meta_dt$side_source2[m_meta_dt$side_source2 == 'NA_pubis'] <- 'pubis'
unique(m_meta_dt$side_source2)


m_lng_dt <- lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
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
ggtree(tree) %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id), size=3, hjust=0) + coord_cartesian(xlim = c(0, 10000))
ggtree(tree) %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id), size=3, hjust=0) + coord_cartesian(xlim = c(0, 65000))
ggtree(tree) %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id), hjust=0) + coord_cartesian(xlim = c(0,100))+
  scale_x_continuous(breaks=seq(0,30,5))+
  theme(axis.text = element_text(size=20))


m_meta_dt <- meta_dt %>% filter(deadbody == 'DB3') %>% rename(taxa = lineage_id) %>% select(taxa, sample_id, Source, Cell_type, Source_side, Source_class,Source_1)
m_meta_dt <- m_meta_dt %>% mutate(source_side_name = paste(Source_side, Source_class, sep='_'))
m_meta_dt$source_side_name[m_meta_dt$source_side_name=='rt_HN'] <- 'HN'
m_meta_dt$source_side_name[m_meta_dt$source_side_name=='lt_HN'] <- 'HN'
m_meta_dt$source_side_name[m_meta_dt$source_side_name=='NA_HN'] <- 'HN'
m_meta_dt$source_side_name[m_meta_dt$source_side_name=='rt_trunk'] <- 'trunk'
m_meta_dt$source_side_name[m_meta_dt$source_side_name=='lt_trunk'] <- 'trunk'
unique(m_meta_dt$source_side_name)
m_lng_dt <- lng_dt %>% filter(deadbody == 'DB3') %>% rename(taxa = lineage_id) %>% select(-deadbody)
m_meta_dt <- left_join(m_lng_dt, m_meta_dt)



p %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id)) + scale_x_continuous(limits = c(0, 70000))
p %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id), hjust=0.3) + coord_cartesian(xlim = c(0,5000)) 
p %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id)) + coord_cartesian(xlim = c(0,100))+ geom_label(aes(x=node, label=sample_id), size=3, hjust=-1)
p %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id)) + coord_cartesian(xlim = c(0,30))+
  theme(axis.text = element_text(size=20))+
  scale_x_continuous(breaks = seq(0,30,2))
p %<+% m_meta_dt + geom_tree(aes(color = source_side_name),alpha=0.7, size=3,show.legend = TRUE)+theme_tree2() +
  geom_tiplab(aes(label = sample_id)) + 
  coord_cartesian(xlim = c(0,30))+ 
  scale_x_continuous(breaks = seq(0,50,2))+
  scale_color_manual(values = source_pal)+
  theme(axis.text = element_text(size=15))

# Tree with CNV presence
m_meta_dt <- m_meta_dt %>% mutate(chrYloss = ifelse(grepl("chrY_loss", CNV), 'Y','N')) 
m_meta_dt <- m_meta_dt %>% mutate(CNVpresence = ifelse(CNV != "none" & is.na(CNV)== F, 'Y',NA))
m_meta_dt$chrYloss[m_meta_dt$chrYloss=='N'] <- NA


p %<+% m_meta_dt + geom_tree(aes(color = CNVpresence),alpha=0.7, size=3,show.legend = TRUE)+theme_tree2() +
  geom_tiplab(aes(label = sample_id)) + 
  coord_cartesian(xlim = c(0,50))+ 
  scale_x_continuous(breaks = seq(0,50,5))+
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

#Trees with CCTTcount
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

# VAF plot
library(tidyverse)
library(cowplot)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/db3/04-3_substitution/re_processing_with_new_batch/perSample/',pattern = 'readc$', full.names = T)
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
res_tbl %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/db3/04-3_substitution/re_processing_with_new_batch/perSample/peakVAF.txt')
length(p)
png('/home/users/sypark/00_Project/06_LineageTracing/db3/04-3_substitution/re_processing_with_new_batch/perSample/peakVAF_plot4.png', width=14, height=7, res=150, units="in")
plot_grid(plotlist=p[37:45], ncol=4, align="hv")
dev.off()


#Extract coverage from sequenza
file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/db3/12_Sequenza/output/', pattern = 'extract.RData$', recursive = T, full.names = T)
dt <- tibble(sample_id = str_split_fixed(str_split_fixed(str_split_fixed(file_list,pattern ='//', n=2)[,2], pattern = '/', n=2)[,1],'\\.',2)[,1] , avg_t_dp = unlist(lapply(file_list, function(x) get(load(x))$avg.depth.tumor)), avg_n_dp = unlist(lapply(file_list, function(x) get(load(x))$avg.depth.normal)))
mean(dt$avg_n_dp)
dt %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/db3/12_Sequenza/output/cov_from_seqz.txt')

#RCS grouping
RCS11 <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db3/04-3_substitution/re_processing_with_new_batch/perSample/3_RCS_002A2.edited.vcf', col_types = cols(`#CHROM`= 'c'), comment='##')
RCS12 <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db3/04-3_substitution/re_processing_with_new_batch/perSample/3_RCS_002H5.edited.vcf', col_types = cols(`#CHROM`= 'c'), comment='##')
RCS2 <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db3/04-3_substitution/re_processing_with_new_batch/perSample/3_RCS_002H3.edited.vcf', col_types = cols(`#CHROM`= 'c'), comment='##')
RCS11 <- RCS11 %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep=':')) %>% select(var_id, `3_RCS_002A2`) %>% separate(`3_RCS_002A2`, c('ref11','var11'), convert = T) %>% mutate(dp11 = ref11+var11) %>% mutate(vaf11 = var11/dp11)
RCS12 <- RCS12 %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep=':')) %>% select(var_id, `3_RCS_002H5`) %>% separate(`3_RCS_002H5`, c('ref12','var12'), convert = T) %>% mutate(dp12 = ref12+var12) %>% mutate(vaf12 = var12/dp12)
RCS2 <- RCS2 %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep=':')) %>% select(var_id, `3_RCS_002H3`) %>% separate(`3_RCS_002H3`, c('ref2','var2'), convert = T) %>% mutate(dp2 = ref2+var2) %>% mutate(vaf2 = var2/dp2)
merge_dt <- full_join(full_join(RCS11, RCS12), RCS2)
ggplot(merge_dt, aes(x=vaf12, y=vaf11))+
  geom_point()+
  geom_density2d()+
  coord_cartesian(xlim=c(0,1), ylim=c(0,1))+
  xlab("3_RCS_002H5_VAF")+ylab("3_RCS_002A2_VAF")
ggplot(subset(merge_dt, is.na(vaf11)==T & is.na(vaf2) ==T))+
  geom_histogram(aes(vaf12), binwidth=0.02)+
  coord_cartesian(xlim=c(0,1))+
  ggtitle("Private subs of 3_RCS_002H5")
ggplot(subset(merge_dt, is.na(vaf12)==T & is.na(vaf2) ==T))+
  geom_histogram(aes(vaf11), binwidth=0.02)+
  coord_cartesian(xlim = c(0,1))+
  ggtitle("Private subs of 3_RCS_002A2")
ggplot(subset(merge_dt, is.na(vaf11)==T & is.na(vaf12) ==T))+
  geom_histogram(aes(vaf2), binwidth=0.02)+
  coord_cartesian(xlim = c(0,1))+
  ggtitle("Private subs of 3_RCS_002H3")


