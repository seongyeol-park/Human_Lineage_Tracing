read_tsv('abcd')


#load library
library(tidyverse)
library(ggsci)
theme_set(theme_bw(base_size = 18)) 
theme_set(theme_gray(base_size = 18)) 

#load dataset
library(ggtree)
tree <- read.tree('/home/users/sypark/00_Project/06_LineageTracing/db2/07_substitution/DB2_Lineage_count_table.txt.nwk')
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Lineage_tracing_summary_200131.txt') %>% filter(deadbody == 'DB2' & current_final_for_lineage == 'Y')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_200131.txt') %>% 


#HM for phylogenetic tree
library(ComplexHeatmap)
library(circlize)
dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db2/07_substitution/DB2_30s_merged.txt.sampleN.sam5.txt.31sCall')
colnames(dt)
dt <- dt %>% mutate(snv_id = paste0(`#CHROM`,':',POS,REF,'>',ALT))
dt <- dt %>% mutate(Nr1=as.integer(str_split_fixed(`dp;Nr1;Nr2;Nr3;Nr3ids`, ';',5)[,2]), Nr3=as.integer(str_split_fixed(`dp;Nr1;Nr2;Nr3;Nr3ids`, ';',5)[,4]))
#dt$total_ids <- paste(str_split_fixed(dt$`dp;Nr1;Nr2;Nr3;Nr3ids`, ';',5)[,5], dt$LOH_ids, sep=';')
#dt$Nr3_loh_n <- lengths(lapply(strsplit(dt$total_ids, ';'), function(x) unique(x[x != 'NA'])))
mVAF <- dt %>% select(ends_with('vafpct'),-`2_blood_3s_merged_vafpct`) %>% apply(., 1, function(x) mean(x[x>0 & is.na(x)==F]))
dt$meanVAF <- mVAF
f_dt <- dt %>% filter(Nr1-Nr3 < 5 & Nr3 - casNo < 5 & Nr1 < 28 & (Nr3 < 5 | meanVAF >=15) & Nr3 > 1)
f_dt <- dt %>% filter(Nr1-Nr3 < 5 & Nr3 - casNo < 5 &  Nr3 > 1)
#f_dt <- f_dt %>% filter(`2_RT3Fb3_001G8_vafpct`==0)
hm_mx <- f_dt %>% select(snv_id, ends_with('vafpct'), -`2_blood_3s_merged_vafpct`, -`2_CT_DP-6_002C3_vafpct`) %>% as.data.frame() %>% column_to_rownames('snv_id') %>% as.matrix()
colnames(hm_mx) <- gsub('_vafpct','',colnames(hm_mx))



dim(hm_mx)
rt_annot = rowAnnotation(blood = f_dt$`2_blood_3s_merged_vafpct`, mVAF=f_dt$meanVAF)
pdf('/home/users/sypark/00_Project/06_LineageTracing/db2/07_substitution/DB2_HM_for_tree_200106.pdf',height=12, width=10)
Heatmap(hm_mx, col = colorRamp2(c(0,50,100), c("blue","yellow","red")),
        right_annotation = rt_annot, row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=12)
        ,cell_fun = function(j,i,x,y,width,height,fill){
          grid.rect(x=x, y=y, width=width, height=height,gp = gpar(col="gray",fill=NA))
        }
)
dev.off()


#check filtered HM
library(ComplexHeatmap)
library(circlize)
dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db2/07_substitution/DB2_30s_merged.txt.sampleN.loh_id.g_fi.fi.sam5.txt.31sCall')
dim(dt)
colnames(dt)
dt <- dt %>% mutate(snv_id = paste0(`#CHROM`,':',POS,REF,'>',ALT))
dt <- dt %>% mutate(Nr1=as.integer(str_split_fixed(`dp;Nr1;Nr2;Nr3;Nr3ids`, ';',5)[,2]), Nr3=as.integer(str_split_fixed(`dp;Nr1;Nr2;Nr3;Nr3ids`, ';',5)[,4]))
#dt$total_ids <- paste(str_split_fixed(dt$`dp;Nr1;Nr2;Nr3;Nr3ids`, ';',5)[,5], dt$LOH_ids, sep=';')
#dt$Nr3_loh_n <- lengths(lapply(strsplit(dt$total_ids, ';'), function(x) unique(x[x != 'NA'])))
mVAF <- dt %>% select(ends_with('vafpct'),-`2_blood_3s_merged_vafpct`,-`2_CT_DP-6_002C3_vafpct`) %>% apply(., 1, function(x) mean(x[x>0 & is.na(x)==F]))
dt$meanVAF <- mVAF
f_dt <- dt
hm_mx <- f_dt %>% select(snv_id, ends_with('vafpct'), -`2_blood_3s_merged_vafpct`, -`2_CT_DP-6_002C3_vafpct`) %>% as.data.frame() %>% column_to_rownames('snv_id') %>% as.matrix()
colnames(hm_mx) <- gsub('_vafpct','',colnames(hm_mx))
dim(hm_mx)
rt_annot = rowAnnotation(blood = f_dt$`2_blood_3s_merged_vafpct`, mVAF=f_dt$meanVAF)
pdf('/home/users/sypark/00_Project/06_LineageTracing/db2/07_substitution/DB2_filtered_HM_200106.pdf',height=12, width=10)
Heatmap(hm_mx, col = colorRamp2(c(0,50,100), c("blue","yellow","red")),
        right_annotation = rt_annot, row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=12)
        ,cell_fun = function(j,i,x,y,width,height,fill){
          grid.rect(x=x, y=y, width=width, height=height,gp = gpar(col="gray",fill=NA))
        }
)
dev.off()


#Extract coverage from sequenza
file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/db2/05_sequenza/', pattern = 'extract.RData$', recursive = T, full.names = T)
dt <- tibble(sample_id = str_split_fixed(str_split_fixed(str_split_fixed(file_list,pattern ='//', n=2)[,2], pattern = '/', n=2)[,1],'\\.',2)[,1] , avg_t_dp = unlist(lapply(file_list, function(x) get(load(x))$avg.depth.tumor)), avg_n_dp = unlist(lapply(file_list, function(x) get(load(x))$avg.depth.normal)))
mean(dt$avg_n_dp)
dt %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/db2/05_sequenza/cov_from_seqz.txt')


###### phylogenetic tree
library(RColorBrewer)
library(colorRamps)
library("scales")

m_meta_dt <- meta_dt %>% rename(taxa = lineage_id) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
m_meta_dt <- m_meta_dt %>% mutate(side_source2 = paste(Source_side, Source2, sep='_'))
m_meta_dt$side_source2[m_meta_dt$side_source2 == 'NA_occiput'] <- 'occiput'
unique(m_meta_dt$side_source2)


m_lng_dt <- lng_dt %>% filter(deadbody == 'DB2') %>% rename(taxa = lineage_id) %>% select(-deadbody)
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

ggtree(tree) %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id), size=3) + scale_x_continuous(limits = c(0, 5000))
ggtree(tree) %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id), hjust=0) + coord_cartesian(xlim = c(0,30))+
  scale_x_continuous(breaks=seq(0,30,5))+
  theme(axis.text = element_text(size=20))


#Trees with anatomical information
#pdf('/home/users/sypark/00_Project/06_LineageTracing/db6/07_substitution/filter3/DB6_phylo_proximal_40_191204.pdf', height=30, width=20)
source_pal = colorRampPalette(brewer.pal(8,"Set1"))(4)
names(source_pal) = unique(m_meta_dt$side_source2)
show_col(source_pal)

p <- ggtree(tree) %<+% m_lng_dt
p + geom_tree(aes(color = side_source2),alpha=0.7, size=3,show.legend = TRUE)+theme_tree2() +
  geom_tiplab(aes(label = sample_id)) + 
  coord_cartesian(xlim = c(0,30))+ 
  scale_x_continuous(breaks = seq(0,50,2))+
  scale_color_manual(values = source_pal)+
  theme(axis.text = element_text(size=15))
p + geom_tree(aes(color = side_source2),alpha=0.7, size=3,show.legend = TRUE)+theme_tree2() +
  geom_tiplab(aes(label = sample_id)) + 
  scale_x_continuous(limits = c(0,5000),breaks = seq(0,5000,1000))+
  scale_color_manual(values = source_pal)+
  theme(axis.text = element_text(size=15))

#dev.off()
plot.new()
plot(x=1)
legend(1,1.4,legend=unique(m_meta_dt$side_source2),
       col=source_pal, lty=1, cex=2, lwd=4)

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


# Draw summary plot of CNV
library(cowplot)
cyto_dt <- read_tsv('/home/users/sypark/03_Tools/annovar/humandb/hg19_cytoBand.txt', col_names = c('chr','start','end','cyto_id','stain'))
acen_dt <- cyto_dt %>% filter(stain == 'acen') 
acen_dt2 <- acen_dt %>% select(start) %>% aggregate(list(acen_dt$chr), min) %>% as.tibble()
acen_dt2$end <- acen_dt %>% select(end) %>% aggregate(list(acen_dt$chr), max) %>% .$end
acen_dt2 <- acen_dt2 %>% mutate(mean = (start+end)/2, chrom = gsub('chr','',Group.1))

file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/db2/05_sequenza/clean_seg/', pattern = 'clean$', full.names = T)
for (fn in file_list){
  print(fn)
  dt <- read_tsv(fn)
  p <- list()
  n=0
  for (chrom in unique(dt$`#CHROM`)){
    n = n+1
    subdt <- dt %>% filter(`#CHROM` == chrom)
    p[[n]] <- ggplot(subdt)+
      geom_segment(aes(x=start_pos, xend=end_pos, y=majCN, yend=majCN), color = "red", size=5, alpha=0.5)+
      geom_segment(aes(x=start_pos, xend=end_pos, y=minCN, yend=minCN), color = "blue",size=5, alpha=0.5 )+
      geom_vline(xintercept = acen_dt2$mean[acen_dt2$chrom == chrom],color = "green", size=2, alpha=0.5)+
      scale_y_continuous(limits=c(0,4), breaks=seq(0,4,1))+
      xlab(chrom)+
      theme_bw() + theme(plot.margin=unit(c(0,0,0,0),"cm"), axis.title.y = element_blank())
  }
  png(paste0(fn, '.png'), width=12, height=2, units = "in", res=100)
  print(plot_grid(plotlist = p, nrow=1))
  dev.off()
}

# VAF plot
library(tidyverse)
library(cowplot)
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/db2/07_substitution/perSample/',pattern = 'readc$', full.names = T)
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
res_tbl %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/db2/07_substitution/perSample/peakVAF.txt')
length(p)
png('/home/users/sypark/00_Project/06_LineageTracing/db2/07_substitution/perSample/peakVAF_plot3.png', width=14, height=7, res=150, units="in")
plot_grid(plotlist=p[25:30], ncol=4, align="hv")
dev.off()


#signature n_subs plot
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Lineage_tracing_summary_200131.txt')
f_meta_dt <- meta_dt %>% filter(current_final_for_lineage =='Y' & deadbody == 'DB2')
f_meta_dt <- f_meta_dt %>% select(sample_id, n_subs, ends_with('_ct')) %>% gather(-sample_id, -n_subs, key='signature', value='count')

ggplot(f_meta_dt, aes(x=sample_id, y=count))+
  geom_bar(aes(fill = signature), stat="identity")+
  xlab('')+ylab('No. of substitutions')+
  theme_bw(base_size=18)+theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
