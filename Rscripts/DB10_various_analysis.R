read_tsv('abcd')

#load library
library(tidyverse)
library(ggsci)
theme_set(theme_bw(base_size = 18)) 
theme_set(theme_gray(base_size = 18)) 

#load dataset
library(ggtree)
tree <- read.tree('/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution/DB10_Lineage_count_table.txt.nwk')
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Lineage_tracing_summary_200209.txt') %>% filter(deadbody == 'DB10' & current_final_for_lineage == 'Y')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_200209.txt') %>% filter(deadbody == 'DB10')


#Draw heatmap for HM
dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution/DB10_50s_merged.txt.sampleN.loh_id.g_fi.fi.sam3.txt.51sCall')
colnames(dt)
dt <- dt %>% mutate(snv_id = paste0(`#CHROM`,':',POS,REF,'>',ALT))
dt <- dt %>% mutate(Nr1=as.integer(str_split_fixed(`dp;Nr1;Nr2;Nr3;Nr3ids`, ';',5)[,2]), Nr3=as.integer(str_split_fixed(`dp;Nr1;Nr2;Nr3;Nr3ids`, ';',5)[,4]))
mVAF <- dt %>% select(ends_with('vafpct'),-c('10_Blood_3s_merged_vafpct','10_ALL13-4_006H8_vafpct','10_ALL13-8_001G5_vafpct','10_ARL10-4_001D4-B6_vafpct','10_ARL7-2_005B2_vafpct',
                      '10_LA14-2_002D1-A5_vafpct','10_LAFb16-3_001E11_vafpct','10_PRL3-1_001B8-E10_vafpct','10_PRLM2-3_001D10_vafpct','10_PLL6-3_002A11_vafpct')) %>% apply(., 1, function(x) mean(x[x>0 & is.na(x)==F]))
dt$meanVAF <- mVAF
f_dt <- dt %>% filter(Nr1-Nr3 < 5 & Nr3 - casNo < 5 & Nr1 <48 & Nr3 > 3)
fs_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution/DB10_false_snv_list.txt', col_names = c('chr','pos','ref','alt'))
fs_dt <- fs_dt %>% mutate(snv_id = paste0(chr,':',pos,ref,'>',alt))
f_dt <- f_dt %>% filter(!snv_id %in% fs_dt$snv_id)
hm_mx <- f_dt %>% select(snv_id, ends_with('vafpct'), -c('10_Blood_3s_merged_vafpct','10_ALL13-4_006H8_vafpct','10_ALL13-8_001G5_vafpct','10_ARL10-4_001D4-B6_vafpct','10_ARL7-2_005B2_vafpct',
                                                         '10_LA14-2_002D1-A5_vafpct','10_LAFb16-3_001E11_vafpct','10_PRL3-1_001B8-E10_vafpct','10_PRLM2-3_001D10_vafpct','10_PLL6-3_002A11_vafpct')) %>% as.data.frame() %>% column_to_rownames('snv_id') %>% as.matrix()
colnames(hm_mx) <- gsub('_vafpct','',colnames(hm_mx))
library(ComplexHeatmap)
library(circlize)
dim(hm_mx)
rt_annot = rowAnnotation(blood = f_dt$`10_Blood_3s_merged_vafpct`, mVAF=f_dt$meanVAF)
pdf('/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution/DB10_HM_for_tree_200108.pdf',height=12, width=10)
Heatmap(hm_mx, col = colorRamp2(c(0,50,100), c("blue","yellow","red")),
        right_annotation = rt_annot, row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=12)
        ,cell_fun = function(j,i,x,y,width,height,fill){
          grid.rect(x=x, y=y, width=width, height=height,gp = gpar(col="gray",fill=NA))
        }
)
dev.off()


f_dt <- dt %>% filter(Nr1-Nr3 < 5 & Nr3 - casNo < 5 & Nr1 <48 & Nr3 == 2)
fs_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution/DB10_false_snv_list.txt', col_names = c('chr','pos','ref','alt'))
fs_dt <- fs_dt %>% mutate(snv_id = paste0(chr,':',pos,ref,'>',alt))
f_dt <- f_dt %>% filter(!snv_id %in% fs_dt$snv_id)
hm_mx <- f_dt %>% select(snv_id, ends_with('vafpct'), -c('10_Blood_3s_merged_vafpct','10_ALL13-4_006H8_vafpct','10_ALL13-8_001G5_vafpct','10_ARL10-4_001D4-B6_vafpct','10_ARL7-2_005B2_vafpct',
                                                         '10_LA14-2_002D1-A5_vafpct','10_LAFb16-3_001E11_vafpct','10_PRL3-1_001B8-E10_vafpct','10_PRLM2-3_001D10_vafpct','10_PLL6-3_002A11_vafpct')) %>% as.data.frame() %>% column_to_rownames('snv_id') %>% as.matrix()
colnames(hm_mx) <- gsub('_vafpct','',colnames(hm_mx))
library(ComplexHeatmap)
library(circlize)
dim(hm_mx)
rt_annot = rowAnnotation(blood = f_dt$`10_Blood_3s_merged_vafpct`, mVAF=f_dt$meanVAF)
pdf('/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution/DB10_HM_for_tree_two_samples_200108.pdf',height=30, width=10)
Heatmap(hm_mx, col = colorRamp2(c(0,50,100), c("blue","yellow","red")),
        right_annotation = rt_annot, row_names_gp=gpar(fontsize=7), column_names_gp=gpar(fontsize=12)
        ,cell_fun = function(j,i,x,y,width,height,fill){
          grid.rect(x=x, y=y, width=width, height=height,gp = gpar(col="gray",fill=NA))
        }
)
dev.off()


dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution/DB10_50s_merged.txt.sampleN.loh_id.g_fi.fi')
colnames(dt)
dt <- dt %>% mutate(snv_id = paste0(`#CHROM`,':',POS,REF,'>',ALT))
dt <- dt %>% mutate(Nr1=as.integer(str_split_fixed(`dp;Nr1;Nr2;Nr3;Nr3ids`, ';',5)[,2]), Nr3=as.integer(str_split_fixed(`dp;Nr1;Nr2;Nr3;Nr3ids`, ';',5)[,4]), Nr3ids =str_split_fixed(`dp;Nr1;Nr2;Nr3;Nr3ids`, ';',5)[,5])
dt %>% group_by(Nr3ids) %>% count() %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution/Nr3ids_count.txt')

#Extract coverage from sequenza
file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/db10/06_sequenza/', pattern = 'extract.RData$', recursive = T, full.names = T)
dt <- tibble(sample_id = str_split_fixed(str_split_fixed(str_split_fixed(file_list,pattern ='//', n=2)[,2], pattern = '/', n=2)[,1],'\\.',2)[,1] , avg_t_dp = unlist(lapply(file_list, function(x) get(load(x))$avg.depth.tumor)), avg_n_dp = unlist(lapply(file_list, function(x) get(load(x))$avg.depth.normal)))
mean(dt$avg_n_dp)
dt %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10/06_sequenza/cov_from_seqz.txt')

# phylogenetic tree
m_meta_dt <- meta_dt %>% rename(taxa = lineage_id) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
m_meta_dt <- m_meta_dt %>% mutate(side_source2 = paste(Source_side, Source2, sep='_'))
m_meta_dt$side_source2[m_meta_dt$side_source2 == 'NA_umbilicus'] <- 'umbilicus'
m_meta_dt$side_source2[m_meta_dt$side_source2 == 'NA_abdomen'] <- 'abdomen'
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
ggtree(tree) %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id), size=3) + coord_cartesian(xlim = c(0, 34000))
ggtree(tree) %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id), hjust=0) + coord_cartesian(xlim = c(0,30))+
  scale_x_continuous(breaks=seq(0,30,5))+
  theme(axis.text = element_text(size=20))




m_meta_dt <- meta_dt %>% rename(taxa = lineage_id) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
m_meta_dt <- m_meta_dt %>% mutate(side_source2 = paste(Source_side, Source2, sep='_'))
m_meta_dt$side_source2[m_meta_dt$side_source2 == 'NA_occiput'] <- 'occiput'
unique(m_meta_dt$side_source2)


m_lng_dt <- lng_dt %>% filter(deadbody == 'DB5') %>% rename(taxa = lineage_id) %>% select(-deadbody)
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
ggtree(tree) %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id), size=3) + coord_cartesian(xlim = c(0, 8000))
ggtree(tree) %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id), hjust=0) + coord_cartesian(xlim = c(0,30))+
  scale_x_continuous(breaks=seq(0,30,5))+
  theme(axis.text = element_text(size=20))



m_meta_dt <- meta_dt %>% rename(taxa = lineage_id) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
m_meta_dt <- m_meta_dt %>% mutate(source_side_name = paste(Source_side, Source_class, sep='_'))
m_meta_dt$source_side_name[m_meta_dt$source_side_name=='rt_HN'] <- 'HN'
m_meta_dt$source_side_name[m_meta_dt$source_side_name=='lt_HN'] <- 'HN'
m_meta_dt$source_side_name[m_meta_dt$source_side_name=='NA_HN'] <- 'HN'
m_meta_dt$source_side_name[m_meta_dt$source_side_name=='rt_trunk'] <- 'trunk'
m_meta_dt$source_side_name[m_meta_dt$source_side_name=='lt_trunk'] <- 'trunk'
unique(m_meta_dt$source_side_name)

m_lng_dt <- lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
m_meta_dt <- left_join(m_meta_dt, m_lng_dt, by='taxa')

#default trees
p %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id), size=3) + scale_x_continuous(limits = c(0, 32000))
p %<+% m_meta_dt + geom_tree()+theme_tree2() +geom_tiplab(aes(label = sample_id), hjust=0) + coord_cartesian(xlim = c(0,30))+
  scale_x_continuous(breaks=seq(0,30,2))+
  theme(axis.text = element_text(size=20))


# draw seqz clean plot
library(cowplot)
cyto_dt <- read_tsv('/home/users/sypark/03_Tools/annovar/humandb/hg19_cytoBand.txt', col_names = c('chr','start','end','cyto_id','stain'))
acen_dt <- cyto_dt %>% filter(stain == 'acen') 
acen_dt2 <- acen_dt %>% select(start) %>% aggregate(list(acen_dt$chr), min) %>% as.tibble()
acen_dt2$end <- acen_dt %>% select(end) %>% aggregate(list(acen_dt$chr), max) %>% .$end
acen_dt2 <- acen_dt2 %>% mutate(mean = (start+end)/2, chrom = gsub('chr','',Group.1))

file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/db10/06_sequenza/seqz_clean/', pattern = 'clean$', full.names = T)
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
  png(paste0(fn, '.png'), width=18, height=3, units = "in", res=100)
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
file_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution/perSample/',pattern = 'readc$', full.names = T)
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
res_tbl %>% write_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution/perSample/peakVAF.txt')
length(p)
png('/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution/perSample/peakVAF_plot4.png', width=14, height=7, res=150, units="in")
plot_grid(plotlist=p[37:50], ncol=4, align="hv")
dev.off()


##### Serial sample mutCN plotting
library(mgsub)
library(gatepoints)
#load input
mutcn_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution/perSample/culture_associated_mutations/10_LA14-2.mother_daughter_union_loci.txt.readc.readc.seqzcn.scF.scF', col_types = cols(`#CHROM`='c'))
sampleid = "10_LA14-2"
mt_nsubs = 9130
dt_nsubs = 9873
#load input2
mutcn_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution/perSample/culture_associated_mutations/10_ARL10-4.mother_daughter_union_loci.txt.readc.readc.seqzcn.scF.scF',col_types = cols(`#CHROM`='c'))
sampleid = "10_ARL10-4"
mt_nsubs = 5614
dt_nsubs = 5908
#load input3
mutcn_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution/perSample/culture_associated_mutations/10_PRL3-1.mother_daughter_union_loci.txt.readc.readc.seqzcn.scF.scF',col_types = cols(`#CHROM`='c'))
sampleid = "10_PRL3-1"
mt_nsubs = 28445
dt_nsubs = 29726
#plotting
mutcn_dt <- mutcn_dt %>% mutate(context = mgsub(paste0(REF,ALT), c('GT','GC','GA','AT','AG','AC'),c('CA','CG','CT','TA','TC','TG'))) 
context_pal = c("#4393c3","black","#d73027","gray","#a6d96a","#f1b6da")
names(context_pal) = c('CA','CG','CT','TA','TC','TG')
#dt_mutCN >= mt_mutCN*1+0.3 
#dt_mutCN < mt_mutCN*(-1)+1.7
#mt_mutCN < 0.3
ca_dt <- mutcn_dt %>% filter(mt_mutCN+0.3 <= dt_mutCN & mt_mutCN < 0.3 & mt_mutCN*(-1)+1.7 > dt_mutCN)
dim(ca_dt)
ggplot(mutcn_dt, aes(x=mt_mutCN, y=dt_mutCN))+
  geom_point(aes(color=context), size=2, alpha=0.3)+
  geom_segment(aes(x=0, y=0.3, xend=0.3, yend=0.7))+
  geom_segment(aes(x=0.3, y=0.7, xend=0.3, yend=1.3))+
  geom_segment(aes(x=0.3, y=1.3, xend=0, yend=1.7))+
  scale_color_manual(values = context_pal)+
  scale_x_continuous(limits = c(0,2.4), breaks=seq(0,2.4,0.2))+
  scale_y_continuous(limits= c(0,2.4), breaks=seq(0,2.4,0.2))+
  ggtitle(paste0(sampleid, ' (total_mt = ',mt_nsubs,', total_dt = ',dt_nsubs,', ca_mt = ',nrow(ca_dt),')'))+
  theme_bw()
ggplot(ca_dt, aes(x=mt_mutCN, y=dt_mutCN))+
  geom_point(aes(color=context), size=2, alpha=0.3)+
  scale_color_manual(values = context_pal)+
  scale_x_continuous(limits = c(0,2), breaks=seq(0,2,0.2))+
  scale_y_continuous(limits = c(0,2), breaks=seq(0,2,0.2))+
  ggtitle(sampleid)+
  theme_bw()
ca_dt %>% write_tsv(paste0('/home/users/sypark/00_Project/06_LineageTracing/db10/05_substitution/perSample/culture_associated_mutations/',sampleid,'.culture_associated_mutations.txt'))

#sig5-sig7 correlation
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Lineage_tracing_summary_200130.txt') %>% filter(deadbody == 'DB10')
colnames(meta_dt)
meta_dt$Memo
m_meta_dt <- meta_dt %>% filter(Memo !='other_person' | is.na(Memo) == T)
m_meta_dt <- meta_dt %>% filter(current_final_for_lineage == 'Y')
dim(m_meta_dt)
library(ggrepel)
ggplot(m_meta_dt, aes(x=seqz_cov, y=Sig1_ct+corSig5_ct_by0.29))+
  geom_point(aes(color=peakVAF, shape=current_final_for_lineage))+
  geom_text_repel(aes(label=Memo))+
  scale_x_continuous(limits=c(15,35))
