#Draw absCN plot
setwd("~/00_Project/06_LineageTracing/db1/11_smoothenedCN/absCN_withDB3LLL-1Fb001G5/")
file_list <- list.files(pattern = "absCN$")
pdf('DB1_AbsCN_vsDB3.genomewide.190225.pdf', height=1, width=30)
for (n in 1:length(file_list)){
  print(file_list[n])
  sampleid<-unlist(strsplit(file_list[n],'\\.'))[1]
  print(sampleid)
  dt <- read_tsv(file_list[n], comment = '#', col_names=F, col_types=cols(X1="c"))
  colnames(dt) <- c('chr1','pos1','Adp','Bdp','absCN')
  chr_list <- c(1:22, "X","Y")
  p<-list()
  i=0
  for (chrom in chr_list){
    i=i+1
    subdt <- dt %>% filter(chr1 == chrom)
    #gender correction for male
    if (chrom =='X' | chrom == 'Y'){
      subdt$absCN <- subdt$absCN/2
    }
    p[[i]] <- ggplot(subdt, aes(x=pos1, y=absCN))+
      geom_point(alpha=0.5)+
      ggtitle(paste0('chr',chrom))+
      scale_y_continuous(limits=c(0,5), breaks=seq(0,5,1))+
      xlab(sampleid)+
      theme(axis.text.x=element_text(angle=45))
  }
  pl<-plot_grid(plotlist=p, nrow=1, align="h")
  print(pl)
}
dev.off()


#Draw vaf heatmap using sampled subs
setwd("~/00_Project/06_LineageTracing/db1/04-2_subs_re")
library(tidyverse)
dt <- read_tsv('DB1_sampled_subs.txt.9sCall')
multiclonals <- c('1_Lt_Wh_Fb_001E11_vafpct')
bloodsamples <-c('1_blood-1_vafpct')
colnames(dt)
dt <- dt %>% select(`#CHROM`,POS, ID, REF, ALT, ends_with('vafpct'))
dt <- dt %>% mutate(snv_id = paste(`#CHROM`,POS,REF,ALT, sep="_")) %>% select(snv_id, ends_with('vafpct'))
dt <- dt %>% select(-one_of(multiclonals), -one_of(bloodsamples))
mx <- dt %>% as.data.frame() %>% column_to_rownames('snv_id') %>% as.matrix()
library(ComplexHeatmap)
Heatmap(mx, na_col = "black", row_names_gp=gpar(fontsize=4))

tmp_mx <- dt %>% filter(`1_Lt_Wh_Fb_001B4_vafpct`==0)%>% as.data.frame() %>% column_to_rownames('snv_id') %>% as.matrix()
Heatmap(tmp_mx, na_col = "black", row_names_gp=gpar(fontsize=8))


#Draw phylogenetic tree using previous nwk
library(tidyverse)
library(ggtree)
setwd("~/00_Project/06_LineageTracing/db1/04-1_subs_previous")
tree <- read.tree('DB1_total_subs.nwk')
ggtree(tree)+theme_tree2()+geom_tiplab()+
  scale_x_continuous(breaks=seq(0,4000,1000), limits=c(0,7000))


