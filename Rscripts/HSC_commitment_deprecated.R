
###########################<Common processing>##############################################
#load library
library(tidyverse)
library(ggsci)
library(ggtree)
library(ggrepel)
library(cowplot)
library(scales)
library(grid)
library(gridExtra)
library(png)

#load dataset
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_200908.txt') %>% filter(current_final_for_lineage == 'Y')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_200908.txt')
DB_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_DB_200918.txt')

# path assign
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

#edit meta_dt
meta_dt <- meta_dt %>% mutate(Source3 = ifelse(is.na(Source_side)==T, Source2, paste0(Source_side,'_',Source2)))
meta_dt$deadbody <- factor(meta_dt$deadbody, levels =c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))

#edit lng_dt
lng_dt$deadbody <- factor(lng_dt$deadbody, levels =c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))

#edit DB_dt
DB_dt$deadbody <- factor(DB_dt$deadbody, levels =c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))

#color setting

#show_col(pal_npg('nrc')(10))
db_pal <- pal_npg('nrc')(10)[c(1:5,9,7)]
names(db_pal) <- unique(meta_dt$deadbody)
cell_pal <- pal_npg('nrc')(6)
names(cell_pal) <- unique(meta_dt$Cell_type)
source2_pal <- pal_d3('category10')(6)
names(source2_pal) <- c('HN','lt_LE','lt_UE','rt_LE','trunk','rt_UE')
timing_pal <- pal_npg('nrc')(10)[c(1,3,4)]
names(timing_pal) <- c('early','early_to_late','late')
pal_combi <- c(pal_npg("nrc")(10), pal_d3('category10')(10), pal_aaas('default')(10))
sig_pal <- pal_combi[c(9,11,1,12,5,17,13)]
names(sig_pal) <- c('v3_1','v3_5','v3_7a','v3_7b','v3_7c','v3_7d','v3_18')


#Fx source
source('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/merged_analysis_functions.R')

#theme setting
theme_syp = theme_bw() + theme(axis.text = element_text(size=15), axis.title = element_text(size=18), strip.text = element_text(size=15), panel.border = element_blank(), axis.line = element_line(), title = element_text(size=18))

########################
this_db = 'DB9'

Draw_blVAF_on_tree <- function(lng_dt, this_db){   #Fig 3d
  library(colorRamps)
  total_break=100
  col_breaks = seq(0, 100, length.out=total_break)
  a <- length(col_breaks[col_breaks < 1] )  #change color at 1
  b <- length(col_breaks[col_breaks < 10])-a  #change color at 10
  rampCol1 <- colorRampPalette(c("gray", "yellow"))(n=a)
  rampCol2 <- colorRampPalette(c("yellow","red"))(n=b)
  rampCol3 <- colorRampPalette(c("red", "blue"))(n=total_break-(a+b))
  my_palette<- c(rampCol1, rampCol2, rampCol3)
  
  file_list <- list.files(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db], 'perLineage/'), pattern='^L.*.pointmts.txt$', full.names = T)
  dt <- do.call(rbind,lapply(file_list, function(x) read_tsv(x, col_types = cols(`#CHROM`='c')) %>% mutate(lineage_id2 = paste0(this_db, '_', gsub('.pointmts.txt','',unlist(strsplit(x,'//'))[2])))))
  dt <- dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep="_")) %>% separate(`blood_dp;blood_var;blood_vafpct`, c('blood_dp','blood_var','blood_vafpct'), sep=';',convert=T)
  dt <- left_join(dt, lng_dt %>% select(lineage_id2, start_mtn, time_of_lineage))
  dt <- dt %>% mutate(adjblVAF = ifelse(time_of_lineage =='early' | blood_var > 1, blood_vafpct,0))   # ignore when the lineage is not early and blood_var == 1
  for (t_lid in unique(dt$lineage_id2)){   #ignore when the lowest 3 blVAFs of the upper stream are 0
    if (dt %>% filter(lineage_id2 == t_lid) %>% arrange(adjblVAF) %>% head(3) %>% pull(adjblVAF) %>% max() == 0){
      dt$adjblVAF[grepl(paste0(t_lid, '-'),dt$lineage_id2)] <- 0
    }
  }
  max_VAF_len=20
  res_tbl <- as.tibble(data.frame(matrix(vector(), 0, max_VAF_len + 1,
                                         dimnames=list(c(), c('lineage_id2',paste0('VAF',1:max_VAF_len)))),
                                  stringsAsFactors=F))
  n=0
  for (lid in unique(dt$lineage_id2)){
    n=n+1
    VAF_v <- dt %>% filter(lineage_id2 == lid) %>% arrange(desc(adjblVAF)) %>% pull(adjblVAF)%>% .[1:max_VAF_len]
    NA_v <- rep(NA, max_VAF_len - length(VAF_v))
    res_tbl[n,] <- c(lid,VAF_v, NA_v)
  }
  res_tbl <- res_tbl %>% mutate_at(vars(starts_with('VAF')), list(as.numeric))
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  m_lng_dt <- lng_dt %>% filter(deadbody == this_db) %>% rename(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt <- left_join(m_lng_dt, res_tbl)
  
  ggtree(tree) %<+% m_lng_dt + theme_tree2()+
    coord_cartesian(xlim = c(-1,35))+
    geom_text2(aes(subset = is.na(Source_class2) == F, x=x-n_pointmt+3, y = y+0.25, label = Source_class2), size=2.5, color="red")+
    geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y, color = VAF1))+
    geom_text2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, y = y, label = VAF1), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2))+
    geom_text2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, y = y, label = VAF2), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3))+
    geom_text2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, y = y, label = VAF3), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4))+
    geom_text2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, y = y, label = VAF4), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, xend = x -n_pointmt +5, y = y, yend = y, color = VAF5))+
    geom_text2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, y = y, label = VAF5), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, xend = x -n_pointmt +6, y = y, yend = y, color = VAF6))+
    geom_text2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, y = y, label = VAF6), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, xend = x -n_pointmt +7, y = y, yend = y, color = VAF7))+
    geom_text2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, y = y, label = VAF7), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, xend = x -n_pointmt +8, y = y, yend = y, color = VAF8))+
    geom_text2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, y = y, label = VAF8), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, xend = x -n_pointmt +9, y = y, yend = y, color = VAF9))+
    geom_text2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, y = y, label = VAF9), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, xend = x -n_pointmt +10, y = y, yend = y, color = VAF10))+
    geom_text2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, y = y, label = VAF10), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, xend = x -n_pointmt +11, y = y, yend = y, color = VAF11))+
    geom_text2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, y = y, label = VAF11), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, xend = x -n_pointmt +12, y = y, yend = y, color = VAF12))+
    geom_text2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, y = y, label = VAF12), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, xend = x -n_pointmt +13, y = y, yend = y, color = VAF13))+
    geom_text2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, y = y, label = VAF13), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, xend = x -n_pointmt +14, y = y, yend = y, color = VAF14))+
    geom_text2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, y = y, label = VAF14), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, xend = x -n_pointmt +15, y = y, yend = y, color = VAF15))+
    geom_text2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, y = y, label = VAF15), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, xend = x -n_pointmt +16, y = y, yend = y, color = VAF16))+
    geom_text2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, y = y, label = VAF16), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, xend = x -n_pointmt +17, y = y, yend = y, color = VAF17))+
    geom_text2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, y = y, label = VAF17), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, xend = x -n_pointmt +18, y = y, yend = y, color = VAF18))+
    geom_text2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, y = y, label = VAF18), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, xend = x -n_pointmt +19, y = y, yend = y, color = VAF19))+
    geom_text2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, y = y, label = VAF19), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, xend = x -n_pointmt +20, y = y, yend = y, color = VAF20))+
    geom_text2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, y = y, label = VAF20), size=2.5)+
    scale_color_gradientn(colors = my_palette, breaks = col_breaks)+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
}





Draw_blVAF_on_tree_simple <- function(lng_dt, this_db){
  library(colorRamps)
  total_break=1000
  col_breaks = seq(0, 100, length.out=total_break)
  brk1=1
  brk2=10
  a <- length(col_breaks[col_breaks < brk1] )  #change color at 1
  b <- length(col_breaks[col_breaks < brk2])-a  #change color at 10
  rampCol1 <- colorRampPalette(c("gray", "yellow"))(n=a)
  rampCol2 <- colorRampPalette(c("yellow","red"))(n=b)
  rampCol3 <- colorRampPalette(c("red", "blue"))(n=total_break-(a+b))
  my_palette<- c(rampCol1, rampCol2, rampCol3)
  
  colfunc <- colorRampPalette(c("red", "gray"))
  #plot(1:a, 1:a, pch = 19, cex=2, col = rev(rampCol1))
  legend_image <- as.raster(matrix(rev(rampCol1), ncol=1))
  plot(c(0,2),c(0,50),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
  rasterImage(legend_image, 0, 0, 0.3,brk1)
  legend_image <- as.raster(matrix(rev(rampCol2), ncol=1))
  rasterImage(legend_image, 0, brk1, 0.3,brk2)
  legend_image <- as.raster(matrix(rev(rampCol3), ncol=1))
  rasterImage(legend_image, 0, brk2, 0.3,50)
  text(x=0.5, y = seq(0,50,10), labels = paste0(seq(0,50,10),'%'), cex=1)
  
  file_list <- list.files(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db], 'perLineage/'), pattern='^L.*.pointmts.txt$', full.names = T)
  dt <- do.call(rbind,lapply(file_list, function(x) read_tsv(x, col_types = cols(`#CHROM`='c')) %>% mutate(lineage_id2 = paste0(this_db, '_', gsub('.pointmts.txt','',unlist(strsplit(x,'//'))[2])))))
  dt <- dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep="_")) %>% separate(`blood_dp;blood_var;blood_vafpct`, c('blood_dp','blood_var','blood_vafpct'), sep=';',convert=T)
  dt <- left_join(dt, lng_dt %>% select(lineage_id2, start_mtn, time_of_lineage))
  dt <- dt %>% mutate(adjblVAF = ifelse(time_of_lineage =='early' | blood_var > 1, blood_vafpct,0))   # ignore when the lineage is not early and blood_var == 1
  for (t_lid in unique(dt$lineage_id2)){   #ignore when the lowest 3 blVAFs of the upper stream are 0
    if (dt %>% filter(lineage_id2 == t_lid) %>% arrange(adjblVAF) %>% head(3) %>% pull(adjblVAF) %>% max() == 0){
      dt$adjblVAF[grepl(paste0(t_lid, '-'),dt$lineage_id2)] <- 0
    }
  }
  max_VAF_len=20
  res_tbl <- as.tibble(data.frame(matrix(vector(), 0, max_VAF_len + 1,
                                         dimnames=list(c(), c('lineage_id2',paste0('VAF',1:max_VAF_len)))),
                                  stringsAsFactors=F))
  n=0
  for (lid in unique(dt$lineage_id2)){
    n=n+1
    VAF_v <- dt %>% filter(lineage_id2 == lid) %>% arrange(desc(adjblVAF)) %>% pull(adjblVAF)%>% .[1:max_VAF_len]
    NA_v <- rep(NA, max_VAF_len - length(VAF_v))
    res_tbl[n,] <- c(lid,VAF_v, NA_v)
  }
  res_tbl <- res_tbl %>% mutate_at(vars(starts_with('VAF')), list(as.numeric))
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  m_lng_dt <- lng_dt %>% filter(deadbody == this_db) %>% rename(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt <- left_join(m_lng_dt, res_tbl)
  
  ggtree(tree) %<+% m_lng_dt + theme_tree2()+
    coord_cartesian(xlim = c(-1,35))+
    geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y, color = VAF1), size=3)+
    geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2), size=3)+
    geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3), size=3)+
    geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4), size=3)+
    geom_segment2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, xend = x -n_pointmt +5, y = y, yend = y, color = VAF5), size=3)+
    geom_segment2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, xend = x -n_pointmt +6, y = y, yend = y, color = VAF6), size=3)+
    geom_segment2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, xend = x -n_pointmt +7, y = y, yend = y, color = VAF7), size=3)+
    geom_segment2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, xend = x -n_pointmt +8, y = y, yend = y, color = VAF8), size=3)+
    geom_segment2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, xend = x -n_pointmt +9, y = y, yend = y, color = VAF9), size=3)+
    geom_segment2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, xend = x -n_pointmt +10, y = y, yend = y, color = VAF10), size=3)+
    geom_segment2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, xend = x -n_pointmt +11, y = y, yend = y, color = VAF11), size=3)+
    geom_segment2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, xend = x -n_pointmt +12, y = y, yend = y, color = VAF12), size=3)+
    geom_segment2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, xend = x -n_pointmt +13, y = y, yend = y, color = VAF13), size=3)+
    geom_segment2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, xend = x -n_pointmt +14, y = y, yend = y, color = VAF14), size=3)+
    geom_segment2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, xend = x -n_pointmt +15, y = y, yend = y, color = VAF15), size=3)+
    geom_segment2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, xend = x -n_pointmt +16, y = y, yend = y, color = VAF16), size=3)+
    geom_segment2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, xend = x -n_pointmt +17, y = y, yend = y, color = VAF17), size=3)+
    geom_segment2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, xend = x -n_pointmt +18, y = y, yend = y, color = VAF18), size=3)+
    geom_segment2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, xend = x -n_pointmt +19, y = y, yend = y, color = VAF19), size=3)+
    geom_segment2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, xend = x -n_pointmt +20, y = y, yend = y, color = VAF20), size=3)+
    scale_color_gradientn(colors = my_palette, breaks = col_breaks)+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
}

Draw_blVAF_on_tree_median <- function(lng_dt, this_db){
  library(colorRamps)
  total_break=1000
  col_breaks = seq(0, 100, length.out=total_break)
  brk1=1
  brk2=10
  a <- length(col_breaks[col_breaks < brk1] )  #change color at 1
  b <- length(col_breaks[col_breaks < brk2])-a  #change color at 10
  rampCol1 <- colorRampPalette(c("gray", "yellow"))(n=a)
  rampCol2 <- colorRampPalette(c("yellow","red"))(n=b)
  rampCol3 <- colorRampPalette(c("red", "blue"))(n=total_break-(a+b))
  my_palette<- c(rampCol1, rampCol2, rampCol3)
  
  colfunc <- colorRampPalette(c("red", "gray"))
  #plot(1:a, 1:a, pch = 19, cex=2, col = rev(rampCol1))
  legend_image <- as.raster(matrix(rev(rampCol1), ncol=1))
  plot(c(0,2),c(0,50),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
  rasterImage(legend_image, 0, 0, 0.3,brk1)
  legend_image <- as.raster(matrix(rev(rampCol2), ncol=1))
  rasterImage(legend_image, 0, brk1, 0.3,brk2)
  legend_image <- as.raster(matrix(rev(rampCol3), ncol=1))
  rasterImage(legend_image, 0, brk2, 0.3,50)
  text(x=0.5, y = seq(0,50,10), labels = paste0(seq(0,50,10),'%'), cex=1)
  
  file_list <- list.files(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db], 'perLineage/'), pattern='^L.*.pointmts.txt$', full.names = T)
  dt <- do.call(rbind,lapply(file_list, function(x) read_tsv(x, col_types = cols(`#CHROM`='c')) %>% mutate(lineage_id2 = paste0(this_db, '_', gsub('.pointmts.txt','',unlist(strsplit(x,'//'))[2])))))
  dt <- dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep="_")) %>% separate(`blood_dp;blood_var;blood_vafpct`, c('blood_dp','blood_var','blood_vafpct'), sep=';',convert=T)
  dt <- left_join(dt, lng_dt %>% select(lineage_id2, start_mtn, time_of_lineage))
  for (t_lid in unique(dt$lineage_id2)){   #ignore when the lowest 3 blVAFs of the upper stream are 0
    if (dt %>% filter(lineage_id2 == t_lid) %>% arrange(blood_vafpct) %>% head(3) %>% pull(blood_vafpct) %>% max() == 0){
      dt$blood_vafpct[grepl(paste0(t_lid, '-'),dt$lineage_id2)] <- 0
    }
  }
  res_tbl <-tibble(lineage_id2 = as.character(), medBlVAF = as.numeric(), time_of_lineage = as.character())
  n=0
  for (t_lid in unique(dt$lineage_id2[dt$time_of_lineage == 'early_to_late'])){
    n=n+1
    med_v <- dt %>% filter(lineage_id2 == t_lid) %>% arrange(desc(blood_vafpct)) %>% head(2) %>% pull(blood_vafpct) %>% median()
    res_tbl[n,] <- list(t_lid, med_v, 'early_to_late')
  }
  res_tbl <- bind_rows(res_tbl, dt %>% filter(time_of_lineage == 'early') %>% group_by(lineage_id2, time_of_lineage) %>% summarise(medBlVAF = median(blood_vafpct)))
  res_tbl <- res_tbl %>% mutate_at(vars(starts_with('VAF')), list(as.numeric))
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  m_lng_dt <- lng_dt %>% filter(deadbody == this_db) %>% rename(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt <- left_join(m_lng_dt, res_tbl)
  ggtree(tree) %<+% m_lng_dt + theme_tree2()+
    coord_cartesian(xlim = c(-1,35))+
    geom_segment2(aes(subset = (time_of_lineage == 'early'), x=x-n_pointmt, xend = x, y = y, yend = y, color = medBlVAF), size=3)+
    geom_segment2(aes(subset = (time_of_lineage == 'early_to_late'), x=x-n_pointmt, xend = x-n_pointmt+2, y = y, yend = y, color = medBlVAF), size=3)+
    scale_color_gradientn(colors = my_palette, breaks = col_breaks)+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
}


Draw_blVAFprop_on_HoriLine <- function(lng_dt, pointmt_path_tbl, this_db, n_earlymt){
  file_list <- list.files(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db], 'perLineage/'), pattern='^L.*.pointmts.txt$', full.names = T)
  dt <- do.call(rbind,lapply(file_list, function(x) read_tsv(x, col_types = cols(`#CHROM`='c')) %>% mutate(lineage_id2 = paste0(this_db, '_', gsub('.pointmts.txt','',unlist(strsplit(x,'//'))[2])))))
  dt <- dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep="_")) %>% separate(`blood_dp;blood_var;blood_vafpct`, c('blood_dp','blood_var','blood_vafpct'), sep=';',convert=T)
  dt <- left_join(dt, lng_dt %>% select(lineage_id2, start_mtn, time_of_lineage))
  dt <- dt %>% mutate(adjblVAF = ifelse(time_of_lineage =='early' | blood_var > 1, blood_vafpct,0))   # ignore when the lineage is not early and blood_var == 1
  for (t_lid in unique(dt$lineage_id2)){   #ignore when the lowest 3 blVAFs of the upper stream are 0
    if (dt %>% filter(lineage_id2 == t_lid) %>% arrange(adjblVAF) %>% head(3) %>% pull(adjblVAF) %>% max() == 0){
      dt$adjblVAF[grepl(paste0(t_lid, '-'),dt$lineage_id2)] <- 0
    }
  }
  dt <- dt %>% filter(time_of_lineage != 'late')
  tmp_dt <- dt %>% group_by(lineage_id2) %>% mutate(rank = rank(adjblVAF, ties.method = "random")) %>% mutate(exp_mtn = start_mtn + rank)
  f_tmp_dt <- tmp_dt %>% filter(exp_mtn == n_earlymt) %>% mutate(adjblVAFprop = adjblVAF/sum(.$adjblVAF)) %>% mutate(exp_prop = 1/nrow(.))
  f_tmp_dt <- f_tmp_dt %>% select(lineage_id2, adjblVAFprop, exp_prop) %>% gather(-lineage_id2, key="group", value="prop")
  f_tmp_dt$lineage_id2 <- factor(f_tmp_dt$lineage_id2, levels = f_tmp_dt %>% filter(group == 'adjblVAFprop') %>% arrange(prop) %>% pull(lineage_id2))
  f_tmp_dt$group <- factor(f_tmp_dt$group, levels = c('exp_prop', 'adjblVAFprop'))
  ggplot(f_tmp_dt, aes(x=group, y=prop, fill = lineage_id2))+ 
    geom_bar(stat="identity")+
    ggtitle(paste0("No. Mts = ",n_earlymt))
}



Draw_obs_exp_ratio_blVAF <- function(lng_dt, pointmt_path_tbl, this_db, color_method){
  require(binom)
  #add expected VAF
  id_list <- lng_dt %>% filter(time_of_lineage !='late' & deadbody == this_db) %>% pull(lineage_id2)
  res_tbl1 = tibble(lineage_id2 = as.character(), expected_VAF = as.numeric())
  n=0
  for (this_lid in id_list){
    n= n+1
    t_db <- unlist(strsplit(this_lid,'_'))[1]
    lng_dt %>% filter(time_of_lineage =='early_to_late' & grepl(this_lid, lineage_id2)) %>% nrow() -> a
    lng_dt %>% filter(time_of_lineage =='early_to_late' & deadbody == t_db ) %>% nrow() -> b
    lower_margin = 0.48*binom.confint(a,b, method="exact")$lower
    upper_margin = 0.48*binom.confint(a,b, method="exact")$upper
    expected_value = 0.48*a/b
    res_tbl1[n,] <- list(this_lid, expected_value)
  }
  m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db), res_tbl1)
  #add median value of first two blVAF
  dt <- tibble(path = as.character(), deadbody = as.character())
  file_list <- list.files(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db], 'perLineage/'), pattern='^L.*.pointmts.txt$', full.names = T)
  dt <- tibble(path= file_list, deadbody = this_db)
  dt <- dt %>% mutate(lineage_id2 = apply(.[c("path","deadbody")], 1, function(x) paste0(x[2], '_', gsub('.pointmts.txt','',unlist(strsplit(x[1],'//')))[2])))
  el_list <- lng_dt %>% filter(time_of_lineage == 'early_to_late') %>% pull(lineage_id2)  
  dt1 <- dt %>% filter(lineage_id2 %in% el_list)
  twoMedVAF <- c()
  for (t_path in dt1$path){  #calculta median value of the two largest VAFs in early_to_late lineages
    print(t_path)
    t_dt <- read_tsv(t_path, col_types = cols(`#CHROM`='c'))
    t_dt <- t_dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep="_")) %>% separate(`blood_dp;blood_var;blood_vafpct`, c('blood_dp','blood_var','blood_vafpct'), sep=';',convert=T) %>%
      filter(blood_dp >= 30) %>% arrange(desc(blood_vafpct))
    t_dt$blood_vafpct[t_dt$blood_var ==1] <- 0
    twoMedVAF <- c(twoMedVAF, t_dt %>% pull(blood_vafpct) %>% .[1:2] %>% median())
  }
  dt1$twoMedVAF <- twoMedVAF
  e_list <- lng_dt %>% filter(time_of_lineage == 'early') %>% pull(lineage_id2)
  dt2 <- dt %>% filter(lineage_id2 %in% e_list)
  for (t_path in dt2$path){  # remove early_to_late lineages if upstreatm lineage have three sequential 0 VAFs
    print(t_path)
    t_lid <- paste0(this_db,'_',gsub('.pointmts.txt','',unlist(strsplit(t_path,'//'))[2]))
    t_dt <- read_tsv(t_path, col_types = cols(`#CHROM`='c'))
    l3_v <- t_dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep="_")) %>% separate(`blood_dp;blood_var;blood_vafpct`, c('blood_dp','blood_var','blood_vafpct'), sep=';',convert=T) %>%
      arrange(blood_vafpct) %>% head(3) %>%  pull(blood_vafpct) 
    if (max(l3_v) == 0){
      dt1$twoMedVAF[grepl(paste0(t_lid, '-'),dt1$lineage_id2)] <- 0
    }
  }
  m_lng_dt <- left_join(m_lng_dt, dt1 %>% select(lineage_id2, twoMedVAF))
  m_lng_dt <- m_lng_dt %>% filter(time_of_lineage != 'late' & lineage_id != 'L0') %>% mutate(medVAF_v = ifelse(time_of_lineage == 'early', median_blood_VAFpct, twoMedVAF)) %>%
    mutate(m_end_mtn = ifelse(time_of_lineage == 'early', end_mtn, start_mtn+2))
  m_lng_dt <- m_lng_dt %>% mutate(first_group = apply(.["lineage_id"], 1, function(x) unlist(strsplit(x, '-'))[1]))
  m_lng_dt <- m_lng_dt %>% mutate(obs_exp_ratio =medVAF_v/expected_VAF/100 )
  con_tbl <- tibble(start_lid = as.character(), end_lid = as.character(), x = as.numeric(), xend = as.numeric(), y=as.numeric(), yend=as.numeric())
  for (lid in m_lng_dt$lineage_id2){
    con_tbl <- bind_rows(con_tbl, m_lng_dt %>% filter(grepl(paste0(lid,'-'), lineage_id2) & step_n == m_lng_dt$step_n[m_lng_dt$lineage_id2 == lid]+1 ) %>% select(lineage_id2, start_mtn, obs_exp_ratio) %>% rename(end_lid = lineage_id2, xend = start_mtn, yend = obs_exp_ratio) %>%
                           mutate(start_lid = lid, x = m_lng_dt$start_mtn[m_lng_dt$lineage_id2 == lid], y = m_lng_dt$obs_exp_ratio[m_lng_dt$lineage_id2 == lid]))
  }
  con_tbl <- con_tbl %>% mutate(group1 = ifelse(yend > y, 'increasing', 'decreasing'))
  max_lid <- con_tbl$end_lid[con_tbl$yend == max(con_tbl$yend)]
  max_lid_list <- unlist(strsplit(max_lid,'-'))
  max_lid_dts <- c()
  for (i in 1:length(max_lid_list)){
    max_lid_dts <- c(max_lid_dts,print(paste(max_lid_list[1:i],collapse='-')))
  }
  con_tbl <- con_tbl %>% mutate(group2 = ifelse(end_lid %in% max_lid_dts, 'max_lineage','other'))
  if (color_method == 'increasing'){
    con_tbl <- con_tbl %>% mutate(color_group = group1)
  } else if (color_method == 'max'){
    con_tbl <- con_tbl %>% mutate(color_group = group2)
  }
  ggplot(m_lng_dt)+
    geom_point(aes(x=start_mtn, y=medVAF_v/expected_VAF/100))+
    geom_segment(data = con_tbl, aes(x=x, xend=xend, y=y, yend=yend, color = color_group))+
    geom_hline(yintercept = 1, linetype = 'dashed', color="blue")+
    scale_color_manual(values = c("increasing" = "red", "decreasing" = "gray", "max_lineage" = "red", "other" = "gray"))+
    coord_cartesian(xlim=c(0,35))+
    ylab("observedVAF/expectedVAF")+xlab("No. of early mutations")+
    ggtitle(this_db)+
    theme_syp
}




#blood VAF on phylogenetic tree
if(F){
  this_db = 'DB9'
  Draw_blVAF_on_tree(lng_dt, this_db)
  Draw_blVAF_on_tree_simple(lng_dt, this_db)
  Draw_blVAF_on_tree_median(lng_dt, this_db)
}

#blood VAF composition change barplot
if(F){
  this_db='DB3'
  n_earlymt=9 # mutation number for proportion (x axis in phylogenetic tree)
  Draw_blVAFprop_on_HoriLine(lng_dt, pointmt_path_tbl, this_db, n_earlymt)
}

#blood observedmedianVAF/expected VAF ratio
if(F){
  this_db = 'DB10'
  color_method="max"  #c("increasing", "max")
  Draw_obs_exp_ratio_blVAF(lng_dt, pointmt_path_tbl, this_db, color_method)
}

#blood specific mutation count
blct_tbl <- tribble(~deadbody, ~numb,~peakVAF,
                    'DB2',4,0.14,
                    'DB3',0,0,
                    'DB5',3,0.14,
                    'DB6',16,0.12,
                    'DB8',4,0.1,
                    'DB9',21,0.18,
                    'DB10',236,0.12)
blct_tbl$deadbody = factor(blct_tbl$deadbody, levels = c('DB2','DB3','DB5','DB6','DB8','DB9','DB10'))
ggplot(blct_tbl, aes(x=numb, y=peakVAF))+
  geom_point(aes(color=deadbody), size=5, alpha=0.5)+
  xlab("No. of blood-specific mutations occured before clonal expansion")+
  ylab("peak VAF of blood-specific mutations")+
  scale_color_manual(values = db_pal)+
  theme_bw()+theme(axis.text = element_text(size=12), axis.title = element_text(size=15))


#DB10 VAF distribution
if(F){
  ch_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/clonal_hematopoiesis/CH_final_list.txt', col_types = cols(`#CHROM`='c'))
  ch_dt <- ch_dt %>% mutate(context = mgsub(paste0(REF,ALT), c("GT","GC","GA","AT","AG","AC"), c("CA","CG","CT","TA","TC","TG")))
  sub_dt <- ch_dt %>% filter(deadbody == 'DB10') %>% mutate(DP = ref_readN + var_readN) #DB10 female -> no adjustment
  colnames(sub_dt)
  sub_dt$VAF %>% summary()
  ggplot(sub_dt)+
    geom_histogram(aes(x=VAF, fill=context))+
    geom_vline(xintercept=0.11)+
    scale_x_continuous(limits=c(0,0.4), breaks=seq(0,0.4,0.1))
  
  ggplot(sub_dt, aes(x=DP, y=VAF, color=context))+
    geom_point()
  
  DP_list <- sub_dt  %>% pull(DP)
  merged_v <- c()
  for (i in 1:1000){
    print(i)
    merged_v = c(merged_v, lapply(DP_list, function(x) { v = rbinom(1, x, 0.11); if (v >=3){v/x}}) %>% unlist())
  }
  merged_v %>% as.tibble() %>% ggplot()+
    geom_histogram(aes(value))+
    geom_vline(xintercept=0.11)+
    scale_x_continuous(limits=c(0,0.4), breaks=seq(0,0.4,0.1))+
    ggtitle("rbinom(DP,0.11) simulation *1000")
}


#DB9 VAF distribution
if(F){
  ch_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/clonal_hematopoiesis/CH_final_list.txt', col_types = cols(`#CHROM`='c'))
  ch_dt <- ch_dt %>% mutate(context = mgsub(paste0(REF,ALT), c("GT","GC","GA","AT","AG","AC"), c("CA","CG","CT","TA","TC","TG")))
  sub_dt <- ch_dt %>% filter(deadbody == 'DB9') %>% mutate(DP = ref_readN + var_readN) #DB10 female -> no adjustment
  colnames(sub_dt)
  sub_dt$VAF %>% summary()
  ggplot(sub_dt)+
    geom_histogram(aes(x=VAF, fill=context))+
    geom_vline(xintercept=0.18)+
    geom_vline(xintercept=0.3)+
    scale_x_continuous(limits=c(0,0.4), breaks=seq(0,0.4,0.1))
  
  ggplot(sub_dt, aes(x=DP, y=VAF, color=context))+
    geom_point()
  DP_list <- sub_dt  %>% pull(DP)
  merged_v=c()
  for (i in 1:1000){
    print(i)
    merged_v = c(merged_v, lapply(DP_list, function(x) { v = rbinom(1, x, 0.18); if (v >=3){v/x}}) %>% unlist())
  }
  merged_v %>% as.tibble() %>% ggplot()+
    geom_histogram(aes(value))+
    geom_vline(xintercept=0.18)+
    scale_x_continuous(limits=c(0,0.4), breaks=seq(0,0.4,0.1))+
    ggtitle("rbinom(DP,0.18) simulation *1000")
  
  merged_v1=c()
  merged_v2=c()
  subDP1 <- sample(DP_list, 11)
  subDP2 <- setdiff(DP_list, subDP1)
  for (i in 1:1000){
    print(i)
    merged_v1 = c(merged_v1, lapply(subDP1, function(x) { v = rbinom(1, x, 0.18); if (v >=3){v/x}}) %>% unlist())
    merged_v2 = c(merged_v2, lapply(subDP2, function(x) { v = rbinom(1, x, 0.3); if (v >=3){v/x}}) %>% unlist())
  }
  merged_v <- c(merged_v1, merged_v2)
  merged_v %>% as.tibble() %>% ggplot()+
    geom_histogram(aes(value))+
    geom_vline(xintercept=0.18)+
    geom_vline(xintercept=0.3)+
    scale_x_continuous(limits=c(0,0.4), breaks=seq(0,0.4,0.1))+
    ggtitle("rbinom(DP,0.18, 0.3) simulation *1000")
}





