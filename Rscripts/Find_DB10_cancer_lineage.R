# DB10 cancer VAF on phylogenetic tree
library(tidyverse)


#load dataset
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_210322.txt') %>% filter(current_final_for_lineage == 'Y')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_210322.txt')
DB_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_DB_210322.txt')

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


#define functions

Color_VAF_on_phylo_2 <- function(input_tbl, nwk_path){  #required columns in in input_tbl: lineage_id, rank, n_vars_same_rank, vaf, n_pointmt
  #color setting
  total_break=100
  start_value=0
  end_value=100
  change_value1=10
  change_value2=30
  col_breaks = seq(start_value, end_value, length.out=total_break)
  a <- length(col_breaks[col_breaks < change_value1] )  #change color 
  b <- length(col_breaks[col_breaks < change_value2])-a  #change color
  rampCol1 <- colorRampPalette(c("gray", "yellow"))(n=a)
  rampCol2 <- colorRampPalette(c("yellow","red"))(n=b)
  rampCol3 <- colorRampPalette(c("red", "blue"))(n=total_break-(a+b))
  my_palette<- c(rampCol1, rampCol2, rampCol3)
  #print legend
  if(T){
    plot.new()
    plot(c(0,2),c(start_value,end_value),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
    colfunc1 <- rev(rampCol1)
    legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
    rasterImage(legend_image1, 0, start_value, 0.3,change_value1)
    colfunc2 <- rev(rampCol2)
    legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
    rasterImage(legend_image2, 0, change_value1, 0.3,change_value2)
    colfunc3 <- rev(rampCol3)
    legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
    rasterImage(legend_image3, 0, change_value2, 0.3,end_value)
    text(x=0.5, y = c(start_value, change_value1, change_value2, end_value), labels =c(start_value, change_value1, change_value2,  end_value), cex=2)
  }
  
  #make vaf tbl
  input_tbl$vaf <- as.numeric(input_tbl$vaf)
  m_input_tbl <- input_tbl %>% mutate(var_g_id = paste0(lineage_id,'_',rank)) %>% filter(is.na(rank)==F) %>% group_by(lineage_id, rank, n_vars_same_rank) %>%
    summarise(mean_vaf = mean(vaf))
  vaf_split_tbl <- tibble(lineage_id = as.character(), vaf_rank = as.character(), mean_vaf=as.numeric())
  for (i in 1:nrow(m_input_tbl)){
    t_id <- m_input_tbl$lineage_id[i]
    t_rank <- m_input_tbl$rank[i]
    t_same <- m_input_tbl$n_vars_same_rank[i]
    t_vaf <- m_input_tbl$mean_vaf[i]
    vaf_ids <- paste0('VAF',seq(t_rank, t_rank + t_same -1))
    t_tbl <- tibble(lineage_id = t_id, vaf_rank = vaf_ids, mean_vaf = t_vaf)
    t_tbl$mean_vaf <- as.numeric(t_tbl$mean_vaf)
    vaf_split_tbl <- bind_rows(vaf_split_tbl, t_tbl)
  }
  annot_tbl <- vaf_split_tbl %>% spread(key=vaf_rank, value=mean_vaf)
  current_cols <- annot_tbl %>% select(starts_with("VAF")) %>% colnames()
  additional_cols <- setdiff(paste0("VAF", 1:20), current_cols)
  annot_tbl <- cbind(annot_tbl, tibble(additional_cols) %>% mutate(value=NA) %>% spread(key=additional_cols, value=value)) %>% as.tibble()
  annot_tbl <- left_join(annot_tbl, input_tbl %>% select(lineage_id, n_pointmt) %>% unique())
  
  #draw colored phylogenetic tree
  library(ggtree)
  tree <- read.tree(nwk_path)
  g <- ggtree(tree) %<+% annot_tbl + theme_tree2()+
    coord_cartesian(xlim = c(-1,30))+
    geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y,color = VAF1) ,size=2)+
    geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2), size=2)+
    geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3), size=2)+
    geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4), size=2)+
    geom_segment2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, xend = x -n_pointmt +5, y = y, yend = y, color = VAF5), size=2)+
    geom_segment2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, xend = x -n_pointmt +6, y = y, yend = y, color = VAF6), size=2)+
    geom_segment2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, xend = x -n_pointmt +7, y = y, yend = y, color = VAF7), size=2)+
    geom_segment2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, xend = x -n_pointmt +8, y = y, yend = y, color = VAF8), size=2)+
    geom_segment2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, xend = x -n_pointmt +9, y = y, yend = y, color = VAF9), size=2)+
    geom_segment2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, xend = x -n_pointmt +10, y = y, yend = y, color = VAF10), size=2)+
    geom_segment2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, xend = x -n_pointmt +11, y = y, yend = y, color = VAF11), size=2)+
    geom_segment2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, xend = x -n_pointmt +12, y = y, yend = y, color = VAF12), size=2)+
    geom_segment2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, xend = x -n_pointmt +13, y = y, yend = y, color = VAF13), size=2)+
    geom_segment2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, xend = x -n_pointmt +14, y = y, yend = y, color = VAF14), size=2)+
    geom_segment2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, xend = x -n_pointmt +15, y = y, yend = y, color = VAF15), size=2)+
    geom_segment2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, xend = x -n_pointmt +16, y = y, yend = y, color = VAF16), size=2)+
    geom_segment2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, xend = x -n_pointmt +17, y = y, yend = y, color = VAF17), size=2)+
    geom_segment2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, xend = x -n_pointmt +18, y = y, yend = y, color = VAF18), size=2)+
    geom_segment2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, xend = x -n_pointmt +19, y = y, yend = y, color = VAF19), size=2)+
    geom_segment2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, xend = x -n_pointmt +20, y = y, yend = y, color = VAF20), size=2)+
    scale_color_gradientn(colors = my_palette, breaks = col_breaks)+
    theme(axis.text = element_text(size=15), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  print(g)
}

Assign_Lineage_ID <- function(ts_dt, this_db){
  merged_tbl = tibble()
  file_list <- list.files(paste0(pointmt_path_tbl$pointmt_path[pointmt_path_tbl$deadbody == this_db],'perLineage/'), pattern = 'L.*\\.pointmts.txt$', full.names = T)
  for (n in 1:length(file_list)){
    dt <- read_tsv(file_list[n], col_types = cols(`#CHROM`='c', REF='c', ALT='c'))
    l_id <- gsub('.pointmts.txt','',unlist(strsplit(file_list[n],'//'))[2])
    dt <- dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_'), lineage_id2 = paste(this_db, l_id, sep='_')) %>%
      separate("blood_dp;blood_var;blood_vafpct", c("blood_dp","blood_var","blood_vafpct"), sep=';') %>%
      mutate(WGS_bloodVAF = as.numeric(blood_vafpct)/100) %>% select(var_id, lineage_id2, WGS_bloodVAF)
    merged_tbl <- bind_rows(merged_tbl, left_join(ts_dt, dt) %>% filter(is.na(lineage_id2) == F))
  }
  return(merged_tbl)
}




#input table (only early mutations) (deprecated)
input_tbl <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10_cancer_WGS/05_pileup_early_mutations/DB10_early_pointmts.txt.readc')
nrow(input_tbl)
input_tbl <- input_tbl %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_')) 
#get f_vaf_tbl to extract rank information
f_vaf_tbl <- readRDS('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Target_seq_analysis_afterTopUp_RDS/f_vaf_tbl_rank.rds')
input_tbl <- left_join(input_tbl, f_vaf_tbl %>% filter(deadbody == 'DB10' & time_of_lineage == 'early') %>% select(var_id, rank, n_vars_same_rank)%>% unique())
input_tbl <- input_tbl %>% separate(`10_Breastcancer_WGS_ref;10_Breastcancer_WGS_var;10_Breastcancer_WGS_ukn;10_Breastcancer_WGS_vafpct`, 
                       c("n_ref","n_var","n_ukn","vaf"), sep=';')
input_tbl <- input_tbl %>% mutate(lineage_id = gsub('DB10_','',lineage_id2))
input_tbl <- left_join(input_tbl, lng_dt %>% filter(deadbody == 'DB10') %>% select(lineage_id, n_pointmt)) #add n_pointmt info from lng_dt
nwk_path <- nwk_tbl$nwk_path[nwk_tbl$deadbody=='DB10']
Color_VAF_on_phylo_2(input_tbl, nwk_path)
                                                    

#input table (total DB10 mutations)
input_tbl <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db10_cancer_WGS/somatic_call_with_blood/05_pileup_early_mutations/DB10_merged_subs_indels.txt.readc', col_types=cols(`#CHROM`='c'))
nrow(input_tbl)
input_tbl <- input_tbl %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_')) 
input_tbl <- input_tbl %>% separate(`10_Breastcancer_WGS_ref;10_Breastcancer_WGS_var;10_Breastcancer_WGS_ukn;10_Breastcancer_WGS_vafpct`, 
                                    c("n_ref","n_var","n_ukn","vaf"), sep=';', convert=T)
input_tbl <- Assign_Lineage_ID(input_tbl, 'DB10')
m_input_tbl <- input_tbl %>% mutate(lineage_id = gsub('DB10_','',lineage_id2))
m_input_tbl <- m_input_tbl %>% select(var_id, n_var, vaf, lineage_id) %>% unique()
m_input_tbl$vaf <- as.numeric(m_input_tbl$vaf)
m_input_tbl <- left_join(m_input_tbl, lng_dt %>% filter(deadbody == 'DB10') %>%  select(lineage_id, time_of_lineage) %>% unique())
f_input_tbl <- m_input_tbl %>% filter(n_var >2 | time_of_lineage == 'early')
f_input_tbl <- f_input_tbl %>% group_by(lineage_id) %>% mutate(rank = rank(-vaf, ties.method='first'))
#f_input_tbl %>% filter(lineage_id == 'L1')
f_input_tbl <- f_input_tbl %>% mutate(n_vars_same_rank = 1)
f_input_tbl <- f_input_tbl %>% filter(rank < 20 )
f_input_tbl <- left_join(f_input_tbl, lng_dt %>% filter(deadbody == 'DB10') %>% select(lineage_id, n_pointmt)) #add n_pointmt info from lng_dt
nwk_path <- nwk_tbl$nwk_path[nwk_tbl$deadbody=='DB10']
pdf( "/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/EDFig8e.pdf", width=8, height=6)
Color_VAF_on_phylo_2(f_input_tbl, nwk_path)
dev.off()

pdf( "/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/EDFig8e_legend.pdf")
#color setting
total_break=100
start_value=0
end_value=100
change_value1=10
change_value2=30
col_breaks = seq(start_value, end_value, length.out=total_break)
a <- length(col_breaks[col_breaks < change_value1] )  #change color 
b <- length(col_breaks[col_breaks < change_value2])-a  #change color
rampCol1 <- colorRampPalette(c("gray", "yellow"))(n=a)
rampCol2 <- colorRampPalette(c("yellow","red"))(n=b)
rampCol3 <- colorRampPalette(c("red", "blue"))(n=total_break-(a+b))
my_palette<- c(rampCol1, rampCol2, rampCol3)
#print legend
if(T){
  plot.new()
  plot(c(0,2),c(start_value,end_value),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
  colfunc1 <- rev(rampCol1)
  legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
  rasterImage(legend_image1, 0, start_value, 0.3,change_value1)
  colfunc2 <- rev(rampCol2)
  legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
  rasterImage(legend_image2, 0, change_value1, 0.3,change_value2)
  colfunc3 <- rev(rampCol3)
  legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
  rasterImage(legend_image3, 0, change_value2, 0.3,end_value)
  text(x=0.5, y = c(start_value, change_value1, change_value2, end_value), labels =c(start_value, change_value1, change_value2,  end_value), cex=2)
}
dev.off()
