
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
meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_210322.txt') %>% filter(current_final_for_lineage == 'Y')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_210322.txt')
DB_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_DB_210322.txt')
rate_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Mutation_rate_simulation/df_forest_plot.txt')
m_rate_dt <- rate_dt %>% select(-ll, -ul) %>%  spread(key=type, value=mean) %>% dplyr::rename(deadbody = ID)

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


#filtered blood VAF on tree
#Fig 3d
if(T){
  Draw_blVAF_on_tree <- function(lng_dt, this_db){   
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
  Draw_blVAF_on_tree(lng_dt, 'DB6')
  Draw_blVAF_on_tree(lng_dt, 'DB9')
}

#updated code
if(T){
  Draw_blVAF_on_tree2 <- function(lng_dt, this_db, legend){   
    #color setting
    total_break=100
    start_value=0
    end_value=50
    change_value1=5
    change_value2=20
    col_breaks = seq(start_value, end_value, length.out=total_break)
    a <- length(col_breaks[col_breaks < change_value1] )  #change color 
    b <- length(col_breaks[col_breaks < change_value2])-a  #change color
    rampCol1 <- colorRampPalette(c("gray", "yellow"))(n=a)
    rampCol2 <- colorRampPalette(c("yellow","red"))(n=b)
    rampCol3 <- colorRampPalette(c("red", "blue"))(n=total_break-(a+b))
    my_palette<- c(rampCol1, rampCol2, rampCol3)
    #print legend
    if(legend =='on'){
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
    }
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
    m_lng_dt <- m_lng_dt %>% mutate(lid = taxa)
    
    g <- ggtree(tree) %<+%  m_lng_dt + theme_tree2()+
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
      #geom_text2(aes(x=x-n_pointmt, y = y+0.5, label = lid), hjust=0,color="red",size=2.5)+
      scale_color_gradientn(limits=c(start_value,end_value),colors = my_palette, breaks = col_breaks)+
      theme(axis.text = element_text(size=15), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))+
      ggtitle(this_db)
    g2 <- ggtree(tree) %<+% m_lng_dt + theme_tree2()+
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
      scale_color_gradientn(limits = c(start_value, end_value), colors = my_palette, breaks = col_breaks)+
      theme(axis.text = element_text(size=15), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
    print(g)
  }
  plot_list=list()
  for (this_db in c('DB2','DB3','DB5','DB6','DB8', 'DB10')){
    plot_list[[this_db]] <- Draw_blVAF_on_tree2(lng_dt, this_db, legend='off')
  }
  #DB10 제외함 tumor contam
}

pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/EDFig8b.pdf', width=12, height=8)
plot_grid(plotlist = plot_list[c(2,4,5,6)], nrow=1, align="h")
dev.off()

pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/EDFig8b.pdf', width=12, height=8)
grid.newpage()
pushViewport(viewport(x=0, y=0.5, width = 0.5, height = 0.5, just = c('left', 'bottom')))
print(plot_list[['DB3']], newpage = F)
popViewport(1)
pushViewport(viewport(x=0.5, y=0, width = 0.5, height = 1, just = c('left', 'bottom')))
print(plot_list[['DB6']], newpage = F)
popViewport(1)
pushViewport(viewport(x=0, y=0, width = 0.5, height = 0.5, just = c('left', 'bottom')))
print(plot_list[['DB8']], newpage = F)
popViewport(1)
dev.off()

pdf('/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/EDFig8b_legend.pdf')
total_break=100
start_value=0
end_value=50
change_value1=5
change_value2=20
col_breaks = seq(start_value, end_value, length.out=total_break)
a <- length(col_breaks[col_breaks < change_value1] )  #change color 
b <- length(col_breaks[col_breaks < change_value2])-a  #change color
rampCol1 <- colorRampPalette(c("gray", "yellow"))(n=a)
rampCol2 <- colorRampPalette(c("yellow","red"))(n=b)
rampCol3 <- colorRampPalette(c("red", "blue"))(n=total_break-(a+b))
my_palette<- c(rampCol1, rampCol2, rampCol3)
#print legend
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
dev.off()




#observed/expected blood VAF ratio
#EDFig 5c
if(T){
  Draw_obs_exp_ratio_blVAF <- function(lng_dt, pointmt_path_tbl, this_db){   
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
      res_tbl1[n,] <- list(this_lid, expected_value*100)
    }
    m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db), res_tbl1)
    
    
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
    res_mx <- res_tbl %>% as.data.frame() %>% column_to_rownames('lineage_id2')
    tmp_v <- apply(res_mx, 1, function(x) mean(x[x >0 & is.na(x) == F]))
    tmp_v2 <- apply(res_mx, 1, function(x) length(x[x >0 & is.na(x) == F]))
    tmp_v2
    m_lng_dt <- left_join(m_lng_dt, tibble(lineage_id2 = names(tmp_v), adj_mean_blVAF = tmp_v, adj_n_blVAF = tmp_v2)) %>% filter(deadbody == this_db)
    m_lng_dt <- m_lng_dt %>% mutate(obs_exp_ratio = adj_mean_blVAF/expected_VAF,
                                    adj_mtn_blVAF = start_mtn + adj_n_blVAF-1)
    
    con_tbl <- tibble(start_lid = as.character(), end_lid = as.character(), x = as.numeric(), xend = as.numeric(), y=as.numeric(), yend=as.numeric())
    for (lid in m_lng_dt$lineage_id2){
      con_tbl <- bind_rows(con_tbl, m_lng_dt %>% filter(grepl(paste0(lid,'-'), lineage_id2) & step_n == m_lng_dt$step_n[m_lng_dt$lineage_id2 == lid]+1 ) %>% select(lineage_id2, adj_mtn_blVAF, obs_exp_ratio) %>% dplyr::rename(end_lid = lineage_id2, xend = adj_mtn_blVAF, yend = obs_exp_ratio) %>%
                             mutate(start_lid = lid, x = m_lng_dt$adj_mtn_blVAF[m_lng_dt$lineage_id2 == lid], y = m_lng_dt$obs_exp_ratio[m_lng_dt$lineage_id2 == lid]))
    }
  
  
    m_lng_dt <- m_lng_dt %>% mutate(vaf_ratio =adj_mean_blVAF/expected_VAF,
                                    adjusted_mtn = start_mtn + adj_n_blVAF -1)
    m_lng_dt %>% arrange(desc(vaf_ratio)) %>% select(deadbody, lineage_id, adjusted_mtn, vaf_ratio) %>% print()
    diverge_point <- m_lng_dt %>% arrange(desc(vaf_ratio)) %>% head(n=1) %>% pull(adjusted_mtn) 
    g <- ggplot(m_lng_dt, aes(x=start_mtn + adj_n_blVAF-1, y=adj_mean_blVAF/expected_VAF))+
      geom_point(size=3, alpha=0.7)+
      coord_cartesian(xlim=c(0,35))+
      geom_segment(data = con_tbl, aes(x=x, xend=xend, y=y, yend=yend))+
      geom_hline(yintercept=1, linetype = 'dashed', color="red")+
      theme_syp+
      xlab("No. of mutations")+ylab("Observed/Expected VAF ratio")+
      ggtitle(this_db)
    print(g)
    #print(diverge_point)
  }
  plot_list <- list()
  n=0
  for (this_db in c('DB3','DB6','DB8','DB9')){
    print(this_db)
    n = n+1
    plot_list[[n]] <- Draw_obs_exp_ratio_blVAF(lng_dt, pointmt_path_tbl, this_db)
  }
  library(cowplot)
  pdf( "/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/EDFig8c.pdf", width=12, height=4)
  plot_grid(plotlist = plot_list, nrow=2, align="hv")
  dev.off()
  #DB3 14 / DB6:9/ DB8: 17/ DB9: 15 / DB10: 11
}

tmp1 <- tibble(deadbody = c('DB3','DB6','DB8','DB9','DB10'), mutation_time = c(14,9,17,15,11))
tmp1 <- left_join(tmp1, m_rate_dt)
tmp1 <- tmp1 %>% mutate(stage = (mutation_time - M0)/M1 +2)
tmp1



#blood-specific mutations, VAF histogram
#EDFig 4d
if(T){
  ch_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/clonal_hematopoiesis/CH_final_list.txt', col_types = cols(`#CHROM`='c'))
  library(mgsub)
  ch_dt <- ch_dt %>% mutate(context = mgsub(paste0(REF,ALT), c("GT","GC","GA","AT","AG","AC"), c("CA","CG","CT","TA","TC","TG")))
  #DB2 F, DB3 M, DB5 M, DB6 F, DB8 F, DB9 M, DB10 F
  male_dbs <- c("DB3","DB5","DB9")
  ch_dt <- ch_dt %>% mutate(adjVAF = ifelse(deadbody %in% male_dbs & `#CHROM` %in% c('X','Y'), VAF/2, VAF))
  ch_dt %>% group_by(deadbody) %>% count()
  plot_list =list()
  n=0
  for (this_db in c('DB6','DB8','DB9')){
    n=n+1
    f_dt <- ch_dt %>% filter(deadbody == this_db)
    g <- ggplot(f_dt, aes(x=adjVAF*100))+
      geom_histogram()+ 
      scale_x_continuous(limits=c(0,100), breaks=seq(0,100,25))+
      xlab("Variant-allele frequency(%)")+ylab("No. of mutations")+
      ggtitle(paste0(this_db, ' (n=', nrow(f_dt), ')'))+
      theme_syp
    plot_list[[n]] <- g
  }
  pdf( "/home/users/sypark/00_Project/06_LineageTracing/Rscripts/Scripts_for_revision1/EDFigs/EDFig8d.pdf", width=4, height=10)
  plot_grid(plotlist=plot_list, ncol=1, align="hv")
  dev.off()
  
  
  ch_dt$deadbody <- factor(ch_dt$deadbody, levels = c('DB2','DB5','DB6','DB8','DB9','DB10'))
  ggplot(subset(ch_dt, deadbody %in% c('DB6','DB8','DB9','DB10')))+
    geom_histogram(aes(x=VAF))+
    theme_syp+
    xlab("Variant-allele frequency")+ylab("No. of blood-specific mutations")+
    facet_wrap(~deadbody, scales = "free_y")
}

#total db VAF and number scatter plot
#Fig 4f
if(T){
  ch_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/merged_data/clonal_hematopoiesis/CH_final_list.txt', col_types = cols(`#CHROM`='c'))
  colnames(ch_dt)
  male_dbs <- c("DB3","DB5","DB9")
  ch_dt <- ch_dt %>% mutate(adjVAF = ifelse(deadbody %in% male_dbs & `#CHROM` %in% c('X','Y'), VAF/2, VAF),
                            VCF = adjVAF/0.48)
  tmp_tbl <- left_join(ch_dt %>% select(var_id, adjVAF, VCF, deadbody) %>% group_by(deadbody) %>% summarise(med_VCF = median(VCF), q1_VCF = quantile(VCF)[2], q2_VCF = quantile(VCF)[4]),
                       ch_dt %>% group_by(deadbody) %>% count())
  tmp_tbl[7,] <-  list('DB3',0,0,0,0)
  tmp_tbl <- tmp_tbl %>% filter(deadbody %in% c('DB3','DB5','DB6','DB8','DB9','DB10'))
  
  ggplot(tmp_tbl, aes(x=n, y=med_VCF, color=deadbody))+
    geom_point(size=3)+
    geom_errorbar(aes(ymin = q1_VCF, ymax = q2_VCF))+
    scale_color_manual(values = db_pal)+
    scale_y_continuous(limits=c(0,0.7))+
    ylab("Fraction of expanded clone in blood")+
    xlab("No. of blood-specific mutations")+
    theme_syp
}

