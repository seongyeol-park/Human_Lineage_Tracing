


#Define Functions
Wilcox_Test_on_Phylo_Ecto <- function(merged_m_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co){
  print(this_db)
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
  fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group = ifelse(dominant_layer != 'ectoderm', dominant_layer, paste0(tissue_type, '_epidermis'))) 
  fm_vaf_tbl <- fm_vaf_tbl %>% mutate(group2 = ifelse(tissue_class2 %in% c('lung', 'liver'), 'endoderm', group))
  f_vaf_tbl <- fm_vaf_tbl %>% filter(group2 %in% c('mesoderm','endoderm','fresh_frozen_epidermis') & is.na(rank)== F)
  f_vaf_tbl <- f_vaf_tbl %>% mutate(var_g_id = paste0(lineage_id2,'_',rank))
  res_tbl = tibble(var_g_id = as.character(), wilcox_p_Eclow =as.numeric(), wilcox_p_Echi = as.numeric(), wilcox_p_both = as.numeric(),
                   ecto_mean = as.numeric(), mesoendo_mean = as.numeric(), ecto_n = as.integer(), mesoendo_n = as.integer())
  n=0
  for(t_varid in unique(f_vaf_tbl$var_g_id)){
    n=n+1
    a <- f_vaf_tbl %>% filter(var_g_id == t_varid & group == 'fresh_frozen_epidermis') %>% pull(VAF)
    b <- f_vaf_tbl %>% filter(var_g_id == t_varid & group != 'fresh_frozen_epidermis') %>% pull(VAF)
    if (length(a) > 0 & length(b) > 0 ){
      res1 <- wilcox.test(a,b, alternative = 'less')
      pval1 <- res1$p.value
      res2 <- wilcox.test(a,b, alternative = 'greater')
      pval2 <- res2$p.value
      res3 <- wilcox.test(a,b, alternative = 'two.sided')
      pval3 <- res3$p.value
      
    } else {
      pval1 <- NA
      pval2 <- NA
      pval3 <- NA
    }
    res_tbl[n,] <- list(t_varid, pval1, pval2, pval3, mean(a), mean(b), length(a), length(b))
  }
  res_tbl <- res_tbl %>% mutate(wilcox_p = ifelse(wilcox_p_Eclow < wilcox_p_Echi, wilcox_p_both, -1*wilcox_p_both), discrete = floor(wilcox_p*100))
  res_tbl2 <- res_tbl %>% select(var_g_id, wilcox_p)
  res_tbl <- res_tbl %>% mutate(discrete = as.factor(discrete)) %>% select(var_g_id, discrete)
  
  #wilcox pvalue on phylogenetic tree
  if(T){ #order in rank (random rank in ties)
    m_res_tbl <- left_join(res_tbl,
                           merged_m_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, rank_ctie_first) %>% unique() %>% 
                             mutate(var_g_id = paste0(lineage_id2, '_', rank))
    )
    m_res_tbl <- m_res_tbl %>% mutate(rank2 = paste0('VAF', rank_ctie_first)) %>% select(-var_id, -var_g_id, -rank_ctie_first, -rank) %>% spread(key=rank2, value=discrete)
    max_VAF_len =20
    for(id in setdiff(paste0('VAF',1:max_VAF_len), colnames(m_res_tbl))){
      m_res_tbl <- m_res_tbl %>% mutate(this = NA)
      colnames(m_res_tbl) <- gsub('this',id,colnames(m_res_tbl))
    }
  }
  
  if(T){ #continuous variable
    library(colorRamps)
    total_break=50
    col_breaks = seq(0, 0.5, length.out=total_break)
    a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
    
    b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
    rampCol1 <- colorRampPalette(c("red", "orange"))(n=a)
    rampCol2 <- colorRampPalette(c("orange","gray"))(n=b)
    rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette1<- c(rampCol1, rampCol2, rampCol3)
    
    #legend..................................
    if(T){
      plot.new()
      plot(c(0,2),c(0,20),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Meso-Endoderm > Ectoderm')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,20)
      text(x=0.5, y = c(0,5,10,20), labels =c(0,0.05, 0.1, 0.2), cex=2)
    }
    
    total_break=50
    col_breaks = seq(0, 0.5, length.out=total_break)
    a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
    
    b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
    rampCol1 <- colorRampPalette(c(pal_combi[4], pal_combi[2]))(n=a)
    rampCol2 <- colorRampPalette(c(pal_combi[2],"gray"))(n=b)
    rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette2<- c(rampCol1, rampCol2, rampCol3)
    
    
    #legend..................................
    if(T){
      plot.new()
      plot(c(0,2),c(0,20),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Meso-Endoderm < Ectoderm')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,20)
      text(x=0.5, y = c(0,5,10,50), labels =c(0,0.05, 0.1, 0.2), cex=2)
    }
    
    my_palette <- c(rev(my_palette2), my_palette1)
    names(my_palette) <- seq(-49,48,1)
  }
  
  #assign point changing to <0.05 and >= 0.05
  if(T){
    m_res_tbl2 <- left_join(res_tbl2,
                            merged_m_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, n_vars_same_rank, start_mtn) %>% unique() %>% 
                              mutate(var_g_id = paste0(lineage_id2, '_', rank)))
    
    el_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early_to_late') %>% pull(lineage_id2)
    changing_lines = tibble()
    for(t_id in el_ids){
      find_v = 'N'
      tmp_v <- unlist(strsplit(t_id, '-'))
      t_id_list <- c()
      for (i in 1:length(tmp_v)){
        t_id_list <- c(t_id_list, paste(tmp_v[1:i], collapse="-"))
      }
      for(n in 1:(length(t_id_list)-1)){
        rank_list <- m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n]) %>% pull(rank) %>% unique() %>% sort()
        if(length(rank_list) == 0){
          next
        }
        for (m in 1:length(rank_list)){
          t_pv <- m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n] & rank == rank_list[m]) %>% pull(wilcox_p) %>% unique() %>% abs()
          if(is.na(t_pv)){
            next
          }
          if (t_pv < pval_co){
            changing_lines <- bind_rows(changing_lines, m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n] & rank == rank_list[m]) %>% select(-var_id) %>% unique())
            find_v = 'Y'
            break
          }
        }
        if (find_v == 'Y'){
          break
        }
      }
    }
    changing_lines <- changing_lines %>% unique()
    changing_lines <- changing_lines %>% mutate(point0.05 = (start_mtn + rank - 1 + start_mtn + rank - 1 + n_vars_same_rank)/2)
  }
  my_func <- function(num){
    if (is.na(num)){
      return (NA)
    } else if(num < -49){
      return (-49)
    } else if ( num > 48){
      return (48) 
    } else {
      return (num)
    }
  }
  m_res_tbl <- m_res_tbl %>% mutate_at(vars(starts_with('VAF')), as.character) %>% mutate_at(vars(starts_with('VAF')), as.numeric) %>%
    rowwise() %>% mutate_at(vars(starts_with('VAF')),my_func) %>% mutate_at(vars(starts_with('VAF')), as.factor)
  m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db) ,m_res_tbl)
  m_lng_dt <- left_join(m_lng_dt, changing_lines %>% select(lineage_id2, point0.05))
  m_lng_dt <- m_lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  g <- ggtree(tree) %<+% m_lng_dt + theme_tree2()+
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
    ggtitle(paste0(this_db, ', Ectoderm-MesoEndoderm divergence'))+
    scale_color_manual(values = my_palette)+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  print(g)
  return(changing_lines)
}
Wilcox_Test_on_Phylo_LR <- function(merged_m_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co, s_range){
  print(this_db)
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
  if(s_range == 'EandT'){
    fm_vaf_tbl <- fm_vaf_tbl %>% mutate(Source_side2 = ifelse(anatomy_class1 %in% c('LE','UE','trunk'), Source_side, 'unknown'))
  } else if (s_range == 'Eonly'){
    fm_vaf_tbl <- fm_vaf_tbl %>% mutate(Source_side2 = ifelse(anatomy_class1 %in% c('LE','UE'), Source_side, 'unknown'))
  } else {
    print('wrong s_range')
    stop()
  }
  print(fm_vaf_tbl %>% group_by(Source_side2, var_id) %>% count() %>% group_by(Source_side2) %>% summarise(mean_n= mean(n)))
  f_vaf_tbl <- fm_vaf_tbl %>% filter(is.na(rank) == F)
  f_vaf_tbl <- f_vaf_tbl %>% mutate(var_g_id = paste0(lineage_id2,'_',rank))
  res_tbl = tibble(var_g_id = as.character(), wilcox_p_rthi =as.numeric(), wilcox_p_lthi = as.numeric(), wilcox_p_both = as.numeric())
  n=0
  for(t_varid in unique(f_vaf_tbl$var_g_id)){
    n=n+1
    a <- f_vaf_tbl %>% filter(var_g_id == t_varid & Source_side2 == 'lt') %>% pull(VAF)
    b <- f_vaf_tbl %>% filter(var_g_id == t_varid & Source_side2 == 'rt') %>% pull(VAF)
    if (length(a) > 0 & length(b) > 0 ){
      res1 <- wilcox.test(a,b, alternative = 'less')
      pval1 <- res1$p.value
      res2 <- wilcox.test(a,b, alternative = 'greater')
      pval2 <- res2$p.value
      res3 <- wilcox.test(a,b, alternative = 'two.sided')
      pval3 <- res3$p.value
      
    } else {
      pval1 <- NA
      pval2 <- NA
      pval3 <- NA
    }
    res_tbl[n,] <- list(t_varid, pval1, pval2, pval3)
  }
  res_tbl <- res_tbl %>% mutate(wilcox_p = ifelse(wilcox_p_rthi < wilcox_p_lthi, wilcox_p_both, -1*wilcox_p_both), discrete = floor(wilcox_p*100))
  res_tbl2 <- res_tbl %>% select(var_g_id, wilcox_p)
  res_tbl <- res_tbl %>% mutate(discrete = as.factor(discrete)) %>% select(var_g_id, discrete)
  
  #wilcox pvalue on phylogenetic tree
  if(T){ #order in rank (random rank in ties)
    m_res_tbl <- left_join(res_tbl,
                           merged_m_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, rank_ctie_first) %>% unique() %>% 
                             mutate(var_g_id = paste0(lineage_id2, '_', rank))
    )
    m_res_tbl <- m_res_tbl %>% mutate(rank2 = paste0('VAF', rank_ctie_first)) %>% select(-var_id, -var_g_id, -rank_ctie_first, -rank) %>% spread(key=rank2, value=discrete)
    max_VAF_len =20
    for(id in setdiff(paste0('VAF',1:max_VAF_len), colnames(m_res_tbl))){
      m_res_tbl <- m_res_tbl %>% mutate(this = NA)
      colnames(m_res_tbl) <- gsub('this',id,colnames(m_res_tbl))
    }
  }
  
  
  if(T){ #continuous variable
    library(colorRamps)
    total_break=50
    col_breaks = seq(0, 0.5, length.out=total_break)
    a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
    
    b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
    rampCol1 <- colorRampPalette(c("red", "orange"))(n=a)
    rampCol2 <- colorRampPalette(c("orange","gray"))(n=b)
    rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette1<- c(rampCol1, rampCol2, rampCol3)
    
    #legend..................................
    if(T){
      plot.new()
      plot(c(0,2),c(0,50),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Rt. VAF > Lt. VAF')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,20)
      text(x=0.5, y = c(0,5,10,20), labels =c(0,0.05, 0.1, 0.2), cex=2)
    }
    
    total_break=50
    col_breaks = seq(0, 0.5, length.out=total_break)
    a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
    
    b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
    rampCol1 <- colorRampPalette(c(pal_combi[4], pal_combi[2]))(n=a)
    rampCol2 <- colorRampPalette(c(pal_combi[2],"gray"))(n=b)
    rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette2<- c(rampCol1, rampCol2, rampCol3)
    
    
    #legend..................................
    if(T){
      plot.new()
      plot(c(0,2),c(0,50),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Rt. VAF < Lt. VAF')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,20)
      text(x=0.5, y = c(0,5,10,20), labels =c(0,0.05, 0.1, 0.2), cex=2)
    }
    
    my_palette <- c(rev(my_palette2), my_palette1)
    names(my_palette) <- seq(-49,48,1)
  }
  
  #assign point changing to <0.05 and >= 0.05
  if(T){
    m_res_tbl2 <- left_join(res_tbl2,
                            merged_m_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, n_vars_same_rank, start_mtn) %>% unique() %>% 
                              mutate(var_g_id = paste0(lineage_id2, '_', rank)))
    
    el_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early_to_late') %>% pull(lineage_id2)
    changing_lines = tibble()
    for(t_id in el_ids){
      find_v = 'N'
      tmp_v <- unlist(strsplit(t_id, '-'))
      t_id_list <- c()
      for (i in 1:length(tmp_v)){
        t_id_list <- c(t_id_list, paste(tmp_v[1:i], collapse="-"))
      }
      for(n in 1:(length(t_id_list)-1)){
        rank_list <- m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n]) %>% pull(rank) %>% unique() %>% sort()
        if(length(rank_list) == 0){
          next
        }
        for (m in 1:length(rank_list)){
          t_pv <- m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n] & rank == rank_list[m]) %>% pull(wilcox_p) %>% unique() %>% abs()
          if(is.na(t_pv)){
            next
          }
          if (t_pv < pval_co){
            changing_lines <- bind_rows(changing_lines, m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n] & rank == rank_list[m]) %>% select(-var_id) %>% unique())
            find_v = 'Y'
            break
          }
        }
        if (find_v == 'Y'){
          break
        }
      }
    }
    changing_lines <- changing_lines %>% unique()
    changing_lines <- changing_lines %>% mutate(point0.05 = (start_mtn + rank - 1 + start_mtn + rank - 1 + n_vars_same_rank)/2)
  }
  my_func <- function(num){
    if (is.na(num)){
      return (NA)
    } else if(num < -49){
      return (-49)
    } else if ( num > 48){
      return (48) 
    } else {
      return (num)
    }
  }
  m_res_tbl <- m_res_tbl %>% mutate_at(vars(starts_with('VAF')), as.character) %>% mutate_at(vars(starts_with('VAF')), as.numeric) %>%
    rowwise() %>% mutate_at(vars(starts_with('VAF')),my_func) %>% mutate_at(vars(starts_with('VAF')), as.factor)
  m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db) ,m_res_tbl)
  m_lng_dt <- left_join(m_lng_dt, changing_lines %>% select(lineage_id2, point0.05))
  m_lng_dt <- m_lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  g <- ggtree(tree) %<+% m_lng_dt + theme_tree2()+
    coord_cartesian(xlim = c(-1,30))+
    geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y,color = VAF1) ,size=2)+
    #geom_text2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, y = y, label = round(VAF1,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2), size=2)+
    #geom_text2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, y = y, label = round(VAF2,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3), size=2)+
    #geom_text2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, y = y, label = round(VAF3,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4), size=2)+
    #geom_text2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, y = y, label = round(VAF4,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, xend = x -n_pointmt +5, y = y, yend = y, color = VAF5), size=2)+
    #geom_text2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, y = y, label = VAF5), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, xend = x -n_pointmt +6, y = y, yend = y, color = VAF6), size=2)+
    #geom_text2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, y = y, label = VAF6), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, xend = x -n_pointmt +7, y = y, yend = y, color = VAF7), size=2)+
    #geom_text2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, y = y, label = VAF7), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, xend = x -n_pointmt +8, y = y, yend = y, color = VAF8), size=2)+
    #geom_text2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, y = y, label = VAF8), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, xend = x -n_pointmt +9, y = y, yend = y, color = VAF9), size=2)+
    #geom_text2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, y = y, label = VAF9), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, xend = x -n_pointmt +10, y = y, yend = y, color = VAF10), size=2)+
    #geom_text2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, y = y, label = VAF10), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, xend = x -n_pointmt +11, y = y, yend = y, color = VAF11), size=2)+
    #geom_text2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, y = y, label = VAF11), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, xend = x -n_pointmt +12, y = y, yend = y, color = VAF12), size=2)+
    #geom_text2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, y = y, label = VAF12), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, xend = x -n_pointmt +13, y = y, yend = y, color = VAF13), size=2)+
    #geom_text2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, y = y, label = VAF13), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, xend = x -n_pointmt +14, y = y, yend = y, color = VAF14), size=2)+
    #geom_text2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, y = y, label = VAF14), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, xend = x -n_pointmt +15, y = y, yend = y, color = VAF15), size=2)+
    #geom_text2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, y = y, label = VAF15), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, xend = x -n_pointmt +16, y = y, yend = y, color = VAF16), size=2)+
    #geom_text2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, y = y, label = VAF16), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, xend = x -n_pointmt +17, y = y, yend = y, color = VAF17), size=2)+
    #geom_text2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, y = y, label = VAF17), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, xend = x -n_pointmt +18, y = y, yend = y, color = VAF18), size=2)+
    #geom_text2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, y = y, label = VAF18), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, xend = x -n_pointmt +19, y = y, yend = y, color = VAF19), size=2)+
    #geom_text2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, y = y, label = VAF19), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, xend = x -n_pointmt +20, y = y, yend = y, color = VAF20), size=2)+
    #geom_text2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, y = y, label = VAF20), size=2.5)+
    ggtitle(paste0(this_db, ', Left-right axis, ',s_range))+
    scale_color_manual(values = my_palette)+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  print(g)
  return(changing_lines)
}
Wilcox_Test_on_Phylo_CC <- function(merged_m_vaf_tbl, this_db, lng_dt, nwk_tbl, pval_co, s_range){
  fm_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db & DP > 100 & time_of_lineage == 'early' & tissue_class != 'cancer_tissue')
  if(s_range == 'HTE'){
    fm_vaf_tbl <- fm_vaf_tbl %>% mutate(Source_side2 = ifelse(anatomy_class1 %in% c('UE','HN') | anatomy_class2 == 'chest', 'cranio',
                                                              ifelse(anatomy_class1 == 'LE' | anatomy_class2 %in% c('abdomen','pubis'), 'caudal','unknown')))
  } else if (s_range == 'HandLE'){
    fm_vaf_tbl <- fm_vaf_tbl %>% mutate(Source_side2 = ifelse(anatomy_class1 %in% c('UE','HN') , 'cranio', 
                                                              ifelse(anatomy_class1 == 'LE', 'caudal','unknown')))
  } else {
    print('wrong s_range')
    stop()
  }
  
  print(this_db)
  print(fm_vaf_tbl %>% group_by(Source_side2, var_id) %>% count() %>% group_by(Source_side2) %>% summarise(mean_n= mean(n)))
  f_vaf_tbl <- fm_vaf_tbl %>% filter(is.na(rank) == F)
  f_vaf_tbl <- f_vaf_tbl %>% mutate(var_g_id = paste0(lineage_id2,'_',rank))
  res_tbl = tibble(var_g_id = as.character(), wilcox_p_crhi =as.numeric(), wilcox_p_cahi = as.numeric(), wilcox_p_both = as.numeric())
  n=0
  for(t_varid in unique(f_vaf_tbl$var_g_id)){
    n=n+1
    a <- f_vaf_tbl %>% filter(var_g_id == t_varid & Source_side2 == 'caudal') %>% pull(VAF)
    b <- f_vaf_tbl %>% filter(var_g_id == t_varid & Source_side2 == 'cranio') %>% pull(VAF)
    if (length(a) > 0 & length(b) > 0 ){
      res1 <- wilcox.test(a,b, alternative = 'less')
      pval1 <- res1$p.value
      res2 <- wilcox.test(a,b, alternative = 'greater')
      pval2 <- res2$p.value
      res3 <- wilcox.test(a,b, alternative = 'two.sided')
      pval3 <- res3$p.value
      
    } else {
      pval1 <- NA
      pval2 <- NA
      pval3 <- NA
    }
    res_tbl[n,] <- list(t_varid, pval1, pval2, pval3)
  }
  res_tbl <- res_tbl %>% mutate(wilcox_p = ifelse(wilcox_p_crhi < wilcox_p_cahi, wilcox_p_both, -1*wilcox_p_both), discrete = floor(wilcox_p*100))
  res_tbl2 <- res_tbl %>% select(var_g_id, wilcox_p)
  res_tbl <- res_tbl %>% mutate(discrete = as.factor(discrete)) %>% select(var_g_id, discrete)
  
  
  #wilcox pvalue on phylogenetic tree
  if(T){ #order in rank (random rank in ties)
    m_res_tbl <- left_join(res_tbl,
                           merged_m_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, rank_ctie_first) %>% unique() %>% 
                             mutate(var_g_id = paste0(lineage_id2, '_', rank))
    )
    m_res_tbl <- m_res_tbl %>% mutate(rank2 = paste0('VAF', rank_ctie_first)) %>% select(-var_id, -var_g_id, -rank_ctie_first, -rank) %>% spread(key=rank2, value=discrete)
    max_VAF_len =20
    for(id in setdiff(paste0('VAF',1:max_VAF_len), colnames(m_res_tbl))){
      m_res_tbl <- m_res_tbl %>% mutate(this = NA)
      colnames(m_res_tbl) <- gsub('this',id,colnames(m_res_tbl))
    }
  }
  
  if(T){ #continuous variable
    library(colorRamps)
    total_break=50
    col_breaks = seq(0, 0.5, length.out=total_break)
    a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
    
    b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
    rampCol1 <- colorRampPalette(c("red", "orange"))(n=a)
    rampCol2 <- colorRampPalette(c("orange","gray"))(n=b)
    rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette1<- c(rampCol1, rampCol2, rampCol3)
    
    #legend..................................
    if(T){
      plot.new()
      plot(c(0,2),c(0,50),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Rt. VAF > Lt. VAF')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,20)
      text(x=0.5, y = c(0,5,10,20), labels =c(0,0.05, 0.1, 0.2), cex=2)
    }
    
    total_break=50
    col_breaks = seq(0, 0.5, length.out=total_break)
    a <- length(col_breaks[col_breaks < 0.05] )  #change color at 1
    
    b <- length(col_breaks[col_breaks < 0.1])-a  #change color at 5
    rampCol1 <- colorRampPalette(c(pal_combi[4], pal_combi[2]))(n=a)
    rampCol2 <- colorRampPalette(c(pal_combi[2],"gray"))(n=b)
    rampCol3 <- colorRampPalette(c("gray", "gray"))(n=total_break-(a+b)-1)
    my_palette2<- c(rampCol1, rampCol2, rampCol3)
    
    
    #legend..................................
    if(T){
      plot.new()
      plot(c(0,2),c(0,50),type = 'n', axes = F,xlab = '', ylab = '', main = 'Pvalues from Rt. VAF < Lt. VAF')
      colfunc1 <- rev(rampCol1)
      legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
      rasterImage(legend_image1, 0, 0, 0.3,5)
      colfunc2 <- rev(rampCol2)
      legend_image2 <- as.raster(matrix(colfunc2, ncol=1))
      rasterImage(legend_image2, 0, 5, 0.3,10)
      colfunc3 <- rev(rampCol3)
      legend_image3 <- as.raster(matrix(colfunc3, ncol=1))
      rasterImage(legend_image3, 0, 10, 0.3,20)
      text(x=0.5, y = c(0,5,10,20), labels =c(0,0.05, 0.1, 0.2), cex=2)
    }
    
    my_palette <- c(rev(my_palette2), my_palette1)
    names(my_palette) <- seq(-49,48,1)
  } 
  
  #assign point changing to <0.05 and >= 0.05
  if(T){
    m_res_tbl2 <- left_join(res_tbl2,
                            merged_m_vaf_tbl %>% filter(deadbody == this_db & time_of_lineage == 'early' & is.na(rank) == F ) %>% select(lineage_id2, var_id, rank, n_vars_same_rank, start_mtn) %>% unique() %>% 
                              mutate(var_g_id = paste0(lineage_id2, '_', rank)))
    
    el_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early_to_late') %>% pull(lineage_id2)
    changing_lines = tibble()
    for(t_id in el_ids){
      find_v = 'N'
      tmp_v <- unlist(strsplit(t_id, '-'))
      t_id_list <- c()
      for (i in 1:length(tmp_v)){
        t_id_list <- c(t_id_list, paste(tmp_v[1:i], collapse="-"))
      }
      for(n in 1:(length(t_id_list)-1)){
        rank_list <- m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n]) %>% pull(rank) %>% unique() %>% sort()
        if(length(rank_list) == 0){
          next
        }
        for (m in 1:length(rank_list)){
          t_pv <- m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n] & rank == rank_list[m]) %>% pull(wilcox_p) %>% unique() %>% abs()
          if(is.na(t_pv)){
            next
          }
          if (t_pv < pval_co){
            changing_lines <- bind_rows(changing_lines, m_res_tbl2 %>% filter(lineage_id2 == t_id_list[n] & rank == rank_list[m]) %>% select(-var_id) %>% unique())
            find_v = 'Y'
            break
          }
        }
        if (find_v == 'Y'){
          break
        }
      }
    }
    changing_lines <- changing_lines %>% unique()
    changing_lines <- changing_lines %>% mutate(point0.05 = (start_mtn + rank - 1 + start_mtn + rank - 1 + n_vars_same_rank)/2)
  }
  my_func <- function(num){
    if (is.na(num)){
      return (NA)
    } else if(num < -49){
      return (-49)
    } else if ( num > 48){
      return (48) 
    } else {
      return (num)
    }
  }
  m_res_tbl <- m_res_tbl %>% mutate_at(vars(starts_with('VAF')), as.character) %>% mutate_at(vars(starts_with('VAF')), as.numeric) %>%
    rowwise() %>% mutate_at(vars(starts_with('VAF')),my_func) %>% mutate_at(vars(starts_with('VAF')), as.factor)
  m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db) ,m_res_tbl)
  m_lng_dt <- left_join(m_lng_dt, changing_lines %>% select(lineage_id2, point0.05))
  m_lng_dt <- m_lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  g <- ggtree(tree) %<+% m_lng_dt + theme_tree2()+
    coord_cartesian(xlim = c(-1,30))+
    geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y,color = VAF1) ,size=2)+
    #geom_text2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, y = y, label = round(VAF1,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2), size=2)+
    #geom_text2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, y = y, label = round(VAF2,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3), size=2)+
    #geom_text2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, y = y, label = round(VAF3,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4), size=2)+
    #geom_text2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, y = y, label = round(VAF4,3)), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, xend = x -n_pointmt +5, y = y, yend = y, color = VAF5), size=2)+
    #geom_text2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt+4, y = y, label = VAF5), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, xend = x -n_pointmt +6, y = y, yend = y, color = VAF6), size=2)+
    #geom_text2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt+5, y = y, label = VAF6), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, xend = x -n_pointmt +7, y = y, yend = y, color = VAF7), size=2)+
    #geom_text2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt+6, y = y, label = VAF7), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, xend = x -n_pointmt +8, y = y, yend = y, color = VAF8), size=2)+
    #geom_text2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt+7, y = y, label = VAF8), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, xend = x -n_pointmt +9, y = y, yend = y, color = VAF9), size=2)+
    #geom_text2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt+8, y = y, label = VAF9), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, xend = x -n_pointmt +10, y = y, yend = y, color = VAF10), size=2)+
    #geom_text2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt+9, y = y, label = VAF10), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, xend = x -n_pointmt +11, y = y, yend = y, color = VAF11), size=2)+
    #geom_text2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt+10, y = y, label = VAF11), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, xend = x -n_pointmt +12, y = y, yend = y, color = VAF12), size=2)+
    #geom_text2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt+11, y = y, label = VAF12), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, xend = x -n_pointmt +13, y = y, yend = y, color = VAF13), size=2)+
    #geom_text2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt+12, y = y, label = VAF13), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, xend = x -n_pointmt +14, y = y, yend = y, color = VAF14), size=2)+
    #geom_text2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt+13, y = y, label = VAF14), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, xend = x -n_pointmt +15, y = y, yend = y, color = VAF15), size=2)+
    #geom_text2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt+14, y = y, label = VAF15), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, xend = x -n_pointmt +16, y = y, yend = y, color = VAF16), size=2)+
    #geom_text2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt+15, y = y, label = VAF16), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, xend = x -n_pointmt +17, y = y, yend = y, color = VAF17), size=2)+
    #geom_text2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt+16, y = y, label = VAF17), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, xend = x -n_pointmt +18, y = y, yend = y, color = VAF18), size=2)+
    #geom_text2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt+17, y = y, label = VAF18), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, xend = x -n_pointmt +19, y = y, yend = y, color = VAF19), size=2)+
    #geom_text2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt+18, y = y, label = VAF19), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, xend = x -n_pointmt +20, y = y, yend = y, color = VAF20), size=2)+
    #geom_text2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt+19, y = y, label = VAF20), size=2.5)+
    ggtitle(paste0(this_db, ', Cranio-caudal axis, ', s_range))+
    scale_color_manual(values = my_palette)+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  print(g)
  return(changing_lines)
}


