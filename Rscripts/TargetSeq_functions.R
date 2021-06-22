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
  merged_tbl <- bind_rows(merged_tbl, ts_dt %>% filter(!var_id %in% merged_tbl$var_id) %>% mutate(lineage_id2 = paste0(this_db,'_GSNP')))
  return(merged_tbl)
}

Make_VAF_tbl <- function(merged_tbl, this_db){
  vaf_tbl1 <- merged_tbl  #copy table
  #extract dp and var information
  for (c_name in col_list){  
    print(c_name)
    vaf_tbl1[[paste0(c_name, '_dp')]] <- apply(vaf_tbl1[c_name], 1, function(x) as.numeric(gsub('dp=','',unlist(strsplit(x,';'))[1])))
    vaf_tbl1[[paste0(c_name, '_var')]] <- apply(vaf_tbl1[c_name], 1, function(x) as.numeric(gsub('var=','',unlist(strsplit(x,';'))[3])))
    vaf_tbl1 <- vaf_tbl1 %>% select(-c_name)
  }
  
  vaf_tbl1 <- vaf_tbl1 %>% select(var_id, lineage_id2, WGS_bloodVAF, ends_with('_dp'), ends_with('_var'))
  #caclulate binom.confint
  library(binom)
  for (c_name in col_list){ 
    print(c_name)
    vaf_tbl1[[paste0(c_name, '_var_low')]] <- round(vaf_tbl1[[paste0(c_name, '_dp')]]*binom.confint(vaf_tbl1[[paste0(c_name, '_var')]],vaf_tbl1[[paste0(c_name, '_dp')]], method="exact")$lower,4)
    vaf_tbl1[[paste0(c_name, '_var_up')]] <- round(vaf_tbl1[[paste0(c_name, '_dp')]]*binom.confint(vaf_tbl1[[paste0(c_name, '_var')]],vaf_tbl1[[paste0(c_name, '_dp')]], method="exact")$upper,4)
    vaf_tbl1[[paste0(c_name, '_vaf')]] <- round(vaf_tbl1[[paste0(c_name, '_var')]]/vaf_tbl1[[paste0(c_name, '_dp')]],4)
    vaf_tbl1[[paste0(c_name, '_vaf_low')]] <- round(vaf_tbl1[[paste0(c_name, '_var_low')]]/vaf_tbl1[[paste0(c_name, '_dp')]],4)
    vaf_tbl1[[paste0(c_name, '_vaf_up')]] <- round(vaf_tbl1[[paste0(c_name, '_var_up')]]/vaf_tbl1[[paste0(c_name, '_dp')]],4)
  }
  #join with lng_dt
  vaf_tbl1 <- left_join(vaf_tbl1, lng_dt %>% filter(deadbody == this_db) %>% select(lineage_id2, time_of_lineage, early_desc_n,step_n,time_of_lineage, start_mtn, end_mtn))
  tmp_tbl <- vaf_tbl1 %>% select(ends_with('_var'))
  vaf_tbl1$nPositive <- apply(tmp_tbl, 1, function(x) sum(x>0))
  vaf_tbl1$meanVarRC <- round(apply(tmp_tbl, 1, function(x) mean(x[x>0])),0)
  tmp_tbl2 <- vaf_tbl1 %>% select(ends_with('_dp'))
  vaf_tbl1$hidp100 <- apply(tmp_tbl2, 1, function(x) sum(x>=100))
  return(vaf_tbl1)
}

Add_Expected_VAF <- function(vaf_tbl1, this_db){
  #consistency of expected VAF and actual VAF of tissues
  id_list <- vaf_tbl1 %>% filter(is.na(time_of_lineage)==F) %>% .$lineage_id2 %>% unique()
  res_tbl1 = tibble(lineage_id2 = as.character(), expected_VAF = as.numeric())
  n=0
  for (this_lid in id_list){
    n= n+1
    lng_dt %>% filter(time_of_lineage =='early_to_late' & grepl(this_lid, lineage_id2)) %>% nrow() -> a
    lng_dt %>% filter(time_of_lineage =='early_to_late' & deadbody == this_db ) %>% nrow() -> b
    lower_margin = 0.48*binom.confint(a,b, method="exact")$lower
    upper_margin = 0.48*binom.confint(a,b, method="exact")$upper
    expected_value = 0.48*a/b
    res_tbl1[n,] <- list(this_lid, expected_value)
  }
  vaf_tbl1 <- left_join(vaf_tbl1, res_tbl1)
  return(vaf_tbl1)
}

Draw_L1plusL2_VAF <- function(merged_m_vaf_tbl){
  GSNP_mVAF <- merged_m_vaf_tbl %>% filter(grepl('GSNP', lineage_id2)) %>% pull(VAF) %>% mean(., na.rm=T)
  my_cols = pal_npg('nrc')(10)[c(1,4)]
  names(my_cols) <- c('L1','L2')
  p_list <- list()
  n=0
  for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
    n=n+1
    m_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db)
    fm_vaf_tbl <- m_vaf_tbl %>% filter(var_id %in% first_variants$var_id[first_variants$deadbody == this_db]) %>% filter(DP > 100)
    tmp_tbl <- fm_vaf_tbl %>% group_by(sample_id, lineage_id2) %>% summarise(meanVAF = mean(VAF)) 
    m_vaf_tbl %>% filter(lineage_id2 == paste0(this_db,'_GSNP')) %>% filter(DP > 100) %>% select(sample_id, lineage_id2, VAF)
    tmp_tbl <- left_join(tmp_tbl,pmeta_dt %>% select(sample_id, dominant_layer, anatomy_class3)) %>% separate(lineage_id2, c('deadbody','lineage_id'), sep="_")
    samples <- tmp_tbl %>% group_by(sample_id) %>% count() %>% filter(n==2) %>% pull(sample_id)
    tmp_tbl <- tmp_tbl %>% filter(sample_id %in% samples)
    #x_order <- tmp_tbl %>% spread(key=lineage_id, value=meanVAF) %>% mutate(L12 = L1+L2) %>% arrange(desc(L12)) %>% pull(sample_id)
    p_list[[n]] <-ggplot(tmp_tbl, aes(x=sample_id, y=meanVAF, fill=lineage_id))+
      geom_bar(stat="identity")+
      geom_hline(yintercept=GSNP_mVAF)+
      #scale_x_discrete(limits=x_order)+
      scale_y_continuous(limits=c(0,0.7), breaks=seq(0,0.6,0.2), expand = c(0,0))+
      scale_fill_manual(values = my_cols)+
      theme_syp+
      theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks.x = element_blank(), legend.position = 'none')+
      ggtitle(this_db)
  }
  library(cowplot)
  plot_grid(plotlist=p_list, nrow=1, align="h")
}


Draw_L1L2_VAF_relative <- function(merged_m_vaf_tbl){
  GSNP_mVAF <- merged_m_vaf_tbl %>% filter(grepl('GSNP', lineage_id2)) %>% pull(VAF) %>% mean(., na.rm=T)
  my_cols = pal_npg('nrc')(10)[c(1,4)]
  names(my_cols) <- c('L1ratio','L2ratio')
  p_list <- list()
  n=0
  for (this_db in c('DB3','DB6','DB8','DB9','DB10')){
    n=n+1
    m_vaf_tbl <- merged_m_vaf_tbl %>% filter(deadbody == this_db)
    fm_vaf_tbl <- m_vaf_tbl %>% filter(var_id %in% first_variants$var_id[first_variants$deadbody == this_db]) %>% filter(DP > 100)
    tmp_tbl <- fm_vaf_tbl %>% group_by(sample_id, lineage_id2) %>% summarise(meanVAF = mean(VAF)) 
    m_vaf_tbl %>% filter(lineage_id2 == paste0(this_db,'_GSNP')) %>% filter(DP > 100) %>% select(sample_id, lineage_id2, VAF)
    tmp_tbl <- left_join(tmp_tbl,pmeta_dt %>% select(sample_id, dominant_layer, anatomy_class3)) %>% separate(lineage_id2, c('deadbody','lineage_id'), sep="_") %>% 
      spread(key=lineage_id, value=meanVAF) %>% mutate(L1ratio=L1/(L1+L2), L2ratio=L2/(L1+L2)) %>% select(-L1, -L2) %>% filter(is.na(L1ratio)==F & is.na(L2ratio) == F) %>%
      gather(-c(sample_id,deadbody,dominant_layer, anatomy_class3), key=lineage_id, value=ratio)
    x_order <- tmp_tbl %>% filter(lineage_id == 'L1ratio') %>% arrange(desc(ratio)) %>% pull(sample_id)
    exp_v1 <- m_vaf_tbl %>% filter(lineage_id2 == paste0(this_db, '_L2')) %>% pull(expected_VAF) %>% unique()
    exp_v2 <- m_vaf_tbl %>% filter(lineage_id2 == paste0(this_db, '_L1')) %>% pull(expected_VAF) %>% unique()
    p_list[[n]] <-ggplot(tmp_tbl, aes(x=sample_id, y=ratio, fill=lineage_id))+
      geom_bar(stat="identity")+
      geom_hline(yintercept=exp_v1/(exp_v1+exp_v2), linetype="dashed")+
      scale_x_discrete(limits=x_order)+
      scale_fill_manual(values = my_cols)+
      theme_syp+
      theme(axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks.x = element_blank(), legend.position = 'none')+
      ggtitle(this_db)
  }
  library(cowplot)
  plot_grid(plotlist=p_list, nrow=2, align="hv")
}



DB_Sample_Plot <- function(fm_vaf_tbl, this_id, this_db){ 
  library(colorRamps)
  total_break=100
  col_breaks = seq(0, 100, length.out=total_break)
  a <- length(col_breaks[col_breaks < 1] )  #change color at 1
  b <- length(col_breaks[col_breaks < 5])-a  #change color at 5
  rampCol1 <- colorRampPalette(c("gray", "yellow"))(n=a)
  rampCol2 <- colorRampPalette(c("yellow","red"))(n=b)
  rampCol3 <- colorRampPalette(c("red", "blue"))(n=total_break-(a+b))
  my_palette<- c(rampCol1, rampCol2, rampCol3)
  
  print(this_id)
  annot <- paste0(unique(fm_vaf_tbl$Source_side[fm_vaf_tbl$sample_id == this_id]),'_',unique(fm_vaf_tbl$anatomy_tissue_class[fm_vaf_tbl$sample_id == this_id]))
  ffm_vaf_tbl <- fm_vaf_tbl %>% filter(sample_id == this_id)
  axis_labels = c(c(0,1,2,3,5),seq(10,70,10))
  g1 <- ggplot(ffm_vaf_tbl, aes(x=log10(expected_VAF*100+1), y=log10(VAF*100+0.1)))+
    geom_point(aes(size = DP), alpha=0.2)+
    geom_abline(slope=1, intercept=0)+
    scale_x_continuous(labels = axis_labels, breaks = log10(axis_labels + 0.1), limits = log10(c(1,70)+0.1))+
    scale_y_continuous(labels = axis_labels, breaks = log10(axis_labels + 0.1))+
    xlab("Expected VAF")+ylab("Actual VAF")+
    theme_syp
  max_VAF_len=20
  res_tbl <- as.tibble(data.frame(matrix(vector(), 0, max_VAF_len + 1,
                                         dimnames=list(c(), c('lineage_id2',paste0('VAF',1:max_VAF_len)))),
                                  stringsAsFactors=F))
  n=0
  id_list <- unique(ffm_vaf_tbl$lineage_id2)[grepl("germline",unique(ffm_vaf_tbl$lineage_id2))==F]
  id_list <-id_list[grepl("GSNP",unique(ffm_vaf_tbl$lineage_id2))==F]
  for (lid in id_list){
    n=n+1
    VAF_v <- ffm_vaf_tbl %>% filter(lineage_id2 == lid) %>% arrange(desc(VAF)) %>% .$VAF
    VAF_v <- round(VAF_v*100,2)
    NA_v <- rep(NA, max_VAF_len - length(VAF_v))
    res_tbl[n,] <- c(lid,VAF_v, NA_v)
  }
  res_tbl <- res_tbl %>% mutate_at(vars(starts_with('VAF')), list(as.numeric))
  early_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early') %>% .$lineage_id2
  res_tbl <- res_tbl %>% filter(lineage_id2 %in% early_ids)
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  m_lng_dt <- lng_dt %>% filter(deadbody == this_db) %>% rename(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt$Source_class2[m_lng_dt$time_of_lineage == 'early'] <- NA
  m_lng_dt <- left_join(m_lng_dt, res_tbl)
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
    ggtitle(paste0(this_id,' (',annot,')'))+
    scale_color_gradientn(colors = my_palette, breaks = col_breaks)+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  library(cowplot)
  print(plot_grid(g1,g2,nrow=1,align="h", rel_widths=c(1,1.5)))
}

Color_VAFratio_on_phylo <- function(input_tbl, this_db, this_id, nwk_tbl){
  #color setting
  library(colorRamps)
  total_break=100
  start_value=-1
  max_value=1
  change_point1 = 0
  col_breaks = seq(start_value, max_value, length.out=total_break)
  a <- length(col_breaks[col_breaks < change_point1] )
  rampCol1 <- colorRampPalette(c("blue", "gray"))(n=a)
  rampCol2 <- colorRampPalette(c("gray", "red"))(n=total_break-a)
  ratio_palette<- c(rampCol1, rampCol2)
  
  
  #print legend
  if(T){
    plot.new()
    plot(c(0,2),c(start_value,max_value),type = 'n', axes = F,xlab = '', ylab = '', main = 'VAF(%)')
    colfunc1 <- rev(my_palette)
    legend_image1 <- as.raster(matrix(colfunc1, ncol=1))
    rasterImage(legend_image1, 0, start_value, 0.3, max_value)
    break_bin=1
    text(x=0.5, y = log10(c(0,0.0001, 0.001, 0.01, 0.1,1)+0.00001), labels = c(0,0.0001, 0.001, 0.01, 0.1,1)*100, cex=2)
  }
  
  #make vaf tbl
  print(this_db)
  f_input_tbl <- input_tbl %>% filter(deadbody == this_db & bulk_id == this_id) 
  f_input_tbl
  f_input_tbl <- f_input_tbl %>% mutate(var_g_id = paste0(lineage_id,'_',rank))
  mf_input_tbl <- f_input_tbl %>% group_by(lineage_id, var_g_id, expected_VAF, rank, n_vars_same_rank) %>% summarise(mean_vaf = mean(vaf))
  mf_input_tbl <- mf_input_tbl %>% mutate(vaf_ratio = mean_vaf/expected_VAF) %>% mutate(log_vaf_ratio = log10(vaf_ratio + 0.01))
  
  ratio_split_tbl <- tibble(lineage_id = as.character(), ratio_rank = as.character(), log_vaf_ratio=as.numeric())
  for (i in 1:nrow(mf_input_tbl)){
    t_id <- mf_input_tbl$lineage_id[i]
    t_rank <- mf_input_tbl$rank[i]
    t_same <- mf_input_tbl$n_vars_same_rank[i]
    t_ratio <- mf_input_tbl$log_vaf_ratio[i]
    ratio_ids <- paste0('RATIO',seq(t_rank, t_rank + t_same -1))
    t_tbl <- tibble(lineage_id = t_id, ratio_rank = ratio_ids, log_vaf_ratio = t_ratio)
    ratio_split_tbl <- bind_rows(ratio_split_tbl, t_tbl)
  }
  annot_tbl <- ratio_split_tbl %>% spread(key=ratio_rank, value=log_vaf_ratio)
  current_cols <- annot_tbl %>% select(starts_with("RATIO")) %>% colnames()
  additional_cols <- setdiff(paste0("RATIO", 1:20), current_cols)
  annot_tbl <- cbind(annot_tbl, tibble(additional_cols) %>% mutate(value=NA) %>% spread(key=additional_cols, value=value)) %>% as.tibble()
  annot_tbl <- left_join(annot_tbl, f_input_tbl %>% select(lineage_id, n_pointmt) %>% unique())
  
  #draw colored phylogenetic tree
  library(ggtree)
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  g <- ggtree(tree) %<+% annot_tbl + theme_tree2()+
    coord_cartesian(xlim = c(-1,30))+
    geom_segment2(aes(subset = is.na(RATIO1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y,color = RATIO1) ,size=2)+
    geom_segment2(aes(subset = is.na(RATIO2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = RATIO2), size=2)+
    geom_segment2(aes(subset = is.na(RATIO3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = RATIO3), size=2)+
    geom_segment2(aes(subset = is.na(RATIO4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = RATIO4), size=2)+
    geom_segment2(aes(subset = is.na(RATIO5) == F, x=x-n_pointmt+4, xend = x -n_pointmt +5, y = y, yend = y, color = RATIO5), size=2)+
    geom_segment2(aes(subset = is.na(RATIO6) == F, x=x-n_pointmt+5, xend = x -n_pointmt +6, y = y, yend = y, color = RATIO6), size=2)+
    geom_segment2(aes(subset = is.na(RATIO7) == F, x=x-n_pointmt+6, xend = x -n_pointmt +7, y = y, yend = y, color = RATIO7), size=2)+
    geom_segment2(aes(subset = is.na(RATIO8) == F, x=x-n_pointmt+7, xend = x -n_pointmt +8, y = y, yend = y, color = RATIO8), size=2)+
    geom_segment2(aes(subset = is.na(RATIO9) == F, x=x-n_pointmt+8, xend = x -n_pointmt +9, y = y, yend = y, color = RATIO9), size=2)+
    geom_segment2(aes(subset = is.na(RATIO10) == F, x=x-n_pointmt+9, xend = x -n_pointmt +10, y = y, yend = y, color = RATIO10), size=2)+
    geom_segment2(aes(subset = is.na(RATIO11) == F, x=x-n_pointmt+10, xend = x -n_pointmt +11, y = y, yend = y, color = RATIO11), size=2)+
    geom_segment2(aes(subset = is.na(RATIO12) == F, x=x-n_pointmt+11, xend = x -n_pointmt +12, y = y, yend = y, color = RATIO12), size=2)+
    geom_segment2(aes(subset = is.na(RATIO13) == F, x=x-n_pointmt+12, xend = x -n_pointmt +13, y = y, yend = y, color = RATIO13), size=2)+
    geom_segment2(aes(subset = is.na(RATIO14) == F, x=x-n_pointmt+13, xend = x -n_pointmt +14, y = y, yend = y, color = RATIO14), size=2)+
    geom_segment2(aes(subset = is.na(RATIO15) == F, x=x-n_pointmt+14, xend = x -n_pointmt +15, y = y, yend = y, color = RATIO15), size=2)+
    geom_segment2(aes(subset = is.na(RATIO16) == F, x=x-n_pointmt+15, xend = x -n_pointmt +16, y = y, yend = y, color = RATIO16), size=2)+
    geom_segment2(aes(subset = is.na(RATIO17) == F, x=x-n_pointmt+16, xend = x -n_pointmt +17, y = y, yend = y, color = RATIO17), size=2)+
    geom_segment2(aes(subset = is.na(RATIO18) == F, x=x-n_pointmt+17, xend = x -n_pointmt +18, y = y, yend = y, color = RATIO18), size=2)+
    geom_segment2(aes(subset = is.na(RATIO19) == F, x=x-n_pointmt+18, xend = x -n_pointmt +19, y = y, yend = y, color = RATIO19), size=2)+
    geom_segment2(aes(subset = is.na(RATIO20) == F, x=x-n_pointmt+19, xend = x -n_pointmt +20, y = y, yend = y, color = RATIO20), size=2)+
    scale_color_gradientn(colors = ratio_palette, breaks = col_breaks)+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  print(g)
  
}



Draw_phylo_tissue <-function(fm_vaf_tbl,lng_dt, this_db){
  max_VAF_len=20
  res_tbl <- as.tibble(data.frame(matrix(vector(), 0, max_VAF_len + 2,
                                         dimnames=list(c(), c('sample_id','lineage_id2',paste0('VAF',1:max_VAF_len)))),
                                  stringsAsFactors=F))
  n=0
  for (sam_id in unique(fm_vaf_tbl$sample_id)){
    print(sam_id)
    id_list <- fm_vaf_tbl %>% filter(sample_id == sam_id) %>% .$lineage_id2 %>% unique()
    id_list <- id_list[grepl('germline',id_list)==F]
    ffm_vaf_tbl <- fm_vaf_tbl %>% filter(sample_id == sam_id)
    for (lid in id_list){
      n=n+1
      VAF_v <- ffm_vaf_tbl %>% filter(lineage_id2 == lid) %>% arrange(desc(VAF)) %>% .$VAF
      VAF_v <- round(VAF_v*100,2)
      NA_v <- rep(NA, max_VAF_len - length(VAF_v))
      res_tbl[n,] <- c(sam_id,lid,VAF_v, NA_v)
    }
  }
  res_tbl <- res_tbl %>% mutate_at(vars(starts_with('VAF')), list(as.numeric)) 
  early_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early') %>% .$lineage_id2
  res_tbl <- res_tbl %>% filter(lineage_id2 %in% early_ids)
  med_tbl <- left_join(res_tbl, pmeta_dt) %>% group_by(tissue_class, lineage_id2) %>% summarise_at(vars(starts_with('VAF')), list(median=median))
  colnames(med_tbl) <- gsub('_median','',colnames(med_tbl))
  library(colorRamps)
  total_break=100
  col_breaks = seq(0, 100, length.out=total_break)
  a <- length(col_breaks[col_breaks < 1] )  #change color at 1
  b <- length(col_breaks[col_breaks < 5])-a  #change color at 5
  rampCol1 <- colorRampPalette(c("gray", "yellow"))(n=a)
  rampCol2 <- colorRampPalette(c("yellow","red"))(n=b)
  rampCol3 <- colorRampPalette(c("red", "blue"))(n=total_break-(a+b))
  my_palette<- c(rampCol1, rampCol2, rampCol3)
  
  p_list <- list()
  n=0
  for (this_tis in unique(med_tbl$tissue_class)){
    n = n+1
    tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
    m_lng_dt <- lng_dt %>% filter(deadbody == this_db) %>% rename(taxa = lineage_id) %>% select(-deadbody)
    m_lng_dt$Source_class2[m_lng_dt$time_of_lineage == 'early'] <- NA
    m_lng_dt <- left_join(m_lng_dt, med_tbl %>% filter(tissue_class == this_tis))
    p_list[[n]] <- ggtree(tree) %<+% m_lng_dt + theme_tree2()+
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
      ggtitle(this_tis)+
      scale_color_gradientn(colors = my_palette, breaks = col_breaks)+
      theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  }
  library(cowplot)
  plot_grid(plotlist=p_list, ncol=3)
}


Draw_Phylo_Epi_Derm <- function(fm_vaf_tbl, lng_dt, this_db){
  max_VAF_len=20
  res_tbl <- as.tibble(data.frame(matrix(vector(), 0, max_VAF_len + 2,
                                         dimnames=list(c(), c('sample_id','lineage_id2',paste0('VAF',1:max_VAF_len)))),
                                  stringsAsFactors=F))
  n=0
  for (sam_id in unique(fm_vaf_tbl$sample_id)){
    print(sam_id)
    id_list <- fm_vaf_tbl %>% filter(sample_id == sam_id) %>% .$lineage_id2 %>% unique()
    id_list <- id_list[grepl('germline',id_list)==F]
    ffm_vaf_tbl <- fm_vaf_tbl %>% filter(sample_id == sam_id)
    for (lid in id_list){
      n=n+1
      VAF_v <- ffm_vaf_tbl %>% filter(lineage_id2 == lid) %>% arrange(desc(VAF)) %>% .$VAF
      VAF_v <- round(VAF_v*100,2)
      NA_v <- rep(NA, max_VAF_len - length(VAF_v))
      res_tbl[n,] <- c(sam_id,lid,VAF_v, NA_v)
    }
  }
  res_tbl <- res_tbl %>% mutate_at(vars(starts_with('VAF')), list(as.numeric)) 
  #res_tbl <- res_tbl %>% select(sample_id, lineage_id2, colnames(res_tbl[3:ncol(res_tbl)])[colSums(is.na(res_tbl[3:ncol(res_tbl)])==F) > 0]) #remove columns with all NAs
  early_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early') %>% .$lineage_id2
  res_tbl <- res_tbl %>% filter(lineage_id2 %in% early_ids)
  med_tbl <- left_join(res_tbl, pmeta_dt) %>% group_by(tissue_class2, lineage_id2) %>% summarise_at(vars(starts_with('VAF')), list(median=median))
  med_tbl$tissue_class2 %>% unique()
  m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db) ,med_tbl %>% filter(tissue_class2 == 'skin_derm'))
  m_lng_dt <- m_lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt$Source_class2[m_lng_dt$time_of_lineage == 'early'] <- NA
  colnames(m_lng_dt) <- gsub('_median','',colnames(m_lng_dt))
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  g1 <- ggtree(tree) %<+% m_lng_dt + theme_tree2()+
    coord_cartesian(xlim = c(-1,35))+
    geom_text2(aes(subset = is.na(Source_class2) == F, x=x-n_pointmt+5, y = y+0.25, label = Source_class2), size=2.5, color="red")+
    geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt.x +1, y = y, yend = y, color = VAF1))+
    geom_text2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt.x, y = y, label = VAF1), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt.x+1, xend = x -n_pointmt.x +2, y = y, yend = y, color = VAF2))+
    geom_text2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt.x+1, y = y, label = VAF2), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt.x+2, xend = x -n_pointmt.x +3, y = y, yend = y, color = VAF3))+
    geom_text2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt.x+2, y = y, label = VAF3), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt.x+3, xend = x -n_pointmt.x +4, y = y, yend = y, color = VAF4))+
    geom_text2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt.x+3, y = y, label = VAF4), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt.x+4, xend = x -n_pointmt.x +5, y = y, yend = y, color = VAF5))+
    geom_text2(aes(subset = is.na(VAF5) == F, x=x-n_pointmt.x+4, y = y, label = VAF5), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt.x+5, xend = x -n_pointmt.x +6, y = y, yend = y, color = VAF6))+
    geom_text2(aes(subset = is.na(VAF6) == F, x=x-n_pointmt.x+5, y = y, label = VAF6), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt.x+6, xend = x -n_pointmt.x +7, y = y, yend = y, color = VAF7))+
    geom_text2(aes(subset = is.na(VAF7) == F, x=x-n_pointmt.x+6, y = y, label = VAF7), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt.x+7, xend = x -n_pointmt.x +8, y = y, yend = y, color = VAF8))+
    geom_text2(aes(subset = is.na(VAF8) == F, x=x-n_pointmt.x+7, y = y, label = VAF8), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt.x+8, xend = x -n_pointmt.x +9, y = y, yend = y, color = VAF9))+
    geom_text2(aes(subset = is.na(VAF9) == F, x=x-n_pointmt.x+8, y = y, label = VAF9), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt.x+9, xend = x -n_pointmt.x +10, y = y, yend = y, color = VAF10))+
    geom_text2(aes(subset = is.na(VAF10) == F, x=x-n_pointmt.x+9, y = y, label = VAF10), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt.x+10, xend = x -n_pointmt.x +11, y = y, yend = y, color = VAF11))+
    geom_text2(aes(subset = is.na(VAF11) == F, x=x-n_pointmt.x+10, y = y, label = VAF11), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt.x+11, xend = x -n_pointmt.x +12, y = y, yend = y, color = VAF12))+
    geom_text2(aes(subset = is.na(VAF12) == F, x=x-n_pointmt.x+11, y = y, label = VAF12), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt.x+12, xend = x -n_pointmt.x +13, y = y, yend = y, color = VAF13))+
    geom_text2(aes(subset = is.na(VAF13) == F, x=x-n_pointmt.x+12, y = y, label = VAF13), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt.x+13, xend = x -n_pointmt.x +14, y = y, yend = y, color = VAF14))+
    geom_text2(aes(subset = is.na(VAF14) == F, x=x-n_pointmt.x+13, y = y, label = VAF14), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt.x+14, xend = x -n_pointmt.x +15, y = y, yend = y, color = VAF15))+
    geom_text2(aes(subset = is.na(VAF15) == F, x=x-n_pointmt.x+14, y = y, label = VAF15), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt.x+15, xend = x -n_pointmt.x +16, y = y, yend = y, color = VAF16))+
    geom_text2(aes(subset = is.na(VAF16) == F, x=x-n_pointmt.x+15, y = y, label = VAF16), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt.x+16, xend = x -n_pointmt.x +17, y = y, yend = y, color = VAF17))+
    geom_text2(aes(subset = is.na(VAF17) == F, x=x-n_pointmt.x+16, y = y, label = VAF17), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt.x+17, xend = x -n_pointmt.x +18, y = y, yend = y, color = VAF18))+
    geom_text2(aes(subset = is.na(VAF18) == F, x=x-n_pointmt.x+17, y = y, label = VAF18), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt.x+18, xend = x -n_pointmt.x +19, y = y, yend = y, color = VAF19))+
    geom_text2(aes(subset = is.na(VAF19) == F, x=x-n_pointmt.x+18, y = y, label = VAF19), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt.x+19, xend = x -n_pointmt.x +20, y = y, yend = y, color = VAF20))+
    geom_text2(aes(subset = is.na(VAF20) == F, x=x-n_pointmt.x+19, y = y, label = VAF20), size=2.5)+
    ggtitle(paste0(this_db, ', median VAF of skin_derm'))+
    scale_color_gradientn(colors = my_palette, breaks = col_breaks)+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db), med_tbl %>% filter(tissue_class2 == 'skin_epi'))
  m_lng_dt <- m_lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt$Source_class2[m_lng_dt$time_of_lineage == 'early'] <- NA
  colnames(m_lng_dt) <- gsub('_median','',colnames(m_lng_dt))
  g2 <- ggtree(tree) %<+% m_lng_dt + theme_tree2()+
    coord_cartesian(xlim = c(-1,35))+
    geom_text2(aes(subset = is.na(Source_class2) == F, x=x-n_pointmt+5, y = y+0.25, label = Source_class2), size=2.5, color="red")+
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
    ggtitle(paste0(this_db, ', median VAF of skin_epi'))+
    scale_color_gradientn(colors = my_palette, breaks = col_breaks)+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  plot_grid(g1,g2, nrow=1)
}


Draw_maxVAF_1v1 <- function(m_vaf_tbl, this_db, lid1,lid2){
  fm_vaf_tbl <- m_vaf_tbl %>% filter(lineage_id2 %in% c(lid1, lid2)) %>% group_by(lineage_id2, sample_id) %>% summarise(maxVAF = max(VAF))
  x_order <- fm_vaf_tbl %>% spread(key=lineage_id2, value=maxVAF) %>% mutate(L12ratio =  get(lid1)/get(lid2)) %>% arrange(desc(L12ratio)) %>% .$sample_id
  c <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early_to_late') %>% nrow()
  median_v <- fm_vaf_tbl %>% filter(lineage_id2 == lid1) %>% .$maxVAF %>% median()
  fm_vaf_tbl$lineage_id2 <- factor(fm_vaf_tbl$lineage_id2, levels = c(lid2,lid1))
  a <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early_to_late'& grepl(lid1, lineage_id2)) %>% nrow()
  b <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early_to_late'& grepl(lid2, lineage_id2)) %>% nrow()
 
  g1 <- ggplot(fm_vaf_tbl, aes(x=sample_id, y=maxVAF))+
    geom_bar(aes(fill = lineage_id2), stat="identity")+
    geom_hline(yintercept= 0.48*a/c, color="red")+
    geom_hline(yintercept= median_v, color="blue")+
    scale_x_discrete(limits = x_order)+
    scale_y_continuous(limits=c(0,0.6))+
    ggtitle(paste0(this_db,", red line = expected value of ",lid1,", blue line = median value of ",lid1))+ylab("maxVAF")+xlab("samples")+
    theme_syp +theme(axis.text.x= element_text(angle=45, hjust=1, size=8))
  mfm_vaf_tbl <- fm_vaf_tbl %>% spread(key=lineage_id2, value=maxVAF) %>% mutate(r1=get(lid1)/(get(lid1)+get(lid2)),
                                                                  r2=get(lid2)/(get(lid1)+get(lid2))) %>%
    gather(-sample_id, -lid1, -lid2, key='lineage',value='ratio') 
  x_order2 <- fm_vaf_tbl %>% spread(key=lineage_id2, value=maxVAF) %>% mutate(r1=get(lid1)/(get(lid1)+get(lid2)),
                                                                  r2=get(lid2)/(get(lid1)+get(lid2))) %>% arrange(desc(r1)) %>% .$sample_id
  median_v2 <- fm_vaf_tbl %>% spread(key=lineage_id2, value=maxVAF) %>% mutate(r1=get(lid1)/(get(lid1)+get(lid2)),
                                                                  r2=get(lid2)/(get(lid1)+get(lid2))) %>% .$r1 %>% median
  mfm_vaf_tbl$lineage[mfm_vaf_tbl$lineage == 'r1'] <- lid1
  mfm_vaf_tbl$lineage[mfm_vaf_tbl$lineage == 'r2'] <- lid2
  mfm_vaf_tbl$lineage = factor(mfm_vaf_tbl$lineage, levels=c(lid2,lid1))
  g2 <- ggplot(mfm_vaf_tbl, aes(x=sample_id, y=ratio, fill=lineage))+
    geom_bar(stat="identity")+
    geom_hline(yintercept=a/(a+b), color="red")+
    geom_hline(yintercept= median_v2, color="blue")+
    scale_x_discrete(limits= x_order2)+
    ggtitle(paste0(this_db,", red line = expected value of ",lid1,", blue line = median value of ",lid1))+ylab("VAFratio")+xlab("samples")+
    theme_syp+theme(axis.text.x = element_text(angle=45, hjust=1,size=8))
  library(cowplot)
  plot_grid(g1,g2,nrow=1)
}


Draw_Lineage_Plot<- function(m_vaf_tbl, db, lid){
  m_lng_dt <- lng_dt %>% filter(deadbody == db) %>% mutate(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt$Source_class2[is.na(m_lng_dt$Source_class2)] <- 'core'
  m_lng_dt <- m_lng_dt %>% mutate(marking = ifelse(lineage_id2 == lid, 'Y', NA))
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
  g1 <- ggtree(tree, aes(color=Source_class2)) %<+% m_lng_dt +theme_tree2() +
    #geom_text2(aes(x=x-n_pointmt, y = y, label =taxa), size=2.5, color="red")+
    geom_text2(aes(x=x-n_pointmt+2, y = y, label =Source3), size=2.5, color="black")+
    geom_segment2(aes(subset = is.na(marking) == F, x=x-n_pointmt, xend = x, y = y, yend = y),color = "red", size=2)+
    coord_cartesian(xlim = c(-1,30))+
    ggtitle(lid)+
    scale_color_manual(values = anatomy_cl3_pal)+
    theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))

  fm_vaf_tbl <- m_vaf_tbl
  ffm_vaf_tbl <- fm_vaf_tbl %>% filter(lineage_id2 == lid) 
  ffm_vaf_tbl$var_id %>% unique() %>% length() -> max_VAF_len
  res_tbl <- as.tibble(data.frame(matrix(vector(), 0, max_VAF_len + 1,
                                         dimnames=list(c(), c('sample_id',paste0('VAF',1:max_VAF_len)))),
                                  stringsAsFactors=F))
  n=0
  id_list <- unique(ffm_vaf_tbl$sample_id)
  for (sid in id_list){
    n=n+1
    VAF_v <- ffm_vaf_tbl %>% filter(sample_id == sid) %>% arrange(desc(VAF)) %>% .$VAF
    VAF_v <- round(VAF_v*100,2)
    NA_v <- rep(NA, max_VAF_len - length(VAF_v))
    res_tbl[n,] <- c(sid,VAF_v, NA_v)
  }
  res_tbl <- res_tbl %>% mutate_at(vars(starts_with('VAF')), list(as.numeric))
  
  m_res_tbl <- res_tbl %>% mutate(VAFmax = ifelse(VAF1 == 0, 1, VAF1))
  m_res_tbl <- m_res_tbl[2:(max_VAF_len+1)]/m_res_tbl$VAFmax
  m_res_tbl$sample_id <- res_tbl$sample_id
  
  tmp_res_tbl <- m_res_tbl
  colnames(tmp_res_tbl)[max_VAF_len] <- 'VAFmax'
  mm_res_tbl <- tmp_res_tbl %>% as.tibble() %>% gather(-sample_id, key='VAF_order', value='VAF') %>% mutate(VAF = as.numeric(VAF))
  mm_res_tbl <- left_join(mm_res_tbl, pmeta_dt)
  mm_res_tbl <- mm_res_tbl %>% mutate(anatomy_class4 = paste(Source_side, anatomy_class2, sep='_'))
  mm_res_tbl$anatomy_tissue_class2 <- gsub('NA_','',paste(mm_res_tbl$Source_side, mm_res_tbl$anatomy_tissue_class, sep="_"))
  mm_res_tbl <- left_join(mm_res_tbl,tmp_res_tbl %>% mutate(VAF_maxR = VAFmax/VAF1) %>% select(sample_id, VAF_maxR))
  mm_res_tbl$VAF_order <- factor(mm_res_tbl$VAF_order, levels = c(paste0('VAF', 1:(max_VAF_len-1)),'VAFmax'))
  
  g2 <- ggplot(mm_res_tbl, aes(x=VAF_order, y=VAF, color = anatomy_class3))+
    geom_point(alpha=0.3)+
    geom_line(aes(group=sample_id), alpha=0.3)+
    #geom_text(data= subset(mm_res_tbl, VAF_maxR > 0), aes(x=max_VAF_len+1, y=VAF_maxR, label=anatomy_tissue_class2), size=3)+
    scale_x_discrete(limits=c(paste0('VAF',1:(max_VAF_len-1)),'VAFmax'))+
    scale_color_manual(values=anatomy_cl3_pal)+
    xlab("Variants (VAF decreasing order)")+ylab("Relative VAF")+
    ggtitle(lid)+
    theme_syp+theme(axis.text.x = element_blank())
  
  #draw group plot
  g_tbl <- left_join(m_res_tbl, pmeta_dt) %>% group_by(anatomy_class3) %>% summarise_at(vars(starts_with('VAF')), 
                                                                                        list(min=min, Q1= ~quantile(.,probs=0.25, na.rm=T), median=median, 
                                                                                             Q3= ~quantile(., probs=0.75,na.rm=T), max=max, mean=mean))
  tmp1 <- g_tbl %>% select(anatomy_class3, ends_with('_median')) %>% gather(-anatomy_class3, key='VAForder', value='VAF_median') %>% mutate(VAForder = gsub('_median','',VAForder))
  tmp2 <- g_tbl %>% select(anatomy_class3, ends_with('_Q1')) %>% gather(-anatomy_class3, key='VAForder', value='VAF_Q1') %>% mutate(VAForder = gsub('_Q1','',VAForder))
  tmp3 <- g_tbl %>% select(anatomy_class3, ends_with('_Q3')) %>% gather(-anatomy_class3, key='VAForder', value='VAF_Q3') %>% mutate(VAForder = gsub('_Q3','',VAForder))
  m_g_tbl <- left_join(left_join(tmp1, tmp2), tmp3) %>% mutate(x_pos = as.numeric(gsub('VAF','', VAForder)))
  pos_adj_v <- seq(-0.2,0.2,length.out=7)
  names(pos_adj_v) <- c('HN', 'internal_organ','trunk','lt_LE','lt_UE','rt_LE','rt_UE')
  m_g_tbl$x_adj <- pos_adj_v[m_g_tbl$anatomy_class3]
  
  g3 <- ggplot(m_g_tbl, aes(x=x_pos+x_adj, y=VAF_median, color= anatomy_class3))+
    geom_point()+
    geom_line(aes(group = anatomy_class3))+
    geom_errorbar(aes(ymin=VAF_Q1, ymax=VAF_Q3))+
    scale_color_manual(values = anatomy_cl3_pal)+
    scale_x_continuous(breaks=seq(1,11,1), labels = c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th'))+
    xlab("Variants (VAF decreasing order)")+ylab("Relative VAF")+
    ggtitle(lid)+
    theme_syp+theme(axis.text.x = element_blank())
  
  g_tbl <- left_join(res_tbl, pmeta_dt) %>% group_by(anatomy_class3) %>% summarise_at(vars(starts_with('VAF')), 
                                                                                        list(min=min, Q1= ~quantile(.,probs=0.25, na.rm=T), median=median, 
                                                                                             Q3= ~quantile(., probs=0.75,na.rm=T), max=max, mean=mean))
  tmp1 <- g_tbl %>% select(anatomy_class3, ends_with('_median')) %>% gather(-anatomy_class3, key='VAForder', value='VAF_median') %>% mutate(VAForder = gsub('_median','',VAForder))
  tmp2 <- g_tbl %>% select(anatomy_class3, ends_with('_Q1')) %>% gather(-anatomy_class3, key='VAForder', value='VAF_Q1') %>% mutate(VAForder = gsub('_Q1','',VAForder))
  tmp3 <- g_tbl %>% select(anatomy_class3, ends_with('_Q3')) %>% gather(-anatomy_class3, key='VAForder', value='VAF_Q3') %>% mutate(VAForder = gsub('_Q3','',VAForder))
  m_g_tbl <- left_join(left_join(tmp1, tmp2), tmp3) %>% mutate(x_pos = as.numeric(gsub('VAF','', VAForder)))
  pos_adj_v <- seq(-0.2,0.2,length.out=7)
  names(pos_adj_v) <- c('HN', 'internal_organ','trunk','lt_LE','lt_UE','rt_LE','rt_UE')
  m_g_tbl$x_adj <- pos_adj_v[m_g_tbl$anatomy_class3]
  g4 <- ggplot(m_g_tbl, aes(x=x_pos+x_adj, y=VAF_median, color= anatomy_class3))+
    geom_point()+
    geom_line(aes(group = anatomy_class3))+
    geom_errorbar(aes(ymin=VAF_Q1, ymax=VAF_Q3))+
    scale_color_manual(values = anatomy_cl3_pal)+
    scale_x_continuous(breaks=seq(1,11,1), labels = c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th'))+
    xlab("Variants (VAF decreasing order)")+ylab("VAF(%)")+
    ggtitle(lid)+
    theme_syp+theme(axis.text.x = element_blank())
  
  
  g_tbl2 <- left_join(res_tbl, pmeta_dt) %>% group_by(tissue_class2) %>% summarise_at(vars(starts_with('VAF')), 
                                                                                        list(min=min, Q1= ~quantile(.,probs=0.25, na.rm=T), median=median, 
                                                                                             Q3= ~quantile(., probs=0.75,na.rm=T), max=max, mean=mean))
  tmp1 <- g_tbl2 %>% select(tissue_class2, ends_with('_median')) %>% gather(-tissue_class2, key='VAForder', value='VAF_median') %>% mutate(VAForder = gsub('_median','',VAForder))
  tmp2 <- g_tbl2 %>% select(tissue_class2, ends_with('_Q1')) %>% gather(-tissue_class2, key='VAForder', value='VAF_Q1') %>% mutate(VAForder = gsub('_Q1','',VAForder))
  tmp3 <- g_tbl2 %>% select(tissue_class2, ends_with('_Q3')) %>% gather(-tissue_class2, key='VAForder', value='VAF_Q3') %>% mutate(VAForder = gsub('_Q3','',VAForder))
  m_g_tbl2 <- left_join(left_join(tmp1, tmp2), tmp3) %>% mutate(x_pos = as.numeric(gsub('VAF','', VAForder)))
  pos_adj_v <- seq(-0.2,0.2,length.out=length(unique(m_g_tbl2$tissue_class2)))
  names(pos_adj_v) <- unique(m_g_tbl2$tissue_class2)
  m_g_tbl2$x_adj <- pos_adj_v[m_g_tbl2$tissue_class2]
  g5 <- ggplot(m_g_tbl2, aes(x=x_pos+x_adj, y=VAF_median, color= tissue_class2))+
    geom_point()+
    geom_line(aes(group = tissue_class2))+
    geom_errorbar(aes(ymin=VAF_Q1, ymax=VAF_Q3))+
    scale_x_continuous(breaks=seq(1,11,1), labels = c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th'))+
    xlab("Variants (VAF decreasing order)")+ylab("VAF(%)")+
    ggtitle(lid)+
    theme_syp+theme(axis.text.x = element_blank())
  g6 <- ggplot(subset(m_g_tbl2, tissue_class2 %in% c("skin_derm","skin_epi")), aes(x=x_pos+x_adj, y=VAF_median, color= tissue_class2))+
    geom_point()+
    geom_line(aes(group = tissue_class2))+
    geom_errorbar(aes(ymin=VAF_Q1, ymax=VAF_Q3))+
    scale_x_continuous(breaks=seq(1,11,1), labels = c('1st','2nd','3rd','4th','5th','6th','7th','8th','9th','10th','11th'))+
    xlab("Variants (VAF decreasing order)")+ylab("VAF(%)")+
    ggtitle(lid)+
    theme_syp+theme(axis.text.x = element_blank())
  
  print(plot_grid(g1,g5,g6,g2,g3,g4, nrow=2))
}

Draw_phylo_timing_skin <-function(fm_vaf_tbl,lng_dt, this_db, cutoff){
    max_VAF_len=20
    res_tbl <- as.tibble(data.frame(matrix(vector(), 0, max_VAF_len + 2,
                                           dimnames=list(c(), c('sample_id','lineage_id2',paste0('VAF',1:max_VAF_len)))),
                                    stringsAsFactors=F))
    n=0
    for (sam_id in unique(fm_vaf_tbl$sample_id)){
      print(sam_id)
      id_list <- fm_vaf_tbl %>% filter(sample_id == sam_id) %>% .$lineage_id2 %>% unique()
      id_list <- id_list[grepl('germline',id_list)==F]
      ffm_vaf_tbl <- fm_vaf_tbl %>% filter(sample_id == sam_id)
      for (lid in id_list){
        n=n+1
        VAF_v <- ffm_vaf_tbl %>% filter(lineage_id2 == lid) %>% arrange(desc(VAF)) %>% .$VAF
        VAF_v <- round(VAF_v*100,2)
        NA_v <- rep(NA, max_VAF_len - length(VAF_v))
        res_tbl[n,] <- c(sam_id,lid,VAF_v, NA_v)
      }
    }
    res_tbl <- res_tbl %>% mutate_at(vars(starts_with('VAF')), list(as.numeric)) 
    #res_tbl <- res_tbl %>% select(sample_id, lineage_id2, colnames(res_tbl[3:ncol(res_tbl)])[colSums(is.na(res_tbl[3:ncol(res_tbl)])==F) > 0]) #remove columns with all NAs
    early_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early') %>% .$lineage_id2
    res_tbl <- res_tbl %>% filter(lineage_id2 %in% early_ids)
    med_tbl <- left_join(res_tbl, pmeta_dt) %>% group_by(tissue_class2, lineage_id2) %>% summarise_at(vars(starts_with('VAF')), list(median=median))
    res_tbl2 <- as.tibble(data.frame(matrix(vector(), 0, max_VAF_len + 1,
                                           dimnames=list(c(), c('lineage_id2',paste0('VAF',1:max_VAF_len)))),
                                    stringsAsFactors=F))
    n=0
    for (lid in unique(med_tbl$lineage_id2)){
      n=n+1
      print(lid)
      res_v=c()
      for (coln in colnames(med_tbl)[3:ncol(med_tbl)]){
        VAFderm <- as.numeric(med_tbl[which(med_tbl$lineage_id2 == lid & med_tbl$tissue_class2 == 'skin_derm'), which(colnames(med_tbl) == coln)])
        VAFepi <- as.numeric(med_tbl[which(med_tbl$lineage_id2 == lid & med_tbl$tissue_class2 == 'skin_epi'), which(colnames(med_tbl) == coln)])
        if (is.na(VAFderm) ==T | is.na(VAFepi) == T){
         VARres=NA 
        }else if(VAFderm  > cutoff & VAFepi > cutoff){
          VARres=2
        }else if (VAFderm > cutoff & VAFepi <= cutoff){
          VARres=1
        }else { VARres=0}
        res_v = c(res_v, VARres)
      }
      res_tbl2[n,] <- c(lid, res_v)
      }
    res_tbl2
    m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db) ,res_tbl2)
    m_lng_dt <- m_lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
    m_lng_dt$Source_class2[m_lng_dt$time_of_lineage == 'early'] <- NA
    colnames(m_lng_dt) <- gsub('_median','',colnames(m_lng_dt))
    tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
    ggtree(tree) %<+% m_lng_dt + theme_tree2()+
      coord_cartesian(xlim = c(-1,35))+
      #geom_text2(aes(subset = is.na(Source_class2) == F, x=x-n_pointmt+5, y = y+0.25, label = Source_class2), size=2.5, color="red")+
      geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y, color = VAF1), size=2)+
      #geom_text2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, y = y, label = VAF1), size=2.5)+
      geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2), size=2)+
      #geom_text2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, y = y, label = VAF2), size=2.5)+
      geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3), size=2)+
      #geom_text2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, y = y, label = VAF3), size=2.5)+
      geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4), size=2)+
      #geom_text2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, y = y, label = VAF4), size=2.5)+
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
      ggtitle(paste0(this_db, ', median VAF > ',cutoff,'%, both derm-epi, only derm, none'))+
      theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
}

Draw_phylo_timing_SampleAlign <-function(fm_vaf_tbl,lng_dt, this_db, cutoff1, cutoff2){
  max_VAF_len=20
  res_tbl <- as.tibble(data.frame(matrix(vector(), 0, max_VAF_len + 2,
                                         dimnames=list(c(), c('sample_id','lineage_id2',paste0('VAF',1:max_VAF_len)))),
                                  stringsAsFactors=F))
  n=0
  for (sam_id in unique(fm_vaf_tbl$sample_id)){
    print(sam_id)
    id_list <- fm_vaf_tbl %>% filter(sample_id == sam_id) %>% .$lineage_id2 %>% unique()
    id_list <- id_list[grepl('germline',id_list)==F]
    ffm_vaf_tbl <- fm_vaf_tbl %>% filter(sample_id == sam_id)
    for (lid in id_list){
      n=n+1
      VAF_v <- ffm_vaf_tbl %>% filter(lineage_id2 == lid) %>% arrange(desc(VAF)) %>% .$VAF
      VAF_v <- round(VAF_v*100,2)
      NA_v <- rep(NA, max_VAF_len - length(VAF_v))
      res_tbl[n,] <- c(sam_id,lid,VAF_v, NA_v)
    }
  }
  res_tbl <- res_tbl %>% mutate_at(vars(starts_with('VAF')), list(as.numeric)) 
  #res_tbl <- res_tbl %>% select(sample_id, lineage_id2, colnames(res_tbl[3:ncol(res_tbl)])[colSums(is.na(res_tbl[3:ncol(res_tbl)])==F) > 0]) #remove columns with all NAs
  early_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early') %>% .$lineage_id2
  res_tbl <- res_tbl %>% filter(lineage_id2 %in% early_ids)
  med_tbl <- left_join(res_tbl, pmeta_dt) %>% group_by(lineage_id2) %>% summarise_at(vars(starts_with('VAF')), list(median=median))
  mx <- med_tbl %>% select(starts_with('VAF'))
  mx[mx < cutoff2] <- 0
  mx[mx >= cutoff2 & mx < cutoff1 ] <- 1
  mx[mx >= cutoff1] <- 2
  m_med_tbl <- bind_cols(med_tbl %>% select(lineage_id2), mx %>% as.tibble())
    colnames(m_med_tbl) <- gsub('_median','',colnames(m_med_tbl))
  m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db),m_med_tbl)
  m_lng_dt <- m_lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt$Source_class2[m_lng_dt$time_of_lineage == 'early'] <- NA
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  ggtree(tree) %<+% m_lng_dt + theme_tree2()+
    coord_cartesian(xlim = c(-1,35))+
    #geom_text2(aes(subset = is.na(Source_class2) == F, x=x-n_pointmt+5, y = y+0.25, label = Source_class2), size=2.5, color="red")+
    geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y, color = VAF1), size=2)+
    #geom_text2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, y = y, label = VAF1), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2), size=2)+
    #geom_text2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, y = y, label = VAF2), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3), size=2)+
    #geom_text2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, y = y, label = VAF3), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4), size=2)+
    #geom_text2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, y = y, label = VAF4), size=2.5)+
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
    ggtitle(paste0(this_db, ', median VAF cutoff1: ',cutoff1,'%, cutoff2: ',cutoff2,'%'))+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
}


Draw_phylo_timing_TotalAlign <-function(fm_vaf_tbl,lng_dt, this_db, cutoff_v){
  max_VAF_len=20
  early_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early') %>% .$lineage_id2 
  fm_vaf_tbl <- fm_vaf_tbl %>% filter(lineage_id2 %in% early_ids)
  fm_vaf_tbl <- fm_vaf_tbl %>% filter(is.na(rank_ctie_first) == F)
  rank_tbl <- fm_vaf_tbl %>% filter(dominant_layer == 'mesoderm') %>% group_by(lineage_id2, var_id, rank_ctie_first) %>% summarise(medianVAF = median(VAF)) %>% 
    mutate(rank2 = paste0('VAF', rank_ctie_first)) %>% ungroup() %>% select(-var_id, -rank_ctie_first) %>% spread(key=rank2, value=medianVAF)
  for(id in setdiff(paste0('VAF',1:max_VAF_len), colnames(rank_tbl))){
    rank_tbl <- rank_tbl %>% mutate(this = NA)
    colnames(rank_tbl) <- gsub('this',id,colnames(rank_tbl))
  }
  change_to_group <- function(numb,cutoff_v){
    if(is.na(numb)){
      return(NA)
    }else{
      n_cutoff_v <- sort(c(numb,cutoff_v), decreasing = T)
      return(paste0('G',min(which(n_cutoff_v == numb))))
    }
  }
  m_rank_tbl <- rank_tbl %>% mutate_at(vars(starts_with('VAF')), list(~ map_chr(., ~ change_to_group(.x, cutoff_v)))) 
  m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db),m_rank_tbl)
  m_lng_dt <- m_lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
  #m_lng_dt$Source_class2[m_lng_dt$time_of_lineage == 'early'] <- NA
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  ggtree(tree) %<+% m_lng_dt + theme_tree2()+
    coord_cartesian(xlim = c(-1,35))+
    #geom_text2(aes(subset = is.na(Source_class2) == F, x=x-n_pointmt+5, y = y+0.25, label = Source_class2), size=2.5, color="red")+
    geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y, color = VAF1), size=2)+
    #geom_text2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, y = y, label = VAF1), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2), size=2)+
    #geom_text2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, y = y, label = VAF2), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3), size=2)+
    #geom_text2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, y = y, label = VAF3), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4), size=2)+
    #geom_text2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, y = y, label = VAF4), size=2.5)+
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
    #ggtitle(paste0(this_db, ', median VAF cutoff1: ',cutoff1,'%, cutoff2: ',cutoff2,'%'))+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
}




Draw_phylo_GroupByVAF <-function(
  fm_vaf_tbl,
  lng_dt,
  this_db,
  cutoff1=NULL,
  cutoff2=NULL){
  max_VAF_len=20
  early_ids <- lng_dt %>% filter(deadbody == this_db & time_of_lineage == 'early') %>% .$lineage_id2 
  fm_vaf_tbl <- fm_vaf_tbl %>% filter(lineage_id2 %in% early_ids)
  rank_tbl <- fm_vaf_tbl %>% group_by(lineage_id2, var_id) %>% summarise(medianVAF = median(VAF)) %>% filter(grepl('germline',lineage_id2)==F) %>% 
    group_by(lineage_id2) %>% mutate(rank = rank(-medianVAF, ties.method = "first")) 
  rank_tbl <- rank_tbl %>% mutate(rank2 = paste0('VAF',rank)) %>% select(-var_id, -rank) %>% spread(key=rank2, value=medianVAF)
  for(id in setdiff(paste0('VAF',1:max_VAF_len), colnames(rank_tbl))){
    rank_tbl <- rank_tbl %>% mutate(this = NA)
    colnames(rank_tbl) <- gsub('this',id,colnames(rank_tbl))
  }
  mx <- rank_tbl %>% ungroup() %>% select(starts_with('VAF'))
  mx = mx*100
  if (is.null(cutoff1)){
    mx[mx < cutoff2] <- 0
    mx[mx >= cutoff2] <- 1
  } else{
    mx[mx < cutoff2] <- 0
    mx[mx >= cutoff2 & mx < cutoff1 ] <- 1
    mx[mx >= cutoff1] <- 2
  }
  

  m_rank_tbl <- bind_cols(rank_tbl %>% select(lineage_id2), mx %>% as.tibble())
   m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db),m_rank_tbl)
  m_lng_dt <- m_lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt$Source_class2[m_lng_dt$time_of_lineage == 'early'] <- NA
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  ggtree(tree) %<+% m_lng_dt + theme_tree2()+
    coord_cartesian(xlim = c(-1,35))+
    #geom_text2(aes(subset = is.na(Source_class2) == F, x=x-n_pointmt+5, y = y+0.25, label = Source_class2), size=2.5, color="red")+
    geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y, color = VAF1), size=2)+
    #geom_text2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, y = y, label = VAF1), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2), size=2)+
    #geom_text2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, y = y, label = VAF2), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3), size=2)+
    #geom_text2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, y = y, label = VAF3), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4), size=2)+
    #geom_text2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, y = y, label = VAF4), size=2.5)+
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
    ggtitle(paste0(this_db, ', median VAF cutoff1: ',cutoff1,'%, cutoff2: ',cutoff2,'%'))+
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
}

Draw_phylo_GroupByManual <-function(
  fm_vaf_tbl,
  lng_dt,
  this_db,
  diff_tbl
  ){
  tmp_tbl <- left_join(fm_vaf_tbl, diff_tbl) 
  tmp_tbl <- tmp_tbl %>% mutate(rank2 = paste0('VAF', rank_ctie_first))
  rank_tbl <- tmp_tbl %>% select(lineage_id2, rank2, diff_group) %>% unique() %>%  spread(key=rank2, value=diff_group)
  max_VAF_len=20
  for(id in setdiff(paste0('VAF',1:max_VAF_len), colnames(rank_tbl))){
    rank_tbl <- rank_tbl %>% mutate(this = NA)
    colnames(rank_tbl) <- gsub('this',id,colnames(rank_tbl))
  }
  m_lng_dt <- left_join(lng_dt %>% filter(deadbody == this_db),rank_tbl)
  m_lng_dt <- m_lng_dt %>% rename(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt$Source_class2[m_lng_dt$time_of_lineage == 'early'] <- NA
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  ggtree(tree) %<+% m_lng_dt + theme_tree2()+
    coord_cartesian(xlim = c(-1,35))+
    #geom_text2(aes(subset = is.na(Source_class2) == F, x=x-n_pointmt+5, y = y+0.25, label = Source_class2), size=2.5, color="red")+
    geom_segment2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, xend = x -n_pointmt +1, y = y, yend = y, color = VAF1), size=2)+
    #geom_text2(aes(subset = is.na(VAF1) == F, x=x-n_pointmt, y = y, label = VAF1), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, xend = x -n_pointmt +2, y = y, yend = y, color = VAF2), size=2)+
    #geom_text2(aes(subset = is.na(VAF2) == F, x=x-n_pointmt+1, y = y, label = VAF2), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, xend = x -n_pointmt +3, y = y, yend = y, color = VAF3), size=2)+
    #geom_text2(aes(subset = is.na(VAF3) == F, x=x-n_pointmt+2, y = y, label = VAF3), size=2.5)+
    geom_segment2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, xend = x -n_pointmt +4, y = y, yend = y, color = VAF4), size=2)+
    #geom_text2(aes(subset = is.na(VAF4) == F, x=x-n_pointmt+3, y = y, label = VAF4), size=2.5)+
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
    theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
}




Draw_Heatmap_Diff <- function(fm_vaf_tbl,
                              diff_tbl
                              )
  {
  #set colors
  library(ggsci)
  pal_combi <- c(pal_npg("nrc")(10), pal_d3('category10')(10), pal_aaas('default')(10))
  anatomy_cl3_pal <- pal_combi[c(24,11,2,8,12,17,3,13,18)]
  names(anatomy_cl3_pal) <- c('lt_LE','lt_UE','lt_trunk','rt_LE','rt_UE','rt_trunk','center_trunk','HN','internal_organ')
  lr_pal = pal_combi[c(8, 11, 10, 18)]
  names(lr_pal) <- c("lt","rt","center","unknown")
  cc_pal = pal_combi[c(24,23,18)]
  names(cc_pal) <- c('cranio','caudal','unknown')

  #make mx
  mx <- fm_vaf_tbl %>% filter(time_of_lineage == 'early') %>% select(var_id, sample_id, VAF) %>% spread(key=sample_id, value=VAF) %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
  e_mx <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'early' & is.na(diff_tbl$diff_group) == F ],]
  e_mx <- e_mx[rowSums(is.na(e_mx)) < ncol(e_mx)/2,]
  e_mx <- t(scale(t(e_mx)))
  l_mx <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'late' & is.na(diff_tbl$diff_group) == F ],]
  l_mx <- l_mx[rowSums(is.na(l_mx)) < ncol(l_mx)/2,]
  l_mx <- t(scale(t(l_mx)))
  if("mid" %in% unique(diff_tbl$diff_group)){
    m_mx <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'mid' & is.na(diff_tbl$diff_group) == F ],]
    m_mx <- m_mx[rowSums(is.na(m_mx)) < ncol(m_mx)/2,]
    m_mx <- t(scale(t(m_mx)))
    mx_list = c('e_mx','m_mx','l_mx')
  }else{
    mx_list = c('e_mx','l_mx')
  }
  name_change = c('early','mid','late')
  names(name_change) <- c('e_mx','m_mx','l_mx')
  #draw heatmap
  annot_dt <- fm_vaf_tbl %>% dplyr::rename(left_right = Source_side, cranio_caudal = Source_side_CC, anatomy = anatomy_class3) %>%
    select(sample_id, left_right, cranio_caudal, anatomy) %>% unique() %>% as.data.frame() %>% column_to_rownames('sample_id')
  anno_colors <- list(left_right = lr_pal, cranio_caudal = cc_pal, anatomy = anatomy_cl3_pal)
  library(pheatmap)
  p_list=list()
  n=1
  for (this_mx in mx_list){
    if(n < length(mx_list)){
      p_list[[n]] <- pheatmap(get(this_mx), annotation_col = annot_dt, show_colnames = F, clustering_method ="ward.D", cluster_rows = F, main=name_change[this_mx], annotation_legend = F, annotation_colors = anno_colors)[[4]]
    } else {
      p_list[[n]] <- pheatmap(get(this_mx), annotation_col = annot_dt, show_colnames = F, clustering_method ="ward.D", cluster_rows = F, main=name_change[this_mx], annotation_legend = T, annotation_colors = anno_colors)[[4]]
    }
    n = n+1
  }
  library(cowplot)
  plot_grid(plotlist=p_list, nrow=1, rel_widths = c(rep(1, length(mx_list)-1),1.25))
}


Draw_PCA_Diff <- function(fm_vaf_tbl, diff_tbl){
  mx <- fm_vaf_tbl %>% filter(time_of_lineage == 'early') %>% select(var_id, sample_id, VAF) %>% spread(key=sample_id, value=VAF) %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
  library(missForest) #NA value imputation using missForest
  Early <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'early' & is.na(diff_tbl$diff_group) == F ],]
  Early <- missForest(t(Early))$ximp
  Mid <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'mid' & is.na(diff_tbl$diff_group) == F ],]
  Mid <- missForest(t(Mid))$ximp
  Late <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'late' & is.na(diff_tbl$diff_group) == F ],]
  Late <- missForest(t(Late))$ximp
  p_list=list()
  n=0
  for (timing in c('Early', 'Mid', 'Late')){
    n = n+1
    df_pca <- prcomp(as.data.frame(get(timing)), scale=T, center = T)
    df_out <- as.data.frame(df_pca$x)
    dt <- left_join(df_out %>% as.data.frame() %>% rownames_to_column('sample_id') %>% as.tibble(), pmeta_dt)
    p_list[[n]] <- ggplot(dt, aes(x=PC1, y=PC2, color=anatomy_class3))+
      geom_point(size=3, alpha=0.7)+
      scale_color_manual(values=anatomy_cl3_pal)+
      theme_syp+
      ggtitle(timing)
  }
  p_list[[1]] <- p_list[[1]]+theme(legend.position = 'none')
  p_list[[2]] <- p_list[[2]]+theme(legend.position = 'none')
  library(cowplot)
  plot_grid(plotlist=p_list, nrow=1, rel_widths = c(1,1,1.25))
  }


Draw_tSNE_Diff <- function(
  fm_vaf_tbl,
  diff_tbl,
  color_type = "LR", #c("LR","CC","AC") LR: left-right, CC: cranio-caudal, AC: anatomy_class
  seed = 1 #the input number of set.seed()
  ){
  set.seed(seed)
  
  #color setting
  library(ggsci)
  pal_combi <- c(pal_npg("nrc")(10), pal_d3('category10')(10), pal_aaas('default')(10))
  anatomy_cl3_pal <- pal_combi[c(24,11,2,8,12,17,3,13,18)]
  names(anatomy_cl3_pal) <- c('lt_LE','lt_UE','lt_trunk','rt_LE','rt_UE','rt_trunk','center_trunk','HN','internal_organ')
  lr_pal = pal_combi[c(8, 11, 10, 18)]
  names(lr_pal) <- c("lt","rt","center","unknown")
  cc_pal = pal_combi[c(24,23,18)]
  names(cc_pal) <- c('cranio','caudal','unknown')
  this_pal = c(anatomy_cl3_pal, lr_pal, cc_pal)
  #make mx in each group
  mx <- fm_vaf_tbl %>% filter(time_of_lineage == 'early') %>% select(var_id, sample_id, VAF) %>% spread(key=sample_id, value=VAF) %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
  mx <- mx[rowSums(mx, na.rm = T)!=0,]
  library(missForest) #NA value imputation using missForest
  Early <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'early' & is.na(diff_tbl$diff_group) == F ],]
  Early <- Early[rowSums(is.na(Early)) < ncol(Early)/2,]
  Early <- missForest(t(Early))$ximp
  Early <- Early[rowSums(Early) >0,]
  Early <- t(scale(t(Early)))
  Late <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'late' & is.na(diff_tbl$diff_group) == F ],]
  Late <- Late[rowSums(is.na(Late)) < ncol(Late)/2,]
  Late <- missForest(t(Late))$ximp
  Late <- Late[rowSums(Late) >0,]
  Late <- t(scale(t(Late)))
  if("mid" %in% diff_tbl$diff_group %>% unique){
    Mid <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'mid' & is.na(diff_tbl$diff_group) == F ],]
    Mid <- Mid[rowSums(is.na(Mid)) < ncol(Mid)/2,]
    Mid <- missForest(t(Mid))$ximp
    Mid <- Mid[rowSums(Mid) >0,]
    Mid <- t(scale(t(Mid)))
  }
  
  #run tsne
  library(Rtsne)
  p_list=list()
  n=0
  if("mid" %in% diff_tbl$diff_group %>% unique){
    timing_group = c("Early","Mid","Late")
  } else{
    timing_group = c("Early","Late")
  }
  for (timing in timing_group){
    print(timing)
    n = n+1
    df_tsne <- Rtsne(get(timing), perplexity=10)
    df_out <- as.data.frame(df_tsne$Y)
    rownames(df_out) <- rownames(get(timing))
    colnames(df_out) <- c('tSNE1','tSNE2')
    dt <- left_join(df_out %>% as.data.frame() %>% rownames_to_column('sample_id') %>% as.tibble(), pmeta_dt)
    if(color_type== "LR"){
      dt <- dt %>% mutate(color_type = Source_side)
    } else if (color_type == 'CC'){
      dt <- dt %>% mutate(color_type = Source_side_CC)
    } else if (color_type == 'AC'){
      dt <- dt %>% mutate(color_type = anatomy_class3)
    }
    if(n < length(timing_group)){
      p_list[[n]] <- ggplot(dt, aes(x=tSNE1, y=tSNE2, color=color_type))+
        geom_point(size=5, alpha=0.7)+
        scale_color_manual(values=this_pal)+
        theme_syp+
        theme(legend.position = 'none')+
        ggtitle(timing)  
    } else{
      p_list[[n]] <- ggplot(dt, aes(x=tSNE1, y=tSNE2, color=color_type))+
        geom_point(size=5, alpha=0.7)+
        scale_color_manual(values=this_pal)+
        theme_syp+
        ggtitle(timing)  
    }
  }
  library(cowplot)
  plot_grid(plotlist=p_list, nrow=1, rel_widths = c(rep(1,length(timing_group)-1),1.25))
}

Draw_tSNE_total <- function(
  fm_vaf_tbl,
  color_type = "LR"  #c("LR","CC","AC") LR: left-right, CC: cranio-caudal, AC: anatomy_class
){
  set.seed(1)
  
  #color setting
  library(ggsci)
  pal_combi <- c(pal_npg("nrc")(10), pal_d3('category10')(10), pal_aaas('default')(10))
  anatomy_cl3_pal <- pal_combi[c(24,11,2,8,12,17,3,13,18)]
  names(anatomy_cl3_pal) <- c('lt_LE','lt_UE','lt_trunk','rt_LE','rt_UE','rt_trunk','center_trunk','HN','internal_organ')
  lr_pal = pal_combi[c(8, 11, 10, 18)]
  names(lr_pal) <- c("lt","rt","center","unknown")
  cc_pal = pal_combi[c(24,23,18)]
  names(cc_pal) <- c('cranio','caudal','unknown')
  this_pal = c(anatomy_cl3_pal, lr_pal, cc_pal)
  #make mx in each group
  mx <- fm_vaf_tbl %>% filter(time_of_lineage == 'early') %>% select(var_id, sample_id, VAF) %>% spread(key=sample_id, value=VAF) %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
  mx <- mx[rowSums(mx, na.rm = T)!=0,]
  mx <- mx[rowSums(is.na(mx)) < ncol(mx)/2,]
  library(missForest) #NA value imputation using missForest
  mx <- missForest(t(mx))$ximp
  mx <- mx[rowSums(mx) >0,]
  mx <- t(scale(t(mx)))
  #run tsne
  library(Rtsne)
  df_tsne <- Rtsne(mx, perplexity=10)
  df_out <- as.data.frame(df_tsne$Y)
  rownames(df_out) <- rownames(mx)
  colnames(df_out) <- c('tSNE1','tSNE2')
  dt <- left_join(df_out %>% as.data.frame() %>% rownames_to_column('sample_id') %>% as.tibble(), pmeta_dt)
  if(color_type== "LR"){
    dt <- dt %>% mutate(color_type = Source_side)
  } else if (color_type == 'CC'){
    dt <- dt %>% mutate(color_type = Source_side_CC)
  } else if (color_type == 'AC'){
    dt <- dt %>% mutate(color_type = anatomy_class3)
  }
  g<- ggplot(dt, aes(x=tSNE1, y=tSNE2, color=color_type))+
        geom_point(size=5, alpha=0.7)+
        scale_color_manual(values=this_pal)+
        theme_syp
  print(g)
}




Draw_Silhouette_Coeff_Anatomy_3group <- function(fm_vaf_tbl, diff_tbl){
  mx <- fm_vaf_tbl %>% filter(time_of_lineage == 'early') %>% select(var_id, sample_id, VAF) %>% spread(key=sample_id, value=VAF) %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
  mx <- mx[rowSums(is.na(mx)) < ncol(mx)-3,] 
  library(missForest) #NA value imputation using missForest
  Early <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'early' & is.na(diff_tbl$diff_group) == F ],]
  Early <- missForest(t(Early))$ximp
  Early <- t(scale(t(Early)))
  Mid <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'mid' & is.na(diff_tbl$diff_group) == F ],]
  Mid <- missForest(t(Mid))$ximp
  Mid <- t(scale(t(Mid)))
  Late <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'late' & is.na(diff_tbl$diff_group) == F ],]
  Late <- missForest(t(Late))$ximp
  Late <- t(scale(t(Late)))
  dist_E <- as.matrix(dist(Early))
  dist_M <- as.matrix(dist(Mid))
  dist_L <- as.matrix(dist(Late))
  res_tbl <- tibble(anatomy_class = as.character(), early_ratio = as.numeric(), mid_ratio = as.numeric(), late_ratio = as.numeric())
  n=0
  for (this_class in unique(fm_vaf_tbl$anatomy_class3)){
    print(this_class)
    n = n+1
    class_list <- fm_vaf_tbl %>% filter(deadbody == this_db & anatomy_class3 == this_class) %>% .$sample_id %>% unique()
    E_within_dist  <- dist_E[class_list,class_list]
    E_within_value  <- mean(E_within_dist[E_within_dist>0])
    M_within_dist  <- dist_M[class_list,class_list]
    M_within_value  <- mean(M_within_dist[M_within_dist>0])
    L_within_dist  <- dist_L[class_list,class_list]
    L_within_value  <- mean(L_within_dist[L_within_dist>0])
    E_other_values = c();M_other_values=c();L_other_values=c()
    for (other_class in setdiff(unique(fm_vaf_tbl$anatomy_class3), this_class)){
      other_class_list <- fm_vaf_tbl %>% filter(deadbody == this_db & anatomy_class3 == other_class) %>% .$sample_id %>% unique()
      E_other_values = c(E_other_values, mean(dist_E[class_list,other_class_list]))
      M_other_values = c(M_other_values, mean(dist_M[class_list,other_class_list]))
      L_other_values = c(L_other_values, mean(dist_L[class_list,other_class_list]))
    }
    e_v = 1-E_within_value/min(E_other_values)
    m_v = 1-M_within_value/min(M_other_values)
    l_v = 1-L_within_value/min(L_other_values)
    res_tbl[n,] <- list(this_class, e_v, m_v, l_v)
  }
  m_res_tbl <- res_tbl %>% gather(-anatomy_class, key='timing',value='Silhouette Coeff') %>% mutate(time = gsub('_ratio','',timing)) 
  m_res_tbl$time <- factor(m_res_tbl$time, levels = c('early','mid','late'))
  ggplot(m_res_tbl, aes(x=time, y=`Silhouette Coeff`, color=anatomy_class))+
    geom_point()+
    geom_line(aes(group = anatomy_class))+
    scale_color_manual(values = anatomy_cl3_pal)+
    theme_syp
}


Draw_Silhouette_Coeff_Anatomy_2group <- function(fm_vaf_tbl, diff_tbl){
  mx <- fm_vaf_tbl %>% filter(time_of_lineage == 'early') %>% select(var_id, sample_id, VAF) %>% spread(key=sample_id, value=VAF) %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
  library(missForest) #NA value imputation using missForest
  Early <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'early' & is.na(diff_tbl$diff_group) == F ],]
  Early <- missForest(t(Early))$ximp
  Early <- t(scale(t(Early)))
  Late <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'late' & is.na(diff_tbl$diff_group) == F ],]
  Late <- missForest(t(Late))$ximp
  Late <- t(scale(t(Late)))
  dist_E <- as.matrix(dist(Early))
  dist_L <- as.matrix(dist(Late))
  res_tbl <- tibble(anatomy_class = as.character(), early_ratio = as.numeric(), late_ratio = as.numeric())
  n=0
  for (this_class in unique(fm_vaf_tbl$anatomy_class3)){
    print(this_class)
    n = n+1
    class_list <- fm_vaf_tbl %>% filter(deadbody == this_db & anatomy_class3 == this_class) %>% .$sample_id %>% unique()
    E_within_dist  <- dist_E[class_list,class_list]
    E_within_value  <- mean(E_within_dist[E_within_dist>0])
    L_within_dist  <- dist_L[class_list,class_list]
    L_within_value  <- mean(L_within_dist[L_within_dist>0])
    E_other_values = c();L_other_values=c()
    for (other_class in setdiff(unique(fm_vaf_tbl$anatomy_class3), this_class)){
      other_class_list <- fm_vaf_tbl %>% filter(deadbody == this_db & anatomy_class3 == other_class) %>% .$sample_id %>% unique()
      E_other_values = c(E_other_values, mean(dist_E[class_list,other_class_list]))
      L_other_values = c(L_other_values, mean(dist_L[class_list,other_class_list]))
    }
    e_v = 1-E_within_value/min(E_other_values)
    l_v = 1-L_within_value/min(L_other_values)
    res_tbl[n,] <- list(this_class, e_v, l_v)
  }
  m_res_tbl <- res_tbl %>% gather(-anatomy_class, key='timing',value='Silhouette Coeff') %>% mutate(time = gsub('_ratio','',timing)) 
  m_res_tbl$time <- factor(m_res_tbl$time, levels = c('early','late'))
  ggplot(m_res_tbl, aes(x=time, y=`Silhouette Coeff`, color=anatomy_class))+
    geom_point()+
    geom_line(aes(group = anatomy_class))+
    scale_color_manual(values = anatomy_cl3_pal)+
    theme_syp
}


Draw_Distance_AllOther <- function(fm_vaf_tbl, diff_tbl){
  mx <- fm_vaf_tbl %>% filter(time_of_lineage == 'early') %>% select(var_id, sample_id, VAF) %>% spread(key=sample_id, value=VAF) %>% as.data.frame() %>% column_to_rownames('var_id') %>% as.matrix()
  library(missForest) #NA value imputation using missForest
  Early <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'early' & is.na(diff_tbl$diff_group) == F ],]
  Early <- missForest(t(Early))$ximp
  Early <- t(scale(t(Early)))
  Mid <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'mid' & is.na(diff_tbl$diff_group) == F ],]
  Mid <- missForest(t(Mid))$ximp
  Mid <- t(scale(t(Mid)))
  Late <- mx[rownames(mx) %in% diff_tbl$var_id[diff_tbl$diff_group == 'late' & is.na(diff_tbl$diff_group) == F ],]
  Late <- missForest(t(Late))$ximp
  Late <- t(scale(t(Late)))
  dist_E <- as.matrix(dist(Early))
  dist_M <- as.matrix(dist(Mid))
  dist_L <- as.matrix(dist(Late))
  res_tbl <- tibble(anatomy_class = as.character(), early_ratio = as.numeric(), mid_ratio = as.numeric(), late_ratio = as.numeric())
  n=0
  for (this_class in unique(fm_vaf_tbl$anatomy_class3)){
    print(this_class)
    n = n+1
    class_list <- fm_vaf_tbl %>% filter(deadbody == this_db & anatomy_class3 == this_class) %>% .$sample_id %>% unique()
    E_within_dist  <- dist_E[class_list,class_list]
    E_within_value  <- mean(E_within_dist[E_within_dist>0])
    M_within_dist  <- dist_M[class_list,class_list]
    M_within_value  <- mean(M_within_dist[M_within_dist>0])
    L_within_dist  <- dist_L[class_list,class_list]
    L_within_value  <- mean(L_within_dist[L_within_dist>0])
    E_other_values = c();M_other_values=c();L_other_values=c()
    other_class_list <- fm_vaf_tbl %>% filter(deadbody == this_db & anatomy_class3 != this_class) %>% .$sample_id %>% unique()
    E_other_value = mean(dist_E[class_list,other_class_list])
    M_other_value = mean(dist_M[class_list,other_class_list])
    L_other_value = mean(dist_L[class_list,other_class_list])
    e_v = 1-E_within_value/E_other_value
    m_v = 1-M_within_value/M_other_value
    l_v = 1-L_within_value/L_other_value
    res_tbl[n,] <- list(this_class, e_v, m_v, l_v)
  }
  m_res_tbl <- res_tbl %>% gather(-anatomy_class, key='timing',value='Clustering Score') %>% mutate(time = gsub('_ratio','',timing)) 
  m_res_tbl$time <- factor(m_res_tbl$time, levels = c('early','mid','late'))
  ggplot(m_res_tbl, aes(x=time, y=`Clustering Score`, color=anatomy_class))+
    geom_point()+
    geom_line(aes(group = anatomy_class))+
    scale_color_manual(values = anatomy_cl3_pal)+
    theme_syp
}


VariantGroupingByTissueVAF <- function(fm_vaf_tbl, lid2, bootstrapR){
  ffm_vaf_tbl <- fm_vaf_tbl %>% filter(lineage_id2 == lid2)
  med_tbl <- ffm_vaf_tbl %>% group_by(var_id) %>% summarise(medVAF = median(VAF)) %>% arrange(desc(medVAF)) 
  var_order <- med_tbl %>% pull(var_id)
  g1 <- ggplot(ffm_vaf_tbl, aes(x=var_id, y=VAF))+
    geom_boxplot(outlier.size=-1)+
    geom_point(alpha=0.1)+
    geom_line(aes(group=sample_id), alpha=0.1)+
    geom_line(data=med_tbl, aes(x=var_id, y=medVAF, group='1'), color='red')+
    scale_x_discrete(limits = var_order)+theme_syp+theme(axis.text.x=element_blank(), axis.title.x = element_blank())+
    ggtitle(lid2)
  
  mycomparison <- combn(unique(ffm_vaf_tbl$var_id),2, simplify = F)
  res_tbl <- tibble(var1= character(), var2= character(), pval = as.numeric())
  for (i in 1:length(mycomparison)){
    tmp_tbl <- ffm_vaf_tbl %>% filter(var_id %in% mycomparison[[i]]) %>% select(sample_id, var_id, VAF) %>%  spread(key=var_id, value=VAF) %>% filter(rowSums(is.na(.))==0)
    if (nrow(tmp_tbl) < 5){
      res_tbl[i,] <- append(mycomparison[[i]],"NA")
    }else{
      res <- t.test(as.numeric(as.matrix(tmp_tbl[2])), as.numeric(as.matrix(tmp_tbl[3])), paired=T)
      res_tbl[i,] <- append(mycomparison[[i]],res$p.value)
    }
  }
  res_tbl$pval = as.numeric(res_tbl$pval)
  res_tbl <- res_tbl %>% filter(is.na(pval)== F) %>% filter(pval >= 0.05)
   if (nrow(res_tbl) == 0){
    g1_2 <- ggplot(res_tbl)
  } else{
    res_tbl$index <- 1:nrow(res_tbl)
    ct_tbl <- fm_vaf_tbl %>% filter(lineage_id2 == lid2) %>% group_by(var_id) %>% count()
    g1_2 <- ggplot(res_tbl)+
      geom_point(aes(x=var1, y=index))+
      geom_point(aes(x=var2, y=index))+
      geom_text(data=ct_tbl, aes(x=var_id, y=0, label = n), color="blue")+
      geom_segment(aes(x = var1, xend = var2, y=index, yend=index), color="red")+
      scale_x_discrete(limits=var_order)+
      theme_syp+theme(axis.text.x= element_blank(), axis.title = element_blank())
    }
  
  #bootstrap
  col_v <- paste0('VAF',1:length(var_order))
  names(col_v) <- var_order
  bt_dt <- ffm_vaf_tbl %>% select(sample_id, var_id, VAF) %>% spread(key=var_id, value=VAF)
  library(boot)
  colnames(bt_dt) <- c('sample_id',col_v[colnames(bt_dt)[2:ncol(bt_dt)]])
  bt_mx <- bt_dt %>% as.data.frame() %>% column_to_rownames('sample_id') %>% as.matrix()
  make_median_list <- function(data, indices){
    d <- data[indices,]
    res <- apply(d, 2, function(x) median(x[is.na(x)==F]))
    return(res)
  }
  bt_res <- boot(data=bt_mx, statistic=make_median_list, R=bootstrapR)
  res_mx <-bt_res$t
  colnames(res_mx) <- colnames(bt_mx)
  rownames(res_mx) <- 1:nrow(res_mx)
  res_dt <- res_mx %>% as.data.frame() %>% rownames_to_column('index') %>% as.tibble()
  m_res_dt <- res_dt%>% gather(-index, key=var_order, value=medianVAF)
  g2 <- ggplot(m_res_dt, aes(x=var_order, y=medianVAF))+
    geom_boxplot()+
    geom_line(aes(group=index), alpha=0.005)+
    scale_x_discrete(limits = col_v)+theme_syp+theme(axis.text.x=element_blank(), axis.title.x = element_blank())
  res_tbl <- tibble(var1 = as.character(), var2=as.character(), var1_hi = as.integer(), var2_hi=as.integer())
  var_list <- paste0('VAF', 1:length(var_order))
  n=0
  for (var1 in var_list){
    for(var2 in setdiff(var_list, var1)){
      n=n+1
      var1_hi <- res_dt %>% filter(get(var1) > get(var2)) %>% nrow()
      var2_hi <- res_dt %>% filter(get(var1) < get(var2)) %>% nrow()
      res_tbl[n,] <- list(var1, var2, var1_hi, var2_hi)
    }
  }
  f_res_tbl <- res_tbl %>% filter(var1_hi < bootstrapR*0.95 & var2_hi < bootstrapR*0.95)
  g_list=list()
  n=0
  for(this_var in var_list){
    n=n+1
    g_list[[n]] <- unique(c(this_var, f_res_tbl %>% filter(var1 == this_var) %>% pull(var2), f_res_tbl %>% filter(var2 == this_var) %>% pull(var1)))
  }
  g_hm_dt <- as.tibble(data.frame(matrix(vector(), length(g_list), length(var_list), dimnames=list(1:length(g_list), var_list))))
  for (i in 1:length(g_list)){
    g_hm_dt[i,][match(g_list[[i]], var_list)] <- i
  }
  g_hm_dt <- g_hm_dt %>% gather(key=vaf_order, value=group) %>% filter(is.na(group) == F)
  g3 <- ggplot(g_hm_dt, aes(x=vaf_order, y=group, fill=group))+
    geom_tile()+
    scale_x_discrete(limits=col_v, labels=names(col_v))+
    theme_syp+theme(axis.text.x=element_text(angle=45, hjust=1), axis.title.x = element_blank(), legend.position = 'none')
  library(cowplot)
  plot_grid(g1,g1_2,g2,g3,ncol=1, align='v')
}

