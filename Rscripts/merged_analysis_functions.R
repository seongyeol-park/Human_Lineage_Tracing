#2021-02-27 lineage_id_cindel -> lineage_id


Draw_Phylo_Full_DB <- function(this_db){
  library(cowplot)
  if (this_db %in% c('DB2','DB5','DB8','DB9','DB10')){
    print(this_db)
    tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
    m_lng_dt <- lng_dt %>% filter(deadbody == this_db) %>% rename(taxa = lineage_id) %>% select(-deadbody) %>% mutate(lid = taxa)
    bin <- ifelse(max(m_lng_dt$end_mtn) > 30000, 10000, ifelse(max(m_lng_dt$end_mtn) > 15000, 5000, 2000))
    p1 <- ggtree(tree) %<+% m_lng_dt + theme_tree2() +
      geom_text2(aes(x=x-n_pointmt, y = y, label = lid), hjust=0,color="red",size=2.5)+
      coord_cartesian(xlim = c(-1,35))+
      ggtitle(this_db)+
      theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
    p2 <- ggtree(tree) %<+% m_lng_dt + geom_tree()+theme_tree2() +
      geom_tiplab(aes(label = paste0(samples, ' (', Source3, ')')), hjust=0)+
      coord_cartesian(xlim = c(max(m_lng_dt$end_mtn)/13, max(m_lng_dt$end_mtn)*1.5))+
      scale_x_continuous(breaks=seq(bin, max(m_lng_dt$end_mtn)*1.5, bin))+
      theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0),"cm"))
    plot_grid(p1,p2, nrow=1, align="h", rel_widths = c(1,3)) %>% print()
  } else if (this_db == 'DB3'){
    tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
    m_lng_dt <- lng_dt %>% filter(deadbody == this_db) %>% rename(taxa = lineage_id) %>% select(-deadbody)%>% mutate(lid = taxa)
    sub_lng_dt <- m_lng_dt %>% filter(end_mtn > 20000)
    db3p1 <- ggtree(tree) %<+% m_lng_dt +theme_tree2() +
      coord_cartesian(xlim = c(-1,35))+
      geom_text2(aes(x=x-n_pointmt, y = y, label = lid), hjust=0,color="red",size=2.5)+
      ggtitle(this_db)+
      theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
    db3p2 <- ggtree(tree) %<+% m_lng_dt + geom_tree()+theme_tree2() +
      coord_cartesian(xlim = c(1000, 15000))+
      geom_tiplab(aes(label = paste0(samples, ' (', Source3, ')')), hjust=0)+
      scale_x_continuous(breaks=seq(2000, 15000,2000))+
      theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0),"cm"))
    db3p3 <- ggtree(tree) %<+% sub_lng_dt + geom_tree()+theme_tree2() +
      coord_cartesian(xlim = c(30000, 120000))+
      geom_tiplab(aes(label = paste0(samples, ' (', Source3, ')')), hjust=0)+
      scale_x_continuous(breaks=seq(20000, 100000,20000), labels = function(x) format(x, scientific = F))+
      theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0),"cm"))
    plot_grid(db3p1,db3p2,db3p3, nrow=1, align="h", rel_widths = c(1,1.5,1.5))
  } else if (this_db == 'DB6'){
    tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
    m_lng_dt <- lng_dt %>% filter(deadbody == this_db) %>% rename(taxa = lineage_id) %>% select(-deadbody)%>% mutate(lid = taxa)
    bin <- ifelse(max(m_lng_dt$end_mtn) > 30000, 10000, ifelse(max(m_lng_dt$end_mtn) > 15000, 5000, 2000))
    p1 <- ggtree(tree) %<+% m_lng_dt +theme_tree2() +
      geom_text2(aes(x=x-n_pointmt, y = y, label = lid), hjust=0,color="red",size=2.5)+
      coord_cartesian(xlim = c(-1,35))+
      ggtitle(this_db)+
      theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
    p2 <- ggtree(tree) %<+% m_lng_dt + geom_tree()+theme_tree2() +
      geom_tiplab(aes(label = paste0(samples, ' (', Source3, ')')), size=3, align=F, linesize=0.1, hjust=-0.3)+
      coord_cartesian(xlim = c(max(m_lng_dt$end_mtn)/13, max(m_lng_dt$end_mtn)*1.5))+
      scale_x_continuous(breaks=seq(bin, max(m_lng_dt$end_mtn)*1.5, bin))+
      theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0),"cm"))
    plot_grid(p1,p2, nrow=1, align="h", rel_widths = c(1.5,3)) 
  }
}

Draw_Phylo_Full_DB_markSample <- function(this_db, t_sid){
  library(cowplot)
  t_age <- meta_dt$age[meta_dt$deadbody == this_db] %>% unique()
  t_gender <- meta_dt$gender[meta_dt$deadbody == this_db] %>% unique()
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  m_lng_dt <- lng_dt %>% filter(deadbody == this_db) %>% rename(taxa = lineage_id) %>% select(-deadbody) %>% mutate(lid = taxa)
  t_lid <- m_lng_dt$lineage_id2[m_lng_dt$samples == t_sid]
  t_lid_list <- unlist(strsplit(t_lid,'\\-'))
  t_lid_line = c()
  for (i in 1:length(t_lid_list)){
    t_lid_line = c(t_lid_line, paste(t_lid_list[1:i], collapse='-'))
  }
  m_lng_dt <- m_lng_dt %>% mutate(group = ifelse(lineage_id2 %in% t_lid_line, 'Y', 'N'))
  if (this_db %in% c('DB2','DB5','DB8','DB9','DB10')){
    print(this_db)
    
    bin <- ifelse(max(m_lng_dt$end_mtn) > 30000, 10000, ifelse(max(m_lng_dt$end_mtn) > 15000, 5000, 2000))
    p1 <- ggtree(tree, aes(color = group)) %<+% m_lng_dt + theme_tree2() +
      #geom_text2(aes(x=x-n_pointmt, y = y, label = lid), hjust=0,color="red",size=2.5)+
      coord_cartesian(xlim = c(-1,35))+
      ggtitle(paste0(this_db, ' (',t_gender,'/',t_age,')'))+
      scale_color_manual(values = c('Y'='red','N'='black'))+
      theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
    p2 <- ggtree(tree,aes(color = group) ) %<+% m_lng_dt + theme_tree2() +
      #geom_tiplab(aes(label = paste0(samples, ' (', Source3, ')')), hjust=0)+
      coord_cartesian(xlim = c(max(m_lng_dt$end_mtn)/13, max(m_lng_dt$end_mtn)*1.5))+
      scale_x_continuous(breaks=seq(bin, max(m_lng_dt$end_mtn)*1.5, bin))+
      scale_color_manual(values = c('Y'='red','N'='black'))+
      theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0),"cm"))
    plot_grid(p1,p2, nrow=1, align="h", rel_widths = c(1,3)) %>% print()
  } else if (this_db == 'DB3'){
    sub_lng_dt <- m_lng_dt %>% filter(end_mtn > 20000)
    db3p1 <-  ggtree(tree, aes(color = group))  %<+% m_lng_dt +theme_tree2() +
      coord_cartesian(xlim = c(-1,35))+
      #geom_text2(aes(x=x-n_pointmt, y = y, label = lid), hjust=0,color="red",size=2.5)+
      scale_color_manual(values = c('Y'='red','N'='black'))+
      ggtitle(this_db)+
      theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
    db3p2 <-  ggtree(tree, aes(color = group))  %<+% m_lng_dt + geom_tree()+theme_tree2() +
      coord_cartesian(xlim = c(1000, 15000))+
      #geom_tiplab(aes(label = paste0(samples, ' (', Source3, ')')), hjust=0)+
      scale_color_manual(values = c('Y'='red','N'='black'))+
      scale_x_continuous(breaks=seq(2000, 15000,2000))+
      theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0),"cm"))
    db3p3 <-  ggtree(tree, aes(color = group))  %<+% sub_lng_dt + geom_tree()+theme_tree2() +
      coord_cartesian(xlim = c(30000, 120000))+
      #geom_tiplab(aes(label = paste0(samples, ' (', Source3, ')')), hjust=0)+
      scale_color_manual(values = c('Y'='red','N'='black'))+
      scale_x_continuous(breaks=seq(20000, 100000,20000), labels = function(x) format(x, scientific = F))+
      theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0),"cm"))
    plot_grid(db3p1,db3p2,db3p3, nrow=1, align="h", rel_widths = c(1,1.5,1.5)) %>% print()
  } else if (this_db == 'DB6'){
    bin <- ifelse(max(m_lng_dt$end_mtn) > 30000, 10000, ifelse(max(m_lng_dt$end_mtn) > 15000, 5000, 2000))
    p1 <-  ggtree(tree, aes(color = group))  %<+% m_lng_dt +theme_tree2() +
      #geom_text2(aes(x=x-n_pointmt, y = y, label = lid), hjust=0,color="red",size=2.5)+
      coord_cartesian(xlim = c(-1,35))+
      ggtitle(this_db)+
      scale_color_manual(values = c('Y'='red','N'='black'))+
      theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
    p2 <-  ggtree(tree, aes(color = group))  %<+% m_lng_dt + geom_tree()+theme_tree2() +
      #geom_tiplab(aes(label = paste0(samples, ' (', Source3, ')')), size=3, align=F, linesize=0.1, hjust=-0.3)+
      coord_cartesian(xlim = c(max(m_lng_dt$end_mtn)/13, max(m_lng_dt$end_mtn)*1.5))+
      scale_x_continuous(breaks=seq(bin, max(m_lng_dt$end_mtn)*1.5, bin))+
      scale_color_manual(values = c('Y'='red','N'='black'))+
      theme(axis.text = element_text(size=15, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0),"cm"))
    plot_grid(p1,p2, nrow=1, align="h", rel_widths = c(1.5,3)) %>% print()
  }
}



Draw_totalTree_Overview <- function(meta_dt, lng_dt, nwk_tbl){
  library(cowplot)
  plist <- list()
  n=0
  for (db in c('DB2','DB5','DB6','DB8','DB9','DB10')){
    print(db)
    n = n+1
    tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
    m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% mutate(l_id = taxa) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
    m_lng_dt <- lng_dt %>% filter(deadbody == db & time_of_lineage != 'germline') %>% rename(taxa = lineage_id) %>% select(-deadbody)
    m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
    bin <- ifelse(max(m_lng_dt$end_mtn) > 30000, 10000, ifelse(max(m_lng_dt$end_mtn) > 15000, 5000, 2000))
    p1 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = time_of_lineage))+theme_tree2() +
      coord_cartesian(xlim = c(-1,35))+
      scale_color_manual(values = timing_pal)+
      ggtitle(db)+
      theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
    p2 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = time_of_lineage))+theme_tree2() +
      #geom_tiplab(aes(label = paste0(samples, ' (', Source3, ')')), hjust=0)+
      coord_cartesian(xlim = c(max(m_lng_dt$end_mtn)/14, max(m_lng_dt$end_mtn)))+
      scale_color_manual(values = timing_pal)+
      scale_x_continuous(breaks=seq(bin, max(m_lng_dt$end_mtn)+500, bin))+
      theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0),"cm"))
    plist[[n]] <- plot_grid(p1,p2, nrow=1, align="h", rel_widths = c(1,2))
  }
  
  db='DB3'
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
  m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% mutate(l_id = taxa) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
  m_lng_dt <- lng_dt %>% filter(deadbody == db& time_of_lineage != 'germline') %>% rename(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
  db3p1 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = time_of_lineage))+theme_tree2() +
    coord_cartesian(xlim = c(-1,35))+
    scale_color_manual(values = timing_pal)+
    ggtitle(db)+
    theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  db3p2 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = time_of_lineage))+theme_tree2() +
    coord_cartesian(xlim = c(500, 7000))+
    scale_color_manual(values = timing_pal)+
    scale_x_continuous(breaks=seq(2000, 8000,2000))+
    theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0),"cm"))
  db3p3 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = time_of_lineage))+theme_tree2() +
    coord_cartesian(xlim = c(10000, max(m_lng_dt$end_mtn)))+
    scale_color_manual(values = timing_pal)+
    scale_x_continuous(breaks=seq(10000, 60000,20000))+
    theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0),"cm"))
  pDB3 <- plot_grid(db3p1,db3p2,db3p3, nrow=1, align="h")
  lplot <- plot_grid(plotlist = list(plist[[1]], pDB3, plist[[2]], plist[[4]], plist[[5]], plist[[6]]), nrow=2)
  plot_grid(plotlist = list(lplot, plist[[3]]), nrow=1, rel_widths = c(1.5,1))
}

Draw_totalEarly_Tree_BranchN <- function(meta_dt, lng_dt, nwk_tbl){
  library(cowplot)
  plist <- list()
  n=0
  for (db in c('DB2','DB5','DB6','DB8','DB9','DB10')){
    print(db)
    n = n+1
    tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
    m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
    m_lng_dt <- lng_dt %>% filter(deadbody == db) %>% mutate(l_id = lineage_id) %>% rename(taxa = lineage_id) %>% select(-deadbody)
    m_lng_dt <- m_lng_dt %>% mutate(multifurcating = ifelse(branch_n >2, branch_n,NA))
    m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
    bin <- ifelse(max(m_lng_dt$end_mtn) > 30000, 10000, ifelse(max(m_lng_dt$end_mtn) > 15000, 5000, 2000))
    plist[[n]] <- ggtree(tree) %<+% m_lng_dt + geom_tree()+theme_tree2() +
      geom_label(aes(label = multifurcating), hjust=1, alpha=0.5, fill="yellow")+
      coord_cartesian(xlim = c(-5,35))+
      scale_color_manual(values=c('Y'='red','N'='black'))+
      ggtitle(db)+
      theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  }
  db='DB3'
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
  m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% mutate(l_id = taxa) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
  m_lng_dt <- lng_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt <- m_lng_dt %>% mutate(multifurcating = ifelse(branch_n >2, branch_n,NA))
  m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
  pDB3 <- ggtree(tree) %<+% m_lng_dt + geom_tree()+theme_tree2() +
    geom_label(aes(label = multifurcating), hjust=1, alpha=0.5, fill="yellow")+
    coord_cartesian(xlim = c(-5,35))+
    ggtitle(db)+
    theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  lplot <- plot_grid(plotlist = list(plist[[1]], pDB3, plist[[2]], plist[[4]], plist[[5]], plist[[6]]), nrow=2)
  plot_grid(plotlist = list(lplot, plist[[3]]), nrow=1, rel_widths = c(1.5,1))
  
}

Draw_totalTree_CtoTpGprop <- function(meta_dt, lng_dt, nwk_tbl){
  library(cowplot)
  plist <- list()
  n=0
  for (db in c('DB2','DB5','DB6','DB8','DB9','DB10')){
    print(db)
    n = n+1
    tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
    m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% mutate(l_id = taxa) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
    m_lng_dt <- lng_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id) %>% select(-deadbody)
    m_lng_dt <- m_lng_dt %>% mutate(CTGprop = ifelse((n_subs >= 300 & sig7_prop >= 5) | n_CtoT == 0, NA, round(n_CtoTG/n_CtoT,2)))
    m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
    bin <- ifelse(max(m_lng_dt$end_mtn) > 30000, 10000, ifelse(max(m_lng_dt$end_mtn) > 15000, 5000, 2000))
    p1 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = CTGprop))+theme_tree2() +
      coord_cartesian(xlim = c(-1,35))+
      ggtitle(db)+
      scale_color_gradient(low ="blue",high = "yellow", breaks = c(0,1))+
      theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
    p2 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = CTGprop))+theme_tree2() +
      #geom_tiplab(aes(label = paste0(samples, ' (', Source3, ')')), hjust=0)+
      coord_cartesian(xlim = c(max(m_lng_dt$end_mtn)/14, max(m_lng_dt$end_mtn)))+
      scale_x_continuous(breaks=seq(bin, max(m_lng_dt$end_mtn)+500, bin))+
      scale_color_gradient(low ="blue",high = "yellow", breaks = c(0,1))+
      theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0),"cm"))
    plist[[n]] <- plot_grid(p1,p2, nrow=1, align="h", rel_widths = c(1,2))
  }
  db='DB3'
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
  m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% mutate(l_id = taxa) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
  m_lng_dt <- lng_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt <- m_lng_dt %>% mutate(CTGprop = ifelse((n_subs >= 300 & sig7_prop >= 5) | n_CtoT == 0, NA, round(n_CtoTG/n_CtoT,2)))
  m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
  db3p1 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = CTGprop))+theme_tree2() +
    coord_cartesian(xlim = c(-1,35))+
    ggtitle(db)+
    scale_color_gradient(low ="blue",high = "yellow", breaks = c(0,1))+
    theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  db3p2 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = CTGprop))+theme_tree2() +
    coord_cartesian(xlim = c(500, 7000))+
    scale_color_gradient(low ="blue",high = "yellow", breaks = c(0,1))+
    scale_x_continuous(breaks=seq(2000, 8000,2000))+
    theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0),"cm"))
  db3p3 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = CTGprop))+theme_tree2() +
    coord_cartesian(xlim = c(10000, max(m_lng_dt$end_mtn)))+
    scale_x_continuous(breaks=seq(10000, 60000,20000))+
    scale_color_gradient(low ="blue",high = "yellow", breaks = c(0,1))+
    theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0),"cm"))
  pDB3 <- plot_grid(db3p1,db3p2,db3p3, nrow=1, align="h")
  lplot <- plot_grid(plotlist = list(plist[[1]], pDB3, plist[[2]], plist[[4]], plist[[5]], plist[[6]]), nrow=2)
  plot_grid(plotlist = list(lplot, plist[[3]]), nrow=1, rel_widths = c(1.5,1))
}


Draw_totalTree_CNV <- function(meta_dt, lng_dt, nwk_tbl){
  library(cowplot)
  plist <- list()
  n=0
  for (db in c('DB2','DB5','DB6','DB8','DB9','DB10')){
    print(db)
    n = n+1
    tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
    m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% mutate(l_id = taxa) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
    m_lng_dt <- lng_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id) %>% select(-deadbody)
    m_lng_dt$CNV_lineage[m_lng_dt$CNV_lineage =='none'] <- NA
    m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
    bin <- ifelse(max(m_lng_dt$end_mtn) > 30000, 10000, ifelse(max(m_lng_dt$end_mtn) > 15000, 5000, 2000))
    p1 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = CNV_lineage))+theme_tree2() +
      geom_label(aes(label=CNV_lineage, fill = CNV_lineage, color=CNV_lineage),hjust=0, alpha=0.2, size=3)+
      coord_cartesian(xlim = c(-1,35))+
      ggtitle(db)+
      theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
    p2 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = CNV_lineage))+theme_tree2() +
      #geom_tiplab(aes(label = paste0(samples, ' (', Source3, ')')), hjust=-0.5, size=2)+
      geom_label(aes(label=CNV_lineage, fill = CNV_lineage, color=CNV_lineage),hjust=0, alpha=0.2, size=3)+
      coord_cartesian(xlim = c(max(m_lng_dt$end_mtn)/14, max(m_lng_dt$end_mtn)*1.3))+
      scale_x_continuous(breaks=seq(bin, max(m_lng_dt$end_mtn)+500, bin))+
      theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0),"cm"))
    plist[[n]] <- plot_grid(p1,p2, nrow=1, align="h", rel_widths = c(1,2))
  }
  db='DB3'
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
  m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% mutate(l_id = taxa) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
  m_lng_dt <- lng_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt$CNV_lineage[m_lng_dt$CNV_lineage =='none'] <- NA
  m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
  db3p1 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = CNV_lineage))+theme_tree2() +
    geom_label(aes(label=CNV_lineage, fill = CNV_lineage, color=CNV_lineage),hjust=0, alpha=0.2, size=3)+
    coord_cartesian(xlim = c(-1,35))+
    ggtitle(db)+
    theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  db3p2 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = CNV_lineage))+theme_tree2() +
    geom_label(aes(label=CNV_lineage, fill = CNV_lineage, color=CNV_lineage),hjust=0, alpha=0.2, size=3)+
    coord_cartesian(xlim = c(500, 7000))+
    scale_x_continuous(breaks=seq(2000, 8000,2000))+
    theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0),"cm"))
  m_lng_dt3 <- m_lng_dt; m_lng_dt3$CNV_lineage[m_lng_dt3$end_mtn < 10000] <- NA
  db3p3 <- ggtree(tree) %<+% m_lng_dt3 + geom_tree(aes(color = CNV_lineage))+theme_tree2() +
    geom_label(aes(label=CNV_lineage, fill = CNV_lineage, color=CNV_lineage),hjust=0, alpha=0.2, size=3)+
    coord_cartesian(xlim = c(15000, max(m_lng_dt$end_mtn)*2))+
    scale_x_continuous(breaks=seq(10000, 60000,20000))+
    theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0),"cm"))
  pDB3 <- plot_grid(db3p1,db3p2,db3p3, nrow=1, align="h", rel_widths = c(1, 1.5, 1.5))
  lplot <- plot_grid(plotlist = list(plist[[1]], pDB3, plist[[2]], plist[[4]], plist[[5]], plist[[6]]), nrow=2)
  plot_grid(plotlist = list(lplot, plist[[3]]), nrow=1, rel_widths = c(1.5,1))
}

Draw_totalTree_SV <- function(meta_dt, lng_dt, nwk_tbl){
  library(cowplot)
  plist <- list()
  n=0
  for (db in c('DB2','DB5','DB6','DB8','DB9','DB10')){
    print(db)
    n = n+1
    tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
    m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% mutate(l_id = taxa) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
    m_lng_dt <- lng_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id) %>% select(-deadbody)
    m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
    m_lng_dt$n_SV[m_lng_dt$n_SV == 0] <- NA
    bin <- ifelse(max(m_lng_dt$end_mtn) > 30000, 10000, ifelse(max(m_lng_dt$end_mtn) > 15000, 5000, 2000))
    p1 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = n_SV))+theme_tree2() +
      geom_label(aes(label=n_SV, fill = n_SV, color=n_SV),hjust=-0.2, alpha=0.2, size=3)+
      coord_cartesian(xlim = c(-1,35))+
      ggtitle(db)+
      theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
    p2 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = n_SV))+theme_tree2() +
      #geom_tiplab(aes(label = paste0(samples, ' (', Source3, ')')), hjust=-0.5, size=2)+
      geom_label(aes(label=n_SV, fill = n_SV, color=n_SV),hjust=-0.2, alpha=0.2, size=3)+
      coord_cartesian(xlim = c(max(m_lng_dt$end_mtn)/14, max(m_lng_dt$end_mtn)*1.3))+
      scale_x_continuous(breaks=seq(bin, max(m_lng_dt$end_mtn)+500, bin))+
      theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0),"cm"))
    plist[[n]] <- plot_grid(p1,p2, nrow=1, align="h", rel_widths = c(1,2))
  }
  db='DB3'
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
  m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% mutate(l_id = taxa) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
  m_lng_dt <- lng_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
  m_lng_dt$n_SV[m_lng_dt$n_SV == 0] <- NA
  db3p1 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = n_SV))+theme_tree2() +
    geom_label(aes(label=n_SV, fill = n_SV, color=n_SV),hjust=-0.2, alpha=0.2, size=3)+
    coord_cartesian(xlim = c(-1,35))+
    ggtitle(db)+
    theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  db3p2 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = n_SV))+theme_tree2() +
    geom_label(aes(label=n_SV, fill = n_SV, color=n_SV),hjust=-0.2, alpha=0.2, size=3)+
    coord_cartesian(xlim = c(500, 7000))+
    scale_x_continuous(breaks=seq(2000, 8000,2000))+
    theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0),"cm"))
  m_lng_dt3 <- m_lng_dt; m_lng_dt3$n_SV[m_lng_dt3$end_mtn < 10000] <- NA
  db3p3 <- ggtree(tree) %<+% m_lng_dt3 + geom_tree(aes(color = n_SV))+theme_tree2() +
    geom_label(aes(label=n_SV, fill = n_SV, color=n_SV),hjust=-0.2, alpha=0.2, size=3)+
    coord_cartesian(xlim = c(15000, max(m_lng_dt$end_mtn)*2))+
    scale_x_continuous(breaks=seq(10000, 60000,20000))+
    theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0),"cm"))
  pDB3 <- plot_grid(db3p1,db3p2,db3p3, nrow=1, align="h", rel_widths = c(1, 1.5, 1.5))
  lplot <- plot_grid(plotlist = list(plist[[1]], pDB3, plist[[2]], plist[[4]], plist[[5]], plist[[6]]), nrow=2)
  plot_grid(plotlist = list(lplot, plist[[3]]), nrow=1, rel_widths = c(1.5,1))
}

Draw_totalTree_Sig7prop <- function(meta_dt, lng_dt, nwk_tbl){
  library(cowplot)
  plist <- list()
  n=0
  for (db in c('DB2','DB5','DB6','DB8','DB9','DB10')){
    print(db)
    n = n+1
    tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
    m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% mutate(l_id = taxa) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
    m_lng_dt <- lng_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id) %>% select(-deadbody)
    m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
    bin <- ifelse(max(m_lng_dt$end_mtn) > 30000, 10000, ifelse(max(m_lng_dt$end_mtn) > 15000, 5000, 2000))
    p1 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = sig7_prop))+theme_tree2() +
      coord_cartesian(xlim = c(-1,35))+
      ggtitle(db)+
      scale_color_gradient(low ="blue",high = "yellow", breaks = c(0,100))+
      theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
    p2 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = sig7_prop))+theme_tree2() +
      #geom_tiplab(aes(label = paste0(samples, ' (', Source3, ')')), hjust=0)+
      coord_cartesian(xlim = c(max(m_lng_dt$end_mtn)/14, max(m_lng_dt$end_mtn)))+
      scale_x_continuous(breaks=seq(bin, max(m_lng_dt$end_mtn)+500, bin))+
      scale_color_gradient(low ="blue",high = "yellow", breaks = c(0,100))+
      theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0),"cm"))
    plist[[n]] <- plot_grid(p1,p2, nrow=1, align="h", rel_widths = c(1,2))
  }
  db='DB3'
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
  m_meta_dt <- meta_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id)%>% mutate(l_id = taxa) %>% select(taxa, setdiff(colnames(meta_dt), c("taxa", "lineage_id")))
  m_lng_dt <- lng_dt %>% filter(deadbody == db) %>% rename(taxa = lineage_id) %>% select(-deadbody)
  m_lng_dt <- left_join(m_lng_dt, m_meta_dt, by='taxa')
  db3p1 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = sig7_prop))+theme_tree2() +
    coord_cartesian(xlim = c(-1,35))+
    ggtitle(db)+
    scale_color_gradient(low ="blue",high = "yellow", breaks = c(0,100))+
    theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0.2),"cm"))
  db3p2 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = sig7_prop))+theme_tree2() +
    coord_cartesian(xlim = c(500, 7000))+
    scale_color_gradient(low ="blue",high = "yellow", breaks = c(0,100))+
    scale_x_continuous(breaks=seq(2000, 8000,2000))+
    theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0,0.2,0),"cm"))
  db3p3 <- ggtree(tree) %<+% m_lng_dt + geom_tree(aes(color = sig7_prop))+theme_tree2() +
    coord_cartesian(xlim = c(10000, max(m_lng_dt$end_mtn)))+
    scale_x_continuous(breaks=seq(10000, 60000,20000))+
    scale_color_gradient(low ="blue",high = "yellow", breaks = c(0,100))+
    theme(axis.text = element_text(size=10, angle=45, hjust=1), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0),"cm"))
  pDB3 <- plot_grid(db3p1,db3p2,db3p3, nrow=1, align="h")
  lplot <- plot_grid(plotlist = list(plist[[1]], pDB3, plist[[2]], plist[[4]], plist[[5]], plist[[6]]), nrow=2)
  plot_grid(plotlist = list(lplot, plist[[3]]), nrow=1, rel_widths = c(1.5,1))
}


Draw_totalTree_AsymRatio <- function(meta_dt, lng_dt, nwk_tbl){
  bifurc_tbl <- lng_dt %>% filter(branch_n == 2 & early_desc_n > 4)
  bifurc_tbl <- bifurc_tbl %>% mutate(sub_desc_n1 = apply(bifurc_tbl[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == paste0(x[2],'-1')][1]),
                                      sub_desc_n2 = apply(bifurc_tbl[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == paste0(x[2],'-2')][1]))
  bifurc_tbl <- bifurc_tbl %>% mutate(sub_desc_n1 = ifelse(lineage_id == 'L0', apply(bifurc_tbl[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == 'L1'][1]), sub_desc_n1),
                                      sub_desc_n2 = ifelse(lineage_id == 'L0', apply(bifurc_tbl[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == 'L2'][1]), sub_desc_n2)) 
  bifurc_tbl <- bifurc_tbl %>% mutate(bifurc_ratio = apply(bifurc_tbl[c("sub_desc_n1","sub_desc_n2")], 1, function(x) max(x[1], x[2])/min(x[1], x[2])))
  bifurc_tbl <- bifurc_tbl %>% mutate(ratio_fill = ifelse(bifurc_ratio > 7, 7, bifurc_ratio))
  
  label_height=c(6,12,10,25,15,8.5,13)
  p <- list()
  n=0
  for (db in c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')){
    n=n+1
    tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==db])
    m_bifurc_tbl <- bifurc_tbl %>% filter(deadbody == db) %>% rename(taxa = lineage_id) %>% select(taxa,sub_desc_n1, sub_desc_n2, bifurc_ratio,ratio_fill)
    #fill data of nodes if information of it's branches is unique
    if (db %in% c('DB3','DB6','DB8','DB9','DB10')){
      First_ratio = bifurc_tbl$bifurc_ratio[bifurc_tbl$deadbody == db & bifurc_tbl$lineage_id == 'L0']
      p[[n]] <- ggtree(tree) %<+% m_bifurc_tbl + theme_tree2() +
        geom_label(aes(label=round(bifurc_ratio,2), fill=ratio_fill),hjust=1.1, alpha=0.9, size=4)+
        geom_label(aes(x=-3,y=10, fill=First_ratio), label=round(First_ratio,2), size=4)+
        scale_fill_gradient2(low= "gray50", mid = "yellow", high="red", limits=c(1.1,7), midpoint=2)+
        coord_cartesian(xlim = c(-5,30))+
        scale_x_continuous(breaks=seq(0,30,5))+
        xlab("No. of substitutions")+ggtitle(db)+
        theme(axis.title = element_text(size=10),axis.text = element_text(size=10), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0.2),"cm"))
    } else {
      p[[n]] <- ggtree(tree) %<+% m_bifurc_tbl + theme_tree2() +
        geom_label(aes(label=round(bifurc_ratio,2), fill=ratio_fill),hjust=1.1, alpha=0.9, size=4)+
        scale_fill_gradient2(low= "gray50", mid = "yellow", high="red", limits=c(1.1,7), midpoint=2)+
        coord_cartesian(xlim = c(-5,30))+
        scale_x_continuous(breaks=seq(0,30,5))+
        xlab("No. of substitutions")+ggtitle(db)+
        theme(axis.title = element_text(size=10),axis.text = element_text(size=10), plot.title = element_text(size=20), plot.margin= unit(c(0.2,0.2,0.2,0.2),"cm"))
    }
  }
  lplot <- plot_grid(plotlist=p[c(1,2,3,5,6,7)], nrow=2, align="hv")
  plot_grid(plotlist=list(lplot, p[[4]]), nrow=1, rel_widths = c(1.5,1))
}

Draw_gradient_legend3 <- function(color_v, start_n, change_n, end_n){
  colfunc <- colorRampPalette(rev(color_v[1:2]))
  #plot(1:10, 1:10, pch = 19, cex=2, col = colfunc(5))
  legend_image <- as.raster(matrix(colfunc(5), ncol=1))
  plot(c(0,3),c(0,10),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
  text(x=0.15, y = c(start_n,change_n), labels = c(start_n,change_n), cex=1.5)
  rasterImage(legend_image, 0, start_n, 0.1,change_n)
  colfunc <- colorRampPalette(rev(color_v[2:3]))
  legend_image <- as.raster(matrix(colfunc(25), ncol=1))
  text(x=0.15, y = end_n, labels = end_n, cex=1.5)
  rasterImage(legend_image, 0, change_n, 0.1,end_n)
}

Draw_Dist_Asym_Obs_Sim <- function(lng_dt, sim_res_path){
  bifurc_tbl <- lng_dt %>% filter(branch_n == 2 & early_desc_n > 4)
  bifurc_tbl <- bifurc_tbl %>% mutate(sub_desc_n1 = apply(bifurc_tbl[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == paste0(x[2],'-1')][1]),
                                      sub_desc_n2 = apply(bifurc_tbl[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == paste0(x[2],'-2')][1]))
  bifurc_tbl <- bifurc_tbl %>% mutate(sub_desc_n1 = ifelse(lineage_id == 'L0', apply(bifurc_tbl[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == 'L1'][1]), sub_desc_n1),
                                      sub_desc_n2 = ifelse(lineage_id == 'L0', apply(bifurc_tbl[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == 'L2'][1]), sub_desc_n2)) 
  bifurc_tbl <- bifurc_tbl %>% mutate(bifurc_ratio = apply(bifurc_tbl[c("sub_desc_n1","sub_desc_n2")], 1, function(x) max(x[1], x[2])/min(x[1], x[2])))
  bifurc_tbl <- bifurc_tbl %>% mutate(ratio_fill = ifelse(bifurc_ratio > 7, 7, bifurc_ratio))
  sim_res <- loadRDS(sim_res_path)
  sim_median_v <- Reduce(c,lapply(sim_res, median)) 
  merged_v <- Reduce(c,sim_res)
  merged_v <- merged_v[is.infinite(merged_v)==F]
  dt <- merged_v %>% as_tibble()
  which.max(density(bifurc_tbl$bifurc_ratio)$y) -> y_idx
  density(bifurc_tbl$bifurc_ratio)$x[y_idx] -> x_peak
  x_median <- median(bifurc_tbl$bifurc_ratio)
  sim_median_v[sim_median_v > x_median]
  sd(sim_median_v)
  ggplot(bifurc_tbl)+
    geom_jitter(aes(x=bifurc_ratio, y=0, color=deadbody), size=5, alpha=0.7, width=0, height=0.005)+
    geom_density(aes(bifurc_ratio), color="black")+
    geom_density(data=dt, aes(value), color="darkgray")+
    geom_vline(xintercept = x_median, linetype="dashed", color="black")+
    geom_vline(xintercept = sim_median, linetype="dashed", color="darkgray")+
    #geom_text(aes(x=4, y=1.0), label=round(x_median,2), color="red")+
    #geom_text(aes(x=2, y=1.0), label=round(sim_median,2), color="blue")+
    scale_x_continuous(breaks=c(seq(1,10,1),seq(10,25,5)))+
    coord_cartesian(xlim=c(1,26))+
    scale_color_manual(values = db_pal)+
    ylab("Density")+xlab("Ratio of descendants at bifurcation")+
    theme_bw(base_size=15)+theme(panel.grid = element_blank())
}

Make_Bifurc_Asymmetry_PNG <- function(lng_dt, save_path){
  bifurc_tbl <- lng_dt %>% filter(branch_n == 2 & early_desc_n > 4)
  bifurc_tbl <- bifurc_tbl %>% mutate(sub_desc_n1 = apply(bifurc_tbl[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == paste0(x[2],'-1')][1]),
                                      sub_desc_n2 = apply(bifurc_tbl[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == paste0(x[2],'-2')][1]))
  bifurc_tbl <- bifurc_tbl %>% mutate(sub_desc_n1 = ifelse(lineage_id == 'L0', apply(bifurc_tbl[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == 'L1'][1]), sub_desc_n1),
                                      sub_desc_n2 = ifelse(lineage_id == 'L0', apply(bifurc_tbl[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == 'L2'][1]), sub_desc_n2)) 
  bifurc_tbl <- bifurc_tbl %>% mutate(bifurc_ratio = apply(bifurc_tbl[c("sub_desc_n1","sub_desc_n2")], 1, function(x) max(x[1], x[2])/min(x[1], x[2])))
  for (i in 1:nrow(bifurc_tbl)){
    asym_v <- c(as.numeric(bifurc_tbl[i,"sub_desc_n1"]),as.numeric(bifurc_tbl[i,"sub_desc_n2"]))
    asym_v <- rev(sort(asym_v))
    lid2 <- as.character(bifurc_tbl[i,"lineage_id2"])
    print(lid2)
    dt <- tibble(y=asym_v, labels = c('n1','n2'))
    png(paste0(save_path,lid2,".png"), bg = "transparent")
    g <-  ggplot(dt, aes(x=factor(0), y=y, fill=labels))+
      geom_bar(stat="identity")+
      coord_polar(theta="y", start=-pi/2, direction=1)+
      scale_fill_manual(values = c('n1'= pal_aaas('default')(10)[2], 'n2'=pal_aaas('default')(10)[1]))+
      theme_bw()+
      theme(axis.ticks=element_blank(),
            axis.text=element_blank(),
            axis.title=element_blank(),
            legend.position="none",
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank(),
            panel.background = element_rect(fill = "transparent", colour = NA), 
            plot.background = element_rect(fill = "transparent", colour = NA))
    print(g)
    dev.off()
  }
}

Draw_indiTree_AsymPNG <- function(meta_dt, lng_dt, nwk_tbl, png_path, this_db, width=1.25, height=0.75){
  bifurc_tbl <- lng_dt %>% filter(branch_n == 2 & early_desc_n > 4)
  bifurc_tbl <- bifurc_tbl %>% mutate(sub_desc_n1 = apply(.[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == paste0(x[2],'-1')][1]),
                                      sub_desc_n2 = apply(.[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == paste0(x[2],'-2')][1])) %>% 
    mutate(sub_desc_n1 = ifelse(lineage_id == 'L0', apply(.[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == 'L1'][1]), sub_desc_n1),
           sub_desc_n2 = ifelse(lineage_id == 'L0', apply(bifurc_tbl[c("deadbody", "lineage_id")],1, function(x) lng_dt$early_desc_n[lng_dt$deadbody == x[1] & lng_dt$lineage_id == 'L2'][1]), sub_desc_n2)) %>% 
    mutate(bifurc_ratio = apply(.[c("sub_desc_n1","sub_desc_n2")], 1, function(x) max(x[1], x[2])/min(x[1], x[2]))) %>%
    mutate(ratio_fill = ifelse(bifurc_ratio > 7, 7, bifurc_ratio)) %>% mutate(png_path = paste0(png_path,lineage_id2,'.png'))
  tree <- read.tree(nwk_tbl$nwk_path[nwk_tbl$deadbody==this_db])
  m_bifurc_tbl <- bifurc_tbl %>% filter(deadbody == this_db) %>% rename(taxa = lineage_id) %>% select(taxa,sub_desc_n1, sub_desc_n2, bifurc_ratio,ratio_fill, png_path)
  g <- ggtree(tree) %<+% m_bifurc_tbl + 
    geom_image(aes(image=png_path))+
    coord_cartesian(xlim=c(-2,30))+
    theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"npc"))
  
  library(grid)    
  dev.off()
  print(g, vp=viewport(angle=-90, width = unit(height, "npc"), height=unit(width, "npc")))
}



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

