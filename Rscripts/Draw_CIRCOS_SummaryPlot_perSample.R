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
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_200827.txt')
DB_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_DB_200522.txt')

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
theme_syp = theme_bw() + theme(axis.text = element_text(size=15), axis.title = element_text(size=18), panel.border = element_blank(), axis.line = element_line(), title = element_text(size=18))

#########################################


library(doMC)
registerDoMC(8)
registerDoMC(1)

######Draw CIRCOS
if(F){
  source('/home/users/sypark/01_Python_files/circos_R/CIRCOS_R_SYPark.R')
  id_list <- meta_dt %>% arrange(deadbody, sample_id) %>% pull(sample_id)
  foreach(t_id = id_list) %dopar% {
    print(t_id)
    png(paste0('/home/users/sypark/00_Project/06_LineageTracing/merged_data/CIRCOS/',t_id,'_CIRCOS.png'),width=12, height=8, units = "in", res = 150)
    t_gender_ori = meta_dt$gender[meta_dt$sample_id == t_id]
    t_gender = gsub('F','XX',gsub('M','XY',t_gender_ori))
    t_gender2 = gsub('M','male',gsub('F','female',t_gender_ori))
    SNV_dt <- read_tsv(paste0('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/',t_id,'.edited.vcf'), comment = '##', col_types = cols(`#CHROM` = 'c')) %>%
      rename(CHR = `#CHROM`)
    SNV_dt <- SNV_dt %>% separate(t_id, c('RD','AD','vaf'), sep=':', convert = T) %>% mutate(vaf = vaf*0.01)
    CNV_dt <- read_tsv(paste0('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Sequenza_results/seqz_seg/',t_id,'.seqz_segments_clean.txt'), col_types = cols(`#CHROM`='c')) %>%
      rename(start.pos = start_pos, end.pos = end_pos, A = majCN, B= minCN) %>% mutate(chromosome = paste0('chr', `#CHROM`))
    SV_dt <- read_tsv(paste0('/home/users/sypark/00_Project/06_LineageTracing/merged_data/structural_variations/SV_final_links/',t_id,'.sv_final.txt'), col_types = cols(`#CHR1`='c', POS1 ='n', CHR2='c', POS2='n')) %>%
      rename(chr1 = `#CHR1`, pos1 = POS1, chr2 = CHR2, pos2 = POS2, ori=terinfo) %>% filter(chr2 != '.')
    load(paste0('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Sequenza_results/seqz_extract_RData/',t_id,'.',t_gender,'.100_20_sequenza_extract.RData'))
    seqz_dt <- get(gsub('RMCS-34','RMCS_34',paste0(t_id,'.',t_gender,'.100_20_sequenza_extract')))
    rm(list=ls(pattern=t_id))
    circos_VCS(
      SNV=SNV_dt, # CHR, POS, REF, ALT, vaf, color(optional), shape(optional(default pch=16)), size(optional(default cex=0.5))
      CNV=CNV_dt, # chromosome ("chr1"...), start.pos, end.pos, A (major CN), B(minor CN) (sequenza segements.txt)
      StV=SV_dt, # chr1, pos1, chr2, pos2, ori("5to5","3to3","3to5","5to3")
      seqz=seqz_dt,   # sequenza_extract.RData
      purity=1.0,
      ploidy=2.0,
      spcs = 'hg19', # hg19 or mm10
      gender=t_gender2,  # female or male
      SNV_plottype = 'both',
      SNV_height=0.4,
      AB_height=0.1,
      CNV_height=0.15,
      dot_transparency=0.3
    )
    title(t_id)
    
    dev.off()
  }
  
}

# CIRCOS for three serial_daughter samples
if(F){
  raw_meta_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_sample_200908.txt')
  id_list <- raw_meta_dt %>% filter(Memo == 'serial_daughter') %>% pull(sample_id)
  for(t_id in id_list[2:3]){
    print(t_id)
    png(paste0('/home/users/sypark/00_Project/06_LineageTracing/merged_data/CIRCOS/',t_id,'_CIRCOS.png'),width=12, height=8, units = "in", res = 150)
    t_gender_ori = meta_dt$gender[raw_meta_dt$sample_id == t_id]
    t_gender = gsub('F','XX',gsub('M','XY',t_gender_ori))
    t_gender2 = gsub('M','male',gsub('F','female',t_gender_ori))
    SNV_dt <- read_tsv(paste0('/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/serial_daughters/',t_id,'.edited.vcf'), comment = '##', col_types = cols(`#CHROM` = 'c')) %>%
      rename(CHR = `#CHROM`)
    SNV_dt <- SNV_dt %>% separate(t_id, c('RD','AD','vaf'), sep=':', convert = T) %>% mutate(vaf = vaf*0.01)
    CNV_dt <- read_tsv(paste0('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Sequenza_results/seqz_seg/serial_daughters/',t_id,'.seqz_segments_clean.txt'), col_types = cols(`#CHROM`='c')) %>%
      rename(start.pos = start_pos, end.pos = end_pos, A = majCN, B= minCN) %>% mutate(chromosome = paste0('chr', `#CHROM`))
    SV_dt <- read_tsv(paste0('/home/users/sypark/00_Project/06_LineageTracing/merged_data/structural_variations/SV_final_links/serial_daughters/',t_id,'.sv_final.txt'), col_types = cols(`#CHR1`='c', POS1 ='n', CHR2='c', POS2='n')) %>%
      rename(chr1 = `#CHR1`, pos1 = POS1, chr2 = CHR2, pos2 = POS2, ori=terinfo) %>% filter(chr2 != '.')
    load(paste0('/home/users/sypark/00_Project/06_LineageTracing/merged_data/Sequenza_results/seqz_extract_RData/serial_daughters/',t_id,'.',t_gender,'.100_20_sequenza_extract.RData'))
    seqz_dt <- get(gsub('RMCS-34','RMCS_34',paste0(t_id,'.',t_gender,'.100_20_sequenza_extract')))
    rm(list=ls(pattern=t_id))
    circos_VCS(
      SNV=SNV_dt, # CHR, POS, REF, ALT, vaf, color(optional), shape(optional(default pch=16)), size(optional(default cex=0.5))
      CNV=CNV_dt, # chromosome ("chr1"...), start.pos, end.pos, A (major CN), B(minor CN) (sequenza segements.txt)
      StV=SV_dt, # chr1, pos1, chr2, pos2, ori("5to5","3to3","3to5","5to3")
      seqz=seqz_dt,   # sequenza_extract.RData
      purity=1.0,
      ploidy=2.0,
      spcs = 'hg19', # hg19 or mm10
      gender=t_gender2,  # female or male
      SNV_plottype = 'both',
      SNV_height=0.4,
      AB_height=0.1,
      CNV_height=0.15,
      dot_transparency=0.3
    )
    title(t_id)
    
    dev.off()
  }
  
}



#make CIRCOS + signature plot per sampe
index_dt <- read_tsv("/home/users/sypark/00_Project/06_LineageTracing/meta_data/indexing_by_ysju.txt", col_names = c('index','sample_id'))
CIRCOS_list <- list.files('/home/users/sypark/00_Project/06_LineageTracing/merged_data/CIRCOS/', pattern = 'png$', full.names = T)
tmp_tbl1 <- tibble(CIRCOS_path = CIRCOS_list) %>% mutate(sample_id = apply(., 1, function(x) unlist(strsplit(unlist(strsplit(x,'//'))[2],'_CIRCOS'))[1]))
png_dt <- index_dt %>% left_join(tmp_tbl1)%>% mutate(deadbody = apply(.["sample_id"], 1, function(x) paste0('DB',unlist(strsplit(x,'_'))[1])))
SNVsig_folder = '/home/users/sypark/00_Project/06_LineageTracing/merged_data/point_mutations/perSample_list_link/sig_SBS15187abcd_a1/'
save_folder = '/home/users/sypark/00_Project/06_LineageTracing/merged_data/Summary_plots/'
foreach (n = 1:nrow(png_dt)) %dopar% {
  t_id = png_dt$sample_id[n]
  print(t_id)
  t_index = png_dt$index[n]
  t_db <- paste0('DB',unlist(strsplit(t_id,'_'))[1])
  png(paste0(save_folder, t_index,'.',t_id,'.summary.png'), width=16, height=9, unit="in", res=150)
  CIRCOSpng <- rasterGrob(readPNG(png_dt$CIRCOS_path[png_dt$sample_id == t_id]))
  spec_dt <- read_tsv(paste0(SNVsig_folder,t_id,'.edited.vcf.anv.snv.signature_spectra.tsv'))
  info_dt <- read_tsv(paste0(SNVsig_folder, t_id, '.edited.vcf.anv.snv.INFO.tsv'), col_names="item")
  info_dt <- info_dt %>% separate(item, c('item','value'), sep=': ', convert = T)
  t_CS <- round(info_dt$value[info_dt$item == "# Cosine similarity"],2)
  chg <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  palette <- c(rgb(30,191,240, max = 255),rgb(5,7,8, max = 255),rgb(230,39,37, max = 255),rgb(203,202,203, max = 255),rgb(161,207,100, max = 255),rgb(237,200,197, max = 255))
  names(palette) <- chg
  spec_dt <- spec_dt %>% mutate(context = paste0(Substitution, '_', Trinucleotide))
  g_ct <- ggplot(spec_dt, aes(x=context, y=Original, fill=Substitution))+
    geom_bar(stat="identity")+
    scale_fill_manual(values = palette)+
    scale_x_discrete(labels = spec_dt$Trinucleotide)+
    theme(axis.text.x = element_text(angle=90, vjust=0.5,hjust=1, size=10), axis.text.y = element_text(size=12), strip.text = element_text(size=12),
          legend.position = 'none', axis.title.x = element_blank())+
    ylab("Count")+
    ggtitle(paste0('Cosine similarity: ',t_CS))+
    facet_grid(cols = vars(Substitution), scales = "free_x")
  ex_dt <- read_tsv(paste0(SNVsig_folder,t_id,'.edited.vcf.anv.snv.signature_exposures.tsv'))
  ex_dt <- ex_dt %>% rename(sig = `Signature #`) %>% filter(sig != 'Unexplained') %>% select(-Proportion) %>% mutate(prop = Exposure/sum(.$Exposure))
  ex_dt$sig <- factor(ex_dt$sig, levels = ex_dt %>% arrange(desc(prop)) %>% pull(sig))
  sig_pal2 <- sig_pal
  names(sig_pal2) <- gsub('v3_','SBS',names(sig_pal))
  g_prop <- ggplot(ex_dt, aes(x=sig, y=prop, fill=sig))+
    geom_bar(stat="identity")+
    geom_text(aes(x=sig, y=prop+0.1, label=round(prop,2)), size=4)+
    scale_fill_manual(values=sig_pal2)+
    ggtitle("Signature proportion")+
    theme(axis.title = element_blank(), legend.position = 'none', axis.text.x = element_text(angle=45, hjust=1, size=10))
  g_phylo <- Draw_Phylo_Full_DB_markSample(t_db, t_id)
  g_tbl <- meta_dt %>% filter(sample_id == t_id) %>% select(sample_id, lineage_id_cindel, Cell_type, Source_side, Source_direction, Source_class, Source2, n_snv, n_indel) %>%
    rename(lineage_id = lineage_id_cindel) %>% as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column('item') %>% rename(value = V1)
  
  grid.newpage()
  pushViewport(viewport(x=0, y=0, width=0.6, height=1, just=c('left','bottom')))
  grid.draw(CIRCOSpng)
  popViewport()
  pushViewport(viewport(x=0.5, y=0, width=0.5, height=0.3, just=c('left','bottom')))
  print(g_ct, newpage =F)
  popViewport()
  pushViewport(viewport(x=0.55, y=0.3, width=0.225, height=0.3, just=c('left','bottom')))
  grid.table(g_tbl, rows=NULL, cols=NULL)
  popViewport()
  pushViewport(viewport(x=0.775, y=0.3, width=0.225, height=0.3, just=c('left','bottom')))
  print(g_prop, newpage =F)
  popViewport()
  pushViewport(viewport(x=0.55, y=0.6, width=0.45, height=0.4, just=c('left','bottom')))
  print(g_phylo, newpage=F)
  popViewport()
  dev.off()
}




