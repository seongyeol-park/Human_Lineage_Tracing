#liver WGS VAF heatmap using early found by fibros
library(tidyverse)
library(circlize)


#VAF heatmap using early mutaitons
mpu_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db3_liver_WGS/08_DB3_previous_call_readc/DB3_early_pointmts.txt.25sCall')
lng_dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/meta_data/Summary_per_lineage_210322.txt')
colnames(lng_dt)
trunk_n <- lng_dt %>% filter(deadbody == 'DB3' & time_of_lineage == 'early_to_late') %>% nrow()
mpu_dt <- mpu_dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep="_"))
mpu_dt$lineage_id <- str_split_fixed(mpu_dt$lineage_id2, '_',2)[,2]
colnames(mpu_dt)



mx <- mpu_dt %>% select(var_id, ends_with('vafpct'), -starts_with('blood')) %>% as.data.frame() %>%
  column_to_rownames('var_id') %>% as.matrix()

annot_rows <- mpu_dt %>% select(var_id, lineage_id) %>% as.data.frame() %>% column_to_rownames('var_id')
left_anno <- ComplexHeatmap::rowAnnotation(df = annot_rows)

library(ComplexHeatmap)
ComplexHeatmap::Heatmap(mx, left_annotation = left_anno,
                        col = colorRamp2(c(0,50,100), c("blue","yellow","red")),
                       row_names_gp=gpar(fontsize=8), column_names_gp=gpar(fontsize=8)
)



#liver WGS VAF heatmap using late nonUV found by fibros
library(tidyverse)

mpu_dt1 <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db3_liver_WGS/08_DB3_previous_call_readc/DB3_lateNonUV_pointmts.txt.25sCall',
                    col_types=cols(`#CHROM`='c'))
mpu_dt2 <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db3_liver_WGS/08_DB3_previous_call_readc/DB3_lateUV_pointmts.txt.25sCall', 
                    col_types = cols(`#CHROM`='c'))
late_dt <- bind_rows(mpu_dt1, mpu_dt2)
late_dt <- late_dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep="_"))
late_dt$var_id_unique <- make.unique(late_dt$var_id)
mx <- late_dt %>% select(var_id_unique, ends_with('vafpct'), -starts_with('blood')) %>% as.data.frame() %>%
  column_to_rownames('var_id_unique') %>% as.matrix()
nrow(mx)
late_dt$max_vaf <- apply(mx, 1, function(x) max(x))
f_dt <- late_dt %>% filter(max_vaf >0)
nrow(f_dt)

f_dt$lineage_id <- str_split_fixed(f_dt$lineage_id2, '_',2)[,2]
s_dt <- f_dt %>% select(var_id, ends_with('vafpct'), -starts_with('blood'))
ids_list <- c()
for (i in 1:nrow(s_dt)){
  print(i)
  sample_list <- s_dt[i,] %>% gather(-var_id, key='sample',value='vafpct') %>% mutate(sample_id = gsub('DB3_liver_','',gsub('_vafpct','',sample))) %>%
    filter(vafpct >= 5) %>% pull(sample_id)
  if(length(sample_list)==0){
    ids = 'None'
  }else{
    ids = paste(sample_list, collapse='_')
  }  
  ids_list <- c(ids_list, ids)
}
ids_list
f_dt$ids_vaf5m <- ids_list
if(F){
  f_dt %>% saveRDS('/home/users/sypark/00_Project/06_LineageTracing/db3_liver_WGS/08_DB3_previous_call_readc/DB3_late_liverLCM_f_dt.rds')  
}
f_dt$ids_vaf5m %>% unique()
this_id = 'R12'
tmp_dt <- f_dt %>% filter(grepl(this_id, ids_vaf5m)) %>% select(var_id, lineage_id, contains(paste0(this_id,'_')))
colnames(tmp_dt) <- gsub(paste0('DB3_liver_',this_id,'_'),'',colnames(tmp_dt))
tmp_dt %>% filter(var>1) %>%  arrange(desc(vafpct))
#tmp_dt %>% filter(var>2) %>% arrange(desc(vafpct))

f_dt %>% filter(var_id == '4_100793395_G_A') %>% View()



tmp_dt %>% group_by(lineage_id) %>% count() %>% arrange(desc(n))
tmp_dt %>% group
ggplot(tmp_dt)+
  geom_density(aes(x=DB3_liver_L1_vafpct, fill=lineage_id))

f_dt <- mpu_dt %>% filter(rowMax > 5)
f_mx <- f_dt %>% select(var_id, ends_with('vafpct'), -starts_with('blood')) %>% as.data.frame() %>%
  column_to_rownames('var_id') %>% as.matrix()
nrow(f_mx)
#annot_rows <- f_dt %>% select(var_id, lineage_id) %>% as.data.frame() %>% column_to_rownames('var_id')
#left_anno <- ComplexHeatmap::rowAnnotation(df = annot_rows)
library(ComplexHeatmap)
ComplexHeatmap::Heatmap(f_mx, left_annotation = left_anno,
                        col = colorRamp2(c(0,50,100), c("blue","yellow","red")),
                        row_names_gp=gpar(fontsize=8), column_names_gp=gpar(fontsize=8)
)




# liver shared mutation filtered

dt <- read_tsv('/home/users/sypark/00_Project/06_LineageTracing/db3_liver_WGS/06_Substitution/DB3_liver_25s_merged.txt.sampleN.loh_id.g_fi.fi.25sCall',
               col_types =cols(`#CHROM`='c'))
dt <- dt %>% mutate(var_id = paste(`#CHROM`, POS, REF, ALT, sep='_'))
dt <- dt %>% mutate(Nr1=as.integer(str_split_fixed(`dp;Nr1;Nr2;Nr3;Nr3ids`, ';',5)[,2]), Nr3=as.integer(str_split_fixed(`dp;Nr1;Nr2;Nr3;Nr3ids`, ';',5)[,4]))
dt <- dt %>% separate(`blood_dp;read;vaf`, c('blood_dp','blood_var','blood_vaf'), sep=';', convert=T)
dt$blood_vaf <- as.numeric(gsub('%','',dt$blood_vaf))
dt$mVAF <- dt %>% select(ends_with('vafpct')) %>% apply(., 1, function(x) mean(x[x>0 & is.na(x)==F]))
mx <- dt %>% select(var_id, ends_with('vafpct')) %>% as.data.frame() %>% column_to_rownames('var_id') %>%
  as.matrix()
#ComplexHeatmap::Heatmap(mx)
dt$sn_mt0 <- rowSums(is.na(mx) ==F & mx >0)
dt$sn_na <- rowSums(is.na(mx))
dt$sn_mt95 <- rowSums(mx > 95 & is.na(mx)==F)
mx2 <- dt %>% select(var_id, ends_with('var')) %>% as.data.frame() %>% column_to_rownames('var_id') %>%
  as.matrix()
dt$sn_v1m <- rowSums(is.na(mx2) ==F & mx2 >=1)
dt$sn_v2m <- rowSums(is.na(mx2) ==F & mx2 >=2)
dt$sn_v3m <- rowSums(is.na(mx2) ==F & mx2 >=2)

ggplot(dt)+
  geom_histogram(aes(x=mVAF))


##filter false positives by yloss
fp_by_yloss <- dt %>% filter((Nr3 ==34 & `#CHROM` %in% c('X','Y'))| sn_mt95 + sn_na > 22) %>% pull(var_id)
length(fp_by_yloss)
f_dt <- dt %>% filter(!var_id %in% fp_by_yloss) 
nrow(f_dt)

#filter by sn1- sn3
ff_dt <- f_dt %>% filter(sn_v1m- sn_v3m <5)
nrow(ff_dt)

#filter by sn_na
ff_dt <- ff_dt %>% filter(sn_na <5)
nrow(ff_dt)

#VAF distribution
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}



id_list <- ff_dt %>% select(ends_with('_vafpct')) %>% colnames()
id_list <- gsub('_vafpct','',id_list)
t_id = 'DB3_liver_C1'
plist=list()
median_tbl <- tibble(sample_id= as.character(), median_vafpct= as.numeric())
for (i in 1:length(id_list)){
  t_id = id_list[i]
  print(t_id)
  tmp <- ff_dt %>% select(var_id, starts_with(t_id)) 
  colnames(tmp) <- gsub(paste0(t_id, '_'), '', colnames(tmp))
  tmp <- tmp %>% filter(var > 1 & is.na(vafpct) == F)
  #mode_v <- getmode(tmp$vafpct)
  median_v <- median(tmp$vafpct)
  plist[[i]] <- ggplot(tmp, aes(x=vafpct))+
    geom_histogram()+
    scale_x_continuous(limits=c(0,100))+
    geom_vline(xintercept=median(tmp$vafpct), color="red")+
    ggtitle(paste0(t_id, ',', round(median(tmp$vafpct),2)))
  median_tbl[i,] <- list(t_id, median_v)
  }
plot_grid(plotlist = plist, nrow=5)
median_tbl %>% filter(grepl('lane', sample_id)) %>% pull(median_vafpct) %>% mean()
median_tbl <- bind_rows(median_tbl %>% filter(!grepl('lane', sample_id)), tibble(sample_id = 'DB3_liver_R10', median_vafpct = 12.1))
median_tbl$median_vafpct %>% summary()



#median vaf tbl



#heatmap
ff_mx <- ff_dt %>% select(var_id, ends_with('vafpct')) %>% as.data.frame() %>% column_to_rownames('var_id') %>%
  as.matrix()
ComplexHeatmap::Heatmap(ff_mx)

#random sample and draw heatmap
tmp <- ff_dt %>% sample_n(100)
#tmp <- ff_dt %>% filter(Nr3 <= 4)
rt_annot = rowAnnotation(blood = tmp$blood_vaf, mVAF=tmp$mVAF)
tmp_mx <- tmp %>% select(var_id, ends_with('vafpct')) %>% as.data.frame() %>% column_to_rownames('var_id') %>%
  as.matrix()
nrow(tmp_mx)
library(circlize)
ComplexHeatmap::Heatmap(tmp_mx, col = colorRamp2(c(0,50,100), c("blue","yellow","red")),
                        right_annotation = rt_annot, row_names_gp=gpar(fontsize=8), column_names_gp=gpar(fontsize=8)
)


#sigle tone heamap
tmp2 <- ff_dt %>% filter(casNo == 1)
rt_annot = rowAnnotation(blood = tmp2$blood_vaf, mVAF=tmp2$mVAF)
tmp_mx2 <- tmp2 %>% select(var_id, ends_with('vafpct')) %>% as.data.frame() %>% column_to_rownames('var_id') %>%
  as.matrix()
nrow(tmp_mx2)
library(circlize)
ComplexHeatmap::Heatmap(tmp_mx2, col = colorRamp2(c(0,50,100), c("blue","yellow","red")),
                        show_row_names = F,
                        right_annotation = rt_annot, column_names_gp=gpar(fontsize=8)
)
#mutli heatmap
tmp2 <- ff_dt %>% filter(casNo > 1)
rt_annot = rowAnnotation(blood = tmp2$blood_vaf, mVAF=tmp2$mVAF)
tmp_mx2 <- tmp2 %>% select(var_id, ends_with('vafpct')) %>% as.data.frame() %>% column_to_rownames('var_id') %>%
  as.matrix()
nrow(tmp_mx2)
library(circlize)
ComplexHeatmap::Heatmap(tmp_mx2, col = colorRamp2(c(0,50,100), c("blue","yellow","red")),
                        show_row_names = F,
                        right_annotation = rt_annot, column_names_gp=gpar(fontsize=8)
                        
)

phm <- pheatmap::pheatmap(tmp_mx2)
var_groups <- sort(cutree(phm$tree_row, k=8))
table(var_groups)
target_vars <- names(var_groups[var_groups!=1])

tmp3 <- tmp2 %>% filter(var_id %in% target_vars)
rt_annot = rowAnnotation(blood = tmp3$blood_vaf, mVAF=tmp3$mVAF)
tmp_mx3 <- tmp3 %>% select(var_id, ends_with('vafpct')) %>% as.data.frame() %>% column_to_rownames('var_id') %>%
  as.matrix()
nrow(tmp_mx3)
library(circlize)
ComplexHeatmap::Heatmap(tmp_mx3, col = colorRamp2(c(0,50,100), c("blue","yellow","red")),
                        row_names_gp=gpar(fontsize=8),
                        right_annotation = rt_annot, column_names_gp=gpar(fontsize=8)
)


