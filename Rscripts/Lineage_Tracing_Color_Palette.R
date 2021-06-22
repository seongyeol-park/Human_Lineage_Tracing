#color code

# Individual
DB2 = '#DF8F44'
DB3 = '#00A1D5'
DB5 = '#B24745'
DB6 = '#79AF97'
DB8 = '#6A6599'
DB9 = '#80796B'
DB10 = '#374E55'
db_pal = c(DB2, DB3, DB5, DB6, DB8, DB9, DB10)
names(db_pal) <- c('DB2','DB3','DB5','DB6','DB8','DB9','DB10')


#Signature
pal_combi <- c(pal_npg("nrc")(10), pal_d3('category10')(10), pal_aaas('default')(10))
sig_pal <- pal_combi[c(9,11,1,12,5,17,13)]
names(sig_pal) <- c('v3_1','v3_5','v3_7a','v3_7b','v3_7c','v3_7d','v3_18')


#germlayer
gl_pal <- c('endoderm' = '#EBB50A',
            'ectoderm' = '#245DB5',
            'mesoderm' = '#C03C3C',
            'meso_endoderm' = "#e3478b" )


#organs (tissue class3)
tc3_pal<- c("non_internal_organ" = 'grey20',
            "lung" = "#8388D6", 
            "bowel" = "#DB9E51",
            "heart" = "#D9CDE1",
            "large_vessel" = "#DDE265",
            "kidney" = "#D5938F",
            "intestine" = "#CCD4AD",
            "pancreas" = "#B250E1",
            "liver" = "#83DA94",
            "spleen" = "#9EDDFF",
            "blood" = "#73B2CF",
            "stomach" = "#E3478B",
            "ovary" = "#81E4D4",
            "uterus"= "#DC89D2")


#left right
lr_pal <- c(
  'lt' = '#DF1919',
  'rt' =  '#023858', #previous color '#367342'
  'center' = "#B09C85FF",
  'unknown' = "#7F7F7FFF"
)


#cranio caudal
cc_pal <- c(
  'cranio' = '#3f007d',
  'caudal' = '#00441b'
)


#total_pal

total_pal <- c(db_pal, sig_pal, gl_pal, tc3_pal, lr_pal)



#color setting

if(F){
  library(randomcoloR)
  this <- randomColor(count=8, hue="green",luminosity = "dark")
  show_col(this)
  db_pal <- pal_npg('nrc')(10)[c(1:5,9,7)]
  names(db_pal) <- unique(meta_dt$deadbody)
  pmeta_dt$anatomy_class3 %>% unique()
  #show_col(c(pal_npg("nrc")(10), pal_d3('category10')(10)))
  pal_combi <- c(pal_npg("nrc")(10), pal_d3('category10')(10), pal_aaas('default')(10))
  #show_col(pal_combi)
  #show_col(pal_npg("nrc")(10)[c(2,4,6,1,5,8,3,6,9)])
  #show_col(pal_d3('category10')(10))
  #show_col(pal_d3('category10')(10)[c(1,10,5,2,3,9,6,8)])
  #anatomy_cl3_pal <-c(pal_d3('category10')(10)[c(1,10,3,9,2,4,5,6)],'black')
  #names(anatomy_cl3_pal) <- c('lt_UE','lt_LE','rt_UE','rt_LE','trunk','internal_organ','HN')
  #anatomy_cl3_pal <- c('#a6cee3','#1f78b4','#33a02c','#fb9a99','#e31a1c','#ff7f00','#cab2d6','#6a3d9a','#ffff99')
  anatomy_cl3_pal <- pal_combi[c(24,11,2,8,12,17,3,13,18)]
  names(anatomy_cl3_pal) <- c('lt_LE','lt_UE','lt_trunk','rt_LE','rt_UE','rt_trunk','center_trunk','HN','internal_organ')
  #show_col(anatomy_cl3_pal)
  source2_pal <- pal_d3('category10')(6)
  names(source2_pal) <- c('HN','lt_LE','lt_UE','rt_LE','trunk','rt_UE')
  
  group_pal <- pal_combi[c(20,22,12, 23,24)]
  names(group_pal) <- c('mixed','fresh_frozen_epidermis','formalin_fixed_epidermis','endoderm','mesoderm')
  dl2_pal <- c(pal_combi[c(8, 4, 13)], colo_ms, colo_ec)
  show_col(dl2_pal)
  names(dl2_pal) <- c('endoderm','mesoderm','ecto_mesoderm', 'meso_endoderm','ectoderm')
  db_pal <- pal_npg('nrc')(10)[c(1:5,9,7)]
  names(db_pal) <- unique(meta_dt$deadbody)
  lr_pal = c('#367342', '#df1919')
  names(lr_pal) = c('rt','lt')
  cc_pal = c(pal_combi[24], pal_combi[23])
  names(cc_pal) = c('cranio','caudal')
  
  
}

