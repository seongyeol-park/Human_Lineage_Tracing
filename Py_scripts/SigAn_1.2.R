#!/bin/Rscript
#
# Written by Joonoh Lim (Ju lab; joonoh.lim@kaist.ac.kr)
#
# First created on 2019-05-21
# Last modified on 2019-08-08
#
# Version.1.0.0 2019-05-22
# Version.1.0.1 2019-05-23: default for signature db is set to its absolute path; VCF including indel calls is supported
# Version 1.0.2 2019-05-25: display help message if no input is supplied
# Version 1.0.3 2019-05-31: chromosomes or contigs other than 1~22, X, and Y can now be properly handled (e.g., GL000211.1 -> chrUn_gl000211)
# Version 1.0.4 2019-06-05: added a new argument for output prefix. Fixed error occrring when no extra chromosome is in in vcf
# Version 1.0.5 2019-06-11: suggested by CJYoon, added a functionality to make use of vaf information to give weights to each mutations (i.e., SNVs)
# Version 1.0.6 2019-06-12: corrected a minor bug regarding output directory
# Version 1.0.7 2019-06-17: resolved a bottleneck for computation speed.
# Version 1.0.8 2019-06-20: added mm10 genome.
# Version 1.0.9 2019-07-04: added new arguments to allow input of a table format
# Version 1.1.0 2019-08-08:
#     1) removed chromosomesExtra, added chromosomes_of_interest and chromosome_prefix, for more robustness.
#     2) now it generates one more output, a INFO file which currently only contains: cosine similarity.
#
# Version 1.2.0 2019-08-19:
#     1) removed the argument '--signature_data' and added a new '--signatureSet' 
#     2) modify the argument '--signature_number'. Previously '1:30' should now be 'all'.
#
# A code to do mutational signature anaysis on a single sample using 30 COSMIC mutationsl signatures (or custom signatures).
# Version M.m.n: M -> Major change of the outlook; m -> major changes, backward incompatible; n -> minor changes, backward compatible
#
.libPaths(c("/home/users/jolim/usr/tools/R/library/3.5.1", .libPaths()))
##
suppressMessages(library(optparse))
option_list = list(
    make_option(c("-i","--input"), type="character", default=NULL, help="Input file (supported: VCF in .vcf, .gz, or .bz2; a tab-delimited text file (TSV) with mutation table (column header: Substitution   Trinucleotide   Count); the code automatically detects the type of input, unless -t or --type is specified)",metavar="character"), # 2019-07-04
    make_option(c("-t","--type"), type="character", default=NULL, help="[Optional] A flag for fail-proof determination of the input data type (either 1 or 2; where 1: VCF, 2: TSV; default: NULL).",metavar="character"), # 2019-07-04
    make_option(c("-r","--header"), type="character", default="Substitution,Trinucleotide,Count", help="[Optional] A column header names, comma-separated, used when the input is of table form instead of vcf (default: 'Substitution,Trinucleotide,Count')",metavar="character"), # 2019-07-04
    make_option(c("-o","--outputDir"), type="character", default="./output/", help="[Optional] Output directory (default: ./output/)",metavar="character"),
    make_option(c("-p","--outputPrefix"), type="character", default="", help="[Optional] a prefix for output file (default: \"\")", metavar="character"),
    make_option(c("-x","--outputSuffix"), type="character", default="", help="[Optional] a suffix for output file (default: \"\")", metavar="character"),
    make_option(c("-n","--signature_number"), type="character", default="all", help="[Optional] Signature numbers to fit, comma-separated or colon-ranged or combination of the two. A special input: 'all'. The default behavior changes accordingly with the choice of COSMIC_v2 or COSMIC_v3. (default: 'all')",metavar="character"), # 2019-08-19
    make_option(c("-g","--refGenome"), type="character", default="hg19", help="[Optional] a reference genome name (supported: hg19, hg38, mm10; default: hg19)",metavar="character"),
    make_option(c("-c","--chromosomes_of_interest"), type="character", default="1:22,X,Y", help="[Optional] chromosome names to include for processing, comma-separated and/or colon-ranged (default: 1:22,X,Y)",metavar="character"), # 2019-08-08
    make_option(c("-d","--chromosome_prefix"), type="character", default="chr", help="[Optional] A chromosome name prefix, e.g. because sometimes it is 1 and sometimes chr1. Case-insensitive (e.g. chr and CHR are treated the same). This doesn't affect the case of no prefix (i.e. chromosome name '1' is not affected by '-d chr'), so just leave it as default when no prefix is the case. (def.: chr)",metavar="character"), # 2019-08-08
    make_option(c("-s","--signatureSet"), type="character", default="COSMIC_v2", help="[Optional] Signature data. Currently supported: either \"COSMIC_v2\" or \"COSMIC_v3\" (default: \"COSMIC_v2\", which uses '/home/users/jolim/Projects/S05_Jongsoo_Yoon/code/db/signatures_probabilities.txt')", metavar="character"), # 2019-08-19"
    make_option(c("-v","--vafVarName"), type="character", default=NULL, help="[Optional] the name of the variable that contains VAF information; this is to give weights (:= vaf * 2) to each mutation. Note that the header name for the column containing the variable is assumed to be 'INFO' (supported: VEP format (e.g. <name>=<value> with field delimiter '|' or ';'); default: NULL)", metavar="character")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
##
if(is.null(opt$input)){parse_args(opt_parser, args = c("--help"))}
##
cat("Initiating...\n")
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10)) # 2019-06-20
suppressMessages(library(pracma))
suppressMessages(library(gtools)) # 2019-08-19
suppressMessages(library(RColorBrewer))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
## Get arguments
if (opt$signatureSet == "COSMIC_v2") {
    sigDF <- fread("/home/users/jolim/Projects/S05_Jongsoo_Yoon/code/db/signatures_probabilities.txt")
} else if (opt$signatureSet == "COSMIC_v3") {
    sigDF <- fread("/home/users/jolim/Projects/S05_Jongsoo_Yoon/code/db/sigProfiler_SBS_signatures_2019_05_22.tsv")
} else {stop("Wrong input for the argument '-s' or '--sigantureSet'. It should be either 'COSMIC_v2' or 'COSMIC_v3' where COSMIC_v2 is the previous version of COSMIC signatures (30 signatures) and COSMIC_v3 is the newer version (released May 2019)")}
# 2019-08-19
inputFile <- opt$input
outputDir <- opt$outputDir
if (opt$outputPrefix=="") {outputPrefix <- opt$outputPrefix} else {outputPrefix <- paste0(opt$outputPrefix,".")}
if (opt$outputSuffix=="") {outputSuffix <- opt$outputSuffix} else {outputSuffix <- paste0(".",opt$outputSuffix)}
inputType <- opt$type
inputHeader <- strsplit(opt$header,split=",")[[1]]

# Helper function
setupTable <- function(mutDF = NULL) {

  # make sure it's DF
  mutDF <- as.data.frame(mutDF)

  # helper functions
  c_ntd <- function(ntd){if(ntd=="A"){return("T")};if(ntd=="C"){return("G")};if(ntd=="G"){return("C")};if(ntd=="T"){return("A")}}
  c_chg <- function(chg){if(substr(chg,1,1)=="A"|substr(chg,1,1)=="G"){out<-paste0(c_ntd(substr(chg,1,1)),">",c_ntd(substr(chg,3,3)))} else {out <- chg};return(out)}
  c_seq <- function(seq){temp <- strsplit(seq,split="")[[1]]; out <- paste0(c_ntd(temp[3]),c_ntd(temp[2]),c_ntd(temp[1])); return(out)}

  # abide by C & T convention
  for (i in 1:nrow(mutDF)) {
    toConvert <- substr(mutDF[i,inputHeader[1]],1,1)
    if (toConvert =="A"| toConvert =="G"){
      mutDF[i,inputHeader[1]] <- c_chg(mutDF[i,inputHeader[1]])
      mutDF[i,inputHeader[2]] <- c_seq(mutDF[i,inputHeader[2]])
    }
  }

  # aggregate
  mutDF <- aggregate(formula(paste(inputHeader[3],"~",inputHeader[1],"+",inputHeader[2],sep=" ")),data=mutDF,sum)

  # arrange
  ntd <- c("A","C","G","T")
  chg <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_pre <- rep(ntd,each=4)
  ctx_suf <- rep(ntd,4)
  type <- rep(chg,each=16) 
  subtype <- paste0(rep(ctx_pre,6),c(rep("C",48),rep("T",48)),rep(ctx_suf,6))

  # return
  out <- left_join(data.frame(type=type,subtype=subtype,stringsAsFactors = FALSE),mutDF,by=c("type"=inputHeader[1],"subtype"=inputHeader[2]))
  
  names(out)[names(out)==inputHeader[3]] <- "orig"
  out$orig[is.na(out$orig)] <- 0

  return(out)
}

inputFileName <- gsub(".*/(.*)","\\1",inputFile)
if (is.null(inputType)) {
  cat("Figuring out the type of input data...\n")
  # Automatic determination of input type
  is.vcf <- !is.null(tryCatch({vcf <- fread(inputFile, skip = "#CHROM", sep = "\t")}, error = function(err){}))
  is.tsv <- !is.null(tryCatch({mutDF <- fread(inputFile, skip = paste(inputHeader,collapse="\t"), sep = "\t")}, error = function(err){}))
  if (is.vcf == TRUE) {
    cat("Input data type: VCF\n")
    inputFileName <- gsub("(.*)[.]vcf","\\1",inputFileName)} else
  if (is.vcf == FALSE & is.tsv == TRUE) {
    cat("Input data type: TSV\n")
    mutDF <- setupTable(mutDF)
    inputFileName <- gsub("(.*)[.]tsv$","\\1",inputFileName)
  } else {stop("Unrecognizable input format. Re-check names of the column header, i.e. -r or --header. Else, try -t or --type to manually specify the format.")}
} else if (inputType == "1") {
  is.vcf <- TRUE
  is.tsv <- FALSE
  cat("Input data type: VCF\n")
  vcf <- fread(inputFile, skip = "#CHROM", sep = "\t")
  inputFileName <- gsub("(.*)[.]vcf","\\1",inputFileName)
} else if (inputType == "2") {
  is.vcf <- FALSE
  is.tsv <- TRUE
  cat("Input data type: VCF\n")
  mutDF <- fread(inputFile, skip = paste(inputHeader,collapse="\t"), sep = "\t")
  mutDF <- setupTable(mutDF)
  inputFileName <- gsub("(.*)[.]tsv$","\\1",inputFileName)
}
# To here, 2019-07-04
refGenome <- opt$refGenome
refGenomes <- list(hg19=BSgenome.Hsapiens.UCSC.hg19, hg38=BSgenome.Hsapiens.UCSC.hg38, mm10=BSgenome.Mmusculus.UCSC.mm10)

vafVarName <- opt$vafVarName # 2019-06-11
## Process signature data
cat("Preparing signature data...\n")
if (opt$signature_number == "all" & opt$signatureSet == "COSMIC_v2") {
    sigN <- names(sigDF)
    sigN <- grep("(Substitution Type|Trinucleotide|Signature .*)",sigN)
    sigDF <- sigDF[,..sigN]
    sigDF <- arrange(sigDF, `Substitution Type`, `Trinucleotide`)
} else if (opt$signature_number == "all" & opt$signatureSet == "COSMIC_v3") {
    sigN <- names(sigDF)
    sigN <- grep("(Type|SubType|SBS.*)",sigN,value=TRUE)
    sigDF <- sigDF[,..sigN]
    sigDF <- arrange(sigDF, `Type`, `SubType`)
} else if (opt$signature_number != "all" & opt$signatureSet == "COSMIC_v2") {
    sigN <- opt$signature_number
    sigN <- strsplit(sigN,split=",")[[1]]
    sigN <- as.numeric(unlist(sapply(sigN, function(x){eval(parse(text=x))})))
    sigN <- c("Substitution Type","Trinucleotide",paste0("Signature ",sort(sigN)))
    sigN <- unique(sigN)
    sigDF <- sigDF[,..sigN]
    sigDF <- arrange(sigDF, `Substitution Type`, `Trinucleotide`)
} else if (opt$signature_number != "all" & opt$signatureSet == "COSMIC_v3") {
    sigN <- opt$signature_number
    sigN <- strsplit(sigN,split=",")[[1]]
    sigN.nonrange <- sigN[!grepl(":",sigN)]
    sigN.range <- sigN[grepl(":",sigN)]
    sigN.range <- unlist(sapply(sigN.range, function(x){eval(parse(text=x))}))
    sigN <- c(sigN.nonrange,sigN.range)
    sigN <- c("Type","SubType",paste0("SBS",mixedsort(sigN)))
    sigN <- unique(sigN)
    sigDF <- sigDF[,..sigN]
    sigDF <- arrange(sigDF, `Type`, `SubType`)
} # 2019-08-19
C <- as.matrix(sigDF[,-c(1,2)])

# 2019-07-04
if (is.vcf == TRUE) {
  ## Process input data
  cat("Processing input data...\n")
  vcf <- filter(vcf, nchar(REF)==1 & nchar(ALT)==1) # Remove indels
  vcf <- filter(vcf, (REF=="A"|REF=="C"|REF=="G"|REF=="T") & (ALT=="A"|ALT=="C"|ALT=="G"|ALT=="T")) # Remove other than A,C,G,T
  if (!is.null(vafVarName)){
    vafInfo <- strsplit(vcf[["INFO"]],split=";|\\|") # 2019-06-11
    vafInfo <- unlist(lapply(vafInfo, function(x){grep(vafVarName,x,value=TRUE)})) # 2019-06-11
    ### Sanity check: Stop if vafInfo is not provided: 2019-06-11
    if (length(vafInfo)==0) { stop("INFO column doesn't have the given varaible name. Do the double check.") } # 2019-06-11
    vaf <- as.numeric(unlist(lapply(strsplit(vafInfo,split="="), function(x){x[2]}))) # 2019-06-11
    vcf$weights <- vaf*2 # 2019-06-11
  }
  COI <- strsplit(opt$chromosomes_of_interest,split=",")[[1]] # 2019-08-08
  COI <- unlist(sapply(COI, function(x){if(grepl("[0-9]+:[0-9]+",x)){eval(parse(text=x))}else{x}}),use.names=FALSE) # 2019-08-08
  chromosome_prefix <- opt$chromosome_prefix # 2019-08-08
  chrs_to_include_idx <- gsub(paste0("^",chromosome_prefix),"",vcf$`#CHROM`,ignore.case=TRUE) %in% COI # 2019-08-08
  vcf <- vcf[chrs_to_include_idx,] # 2019-08-08
  #
  str <- as.vector(sapply(vcf[["REF"]], function(x){if (x == "A" | x == "G") {"-"} else {"+"}}))
  if (grep("chr",vcf[["#CHROM"]]) == 1 && unique(substr(vcf[["#CHROM"]],1,3))=="chr") {gr <- GRanges(seqnames = vcf[["#CHROM"]], ranges = IRanges(start = vcf[["POS"]], width = 1), strand = str)} else {
      gr <- GRanges(seqnames = paste0("chr",vcf[["#CHROM"]]), ranges = IRanges(start = vcf[["POS"]], width = 1), strand = str)}
  df <- data.frame(REF=vcf[["REF"]],ALT=vcf[["ALT"]])
  c_ntd_vec <- function(x){A_i<-x=="A";C_i<-x=="C";G_i<-x=="G";T_i<-x=="T";x[A_i]<-"T";x[C_i]<-"G";x[G_i]<-"C";x[T_i]<-"A";return(x)} # 2019-06-17
  
  # 2019-06-17
  REF_A_idx <- df$REF == "A";
  REF_G_idx <- df$REF == "G";
  df$REF[REF_A_idx] <- "T"; df$ALT[REF_A_idx] <- c_ntd_vec(df$ALT[REF_A_idx]); df$REF[REF_G_idx] <- "C"; df$ALT[REF_G_idx] <- c_ntd_vec(df$ALT[REF_G_idx])
  #
  tryCatch({mutCtx <- as.character(getSeq(refGenomes[[refGenome]], gr+1))},
           error = function(err){cat(paste0("ERROR: ",err));cat(paste0("It seems that a wrong reference genome has been chosen.\nCurrent choice: ",refGenome,"\n"));cat("Terminating...\n")
             opt <- options(show.error.messages = FALSE);on.exit(options(opt));stop()},
           warning = function(war){},
           finally = {}
           )
  ### Sanity check
  cat("Checking if the reference genome has been chosen properly...")
  if (all(substr(mutCtx,2,2)==df$REF)) {cat("OK\n")} else {cat("NOT OK\nThe reference genome does not match with REF in the VCF."); stop()}
  ##
  ntd <- c("A","C","G","T")
  chg <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_pre <- rep(ntd,each=4)
  ctx_suf <- rep(ntd,4)
  type <- rep(chg,each=16) 
  subtype <- paste0(rep(ctx_pre,6),c(rep("C",48),rep("T",48)),rep(ctx_suf,6))
  mutDF <- data.frame(type=type,subtype=subtype,stringsAsFactors = FALSE)
  # 2019-06-11: modified to include vaf info
  if (!is.null(vafVarName)) {
    mutDT <- data.table(feature=paste0(df$REF,">",df$ALT,":",mutCtx))
    mutDT$vaf <- vcf$weights
    mutDT <- mutDT[, .(count =.N, weightSum = sum(vaf)), by = feature ]
    mutDT$Freq <- mutDT$count * mutDT$weightSum # Freq column is weighted count
    mutDT <- mutDT[,c("feature","Freq")]
    mutDT <- separate(mutDT,col="feature",into=c("type","subtype"),sep=":")
    mutDF <- left_join(mutDF,mutDT,by=c("type","subtype"))
  } else {
    mutDF <- left_join(mutDF,separate(as.data.frame(table(paste0(df$REF,">",df$ALT,":",mutCtx))),col="Var1",into=c("type","subtype"),sep=":"),by=c("type","subtype"))
  }
  #
  names(mutDF)[names(mutDF)=="Freq"] <- "orig"
  mutDF$orig[is.na(mutDF$orig)] <- 0
}

## Signature analysis
sigNames <- sigN[-c(1,2)]
cat("Running signature extraction...\n")
res <- lsqnonneg(C,mutDF$orig)
exp <- data.frame(Signature = sigNames, exp = res$x, prop = trunc(res$x/sum(res$x)*10^3)/10)
exp <- arrange(exp, desc(exp))
exp$Signature <- factor(exp$Signature, levels = exp$Signature)#, ordered = TRUE)
mutDF$recon <- C %*% res$x
mutDF$resid <- mutDF$orig - mutDF$recon
## Plot
cat("Saving files...\n")
chg <- c("C>A","C>G","C>T","T>A","T>C","T>G") # 2019-07-04
plotMutationalSpectrum <- function(df, features = c("type","subtype"), bar_width = 0.8, alpha = 0.01) {
  # Colors
  palette <- c(rgb(30,191,240, max = 255),rgb(5,7,8, max = 255),rgb(230,39,37, max = 255),rgb(203,202,203, max = 255),rgb(161,207,100, max = 255),rgb(237,200,197, max = 255))
  names(palette) <- chg
  # Plot
  names(df) <- c("type","subtype","Original","Reconstructed","Residual")
  p <- ggplot(data = melt(df, id.vars = features),aes_string(x = features[ features != "type" ], y = "value", fill = "type")) +
    geom_bar(stat = "identity", color = "black", width = bar_width, size = 0.5) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = alpha) +
    scale_fill_manual("Type", values = palette) +
    theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1),strip.text.x = element_text(size = 12, face = "bold"),strip.text.y = element_text(size = 12, face = "bold"),
          axis.title.x=element_blank()) +
    facet_grid(reformulate("type","variable"), scales = "free_x") +
    ylab("Counts")
  return(p)
}
exp.nonzero <- filter(exp, prop > 0)
plotExposure <- function(df, mode = NULL){
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  colors <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  if (mode == "pie"){
    p <- ggplot(df, aes(x = "", y = prop, fill = Signature)) +
      geom_bar(width = 0.5, stat = "identity", color = "black", size = 0.2) +
      coord_polar(theta = "y", start = 0, direction = -1) +
      scale_fill_manual(values = colors[1:nrow(df)], labels = df$Signature) +
      theme_void()
  }
  if (mode == "bar"){
    df <- arrange(df,exp)
    df$Signature <- factor(df$Signature, levels = df$Signature)
    p <- ggplot(df, aes(x=Signature, y=exp, fill=Signature)) +
      geom_bar(stat='identity', show.legend = FALSE, color = "black", size = 0.2) + 
      scale_fill_manual(values = rev(colors[1:nrow(df)])) +
      xlab("Signature #") + ylab("Exposure (Counts)") +
      theme(legend.position="none") +
      scale_x_discrete(limits=df$Signature, labels=gsub("Signature ","",df$Signature)) +
      annotate("label", x = df$Signature, y = df$exp, label = c(paste0(df$prop[df$prop > 0]," %"),rep(NA,sum(df$prop == 0))), size = 2.5, hjust=-0.1) +
      coord_flip(clip = 'off')
  }
  return(p)
}
plotCosineSim <- function(df) {
  A <- as.matrix(mutDF$orig); B <- as.matrix(mutDF$recon); cosSim <- (t(A) %*% B)/sqrt(t(A) %*% A)/sqrt(t(B) %*% B) # 2019-08-08
  text <- paste0("Cosine similarity: ",trunc(cosSim*10^3)/10^3) # 2019-08-08
  p <- ggplot() +
    annotate("text", x = 4, y = 25, size=5, label = text) +
    theme_void()
  return(list(p=p,cosSim=cosSim)) # 2019-08-08
}
plotNucleotideChanges <- function(df){
  df <- aggregate(df$orig, by=list(df$type), sum)
  names(df) <- c("Substitution","Count")
  df$Prop <- trunc(df$Count/sum(df$Count)*10^3)/10
  df$Pos <- cumsum(rev(df$Count)) - rev(df$Count)/2
  palette <- c(rgb(30,191,240, max = 255),rgb(5,7,8, max = 255),rgb(230,39,37, max = 255),rgb(203,202,203, max = 255),rgb(161,207,100, max = 255),rgb(237,200,197, max = 255))
  p <- ggplot(df, aes(x = "", y = Count, fill = Substitution)) +
    geom_bar(width = 0.5, stat = "identity", color = "black", size = 0.2) +
    coord_polar(theta = "y", start = 0, direction = -1) +
    scale_fill_manual(values = palette, labels = df$Substitution) +
    annotate("label", x = 1.25, y = df$Pos, label = rev(paste0(df$Prop," %")), size = 2.5) +
    theme_void()
  return(p)
}

suppressWarnings(p.m <- plotMutationalSpectrum(mutDF))
p.ep <- plotExposure(exp.nonzero, mode = "pie")
suppressWarnings(p.eb <- plotExposure(exp.nonzero, mode = "bar"))
plotCosineSim.out <- plotCosineSim(mutDF) # 2019-08-08
cosSim <- plotCosineSim.out$cosSim # 2019-08-08
p.c <- plotCosineSim.out$p # 2019-08-08
p.chg <- plotNucleotideChanges(mutDF)
pdf(NULL)
p <- grid.arrange(grobs=list(p.chg,p.m,p.eb,p.c),layout_matrix=rbind(c(1,1,2,2,2,2,2,2),c(3,3,2,2,2,2,2,2),c(4,4,2,2,2,2,2,2)),
                  top=textGrob(inputFileName, gp=gpar(fontsize=20,font=2)))
# Save to files
system(paste0("if [ ! -d ",outputDir," ];then mkdir ",outputDir,"; fi")) # 2019-06-12
suppressMessages(ggsave(filename = paste0(outputDir,"/",outputPrefix,inputFileName,".signature_analysis",outputSuffix,".png"), plot = p,
                        width = 45, height = 15, units = "cm",
                        dpi = 300))
names(mutDF) <- c("Substitution","Trinucleotide","Original","Reconstructed","Residual")
fwrite(mutDF,file = paste0(outputDir,"/",outputPrefix,inputFileName,".signature_spectra",outputSuffix,".tsv"),sep="\t")
names(exp) <- c("Signature #","Exposure","Proportion")
fwrite(exp, file = paste0(outputDir,"/",outputPrefix,inputFileName,".signature_exposures",outputSuffix,".tsv"),sep="\t")
# 2019-08-08, INFO file
INFO.file <- paste0(outputDir,"/",outputPrefix,inputFileName,".INFO",outputSuffix,".tsv")
system(paste0("touch ",INFO.file))
info <- paste0("# Cosine similarity: ",cosSim)
system(paste0("echo \"",info,"\" > ",INFO.file))
#
cat("Job completed.\n")

  
