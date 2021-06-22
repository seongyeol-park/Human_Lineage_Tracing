#!/bin/Rscript
#
# Written by Joonoh Lim (Ju lab; joonoh.lim@kaist.ac.kr)
#
# First created on 2019-05-21
# Last modified on 2019-10-21
#
# Version.1.0.0 2019-05-22: the first code has been implemented
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
# Version 1.2.1 2019-08-22: contained the signature data inside the code. No need to care about the file paths.
# Version 1.2.1 2019-08-26: corrected a small bug
# Version 1.3.1 2019-10-21: added a new option "COSMIC_v2v3" to -s --signatureSet (requested by JHPark).
# Version 1.3.2 2019-12-09: The code has been re-factored.
#     1) New argument '-m or --mode' takes either \"SISO\" or \"MISO\" which stand for 'single-input single-output' and 'multiple-input single-output', respectively.
#        The MISO mode takes in multiple files (either VCFs or TSVs or mixed (NOTE: mixed is only allowed when -t --type is NOT specified, i.e., NULL)), merge the counts across input data, and analyze mutational spectra.
#     2) New argument '-z or --log' displays this update log message.
#     3) New argument '-b or --outputBase', if set, overrides the input filename (suggested by CJYoon)
#     4) A small problem with naming output file has been resolved. (reported by CJYoon)
#     5) A small problem with specifying signature_name has been resolved. (reported by CJYoon)
# Version 1.3.3 2019-12-10:
#     1) New argument -e --regions_of_interest has been added (requested by CJYoon)
#     2) Some minor bugs have been fixed.
# Version 1.3.4 2019-12-16:
#     1) A major problem with handling missing SNVs in ROI has ben resolved. (reported by CJYoon)
#     2) A new flag -l or --list has been added, which when ON lets -i or --input takes a text file that contains a list of input file (e.g., vcf) paths. This for the case of arbitrarily many multiple inputs.
#     3) For displaying "Missing input file ...", changed the coloring from \\e[35;5;1m to \\e[38;5;1m
# Version 1.3.5 2019-12-16:
#     1) A minor problem with saving the output into file when SNVs are missing from all samples has been resolved.
# Version 1.3.6 2020-1-29:
#     1) A custom signature in the form of VCF or table (currently only VCF format is supported) (requested by SYPark)
#        NOTE: it should be a table of two columns (signature_name, vcf file path) WITHOUT a header
# Version 2.0.0 2020-1-29:
#     1) Now it generates a new file that contains the input signature data (e.g., "<some names>.signature_individual_spectra.tsv")
#     2) Now the analysis figure (i.e., PNG) contains a new information: total mutation count. (requested by YHAn)
#
# A code to do mutational signature anaysis on a single sample using 30 COSMIC mutationsl signatures (or custom signatures).
# Version M.m.n: M -> Major change of the outlook; m -> major changes, backward incompatible; n -> major or minor changes, backward compatible
#
.libPaths(c("/home/users/jolim/usr/tools/R/library/3.5.1", .libPaths()))
##
suppressMessages(library(optparse))
option_list = list(
    make_option(c("-i","--input"), type="character", default=NULL, help="Input file (supported: VCF in .vcf, .gz, or .bz2; a tab-delimited text file (TSV) with mutation table (column header: Substitution   Trinucleotide   Count); the code automatically detects the type of input, unless -t or --type is specified)",metavar="character"), # 2019-07-04
    make_option(c("-m","--mode"), type="character", default="SISO", help="A mode of taking input(s) (supported: SISO, MISO) (default: SISO)"), # 2019-12-09
    make_option(c("-t","--type"), type="character", default=NULL, help="[Optional] A flag for fail-proof determination of the input data type (either 1 or 2; where 1: VCF, 2: TSV; default: NULL).",metavar="character"), # 2019-07-04
    make_option(c("-r","--header"), type="character", default="Substitution,Trinucleotide,Count", help="[Optional] A column header names, comma-separated, used when the input is of table form instead of vcf (default: 'Substitution,Trinucleotide,Count')",metavar="character"), # 2019-07-04
    make_option(c("-o","--outputDir"), type="character", default="./output/", help="[Optional] Output directory (default: ./output/)",metavar="character"),
    make_option(c("-b","--outputBase"), type="character", default=NULL, help="[Optional] a base name for output file. If NULL, the input file name is used as a base name for the output file. (default: NULL)", metavar="character"),    
    make_option(c("-p","--outputPrefix"), type="character", default="", help="[Optional] a prefix for output file (default: \"\")", metavar="character"),
    make_option(c("-x","--outputSuffix"), type="character", default="", help="[Optional] a suffix for output file (default: \"\")", metavar="character"),
#    make_option(c("-n","--signature_number"), type="character", default="1:30", help="[Optional] Signature numbers to fit, comma-separated or colon-ranged or combination of the two. (default: 1:30)",metavar="character"),
    make_option(c("-s","--signatureSet"), type="character", default="COSMIC_v2", help="[Optional] Signature data. Currently supported: either \"COSMIC_v2\" or \"COSMIC_v3\" or \"COSMIC_v2v3\" (default: \"COSMIC_v2\")", metavar="character"), # 2019-08-19; 2019-10-21
    make_option(c("-n","--signature_number"), type="character", default="all", help="[Optional] Signature numbers to fit, comma-separated or colon-ranged or combination of the two. A special input: 'all'. The default behavior changes accordingly with the choice of COSMIC_v2 or COSMIC_v3 or COSMIC_v2v3. (default: 'all'). If the signature set is set to COSMIC_v2v3, the special input \"all\" is NOT ALLOWED and the syntax is \'v2_<signature_numbers>__v3_<signature_numbers>\' where __ is double underscore (e.g. v2_1,2:5__v3_1:4,7a,7b ). ",metavar="character"), # 2019-08-19
    make_option(c("-g","--refGenome"), type="character", default="hg19", help="[Optional] a reference genome name (supported: hg19, hg38, mm10; default: hg19)",metavar="character"),
    make_option(c("-c","--chromosomes_of_interest"), type="character", default="1:22,X,Y", help="[Optional] chromosome names to include for processing, comma-separated and/or colon-ranged (default: 1:22,X,Y)",metavar="character"), # 2019-08-08
    make_option(c("-e","--regions_of_interest"), type="character", default=NULL, help="[Optional] regions to include for processing, comma-separated and/or colon-ranged (e.g., 1:100000,200000:300000) (default: NULL). NOTE: chromosomes_of_interest should be specified and single for this argument to work properly.",metavar="character"), # 2019-12-10
    make_option(c("-d","--chromosome_prefix"), type="character", default="chr", help="[Optional] A chromosome name prefix, e.g. because sometimes it is 1 and sometimes chr1. Case-insensitive (e.g. chr and CHR are treated the same). This doesn't affect the case of no prefix (i.e. chromosome name '1' is not affected by '-d chr'), so just leave it as default when no prefix is the case. (def.: chr)",metavar="character"), # 2019-08-08
#    make_option(c("-c","--chromosomesExtra"), type="character", default="0", help="[Optional] a flag to include chromosomes or contigs other than 1-22, X, and Y. 1 to include, 0 to ignore. (default: 0)",metavar="character"),
#    make_option(c("-s","--signature_data"), type="character", default="/home/users/jolim/Projects/S05_Jongsoo_Yoon/code/db/signatures_probabilities.txt", help="[Optional] Signature data (default: /home/users/jolim/Projects/S05_Jongsoo_Yoon/code/db/signatures_probabilities.txt)", metavar="character"),
    make_option(c("-u","--custom_signatures"), type="character", default=NULL, help="[Optional] a file containing a table of two columns (signature_name, vcf_file_path). NOTE that names for custom signatures should not overlap with those of COSMIC signatures (currently supported: VCF) (def.: NULL)"), # 2020-1-29
    make_option(c("-v","--vafVarName"), type="character", default=NULL, help="[Optional] the name of the variable that contains VAF information; this is to give weights (:= vaf * 2) to each mutation. Note that the header name for the column containing the variable is assumed to be 'INFO' (supported: VEP format (e.g. <name>=<value> with field delimiter '|' or ';'); default: NULL)", metavar="character"),
    make_option(c("-l","--list"), action="store_true", default=FALSE, help="[Optional] If this flag is used, -i (--input) takes a text file as an input that contains a list of input file paths (This is for the case of arbitrarily many multiple inputs)"),
    make_option(c("-z","--log"), action="store_true", default=FALSE, help="[Optional] Display update logs")
)
infoText <- "\t\tVERSION\t\tDATE\t\tCHANGES\n\t\t1.3.1\t\t19.10.21\t\"COSMIC_v2v3\" is now available for -s --signatureSet. Refer to the above for its usage along with -n --signature_number.\n" # 2019-10-21

opt_parser = OptionParser(option_list = option_list, epilogue = infoText)
opt = parse_args(opt_parser)

# HELPER FUNCTION: stop_quietly
stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}
##
if(opt$log){
  cat(
  "
  #=================#
  # VERSION HISTORY #
  #=================#

  Version.1.0.0 2019-05-22: the first code has been implemented
  Version.1.0.1 2019-05-23: default for signature db is set to its absolute path; VCF including indel calls is supported
  Version 1.0.2 2019-05-25: display help message if no input is supplied
  Version 1.0.3 2019-05-31: chromosomes or contigs other than 1~22, X, and Y can now be properly handled (e.g., GL000211.1 -> chrUn_gl000211)
  Version 1.0.4 2019-06-05: added a new argument for output prefix. Fixed error occrring when no extra chromosome is in in vcf
  Version 1.0.5 2019-06-11: suggested by CJYoon, added a functionality to make use of vaf information to give weights to each mutations (i.e., SNVs)
  Version 1.0.6 2019-06-12: corrected a minor bug regarding output directory
  Version 1.0.7 2019-06-17: resolved a bottleneck for computation speed.
  Version 1.0.8 2019-06-20: added mm10 genome.
  Version 1.0.9 2019-07-04: added new arguments to allow input of a table format
  Version 1.1.0 2019-08-08:
      1) removed chromosomesExtra, added chromosomes_of_interest and chromosome_prefix, for more robustness.
      2) now it generates one more output, a INFO file which currently only contains: cosine similarity.

  Version 1.2.0 2019-08-19:
      1) removed the argument '--signature_data' and added a new '--signatureSet'
      2) modify the argument '--signature_number'. Previously '1:30' should now be 'all'.

  Version 1.2.1 2019-08-22: contained the signature data inside the code. No need to care about the file paths.
  Version 1.2.1 2019-08-26: corrected a small bug
  Version 1.3.1 2019-10-21: added a new option \"COSMIC_v2v3\" to -s --signatureSet (requested by JHPark).
  Version 1.3.2 2019-12-09: The code has been re-factored.
      1) New argument '-m or --mode' takes either \"SISO\" or \"MISO\" which stand for 'single-input single-output' and 'multiple-input single-output', respectively.
         The MISO mode takes in multiple files (either VCFs or TSVs or mixed (NOTE: mixed is only allowed when -t --type is NOT specified, i.e., NULL)), merge the counts across input data, and analyze mutational spectra.
      2) New argument '-z or --log' displays this update log message.
      3) New argument '-b or --outputBase', if set, overrides the input filename (suggested by CJYoon)
      4) A small problem with naming output file has been resolved. (reported by CJYoon)
      5) A small problem with specifying signature_name has been resolved. (reported by CJYoon)
  Version 1.3.3 2019-12-10:
      1) New argument -e --regions_of_interest has been added (requested by CJYoon)
      2) Some minor bugs have been fixed.

  Version 1.3.4 2019-12-16:
      1) A major problem with handling missing SNVs in ROI has ben resolved. (reported by CJYoon)
      2) A new flag -l or --list has been added, which when ON lets -i or --input takes a text file that contains a list of input file (e.g., vcf) paths. This for the case of arbitrarily many multiple inputs.
      3) For displaying \"Missing input file ...\", changed the coloring from \\e[35;5;1m to \\e[38;5;1m
  Version 1.3.5 2019-12-16:
      1) A minor problem with saving the output into file when SNVs are missing from all samples has been resolved.
  Version 1.3.6 2020-1-29:
      1) A custom signature in the form of VCF or table (currently only VCF format is supported) (requested by SYPark)
      NOTE: it should be a table of two columns (signature_name, vcf file path) WITHOUT a header
  Version 2.0.0 2020-1-29:
      1) Now it generates a new file that contains the input signature data (e.g., \"<some names>.signature_individual_spectra.tsv\")
      2) Now the analysis figure (i.e., PNG) contains a new information: total mutation count. (requested by YHAn)

  A code to do mutational signature anaysis on a single sample using 30 COSMIC mutationsl signatures (or custom signatures).
  Version M.m.n: M -> Major change of the outlook; m -> major changes, backward incompatible; n -> major or minor changes, backward compatible
  \n"
  )
  } else # display update log
if(is.null(opt$input)){parse_args(opt_parser, args = c("--help"))} else { # Show help message if input is missing
##
run <- function(){
  
  system(paste0("echo -e \"\\e[01m\\e[38;5;15mMODE: ",opt$mode,"\\e[0m\"")) 
  system(paste0("echo -e \"\\e[01m\\e[38;5;15mINPUT: ",opt$input,"\\e[0m\"")) 
  system("echo -e \"\\e[38;5;15mInitiating...\\e[0m\"")
  suppressMessages(library(tidyverse))
  suppressMessages(library(data.table))
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
  suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10)) 
  suppressMessages(library(pracma))
  suppressMessages(library(gtools)) 
  suppressMessages(library(RColorBrewer))
  suppressMessages(library(grid))
  suppressMessages(library(gridExtra))
  
  ## M0 : Argument handling ----------------------------------------------------
  inputMode <- opt$mode 
  # Sanity checks
  if (inputMode == "MISO" & is.null(opt$outputBase)){
    system("echo -e \"\\e[38;5;1mERROR: If inputMode is \"MISO\", --outputBase (-b) should be specified.\\e[0m\"")
    stop_quietly()
  }
  if (inputMode == "SISO" & (length(strsplit(opt$input,split=",")[[1]]) > 1)){
    system("echo -e \"\\e[38;5;1mERROR: The inputMode is \"SISO\", but multiple files have been input.\\e[0m\"")
    stop_quietly()
  }
  
  
  outputDir <- opt$outputDir
  if (opt$outputPrefix=="") {outputPrefix <- opt$outputPrefix} else {outputPrefix <- paste0(opt$outputPrefix,".")}
  if (opt$outputSuffix=="") {outputSuffix <- opt$outputSuffix} else {outputSuffix <- paste0(".",opt$outputSuffix)}
  
  signature_number <- opt$signature_number
  tmp.signature_number <- unlist(strsplit(signature_number,split="__"))
  if (length(tmp.signature_number)>2){system("echo -e \"\\e[49mUnrecoginzable input for -n --signature_number\\e[0m\""); stop_quietly()} # sanity check
  tmp.v2 <- any(grepl("^v2_",tmp.signature_number))
  tmp.v3 <- any(grepl("^v3_",tmp.signature_number))
  if (tmp.v2 & tmp.v3) {
    if(opt$signatureSet != "COSMIC_v2v3"){system("echo -e \"\\e[38;5;214mWARNING: There is a mismatch between --signatureSet and --signature_number. The latter takes precedence and signatureSet is set to COSMIC_v2v3.\\e[0m\"")} # \e[38;5;214 is a orange color
    signatureSet <- "COSMIC_v2v3"
  } else if (tmp.v2) {
    if(opt$signatureSet != "COSMIC_v2"){system("echo -e \"\\e[38;5;214mWARNING: There is a mismatch between --signatureSet and --signature_number. The latter takes precedence and signatureSet is set to COSMIC_v2.\\e[0m\"")}
    signatureSet <- "COSMIC_v2"
    signature_number <- gsub("^v2_","",signature_number)
  } else if (tmp.v3) {
    if(opt$signatureSet != "COSMIC_v3"){system("echo -e \"\\e[38;5;214mWARNING: There is a mismatch between --signatureSet and --signature_number. The latter takes precedence and signatureSet is set to COSMIC_v3.\\e[0m\"")}
    signatureSet <- "COSMIC_v3"
    signature_number <- gsub("^v3_","",signature_number)
  } else {
    signatureSet <- opt$signatureSet
  } 
  
  if (signatureSet == "COSMIC_v2") {
    sigDF <- as.data.table(read.csv(textConnection(paste0(COSMIC_v2,collapse="\n")),sep="\t",check.names=FALSE))
  } else if (signatureSet == "COSMIC_v3") {
    sigDF <- as.data.table(read.csv(textConnection(paste0(COSMIC_v3,collapse="\n")),sep="\t",check.names=FALSE))
  } else if (signatureSet == "COSMIC_v2v3") {
    sigDF.v2 <- as.data.table(read.csv(textConnection(paste0(COSMIC_v2,collapse="\n")),sep="\t",check.names=FALSE))
    sigDF.v3 <- as.data.table(read.csv(textConnection(paste0(COSMIC_v3,collapse="\n")),sep="\t",check.names=FALSE))
  } else {
    system("echo -e \"\\e[35;5;1mWrong input for the argument '-s' or '--sigantureSet'. It should be either 'COSMIC_v2' or 'COSMIC_v3' where COSMIC_v2 is the previous version of COSMIC signatures (30 signatures) and COSMIC_v3 is the newer version (released May 2019)\\e[0m\"")
    stop_quietly()
    }
  
  # OTHER ARGUMENTS
  refGenome <- opt$refGenome
  refGenomes <- list(hg19=BSgenome.Hsapiens.UCSC.hg19, hg38=BSgenome.Hsapiens.UCSC.hg38, mm10=BSgenome.Mmusculus.UCSC.mm10)
  vafVarName <- opt$vafVarName
  
  ## M1 : Data pre-processing --------------------------------------------------
  # HELPER FUNCTION :: setupTable
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
  # HELPER FUNCTION :: processVCF
  processVCF <- function(vcf = NULL, filtering = TRUE){
    # vcf must be data.table since it is from fread, but to make sure
    if (!is.data.table(vcf)){vcf <- as.data.table(vcf)}
    
    ## Process input data
    system("echo -e \"\\e[38;5;15mProcessing input data...\\e[0m\"")
    vcf <- vcf[nchar(REF)==1 & nchar(ALT)==1 & (REF %in% c("A","C","G","T")) & (ALT %in% c("A","C","G","T")),]

    if (!is.null(vafVarName)){
      vafInfo <- strsplit(vcf[["INFO"]],split=";|\\|")
      vafInfo <- unlist(lapply(vafInfo, function(x){grep(vafVarName,x,value=TRUE)}))
      ### Sanity check: Stop if vafInfo is not provided
      if (length(vafInfo)==0) { stop("INFO column doesn't have the given varaible name. Do the double check.") }
      vaf <- as.numeric(unlist(lapply(strsplit(vafInfo,split="="), function(x){x[2]})))
      vcf$weights <- vaf*2 
    }

    # FILTERING STEP
    if(filtering){

      ## Filtering by chromosome (i.e., #CHROM)
      COI <- strsplit(opt$chromosomes_of_interest,split=",")[[1]] 
      COI <- unlist(sapply(COI, function(x){if(grepl("[0-9]+:[0-9]+",x)){eval(parse(text=x))}else{x}}),use.names=FALSE) 
      chromosome_prefix <- opt$chromosome_prefix
      chrs_to_include_idx <- gsub(paste0("^",chromosome_prefix),"",vcf$`#CHROM`,ignore.case=TRUE) %in% COI 
      vcf <- vcf[chrs_to_include_idx,]
  
      ## Filtering by region (i.e., POS)
      if (!is.null(opt$regions_of_interest)){
          # Sanity check: whether a single COI is specified
          if(length(COI)!=1){stop("-c (--chromosomes_of_interest) must be single when -e (--regions_of_interest) is specified.")}
          #
          ROI <- strsplit(opt$regions_of_interest,split=",")[[1]] 
          lapply(ROI,function(roi){
              tmp.range <- strsplit(roi,split="(:|-)")[[1]]
              if(length(tmp.range)!=2){stop("Invalid input for -e or --regions_of_interest")} # sanity check: whether the format is num:num or num-num
              tmp.start <- tmp.range[1]
              tmp.end <- tmp.range[2]
              fc <- paste0("(",tmp.start," < as.numeric(POS) & as.numeric(POS) < ",tmp.end,")")
          }) -> tmp.ROI.fc

          ROI.fc <- paste0(tmp.ROI.fc,collapse=" | ")

          vcf <- eval(parse(text=paste0("vcf[",ROI.fc,",]")))
      }

    }

    # EXCEPTION HANDLING EXCEPTION HANDLING EXCEPTION HANDLING, 
    if(nrow(vcf)==0){
      ntd <- c("A","C","G","T")
      chg <- c("C>A","C>G","C>T","T>A","T>C","T>G")
      ctx_pre <- rep(ntd,each=4)
      ctx_suf <- rep(ntd,4)
      type <- rep(chg,each=16) 
      subtype <- paste0(rep(ctx_pre,6),c(rep("C",48),rep("T",48)),rep(ctx_suf,6))
      mutDF <- data.frame(type=type,subtype=subtype,stringsAsFactors = FALSE)
      mutDF[["orig"]] <- 0
      return(mutDF)
    } # otherwise, proceed
    # EXCEPTION HANDLING EXCEPTION HANDLING EXCEPTION HANDLING,

    str <- as.vector(sapply(vcf[["REF"]], function(x){if (x == "A" | x == "G") {"-"} else {"+"}}))
    if (grep("chr",vcf[["#CHROM"]]) == 1 && unique(substr(vcf[["#CHROM"]],1,3))=="chr") {gr <- GRanges(seqnames = vcf[["#CHROM"]], ranges = IRanges(start = vcf[["POS"]], width = 1), strand = str)} else {
      gr <- GRanges(seqnames = paste0("chr",vcf[["#CHROM"]]), ranges = IRanges(start = vcf[["POS"]], width = 1), strand = str)}
    df <- data.frame(REF=vcf[["REF"]],ALT=vcf[["ALT"]])
    c_ntd_vec <- function(x){A_i<-x=="A";C_i<-x=="C";G_i<-x=="G";T_i<-x=="T";x[A_i]<-"T";x[C_i]<-"G";x[G_i]<-"C";x[T_i]<-"A";return(x)} 
    
    REF_A_idx <- df$REF == "A";
    REF_G_idx <- df$REF == "G";
    df$REF[REF_A_idx] <- "T"; df$ALT[REF_A_idx] <- c_ntd_vec(df$ALT[REF_A_idx]); df$REF[REF_G_idx] <- "C"; df$ALT[REF_G_idx] <- c_ntd_vec(df$ALT[REF_G_idx])
    #
    tryCatch({mutCtx <- as.character(getSeq(refGenomes[[refGenome]], gr+1))},
             error = function(err){system(paste0("echo -e \"\\e[41mERROR: ",err));system(paste0("echo -e \"It seems that a wrong reference genome has been chosen.\nCurrent choice: ",refGenome,"\\e[0m\""));system("echo -e \"Terminating...\"")
               opt <- options(show.error.messages = FALSE);on.exit(options(opt));stop_quietly()},
             warning = function(war){},
             finally = {}
    )
    ### Sanity check
    system("echo -e \"\\e[38;5;15mChecking if the reference genome has been chosen properly...\\c\\e[0m\"") # \c --> echo without newline
    if (all(substr(mutCtx,2,2)==df$REF)) {system("echo \"OK\"")} else {system("echo -e \"NOT OK\nThe reference genome does not match with REF in the VCF.\""); stop_quietly()}
    ##
    ntd <- c("A","C","G","T")
    chg <- c("C>A","C>G","C>T","T>A","T>C","T>G")
    ctx_pre <- rep(ntd,each=4)
    ctx_suf <- rep(ntd,4)
    type <- rep(chg,each=16) 
    subtype <- paste0(rep(ctx_pre,6),c(rep("C",48),rep("T",48)),rep(ctx_suf,6))
    mutDF <- data.frame(type=type,subtype=subtype,stringsAsFactors = FALSE)
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
    names(mutDF)[names(mutDF)=="Freq"] <- "orig"
    mutDF$orig[is.na(mutDF$orig)] <- 0
    
    return(mutDF)
  }
  
  # CUSTOM SIGNATURE Pt.1: START
  if(!is.null(opt$custom_signatures)){
    custom_signatures_list_path <- opt$custom_signatures                                                                                                                                                                          
    custom_signatures_vcf_list <- fread(custom_signatures_list_path, header=FALSE)                                                                                                                                                
    if(ncol(custom_signatures_vcf_list)!=2){                                                                                                                                                                                      
      custom_signatures_vcf_list <- fread(custom_signatures_list_path, header=FALSE, sep="\t")                                                                                                                                    
    }                                                                                                                                                                                                                             
    if(ncol(custom_signatures_vcf_list)!=2){                                                                                                                                                                                      
      custom_signatures_vcf_list <- fread(custom_signatures_list_path, header=FALSE, sep=",")                                                                                                                                     
    }                                                                                                                                                                                                                             
    custom_signatures_vcf_list <- set_names(custom_signatures_vcf_list, c("SIGNAME", "FILE"))      
    lapply(1:nrow(custom_signatures_vcf_list), function(i){
        tmp.vcf <- fread(custom_signatures_vcf_list$FILE[i])
        tmp.df <- processVCF(vcf = tmp.vcf, filtering = FALSE)
        tmp.df <- set_names(tmp.df, c("type","subtype",custom_signatures_vcf_list$SIGNAME[i]))
        if (i==1) {tmp.df} else {tmp.df[,3]}
    }) -> custom_signatures_table_list
    custom_signatures_table <- do.call(cbind, custom_signatures_table_list)
    if (ncol(custom_signatures_table)==3) { 
      # i.e., if only one custom signature
      custom_signatures_table[,3] <- custom_signatures_table[,3]/sum(custom_signatures_table[,3])
    } else {
      custom_signatures_table[,3:ncol(custom_signatures_table)] <- custom_signatures_table[,3:ncol(custom_signatures_table)]/colSums(custom_signatures_table[,3:ncol(custom_signatures_table)])
    }
  }
  # CUSTOM SIGNATURE Pt.1: END


  # Submodule 1 : SISO
  S1_dataProcessing_SISO <- function(inputFile = NULL, inputType = NULL){
    
    ## Get arguments
    inputHeader <- strsplit(opt$header,split=",")[[1]]
    
    inputFileName <- gsub(".*/(.*)","\\1",inputFile)
    if (is.null(inputType)) {
      system("echo -e \"\\e[38;5;15mFiguring out the type of input data...\\e[0m\"")
      # Automatic determination of input type
      is.vcf <- !is.null(tryCatch({vcf <- fread(inputFile, skip = "#CHROM", sep = "\t")}, error = function(err){}))
      is.tsv <- !is.null(tryCatch({mutDF <- fread(inputFile, skip = paste(inputHeader,collapse="\t"), sep = "\t")}, error = function(err){}))
      
      if (is.vcf == TRUE) {
        system("echo -e \"\\e[38;5;15mInput data type: VCF\\e[0m\"")
        inputFileName <- gsub("(.*)[.]vcf$","\\1",gsub("(.*)[.]gz$","\\1",inputFileName))
        
      } else if (is.vcf == FALSE & is.tsv == TRUE) {
        system("echo -e \"\\e[38;5;15mInput data type: TSV\\e[0m\"")
        mutDF <- setupTable(mutDF)
        inputFileName <- gsub("(.*)[.]tsv$","\\1",inputFileName)
        
      } else {
        system("echo -e \"\\e[38;5;1mMissing input file or unrecognizable input format. In case of the latter, re-check names of the column header, i.e. -r or --header. Else, try -t or --type to manually specify the format.\\e[0m\"")
        stop_quietly()}
    } else if (inputType == "1") {
      is.vcf <- TRUE
      is.tsv <- FALSE 
      system("echo -e \"\\e[38;5;15mInput data type: VCF\\e[0m\"")
      vcf <- fread(inputFile, skip = "#CHROM", sep = "\t")
      inputFileName <- gsub("(.*)[.]vcf$","\\1",gsub("(.*)[.]gz$","\\1",inputFileName))
      
    } else if (inputType == "2") {
      is.vcf <- FALSE
      is.tsv <- TRUE
      system("echo -e \"\\e[38;5;15mInput data type: VCF\\e[0m\"")
      mutDF <- fread(inputFile, skip = paste(inputHeader,collapse="\t"), sep = "\t")
      mutDF <- setupTable(mutDF)
      inputFileName <- gsub("(.*)[.]tsv$","\\1",inputFileName)
      
    }
    
    if (is.null(opt$outputBase)) {
      outputBase <- inputFileName
    } else {
      outputBase <- opt$outputBase
    }
    
    if (is.vcf == TRUE) {mutDF <- processVCF(vcf = vcf)}
    
    # RETURN
    return(list(outputBase = outputBase,
                mutDF = mutDF))
    
  } 

  # Submodule 2 : MISO
  S2_dataProcessing_MISO <- function(inputFiles = NULL, inputType = NULL){
    
    ## Get arguments
    inputHeader <- strsplit(opt$header,split=",")[[1]]
    
    lapply(1:length(inputFiles), function(i){
      system(paste0("echo -e \"\\e[38;5;40mProcessing input files: ",i,"/",length(inputFiles),"...\\e[0m\""))
      tmp <- S1_dataProcessing_SISO(inputFile = inputFiles[i], inputType = inputType)
      if(i==1){outputBase <<- tmp[["outputBase"]]}
      return(tmp[["mutDF"]])
    }) -> tmp
    mutDF <- Reduce(function(a,b){aggregate(orig ~ type + subtype, rbind(a,b), sum)}, tmp)
    mutDF <- arrange(mutDF,type,subtype)
    
    # RETURN
    return(list(
      outputBase = outputBase,
      mutDF = mutDF
    ))
    
  } 	
  # M1's Main, 2019-12-09
  # NOTE: if inputType is NULL, each inputFile is assessed separately whether it is vcf or tsv, meaning it is possible to mix the input formats.
  #       But, it inputType is set to either one, all the files in inputFiles are assumed to be in the specified format.
  inputType <- opt$type
  if (inputMode == "SISO") {
    
    if(opt$list){
      system("echo -e \"\\e[38;5;214mWARNING: The flag -l --list is on, but the mode is SISO. Consider not using it when the mode is SISO but using a usual single data input (i.e., vcf(.gz) or tsv).\\e[0m\"")
      inputFiles <- fread(opt$input, header=FALSE)[[1]]
      ## Sanity check
      if(length(inputFiles)>1){
        system("echo -e \"\\e[38;5;214mWARNING: The mode is SISO, but the number of input file is greater than 1. Launching the MISO procedure instead...\\e[0m\"")
        tmp <- S2_dataProcessing_MISO(inputFiles = inputFiles, inputType = inputType)
      } else {
        tmp <- S1_dataProcessing_SISO(inputFile = inputFiles, inputType = inputType)
      }
    } 

    if(!opt$list){
      inputFile <- opt$input
      inputType <- opt$type
      tmp <- S1_dataProcessing_SISO(inputFile = inputFile, inputType = inputType)
    }  

  } else 
      
  if (inputMode == "MISO") {
   
    if (opt$list){
      inputFiles <- fread(opt$input, header=FALSE)[[1]]
    } 
    if (!opt$list){
      inputFiles <- strsplit(opt$input,split=",")[[1]]
    } 
    
    ## sanity check
    if(length(inputFiles)<2){
      system("echo -e \"\\e[38;5;214mWARNING: There is a mismatch between --input and --mode. The mode is MISO, but the number of input file is 1. Launching the SISO procedure instead...\\e[0m\"")
      tmp <- S1_dataProcessing_SISO(inputFile = inputFiles, inputType = inputType)
    } else {
      tmp <- S2_dataProcessing_MISO(inputFiles = inputFiles, inputType = inputType) 
    }

  } else {system("echo -e \"\\e[49mInvalid input for for --inputMode. Currently supported are \"SISO\" and \"MISO\".\"");stop_quietly()}

  outputBase <- tmp[["outputBase"]]
  mutDF <- tmp[["mutDF"]]
	
	## M2 : Signature data processing ----------------------------------------------------
	system("echo -e \"\\e[38;5;15mPreparing signature data...\\e[0m\"")
 
	if (signature_number == "all" & signatureSet == "COSMIC_v2") {
	    sigN <- names(sigDF)
	    sigN <- grep("(Substitution Type|Trinucleotide|Signature .*)",sigN,value=TRUE) # 2019-08-26
	    sigDF <- sigDF[,..sigN]
	    sigDF <- arrange(sigDF, `Substitution Type`, `Trinucleotide`)
	} else if (signature_number == "all" & signatureSet == "COSMIC_v3") {
	    sigN <- names(sigDF)
	    sigN <- grep("(Type|SubType|SBS.*)",sigN,value=TRUE)
	    sigDF <- sigDF[,..sigN]
	    sigDF <- arrange(sigDF, `Type`, `SubType`)
	} else if (signature_number != "all" & signatureSet == "COSMIC_v2") {
	    sigN <- signature_number
	    sigN <- strsplit(sigN,split=",")[[1]]
	    sigN <- as.numeric(unlist(sapply(sigN, function(x){eval(parse(text=x))})))
	    sigN <- c("Substitution Type","Trinucleotide",paste0("Signature ",sort(sigN)))
	    sigN <- unique(sigN)
	    sigDF <- sigDF[,..sigN]
	    sigDF <- arrange(sigDF, `Substitution Type`, `Trinucleotide`)
	} else if (signature_number != "all" & signatureSet == "COSMIC_v3") {
	    sigN <- signature_number
	    sigN <- strsplit(sigN,split=",")[[1]]
	    sigN.nonrange <- sigN[!grepl(":",sigN)]
	    sigN.range <- sigN[grepl(":",sigN)]
	    sigN.range <- unlist(sapply(sigN.range, function(x){eval(parse(text=x))}))
	    sigN <- c(sigN.nonrange,sigN.range)
	    sigN <- c("Type","SubType",paste0("SBS",mixedsort(sigN)))
	    sigN <- unique(sigN)
	    sigDF <- sigDF[,..sigN]
	    sigDF <- arrange(sigDF, `Type`, `SubType`)
    } else if (signature_number == "all" & signatureSet == "COSMIC_v2v3") {
    # NOTE: (signature_number == "all" & opt$signatureSet == "COSMIC_v2v3") is NOT ALLOWED
      system("echo -e \"\\e[35;5;1mWhen signatureSet is COSMIC_v2v3, the special input \"all\" for signature_number is NOT allowed.\\e[0m\"")
      stop_quietly()
    } else if (signature_number != "all" & signatureSet == "COSMIC_v2v3") {
        sigN <- signature_number
        sigN <- sort(strsplit(sigN,split="__")[[1]]) # "v2_<>" "v3_<>"
        # Sanity check
        if(!(grepl("^v2_",sigN[1]) & grepl("^v3_",sigN[2]))) {
          system("echo -e \"\\e[35;5;1mWrong input for signature_number.\\e[0m\"")
          stop_quietly()
          }
        # v2
        sigN.v2 <- gsub("^v2_(.*)","\\1",sigN[1])
	    sigN.v2 <- strsplit(sigN.v2,split=",")[[1]]
	    sigN.v2 <- as.numeric(unlist(sapply(sigN.v2, function(x){eval(parse(text=x))})))
	    sigN.v2 <- c("Substitution Type","Trinucleotide",paste0("Signature ",sort(sigN.v2)))
	    sigN.v2 <- unique(sigN.v2)
	    sigDF.v2 <- sigDF.v2[,..sigN.v2]
	    sigDF.v2 <- arrange(sigDF.v2, `Substitution Type`, `Trinucleotide`)
        # v3
        sigN.v3 <- gsub("^v3_(.*)","\\1",sigN[2])
	    sigN.v3 <- strsplit(sigN.v3,split=",")[[1]]
	    sigN.v3.nonrange <- sigN.v3[!grepl(":",sigN.v3)]
	    sigN.v3.range <- sigN.v3[grepl(":",sigN.v3)]
	    sigN.v3.range <- unlist(sapply(sigN.v3.range, function(x){eval(parse(text=x))}))
	    sigN.v3 <- c(sigN.v3.nonrange,sigN.v3.range)
	    sigN.v3 <- c("Type","SubType",paste0("SBS",mixedsort(sigN.v3)))
	    sigN.v3 <- unique(sigN.v3)
	    sigDF.v3 <- sigDF.v3[,..sigN.v3]
	    sigDF.v3 <- arrange(sigDF.v3, `Type`, `SubType`)
        # v2 + v3
        sigDF <- cbind(sigDF.v3[,c(1,2)],sigDF.v2[,-c(1,2)],sigDF.v3[,-c(1,2)])
        # Update sigN
        sigN <- names(sigDF)
    } else {
      system("echo -e \"\\e[35;5;1mWrong inputs or combination for signature set and signature numbers.\\e[0m\"")
      stop_quietly()
    } 
	

	# CUSTOM SIGNATURE Pt.2: START
    if(!is.null(opt$custom_signatures)){
      sigN <- c(names(custom_signatures_table), names(sigDF)[-c(1,2)])
      sigDF <- cbind(custom_signatures_table, sigDF[,-c(1,2)])
      sigDF <- set_names(sigDF, sigN)
	}
	# CUSTOM SIGNATURE Pt.2: END

	C <- as.matrix(sigDF[,-c(1,2)])
	
	## M3 : Signature analysis --------------------------------------------------------
	sigNames <- sigN[-c(1,2)]
	system("echo -e \"\\e[38;5;15mRunning signature extraction...\\e[0m\"")
	res <- lsqnonneg(C,mutDF$orig)

    exp <- data.frame(Signature = sigNames, exp = res$x, prop = if(sum(res$x)!=0){trunc(res$x/sum(res$x)*10^3)/10}else{NA}) # 2019-12-16
	exp <- arrange(exp, desc(exp))
	exp$Signature <- factor(exp$Signature, levels = exp$Signature)#, ordered = TRUE)
	mutDF$recon <- C %*% res$x
	mutDF$resid <- mutDF$orig - mutDF$recon

	## Plot
    system("echo -e \"\\e[38;5;15mSaving files...\"")
	chg <- c("C>A","C>G","C>T","T>A","T>C","T>G") 
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
	      coord_flip(clip = 'off')

          # EXCEPTION HANDLING EXCEPTION HANDLING EXCEPTION HANDLING 
	      if(!is.na(sum(df$prop))){
            p <- p + annotate("label", x = df$Signature, y = df$exp, label = c(paste0(df$prop[df$prop > 0]," %"),rep(NA,sum(df$prop == 0))), size = 2.5, hjust=-0.1)
          } else {}
          # EXCEPTION HANDLING EXCEPTION HANDLING EXCEPTION HANDLING 

	  }
	  return(p)
	}
	plotNucleotideChanges <- function(df){
	  df <- aggregate(df$orig, by=list(df$type), sum)
	  names(df) <- c("Substitution","Count")
      # EXCEPTION HANDLING EXCEPTION HANDLING EXCEPTION HANDLING 
      if(sum(df$Count)!=0){
        df$Prop <- trunc(df$Count/sum(df$Count)*10^3)/10
      } else {df$Prop <- NA}
      # EXCEPTION HANDLING EXCEPTION HANDLING EXCEPTION HANDLING 
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
	plotCosineSim <- function(df) {
	  A <- as.matrix(mutDF$orig)
      B <- as.matrix(mutDF$recon)
      # EXCEPTION HANDLING EXCEPTION HANDLING EXCEPTION HANDLING, 
      if(sum(A)==0|sum(B)==0){
        cosSim <- NA
      } else {
        cosSim <- (t(A) %*% B)/sqrt(t(A) %*% A)/sqrt(t(B) %*% B) 
      }
      # EXCEPTION HANDLING EXCEPTION HANDLING EXCEPTION HANDLING, 
	  text <- paste0("Cosine similarity: ",trunc(cosSim*10^3)/10^3) 
	  p <- ggplot() +
	    annotate("text", x = 0, y = 0, size=5, label = text) +
	    theme_void()
	  return(list(p=p,cosSim=cosSim)) 
	}
    plotTotalCounts <- function(df){
       tot <- sum(mutDF$orig)
       text <- paste0("Total mutation counts: ",tot)
       p <- ggplot() +
         annotate("text", x=0, y=0, size=5, label=text) +
         theme_void()
       return(list(p=p,totCount=tot))
    }
	
	suppressWarnings(p.m <- plotMutationalSpectrum(mutDF))
    # EXCEPTION HANDLING EXCEPTION HANDLING EXCEPTION HANDLING, 
	if(!is.na(sum(exp$prop))){exp.nonzero <- filter(exp, prop > 0)}
    if(is.na(sum(exp$prop))){exp.nonzero <- exp} # Note that it is not actually exp.nonzero, but regarded as if for convenience
    # EXCEPTION HANDLING EXCEPTION HANDLING EXCEPTION HANDLING, 
    
    p.ep <- plotExposure(exp.nonzero, mode = "pie")
    
    suppressWarnings(p.eb <- plotExposure(exp.nonzero, mode = "bar"))
    
    plotCosineSim.out <- plotCosineSim(mutDF) 
	cosSim <- plotCosineSim.out$cosSim 
	p.c <- plotCosineSim.out$p 

	p.chg <- plotNucleotideChanges(mutDF)

    plotTotalCounts.out <- plotTotalCounts(mutDF) 
    p.tot <- plotTotalCounts.out$p
    totCount <- plotTotalCounts.out$totCount

	pdf(NULL)
	p <- grid.arrange(grobs=list(p.chg,p.m,p.eb,p.tot,p.c),
                      layout_matrix=rbind(c(1,1,2,2,2,2,2,2),
                                          c(1,1,2,2,2,2,2,2),
                                          c(1,1,2,2,2,2,2,2),
                                          c(1,1,2,2,2,2,2,2),
                                          c(3,3,2,2,2,2,2,2),
                                          c(3,3,2,2,2,2,2,2),
                                          c(3,3,2,2,2,2,2,2),
                                          c(3,3,2,2,2,2,2,2),
                                          c(NA,NA,2,2,2,2,2,2),
                                          c(4,4,2,2,2,2,2,2),
                                          c(5,5,2,2,2,2,2,2),
                                          c(NA,NA,2,2,2,2,2,2)
                                          ),
	                  top=textGrob(outputBase, gp=gpar(fontsize=20,font=2)))
	# Save to files
	system(paste0("if [ ! -d ",outputDir," ];then mkdir -p ",outputDir,"; fi")) 
	suppressMessages(ggsave(filename = paste0(outputDir,"/",outputPrefix,outputBase,".signature_analysis",outputSuffix,".png"), plot = p,
	                        width = 45, height = 15, units = "cm",
	                        dpi = 300))
    names(sigDF) <- c("Substitution","Trinucleotide",sigNames) 
    fwrite(sigDF,file = paste0(outputDir,"/",outputPrefix,outputBase,".signature_individual_spectra",outputSuffix,".tsv"),sep="\t") 
	names(mutDF) <- c("Substitution","Trinucleotide","Original","Reconstructed","Residual")
	fwrite(mutDF,file = paste0(outputDir,"/",outputPrefix,outputBase,".signature_spectra",outputSuffix,".tsv"),sep="\t")
	names(exp) <- c("Signature #","Exposure","Proportion")
	fwrite(exp, file = paste0(outputDir,"/",outputPrefix,outputBase,".signature_exposures",outputSuffix,".tsv"),sep="\t")
	INFO.file <- paste0(outputDir,"/",outputPrefix,outputBase,".INFO",outputSuffix,".tsv")
	system(paste0("touch ",INFO.file))
	info <- paste0("# Cosine similarity: ",cosSim,"\n",
                   "# Total mutation count: ",totCount)
	system(paste0("echo \"",info,"\" > ",INFO.file))
	#
	system("echo -e \"\\e[42mJob completed.\\e[0m\"")
    #}}}

}
# Signature Data :: COSMIC_v2
# {{{
COSMIC_v2 <- c(
"Substitution Type\tTrinucleotide\tSomatic Mutation Type\tSignature 1\tSignature 2\tSignature 3\tSignature 4\tSignature 5\tSignature 6\tSignature 7\tSignature 8\tSignature 9\tSignature 10\tSignature 11\tSignature 12\tSignature 13\tSignature 14\tSignature 15\tSignature 16\tSignature 17\tSignature 18\tSignature 19\tSignature 20\tSignature 21\tSignature 22\tSignature 23\tSignature 24\tSignature 25\tSignature 26\tSignature 27\tSignature 28\tSignature 29\tSignature 30\t\t\t\t\t\t\t\r",
"C>A\tACA\tA[C>A]A\t0.011098326166\t0.000682708227\t0.022172306775\t0.036500000000\t0.014941547714\t0.001700000000\t0.000400000000\t0.036718003777\t0.012000000000\t0.000700000000\t0.000200000000\t0.007700000000\t0.000334757187\t0.000100000000\t0.001300000000\t0.016100000000\t0.001832019172\t0.050536418625\t0.010700000000\t0.001179961596\t0.000100000000\t0.001504070375\t0.000453360683\t0.028645992536\t0.009896768001\t0.002039772910\t0.005205626909\t0.001397438842\t0.069981987280\t0.000000000000\t\t\t\t\t\t\t\r",
"C>A\tACC\tA[C>A]C\t0.009149340734\t0.000619107232\t0.017871675353\t0.030900000000\t0.008960917997\t0.002800000000\t0.000500000000\t0.033245722246\t0.006700000000\t0.001000000000\t0.001000000000\t0.004700000000\t0.000648736139\t0.004200000000\t0.004000000000\t0.009700000000\t0.000342235622\t0.010939824789\t0.007400000000\t0.002211505148\t0.000700000000\t0.002451011015\t0.000366800484\t0.020214638433\t0.006998928770\t0.001487162284\t0.004738227415\t0.000917187665\t0.055152357207\t0.000000000000\t\t\t\t\t\t\t\r",
"C>A\tACG\tA[C>A]G\t0.001490070468\t0.000099278956\t0.002138339617\t0.018300000000\t0.002207846012\t0.000500000000\t0.000000000000\t0.002525311256\t0.000500000000\t0.000300000000\t0.000000000000\t0.001700000000\t0.000038144594\t0.000500000000\t0.000000000000\t0.002200000000\t0.000001576225\t0.002288072697\t0.000500000000\t0.000000161691\t0.000000000000\t0.000000000000\t0.000000000000\t0.020478996485\t0.001448442998\t0.000283945614\t0.000782697937\t0.000000000000\t0.017846984044\t0.001967300000\t\t\t\t\t\t\t\r",
"C>A\tACT\tA[C>A]T\t0.006233885236\t0.000323891363\t0.016265145590\t0.024300000000\t0.009206905292\t0.001900000000\t0.000400000000\t0.033598549516\t0.006800000000\t0.009200000000\t0.000200000000\t0.004600000000\t0.000846658460\t0.029600000000\t0.005700000000\t0.008800000000\t0.003179664772\t0.019424091364\t0.007400000000\t0.003008009987\t0.000600000000\t0.000922452518\t0.000000000000\t0.024600145425\t0.004966565004\t0.000597865638\t0.002718242455\t0.000513409971\t0.026804715973\t0.000000000000\t\t\t\t\t\t\t\r",
"C>G\tACA\tA[C>G]A\t0.001801068369\t0.000263481044\t0.024002615202\t0.009700000000\t0.011671021651\t0.001300000000\t0.000000000000\t0.008356844806\t0.004800000000\t0.000500000000\t0.000700000000\t0.003100000000\t0.003775165298\t0.000100000000\t0.001100000000\t0.004800000000\t0.001661704195\t0.001516887580\t0.005800000000\t0.000697041070\t0.000500000000\t0.000527644040\t0.000000000000\t0.011993038821\t0.008032567646\t0.001272881297\t0.001324163246\t0.000255194996\t0.009301971097\t0.000000000000\t\t\t\t\t\t\t\r",
"C>G\tACC\tA[C>G]C\t0.002580908520\t0.000269866049\t0.012160303668\t0.005400000000\t0.007292088583\t0.001200000000\t0.000000000000\t0.004306380304\t0.002300000000\t0.000300000000\t0.000300000000\t0.001500000000\t0.000920824828\t0.000000000000\t0.000100000000\t0.002400000000\t0.001626898399\t0.002498784527\t0.001900000000\t0.002059310006\t0.000800000000\t0.000000000000\t0.000000000000\t0.008425119792\t0.001635970413\t0.001528194862\t0.001771084129\t0.000268898964\t0.003479048871\t0.000000000000\t\t\t\t\t\t\t\r",
"C>G\tACG\tA[C>G]G\t0.000592548022\t0.000219233911\t0.005275419458\t0.003100000000\t0.002303839151\t0.000000000000\t0.000000000000\t0.000584415303\t0.000000000000\t0.000000000000\t0.000000000000\t0.000000000000\t0.000019890489\t0.000000000000\t0.000600000000\t0.000000000000\t0.000025801897\t0.002614508983\t0.000000000000\t0.000012734958\t0.000000000000\t0.000000000000\t0.000000000000\t0.000000000000\t0.000000000000\t0.000307246308\t0.000000000000\t0.000000000000\t0.000154433523\t0.004820000000\t\t\t\t\t\t\t\r",
"C>G\tACT\tA[C>G]T\t0.002963986287\t0.000610973492\t0.023277656342\t0.005400000000\t0.011696245747\t0.001800000000\t0.000100000000\t0.008634906912\t0.003800000000\t0.000200000000\t0.000900000000\t0.002500000000\t0.003860632488\t0.000100000000\t0.001000000000\t0.007300000000\t0.001328507581\t0.003983011170\t0.007200000000\t0.000848587116\t0.001800000000\t0.000299454735\t0.000000000000\t0.003881383185\t0.003428449701\t0.002498251826\t0.001322719819\t0.000307698016\t0.003976789927\t0.000000000000\t\t\t\t\t\t\t\r",
"C>T\tACA\tA[C>T]A\t0.029514532745\t0.007441556817\t0.017872170732\t0.012000000000\t0.021839187574\t0.031200000000\t0.000000000000\t0.018066687156\t0.009300000000\t0.000000000000\t0.022500000000\t0.012100000000\t0.001480180764\t0.029300000000\t0.011700000000\t0.013500000000\t0.009350488403\t0.003825623912\t0.022100000000\t0.032858047272\t0.005100000000\t0.003113389695\t0.019767154372\t0.006315518136\t0.020987618053\t0.005907232039\t0.013723285787\t0.005433792471\t0.005197221326\t0.065119000000\t\t\t\t\t\t\t\r",
"C>T\tACC\tA[C>T]C\t0.014322747041\t0.002726312428\t0.008896034323\t0.007500000000\t0.012756105638\t0.016300000000\t0.019700000000\t0.005650015029\t0.005600000000\t0.003200000000\t0.109900000000\t0.005400000000\t0.000040508707\t0.032100000000\t0.016900000000\t0.007600000000\t0.004224289013\t0.002465214261\t0.020300000000\t0.022196046913\t0.006200000000\t0.001524860290\t0.081880579925\t0.003844326917\t0.013141283697\t0.010626457816\t0.005424306001\t0.001815415488\t0.008976469990\t0.054397000000\t\t\t\t\t\t\t\r",
"C>T\tACG\tA[C>T]G\t0.171646931305\t0.003322083281\t0.003572708350\t0.002800000000\t0.016760177116\t0.090800000000\t0.000100000000\t0.019265105465\t0.012500000000\t0.012600000000\t0.007200000000\t0.001900000000\t0.000316594671\t0.049600000000\t0.030900000000\t0.003500000000\t0.009229213313\t0.008463956886\t0.036800000000\t0.034572014703\t0.025600000000\t0.003625641201\t0.009248715085\t0.003110664425\t0.019659894285\t0.019930323239\t0.006706429454\t0.005432039331\t0.031638996360\t0.020460400000\t\t\t\t\t\t\t\r",
"C>T\tACT\tA[C>T]T\t0.012623763151\t0.003326528417\t0.014797612209\t0.005900000000\t0.016477954692\t0.014900000000\t0.004300000000\t0.020805967497\t0.007600000000\t0.004700000000\t0.063900000000\t0.006700000000\t0.000000000000\t0.027700000000\t0.006600000000\t0.010100000000\t0.005735961173\t0.005605739927\t0.026600000000\t0.018819559611\t0.003000000000\t0.000000000000\t0.033209652209\t0.003420418673\t0.005728519312\t0.011334530552\t0.004475614034\t0.001629213189\t0.004769666594\t0.021935900000\t\t\t\t\t\t\t\r",
"T>A\tATA\tA[T>A]A\t0.004021520333\t0.000000132377\t0.008428563566\t0.004800000000\t0.008902904323\t0.000600000000\t0.001000000000\t0.013365082811\t0.012100000000\t0.000000000000\t0.000200000000\t0.005800000000\t0.001243618346\t0.000000000000\t0.002100000000\t0.008100000000\t0.000613002450\t0.003037210804\t0.001100000000\t0.000574401358\t0.008300000000\t0.049619229712\t0.000000000000\t0.001166747126\t0.020430527125\t0.004459339879\t0.143076389009\t0.000964030891\t0.000000000000\t0.007574300000\t\t\t\t\t\t\t\r",
"T>A\tATC\tA[T>A]C\t0.002371144163\t0.000113069562\t0.007373096670\t0.003900000000\t0.007399202496\t0.003300000000\t0.000800000000\t0.012431202077\t0.004200000000\t0.000200000000\t0.000100000000\t0.006500000000\t0.001060001878\t0.000800000000\t0.009000000000\t0.007300000000\t0.001878587472\t0.000811167214\t0.005300000000\t0.001971462059\t0.029200000000\t0.011666740464\t0.000000000000\t0.000015173057\t0.008892560661\t0.012822242916\t0.001464085387\t0.003365667809\t0.003049595749\t0.003738000000\t\t\t\t\t\t\t\r",
"T>A\tATG\tA[T>A]G\t0.002810909959\t0.000532929416\t0.007357300020\t0.010000000000\t0.011508051457\t0.000000000000\t0.000900000000\t0.014036657067\t0.006800000000\t0.000000000000\t0.000700000000\t0.004600000000\t0.000102554360\t0.000000000000\t0.000900000000\t0.008000000000\t0.001087349685\t0.004503207138\t0.004500000000\t0.000190586616\t0.002700000000\t0.067915454213\t0.000237418521\t0.006332455573\t0.006241948480\t0.001172344272\t0.001226280897\t0.001785753650\t0.003909981984\t0.006590600000\t\t\t\t\t\t\t\r",
"T>A\tATT\tA[T>A]T\t0.008360909345\t0.000149110959\t0.008753938991\t0.003000000000\t0.011193623195\t0.005300000000\t0.003500000000\t0.024109916596\t0.018500000000\t0.001200000000\t0.000300000000\t0.006000000000\t0.000971812099\t0.000900000000\t0.007900000000\t0.008400000000\t0.002045493651\t0.002736523185\t0.003300000000\t0.007421103528\t0.004300000000\t0.014599917580\t0.000000000000\t0.000004967216\t0.014850877622\t0.003993137549\t0.002535197951\t0.005653089100\t0.005105527437\t0.009344700000\t\t\t\t\t\t\t\r",
"T>C\tATA\tA[T>C]A\t0.013915773030\t0.001303509558\t0.013035709432\t0.008400000000\t0.035366735114\t0.007500000000\t0.000500000000\t0.016206510101\t0.025200000000\t0.003100000000\t0.000600000000\t0.046600000000\t0.000549454166\t0.002600000000\t0.011500000000\t0.071700000000\t0.002697314025\t0.002933053644\t0.014500000000\t0.016395062964\t0.026500000000\t0.016839930949\t0.001068448820\t0.004726009142\t0.020663282773\t0.055028652133\t0.005047301081\t0.004092746888\t0.010119181163\t0.000000000000\t\t\t\t\t\t\t\r",
"T>C\tATC\tA[T>C]C\t0.006274960600\t0.000425576677\t0.009186289690\t0.002000000000\t0.013771102129\t0.005600000000\t0.000100000000\t0.007787802820\t0.010500000000\t0.005200000000\t0.001500000000\t0.031700000000\t0.000145506661\t0.004600000000\t0.007300000000\t0.017100000000\t0.008078150072\t0.001513284763\t0.004800000000\t0.019251793639\t0.016500000000\t0.001108262382\t0.000221573159\t0.000627949341\t0.009502478424\t0.027594604208\t0.001845691228\t0.007880729211\t0.006349682655\t0.008754700000\t\t\t\t\t\t\t\r",
"T>C\tATG\tA[T>C]G\t0.010137636154\t0.000575473007\t0.011716920077\t0.008100000000\t0.028449051657\t0.021700000000\t0.000700000000\t0.009852583553\t0.017600000000\t0.002600000000\t0.000800000000\t0.048700000000\t0.000279178024\t0.003500000000\t0.007300000000\t0.034800000000\t0.004432870698\t0.003329668955\t0.013500000000\t0.054784918009\t0.022400000000\t0.004913616405\t0.000498801514\t0.002723758627\t0.019911492456\t0.051791044050\t0.005688813879\t0.005520833604\t0.014659610308\t0.009049800000\t\t\t\t\t\t\t\r",
"T>C\tATT\tA[T>C]T\t0.009256316389\t0.001488164093\t0.016978870798\t0.003600000000\t0.027302974158\t0.002300000000\t0.001100000000\t0.021829463826\t0.024700000000\t0.004000000000\t0.000300000000\t0.025000000000\t0.000267461432\t0.002500000000\t0.006300000000\t0.049700000000\t0.003140370209\t0.005491317771\t0.016000000000\t0.007341269916\t0.018900000000\t0.002630523933\t0.000313240541\t0.007122368470\t0.014488390255\t0.039072244419\t0.002592754577\t0.003232484000\t0.007892413403\t0.000000000000\t\t\t\t\t\t\t\r",
"T>G\tATA\tA[T>G]A\t0.001587636423\t0.000033718411\t0.002351368076\t0.000000000000\t0.003462151975\t0.000000000000\t0.000000000000\t0.004219719284\t0.033500000000\t0.000000000000\t0.000000000000\t0.001200000000\t0.000386199790\t0.000000000000\t0.000000000000\t0.005200000000\t0.000053139249\t0.000386138197\t0.001400000000\t0.000001765291\t0.000000000000\t0.001651998615\t0.000000000000\t0.000000000000\t0.005281679682\t0.000080653426\t0.001136337392\t0.006364675555\t0.000000000000\t0.005311800000\t\t\t\t\t\t\t\r",
"T>G\tATC\tA[T>G]C\t0.001784091288\t0.000024844826\t0.001464231129\t0.000200000000\t0.002246592659\t0.001700000000\t0.000000000000\t0.000793671329\t0.005100000000\t0.004000000000\t0.000000000000\t0.001000000000\t0.000112383214\t0.002000000000\t0.000600000000\t0.004200000000\t0.006468661271\t0.001254982080\t0.000000000000\t0.000381387913\t0.003100000000\t0.000000000000\t0.000000000000\t0.001463950688\t0.000000000000\t0.000163171210\t0.000958976380\t0.016610905647\t0.000000000000\t0.002459200000\t\t\t\t\t\t\t\r",
"T>G\tATG\tA[T>G]G\t0.001385830552\t0.000273454078\t0.009053777642\t0.001500000000\t0.005490216837\t0.000700000000\t0.000900000000\t0.006865800374\t0.005200000000\t0.001500000000\t0.000100000000\t0.001800000000\t0.000036000696\t0.000700000000\t0.001200000000\t0.006100000000\t0.003481549012\t0.002369216073\t0.001000000000\t0.000658758087\t0.001600000000\t0.001086872758\t0.000000000000\t0.000847612530\t0.007136421644\t0.001255079507\t0.000067660610\t0.006592842245\t0.002925171448\t0.006295500000\t\t\t\t\t\t\t\r",
"T>G\tATT\tA[T>G]T\t0.003158539312\t0.000217694872\t0.007030976419\t0.000200000000\t0.003821226712\t0.002900000000\t0.000000000000\t0.002939922259\t0.019300000000\t0.013700000000\t0.000000000000\t0.003000000000\t0.000192493754\t0.005600000000\t0.002900000000\t0.007600000000\t0.030926375169\t0.000338602881\t0.001700000000\t0.000532125222\t0.001800000000\t0.000000000000\t0.000000000000\t0.001075455294\t0.000000000000\t0.001850383492\t0.000256027748\t0.086392897435\t0.000000000000\t0.009541600000\t\t\t\t\t\t\t\r",
"C>A\tCCA\tC[C>A]A\t0.006595870109\t0.000677444988\t0.018781725573\t0.046100000000\t0.009674904330\t0.010100000000\t0.001200000000\t0.031723756590\t0.009800000000\t0.003100000000\t0.000700000000\t0.013500000000\t0.001710089561\t0.005600000000\t0.010600000000\t0.015900000000\t0.001032430218\t0.088768108819\t0.011200000000\t0.017377110647\t0.002000000000\t0.004549692944\t0.000164739382\t0.063559283793\t0.014832947875\t0.003705850083\t0.005065073268\t0.001168515637\t0.051410211701\t0.000000000000\t\t\t\t\t\t\t\r",
"C>A\tCCC\tC[C>A]C\t0.007342367815\t0.000213681013\t0.015760457842\t0.061400000000\t0.004952300642\t0.024100000000\t0.000600000000\t0.025505407135\t0.005700000000\t0.000900000000\t0.001700000000\t0.011200000000\t0.001159256581\t0.010200000000\t0.008400000000\t0.010000000000\t0.000421880127\t0.020641390557\t0.015900000000\t0.036502462984\t0.001400000000\t0.003764473854\t0.000736874802\t0.033757004716\t0.007822175302\t0.003980723373\t0.002234153344\t0.000334291811\t0.025825650831\t0.000000000000\t\t\t\t\t\t\t\r",
"C>A\tCCG\tC[C>A]G\t0.000892840370\t0.000006770460\t0.001963389846\t0.008800000000\t0.002800627307\t0.009100000000\t0.000000000000\t0.001159624303\t0.000000000000\t0.000700000000\t0.001000000000\t0.002800000000\t0.000244166533\t0.000900000000\t0.001500000000\t0.002200000000\t0.000297462769\t0.017178402480\t0.001800000000\t0.012482587542\t0.002700000000\t0.000900163331\t0.000163953661\t0.022428985806\t0.001276976732\t0.000811741970\t0.000266312161\t0.000053651998\t0.014496183324\t0.002262400000\t\t\t\t\t\t\t\r",
"C>A\tCCT\tC[C>A]T\t0.007186581642\t0.000416332906\t0.014722861126\t0.043200000000\t0.011013465824\t0.057100000000\t0.001300000000\t0.028791172971\t0.009100000000\t0.016000000000\t0.001400000000\t0.007100000000\t0.001256768244\t0.125700000000\t0.022800000000\t0.008400000000\t0.000031479429\t0.037676958929\t0.009600000000\t0.103401226249\t0.005600000000\t0.004439846231\t0.000722731835\t0.020086515438\t0.012563654719\t0.019038431345\t0.003100570009\t0.000186671902\t0.040355074099\t0.000000000000\t\t\t\t\t\t\t\r",
"C>G\tCCA\tC[C>G]A\t0.001284983446\t0.000027772077\t0.016832571989\t0.010500000000\t0.007537523174\t0.000000000000\t0.000100000000\t0.006619077216\t0.001800000000\t0.000000000000\t0.000000000000\t0.001900000000\t0.005333651375\t0.000000000000\t0.000900000000\t0.004000000000\t0.001924978751\t0.000458463822\t0.000200000000\t0.000013316260\t0.000100000000\t0.000582017663\t0.000199049174\t0.013476821609\t0.007200790410\t0.001279106069\t0.004730649426\t0.000103366092\t0.002618863169\t0.000000000000\t\t\t\t\t\t\t\r",
"C>G\tCCC\tC[C>G]C\t0.000702134818\t0.000279643862\t0.013531441527\t0.009700000000\t0.007633226665\t0.000000000000\t0.000400000000\t0.006077702895\t0.002000000000\t0.000000000000\t0.000200000000\t0.001900000000\t0.000912303471\t0.000000000000\t0.000300000000\t0.005000000000\t0.000709367774\t0.003654228330\t0.002800000000\t0.000548144708\t0.000000000000\t0.000825299647\t0.000283776026\t0.016105687737\t0.004598473117\t0.001214977303\t0.001261644842\t0.000144310792\t0.001111512270\t0.000000000000\t\t\t\t\t\t\t\r",
"C>G\tCCG\tC[C>G]G\t0.000506289594\t0.000019161576\t0.004176457524\t0.006300000000\t0.002613760399\t0.000000000000\t0.000600000000\t0.000656033433\t0.000000000000\t0.000000000000\t0.000800000000\t0.000000000000\t0.000285467331\t0.000000000000\t0.000400000000\t0.002800000000\t0.000003517869\t0.005405774245\t0.000200000000\t0.000000000000\t0.000000000000\t0.000382114636\t0.000000000000\t0.000797471998\t0.000823709283\t0.000207832588\t0.000360856587\t0.000000000000\t0.001520046950\t0.003246100000\t\t\t\t\t\t\t\r",
"C>G\tCCT\tC[C>G]T\t0.001381542722\t0.000312781573\t0.024046390424\t0.009400000000\t0.009417097876\t0.000200000000\t0.000300000000\t0.007805026957\t0.003900000000\t0.000000000000\t0.000200000000\t0.001900000000\t0.006383097568\t0.000000000000\t0.001700000000\t0.008000000000\t0.000953776334\t0.005519604820\t0.000800000000\t0.000120730556\t0.000300000000\t0.001760645911\t0.000195775338\t0.010388665882\t0.000196360325\t0.002297026211\t0.000268026230\t0.000290999509\t0.000139750691\t0.000000000000\t\t\t\t\t\t\t\r",
"C>T\tCCA\tC[C>T]A\t0.020896446965\t0.015019506931\t0.014395059491\t0.021000000000\t0.022769209355\t0.008500000000\t0.075400000000\t0.004865810021\t0.009800000000\t0.001200000000\t0.031000000000\t0.018700000000\t0.003647278451\t0.006000000000\t0.005100000000\t0.021000000000\t0.021265205902\t0.009767511475\t0.043800000000\t0.011145185306\t0.003400000000\t0.007731349522\t0.060572627805\t0.006035447394\t0.014333994691\t0.006594444885\t0.017438484778\t0.004391532554\t0.009126471767\t0.069447200000\t\t\t\t\t\t\t\r",
"C>T\tCCC\tC[C>T]C\t0.018501704770\t0.003516918098\t0.008544781192\t0.014400000000\t0.017509476305\t0.009900000000\t0.100700000000\t0.003980150431\t0.006900000000\t0.002400000000\t0.151800000000\t0.009400000000\t0.001004229295\t0.016700000000\t0.008800000000\t0.009800000000\t0.003847235312\t0.005654828462\t0.075300000000\t0.018152746866\t0.009100000000\t0.007110850526\t0.186375634076\t0.000000000000\t0.005936643659\t0.006510631971\t0.007144870207\t0.002647660966\t0.007067603037\t0.063840300000\t\t\t\t\t\t\t\r",
"C>T\tCCG\tC[C>T]G\t0.095577217268\t0.004979275715\t0.003518465840\t0.007600000000\t0.012862237589\t0.090100000000\t0.020800000000\t0.008338639182\t0.007600000000\t0.010900000000\t0.013500000000\t0.004600000000\t0.000129202148\t0.029300000000\t0.015500000000\t0.001600000000\t0.013176112314\t0.016183222964\t0.023800000000\t0.044384942487\t0.031800000000\t0.005502030965\t0.015149082593\t0.002924419074\t0.011559391627\t0.011904935324\t0.007503561655\t0.002284563551\t0.019875180409\t0.017312600000\t\t\t\t\t\t\t\r",
"C>T\tCCT\tC[C>T]T\t0.017113307642\t0.008956552775\t0.016075545712\t0.020100000000\t0.020429513385\t0.008700000000\t0.078800000000\t0.018843627434\t0.009700000000\t0.004300000000\t0.081600000000\t0.020500000000\t0.000611769970\t0.018600000000\t0.009600000000\t0.021100000000\t0.009648025316\t0.005153318613\t0.115900000000\t0.013847136214\t0.004000000000\t0.004169777241\t0.156015396977\t0.003971201575\t0.013632596885\t0.006238556327\t0.007148027702\t0.002423984267\t0.001634058687\t0.034231800000\t\t\t\t\t\t\t\r",
"T>A\tCTA\tC[T>A]A\t0.001182587416\t0.000154598981\t0.007571320260\t0.007500000000\t0.005043163511\t0.000100000000\t0.000000000000\t0.012024217522\t0.005500000000\t0.000000000000\t0.000400000000\t0.001700000000\t0.000285158991\t0.000000000000\t0.001600000000\t0.005400000000\t0.000030079241\t0.000412346353\t0.004300000000\t0.000045042364\t0.000000000000\t0.086898545516\t0.000594528453\t0.002766634838\t0.026296323089\t0.003561422681\t0.330601017471\t0.002566579967\t0.001765447503\t0.005311800000\t\t\t\t\t\t\t\r",
"T>A\tCTC\tC[T>A]C\t0.001903166857\t0.000464019028\t0.012725463528\t0.011100000000\t0.006208843676\t0.002600000000\t0.001400000000\t0.017881092423\t0.004900000000\t0.000100000000\t0.001100000000\t0.008600000000\t0.000222287328\t0.003400000000\t0.004300000000\t0.012800000000\t0.010223747696\t0.002491411310\t0.005300000000\t0.004389055645\t0.003700000000\t0.053706846809\t0.000000000000\t0.002694462454\t0.015799883085\t0.003901610981\t0.002813598808\t0.019610897429\t0.000383359117\t0.006787300000\t\t\t\t\t\t\t\r",
"T>A\tCTG\tC[T>A]G\t0.001487960630\t0.000230409809\t0.011508561865\t0.034200000000\t0.010507568970\t0.000800000000\t0.000700000000\t0.016356774533\t0.005100000000\t0.000300000000\t0.000600000000\t0.003100000000\t0.000286951828\t0.000200000000\t0.002700000000\t0.013000000000\t0.004467036232\t0.005523045713\t0.002100000000\t0.001127835842\t0.002500000000\t0.213220866558\t0.000173120400\t0.006911164721\t0.062860606381\t0.002390304180\t0.004396857083\t0.008833600032\t0.003556845215\t0.008262800000\t\t\t\t\t\t\t\r",
"T>A\tCTT\tC[T>A]T\t0.002179344412\t0.000574885587\t0.016456177696\t0.011500000000\t0.009784640746\t0.001100000000\t0.001800000000\t0.026247801917\t0.006900000000\t0.000900000000\t0.000400000000\t0.005400000000\t0.000000800446\t0.002000000000\t0.000700000000\t0.009600000000\t0.031244883388\t0.003313713273\t0.002200000000\t0.002385541321\t0.001000000000\t0.061827507476\t0.000457813099\t0.002915200165\t0.036622351112\t0.001635922667\t0.005430891633\t0.010942609820\t0.002232755492\t0.000000000000\t\t\t\t\t\t\t\r",
"T>C\tCTA\tC[T>C]A\t0.004176674882\t0.000546619204\t0.007895012611\t0.005200000000\t0.014238629786\t0.006200000000\t0.001000000000\t0.005173842417\t0.015400000000\t0.000500000000\t0.000100000000\t0.031200000000\t0.000042751466\t0.000600000000\t0.003200000000\t0.014600000000\t0.007426999646\t0.000192169320\t0.006600000000\t0.014867696528\t0.021900000000\t0.004027346334\t0.000000000000\t0.001848876535\t0.009028014894\t0.037889377691\t0.004899079238\t0.020738003352\t0.000000000000\t0.000000000000\t\t\t\t\t\t\t\r",
"T>C\tCTC\tC[T>C]C\t0.005252593331\t0.000392332573\t0.014431170362\t0.002600000000\t0.012489500440\t0.004000000000\t0.000800000000\t0.010936893768\t0.011800000000\t0.001900000000\t0.000400000000\t0.038500000000\t0.000043307542\t0.005300000000\t0.005000000000\t0.019400000000\t0.039657001207\t0.003050192581\t0.007000000000\t0.018356161025\t0.010000000000\t0.001352443929\t0.000000000000\t0.004453689857\t0.010771632420\t0.021740734953\t0.005724448467\t0.043283797246\t0.002774953681\t0.000000000000\t\t\t\t\t\t\t\r",
"T>C\tCTG\tC[T>C]G\t0.007013225310\t0.000362024508\t0.008422972144\t0.010000000000\t0.018741631229\t0.027000000000\t0.001000000000\t0.005657773602\t0.015400000000\t0.001500000000\t0.000100000000\t0.041500000000\t0.000054272534\t0.001600000000\t0.012800000000\t0.019000000000\t0.028609696263\t0.002596464193\t0.004000000000\t0.074150562563\t0.021100000000\t0.002749666137\t0.000320704886\t0.000000000000\t0.013350661310\t0.047240103343\t0.005661478992\t0.037088261728\t0.006902414709\t0.000000000000\t\t\t\t\t\t\t\r",
"T>C\tCTT\tC[T>C]T\t0.006713813119\t0.000560900126\t0.011932430465\t0.005400000000\t0.019054082944\t0.003300000000\t0.004500000000\t0.011120389322\t0.027500000000\t0.002000000000\t0.000900000000\t0.032400000000\t0.000000763944\t0.002700000000\t0.003500000000\t0.026300000000\t0.101869657701\t0.001310028744\t0.005300000000\t0.011001118253\t0.020400000000\t0.002363353537\t0.000000000000\t0.001597389007\t0.008961029966\t0.020740621431\t0.003900589062\t0.005052716018\t0.003041424895\t0.000000000000\t\t\t\t\t\t\t\r",
"T>G\tCTA\tC[T>G]A\t0.000302691186\t0.000113784603\t0.001974181146\t0.000000000000\t0.002041648941\t0.000100000000\t0.000000000000\t0.002424050552\t0.021600000000\t0.002200000000\t0.000000000000\t0.001300000000\t0.000015094208\t0.000000000000\t0.000000000000\t0.000200000000\t0.000091469787\t0.000061002548\t0.000000000000\t0.000023264084\t0.000000000000\t0.001120356515\t0.000000000000\t0.000000000000\t0.000000000000\t0.000095396260\t0.000000000000\t0.009589841941\t0.001116854814\t0.004623300000\t\t\t\t\t\t\t\r",
"T>G\tCTC\tC[T>G]C\t0.002098502440\t0.000022095087\t0.005824295545\t0.001300000000\t0.003478502118\t0.004000000000\t0.000800000000\t0.002096210616\t0.006400000000\t0.001800000000\t0.000100000000\t0.004900000000\t0.000147624840\t0.002000000000\t0.001900000000\t0.008200000000\t0.019814293551\t0.002015696913\t0.001100000000\t0.008476494525\t0.000500000000\t0.001064583571\t0.000190144342\t0.002185437348\t0.002886454337\t0.003919455100\t0.000381244984\t0.033913223722\t0.000000000000\t0.006000400000\t\t\t\t\t\t\t\r",
"T>G\tCTG\tC[T>G]G\t0.001599548500\t0.000228245948\t0.010464654525\t0.004600000000\t0.007147456226\t0.005000000000\t0.000900000000\t0.006604046294\t0.012600000000\t0.003700000000\t0.001100000000\t0.004500000000\t0.000045987300\t0.000500000000\t0.000700000000\t0.006700000000\t0.013449837633\t0.005311967254\t0.001300000000\t0.015383776115\t0.000400000000\t0.004188568126\t0.000000000000\t0.000113909674\t0.011769492499\t0.005497357314\t0.000403257236\t0.018516481606\t0.002169693346\t0.007377500000\t\t\t\t\t\t\t\r",
"T>G\tCTT\tC[T>G]T\t0.002758537619\t0.000067111340\t0.008724387302\t0.001200000000\t0.011486811509\t0.008600000000\t0.001300000000\t0.004866713859\t0.050900000000\t0.018200000000\t0.000900000000\t0.006300000000\t0.000463714722\t0.006300000000\t0.004500000000\t0.018600000000\t0.261456614079\t0.002341629178\t0.001900000000\t0.002185394659\t0.003200000000\t0.001629409573\t0.000210834980\t0.002533718333\t0.008550807482\t0.006788718691\t0.001439096069\t0.118967102954\t0.001143663251\t0.009148100000\t\t\t\t\t\t\t\r",
"C>A\tGCA\tG[C>A]A\t0.008232603989\t0.000352013353\t0.009696539676\t0.037600000000\t0.011892168979\t0.002400000000\t0.000300000000\t0.023682328869\t0.011800000000\t0.001400000000\t0.000400000000\t0.006200000000\t0.000132109643\t0.001800000000\t0.002400000000\t0.009600000000\t0.006535404875\t0.128724158080\t0.003200000000\t0.001116123789\t0.000100000000\t0.001298370161\t0.000349907495\t0.054676448679\t0.013465295053\t0.001375311846\t0.010755871858\t0.002136629077\t0.078046610108\t0.008853000000\t\t\t\t\t\t\t\r",
"C>A\tGCC\tG[C>A]C\t0.005758021410\t0.000133816879\t0.010843341116\t0.039900000000\t0.009247857495\t0.005800000000\t0.000100000000\t0.015821896445\t0.009200000000\t0.002200000000\t0.001000000000\t0.005600000000\t0.000754244011\t0.011400000000\t0.009900000000\t0.009400000000\t0.001293803953\t0.016092878663\t0.003100000000\t0.003271428605\t0.003800000000\t0.001809622152\t0.000568337772\t0.109945988997\t0.006468298058\t0.001961639088\t0.012442335117\t0.000429208379\t0.073460914669\t0.009344900000\t\t\t\t\t\t\t\r",
"C>A\tGCG\tG[C>A]G\t0.000616335232\t0.000178441660\t0.000929140538\t0.022700000000\t0.002809188466\t0.002100000000\t0.000000000000\t0.000850913754\t0.000000000000\t0.000200000000\t0.000000000000\t0.001500000000\t0.000096820775\t0.001100000000\t0.001300000000\t0.003600000000\t0.001283663378\t0.009271416762\t0.000000000000\t0.003751612508\t0.000300000000\t0.000333838055\t0.000149417834\t0.046524869010\t0.000658103089\t0.000013229664\t0.001629628347\t0.000000000000\t0.017383296915\t0.000885300000\t\t\t\t\t\t\t\r",
"C>A\tGCT\tG[C>A]T\t0.004459080311\t0.000122832046\t0.012215382639\t0.025800000000\t0.010301267456\t0.008700000000\t0.000100000000\t0.021061235802\t0.008500000000\t0.008800000000\t0.000600000000\t0.002700000000\t0.000551103672\t0.073600000000\t0.030900000000\t0.006300000000\t0.003103766672\t0.072017070126\t0.007200000000\t0.012937194450\t0.002300000000\t0.000519148161\t0.000000000000\t0.050428709482\t0.010473080347\t0.001934996293\t0.007412355152\t0.001318507256\t0.059602472118\t0.008164500000\t\t\t\t\t\t\t\r",
"C>G\tGCA\tG[C>G]A\t0.000602122711\t0.000044838514\t0.011916844569\t0.007000000000\t0.005559423089\t0.000000000000\t0.000000000000\t0.003724567946\t0.001100000000\t0.000000000000\t0.000000000000\t0.001200000000\t0.001349419584\t0.000200000000\t0.002200000000\t0.003100000000\t0.001296324191\t0.000923879703\t0.002000000000\t0.001339769321\t0.000400000000\t0.000000000000\t0.000000000000\t0.004945060708\t0.004169655524\t0.001320607475\t0.001418166387\t0.000081591411\t0.000265759858\t0.007377500000\t\t\t\t\t\t\t\r",
"C>G\tGCC\tG[C>G]C\t0.002393352209\t0.000014520103\t0.009823653397\t0.009100000000\t0.005388801221\t0.003000000000\t0.000400000000\t0.003046596953\t0.002900000000\t0.000200000000\t0.000200000000\t0.001100000000\t0.000788807743\t0.000100000000\t0.000500000000\t0.003200000000\t0.000260760911\t0.000648647922\t0.000000000000\t0.001811151475\t0.004200000000\t0.000000000000\t0.000137501074\t0.009369875774\t0.002392804944\t0.001845501071\t0.001346716782\t0.000105219031\t0.000945012884\t0.006590600000\t\t\t\t\t\t\t\r",
"C>G\tGCG\tG[C>G]G\t0.000000248534\t0.000040665889\t0.001671054250\t0.006200000000\t0.001100589633\t0.000000000000\t0.000000000000\t0.000321212793\t0.000000000000\t0.000000000000\t0.000000000000\t0.000000000000\t0.000081791214\t0.000000000000\t0.000600000000\t0.000000000000\t0.000004377552\t0.000542839968\t0.000000000000\t0.000000637538\t0.000000000000\t0.000214695851\t0.000000000000\t0.005625726540\t0.000418193916\t0.000205060777\t0.000000000000\t0.000000000000\t0.000000000000\t0.001573900000\t\t\t\t\t\t\t\r",
"C>G\tGCT\tG[C>G]T\t0.000890080731\t0.000268429459\t0.017914334306\t0.006000000000\t0.006041247837\t0.001700000000\t0.000000000000\t0.005775663598\t0.004400000000\t0.000400000000\t0.000000000000\t0.000900000000\t0.002680123244\t0.000300000000\t0.001800000000\t0.002900000000\t0.001820392786\t0.001376724888\t0.000000000000\t0.001238264182\t0.000400000000\t0.000538438803\t0.000000000000\t0.003055770064\t0.001994012217\t0.001225532856\t0.004560685974\t0.000159679354\t0.004796744744\t0.000000000000\t\t\t\t\t\t\t\r",
"C>T\tGCA\tG[C>T]A\t0.024943814154\t0.006390807876\t0.016127273332\t0.008700000000\t0.020038258645\t0.065300000000\t0.000000000000\t0.002261262686\t0.006200000000\t0.000000000000\t0.026100000000\t0.012400000000\t0.002934137577\t0.055700000000\t0.047600000000\t0.012500000000\t0.002068754725\t0.002703068727\t0.037400000000\t0.056793947985\t0.016300000000\t0.001915370854\t0.053282058965\t0.023679051295\t0.010038112013\t0.009606515027\t0.014215313743\t0.003384240869\t0.011638793032\t0.048593400000\t\t\t\t\t\t\t\r",
"C>T\tGCC\tG[C>T]C\t0.027161494035\t0.001995816784\t0.008208632437\t0.008000000000\t0.018022312238\t0.077300000000\t0.021000000000\t0.001616649187\t0.006900000000\t0.013400000000\t0.097500000000\t0.009200000000\t0.002515171146\t0.061800000000\t0.134500000000\t0.009200000000\t0.000325483399\t0.004199307536\t0.061400000000\t0.054429857270\t0.039200000000\t0.003229133547\t0.162918343743\t0.015919201840\t0.006953875155\t0.019507426720\t0.020143646182\t0.001687414029\t0.003781751661\t0.049478700000\t\t\t\t\t\t\t\r",
"C>T\tGCG\tG[C>T]G\t0.103570762296\t0.000303020977\t0.001212922866\t0.002300000000\t0.013193817648\t0.133900000000\t0.000200000000\t0.011606711600\t0.008800000000\t0.027000000000\t0.009000000000\t0.002300000000\t0.003205264550\t0.097800000000\t0.182900000000\t0.000900000000\t0.007511101358\t0.011897360825\t0.027800000000\t0.004592688609\t0.001900000000\t0.003300598879\t0.013046494744\t0.012648077606\t0.010468363027\t0.022502717383\t0.009368468496\t0.002189308960\t0.019684388045\t0.015738700000\t\t\t\t\t\t\t\r",
"C>T\tGCT\tG[C>T]T\t0.017689854381\t0.003265818188\t0.010611649208\t0.008200000000\t0.019502915335\t0.052400000000\t0.016100000000\t0.006284928996\t0.010100000000\t0.015200000000\t0.052200000000\t0.011500000000\t0.001687730545\t0.055500000000\t0.080000000000\t0.011600000000\t0.000925697167\t0.001913380905\t0.091700000000\t0.049580863753\t0.016600000000\t0.001114059570\t0.107590399754\t0.018129778378\t0.009823820735\t0.017307472184\t0.007589625950\t0.001043993631\t0.005761336823\t0.018886500000\t\t\t\t\t\t\t\r",
"T>A\tGTA\tG[T>A]A\t0.000689289439\t0.000114994210\t0.004434782559\t0.006900000000\t0.006740387717\t0.000000000000\t0.000000000000\t0.008139178608\t0.003000000000\t0.000000000000\t0.000300000000\t0.001500000000\t0.000137960405\t0.000000000000\t0.000000000000\t0.003000000000\t0.000008545689\t0.001875434081\t0.002600000000\t0.000163129836\t0.000000000000\t0.053709145693\t0.000000000000\t0.002225127336\t0.036217800837\t0.002242886513\t0.109851871980\t0.000551460783\t0.002930098053\t0.003738000000\t\t\t\t\t\t\t\r",
"T>A\tGTC\tG[T>A]C\t0.000552409528\t0.000293862060\t0.005615432411\t0.005200000000\t0.004021596480\t0.002800000000\t0.000000000000\t0.007910080318\t0.002900000000\t0.000100000000\t0.000300000000\t0.003900000000\t0.000252867405\t0.002500000000\t0.003800000000\t0.004000000000\t0.002573466959\t0.001510309757\t0.000000000000\t0.003799179503\t0.014200000000\t0.009489396694\t0.000183203812\t0.004853149157\t0.005847109647\t0.005207287095\t0.000023455678\t0.003163987597\t0.000000000000\t0.003344500000\t\t\t\t\t\t\t\r",
"T>A\tGTG\tG[T>A]G\t0.001200228847\t0.000088721344\t0.008070136472\t0.013300000000\t0.006814352288\t0.000600000000\t0.000600000000\t0.009057263056\t0.002800000000\t0.000000000000\t0.000400000000\t0.001900000000\t0.000153517461\t0.000000000000\t0.002600000000\t0.005100000000\t0.003538822593\t0.001505368776\t0.001600000000\t0.000828333700\t0.003100000000\t0.050293902428\t0.000000000000\t0.000195758628\t0.038241179315\t0.001358121778\t0.001565937159\t0.002152847013\t0.000533632401\t0.005016700000\t\t\t\t\t\t\t\r",
"T>A\tGTT\tG[T>A]T\t0.002107136837\t0.000215755466\t0.008679121744\t0.004500000000\t0.005101032909\t0.002100000000\t0.001000000000\t0.017812238647\t0.003600000000\t0.000400000000\t0.000100000000\t0.003100000000\t0.000122745552\t0.000900000000\t0.004300000000\t0.003300000000\t0.003879009961\t0.003611025202\t0.005200000000\t0.001900746587\t0.004100000000\t0.010601457233\t0.000000000000\t0.000000000000\t0.020167954419\t0.002513112942\t0.000000000000\t0.002007733244\t0.000000000000\t0.006393900000\t\t\t\t\t\t\t\r",
"T>C\tGTA\tG[T>C]A\t0.011247835116\t0.000004673852\t0.006850055574\t0.006100000000\t0.016294982792\t0.012200000000\t0.000200000000\t0.005720778549\t0.009100000000\t0.006400000000\t0.001000000000\t0.057300000000\t0.001054081222\t0.007800000000\t0.019500000000\t0.016100000000\t0.000011401335\t0.000529126468\t0.003900000000\t0.002469317074\t0.142600000000\t0.004217454113\t0.000251168628\t0.001269749840\t0.017489823879\t0.098052824719\t0.004777380354\t0.003655048127\t0.007090803763\t0.000000000000\t\t\t\t\t\t\t\r",
"T>C\tGTC\tG[T>C]C\t0.006999724257\t0.000186065927\t0.006260749014\t0.001600000000\t0.009574560593\t0.005900000000\t0.000100000000\t0.005181051502\t0.009700000000\t0.010100000000\t0.001600000000\t0.027400000000\t0.000529899747\t0.012700000000\t0.049900000000\t0.009900000000\t0.009657713318\t0.002279130829\t0.003000000000\t0.001565632906\t0.066500000000\t0.000000000000\t0.000538218489\t0.006710339996\t0.004045612446\t0.040226363963\t0.001964954330\t0.007071175327\t0.007111327534\t0.006590600000\t\t\t\t\t\t\t\r",
"T>C\tGTG\tG[T>C]G\t0.004977592617\t0.000000495044\t0.006098762993\t0.004200000000\t0.014103929348\t0.011500000000\t0.000300000000\t0.005040707471\t0.009400000000\t0.005400000000\t0.001900000000\t0.041000000000\t0.000630095718\t0.008200000000\t0.015400000000\t0.014800000000\t0.004654129131\t0.005871056581\t0.006000000000\t0.027341755021\t0.088400000000\t0.001474884534\t0.000000000000\t0.000277684976\t0.013402246274\t0.044620700983\t0.001229167749\t0.008232238440\t0.005643063554\t0.007869400000\t\t\t\t\t\t\t\r",
"T>C\tGTT\tG[T>C]T\t0.010667406133\t0.000578598104\t0.007509343859\t0.002400000000\t0.015666874692\t0.004200000000\t0.008100000000\t0.008488220258\t0.012600000000\t0.012400000000\t0.002300000000\t0.040500000000\t0.000676187888\t0.008600000000\t0.009400000000\t0.014200000000\t0.008149234601\t0.004771565739\t0.009900000000\t0.001498264796\t0.100300000000\t0.001561042738\t0.000995769681\t0.004238477669\t0.014159693736\t0.055528068280\t0.003747225013\t0.000407884510\t0.011237674345\t0.008459600000\t\t\t\t\t\t\t\r",
"T>G\tGTA\tG[T>G]A\t0.000099045003\t0.000095552392\t0.004144488018\t0.000000000000\t0.001627664547\t0.000000000000\t0.000000000000\t0.000974578665\t0.007200000000\t0.000000000000\t0.000000000000\t0.000000000000\t0.000018489857\t0.000000000000\t0.000000000000\t0.000000000000\t0.000039193821\t0.001341532464\t0.000000000000\t0.000003579149\t0.000000000000\t0.000000000000\t0.000000000000\t0.000000000000\t0.005546226924\t0.000037393232\t0.000000000000\t0.000000000000\t0.000000000000\t0.003246100000\t\t\t\t\t\t\t\r",
"T>G\tGTC\tG[T>G]C\t0.000202365646\t0.000047002381\t0.004501985313\t0.000000000000\t0.000327734937\t0.001600000000\t0.000000000000\t0.000524821555\t0.000600000000\t0.002000000000\t0.000000000000\t0.000400000000\t0.000093373513\t0.001900000000\t0.002200000000\t0.003200000000\t0.009078461515\t0.000000920431\t0.000000000000\t0.002830575073\t0.000600000000\t0.000000000000\t0.000000000000\t0.001073632302\t0.000000000000\t0.002460704995\t0.000050790565\t0.000000000000\t0.000000000000\t0.001869000000\t\t\t\t\t\t\t\r",
"T>G\tGTG\tG[T>G]G\t0.001188353185\t0.000109925703\t0.016391452608\t0.001800000000\t0.005948879802\t0.001000000000\t0.001700000000\t0.006087753467\t0.005000000000\t0.000900000000\t0.001000000000\t0.001100000000\t0.000010579194\t0.001200000000\t0.003700000000\t0.000800000000\t0.004782775461\t0.009694999981\t0.004300000000\t0.002733219405\t0.000400000000\t0.001877589181\t0.000187132414\t0.002924696996\t0.008685244444\t0.000817201563\t0.004742918550\t0.004217808299\t0.003640593813\t0.003344500000\t\t\t\t\t\t\t\r",
"T>G\tGTT\tG[T>G]T\t0.000800723342\t0.000086477180\t0.007067236643\t0.000200000000\t0.003307466591\t0.003500000000\t0.000900000000\t0.005427433785\t0.018500000000\t0.003000000000\t0.000000000000\t0.003200000000\t0.000053741207\t0.003800000000\t0.019000000000\t0.004400000000\t0.063497756582\t0.002890553519\t0.002500000000\t0.003558558565\t0.000800000000\t0.000000000000\t0.000000000000\t0.000825098263\t0.002768571683\t0.007833561271\t0.002298476030\t0.031648997321\t0.006181023776\t0.005606900000\t\t\t\t\t\t\t\r",
"C>A\tTCA\tT[C>A]A\t0.012250063666\t0.015127427519\t0.011653224909\t0.033000000000\t0.014774016248\t0.001700000000\t0.001000000000\t0.027032383039\t0.022200000000\t0.037400000000\t0.000900000000\t0.006600000000\t0.048166047937\t0.000500000000\t0.005700000000\t0.012200000000\t0.001824558257\t0.067080643529\t0.005100000000\t0.000688299692\t0.000100000000\t0.004051634552\t0.002071289985\t0.031451263182\t0.006938146585\t0.002680062424\t0.007079284522\t0.000775902768\t0.026894406511\t0.000000000000\t\t\t\t\t\t\t\r",
"C>A\tTCC\tT[C>A]C\t0.011162229329\t0.006532492454\t0.016606775011\t0.053800000000\t0.012043464964\t0.002900000000\t0.002000000000\t0.018089773290\t0.004300000000\t0.010300000000\t0.002500000000\t0.009900000000\t0.017329557089\t0.016600000000\t0.006200000000\t0.014500000000\t0.001286518338\t0.041869009693\t0.006300000000\t0.004058229060\t0.000500000000\t0.002739770937\t0.000696410201\t0.062538293699\t0.010246607400\t0.002031643866\t0.006268800628\t0.000829292467\t0.031608136996\t0.000000000000\t\t\t\t\t\t\t\r",
"C>A\tTCG\tT[C>A]G\t0.002275495686\t0.001656455398\t0.001357239439\t0.010400000000\t0.003902362388\t0.001100000000\t0.000200000000\t0.001694875197\t0.000000000000\t0.003100000000\t0.000100000000\t0.001900000000\t0.002293157306\t0.000000000000\t0.002500000000\t0.004000000000\t0.002172941968\t0.015162379290\t0.001800000000\t0.001001102231\t0.000000000000\t0.000513151070\t0.000099393633\t0.012802590728\t0.000000000000\t0.000265071700\t0.001446403415\t0.000000000000\t0.012500468152\t0.001672200000\t\t\t\t\t\t\t\r",
"C>A\tTCT\tT[C>A]T\t0.015259102491\t0.012394610714\t0.016328075975\t0.037000000000\t0.018243395705\t0.005800000000\t0.001300000000\t0.038141330868\t0.032200000000\t0.308300000000\t0.000400000000\t0.004900000000\t0.018634581512\t0.051600000000\t0.014500000000\t0.015700000000\t0.005012316394\t0.117877846666\t0.008700000000\t0.012229104507\t0.000400000000\t0.003605650899\t0.000207430191\t0.030297679613\t0.024625787226\t0.003017187457\t0.009806277750\t0.001998073667\t0.041207176390\t0.000000000000\t\t\t\t\t\t\t\r",
"C>G\tTCA\tT[C>G]A\t0.001874853199\t0.037241853149\t0.016041405438\t0.003200000000\t0.002681056961\t0.000000000000\t0.000200000000\t0.004117878837\t0.003300000000\t0.000000000000\t0.000100000000\t0.002400000000\t0.280234946548\t0.000000000000\t0.000900000000\t0.005000000000\t0.000392684657\t0.001802969051\t0.000500000000\t0.000024298469\t0.000200000000\t0.000225490615\t0.000230608944\t0.010220857436\t0.009798724567\t0.004201896751\t0.002870433721\t0.000521370368\t0.000820444889\t0.000000000000\t\t\t\t\t\t\t\r",
"C>G\tTCC\tT[C>G]C\t0.002067418791\t0.000019413411\t0.020149920428\t0.010500000000\t0.007924048307\t0.000200000000\t0.001000000000\t0.003799649592\t0.002500000000\t0.000000000000\t0.001100000000\t0.002700000000\t0.063885286587\t0.000000000000\t0.002600000000\t0.006400000000\t0.000655602418\t0.007863006794\t0.001100000000\t0.000151998550\t0.000300000000\t0.000628695020\t0.000165787009\t0.007578365252\t0.003266507150\t0.002807761534\t0.003107967569\t0.000667325890\t0.000444119287\t0.000000000000\t\t\t\t\t\t\t\r",
"C>G\tTCG\tT[C>G]G\t0.000304897004\t0.001625465457\t0.002527910936\t0.003100000000\t0.001319075649\t0.000000000000\t0.000300000000\t0.000025653018\t0.000000000000\t0.000000000000\t0.000300000000\t0.000000000000\t0.009528486008\t0.000000000000\t0.000000000000\t0.000000000000\t0.000013641541\t0.001513223334\t0.000000000000\t0.000000000000\t0.000000000000\t0.000395807994\t0.000000000000\t0.001592714217\t0.002057436753\t0.000000006197\t0.000012810409\t0.000000000000\t0.000143297155\t0.001967300000\t\t\t\t\t\t\t\r",
"C>G\tTCT\tT[C>G]T\t0.003151574460\t0.066879899684\t0.032673640762\t0.005000000000\t0.006644595708\t0.000100000000\t0.000700000000\t0.006247553738\t0.004900000000\t0.000100000000\t0.000800000000\t0.001500000000\t0.330225481745\t0.000100000000\t0.001500000000\t0.016900000000\t0.002852380395\t0.003567102065\t0.010500000000\t0.000683059976\t0.000300000000\t0.000257375147\t0.000438955809\t0.006100727190\t0.016400916850\t0.008019382986\t0.002130587503\t0.001245779535\t0.001369869047\t0.000000000000\t\t\t\t\t\t\t\r",
"C>T\tTCA\tT[C>T]A\t0.014492099634\t0.419941399594\t0.008880209061\t0.003500000000\t0.010997560324\t0.007400000000\t0.120200000000\t0.007180668540\t0.005000000000\t0.003700000000\t0.030200000000\t0.013800000000\t0.113842093163\t0.009200000000\t0.011200000000\t0.017200000000\t0.016100784572\t0.014923391701\t0.021800000000\t0.010710322070\t0.006400000000\t0.002016721688\t0.014558613696\t0.026571894346\t0.011077520800\t0.011303288768\t0.018757776459\t0.007622257672\t0.012076583098\t0.084989200000\t\t\t\t\t\t\t\r",
"C>T\tTCC\tT[C>T]C\t0.017680775357\t0.081972496144\t0.013529572798\t0.007000000000\t0.020644964803\t0.006700000000\t0.288700000000\t0.004061032938\t0.008400000000\t0.021100000000\t0.158900000000\t0.009500000000\t0.015024839001\t0.028300000000\t0.015600000000\t0.013100000000\t0.016779984064\t0.009962744310\t0.033000000000\t0.008175509770\t0.008600000000\t0.001990834246\t0.043510446939\t0.012055691082\t0.006448787540\t0.010807526975\t0.007989184406\t0.005116619241\t0.000000000000\t0.090301000000\t\t\t\t\t\t\t\r",
"C>T\tTCG\tT[C>T]G\t0.076002221712\t0.047720186368\t0.001705408858\t0.001100000000\t0.007534489574\t0.039100000000\t0.099200000000\t0.005535121520\t0.004700000000\t0.214100000000\t0.008000000000\t0.001000000000\t0.006102097956\t0.009400000000\t0.010400000000\t0.003600000000\t0.016272095596\t0.011549715126\t0.009000000000\t0.023163899202\t0.012600000000\t0.001498473091\t0.002752640544\t0.005363359680\t0.014083416181\t0.010364329328\t0.004609130971\t0.001918028110\t0.015830196802\t0.015148500000\t\t\t\t\t\t\t\r",
"C>T\tTCT\tT[C>T]T\t0.013761704021\t0.228674918306\t0.010304416928\t0.007700000000\t0.011787461775\t0.004700000000\t0.084400000000\t0.012209070497\t0.009600000000\t0.039200000000\t0.095400000000\t0.009000000000\t0.028140661233\t0.022900000000\t0.007600000000\t0.019300000000\t0.019827417550\t0.005398013994\t0.026700000000\t0.009540087238\t0.006400000000\t0.001182226502\t0.024328785234\t0.005532499180\t0.012418585885\t0.007249145360\t0.011880932482\t0.006193158626\t0.013462137300\t0.045544000000\t\t\t\t\t\t\t\r",
"T>A\tTTA\tT[T>A]A\t0.005600155423\t0.000080752486\t0.007132681354\t0.004500000000\t0.009205637343\t0.000200000000\t0.001000000000\t0.015139996236\t0.021500000000\t0.000000000000\t0.000300000000\t0.005700000000\t0.000403261316\t0.000000000000\t0.000800000000\t0.008500000000\t0.000698812311\t0.004663278043\t0.003900000000\t0.000924547968\t0.002700000000\t0.040810502825\t0.000000000000\t0.000000000000\t0.021830329049\t0.000628145256\t0.052813887494\t0.000203237563\t0.001191207434\t0.007672600000\t\t\t\t\t\t\t\r",
"T>A\tTTC\tT[T>A]C\t0.001999079260\t0.000004802898\t0.009102797827\t0.004600000000\t0.006835344247\t0.000800000000\t0.001500000000\t0.012093554606\t0.002700000000\t0.000800000000\t0.000000000000\t0.013200000000\t0.000771813274\t0.001600000000\t0.001900000000\t0.008700000000\t0.002579358604\t0.002578155543\t0.004800000000\t0.000708206814\t0.004500000000\t0.025570395801\t0.000000000000\t0.007994129155\t0.003405686728\t0.005073753861\t0.000047091785\t0.010180964868\t0.001682664015\t0.006492200000\t\t\t\t\t\t\t\r",
"T>A\tTTG\tT[T>A]G\t0.001090065693\t0.000066605043\t0.006565863041\t0.008200000000\t0.007144208769\t0.000000000000\t0.001000000000\t0.008362546590\t0.014000000000\t0.000000000000\t0.000000000000\t0.003700000000\t0.000245194218\t0.000000000000\t0.000600000000\t0.006300000000\t0.000021133639\t0.002016738361\t0.005200000000\t0.000097003127\t0.000000000000\t0.069845218085\t0.000231525618\t0.011202143591\t0.010191520401\t0.000020138391\t0.000000000000\t0.003250427227\t0.002438756000\t0.003639600000\t\t\t\t\t\t\t\r",
"T>A\tTTT\tT[T>A]T\t0.003981022761\t0.000276315891\t0.014712135403\t0.004500000000\t0.010240850413\t0.000700000000\t0.005000000000\t0.027653296254\t0.013900000000\t0.002800000000\t0.000800000000\t0.009400000000\t0.000160345975\t0.001200000000\t0.000300000000\t0.008000000000\t0.003883130589\t0.002818679447\t0.008500000000\t0.002313489060\t0.000800000000\t0.023579061797\t0.000000000000\t0.000270578768\t0.014527539313\t0.001235538470\t0.001929770813\t0.016492812256\t0.000000000000\t0.000000000000\t\t\t\t\t\t\t\r",
"T>C\tTTA\tT[T>C]A\t0.008073616351\t0.000101932352\t0.009115495077\t0.002800000000\t0.017059640976\t0.005900000000\t0.002300000000\t0.004901072379\t0.015600000000\t0.003000000000\t0.000200000000\t0.050000000000\t0.000784056118\t0.001800000000\t0.006800000000\t0.024600000000\t0.001708586711\t0.000302303565\t0.005200000000\t0.005734229562\t0.038900000000\t0.004803669740\t0.000000000000\t0.001682291048\t0.014515050090\t0.028789680542\t0.001721015277\t0.009929552068\t0.011498665964\t0.000000000000\t\t\t\t\t\t\t\r",
"T>C\tTTC\tT[T>C]C\t0.004857381178\t0.000470343659\t0.010953790932\t0.001600000000\t0.014195709147\t0.003500000000\t0.001800000000\t0.006081011299\t0.013700000000\t0.009700000000\t0.000900000000\t0.034200000000\t0.000132974888\t0.005600000000\t0.009900000000\t0.018300000000\t0.011074152072\t0.006260791115\t0.002000000000\t0.011234254668\t0.031300000000\t0.000912157512\t0.000230608944\t0.002344527656\t0.007406083457\t0.036639492442\t0.003489212553\t0.029010563206\t0.003584215038\t0.009148100000\t\t\t\t\t\t\t\r",
"T>C\tTTG\tT[T>C]G\t0.008325454207\t0.000192354203\t0.006113033219\t0.003600000000\t0.012596501722\t0.010600000000\t0.001900000000\t0.001712484699\t0.009800000000\t0.006500000000\t0.000700000000\t0.032000000000\t0.000037402404\t0.003400000000\t0.003900000000\t0.016500000000\t0.005932002763\t0.001541265329\t0.005700000000\t0.023767611345\t0.024900000000\t0.001781735680\t0.000239513775\t0.000481274164\t0.004892531217\t0.019143694280\t0.001152846581\t0.015683218557\t0.001694296670\t0.006000400000\t\t\t\t\t\t\t\r",
"T>C\tTTT\tT[T>C]T\t0.006257105605\t0.000585330275\t0.010774336259\t0.002200000000\t0.017375139145\t0.002900000000\t0.002400000000\t0.010003441830\t0.030900000000\t0.009900000000\t0.000500000000\t0.037100000000\t0.000179702846\t0.002600000000\t0.004700000000\t0.017400000000\t0.007208218396\t0.001357797089\t0.003800000000\t0.005395186301\t0.027400000000\t0.000000000000\t0.000000000000\t0.000358148002\t0.017014780759\t0.023071688831\t0.001527957003\t0.017291563618\t0.003799705772\t0.000000000000\t\t\t\t\t\t\t\r",
"T>G\tTTA\tT[T>G]A\t0.001397553749\t0.000071736950\t0.005427184241\t0.000000000000\t0.005202874420\t0.000900000000\t0.000000000000\t0.001743221362\t0.050200000000\t0.005000000000\t0.000000000000\t0.001900000000\t0.000546581847\t0.000000000000\t0.000000000000\t0.006800000000\t0.000133410642\t0.000036318404\t0.003200000000\t0.000000269147\t0.000000000000\t0.001651698761\t0.000000000000\t0.000000000000\t0.002081479927\t0.000008534935\t0.001189022454\t0.009389658706\t0.000000000000\t0.008656300000\t\t\t\t\t\t\t\r",
"T>G\tTTC\tT[T>G]C\t0.001291736985\t0.000014281456\t0.006160250423\t0.000300000000\t0.005131607914\t0.001900000000\t0.001000000000\t0.002549838329\t0.008100000000\t0.009200000000\t0.000000000000\t0.002700000000\t0.000235347085\t0.001500000000\t0.000400000000\t0.006900000000\t0.009613345245\t0.003233837990\t0.001800000000\t0.000377234437\t0.001800000000\t0.000000000000\t0.000000000000\t0.000000000000\t0.000578978775\t0.002718509845\t0.000280295354\t0.030117076984\t0.000000000000\t0.004328200000\t\t\t\t\t\t\t\r",
"T>G\tTTG\tT[T>G]G\t0.002031076880\t0.000206615168\t0.011076526261\t0.003000000000\t0.006055254100\t0.001100000000\t0.001000000000\t0.006030395185\t0.008800000000\t0.002200000000\t0.000300000000\t0.001100000000\t0.000000479022\t0.000200000000\t0.000900000000\t0.004900000000\t0.004522462318\t0.000754601773\t0.001100000000\t0.000515421599\t0.000300000000\t0.002572751958\t0.000247501933\t0.001360504850\t0.009429195885\t0.001369161173\t0.002353055589\t0.012698750811\t0.000353695951\t0.008262800000\t\t\t\t\t\t\t\r",
"T>G\tTTT\tT[T>G]T\t0.004030128160\t0.000023598204\t0.013000984209\t0.001100000000\t0.013369935834\t0.007200000000\t0.001400000000\t0.007223998893\t0.054500000000\t0.063300000000\t0.000300000000\t0.003200000000\t0.000670588327\t0.002500000000\t0.003300000000\t0.016300000000\t0.058040407760\t0.002126441532\t0.001300000000\t0.000615656691\t0.000300000000\t0.000000000000\t0.000000000000\t0.000069515778\t0.007869671577\t0.002568076726\t0.000139561285\t0.233659783276\t0.006104834136\t0.000000000000\t\t\t\t\t\t\t"
)

# Signature Data :: COSMIC_v3
COSMIC_v3 <- c(
 "Type\tSubType\tSBS1\tSBS2\tSBS3\tSBS4\tSBS5\tSBS6\tSBS7a\tSBS7b\tSBS7c\tSBS7d\tSBS8\tSBS9\tSBS10a\tSBS10b\tSBS11\tSBS12\tSBS13\tSBS14\tSBS15\tSBS16\tSBS17a\tSBS17b\tSBS18\tSBS19\tSBS20\tSBS21\tSBS22\tSBS23\tSBS24\tSBS25\tSBS26\tSBS27\tSBS28\tSBS29\tSBS30\tSBS31\tSBS32\tSBS33\tSBS34\tSBS35\tSBS36\tSBS37\tSBS38\tSBS39\tSBS40\tSBS41\tSBS42\tSBS43\tSBS44\tSBS45\tSBS46\tSBS47\tSBS48\tSBS49\tSBS50\tSBS51\tSBS52\tSBS53\tSBS54\tSBS55\tSBS56\tSBS57\tSBS58\tSBS59\tSBS60\tSBS84\tSBS85\r",
 "C>A\tACA\t8.86E-04\t5.80E-07\t2.08E-02\t4.22E-02\t1.20E-02\t4.25E-04\t6.70E-05\t2.33E-03\t4.83E-03\t4.04E-05\t4.41E-02\t5.58E-04\t2.19E-03\t1.82E-04\t1.46E-04\t4.52E-03\t1.82E-03\t1.12E-03\t9.44E-04\t1.60E-02\t2.07E-03\t6.08E-04\t5.15E-02\t1.27E-03\t6.19E-04\t1.57E-04\t6.01E-03\t8.35E-04\t3.64E-02\t9.90E-03\t8.73E-04\t5.21E-03\t7.84E-04\t6.32E-02\t1.80E-03\t9.54E-03\t2.23E-02\t3.11E-03\t4.87E-03\t8.83E-03\t2.52E-02\t3.95E-03\t1.28E-02\t1.17E-02\t2.82E-02\t2.11E-03\t1.16E-03\t2.95E-02\t7.68E-18\t9.11E-03\t4.40E-03\t6.78E-02\t8.55E-04\t2.51E-02\t1.19E-01\t1.41E-01\t1.52E-02\t5.38E-03\t2.16E-03\t5.88E-03\t1.26E-02\t1.23E-02\t5.89E-02\t3.59E-03\t6.15E-03\t0.003471994832\t0.006080257390\r",
 "C>A\tACC\t2.28E-03\t1.48E-04\t1.65E-02\t3.33E-02\t9.44E-03\t5.24E-04\t1.79E-04\t4.61E-04\t1.15E-03\t7.65E-04\t4.78E-02\t4.09E-03\t1.77E-03\t6.54E-03\t5.52E-04\t1.13E-03\t7.21E-04\t1.31E-02\t4.97E-04\t2.90E-03\t9.18E-04\t1.29E-04\t1.58E-02\t6.41E-04\t1.39E-03\t2.36E-03\t9.54E-05\t3.99E-04\t2.65E-02\t7.00E-03\t5.28E-04\t4.74E-03\t2.53E-03\t5.12E-02\t5.06E-04\t1.85E-02\t1.88E-02\t2.23E-03\t6.96E-03\t4.62E-02\t8.32E-03\t1.45E-03\t1.01E-02\t7.15E-03\t1.34E-02\t1.22E-03\t2.06E-02\t5.90E-03\t1.50E-04\t2.85E-03\t4.70E-03\t2.97E-02\t8.41E-04\t7.03E-03\t1.27E-01\t1.69E-03\t6.54E-03\t1.96E-03\t7.96E-04\t2.05E-03\t1.57E-02\t1.43E-03\t6.75E-03\t2.37E-03\t7.79E-04\t0.005007308997\t0.000879936873\r",
 "C>A\tACG\t1.77E-04\t5.23E-05\t1.75E-03\t1.56E-02\t1.85E-03\t5.20E-05\t7.12E-05\t1.86E-04\t3.77E-04\t2.50E-04\t4.62E-03\t4.26E-04\t1.50E-04\t5.35E-05\t9.41E-05\t5.40E-04\t2.64E-04\t4.13E-04\t4.61E-05\t1.02E-03\t4.76E-05\t5.82E-05\t2.43E-03\t2.46E-04\t2.18E-05\t2.94E-04\t7.65E-04\t9.85E-08\t1.48E-02\t1.45E-03\t1.14E-04\t7.83E-04\t3.53E-04\t1.92E-02\t9.13E-05\t1.66E-03\t4.46E-03\t4.14E-04\t5.20E-05\t1.39E-03\t2.24E-03\t1.06E-03\t1.90E-03\t2.67E-03\t2.92E-03\t6.14E-05\t3.34E-05\t6.26E-04\t9.16E-07\t1.66E-03\t3.00E-04\t1.38E-03\t2.58E-04\t7.11E-03\t5.77E-03\t4.80E-03\t4.14E-03\t3.81E-02\t1.64E-03\t4.49E-05\t2.06E-04\t3.33E-03\t8.25E-04\t1.42E-04\t2.23E-16\t0.000452641972\t0.000306011349\r",
 "C>A\tACT\t1.28E-03\t9.78E-05\t1.22E-02\t2.95E-02\t6.61E-03\t1.80E-04\t2.48E-04\t7.10E-04\t1.96E-03\t4.05E-03\t4.70E-02\t3.05E-03\t1.70E-02\t1.63E-05\t2.66E-04\t1.22E-03\t3.48E-04\t8.26E-02\t1.11E-03\t1.06E-02\t6.18E-05\t4.56E-04\t2.14E-02\t5.71E-04\t1.24E-03\t6.21E-04\t1.72E-03\t5.61E-18\t2.20E-02\t4.97E-03\t6.19E-04\t2.72E-03\t3.97E-03\t3.56E-02\t5.56E-04\t6.28E-03\t1.46E-02\t1.90E-03\t1.23E-03\t2.16E-02\t1.79E-02\t1.85E-03\t8.85E-03\t7.40E-03\t1.48E-02\t1.33E-03\t7.97E-03\t2.07E-03\t5.78E-03\t9.64E-03\t3.40E-03\t1.19E-02\t3.89E-05\t7.97E-03\t8.32E-02\t2.12E-02\t9.24E-03\t1.71E-05\t4.14E-04\t1.17E-03\t2.30E-02\t8.64E-03\t5.58E-03\t1.45E-02\t4.44E-04\t0.009295948365\t0.002717279233\r",
 "C>A\tCCA\t3.12E-04\t2.08E-04\t2.25E-02\t8.07E-02\t7.43E-03\t1.82E-03\t4.55E-04\t1.14E-03\t1.09E-04\t1.45E-02\t4.01E-02\t4.80E-03\t3.18E-03\t5.21E-04\t6.18E-04\t5.42E-03\t1.40E-03\t1.63E-02\t4.18E-02\t1.04E-02\t2.95E-04\t2.71E-04\t7.40E-02\t1.32E-03\t3.76E-02\t2.22E-16\t1.83E-03\t1.51E-04\t2.54E-02\t1.48E-02\t1.41E-03\t5.07E-03\t2.03E-03\t5.42E-02\t7.68E-04\t1.07E-02\t1.24E-03\t2.70E-03\t1.95E-03\t2.99E-02\t5.58E-02\t3.22E-03\t3.81E-01\t1.01E-02\t2.08E-02\t2.38E-03\t1.16E-02\t6.63E-03\t2.84E-05\t2.38E-01\t1.56E-02\t2.46E-02\t1.18E-03\t2.29E-01\t2.25E-02\t6.98E-03\t5.74E-02\t3.34E-01\t6.32E-05\t6.25E-03\t2.48E-03\t8.70E-03\t6.18E-03\t2.23E-03\t1.18E-03\t0.000000000000\t0.004264235123\r",
 "C>A\tCCC\t1.79E-03\t9.53E-05\t2.53E-02\t7.97E-02\t6.14E-03\t4.09E-03\t2.40E-04\t1.55E-03\t3.09E-03\t3.22E-03\t3.88E-02\t1.92E-03\t2.02E-04\t3.18E-03\t9.86E-04\t4.28E-03\t9.68E-04\t1.97E-02\t5.05E-03\t9.03E-03\t1.16E-03\t4.11E-05\t1.96E-02\t7.67E-04\t9.89E-02\t2.60E-03\t4.43E-03\t6.30E-04\t3.49E-02\t7.82E-03\t1.51E-03\t2.23E-03\t1.09E-04\t3.69E-02\t3.88E-04\t1.16E-02\t4.92E-03\t1.91E-03\t3.67E-05\t3.30E-02\t2.02E-02\t8.14E-04\t6.95E-02\t9.06E-03\t1.47E-02\t9.19E-04\t3.28E-02\t8.20E-05\t3.03E-02\t1.51E-01\t7.20E-03\t5.27E-03\t3.57E-05\t2.86E-02\t2.89E-03\t3.17E-03\t1.90E-02\t3.55E-03\t2.57E-03\t3.75E-03\t6.68E-03\t3.94E-03\t1.83E-03\t5.93E-03\t2.53E-04\t0.005218571801\t0.000306234865\r",
 "C>A\tCCG\t9.32E-05\t2.23E-16\t2.51E-03\t2.45E-02\t3.46E-03\t1.65E-03\t2.23E-16\t4.08E-04\t1.91E-03\t8.51E-05\t3.41E-03\t1.37E-03\t4.64E-05\t8.38E-05\t9.20E-06\t1.17E-03\t4.33E-06\t1.38E-03\t6.62E-05\t1.74E-03\t2.22E-16\t2.23E-16\t1.20E-02\t2.36E-04\t7.43E-03\t1.98E-03\t1.36E-04\t1.30E-04\t1.60E-02\t1.28E-03\t3.62E-04\t2.66E-04\t2.22E-16\t1.48E-02\t4.07E-04\t3.29E-04\t1.66E-05\t3.71E-04\t2.98E-04\t2.20E-03\t7.10E-03\t7.47E-04\t5.10E-02\t5.74E-03\t2.54E-03\t2.61E-04\t1.88E-03\t2.56E-05\t3.74E-03\t6.79E-02\t4.50E-03\t2.71E-04\t1.37E-01\t5.15E-01\t5.18E-18\t4.30E-03\t1.36E-02\t6.77E-02\t5.93E-03\t6.41E-04\t6.72E-18\t1.07E-03\t1.46E-03\t8.15E-04\t3.79E-04\t0.000069231298\t0.000157256840\r",
 "C>A\tCCT\t2.23E-16\t4.21E-04\t1.54E-02\t6.90E-02\t6.49E-03\t9.60E-03\t3.72E-04\t1.49E-03\t1.35E-02\t1.17E-04\t3.30E-02\t3.06E-03\t1.73E-02\t8.33E-03\t1.50E-04\t1.10E-03\t2.22E-16\t2.60E-01\t2.87E-02\t2.49E-03\t1.37E-04\t1.82E-04\t3.63E-02\t3.44E-05\t2.52E-01\t1.91E-02\t1.10E-04\t5.80E-04\t1.71E-02\t1.26E-02\t7.14E-03\t3.10E-03\t2.03E-04\t5.19E-02\t2.97E-04\t2.79E-02\t2.67E-03\t2.53E-03\t6.20E-04\t8.25E-02\t6.05E-02\t1.77E-03\t1.84E-01\t6.03E-03\t1.47E-02\t2.65E-03\t2.11E-02\t5.22E-18\t8.67E-02\t1.94E-01\t8.30E-03\t6.82E-04\t8.45E-04\t4.90E-02\t1.15E-03\t5.07E-03\t1.88E-02\t4.92E-03\t5.53E-03\t2.84E-05\t7.13E-02\t8.21E-03\t1.96E-03\t2.02E-03\t5.41E-05\t0.005640381180\t0.000860850681\r",
 "C>A\tGCA\t1.58E-03\t8.62E-05\t7.13E-03\t3.16E-02\t1.02E-02\t4.27E-04\t4.06E-05\t7.46E-04\t1.05E-03\t1.68E-04\t2.42E-02\t7.72E-03\t1.21E-03\t7.31E-04\t1.65E-04\t1.98E-03\t9.49E-04\t2.08E-03\t2.29E-03\t8.61E-03\t1.78E-03\t9.98E-04\t1.09E-01\t8.14E-04\t1.08E-02\t2.14E-03\t2.23E-03\t3.37E-04\t6.22E-02\t1.35E-02\t4.96E-04\t1.08E-02\t1.56E-04\t9.71E-02\t9.43E-04\t1.31E-02\t4.51E-03\t4.41E-03\t2.81E-03\t9.87E-03\t6.77E-02\t2.04E-03\t1.20E-02\t6.88E-03\t1.42E-02\t3.09E-03\t7.71E-03\t9.30E-03\t7.68E-18\t7.55E-03\t7.20E-03\t1.11E-02\t7.52E-04\t7.88E-03\t1.42E-02\t3.94E-03\t8.24E-03\t4.32E-04\t5.19E-03\t8.55E-04\t2.17E-04\t1.89E-03\t2.74E-03\t1.03E-01\t5.17E-04\t0.002703472389\t0.001439703098\r",
 "C>A\tGCC\t3.39E-04\t2.23E-16\t1.09E-02\t3.48E-02\t7.65E-03\t9.13E-04\t1.38E-04\t2.15E-04\t7.42E-04\t2.22E-16\t2.62E-02\t6.68E-03\t4.41E-04\t2.55E-03\t4.05E-04\t4.17E-04\t4.51E-04\t2.24E-02\t2.49E-03\t4.25E-03\t1.28E-03\t2.78E-06\t1.73E-02\t2.72E-04\t5.65E-03\t7.07E-03\t1.73E-03\t5.41E-04\t1.28E-01\t6.47E-03\t6.66E-04\t1.24E-02\t1.02E-03\t9.03E-02\t4.94E-04\t2.15E-02\t3.62E-03\t9.68E-04\t4.68E-04\t5.02E-02\t1.92E-02\t6.29E-04\t6.28E-03\t6.95E-03\t1.02E-02\t7.23E-04\t2.68E-02\t1.69E-03\t4.67E-03\t5.73E-03\t3.00E-03\t3.48E-03\t3.46E-04\t2.48E-03\t7.58E-03\t3.28E-03\t7.32E-03\t3.44E-03\t9.13E-03\t1.92E-04\t7.34E-03\t7.52E-04\t1.08E-03\t1.29E-01\t8.66E-04\t0.010228858800\t0.000504268015\r",
 "C>A\tGCG\t5.87E-04\t1.39E-05\t1.43E-03\t1.53E-02\t2.34E-03\t3.57E-04\t4.25E-05\t7.29E-05\t2.22E-16\t1.82E-05\t2.74E-03\t1.71E-04\t1.05E-05\t8.52E-05\t1.39E-05\t4.47E-04\t1.64E-04\t1.11E-03\t4.04E-04\t1.36E-03\t1.06E-04\t9.24E-05\t7.56E-03\t2.05E-04\t1.57E-03\t1.69E-03\t5.15E-04\t1.19E-04\t6.59E-02\t6.58E-04\t1.29E-04\t1.63E-03\t3.95E-05\t2.30E-02\t2.22E-16\t9.72E-04\t8.77E-04\t3.42E-04\t2.27E-04\t5.07E-04\t4.83E-03\t5.06E-04\t1.37E-03\t4.53E-03\t2.55E-03\t2.20E-04\t1.06E-03\t1.77E-04\t6.05E-04\t2.09E-03\t4.00E-04\t7.83E-04\t7.39E-04\t6.83E-03\t2.44E-04\t4.28E-03\t4.04E-03\t5.08E-02\t5.64E-03\t7.72E-05\t9.79E-06\t6.57E-18\t2.73E-03\t1.21E-02\t3.90E-04\t0.000763601342\t0.000362282705\r",
 "C>A\tGCT\t2.23E-16\t5.12E-05\t1.00E-02\t2.12E-02\t6.82E-03\t1.23E-03\t2.23E-16\t1.47E-04\t3.63E-04\t7.40E-05\t2.69E-02\t8.56E-03\t1.11E-02\t2.80E-03\t2.22E-16\t2.75E-05\t6.15E-05\t1.26E-01\t1.82E-02\t5.80E-03\t4.70E-05\t5.34E-04\t6.23E-02\t2.73E-04\t2.49E-03\t8.90E-03\t5.51E-04\t3.09E-05\t5.44E-02\t1.05E-02\t8.17E-04\t7.41E-03\t1.58E-04\t8.46E-02\t8.62E-04\t3.51E-03\t3.27E-03\t3.19E-03\t1.35E-05\t2.29E-02\t6.90E-02\t1.29E-03\t9.09E-03\t5.89E-03\t1.08E-02\t1.95E-03\t1.81E-02\t1.39E-03\t1.13E-02\t6.99E-03\t2.40E-03\t1.41E-03\t1.99E-05\t5.84E-03\t2.04E-03\t3.45E-03\t4.21E-03\t1.14E-03\t6.73E-03\t7.38E-04\t1.57E-02\t1.95E-03\t3.50E-04\t4.99E-01\t1.33E-04\t0.010179165989\t0.000762947056\r",
 "C>A\tTCA\t6.58E-05\t6.39E-04\t8.41E-03\t3.98E-02\t7.86E-03\t3.56E-04\t3.90E-04\t2.90E-03\t2.52E-03\t1.14E-02\t2.98E-02\t1.70E-02\t9.40E-02\t4.46E-04\t1.80E-04\t3.48E-03\t7.90E-02\t3.47E-04\t1.72E-02\t1.08E-02\t1.90E-04\t1.24E-05\t7.38E-02\t1.38E-03\t1.71E-03\t1.46E-03\t7.84E-04\t2.13E-03\t2.98E-02\t6.94E-03\t9.38E-04\t7.08E-03\t1.43E-02\t4.52E-02\t6.96E-04\t1.86E-03\t9.73E-03\t1.19E-02\t4.27E-03\t5.27E-03\t1.41E-01\t2.00E-03\t1.77E-02\t5.98E-03\t2.60E-02\t3.31E-03\t5.10E-03\t6.30E-03\t8.53E-04\t6.40E-02\t1.23E-02\t1.69E-02\t1.24E-01\t3.29E-02\t1.90E-02\t3.57E-02\t4.64E-01\t2.67E-03\t6.69E-04\t1.74E-03\t2.19E-01\t1.43E-02\t7.37E-03\t2.37E-05\t2.02E-03\t0.002083564327\t0.001944387330\r",
 "C>A\tTCC\t2.53E-03\t1.70E-04\t1.68E-02\t4.69E-02\t9.17E-03\t4.98E-04\t6.36E-04\t3.52E-03\t1.42E-02\t7.76E-03\t2.51E-02\t2.46E-03\t9.97E-03\t3.80E-03\t9.38E-04\t2.88E-03\t2.94E-02\t2.43E-02\t5.81E-03\t5.87E-03\t2.59E-05\t2.23E-16\t4.36E-02\t1.13E-03\t3.37E-03\t2.36E-04\t1.74E-03\t6.64E-04\t5.19E-02\t1.02E-02\t6.30E-04\t6.27E-03\t4.52E-04\t4.11E-02\t2.43E-04\t9.38E-03\t9.29E-03\t5.11E-03\t1.22E-03\t2.31E-02\t7.72E-02\t1.59E-03\t1.15E-02\t7.86E-03\t2.10E-02\t1.72E-03\t5.40E-02\t5.37E-06\t2.22E-03\t6.42E-02\t7.10E-03\t8.11E-05\t3.79E-03\t5.11E-03\t8.86E-04\t2.86E-03\t4.68E-02\t3.64E-03\t4.05E-05\t8.28E-18\t1.32E-02\t7.09E-04\t4.29E-04\t2.02E-02\t1.63E-04\t0.000009064386\t0.000169427568\r",
 "C>A\tTCG\t2.23E-16\t8.53E-05\t1.41E-03\t1.09E-02\t2.30E-03\t2.00E-04\t3.44E-05\t3.00E-04\t2.22E-16\t2.22E-16\t1.90E-03\t4.45E-04\t1.16E-03\t5.73E-04\t2.22E-16\t5.05E-04\t1.21E-02\t7.08E-06\t7.69E-05\t1.49E-03\t2.22E-16\t9.51E-05\t1.27E-02\t2.85E-04\t2.22E-16\t6.83E-05\t2.22E-16\t9.12E-05\t9.79E-03\t0.00E+00\t1.23E-04\t1.45E-03\t1.29E-04\t1.74E-02\t8.26E-05\t9.55E-04\t2.38E-03\t1.37E-03\t2.22E-16\t9.47E-05\t1.47E-02\t1.12E-03\t2.38E-03\t3.01E-03\t3.17E-03\t2.07E-04\t9.04E-04\t3.07E-18\t3.15E-07\t1.83E-02\t2.50E-03\t2.39E-05\t7.09E-01\t4.05E-02\t4.19E-05\t4.54E-03\t1.34E-01\t4.47E-04\t8.28E-04\t9.91E-05\t5.46E-03\t2.23E-04\t2.73E-03\t6.34E-04\t2.05E-04\t0.000000000000\t0.000000003989\r",
 "C>A\tTCT\t5.88E-06\t4.40E-04\t9.97E-03\t4.13E-02\t1.27E-02\t7.47E-04\t6.13E-04\t1.41E-03\t3.73E-03\t1.91E-04\t3.19E-02\t2.73E-02\t6.70E-01\t8.16E-02\t4.73E-05\t2.35E-03\t3.55E-02\t1.51E-01\t8.05E-03\t9.25E-03\t2.88E-04\t2.57E-03\t1.22E-01\t1.38E-03\t9.04E-04\t2.22E-16\t3.96E-04\t2.40E-04\t2.28E-02\t2.46E-02\t1.40E-03\t9.81E-03\t5.57E-04\t4.11E-02\t6.65E-04\t4.49E-03\t9.64E-03\t1.28E-02\t3.20E-05\t1.53E-02\t2.08E-01\t2.44E-03\t2.78E-02\t6.55E-03\t2.41E-02\t5.90E-03\t2.73E-02\t6.52E-18\t1.66E-02\t6.20E-02\t3.00E-03\t1.92E-04\t1.87E-04\t8.75E-03\t1.85E-03\t3.70E-04\t6.10E-02\t8.71E-04\t7.20E-03\t4.06E-03\t4.02E-01\t6.19E-03\t1.47E-02\t3.83E-03\t4.06E-04\t0.001203779266\t0.002057492391\r",
 "C>G\tACA\t1.86E-03\t2.23E-16\t1.97E-02\t6.89E-03\t1.01E-02\t4.71E-04\t6.49E-05\t8.56E-06\t1.12E-03\t1.18E-03\t4.33E-03\t4.84E-03\t2.23E-16\t3.05E-04\t1.70E-04\t1.22E-03\t3.87E-03\t2.22E-16\t1.13E-04\t1.74E-03\t1.01E-03\t1.46E-04\t1.73E-03\t3.19E-03\t8.77E-03\t1.69E-05\t1.66E-03\t5.74E-18\t3.34E-03\t8.03E-03\t4.27E-04\t1.32E-03\t1.54E-04\t3.14E-03\t1.07E-03\t8.32E-03\t1.99E-03\t1.41E-03\t7.62E-04\t3.58E-03\t1.84E-03\t3.43E-02\t2.69E-03\t4.73E-02\t1.22E-02\t5.33E-03\t4.06E-04\t3.61E-04\t3.18E-03\t3.11E-03\t1.20E-03\t1.03E-03\t1.98E-04\t2.88E-04\t3.85E-03\t3.96E-03\t1.73E-03\t7.03E-03\t1.55E-03\t2.85E-03\t4.18E-04\t1.54E-02\t2.06E-02\t2.94E-03\t3.77E-04\t0.006758077485\t0.007234927130\r",
 "C>G\tACC\t1.22E-03\t1.33E-04\t1.17E-02\t2.84E-03\t5.70E-03\t2.89E-04\t2.16E-05\t1.89E-04\t3.19E-04\t2.55E-03\t2.95E-03\t1.96E-03\t2.43E-05\t1.34E-04\t2.90E-04\t6.42E-04\t9.51E-04\t5.78E-06\t2.69E-04\t2.01E-03\t5.69E-04\t4.34E-05\t2.59E-03\t1.79E-03\t1.51E-03\t7.19E-04\t1.95E-04\t4.08E-05\t3.67E-03\t1.64E-03\t3.78E-04\t1.77E-03\t2.24E-04\t2.63E-03\t4.67E-04\t3.16E-03\t2.03E-03\t1.30E-03\t1.13E-03\t4.55E-03\t1.67E-03\t1.12E-02\t1.58E-03\t2.25E-02\t8.68E-03\t3.16E-03\t5.06E-04\t3.13E-05\t2.33E-04\t1.93E-03\t1.50E-03\t6.17E-04\t1.14E-04\t1.35E-04\t1.10E-03\t4.54E-03\t1.52E-03\t5.46E-03\t4.96E-04\t1.68E-04\t1.66E-04\t3.13E-03\t1.74E-03\t2.44E-03\t5.88E-04\t0.010811956859\t0.002551355408\r",
 "C>G\tACG\t1.15E-04\t1.52E-05\t2.53E-04\t1.28E-03\t1.72E-03\t2.22E-16\t5.08E-05\t5.05E-04\t2.22E-16\t9.06E-04\t2.86E-04\t1.36E-03\t2.23E-16\t2.22E-16\t5.96E-05\t7.51E-05\t2.64E-04\t3.91E-05\t3.07E-05\t8.47E-04\t1.52E-04\t2.23E-16\t1.92E-03\t6.16E-04\t4.69E-04\t2.22E-16\t2.31E-04\t2.86E-18\t9.20E-04\t0.00E+00\t6.83E-05\t0.00E+00\t2.22E-16\t8.96E-04\t1.31E-04\t3.00E-03\t1.13E-03\t2.53E-04\t2.22E-16\t3.46E-04\t5.00E-04\t5.35E-03\t6.36E-04\t1.36E-02\t2.25E-03\t1.21E-03\t1.77E-05\t4.43E-03\t3.14E-05\t2.70E-04\t0.00E+00\t3.28E-04\t8.03E-05\t2.04E-04\t8.27E-04\t1.97E-03\t8.18E-04\t2.51E-03\t1.70E-03\t5.10E-04\t4.28E-04\t2.77E-03\t1.69E-03\t8.30E-04\t2.76E-04\t0.000045986523\t0.000422196736\r",
 "C>G\tACT\t1.14E-03\t9.12E-05\t1.74E-02\t3.55E-03\t1.01E-02\t2.54E-04\t9.56E-05\t2.85E-04\t1.93E-03\t2.91E-03\t5.71E-03\t6.41E-03\t2.23E-16\t2.22E-16\t5.35E-04\t1.44E-03\t4.00E-03\t2.04E-04\t3.25E-04\t1.15E-02\t1.96E-05\t3.19E-04\t4.08E-03\t2.14E-03\t4.42E-03\t6.19E-05\t6.20E-05\t5.61E-18\t2.85E-03\t3.43E-03\t7.28E-04\t1.32E-03\t8.33E-05\t2.49E-03\t7.31E-04\t3.78E-03\t1.75E-03\t1.30E-03\t2.03E-03\t9.48E-03\t1.91E-03\t1.88E-02\t2.44E-03\t4.34E-02\t1.46E-02\t8.27E-03\t4.31E-04\t3.54E-04\t3.29E-03\t1.43E-03\t2.20E-03\t1.72E-03\t6.86E-04\t1.76E-04\t1.97E-03\t9.00E-04\t1.71E-03\t3.04E-03\t5.18E-03\t1.27E-03\t1.23E-04\t5.07E-03\t2.03E-03\t3.83E-03\t2.50E-04\t0.029594216358\t0.013019760586\r",
 "C>G\tCCA\t2.41E-05\t2.23E-16\t1.93E-02\t7.83E-03\t6.96E-03\t4.11E-05\t2.23E-16\t4.58E-04\t2.06E-04\t1.32E-03\t5.00E-03\t8.93E-04\t4.57E-05\t2.99E-05\t2.61E-04\t1.34E-03\t5.61E-03\t8.55E-05\t8.36E-05\t3.55E-03\t5.53E-04\t1.94E-04\t6.13E-04\t1.77E-03\t7.69E-05\t1.97E-05\t1.05E-03\t1.74E-04\t2.61E-03\t7.20E-03\t3.64E-04\t4.73E-03\t1.44E-05\t3.29E-03\t1.27E-03\t2.05E-03\t1.15E-03\t1.50E-03\t1.19E-04\t3.86E-03\t1.24E-03\t4.42E-04\t1.69E-03\t2.88E-02\t7.56E-03\t2.08E-03\t4.26E-03\t2.09E-05\t5.36E-05\t1.43E-03\t5.00E-04\t9.56E-04\t9.93E-05\t2.16E-04\t6.56E-05\t2.85E-03\t1.18E-03\t2.72E-03\t1.13E-03\t2.94E-03\t2.04E-04\t8.88E-04\t7.94E-04\t2.08E-03\t4.71E-05\t0.000007117033\t0.001223945047\r",
 "C>G\tCCC\t7.68E-05\t3.51E-04\t1.41E-02\t6.87E-03\t9.07E-03\t3.72E-05\t1.31E-04\t9.60E-04\t2.26E-03\t2.86E-04\t5.46E-03\t1.54E-03\t2.23E-16\t2.22E-16\t1.77E-04\t8.91E-04\t6.48E-04\t7.86E-05\t1.70E-05\t4.54E-03\t2.35E-04\t4.66E-05\t2.26E-03\t1.81E-03\t4.38E-06\t2.73E-04\t7.85E-04\t2.55E-04\t5.22E-03\t4.60E-03\t3.22E-04\t1.26E-03\t7.55E-05\t6.34E-03\t7.93E-04\t1.86E-03\t8.24E-04\t4.10E-04\t1.18E-04\t6.30E-03\t1.37E-03\t3.13E-04\t1.26E-03\t3.77E-02\t9.35E-03\t1.20E-03\t6.62E-03\t2.12E-04\t2.49E-05\t8.00E-04\t3.00E-04\t2.37E-04\t1.22E-04\t2.56E-04\t2.89E-04\t1.20E-04\t8.35E-04\t6.03E-03\t8.83E-04\t2.86E-03\t2.25E-04\t2.87E-03\t1.28E-03\t1.73E-03\t1.42E-04\t0.003440786012\t0.000210178860\r",
 "C>G\tCCG\t3.52E-04\t4.87E-05\t1.50E-03\t3.40E-03\t2.49E-03\t8.43E-05\t7.75E-05\t7.91E-04\t2.62E-04\t1.32E-03\t1.05E-05\t4.37E-04\t2.01E-05\t2.22E-16\t1.12E-05\t1.73E-04\t4.88E-04\t1.91E-05\t5.02E-05\t1.46E-03\t2.32E-04\t2.23E-16\t3.58E-03\t4.81E-05\t1.88E-04\t1.01E-05\t2.39E-04\t1.86E-18\t9.45E-04\t8.24E-04\t4.17E-05\t3.61E-04\t2.22E-16\t1.30E-03\t1.20E-04\t8.54E-04\t5.37E-04\t3.18E-04\t1.90E-04\t8.70E-04\t4.19E-04\t3.21E-04\t9.57E-04\t1.23E-02\t1.68E-03\t4.11E-04\t4.70E-04\t2.54E-03\t7.68E-18\t6.43E-04\t1.00E-04\t6.70E-04\t3.63E-05\t2.38E-04\t1.42E-03\t2.03E-03\t5.68E-04\t1.62E-03\t4.34E-04\t2.64E-04\t6.66E-05\t7.63E-04\t5.11E-05\t6.90E-04\t3.37E-04\t0.000002195949\t0.000040247509\r",
 "C>G\tCCT\t2.23E-16\t6.36E-05\t2.23E-02\t6.82E-03\t9.53E-03\t1.44E-05\t9.56E-05\t8.90E-04\t4.24E-03\t2.22E-16\t5.49E-03\t3.03E-03\t3.01E-05\t2.22E-16\t9.64E-06\t1.56E-03\t6.40E-03\t1.78E-04\t8.45E-05\t9.62E-03\t3.90E-05\t1.93E-04\t3.79E-03\t3.22E-03\t1.01E-05\t1.73E-04\t3.20E-04\t2.27E-04\t4.71E-03\t1.96E-04\t6.25E-04\t2.68E-04\t2.22E-16\t6.69E-03\t1.02E-03\t3.42E-03\t4.95E-04\t6.62E-04\t7.36E-05\t1.62E-02\t2.09E-03\t5.81E-04\t1.27E-03\t4.02E-02\t1.39E-02\t5.03E-03\t6.84E-03\t8.53E-04\t8.64E-05\t1.91E-03\t1.50E-03\t1.37E-04\t2.45E-04\t1.91E-04\t5.12E-04\t2.61E-03\t1.20E-03\t5.08E-04\t1.33E-03\t3.31E-03\t7.69E-04\t2.87E-03\t1.25E-03\t3.46E-03\t1.46E-04\t0.013308248775\t0.000484083725\r",
 "C>G\tGCA\t9.59E-06\t6.54E-05\t1.27E-02\t5.29E-03\t4.68E-03\t1.94E-05\t1.96E-05\t5.11E-05\t1.82E-04\t1.96E-04\t2.58E-03\t1.36E-03\t1.50E-05\t1.72E-05\t4.67E-05\t5.64E-04\t1.50E-03\t2.22E-16\t9.94E-04\t1.05E-03\t8.78E-04\t5.99E-05\t6.35E-04\t1.31E-03\t6.69E-04\t3.88E-04\t1.40E-03\t4.32E-18\t2.72E-03\t4.17E-03\t3.95E-04\t1.42E-03\t1.83E-04\t2.27E-03\t5.15E-04\t1.83E-03\t8.76E-04\t7.71E-04\t3.34E-05\t3.38E-03\t8.86E-04\t1.32E-03\t1.20E-03\t2.40E-02\t5.06E-03\t2.04E-03\t2.40E-03\t1.36E-02\t2.19E-04\t3.33E-04\t1.00E-04\t1.68E-03\t7.84E-05\t1.52E-04\t2.67E-02\t3.41E-03\t1.38E-03\t1.38E-03\t2.19E-03\t8.51E-04\t1.87E-04\t3.79E-03\t2.12E-03\t1.44E-03\t1.63E-04\t0.011151059622\t0.001595924101\r",
 "C>G\tGCC\t1.64E-04\t2.61E-04\t8.63E-03\t5.75E-03\t5.08E-03\t6.85E-04\t3.32E-05\t4.52E-04\t4.73E-04\t2.22E-16\t2.39E-03\t2.83E-03\t8.31E-05\t1.36E-04\t1.14E-04\t2.89E-04\t4.79E-04\t4.33E-05\t6.08E-04\t3.09E-03\t8.35E-04\t5.55E-05\t1.34E-03\t1.29E-03\t2.86E-03\t1.42E-03\t9.30E-04\t1.27E-04\t5.79E-03\t2.39E-03\t5.79E-04\t1.35E-03\t2.22E-16\t1.64E-03\t2.56E-04\t7.51E-03\t9.00E-04\t7.76E-04\t1.93E-04\t2.14E-02\t1.08E-03\t9.38E-04\t1.08E-03\t2.13E-02\t5.37E-03\t2.69E-03\t1.42E-03\t1.52E-02\t1.19E-03\t1.01E-03\t8.00E-04\t2.28E-03\t2.17E-04\t1.10E-04\t3.30E-02\t1.14E-03\t1.45E-03\t8.47E-02\t2.20E-04\t1.41E-03\t9.83E-05\t2.57E-03\t1.49E-03\t1.72E-03\t1.35E-04\t0.015758341831\t0.005060362002\r",
 "C>G\tGCG\t1.66E-04\t2.23E-16\t2.42E-03\t3.30E-03\t1.51E-03\t2.22E-16\t4.51E-05\t3.37E-04\t1.80E-04\t5.45E-04\t4.88E-04\t5.21E-04\t2.23E-16\t2.22E-16\t3.24E-06\t9.54E-05\t2.41E-04\t2.22E-16\t1.01E-05\t8.13E-05\t2.37E-05\t2.23E-16\t6.39E-04\t5.97E-04\t6.92E-04\t1.17E-04\t2.11E-04\t1.83E-18\t3.14E-03\t4.18E-04\t7.38E-05\t0.00E+00\t2.22E-16\t6.88E-06\t1.17E-04\t4.57E-04\t6.86E-05\t4.66E-05\t2.22E-16\t5.91E-04\t1.72E-04\t5.77E-04\t1.02E-03\t9.86E-03\t1.19E-03\t1.56E-04\t2.23E-04\t1.28E-02\t7.68E-18\t6.95E-04\t0.00E+00\t4.65E-04\t4.31E-05\t3.67E-05\t2.80E-03\t2.67E-04\t1.01E-03\t3.00E-03\t1.59E-03\t5.72E-03\t5.62E-06\t1.86E-03\t1.71E-03\t3.59E-04\t2.02E-04\t0.001697967239\t0.000026827588\r",
 "C>G\tGCT\t2.23E-16\t1.34E-04\t1.47E-02\t3.04E-03\t6.72E-03\t3.21E-04\t4.20E-05\t2.77E-04\t1.28E-03\t8.82E-04\t4.86E-03\t8.57E-03\t3.54E-05\t1.13E-04\t1.86E-04\t5.04E-04\t2.70E-03\t8.03E-05\t4.51E-04\t3.23E-03\t7.87E-04\t2.23E-16\t1.23E-03\t3.20E-04\t3.02E-03\t2.07E-04\t4.55E-04\t4.04E-18\t3.38E-03\t1.99E-03\t3.05E-04\t4.56E-03\t2.22E-16\t3.21E-03\t3.99E-04\t3.36E-03\t2.28E-04\t9.31E-04\t2.69E-05\t8.74E-03\t1.05E-03\t1.27E-03\t9.17E-04\t3.24E-02\t6.98E-03\t3.56E-03\t2.77E-03\t1.80E-02\t1.85E-03\t9.32E-04\t6.00E-04\t2.22E-03\t6.53E-05\t2.08E-04\t1.43E-02\t2.85E-03\t1.02E-03\t5.73E-04\t4.94E-03\t6.82E-05\t1.75E-06\t4.36E-04\t1.08E-03\t2.26E-03\t6.43E-05\t0.076057560709\t0.003186103940\r",
 "C>G\tTCA\t2.23E-16\t2.23E-16\t1.33E-02\t3.58E-03\t7.81E-03\t4.97E-05\t2.23E-16\t2.54E-04\t1.94E-03\t3.74E-03\t1.47E-03\t2.69E-03\t5.55E-05\t2.22E-16\t1.56E-05\t9.69E-04\t3.17E-01\t1.97E-05\t2.14E-04\t2.51E-03\t3.75E-03\t1.60E-04\t2.14E-03\t1.75E-03\t1.29E-04\t6.54E-04\t3.56E-03\t2.11E-04\t3.75E-03\t9.80E-03\t1.31E-03\t2.87E-03\t1.15E-04\t1.07E-04\t1.27E-04\t3.96E-04\t4.71E-18\t1.30E-02\t6.18E-07\t6.63E-03\t3.98E-03\t1.79E-03\t2.66E-03\t3.47E-02\t1.25E-02\t1.99E-02\t5.59E-04\t4.42E-06\t6.73E-04\t3.30E-03\t2.80E-03\t1.38E-04\t1.37E-04\t4.80E-05\t9.32E-04\t2.05E-03\t1.98E-03\t3.74E-03\t4.69E-03\t7.06E-03\t1.13E-03\t2.28E-03\t2.13E-03\t2.35E-06\t5.29E-04\t0.000000000000\t0.000947995470\r",
 "C>G\tTCC\t2.03E-04\t9.87E-05\t2.03E-02\t5.95E-03\t7.60E-03\t1.13E-04\t5.09E-05\t1.62E-03\t4.61E-03\t5.88E-04\t1.34E-03\t2.00E-03\t2.23E-16\t2.22E-16\t2.24E-04\t1.55E-03\t6.35E-02\t1.50E-04\t6.12E-05\t5.23E-03\t7.23E-04\t6.18E-05\t4.64E-03\t3.12E-03\t2.82E-05\t6.31E-04\t1.52E-03\t2.66E-04\t5.58E-03\t3.27E-03\t7.45E-04\t3.11E-03\t1.36E-04\t6.83E-03\t9.86E-04\t2.92E-03\t1.10E-03\t2.51E-03\t1.40E-04\t5.63E-03\t2.41E-03\t1.11E-03\t2.18E-03\t3.67E-02\t1.41E-02\t1.42E-02\t1.87E-03\t1.03E-03\t9.40E-06\t2.01E-03\t2.40E-03\t1.19E-03\t3.50E-04\t1.53E-04\t4.93E-04\t3.61E-03\t1.44E-03\t4.99E-03\t5.17E-04\t3.28E-03\t2.94E-04\t2.80E-03\t1.80E-03\t2.09E-03\t3.48E-04\t0.000113943392\t0.001869404519\r",
 "C>G\tTCG\t6.94E-05\t4.12E-05\t5.20E-04\t1.83E-03\t1.71E-03\t1.95E-06\t6.35E-06\t3.64E-04\t8.75E-06\t8.73E-04\t7.82E-05\t4.04E-04\t2.23E-16\t2.22E-16\t2.22E-16\t2.14E-04\t1.38E-02\t1.83E-05\t5.84E-05\t4.39E-04\t3.80E-04\t1.55E-05\t8.99E-04\t4.41E-04\t1.24E-04\t6.27E-05\t3.81E-04\t2.56E-18\t1.35E-03\t2.06E-03\t9.41E-05\t1.28E-05\t2.46E-05\t8.97E-04\t1.06E-04\t4.16E-04\t9.97E-04\t8.04E-05\t2.22E-16\t1.55E-04\t3.02E-04\t8.60E-04\t1.10E-03\t9.62E-03\t1.74E-03\t1.93E-03\t2.07E-04\t1.47E-03\t7.68E-18\t1.78E-04\t1.00E-04\t0.00E+00\t2.82E-05\t6.33E-05\t1.03E-04\t2.16E-03\t5.81E-04\t1.87E-03\t1.23E-03\t2.29E-04\t8.39E-05\t1.33E-03\t6.30E-04\t4.05E-04\t6.59E-05\t0.000000720571\t0.000169996975\r",
 "C>G\tTCT\t1.56E-03\t1.18E-03\t2.35E-02\t3.38E-03\t1.25E-02\t1.28E-04\t2.81E-04\t7.14E-04\t1.23E-02\t8.53E-05\t2.49E-04\t5.63E-03\t5.45E-05\t4.34E-05\t1.85E-04\t3.05E-03\t3.69E-01\t3.21E-04\t1.74E-04\t2.85E-02\t2.83E-04\t1.50E-03\t3.18E-03\t1.61E-03\t1.12E-03\t6.87E-04\t3.88E-04\t4.66E-04\t4.07E-03\t1.64E-02\t2.21E-03\t2.13E-03\t1.43E-04\t3.36E-04\t1.09E-04\t4.24E-03\t6.92E-04\t2.01E-02\t2.15E-03\t2.12E-02\t7.59E-03\t2.55E-03\t3.95E-03\t5.07E-02\t2.37E-02\t4.46E-02\t3.81E-04\t2.51E-04\t2.59E-03\t3.29E-03\t4.90E-03\t9.24E-04\t1.48E-04\t5.44E-04\t4.21E-04\t3.33E-03\t2.42E-03\t3.74E-03\t3.11E-04\t2.11E-03\t4.11E-04\t1.29E-02\t1.81E-02\t1.52E-03\t7.61E-04\t0.004953383455\t0.006244142553\r",
 "C>T\tACA\t2.50E-02\t6.11E-05\t1.42E-02\t8.70E-03\t3.26E-02\t6.41E-02\t1.22E-04\t1.41E-03\t1.33E-03\t3.33E-03\t6.81E-03\t8.10E-03\t1.04E-03\t5.49E-04\t2.58E-02\t5.51E-03\t7.18E-03\t1.24E-02\t5.05E-03\t1.76E-03\t1.22E-03\t5.71E-04\t9.20E-03\t2.80E-02\t4.93E-02\t2.17E-03\t1.76E-03\t2.26E-02\t1.18E-02\t2.10E-02\t2.66E-03\t1.37E-02\t7.77E-04\t5.09E-03\t1.06E-01\t1.77E-02\t8.64E-02\t5.31E-03\t3.52E-03\t2.21E-03\t5.96E-03\t6.24E-03\t3.45E-03\t9.60E-03\t1.97E-02\t9.88E-03\t2.24E-02\t2.90E-05\t8.99E-02\t1.33E-03\t1.52E-02\t4.27E-03\t2.87E-04\t5.75E-04\t3.52E-03\t3.58E-02\t2.69E-03\t7.22E-03\t1.48E-03\t5.11E-03\t7.40E-03\t3.26E-02\t8.05E-02\t3.92E-03\t2.76E-03\t0.053885768706\t0.001374780580\r",
 "C>T\tACC\t6.32E-03\t1.38E-03\t1.24E-02\t4.18E-03\t1.79E-02\t2.29E-02\t1.77E-02\t3.81E-02\t7.68E-03\t1.40E-03\t2.82E-03\t5.42E-03\t8.34E-04\t5.25E-03\t1.47E-01\t2.12E-03\t2.22E-16\t1.01E-02\t6.91E-03\t3.73E-03\t1.71E-03\t7.03E-05\t4.65E-03\t1.41E-02\t1.30E-02\t3.92E-03\t2.20E-04\t8.09E-02\t1.12E-02\t1.31E-02\t2.54E-03\t5.42E-03\t3.38E-03\t4.31E-03\t9.04E-02\t9.28E-03\t1.23E-01\t2.13E-03\t3.71E-03\t3.99E-03\t3.28E-03\t3.06E-03\t2.15E-03\t4.04E-03\t1.21E-02\t3.36E-03\t3.98E-02\t1.28E-04\t1.69E-02\t6.01E-04\t8.90E-03\t9.40E-04\t4.38E-04\t3.03E-04\t1.80E-03\t5.65E-03\t2.10E-03\t1.96E-03\t5.29E-04\t2.33E-04\t9.82E-04\t1.29E-02\t6.76E-03\t2.86E-03\t1.34E-03\t0.042106831614\t0.001478441742\r",
 "C>T\tACG\t3.65E-01\t3.27E-05\t2.57E-03\t7.83E-04\t6.18E-03\t1.20E-01\t2.23E-16\t5.43E-04\t9.27E-04\t1.21E-03\t3.23E-03\t4.38E-03\t5.57E-03\t4.14E-04\t7.58E-03\t1.18E-04\t3.03E-04\t2.65E-02\t7.45E-02\t3.27E-03\t3.95E-05\t2.87E-04\t1.22E-02\t1.55E-03\t1.33E-02\t7.40E-03\t1.36E-03\t8.72E-03\t5.14E-04\t1.97E-02\t5.33E-03\t6.71E-03\t7.05E-03\t2.01E-02\t1.44E-02\t7.59E-03\t3.66E-02\t7.16E-03\t1.18E-04\t5.80E-04\t7.07E-03\t6.33E-03\t5.18E-03\t6.00E-03\t2.67E-04\t3.43E-03\t4.40E-04\t1.22E-03\t1.43E-02\t5.41E-03\t4.60E-03\t3.05E-03\t1.21E-03\t2.01E-03\t1.39E-02\t5.65E-03\t4.08E-03\t6.95E-03\t3.43E-04\t5.32E-03\t1.52E-03\t1.76E-02\t8.63E-03\t2.56E-03\t4.21E-03\t0.002540673114\t0.000443775624\r",
 "C>T\tACT\t9.58E-03\t1.86E-03\t1.21E-02\t4.25E-03\t2.20E-02\t2.13E-02\t7.42E-03\t1.08E-02\t7.48E-03\t3.29E-02\t1.25E-02\t6.49E-03\t2.99E-03\t7.93E-03\t1.09E-01\t2.38E-03\t3.18E-04\t1.61E-02\t1.78E-03\t6.00E-03\t7.37E-04\t1.47E-04\t7.97E-03\t3.68E-02\t1.30E-02\t4.17E-05\t6.16E-04\t3.50E-02\t9.15E-03\t5.73E-03\t2.79E-03\t4.48E-03\t5.06E-03\t1.83E-03\t2.93E-02\t1.01E-02\t1.52E-01\t3.07E-03\t3.87E-03\t9.00E-04\t4.31E-03\t5.65E-03\t2.44E-03\t1.28E-02\t1.69E-02\t6.35E-03\t3.52E-02\t9.10E-04\t3.03E-02\t1.22E-03\t6.80E-03\t2.49E-03\t2.35E-04\t4.15E-04\t6.92E-03\t3.20E-03\t1.76E-03\t1.94E-03\t4.99E-03\t7.76E-04\t2.84E-03\t1.54E-02\t7.31E-03\t3.39E-03\t2.98E-04\t0.054345874146\t0.004563583260\r",
 "C>T\tCCA\t2.00E-03\t4.32E-03\t1.61E-02\t1.27E-02\t1.94E-02\t1.36E-02\t5.00E-02\t1.07E-01\t2.71E-03\t3.56E-02\t2.67E-03\t7.24E-03\t3.68E-05\t1.51E-04\t4.07E-03\t6.59E-03\t1.93E-05\t2.73E-03\t1.08E-02\t7.17E-03\t8.21E-04\t4.57E-04\t1.17E-02\t1.19E-01\t1.06E-02\t1.35E-03\t1.21E-03\t6.22E-02\t5.38E-03\t1.43E-02\t1.69E-03\t1.74E-02\t2.92E-04\t7.61E-03\t1.10E-01\t4.12E-02\t1.69E-02\t4.59E-03\t4.75E-03\t5.10E-03\t6.75E-03\t3.61E-03\t4.24E-03\t9.11E-03\t1.96E-02\t6.91E-03\t1.87E-02\t7.62E-06\t3.57E-03\t2.14E-03\t2.84E-02\t8.64E-04\t1.82E-05\t4.70E-04\t3.84E-03\t3.23E-03\t3.27E-03\t1.46E-02\t3.37E-03\t7.25E-03\t2.98E-04\t1.49E-02\t8.09E-03\t2.87E-03\t1.52E-03\t0.018028129532\t0.001636999066\r",
 "C>T\tCCC\t2.70E-04\t1.70E-04\t2.02E-02\t1.15E-02\t1.90E-02\t1.48E-02\t5.53E-02\t1.82E-01\t5.97E-03\t2.60E-03\t4.02E-03\t8.76E-03\t2.23E-16\t8.44E-05\t1.26E-01\t2.54E-04\t1.35E-03\t1.38E-03\t1.11E-02\t2.84E-03\t1.72E-03\t5.89E-05\t6.76E-03\t9.20E-02\t1.30E-02\t5.16E-03\t3.19E-03\t1.89E-01\t8.85E-03\t5.94E-03\t2.83E-03\t7.14E-03\t2.10E-04\t1.64E-03\t1.03E-01\t1.64E-01\t2.84E-02\t1.53E-03\t5.14E-04\t8.32E-02\t6.68E-03\t2.10E-03\t4.73E-03\t3.86E-03\t1.73E-02\t1.85E-03\t4.81E-02\t4.35E-18\t2.86E-03\t5.07E-05\t1.33E-02\t4.61E-03\t3.01E-04\t4.93E-04\t4.25E-03\t4.70E-03\t1.82E-03\t4.60E-03\t2.26E-04\t4.00E-03\t5.31E-04\t1.21E-02\t9.92E-03\t4.09E-04\t1.88E-03\t0.040621303139\t0.000369373076\r",
 "C>T\tCCG\t1.96E-01\t2.77E-03\t2.24E-03\t3.31E-03\t1.75E-02\t1.22E-01\t6.71E-03\t2.21E-02\t8.94E-03\t2.64E-03\t1.49E-03\t6.51E-03\t2.23E-16\t3.70E-03\t9.20E-04\t4.09E-05\t1.99E-03\t1.16E-02\t2.48E-02\t6.91E-04\t1.90E-03\t5.43E-04\t1.68E-02\t1.57E-02\t1.65E-02\t6.22E-03\t4.26E-04\t1.42E-02\t6.35E-03\t1.16E-02\t3.34E-03\t7.50E-03\t2.55E-03\t1.53E-02\t1.68E-02\t8.20E-03\t1.12E-02\t5.88E-03\t9.45E-03\t2.67E-03\t5.33E-03\t3.27E-03\t3.50E-03\t7.57E-03\t2.36E-03\t1.01E-02\t1.92E-03\t9.72E-05\t1.55E-02\t2.20E-03\t8.30E-03\t1.52E-03\t1.87E-05\t1.38E-03\t1.08E-03\t2.07E-02\t2.38E-03\t1.37E-02\t2.83E-02\t6.59E-05\t7.95E-03\t6.40E-03\t7.99E-03\t1.10E-03\t4.06E-03\t0.000000000000\t0.000024057242\r",
 "C>T\tCCT\t1.96E-04\t1.45E-03\t2.30E-02\t1.60E-02\t2.31E-02\t1.38E-02\t5.66E-02\t1.26E-01\t2.06E-02\t1.38E-04\t1.45E-02\t1.11E-02\t1.03E-04\t5.95E-04\t6.51E-02\t3.38E-03\t1.45E-03\t1.36E-03\t1.56E-02\t1.93E-02\t1.74E-03\t1.36E-04\t6.65E-03\t2.58E-01\t1.17E-02\t3.25E-04\t1.60E-03\t1.56E-01\t4.29E-03\t1.36E-02\t1.65E-03\t7.15E-03\t3.32E-03\t6.73E-04\t4.62E-02\t1.45E-01\t4.28E-02\t2.70E-03\t3.53E-05\t4.70E-02\t5.75E-03\t6.06E-03\t4.10E-03\t6.99E-03\t2.16E-02\t7.73E-03\t3.29E-02\t7.44E-04\t5.04E-03\t2.73E-03\t6.40E-03\t3.00E-03\t4.79E-04\t5.19E-04\t1.19E-02\t4.06E-03\t2.57E-03\t1.60E-02\t3.71E-03\t7.19E-04\t5.08E-04\t3.20E-02\t2.77E-02\t2.42E-03\t2.51E-03\t0.042237294613\t0.002350306147\r",
 "C>T\tGCA\t4.44E-03\t2.23E-16\t1.64E-02\t4.09E-03\t1.94E-02\t9.36E-02\t8.49E-05\t8.55E-04\t1.57E-04\t3.65E-04\t2.76E-04\t6.98E-03\t4.22E-05\t3.63E-04\t1.97E-02\t4.48E-03\t7.18E-03\t1.71E-02\t7.37E-02\t6.50E-03\t1.64E-03\t6.11E-05\t9.68E-03\t6.10E-02\t4.22E-02\t9.32E-03\t2.38E-03\t5.37E-02\t2.99E-02\t1.00E-02\t2.75E-03\t1.42E-02\t3.25E-03\t1.04E-02\t4.97E-02\t6.80E-03\t2.26E-02\t5.35E-03\t2.33E-04\t7.29E-03\t3.72E-03\t2.55E-03\t2.37E-03\t8.03E-03\t1.43E-02\t6.01E-03\t3.36E-02\t7.66E-04\t9.55E-02\t1.00E-03\t1.26E-02\t1.06E-03\t3.45E-04\t1.84E-04\t2.42E-03\t4.47E-03\t1.78E-03\t9.18E-03\t2.22E-04\t6.50E-03\t9.49E-06\t5.14E-03\t6.78E-03\t2.44E-03\t1.24E-03\t0.089177561445\t0.000887033367\r",
 "C>T\tGCC\t9.28E-05\t1.29E-07\t1.34E-02\t4.46E-03\t2.16E-02\t1.11E-01\t1.17E-02\t2.63E-02\t9.51E-06\t1.81E-04\t2.45E-04\t6.53E-03\t1.25E-02\t1.51E-02\t1.15E-01\t1.82E-03\t9.21E-04\t8.16E-03\t1.29E-01\t5.29E-03\t2.22E-16\t1.57E-04\t8.17E-03\t3.97E-02\t5.32E-02\t1.96E-02\t1.87E-03\t1.54E-01\t4.37E-02\t6.95E-03\t5.45E-03\t2.01E-02\t1.88E-04\t7.96E-03\t6.10E-02\t1.14E-02\t3.79E-02\t5.14E-03\t2.22E-16\t1.78E-02\t3.40E-03\t2.31E-03\t2.96E-03\t4.50E-03\t1.16E-02\t8.64E-03\t1.52E-01\t4.10E-18\t8.14E-02\t9.33E-04\t4.50E-03\t1.71E-03\t7.39E-06\t1.18E-04\t1.58E-03\t3.09E-03\t1.22E-03\t1.19E-02\t1.31E-02\t1.51E-03\t5.06E-04\t6.93E-03\t3.15E-03\t2.35E-03\t1.44E-03\t0.111629551116\t0.000484366313\r",
 "C>T\tGCG\t2.18E-01\t2.23E-16\t5.28E-04\t4.82E-04\t1.23E-02\t1.78E-01\t5.52E-06\t1.76E-04\t3.43E-05\t1.57E-04\t1.91E-03\t4.76E-03\t1.49E-02\t1.26E-02\t1.78E-03\t7.46E-05\t1.78E-04\t1.92E-05\t2.77E-01\t1.17E-03\t1.61E-03\t8.55E-04\t1.48E-02\t3.21E-03\t2.05E-02\t4.34E-05\t2.22E-16\t1.23E-02\t2.54E-02\t1.05E-02\t6.68E-03\t9.37E-03\t2.15E-04\t1.83E-02\t5.24E-03\t1.97E-03\t1.04E-02\t5.84E-03\t1.67E-02\t2.53E-03\t5.70E-03\t4.32E-03\t4.48E-03\t5.48E-03\t3.93E-18\t1.71E-02\t4.55E-03\t9.66E-06\t4.96E-02\t1.51E-03\t1.50E-03\t1.72E-03\t3.03E-04\t1.04E-03\t5.50E-03\t4.28E-03\t2.83E-03\t5.99E-03\t1.91E-02\t1.19E-02\t2.06E-03\t8.44E-03\t7.36E-03\t3.24E-18\t3.36E-03\t0.009223925314\t0.000116423626\r",
 "C>T\tGCT\t3.84E-05\t5.04E-05\t1.13E-02\t3.81E-03\t1.86E-02\t7.44E-02\t1.61E-02\t5.39E-03\t7.58E-04\t4.21E-02\t4.61E-03\t1.07E-02\t1.35E-02\t3.19E-02\t6.96E-02\t2.67E-03\t1.27E-03\t1.05E-02\t7.72E-02\t6.72E-03\t6.60E-05\t2.23E-16\t6.69E-03\t1.09E-01\t3.31E-02\t8.05E-03\t1.76E-03\t1.07E-01\t2.03E-02\t9.82E-03\t4.98E-03\t7.59E-03\t1.02E-04\t5.12E-03\t1.77E-02\t5.70E-03\t5.32E-02\t3.93E-03\t5.17E-04\t6.03E-03\t3.63E-03\t2.32E-03\t3.20E-03\t8.61E-03\t1.27E-02\t8.41E-03\t8.65E-02\t1.86E-03\t6.66E-02\t1.99E-03\t4.00E-03\t4.04E-03\t1.10E-04\t3.08E-04\t1.63E-02\t5.95E-03\t1.02E-03\t3.10E-03\t2.48E-04\t5.29E-03\t8.30E-04\t6.18E-03\t5.03E-03\t2.96E-03\t8.94E-04\t0.180440456121\t0.004167105016\r",
 "C>T\tTCA\t1.11E-03\t5.36E-01\t5.76E-03\t3.74E-03\t2.03E-02\t1.20E-02\t2.38E-01\t6.61E-03\t1.73E-02\t2.07E-04\t5.04E-03\t4.27E-03\t3.66E-04\t2.73E-03\t9.78E-03\t5.52E-03\t6.93E-03\t3.54E-03\t2.65E-02\t6.36E-03\t5.32E-04\t1.08E-04\t1.14E-02\t4.39E-02\t8.07E-03\t3.60E-03\t1.59E-03\t1.53E-02\t2.40E-03\t1.11E-02\t3.29E-03\t1.88E-02\t4.84E-03\t6.39E-03\t9.80E-02\t7.48E-03\t2.48E-02\t1.47E-02\t1.86E-03\t2.96E-03\t1.07E-02\t4.02E-03\t6.45E-03\t3.88E-03\t3.00E-02\t1.17E-02\t4.82E-03\t3.72E-05\t1.40E-03\t2.40E-03\t1.85E-02\t7.24E-03\t2.14E-04\t1.87E-04\t1.94E-02\t1.70E-03\t1.47E-03\t5.73E-03\t7.45E-18\t5.82E-11\t2.00E-04\t8.14E-03\t4.29E-03\t1.61E-04\t4.07E-05\t0.000070752091\t0.003823610731\r",
 "C>T\tTCC\t3.73E-05\t9.73E-02\t1.33E-02\t6.99E-03\t2.30E-02\t1.05E-02\t3.31E-01\t2.11E-01\t1.37E-02\t5.65E-03\t2.52E-04\t7.33E-03\t2.76E-03\t4.85E-02\t1.58E-01\t2.64E-03\t4.48E-03\t1.37E-04\t9.56E-03\t3.68E-03\t9.06E-05\t6.04E-04\t1.14E-02\t2.50E-02\t9.19E-03\t8.07E-05\t1.32E-03\t4.36E-02\t1.22E-02\t6.45E-03\t2.53E-03\t7.99E-03\t4.13E-03\t3.61E-18\t1.34E-01\t2.46E-02\t3.87E-02\t3.17E-03\t2.89E-04\t1.26E-02\t6.86E-03\t2.50E-03\t7.23E-03\t5.28E-03\t2.31E-02\t6.47E-03\t3.27E-02\t5.06E-18\t2.25E-03\t2.02E-03\t1.26E-02\t3.71E-02\t3.17E-04\t5.82E-04\t1.04E-01\t1.77E-03\t1.53E-03\t5.77E-03\t6.61E-04\t4.11E-03\t1.84E-04\t1.73E-02\t1.53E-02\t2.80E-04\t9.58E-04\t0.026909334242\t0.002998120759\r",
 "C>T\tTCG\t1.10E-01\t4.42E-02\t9.89E-04\t2.81E-04\t1.59E-02\t5.20E-02\t7.41E-02\t5.74E-02\t4.64E-05\t1.75E-03\t4.17E-05\t3.26E-03\t2.79E-03\t4.37E-01\t2.94E-07\t8.85E-05\t1.19E-03\t5.32E-03\t2.12E-02\t3.02E-03\t1.77E-04\t6.24E-05\t1.18E-02\t4.32E-03\t2.64E-03\t3.17E-04\t3.77E-05\t2.87E-03\t2.69E-03\t1.41E-02\t2.63E-03\t4.61E-03\t1.74E-03\t1.44E-02\t1.47E-02\t8.41E-03\t9.85E-03\t3.85E-03\t3.00E-04\t7.65E-04\t3.89E-03\t3.43E-03\t3.47E-03\t4.86E-03\t2.12E-03\t4.92E-03\t3.64E-04\t3.07E-18\t7.68E-18\t1.75E-03\t3.70E-03\t1.94E-03\t4.02E-04\t9.34E-04\t3.97E-03\t6.49E-03\t2.19E-03\t1.02E-02\t1.19E-13\t4.81E-03\t1.54E-04\t4.60E-03\t3.34E-03\t4.54E-18\t1.73E-03\t0.000006923779\t0.000157958968\r",
 "C>T\tTCT\t2.23E-16\t3.00E-01\t8.17E-03\t7.22E-03\t2.10E-02\t6.99E-03\t1.07E-01\t1.02E-01\t2.32E-03\t5.00E-02\t5.27E-03\t5.46E-03\t2.28E-02\t1.28E-01\t1.19E-01\t6.12E-03\t1.42E-03\t3.59E-04\t1.74E-03\t2.21E-02\t3.48E-04\t8.95E-04\t3.75E-03\t6.29E-02\t5.00E-05\t7.31E-05\t2.22E-16\t2.63E-02\t4.71E-03\t1.24E-02\t2.13E-03\t1.19E-02\t1.99E-03\t4.28E-03\t5.25E-02\t1.95E-02\t6.70E-02\t9.04E-03\t8.08E-04\t8.04E-03\t8.33E-03\t5.84E-03\t6.58E-03\t4.23E-03\t2.24E-02\t1.15E-02\t1.51E-02\t8.36E-04\t3.02E-03\t2.95E-03\t4.40E-03\t8.09E-02\t5.85E-06\t3.12E-04\t8.18E-02\t4.24E-02\t1.65E-03\t7.71E-03\t2.36E-04\t3.14E-02\t8.20E-05\t3.72E-02\t1.79E-01\t2.09E-03\t3.91E-03\t0.016839668838\t0.038650889984\r",
 "T>A\tATA\t8.00E-04\t9.59E-05\t5.49E-03\t9.42E-03\t7.38E-03\t9.82E-05\t1.94E-03\t2.49E-03\t1.24E-02\t8.89E-05\t1.95E-02\t1.27E-02\t5.07E-05\t1.82E-03\t4.17E-04\t5.17E-03\t2.22E-16\t6.13E-04\t2.04E-03\t2.31E-02\t1.01E-03\t6.26E-05\t2.49E-03\t7.86E-04\t5.11E-04\t3.46E-03\t6.39E-02\t8.28E-18\t1.58E-03\t2.04E-02\t1.52E-03\t1.43E-01\t1.72E-04\t6.19E-04\t3.63E-04\t7.86E-03\t7.13E-04\t8.68E-04\t1.03E-01\t4.89E-03\t2.44E-03\t2.25E-03\t1.99E-03\t5.31E-03\t6.10E-03\t7.46E-03\t1.33E-02\t1.53E-02\t4.78E-03\t5.25E-04\t2.70E-03\t1.38E-01\t6.60E-04\t1.29E-04\t5.36E-03\t1.13E-01\t1.77E-03\t5.89E-03\t9.56E-04\t4.26E-03\t5.05E-06\t3.28E-03\t2.45E-02\t2.00E-03\t3.97E-03\t0.001562338887\t0.111858717410\r",
 "T>A\tATC\t2.23E-03\t8.78E-04\t7.21E-03\t3.53E-03\t6.93E-03\t5.46E-04\t6.70E-04\t4.67E-04\t3.42E-03\t1.08E-03\t2.13E-02\t2.82E-03\t1.94E-04\t3.37E-04\t4.41E-04\t2.71E-03\t6.52E-05\t3.72E-04\t8.41E-04\t7.74E-03\t9.97E-04\t2.89E-04\t8.78E-04\t2.55E-04\t5.19E-03\t5.65E-03\t7.01E-03\t5.56E-18\t2.25E-03\t8.89E-03\t3.52E-03\t1.46E-03\t4.78E-04\t2.00E-03\t3.93E-04\t6.07E-03\t8.87E-04\t1.14E-03\t4.44E-02\t7.36E-03\t1.43E-03\t1.66E-03\t1.94E-03\t2.88E-03\t4.44E-03\t4.30E-03\t4.01E-03\t3.78E-03\t2.60E-03\t4.97E-04\t1.80E-03\t7.16E-03\t2.66E-04\t7.78E-05\t3.94E-03\t1.35E-03\t1.18E-03\t2.59E-03\t3.19E-03\t6.18E-04\t4.29E-03\t6.43E-03\t1.46E-03\t2.50E-03\t2.23E-04\t0.000397461148\t0.013051252824\r",
 "T>A\tATG\t1.14E-03\t2.55E-05\t9.64E-03\t1.24E-02\t9.12E-03\t2.22E-16\t5.87E-04\t1.56E-03\t2.11E-04\t3.70E-03\t2.28E-02\t2.76E-03\t2.23E-16\t2.34E-05\t3.93E-05\t3.12E-03\t5.90E-04\t3.54E-04\t5.40E-04\t8.81E-03\t2.00E-04\t6.57E-05\t2.60E-03\t9.52E-04\t3.03E-04\t2.48E-04\t6.06E-02\t3.54E-04\t2.00E-03\t6.24E-03\t6.25E-04\t1.23E-03\t1.52E-04\t2.48E-03\t6.38E-04\t5.03E-03\t7.12E-04\t6.96E-04\t3.46E-02\t1.33E-03\t2.53E-03\t2.49E-03\t1.08E-03\t6.94E-03\t6.93E-03\t3.03E-03\t5.01E-03\t1.43E-03\t1.08E-06\t2.13E-03\t5.00E-04\t7.45E-03\t4.65E-04\t2.41E-04\t5.73E-04\t3.68E-02\t1.19E-03\t2.87E-03\t4.23E-03\t8.10E-03\t5.13E-04\t5.80E-03\t4.69E-03\t2.19E-03\t3.52E-04\t0.000070829550\t0.021107255813\r",
 "T>A\tATT\t1.83E-04\t2.23E-16\t6.11E-03\t4.21E-03\t6.20E-03\t1.66E-03\t3.71E-03\t8.29E-03\t5.07E-02\t4.11E-02\t2.50E-02\t1.89E-02\t6.66E-03\t2.85E-03\t3.40E-06\t2.82E-03\t6.06E-04\t1.85E-03\t9.54E-03\t1.10E-02\t1.46E-03\t6.33E-04\t7.84E-03\t6.54E-04\t2.56E-02\t1.60E-03\t2.51E-03\t1.61E-04\t1.80E-03\t1.49E-02\t3.07E-03\t2.54E-03\t1.70E-03\t1.96E-03\t4.90E-04\t4.16E-03\t8.44E-04\t5.07E-03\t3.16E-01\t7.37E-03\t3.07E-03\t2.77E-03\t2.72E-03\t6.12E-03\t1.07E-02\t6.96E-02\t6.85E-03\t4.28E-03\t3.08E-02\t6.05E-04\t2.00E-03\t6.76E-02\t6.75E-04\t3.54E-04\t2.41E-03\t2.95E-02\t1.19E-03\t4.63E-03\t3.01E-04\t3.92E-03\t2.90E-02\t6.31E-03\t8.77E-03\t2.80E-03\t7.28E-04\t0.010652827921\t0.070288409275\r",
 "T>A\tCTA\t4.30E-05\t2.43E-04\t7.94E-03\t1.44E-02\t3.54E-03\t6.66E-05\t1.77E-04\t1.12E-03\t3.47E-03\t1.61E-03\t1.90E-02\t4.79E-03\t2.23E-16\t7.09E-05\t1.84E-04\t5.46E-04\t2.09E-04\t1.55E-04\t3.48E-04\t3.21E-03\t4.34E-04\t9.36E-05\t5.42E-04\t9.60E-05\t1.44E-04\t6.19E-06\t1.28E-01\t4.62E-04\t1.32E-03\t2.63E-02\t9.44E-04\t3.31E-01\t1.40E-04\t1.95E-03\t2.14E-04\t2.06E-02\t3.25E-04\t7.46E-04\t2.37E-02\t5.26E-02\t1.26E-03\t1.57E-03\t8.97E-04\t4.94E-03\t1.09E-03\t3.62E-03\t1.26E-02\t1.36E-03\t1.24E-05\t6.97E-04\t4.20E-03\t1.40E-02\t8.62E-05\t7.41E-05\t1.93E-04\t6.11E-03\t1.22E-03\t3.69E-03\t7.45E-18\t7.39E-05\t1.19E-03\t6.39E-03\t2.48E-03\t7.03E-05\t6.64E-05\t0.000136550096\t0.009896320186\r",
 "T>A\tCTC\t3.93E-04\t2.44E-04\t1.78E-02\t1.22E-02\t5.19E-03\t5.51E-04\t2.61E-04\t9.58E-04\t1.34E-02\t8.13E-05\t2.92E-02\t3.13E-03\t2.23E-16\t1.06E-04\t2.12E-04\t5.60E-03\t3.42E-04\t2.95E-03\t2.44E-04\t5.70E-03\t4.05E-02\t2.24E-03\t3.06E-03\t4.03E-04\t6.48E-03\t6.88E-04\t5.51E-02\t4.47E-18\t1.87E-03\t1.58E-02\t1.45E-03\t2.81E-03\t1.35E-03\t1.41E-03\t7.77E-04\t2.11E-02\t6.50E-04\t9.92E-04\t4.09E-04\t3.30E-02\t1.46E-03\t1.16E-03\t2.02E-03\t6.44E-03\t7.62E-03\t3.19E-03\t1.72E-03\t1.21E-05\t3.25E-03\t1.04E-03\t3.30E-03\t5.71E-04\t1.02E-04\t5.14E-05\t5.35E-04\t1.49E-03\t2.31E-04\t4.50E-03\t7.32E-04\t1.10E-06\t9.74E-04\t6.59E-03\t2.23E-02\t1.64E-03\t5.14E-04\t0.000000000000\t0.001136017362\r",
 "T>A\tCTG\t3.24E-04\t1.49E-04\t1.32E-02\t3.25E-02\t5.01E-03\t1.15E-04\t1.12E-04\t8.83E-04\t9.61E-04\t1.20E-03\t2.76E-02\t5.31E-03\t2.23E-16\t2.22E-16\t3.05E-04\t1.41E-03\t2.26E-04\t3.17E-04\t8.20E-05\t1.03E-03\t1.24E-02\t7.17E-04\t3.21E-03\t2.94E-04\t1.23E-04\t1.25E-04\t1.68E-01\t1.32E-04\t5.03E-03\t6.29E-02\t1.14E-03\t4.40E-03\t8.06E-04\t1.98E-03\t5.87E-04\t1.87E-02\t1.41E-03\t4.79E-04\t1.17E-03\t5.82E-02\t2.14E-03\t1.60E-03\t1.78E-03\t6.80E-03\t3.46E-03\t3.12E-03\t1.49E-02\t2.96E-04\t5.82E-04\t3.41E-03\t4.30E-03\t5.82E-05\t3.10E-04\t2.32E-04\t5.16E-04\t6.52E-03\t1.23E-03\t3.30E-03\t1.92E-03\t1.60E-02\t1.59E-03\t1.29E-03\t8.58E-04\t1.39E-04\t3.81E-07\t0.000125175761\t0.000178582568\r",
 "T>A\tCTT\t2.60E-04\t2.23E-16\t1.27E-02\t1.18E-02\t4.30E-03\t1.79E-04\t3.71E-04\t1.12E-04\t4.87E-02\t1.15E-03\t4.18E-02\t4.44E-03\t3.76E-04\t2.34E-04\t1.62E-04\t1.50E-03\t5.23E-04\t1.38E-03\t1.05E-03\t6.53E-03\t7.23E-02\t7.79E-03\t2.26E-03\t2.18E-04\t1.02E-04\t2.22E-16\t6.03E-02\t3.75E-04\t4.80E-03\t3.66E-02\t6.42E-04\t5.43E-03\t1.03E-03\t1.98E-03\t6.12E-04\t3.10E-02\t6.63E-04\t1.12E-03\t2.68E-04\t6.64E-02\t1.96E-03\t2.06E-03\t2.36E-03\t8.27E-03\t9.28E-03\t9.01E-03\t7.57E-03\t2.03E-04\t1.38E-05\t1.34E-03\t3.60E-03\t7.60E-05\t3.62E-04\t1.06E-04\t2.09E-04\t3.62E-03\t8.56E-04\t8.61E-03\t6.66E-05\t5.23E-03\t6.09E-06\t1.07E-02\t1.30E-03\t3.71E-03\t1.32E-04\t0.000000000000\t0.010551155388\r",
 "T>A\tGTA\t8.12E-05\t2.23E-16\t7.06E-03\t9.22E-03\t3.67E-03\t5.10E-05\t4.01E-04\t5.26E-04\t1.73E-03\t2.37E-05\t1.28E-02\t2.65E-03\t3.31E-06\t2.22E-16\t1.54E-06\t8.26E-04\t6.70E-06\t5.77E-05\t5.11E-04\t9.40E-04\t6.46E-05\t9.73E-05\t1.60E-03\t1.89E-04\t5.56E-05\t1.17E-03\t5.37E-02\t6.26E-18\t1.57E-03\t3.62E-02\t7.83E-04\t1.10E-01\t1.16E-04\t3.43E-03\t3.29E-04\t1.90E-03\t8.79E-03\t8.10E-04\t1.90E-02\t7.17E-04\t8.80E-04\t1.15E-03\t8.93E-04\t4.96E-03\t2.65E-03\t2.33E-03\t1.20E-02\t6.18E-03\t1.16E-03\t3.00E-04\t6.10E-03\t1.60E-02\t1.30E-04\t5.04E-05\t7.32E-04\t2.04E-02\t5.24E-04\t4.45E-03\t2.33E-03\t2.43E-03\t4.20E-04\t3.33E-03\t2.07E-03\t2.55E-04\t1.09E-04\t0.001098638590\t0.017239232975\r",
 "T>A\tGTC\t1.29E-04\t1.31E-04\t7.48E-03\t3.60E-03\t2.91E-03\t4.89E-04\t2.24E-04\t2.56E-04\t2.27E-03\t1.12E-04\t1.35E-02\t2.12E-03\t2.23E-16\t1.23E-04\t3.75E-04\t1.81E-03\t1.30E-04\t2.73E-03\t5.89E-04\t6.08E-04\t3.40E-03\t3.47E-04\t1.95E-03\t7.78E-05\t2.56E-03\t3.26E-03\t1.02E-02\t1.40E-04\t2.35E-03\t5.85E-03\t1.15E-03\t2.35E-05\t1.18E-04\t7.89E-04\t1.82E-04\t2.91E-03\t8.22E-03\t3.62E-04\t2.42E-03\t5.56E-03\t6.31E-04\t7.74E-04\t1.12E-03\t3.52E-03\t3.43E-03\t2.29E-03\t1.84E-03\t1.33E-03\t2.36E-03\t6.38E-04\t1.00E-03\t2.28E-04\t1.24E-04\t9.10E-05\t1.05E-03\t3.85E-03\t3.99E-04\t2.70E-03\t8.99E-03\t1.85E-04\t2.37E-03\t5.97E-04\t6.26E-04\t1.98E-03\t1.10E-04\t0.000065576983\t0.001005109264\r",
 "T>A\tGTG\t2.46E-04\t6.38E-05\t1.16E-02\t1.31E-02\t4.11E-03\t1.95E-04\t2.90E-04\t6.56E-04\t2.22E-16\t1.35E-03\t1.56E-02\t7.88E-04\t7.30E-05\t2.22E-16\t1.19E-04\t2.51E-03\t4.70E-05\t1.20E-04\t1.87E-04\t5.55E-03\t6.95E-04\t2.41E-04\t1.83E-03\t4.86E-04\t2.13E-03\t7.42E-04\t3.75E-02\t4.68E-18\t2.49E-03\t3.82E-02\t6.08E-04\t1.57E-03\t7.17E-05\t3.66E-03\t9.86E-05\t5.77E-03\t1.31E-02\t8.27E-04\t6.53E-03\t1.12E-03\t1.17E-03\t1.08E-03\t1.21E-03\t6.67E-03\t4.41E-03\t1.43E-03\t7.77E-03\t4.50E-03\t2.21E-03\t1.03E-03\t2.40E-03\t1.31E-05\t1.43E-04\t1.10E-04\t1.99E-04\t4.29E-03\t1.29E-03\t3.96E-03\t2.64E-03\t1.25E-02\t4.36E-04\t2.46E-03\t2.04E-02\t5.05E-04\t3.60E-04\t0.000337380913\t0.003520166086\r",
 "T>A\tGTT\t2.58E-04\t6.28E-05\t1.23E-02\t5.52E-03\t3.22E-03\t2.48E-04\t3.06E-04\t5.59E-05\t2.03E-02\t8.76E-03\t2.81E-02\t2.56E-03\t4.26E-04\t1.08E-04\t1.64E-04\t9.15E-04\t3.58E-04\t8.57E-04\t8.40E-04\t9.66E-04\t4.80E-03\t1.58E-03\t3.46E-03\t2.83E-04\t4.30E-03\t9.37E-04\t9.32E-03\t5.78E-18\t1.31E-03\t2.02E-02\t8.34E-04\t0.00E+00\t3.24E-04\t1.07E-03\t1.91E-04\t7.89E-04\t6.16E-03\t1.76E-03\t7.48E-03\t2.86E-03\t1.12E-03\t1.68E-03\t1.84E-03\t4.09E-03\t5.77E-03\t4.52E-03\t3.25E-03\t1.31E-03\t3.18E-03\t5.10E-04\t1.20E-03\t3.28E-06\t2.31E-04\t1.31E-04\t7.86E-04\t1.46E-04\t8.11E-04\t3.11E-03\t1.85E-03\t5.12E-03\t5.37E-03\t7.93E-03\t2.77E-03\t1.51E-03\t6.37E-04\t0.000135756563\t0.011658501206\r",
 "T>A\tTTA\t6.72E-03\t4.22E-04\t6.51E-03\t8.88E-03\t6.31E-03\t4.54E-04\t3.56E-04\t3.24E-03\t1.54E-02\t1.67E-02\t1.15E-02\t2.60E-02\t2.04E-03\t3.06E-03\t1.86E-04\t3.75E-03\t2.96E-04\t2.00E-04\t1.09E-03\t1.23E-02\t1.19E-04\t1.19E-03\t1.08E-02\t4.94E-04\t1.20E-02\t2.57E-04\t7.18E-02\t7.77E-18\t4.12E-03\t2.18E-02\t1.13E-03\t5.28E-02\t1.88E-04\t7.64E-04\t2.82E-04\t6.70E-03\t5.12E-04\t4.10E-03\t2.72E-01\t4.30E-03\t2.46E-03\t2.45E-03\t2.44E-03\t6.69E-03\t7.13E-03\t8.06E-02\t8.50E-03\t7.63E-05\t6.96E-03\t4.07E-04\t4.50E-03\t1.56E-01\t2.35E-05\t4.20E-04\t1.69E-03\t2.01E-02\t1.05E-03\t3.10E-03\t1.55E-03\t3.61E-03\t9.42E-03\t3.02E-03\t1.19E-02\t2.09E-03\t1.13E-03\t0.006569907486\t0.060179901809\r",
 "T>A\tTTC\t2.23E-16\t1.92E-04\t1.05E-02\t4.31E-03\t4.14E-03\t1.14E-04\t3.89E-04\t2.24E-03\t3.42E-02\t1.23E-03\t1.95E-02\t1.85E-03\t3.60E-04\t6.65E-04\t2.02E-04\t6.45E-03\t2.32E-04\t1.73E-03\t1.57E-04\t4.67E-03\t5.70E-03\t2.42E-04\t2.73E-03\t4.75E-04\t3.18E-04\t7.95E-04\t2.36E-02\t5.27E-18\t3.27E-03\t3.41E-03\t1.26E-03\t4.71E-05\t1.25E-03\t1.93E-03\t2.41E-04\t3.75E-03\t2.45E-03\t2.17E-04\t2.22E-16\t5.30E-03\t1.02E-03\t1.13E-03\t9.28E-04\t3.69E-03\t1.03E-02\t1.43E-02\t2.26E-03\t9.38E-04\t2.07E-05\t2.62E-04\t4.30E-03\t1.23E-05\t7.64E-06\t5.70E-05\t6.64E-04\t2.58E-03\t4.69E-04\t3.56E-03\t1.79E-04\t6.95E-05\t8.05E-04\t9.23E-03\t8.98E-04\t1.86E-03\t2.23E-16\t0.000864507132\t0.007794895232\r",
 "T>A\tTTG\t2.80E-05\t5.01E-05\t6.71E-03\t9.54E-03\t4.79E-03\t2.95E-05\t2.23E-16\t1.70E-03\t1.22E-02\t5.60E-03\t1.29E-02\t1.67E-02\t6.93E-04\t7.06E-05\t2.01E-04\t2.77E-03\t1.88E-04\t3.57E-04\t1.49E-04\t5.18E-03\t1.83E-03\t2.23E-16\t1.79E-03\t3.73E-04\t1.76E-04\t9.85E-05\t7.32E-02\t1.85E-04\t1.99E-03\t1.02E-02\t2.78E-04\t0.00E+00\t3.21E-04\t9.24E-04\t1.90E-04\t3.27E-03\t1.76E-03\t6.60E-04\t7.77E-03\t1.38E-03\t8.72E-04\t1.28E-03\t8.04E-04\t5.46E-03\t5.02E-03\t9.29E-03\t5.33E-03\t5.14E-05\t5.02E-06\t6.54E-04\t5.00E-04\t4.52E-04\t8.32E-06\t7.53E-05\t7.54E-05\t4.03E-03\t7.78E-04\t1.60E-03\t3.29E-04\t9.48E-03\t6.36E-03\t4.40E-03\t8.29E-04\t4.53E-04\t4.24E-05\t0.000784032307\t0.006733400083\r",
 "T>A\tTTT\t2.25E-03\t2.05E-04\t1.30E-02\t5.62E-03\t8.75E-03\t4.26E-04\t7.98E-04\t2.22E-16\t2.79E-01\t2.73E-02\t3.79E-02\t1.70E-02\t5.29E-03\t1.52E-03\t6.74E-04\t5.25E-03\t6.20E-05\t1.82E-03\t4.75E-04\t7.16E-03\t1.98E-03\t1.31E-03\t7.10E-03\t6.73E-04\t6.71E-03\t7.64E-05\t3.01E-02\t7.32E-18\t2.19E-03\t1.45E-02\t1.31E-03\t1.93E-03\t5.72E-03\t1.20E-03\t7.01E-04\t8.88E-03\t7.11E-04\t2.72E-03\t4.60E-02\t5.41E-03\t2.70E-03\t1.95E-03\t2.50E-03\t1.07E-02\t2.10E-02\t5.35E-02\t3.52E-03\t8.77E-18\t1.39E-02\t1.39E-03\t3.40E-03\t4.35E-03\t4.26E-04\t3.31E-04\t2.95E-03\t6.38E-03\t1.34E-03\t5.08E-03\t7.45E-18\t2.10E-03\t7.44E-04\t2.89E-02\t1.46E-02\t2.99E-03\t5.30E-04\t0.007059955557\t0.062798279518\r",
 "T>C\tATA\t1.09E-03\t6.17E-05\t1.65E-02\t7.95E-03\t4.62E-02\t2.21E-03\t4.90E-04\t1.02E-03\t2.88E-04\t4.36E-03\t9.77E-03\t2.03E-02\t7.45E-04\t9.35E-03\t1.14E-03\t7.36E-02\t1.91E-03\t1.26E-02\t3.89E-03\t2.42E-01\t9.65E-05\t9.73E-05\t2.02E-03\t2.59E-03\t1.07E-02\t1.67E-02\t2.50E-02\t8.86E-04\t6.40E-03\t2.07E-02\t8.48E-02\t5.05E-03\t5.49E-04\t1.96E-03\t6.11E-04\t2.04E-02\t9.36E-04\t1.25E-02\t4.62E-03\t2.04E-03\t5.25E-03\t7.33E-02\t2.79E-03\t1.12E-02\t1.87E-02\t1.44E-02\t1.16E-02\t8.88E-06\t1.93E-02\t1.37E-03\t1.85E-02\t2.83E-03\t1.10E-03\t6.13E-04\t2.63E-03\t2.59E-03\t3.75E-03\t1.13E-02\t2.92E-02\t4.35E-03\t1.25E-02\t2.66E-02\t8.05E-02\t9.95E-03\t1.29E-03\t0.000000000000\t0.029567020575\r",
 "T>C\tATC\t3.04E-03\t2.15E-05\t7.76E-03\t7.15E-04\t1.34E-02\t1.13E-03\t2.27E-04\t4.56E-04\t5.33E-04\t3.93E-03\t5.01E-03\t1.05E-02\t3.62E-03\t5.10E-04\t9.38E-04\t4.03E-02\t1.42E-04\t6.30E-03\t1.96E-03\t1.29E-02\t1.27E-02\t2.52E-04\t2.09E-03\t5.94E-04\t4.55E-03\t1.00E-02\t8.99E-04\t1.78E-04\t2.18E-03\t9.50E-03\t2.61E-02\t1.85E-03\t4.41E-04\t2.29E-03\t4.58E-04\t7.40E-03\t7.01E-04\t3.75E-03\t1.23E-03\t8.21E-04\t1.69E-03\t3.84E-02\t1.64E-03\t5.64E-03\t7.34E-03\t5.79E-03\t3.51E-03\t3.13E-02\t3.24E-03\t4.18E-04\t4.78E-02\t5.67E-03\t1.89E-04\t8.14E-05\t9.28E-03\t1.93E-02\t1.88E-03\t3.70E-03\t7.75E-03\t3.36E-06\t6.16E-03\t9.66E-03\t4.93E-03\t2.61E-03\t2.77E-04\t0.000323133451\t0.020410486307\r",
 "T>C\tATG\t1.06E-04\t1.32E-05\t1.23E-02\t5.00E-03\t3.82E-02\t4.24E-03\t1.87E-04\t1.50E-03\t1.07E-03\t9.19E-03\t6.98E-03\t1.43E-02\t2.14E-03\t5.31E-04\t3.20E-04\t7.05E-02\t4.28E-04\t1.24E-02\t2.40E-04\t6.39E-02\t3.76E-04\t2.23E-16\t2.08E-03\t2.31E-03\t1.86E-02\t3.88E-03\t7.75E-03\t4.70E-04\t1.76E-03\t1.99E-02\t7.85E-02\t5.69E-03\t4.05E-04\t2.75E-03\t2.77E-04\t8.12E-03\t1.09E-03\t7.34E-02\t1.95E-05\t9.20E-04\t3.84E-03\t1.68E-02\t2.17E-03\t9.64E-03\t1.55E-02\t8.04E-03\t8.56E-03\t5.86E-04\t2.72E-02\t2.07E-03\t3.29E-02\t2.09E-03\t2.85E-04\t2.00E-04\t1.13E-03\t3.95E-03\t2.62E-03\t8.25E-03\t2.02E-01\t4.02E-03\t9.55E-03\t1.03E-02\t6.75E-03\t3.52E-03\t3.52E-06\t0.000011504801\t0.005070940432\r",
 "T>C\tATT\t5.74E-03\t1.55E-04\t1.73E-02\t1.56E-03\t3.83E-02\t5.02E-04\t6.15E-04\t2.15E-03\t2.65E-03\t5.10E-02\t1.49E-02\t2.90E-02\t3.36E-03\t1.50E-04\t1.19E-03\t6.38E-02\t9.59E-04\t5.97E-03\t9.29E-04\t1.36E-01\t7.79E-03\t5.44E-04\t4.24E-03\t1.74E-03\t1.03E-04\t2.56E-02\t7.55E-04\t2.61E-04\t2.44E-03\t1.45E-02\t3.72E-02\t2.59E-03\t5.95E-04\t1.15E-03\t5.45E-04\t8.47E-03\t6.10E-04\t9.21E-03\t2.34E-03\t2.03E-03\t5.12E-03\t4.80E-02\t3.49E-03\t1.36E-02\t1.94E-02\t1.77E-02\t2.70E-03\t8.80E-03\t2.67E-04\t1.46E-03\t3.06E-02\t1.30E-03\t1.11E-03\t4.82E-04\t4.92E-03\t1.87E-03\t4.92E-03\t7.18E-03\t7.83E-03\t1.37E-02\t1.24E-02\t1.79E-02\t8.16E-03\t1.08E-02\t3.63E-05\t0.009382160052\t0.059234079950\r",
 "T>C\tCTA\t2.23E-16\t3.02E-04\t8.83E-03\t5.18E-03\t1.33E-02\t1.66E-03\t4.51E-04\t1.53E-03\t7.57E-03\t1.02E-02\t2.39E-03\t1.37E-02\t1.22E-04\t1.32E-03\t3.28E-04\t4.54E-02\t4.74E-05\t3.16E-03\t7.50E-04\t1.07E-02\t1.11E-02\t6.47E-04\t3.86E-04\t4.52E-04\t1.17E-02\t1.15E-02\t6.74E-03\t6.86E-18\t2.27E-03\t9.03E-03\t5.20E-02\t4.90E-03\t1.32E-03\t1.18E-03\t1.27E-03\t8.15E-03\t5.93E-04\t1.51E-02\t5.94E-04\t1.33E-03\t1.50E-03\t5.00E-02\t1.95E-03\t7.34E-03\t1.09E-02\t1.03E-02\t3.54E-03\t2.72E-03\t1.46E-02\t1.11E-03\t2.23E-02\t1.91E-03\t1.45E-04\t1.21E-04\t8.39E-04\t4.18E-03\t2.39E-03\t1.01E-02\t1.14E-02\t4.18E-03\t3.37E-03\t1.20E-02\t3.92E-03\t2.20E-03\t8.79E-04\t0.002269665991\t0.007184033893\r",
 "T>C\tCTC\t2.50E-03\t2.85E-05\t1.45E-02\t2.36E-03\t1.04E-02\t1.03E-03\t2.59E-04\t8.40E-04\t3.58E-03\t8.89E-05\t7.00E-03\t1.53E-02\t4.12E-05\t1.24E-04\t7.13E-04\t4.97E-02\t1.22E-04\t5.34E-03\t2.15E-03\t3.17E-03\t1.73E-01\t6.46E-03\t2.46E-03\t8.14E-05\t4.83E-03\t2.82E-03\t1.34E-03\t4.47E-18\t2.24E-03\t1.08E-02\t3.12E-02\t5.72E-03\t2.62E-03\t1.84E-03\t1.38E-03\t9.91E-03\t9.27E-04\t3.01E-02\t4.04E-05\t2.02E-03\t1.88E-03\t2.18E-02\t2.32E-03\t7.91E-03\t9.87E-03\t5.25E-03\t3.13E-03\t1.05E-01\t3.20E-03\t2.87E-03\t8.34E-02\t6.28E-03\t4.50E-04\t1.78E-04\t1.77E-02\t7.11E-02\t1.31E-03\t4.07E-03\t1.70E-02\t4.98E-03\t2.64E-03\t4.56E-03\t9.50E-03\t5.73E-03\t2.79E-03\t0.000003247249\t0.006209066076\r",
 "T>C\tCTG\t3.60E-04\t2.23E-16\t1.01E-02\t5.50E-03\t2.08E-02\t5.32E-03\t1.04E-04\t7.30E-04\t2.17E-03\t9.35E-04\t1.26E-03\t1.54E-02\t1.05E-03\t2.53E-05\t1.31E-04\t3.74E-02\t2.22E-16\t4.47E-03\t1.39E-03\t9.61E-03\t8.73E-02\t3.92E-03\t1.62E-03\t4.98E-04\t2.61E-02\t1.44E-02\t2.03E-03\t3.89E-04\t4.22E-03\t1.34E-02\t7.87E-02\t5.66E-03\t2.49E-03\t3.91E-03\t1.27E-03\t4.46E-03\t1.84E-03\t4.00E-02\t8.11E-04\t1.56E-03\t1.85E-03\t9.95E-03\t2.54E-03\t7.77E-03\t1.05E-02\t5.57E-03\t3.56E-03\t9.25E-03\t3.38E-02\t1.56E-03\t5.04E-02\t3.10E-03\t3.21E-05\t1.08E-04\t2.42E-03\t2.93E-03\t1.23E-03\t8.49E-03\t1.93E-01\t9.49E-04\t4.81E-04\t8.30E-03\t5.93E-03\t1.64E-03\t5.80E-04\t0.003628314755\t0.000839951905\r",
 "T>C\tCTT\t4.26E-05\t1.86E-04\t1.58E-02\t2.56E-03\t1.49E-02\t1.27E-03\t1.82E-03\t5.73E-03\t4.31E-02\t3.97E-02\t6.46E-03\t3.29E-02\t2.44E-06\t8.15E-04\t8.45E-05\t4.59E-02\t1.12E-04\t2.06E-03\t1.75E-03\t1.55E-02\t4.26E-01\t1.83E-02\t2.56E-03\t3.91E-04\t6.90E-03\t2.17E-02\t5.36E-04\t4.71E-05\t3.63E-03\t8.96E-03\t5.55E-02\t3.90E-03\t3.67E-04\t1.78E-03\t4.10E-03\t1.38E-02\t1.09E-03\t3.07E-03\t1.59E-03\t3.18E-03\t2.16E-03\t5.04E-02\t4.27E-03\t8.59E-03\t1.36E-02\t7.87E-03\t3.24E-03\t9.45E-03\t7.71E-03\t6.02E-04\t9.58E-02\t5.45E-03\t6.77E-05\t9.74E-05\t1.38E-03\t6.64E-03\t2.34E-03\t5.42E-03\t2.69E-02\t6.84E-04\t2.06E-04\t2.43E-02\t2.82E-02\t1.22E-02\t2.10E-03\t0.000918581801\t0.016913140039\r",
 "T>C\tGTA\t1.05E-03\t2.23E-16\t9.98E-03\t2.90E-03\t1.41E-02\t2.57E-03\t1.01E-04\t2.66E-04\t2.22E-16\t8.00E-04\t2.60E-03\t9.52E-03\t2.20E-03\t1.47E-02\t8.09E-04\t7.11E-02\t2.87E-04\t1.10E-02\t1.29E-02\t2.06E-03\t1.23E-04\t5.95E-05\t7.08E-04\t9.54E-04\t2.76E-03\t2.23E-01\t7.21E-04\t2.18E-04\t3.70E-03\t1.75E-02\t8.35E-02\t4.78E-03\t4.91E-04\t1.08E-03\t9.35E-04\t4.53E-03\t9.14E-03\t1.22E-02\t1.09E-04\t5.65E-04\t2.34E-03\t3.07E-02\t1.70E-03\t7.47E-03\t6.01E-03\t6.06E-03\t5.82E-03\t7.50E-18\t2.70E-02\t6.54E-04\t2.21E-02\t4.17E-04\t1.57E-04\t6.59E-05\t1.38E-03\t1.18E-03\t2.42E-03\t5.61E-03\t1.14E-02\t2.61E-03\t7.30E-04\t3.86E-03\t2.18E-03\t7.10E-04\t7.61E-04\t0.000000000000\t0.021590706956\r",
 "T>C\tGTC\t1.90E-03\t2.74E-05\t5.64E-03\t1.18E-03\t7.94E-03\t1.11E-03\t7.41E-05\t2.22E-16\t4.23E-04\t2.32E-03\t2.62E-03\t9.78E-03\t6.28E-03\t5.33E-03\t8.60E-04\t3.04E-02\t1.82E-04\t6.54E-03\t1.16E-02\t1.19E-03\t1.49E-02\t1.50E-03\t2.26E-03\t3.71E-04\t2.26E-03\t6.63E-02\t1.24E-03\t4.18E-04\t3.65E-03\t4.05E-03\t2.48E-02\t1.96E-03\t5.87E-04\t3.59E-03\t5.55E-04\t6.57E-03\t7.67E-03\t7.94E-03\t6.46E-04\t2.03E-03\t1.40E-03\t1.82E-02\t2.16E-03\t6.89E-03\t4.58E-03\t3.74E-03\t3.31E-03\t2.19E-02\t1.16E-03\t5.70E-04\t5.76E-02\t1.09E-03\t1.29E-04\t1.51E-04\t5.38E-03\t2.27E-03\t2.05E-03\t6.68E-03\t4.43E-03\t3.17E-03\t2.73E-05\t2.29E-03\t7.80E-04\t3.18E-03\t6.35E-04\t0.000865945066\t0.011606134984\r",
 "T>C\tGTG\t1.17E-03\t1.07E-04\t1.09E-02\t2.37E-03\t1.50E-02\t2.45E-03\t7.25E-05\t5.29E-04\t2.32E-04\t2.57E-03\t2.81E-03\t6.31E-03\t2.23E-16\t1.38E-03\t5.61E-04\t4.15E-02\t1.95E-04\t1.10E-02\t5.95E-03\t4.28E-03\t3.57E-03\t1.01E-03\t3.98E-03\t9.53E-04\t4.06E-03\t1.17E-01\t5.81E-04\t4.68E-18\t2.67E-03\t1.34E-02\t7.12E-02\t1.23E-03\t7.31E-04\t2.99E-03\t5.01E-04\t5.13E-03\t9.62E-03\t2.07E-01\t4.40E-04\t1.91E-04\t1.84E-03\t9.83E-03\t1.87E-03\t5.74E-03\t5.65E-03\t3.51E-03\t3.21E-03\t1.20E-03\t1.73E-02\t9.07E-04\t3.53E-02\t3.10E-03\t4.67E-05\t6.53E-05\t3.71E-03\t4.56E-03\t9.14E-04\t2.56E-03\t1.21E-01\t6.01E-03\t1.65E-03\t4.72E-03\t1.11E-02\t9.86E-04\t5.00E-04\t0.000197907041\t0.006462687955\r",
 "T>C\tGTT\t7.13E-05\t1.25E-05\t1.02E-02\t7.55E-04\t1.37E-02\t9.06E-04\t8.11E-04\t6.85E-03\t7.07E-03\t3.18E-01\t6.29E-03\t9.86E-03\t5.29E-03\t1.15E-02\t2.37E-04\t4.63E-02\t2.22E-16\t1.12E-02\t5.86E-03\t2.40E-03\t5.24E-03\t1.47E-03\t3.15E-03\t2.36E-03\t9.28E-04\t1.33E-01\t8.34E-04\t7.99E-04\t3.99E-03\t1.42E-02\t5.52E-02\t3.75E-03\t8.51E-04\t1.55E-03\t4.45E-04\t4.69E-03\t1.29E-02\t9.71E-03\t2.22E-16\t4.21E-04\t2.65E-03\t4.00E-02\t3.05E-03\t7.51E-03\t7.31E-03\t3.89E-03\t2.39E-03\t5.18E-04\t1.11E-02\t1.39E-03\t6.07E-02\t3.19E-03\t1.10E-04\t1.75E-04\t1.17E-03\t4.01E-03\t3.11E-03\t6.99E-03\t5.41E-03\t1.26E-02\t2.58E-03\t5.36E-03\t2.02E-03\t1.19E-03\t6.96E-04\t0.003662473198\t0.054089001127\r",
 "T>C\tTTA\t2.55E-04\t2.23E-16\t1.52E-02\t2.41E-03\t1.80E-02\t1.48E-03\t2.15E-03\t6.00E-03\t3.02E-02\t8.12E-02\t2.08E-03\t1.40E-02\t7.30E-04\t8.86E-03\t1.13E-04\t8.61E-02\t2.21E-05\t6.27E-03\t3.42E-03\t1.30E-02\t2.09E-03\t1.06E-03\t1.80E-03\t1.24E-03\t1.38E-02\t8.42E-02\t4.63E-03\t7.77E-18\t1.46E-03\t1.45E-02\t4.23E-02\t1.72E-03\t1.31E-03\t2.59E-03\t5.44E-04\t4.59E-03\t8.97E-04\t1.10E-02\t1.40E-03\t3.68E-04\t2.97E-03\t4.77E-02\t2.11E-03\t5.94E-03\t1.02E-02\t6.36E-02\t6.82E-03\t1.37E-04\t1.99E-02\t1.28E-03\t9.30E-03\t1.37E-03\t6.05E-04\t1.73E-04\t1.09E-03\t4.46E-03\t2.77E-03\t7.11E-03\t4.73E-03\t4.88E-03\t5.68E-04\t4.90E-03\t3.31E-03\t2.05E-03\t2.40E-04\t0.000432939443\t0.012963549217\r",
 "T>C\tTTC\t3.39E-03\t4.35E-05\t8.99E-03\t5.03E-04\t9.71E-03\t1.09E-03\t5.96E-04\t1.73E-03\t2.09E-02\t1.43E-02\t2.64E-03\t1.36E-02\t5.17E-03\t8.10E-05\t4.67E-04\t4.63E-02\t2.49E-05\t2.88E-03\t5.92E-03\t3.67E-03\t1.53E-02\t1.10E-03\t6.18E-03\t5.53E-04\t1.66E-02\t5.48E-02\t2.22E-16\t1.90E-04\t3.00E-03\t7.41E-03\t3.14E-02\t3.49E-03\t2.05E-03\t1.50E-03\t2.86E-04\t1.03E-02\t1.22E-03\t1.25E-02\t6.24E-03\t8.22E-04\t2.28E-03\t4.51E-02\t2.19E-03\t6.44E-03\t7.55E-03\t2.68E-02\t3.03E-03\t3.59E-02\t1.05E-02\t1.10E-03\t5.43E-02\t1.75E-02\t4.04E-04\t8.47E-05\t2.90E-02\t3.93E-03\t2.18E-03\t5.03E-03\t8.19E-03\t9.80E-04\t1.32E-03\t2.67E-02\t1.96E-02\t2.24E-03\t1.09E-03\t0.000538464673\t0.023324448028\r",
 "T>C\tTTG\t4.16E-04\t1.17E-04\t6.93E-03\t1.74E-03\t1.27E-02\t2.13E-03\t3.19E-04\t3.81E-03\t1.29E-02\t2.28E-02\t2.22E-16\t9.31E-03\t3.63E-04\t6.68E-04\t1.55E-04\t4.00E-02\t9.98E-05\t4.31E-03\t1.32E-03\t7.23E-03\t2.35E-03\t1.21E-03\t3.26E-03\t1.08E-03\t1.20E-02\t2.76E-02\t1.22E-03\t3.02E-04\t2.03E-03\t4.89E-03\t3.03E-02\t1.15E-03\t1.41E-03\t3.79E-04\t2.03E-04\t4.01E-03\t1.69E-03\t3.20E-01\t1.10E-03\t4.38E-04\t1.38E-03\t1.00E-02\t1.92E-03\t4.98E-03\t6.01E-03\t2.05E-02\t3.41E-03\t2.78E-05\t1.91E-02\t1.24E-03\t2.37E-02\t1.58E-03\t2.97E-04\t6.31E-05\t1.02E-03\t2.47E-03\t1.07E-03\t4.72E-03\t8.84E-02\t4.16E-03\t1.85E-03\t2.96E-03\t2.16E-03\t1.35E-03\t1.09E-04\t0.000022365153\t0.003823534542\r",
 "T>C\tTTT\t4.33E-03\t3.58E-05\t1.39E-02\t6.36E-04\t1.85E-02\t1.74E-03\t8.12E-04\t2.97E-03\t8.81E-02\t4.31E-02\t4.57E-03\t3.84E-02\t3.25E-03\t2.51E-05\t1.06E-04\t5.97E-02\t5.31E-04\t2.92E-03\t1.35E-03\t8.39E-03\t2.10E-02\t1.00E-03\t3.37E-03\t1.47E-03\t1.93E-02\t4.07E-02\t2.22E-16\t7.32E-18\t3.67E-03\t1.70E-02\t5.74E-02\t1.53E-03\t2.06E-03\t8.22E-04\t4.45E-04\t1.21E-02\t6.84E-04\t9.12E-03\t7.08E-03\t1.06E-03\t2.80E-03\t6.14E-02\t3.06E-03\t9.48E-03\t1.02E-02\t4.61E-02\t2.36E-03\t1.84E-03\t3.71E-02\t5.45E-04\t3.65E-02\t1.15E-02\t5.32E-04\t2.08E-04\t5.11E-03\t2.13E-02\t2.00E-03\t2.19E-03\t1.25E-02\t8.87E-04\t2.85E-04\t1.95E-01\t3.20E-02\t2.61E-03\t1.32E-03\t0.009820763708\t0.094054377313\r",
 "T>G\tATA\t1.72E-04\t2.38E-04\t3.95E-03\t1.05E-03\t3.86E-03\t2.57E-05\t4.73E-04\t8.44E-05\t1.19E-02\t3.70E-03\t1.64E-03\t4.51E-02\t1.92E-03\t3.47E-05\t2.18E-04\t2.01E-03\t6.14E-05\t2.12E-04\t4.68E-04\t1.12E-02\t8.51E-04\t1.41E-03\t7.21E-04\t2.33E-03\t1.05E-05\t3.92E-04\t2.33E-03\t8.28E-18\t3.21E-03\t5.28E-03\t2.92E-03\t1.14E-03\t9.87E-03\t4.07E-04\t3.90E-05\t3.68E-03\t3.86E-04\t9.79E-04\t2.82E-03\t3.05E-03\t1.09E-03\t2.17E-02\t1.15E-03\t3.88E-03\t7.23E-03\t7.53E-03\t1.20E-03\t2.63E-03\t2.40E-08\t1.17E-04\t1.00E-04\t1.71E-04\t5.17E-04\t1.96E-04\t7.43E-04\t7.23E-04\t1.73E-03\t4.59E-03\t3.87E-05\t2.62E-04\t1.21E-04\t7.06E-03\t1.35E-02\t1.79E-03\t2.50E-02\t0.000127137772\t0.023838221813\r",
 "T>G\tATC\t2.07E-04\t7.46E-05\t2.60E-03\t1.59E-04\t2.60E-03\t2.04E-04\t2.47E-04\t9.32E-04\t1.05E-03\t2.66E-03\t6.68E-04\t6.93E-03\t1.24E-03\t5.33E-05\t4.67E-05\t6.64E-04\t2.22E-16\t2.97E-03\t4.15E-04\t4.41E-03\t7.00E-04\t4.48E-03\t8.58E-04\t1.89E-03\t2.13E-05\t4.37E-04\t5.09E-05\t5.43E-18\t1.40E-03\t0.00E+00\t2.53E-03\t9.59E-04\t1.73E-02\t8.17E-04\t1.10E-04\t2.48E-03\t4.29E-04\t3.02E-04\t4.62E-04\t9.66E-04\t6.23E-04\t1.06E-02\t1.07E-03\t2.97E-03\t7.05E-03\t1.64E-03\t4.05E-04\t4.07E-03\t3.26E-05\t5.84E-04\t8.00E-04\t1.72E-03\t4.10E-05\t8.28E-05\t6.07E-04\t4.09E-03\t5.29E-04\t2.90E-03\t1.58E-03\t9.20E-04\t7.99E-18\t7.89E-04\t7.88E-04\t1.56E-03\t1.05E-02\t0.001270955861\t0.003434920439\r",
 "T>G\tATG\t2.68E-04\t2.05E-06\t6.30E-03\t1.42E-03\t7.94E-03\t1.07E-04\t2.76E-04\t1.47E-03\t6.04E-04\t1.06E-02\t4.77E-03\t5.33E-03\t3.20E-04\t5.23E-04\t1.79E-04\t8.08E-04\t8.25E-05\t7.09E-04\t1.85E-04\t4.62E-03\t7.19E-04\t9.90E-04\t1.99E-03\t1.77E-03\t2.69E-05\t1.87E-04\t2.95E-03\t5.68E-18\t2.06E-03\t7.14E-03\t2.29E-03\t6.77E-05\t7.68E-03\t4.20E-04\t4.10E-04\t4.02E-03\t1.27E-03\t4.75E-04\t1.97E-04\t4.48E-03\t9.36E-04\t1.68E-02\t2.07E-03\t7.66E-03\t8.06E-03\t2.70E-03\t6.37E-04\t3.49E-02\t7.47E-06\t9.84E-04\t1.60E-03\t3.68E-03\t3.51E-04\t1.26E-04\t3.94E-03\t2.78E-02\t1.19E-03\t5.06E-03\t1.15E-03\t1.45E-01\t1.27E-05\t1.23E-03\t2.71E-03\t1.70E-03\t2.61E-02\t0.000952698782\t0.004485500819\r",
 "T>G\tATT\t1.12E-04\t3.76E-06\t3.97E-03\t2.16E-04\t4.83E-03\t3.86E-04\t3.09E-04\t7.33E-04\t2.22E-16\t1.43E-02\t7.40E-04\t2.39E-02\t9.22E-03\t2.22E-16\t2.85E-04\t1.28E-03\t3.20E-04\t1.11E-02\t3.63E-03\t7.96E-03\t1.25E-02\t6.91E-02\t1.84E-03\t1.62E-03\t5.82E-05\t6.95E-05\t4.18E-05\t7.43E-18\t4.71E-03\t0.00E+00\t5.93E-03\t2.56E-04\t1.15E-01\t8.10E-04\t7.94E-04\t7.83E-04\t6.45E-18\t7.58E-04\t1.15E-03\t2.63E-03\t9.84E-04\t2.47E-02\t2.86E-03\t6.42E-03\t1.80E-02\t7.61E-03\t9.11E-04\t5.64E-03\t1.42E-03\t1.66E-04\t6.00E-04\t8.08E-03\t1.99E-04\t1.26E-04\t1.03E-04\t2.13E-03\t1.86E-03\t2.06E-03\t3.40E-04\t2.28E-03\t5.75E-03\t5.97E-03\t2.08E-03\t3.44E-03\t1.64E-02\t0.001764975843\t0.003178596391\r",
 "T>G\tCTA\t3.55E-05\t1.53E-04\t4.37E-03\t8.70E-04\t2.59E-03\t4.95E-05\t8.42E-05\t3.59E-04\t3.28E-03\t1.48E-03\t1.22E-03\t2.60E-02\t1.55E-04\t2.27E-04\t2.22E-16\t8.93E-04\t4.14E-05\t4.85E-05\t2.06E-04\t3.96E-03\t3.10E-04\t2.05E-03\t1.11E-04\t7.14E-04\t9.63E-04\t4.35E-06\t4.91E-03\t6.85E-18\t1.43E-03\t0.00E+00\t1.42E-03\t0.00E+00\t1.46E-02\t7.81E-04\t1.98E-04\t1.97E-03\t5.81E-05\t7.28E-04\t1.06E-03\t3.71E-03\t4.78E-04\t2.60E-03\t9.83E-04\t5.30E-03\t3.35E-03\t3.47E-03\t8.70E-04\t2.73E-03\t2.34E-04\t1.02E-04\t0.00E+00\t1.37E-04\t1.34E-04\t1.06E-04\t8.94E-04\t2.78E-03\t8.03E-04\t1.49E-03\t4.36E-03\t1.65E-03\t5.82E-03\t1.40E-04\t5.03E-04\t1.01E-03\t3.23E-03\t0.000897193499\t0.003736474707\r",
 "T>G\tCTC\t2.12E-04\t1.55E-04\t7.04E-03\t4.41E-04\t4.61E-03\t7.84E-04\t1.72E-04\t7.92E-04\t5.58E-03\t4.98E-05\t1.95E-03\t7.41E-03\t2.23E-16\t1.15E-03\t2.22E-16\t1.14E-03\t2.82E-06\t1.21E-03\t1.57E-04\t1.08E-03\t1.04E-03\t3.14E-02\t1.99E-03\t2.09E-03\t1.60E-03\t5.44E-05\t8.80E-04\t1.48E-04\t8.80E-04\t2.89E-03\t4.94E-03\t3.81E-04\t2.21E-02\t5.69E-06\t2.15E-04\t1.43E-03\t4.50E-04\t9.49E-04\t2.22E-16\t5.18E-04\t5.66E-04\t2.38E-03\t1.41E-03\t7.39E-03\t5.67E-03\t2.41E-03\t7.53E-04\t7.26E-03\t2.53E-03\t2.81E-04\t1.80E-03\t5.15E-05\t9.28E-05\t2.67E-05\t6.18E-05\t2.56E-03\t7.75E-04\t4.24E-03\t1.89E-02\t4.66E-04\t3.93E-07\t8.47E-04\t1.04E-03\t2.62E-03\t4.16E-03\t0.002148873221\t0.002478170057\r",
 "T>G\tCTG\t1.28E-04\t1.21E-04\t1.07E-02\t3.79E-03\t6.13E-03\t9.24E-04\t2.23E-16\t1.29E-03\t1.19E-03\t2.58E-05\t4.70E-03\t1.42E-02\t5.35E-05\t4.79E-04\t2.22E-16\t5.93E-04\t2.22E-16\t4.35E-04\t2.72E-04\t2.24E-03\t1.94E-06\t1.37E-02\t2.93E-03\t1.33E-03\t2.71E-03\t2.01E-05\t8.03E-03\t6.52E-05\t3.04E-03\t1.18E-02\t5.56E-03\t4.03E-04\t1.55E-02\t3.32E-04\t7.85E-05\t1.50E-03\t1.09E-03\t1.07E-03\t4.36E-04\t6.36E-03\t8.87E-04\t3.77E-03\t1.59E-03\t9.83E-03\t5.59E-03\t3.04E-03\t5.68E-03\t4.65E-02\t3.47E-03\t6.49E-04\t1.60E-03\t3.67E-03\t1.04E-04\t1.12E-04\t3.66E-03\t5.93E-03\t7.96E-04\t2.98E-03\t2.21E-02\t1.13E-01\t4.45E-03\t2.19E-03\t6.99E-04\t2.02E-03\t7.11E-03\t0.000788671403\t0.001120299176\r",
 "T>G\tCTT\t1.71E-04\t2.03E-04\t7.03E-03\t1.20E-03\t7.30E-03\t1.30E-03\t2.65E-05\t3.68E-06\t3.40E-03\t1.56E-04\t5.74E-04\t5.65E-02\t1.30E-03\t5.23E-03\t1.33E-04\t3.44E-04\t5.55E-05\t9.62E-04\t2.49E-03\t2.54E-03\t5.31E-03\t5.50E-01\t4.14E-03\t2.20E-03\t5.11E-03\t1.85E-04\t1.19E-03\t2.47E-04\t4.71E-03\t8.55E-03\t7.80E-03\t1.44E-03\t1.01E-01\t7.67E-04\t2.65E-04\t9.49E-04\t6.76E-04\t1.92E-04\t6.55E-07\t2.71E-03\t1.07E-03\t4.11E-03\t6.88E-03\t4.94E-03\t1.31E-02\t1.05E-02\t3.48E-03\t6.26E-18\t3.84E-04\t7.28E-04\t2.00E-03\t1.07E-04\t6.57E-05\t9.29E-05\t4.71E-04\t2.30E-03\t8.63E-04\t2.72E-03\t8.86E-04\t1.18E-03\t6.78E-18\t5.16E-03\t4.73E-04\t2.12E-02\t1.80E-05\t0.000000000000\t0.001756870383\r",
 "T>G\tGTA\t2.23E-16\t2.09E-04\t4.42E-03\t5.98E-04\t2.52E-03\t4.35E-05\t9.06E-05\t3.23E-04\t2.48E-03\t7.48E-05\t1.31E-03\t8.30E-03\t2.04E-04\t4.81E-04\t9.57E-05\t4.12E-04\t8.77E-05\t7.79E-05\t2.03E-04\t9.02E-04\t4.95E-04\t2.81E-04\t8.70E-04\t6.36E-04\t2.22E-16\t2.86E-04\t1.59E-03\t6.26E-18\t1.07E-03\t5.55E-03\t1.37E-03\t0.00E+00\t1.16E-03\t4.49E-05\t1.88E-04\t3.04E-04\t1.15E-03\t5.81E-04\t3.65E-04\t1.76E-03\t1.85E-04\t1.18E-02\t5.04E-04\t3.33E-03\t2.49E-03\t1.80E-03\t1.14E-03\t6.11E-02\t7.68E-18\t4.40E-04\t0.00E+00\t2.20E-02\t7.64E-06\t6.61E-05\t8.67E-03\t2.19E-03\t1.13E-03\t4.38E-03\t1.70E-04\t5.67E-03\t6.72E-18\t1.40E-03\t8.07E-03\t1.23E-03\t2.82E-02\t0.000018879250\t0.007699623603\r",
 "T>G\tGTC\t2.23E-16\t1.33E-04\t2.34E-03\t2.53E-04\t1.71E-03\t2.13E-04\t1.21E-04\t3.21E-04\t1.66E-03\t2.23E-05\t6.56E-04\t5.85E-03\t8.33E-04\t9.42E-04\t2.22E-16\t2.13E-04\t2.22E-16\t2.10E-03\t5.80E-04\t9.01E-04\t2.04E-04\t6.25E-03\t3.70E-05\t7.03E-04\t6.65E-05\t5.33E-04\t2.36E-04\t8.69E-05\t1.53E-03\t0.00E+00\t1.21E-03\t5.08E-05\t1.93E-03\t9.17E-05\t3.53E-04\t8.72E-04\t7.40E-04\t6.37E-04\t2.22E-16\t5.07E-03\t3.13E-04\t6.86E-03\t9.39E-04\t5.80E-03\t2.24E-03\t9.50E-04\t7.05E-05\t7.56E-02\t1.83E-05\t5.20E-04\t1.00E-04\t1.04E-02\t1.35E-04\t4.37E-05\t9.59E-03\t5.30E-03\t5.61E-04\t1.74E-03\t1.53E-03\t3.20E-03\t9.19E-07\t1.71E-04\t2.69E-03\t1.65E-03\t2.88E-02\t0.001229435537\t0.000800067973\r",
 "T>G\tGTG\t3.48E-04\t3.71E-05\t1.09E-02\t2.43E-03\t5.30E-03\t2.74E-04\t3.25E-04\t2.41E-03\t6.03E-03\t3.80E-03\t3.97E-03\t3.97E-03\t2.23E-16\t1.73E-03\t1.30E-04\t4.96E-04\t4.26E-04\t7.81E-04\t3.47E-04\t6.75E-04\t8.35E-04\t1.48E-03\t5.83E-03\t5.10E-04\t8.70E-05\t2.13E-04\t4.50E-04\t1.57E-04\t3.13E-03\t8.69E-03\t2.35E-03\t4.74E-03\t3.13E-03\t7.70E-04\t1.38E-03\t2.55E-03\t3.01E-03\t5.42E-03\t5.65E-05\t4.04E-03\t1.53E-03\t9.83E-03\t2.16E-03\t7.72E-03\t4.45E-03\t1.20E-03\t4.10E-04\t2.41E-01\t1.03E-04\t1.79E-04\t1.20E-03\t2.72E-02\t6.63E-05\t1.47E-04\t3.51E-02\t7.62E-02\t1.37E-03\t4.73E-03\t7.45E-18\t2.99E-01\t6.72E-18\t8.61E-04\t2.34E-02\t2.29E-03\t6.62E-01\t0.000212751233\t0.001946513805\r",
 "T>G\tGTT\t1.46E-05\t2.23E-16\t5.83E-03\t2.52E-04\t2.35E-03\t7.87E-04\t8.35E-04\t1.83E-03\t6.75E-03\t6.84E-03\t1.29E-03\t1.94E-02\t2.21E-03\t5.16E-03\t1.54E-04\t1.65E-04\t2.89E-06\t2.25E-03\t1.39E-03\t2.39E-03\t7.24E-03\t1.14E-01\t2.14E-03\t3.41E-04\t3.56E-03\t2.59E-03\t1.96E-05\t5.78E-18\t2.68E-03\t2.77E-03\t6.05E-03\t2.30E-03\t2.06E-02\t1.15E-03\t1.60E-04\t2.81E-03\t7.85E-04\t4.51E-04\t4.17E-04\t3.35E-03\t4.89E-04\t1.71E-02\t3.34E-03\t5.94E-03\t5.48E-03\t3.47E-03\t9.20E-04\t7.46E-02\t2.80E-03\t4.96E-04\t6.00E-04\t8.84E-02\t1.15E-05\t1.02E-04\t3.65E-02\t2.09E-03\t8.64E-04\t2.03E-03\t1.21E-05\t1.44E-03\t1.14E-04\t2.83E-03\t1.69E-02\t7.45E-03\t4.74E-02\t0.001340078646\t0.002927905616\r",
 "T>G\tTTA\t2.23E-16\t1.67E-05\t7.25E-03\t3.77E-04\t5.22E-03\t1.05E-04\t1.28E-04\t9.55E-04\t1.93E-02\t2.11E-04\t2.22E-16\t6.57E-02\t2.69E-03\t1.13E-02\t8.51E-05\t1.63E-03\t1.84E-04\t9.92E-05\t2.31E-04\t9.78E-03\t2.22E-16\t1.22E-05\t6.96E-04\t1.66E-03\t2.22E-16\t4.65E-04\t5.06E-03\t1.86E-04\t1.23E-04\t2.08E-03\t1.43E-03\t1.19E-03\t4.23E-02\t4.67E-04\t1.18E-04\t2.25E-03\t2.26E-05\t4.68E-04\t2.22E-16\t4.29E-03\t7.88E-04\t1.83E-02\t1.17E-03\t8.39E-03\t8.46E-03\t4.19E-02\t1.15E-03\t2.04E-04\t7.68E-18\t2.85E-04\t7.00E-04\t2.41E-03\t2.50E-04\t1.50E-04\t3.27E-04\t2.55E-03\t1.43E-03\t1.30E-03\t3.24E-05\t6.06E-03\t9.46E-03\t3.22E-03\t7.98E-04\t1.60E-03\t1.61E-03\t0.001875698241\t0.007106411217\r",
 "T>G\tTTC\t5.51E-05\t7.04E-05\t6.28E-03\t1.74E-04\t6.56E-03\t2.87E-04\t1.16E-04\t1.55E-03\t1.74E-02\t1.15E-04\t1.14E-03\t8.62E-03\t2.23E-16\t5.55E-03\t6.73E-05\t1.13E-03\t2.22E-16\t1.20E-03\t2.94E-04\t5.51E-03\t1.16E-04\t8.75E-03\t2.10E-03\t2.02E-03\t2.22E-16\t3.79E-04\t7.52E-04\t5.27E-18\t2.64E-03\t5.79E-04\t1.72E-03\t2.80E-04\t3.57E-02\t9.93E-05\t9.64E-05\t2.38E-03\t2.24E-04\t1.64E-04\t2.22E-16\t1.94E-03\t7.44E-04\t1.95E-02\t9.03E-04\t5.08E-03\t1.14E-02\t1.55E-02\t2.74E-04\t8.25E-03\t7.68E-18\t1.35E-04\t1.00E-03\t2.25E-05\t2.22E-04\t1.38E-04\t6.55E-04\t3.47E-03\t1.12E-03\t4.36E-03\t7.06E-04\t5.65E-03\t1.00E-06\t3.96E-03\t1.59E-03\t2.33E-03\t1.06E-02\t0.000737916543\t0.006404453857\r",
 "T>G\tTTG\t5.83E-04\t9.54E-05\t8.05E-03\t2.32E-03\t6.94E-03\t3.24E-04\t2.23E-16\t1.35E-03\t7.64E-03\t1.25E-04\t3.09E-03\t1.09E-02\t2.16E-05\t2.76E-03\t1.01E-04\t8.75E-04\t3.67E-05\t1.14E-04\t1.48E-04\t5.24E-03\t9.20E-04\t4.78E-03\t1.45E-03\t2.07E-03\t2.79E-03\t1.31E-04\t5.57E-03\t1.82E-04\t3.04E-03\t9.43E-03\t2.84E-03\t2.35E-03\t1.40E-02\t1.80E-04\t8.13E-04\t5.35E-03\t8.18E-04\t9.59E-04\t1.52E-03\t4.82E-03\t9.02E-04\t3.02E-02\t1.90E-03\t7.28E-03\t8.66E-03\t1.94E-02\t2.93E-03\t2.28E-02\t9.17E-04\t7.08E-04\t9.00E-04\t1.54E-03\t3.50E-04\t1.12E-04\t8.72E-03\t2.76E-03\t1.48E-03\t6.15E-03\t7.69E-03\t7.95E-02\t1.01E-04\t5.40E-03\t6.74E-03\t1.96E-03\t2.39E-02\t0.000000000000\t0.003409999931\r",
 "T>G\tTTT\t2.23E-16\t2.23E-16\t1.05E-02\t5.68E-04\t1.35E-02\t1.01E-03\t8.29E-05\t1.77E-03\t2.17E-02\t1.66E-04\t9.96E-04\t6.39E-02\t1.89E-02\t9.09E-02\t5.55E-05\t2.21E-03\t1.89E-05\t6.05E-04\t5.99E-03\t1.58E-02\t4.58E-03\t1.22E-01\t5.16E-03\t4.23E-03\t1.51E-03\t2.22E-16\t5.95E-04\t7.32E-18\t3.96E-03\t7.87E-03\t9.45E-03\t1.40E-04\t4.79E-01\t9.28E-04\t8.89E-03\t7.30E-03\t6.35E-04\t2.25E-03\t2.32E-04\t3.62E-03\t1.81E-03\t2.91E-02\t5.29E-03\t1.09E-02\t2.50E-02\t8.79E-02\t1.74E-03\t1.32E-03\t9.00E-03\t7.24E-04\t1.00E-03\t3.62E-03\t3.20E-07\t3.64E-04\t1.60E-03\t2.91E-03\t2.53E-03\t7.05E-03\t4.99E-04\t4.12E-02\t2.90E-02\t1.19E-01\t2.17E-02\t1.30E-17\t3.18E-02\t0.006603134330\t0.012708761755\r")
#}}}

# MAIN
run()
}
