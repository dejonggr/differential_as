#!/usr/bin/env Rscript

arsenal <- read.table("~/Downloads/John_AS_Bna_test2.tsv", sep = '\t', header = TRUE); arsenal$species<-NULL
sampleinfo<- read.table("/Users/grantdejong/Desktop/sampleinfojon2.txt", sep = '\t', header = T)
dea <- read.csv("/Users/grantdejong/Desktop/John_DAS_list_Bna_test2.csv", header = T)


# Add args to include: c("PSI", "PSI-SUPPA", "read counts") inputs
# Edit --help docs

#####################################################################################
### Calculates differential splicing based on read counts of alt/cons transcripts ###
#####################################################################################
###################################
### Required packages and setup ###
###################################
rm( list= ls())
args <- commandArgs(trailingOnly = TRUE)

## Libraries
if(length(args) > 1) {
  # Some of these are redundant
  require(magrittr, quietly = T, warn.conflicts = FALSE)
  require(devtools, quietly = T, warn.conflicts = FALSE)
  require(plyr, quietly = T, warn.conflicts = FALSE)
  require(dplyr, quietly = T, warn.conflicts = FALSE)
  require(tidyr, quietly = T, warn.conflicts = FALSE)
  require(reshape2, quietly = T, warn.conflicts = FALSE)
  require(stats, quietly = T, warn.conflicts = FALSE)
  require(car, quietly = T, warn.conflicts = FALSE)
  require(betareg, quietly = T, warn.conflicts = FALSE)
}

## Help section
if(length(args) < 1) {
  args <- c("--help")
}
if("--help" %in% args) {
  cat("
      Differential Alternative Splicing Analysis
      
      Arguments:
      --arg1=junction_counts.csv  - Table of constitutive vs alternative read counts per 
      junction. Must be csv.
      --arg2=sampleinfo.tsv       - Sample info. Must be tsv. Two columns: 
      (1) sample_data - a unique rep ID (2) conditions - 
      condition for each replicate.
      OPTIONAL: 
      SUPPA                       - If you want to run this with SUPPA, include this string
      out=/output/filepath/       - out=(insert filepath here)
      \n\n\n")
  
  q(save="no")
}
if(length(args) > 3) {
  cat("Too many args. You need (1) junction read counts (2) sample info.")
}

## Loading the data
print("Loading the data")
arsenal <- read.csv(args[1], row.names = 1); arsenal$species<-NULL
sampleinfo <- read.table(args[2], sep = '\t', header = T)

## Setting up reference level for trt vs ctrl experimental conditions
sampleinfo <- if(any(grepl("ref=", args))){
  reference=args[grep("ref=", args)] %>% as.character(.) %>% gsub("ref=(.*?)", "\\1", .)
  sampleinfo %>% mutate(condition = relevel(.$condition, reference))
} else if("control" %in% sampleinfo$condition) {
  sampleinfo %>% mutate(condition = relevel(.$condition, "control"))
} else {
  return(sampleinfo)
}

#################
### Functions ###
#################
## Format tack output
format_tack<-function(x){
  x = unite(x, Event, c("gene_name", "junct", "class", "unique_code")) %>%
    filter_at(vars(contains("alt")), all_vars(.>0)) %>% # Maybe filter with rowSums?
    filter_at(vars(contains("nrm")), all_vars(.>0)) %>%
    na.omit(.)
}
# Remove events where alt>nrm

## Differential Splicing: Logisitic regression
BLR_format<-function(x){
  alt = x %>% select(-contains("nrm")) %>% gather(key = "read_type", value = "alt_counts", contains("alt"), factor_key = T) %>%
    mutate(sample_info=gsub("(.*?)_a.*$", "\\1",.$read_type)) %>% select(-read_type)
  nrm = x %>% select(-contains("alt")) %>% gather(key = "read_type", value = "nrm_counts", contains("nrm"), factor_key = T) %>%
    mutate(sample_info=gsub("(.*?)_n.*$", "\\1",.$read_type)) %>% select(-read_type)
  p = merge(alt,nrm,by=c("Event", "sample_info")) %>% 
    merge(., sampleinfo, by="sample_info")
  dv1 = cbind(p,sapply(levels(p$condition), function(x) as.integer(x == p$condition)))
  if("control" %in% sampleinfo$condition){ # i.e. control = default ref
    x = cbind(p,sapply(levels(p$condition), function(x){dv1[x]=ifelse(dv1$control==1, -1, dv1[,x])}))
    x = x %>%  mutate(condition = relevel(.$condition, "control"))
    invisible(return(x))
  }
  else{invisible(return(dv1))}
}
BLR_fit<-function(x){
  form = noquote(args[grep("formula=", args)] %>% gsub("formula=(.*?)", "\\1", .))
  fit=with(x, 
           by(x, x[,"Event"],
              function(x){
                fit=glm(noquote(paste0("cbind(nrm_counts,alt_counts) ~ ", form)), family= binomial(link="logit"), data = x)}))
}
BLR_test<-function(x){
  test=lapply(x, function(x){
    p=Anova(x, type = 2, error.estimate = "dispersion")
    t(p)[3,] }) %>% do.call(rbind,.)
  names = rownames(test)
  test = mapply(function(x){p.adjust(x, method = "BH")}, as.data.frame(test), SIMPLIFY = T, USE.NAMES = T) %>% 
    `rownames<-`(names)
  if(any(grepl("ref", args)) || any("control" %in% levels(sampleinfo$condition))){colnames(test) = levels(sampleinfo$condition)[-1]}
  invisible(return(test))
}

## PSI LFC estimation
# It may be good to add dispersion shrinkage a la apeglm but ensure it's useful for these data
# Include reference group?
fold_change<-function(x){
  ref=levels(sampleinfo$condition)[1]
  for(i in levels(sampleinfo$condition)){
    a=as.numeric(x[i])
    r=as.numeric(x[ref])
    if(!i==ref){
      x=log2(a/r)
      names(x)=paste0(i,":",ref,"_log2fc")
      return(x)
      }
  }
}
LFC_standard<-function(x){
  x = x %>% mutate(PSI=alt_counts/(nrm_counts+alt_counts))
  x = with(x, 
           by(x, x[,"Event"], 
              function(x){
                by(x, x["condition"], function(x) as.data.frame(mean(x$PSI))) # Calculate means
              }))
  lapply(x, function(x){fold_change(x)})
}
LFC_shrink<-function(x){
  
} # WIP


#############################
### Differential splicing ###
#############################

## FILTERING: I added code to filter PSI outside of (0,1) range - these values are biologically impossible or unrealisitc
# If you would prefer to include them, do not round them to 0, otherwise, add a fake near-infinite value.
# Alternatively, you could add a "qualitative differentiation" metric that shows conditional splicing asymmetry
Juncs = format_tack(arsenal)

## Diff. splicing analysis
message("Fitting model and estimating significance...")
Juncs = BLR_format(Juncs)
d_PSI = BLR_fit(Juncs)

message("Calculating fold-changes...")
LFC = LFC_standard(Juncs) %>% do.call(rbind,.)

message("Estimating significance and mutliple test corrections...")
out=BLR_test(d_PSI)

dpsi_output = merge(out, LFC, by = "row.names") %>% rename(Event = Row.names)

message("Saving output. If you provided none, check your R working directory.")
if(any(grepl("out=", args))){
  output=args[grep("out=", args)] %>% as.character(.) %>% gsub("out=(.*?)", "\\1", .)
  write.table(dpsi_output, file = paste0(output,"differential_splicing.txt"), row.names = T, col.names = T, quote = F, sep = "\t")
}else{write.table(dpsi_output, row.names = T, col.names = T, sep = "\t", quote = F)}

message("Done!")

#################
### Plot test ###
#################
if(any(grepl("plot", args))){

## Plotting functions ##
  #Plot
  violin_plot<-function(x){
    ggplot(x, aes(x=factor(condition)))+
      geom_point(aes(y=PSI))+
      ylim(0,1)+
      xlab("Conditions")+
      ylab("PSI")+
      guides(fill=FALSE)+
      theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank())+
      theme(axis.text.x=element_text(angle=90,hjust=1))}
  # Create plot-friendly dataframe
  filter_df=function(x,y){
    x = 
      sig_names = rownames(x)
    y = y %>% mutate(PSI=alt_counts/(nrm_counts+alt_counts))
    y = y[which(y$Event %in% sig_names),]
    with(y, 
         by(y, y[,"Event"], 
            function(x) x))
  }
  
dpsi_output %>% filter_at(vars(-contains("log")), any_vars(.<0.05))

# Make data set
  v = filter_df(dpsi_output)

# Generate the plots
  v_plots=lapply(head(v,100), violin_plot) # for each plot, only use sample
  
# Save the plots
  if(any(grepl("plot", args))){
    plot_output=args[grep("plot=", args)] %>% as.character(.) %>% gsub("out=(.*?)", "\\1", .)
    lapply(names(v_plots), function(x){ggsave(filename = paste0(plot_output,"/",x,".png"),  
                                          plot=v_plots[[x]],
                                          width = 3, height = 3, 
                                          units = "in",
                                          dpi = 800)})
  } else {
    lapply(names(v_plots), function(x){ggsave(filename = paste0(output,x,".png"),  
                                              plot=v_plots[[x]],
                                              width = 3, height = 3, 
                                              units = "in",
                                              dpi = 800)})
}
}
