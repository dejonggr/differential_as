#!/usr/bin/env Rscript

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
  suppressWarnings(require(magrittr, quietly = T, warn.conflicts = FALSE))
  suppressWarnings(require(devtools, quietly = T, warn.conflicts = FALSE))
  suppressWarnings(require(plyr, quietly = T, warn.conflicts = FALSE))
  suppressWarnings(require(dplyr, quietly = T, warn.conflicts = FALSE))
  suppressWarnings(require(bindrcpp, quietly = T, warn.conflicts = FALSE))
  suppressWarnings(require(tidyr, quietly = T, warn.conflicts = FALSE))
  suppressWarnings(require(reshape2, quietly = T, warn.conflicts = FALSE))
  suppressWarnings(require(stats, quietly = T, warn.conflicts = FALSE))
  suppressWarnings(require(car, quietly = T, warn.conflicts = FALSE))
  suppressWarnings(require(betareg, quietly = T, warn.conflicts = FALSE))
  suppressWarnings(require(brglm2, quietly = T, warn.conflicts = FALSE))
}


## Help section
if(length(args) < 1) {
  args <- c("--help")
}
if("--help" %in% args) {
  cat("
      Differential Alternative Splicing Analysis
      
      Arguments:
      --arg1=junction_counts.tsv  	- Table of constitutive vs alternative read counts per 
        junction. Must be tsv.
      OPTIONAL: 
       sample = sampleinfo 		- string of condition information: uses this order your factors before DAS analysis.
       plot =/plots/output/file/path  	- This will output a sample of 100 plots so choose
          the directory wisely!
       plot                          	- Will output the plot by default to your other output folder
        (WIP) SUPPA                   	- If you want to run this with SUPPA, include this string
       out=/output/file/path/        	- out=(insert filepath here)
        \n\n\n")
  
  q(save="no")
}
if(length(args) > 3) {
  cat("Too many args. You need (1) junction read counts (2) sample info.")
}

## Loading the data
message("Loading the data")
arsenal <- read.table(args[1], sep = '\t', header = T); arsenal$species<-NULL

# Data.frame of 
if(!any(grepl("sample", args))){
  # take column names as default
  message("No sampleinfo vector provided; using 'M' and 'I' as default conditions.")
  sampleinfo = colnames(arsenal)[grep("_alt", colnames(arsenal))] %>% gsub("[0-9].*$", "", .) %>% unique
  sampleinfo = factor(sampleinfo, levels = c(sampleinfo[grep("M", sampleinfo)], sampleinfo[grep("I", sampleinfo)]))
} else if(any(grepl("sample", args))){
sampleinfo = read.table(args[2], sep = '\t', header = T)
}


#################
### Functions ###
#################

# only accept moderlately support events; default is 8 read in total
filter_events = function(x,n=8){
  cols_match = c("alt", levels(factor(gsub("(.*?)[0-9]_.*", "\\1", colnames(x)[grepl("_alt", colnames(x)) | grepl("_nrm", colnames(x))]))))
  
  x %>% filter_at(vars(contains("alt")), all_vars(.>0)) %>%
    filter_at(vars(contains("nrm")), all_vars(.>0)) %>%
    mutate(alt_sum = rowSums(.[,grepl(cols_match[1], colnames(.))])) %>%
    filter(alt_sum>n) %>%
    select(-alt_sum)
}

# format for regression fit
BLR_format = function(x){
  alt = x %>% select(-contains("nrm")) %>% gather(key = "read_type", value = "alt_counts", contains("alt"), factor_key = T) %>%
    mutate(sample_data=gsub("_alt", "",.$read_type)) %>% select(-read_type) 
  nrm = x %>% select(-contains("alt")) %>% gather(key = "read_type", value = "nrm_counts", contains("nrm"), factor_key = T) %>%
    mutate(sample_data=gsub("_nrm", "\\1",.$read_type)) %>% select(-read_type)
  df = cbind(alt,nrm)
  if(all(df[,1]==df[,4]) & all(df[,3]==df[,6])){
  out = df[,unique(colnames(df))] %>%
    mutate(condition_name = gsub("[0-9]+","",.$sample_data)) %>%
    mutate(condition_name = factor(condition_name, levels = levels(sampleinfo))) %>%
    mutate(condition = ifelse(condition_name == levels(factor(condition_name))[1], 1, 0)) %>%
    mutate(condition = as.factor(condition))
  } else {
    message("There was a problem with the event merge!")
  }
  invisible(out)
} 
BLR_fit = function(e){
 # benchmark placed it at 46% the run time of lapply
  fit=with(e, 
           by(e, e[,"event_id"],
              function(x){
                fit=glm(cbind(nrm_counts,alt_counts) ~ condition, family= binomial(link="logit"), data = x, method = "brglmFit")}))
}
BLR_test = function(x){
  test=lapply(x, function(x){
    p=Anova(x, type = 3)
    p$`Pr(>Chisq)`[1]})
  p.adjust(test, method = "BH")
}

# Calculate FC and delta PSI for AS
LFC_standard = function(x){
  
  # Calculate PSI
  x = mutate(x, PSI = alt_counts/(nrm_counts+alt_counts))
  
  df = with(x, 
           by(x, x[,"event_id"], 
              function(x){
                by(x, x["condition"], function(x) as.data.frame(mean(x$PSI))) # Calculate means
              }))
  fold_change(df, n)
}
fold_change = function(x,n){
  lapply(x, function(x){
    name_fc = paste0(n, "_log2fc")
    name_dpsi = paste0(n, "_dPSI")
   object = data.frame(log2(x[[1]]/x[[2]]), # fold change
               x[[1]]-x[[2]])  # dpsi
  colnames(object) = c(name_fc, name_dpsi)
  invisible(object)
  }) %>% do.call(rbind,.)
}

#############################
### Differential splicing ###
#############################

Juncs = filter_events(arsenal) %>% 
  # clean up df:
  unite(event_id, c("gene_name", "junct", "class", "unique_code")) %>%
  BLR_format(.)

message("Fitting models...")
rPSI = BLR_fit(Juncs)
message("Calculating LFCs and delta PSI")

message("Calculating LFCs...")
n = paste0(levels(factor(Juncs[,"condition_name"])), collapse = "|")
LFC = LFC_standard(Juncs)

message("Estimating significance and mutliple test corrections...")
out = data.frame(BLR_test(rPSI))
colnames(out) = paste0(n, "_padj")

message("Saving results. If you provided none, check your R working directory.")
dpsi_output = merge(out, LFC, by = "row.names") %>% rename(Event = Row.names)


if(any(grepl("out=", args))){
  output=args[grep("out=", args)] %>% as.character(.) %>% gsub("out=(.*?)", "\\1", .)
  write.table(dpsi_output, file = paste0(output,"differential_splicing.txt"), row.names = T, col.names = T, quote = F, sep = "\t")
}else{write.table(dpsi_output, row.names = T, col.names = T, sep = "\t", quote = F)}

message("Done!")

#################
### Plot test ###
#################
if(any(grepl("plot", args))){

  sample_plots = sample(unique(dpsi_output[dpsi_output[,grep("_padj", colnames(dpsi_output))] < 0.01,]$Event),10)
  
  sample_juncs = Juncs[Juncs$event_id %in% sample_plots,] %>%
    mutate(PSI = alt_counts/(alt_counts+nrm_counts))
  
## Plotting functions ##
  #Plot
  point_plot=function(x){
    ggplot(x, aes(x=factor(condition_name)))+
      geom_point(aes(y=PSI))+
      ylim(0,1)+
      xlab("Conditions")+
      ylab("PSI")+
      guides(fill=FALSE)+
      theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank())+
      theme(axis.text.x=element_text(angle=90,hjust=1))}
  
  # Generate the plots
  p_plots = lapply(sample_plots, function(x){
    p = sample_juncs[sample_juncs$event_id==x,]
    point_plot(p)
  })

# Save the plots
  if(any(grepl("plot=", args))){
    plot_output=args[grep("plot=", args)] %>% as.character(.) %>% gsub("plot=(.*?)", "\\1", .)
    lapply(names(p_plots), function(x){ggsave(filename = paste0(plot_output,"/",x,".png"),  
                                          plot=v_plots[[x]],
                                          width = 3, height = 3, 
                                          units = "in",
                                          dpi = 800)})
  } else {
    lapply(names(p_plots), function(x){ggsave(filename = paste0(output,x,".png"),  
                                              plot=v_plots[[x]],
                                              width = 3, height = 3, 
                                              units = "in",
                                              dpi = 800)})
}
}

###############
### Testing ###
###############

#arsenal <- read.csv("/Users/grantdejong/Desktop/BNapusSenseAS_crit5_minf0.05_maxf0.csv"); arsenal$species<-NULL
#arsenal <- read.table("", sep = '\t', header = TRUE); arsenal$species<-NULL
#sampleinfo<- read.table("/Users/grantdejong/Desktop/sampleinfo-jon.tsv", sep = '\t', header = T)
#dea <- read.csv("/Users/grantdejong/Desktop/John_DAS_list_Bna_test2.csv", header = T)

