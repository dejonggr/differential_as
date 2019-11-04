#!/usr/bin/env Rscript
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

#################
### Functions ###
#################
## Filtering PSI
filter_psi<-function(x){
  colnames(x) = paste0(colnames(x), "_PSI")
  x$Event = as.character(rownames(x))
  if(!("SUPPA" %in% args)){
    x %>% filter_at(vars(contains("PSI")), all_vars(.>1)) %>% 
      filter_at(vars(contains("PSI")), all_vars(.<100)) %>% na.omit(.)
  }else{
    x %>% filter_at(vars(contains("PSI")), all_vars(.>.01)) %>% 
      filter_at(vars(contains("PSI")), all_vars(.<1)) %>% na.omit(.)
  }
}

## Differential Splicing: Beta-regression

BR_format<-function(x){
  p = x %>% gather(key = "sample_info", value = "PSI", contains("PSI"), factor_key = T) %>%
    mutate(sample_info = gsub("(.*?)_PSI", "\\1", .$sample_info)) %>%
    merge(., sampleinfo, by = "sample_info")
  dv1 = cbind(p,sapply(levels(p$condition), function(x) as.integer(x == p$condition)))
  x = cbind(p,sapply(levels(p$condition), function(x){dv1[x]=ifelse(dv1$control==1, -1, dv1[,x])}))
  if("control" %in% sampleinfo$condition){
  x = x %>%  mutate(condition = relevel(.$condition, "control"))
  }
  if(!("SUPPA" %in% args)){
      x= x %>% mutate(PSI=.$PSI/100)
      invisible(return(x))
  }else(invisible(return(x)))
}

BR_fit<-function(x){
  fit=with(x, 
           by(x, x[,"Event"], 
              function(x){fit=betareg(PSI ~ cold + heat + drought, data = x)}))
}

BR_test<-function(x){
  test=lapply(x, function(x){
    Anova(x, type = 2, error.estimate = "dispersion")
  })
}

# fit=betareg(PSI ~ control:cold + control:heat + control:drought

#BR_test<-function(x){
#  test=lapply(x, function(x){
#    summary(x)$coef$mean[,4]})
#  test = test %>% do.call(rbind,.) %>% 
#    as.data.frame(.) %>% 
#    select(-`(Intercept)`)
#  names=row.names(test)
#  mapply(function(x){p.adjust(x, method = "BH")}, test, SIMPLIFY = T, USE.NAMES = T) %>%
#    `rownames<-`(names)
#}

#############################
### Differential splicing ###
#############################

## FILTERING: I added code to filter PSI outside of (0,1) range - these values are biologically impossible or unrealisitc
# If you would prefer to include them, do not round them to 0, otherwise, add a fake near-infinite value.
# Alternatively, you could add a "qualitative differentiation" metric that shows conditional splicing asymmetry
arsenal = filter_psi(arsenal)

## Diff. splicing analysis
print("Fitting model and estimating significance...")
PSI = BR_format(arsenal)
d_PSI = BR_fit(PSI)

print("Estimating significance and mutliple test corrections...")
out=BR_test(d_PSI)

print("Saving output. If not output given, check your R working directory.")
if(any(grepl("out=", args))){
  output=args[grep("out=", args)] %>% as.character(.) %>% gsub("out=(.*?)", "\\1", .)
  write.table(out, file = paste0(output,"differential_splicing.txt"), row.names = T, col.names = T, quote = F, sep = "\t")
}else{write.table(out, row.names = T, col.names = T, sep = "\t", quote = F)}

print("Done!")

#################
### Plot test ###
#################

# No FDR
z=lapply(out, function(x) t(x)[3,]) %>% do.call(rbind,.)

# FDR
names=rownames(z)
z2=mapply(function(x){p.adjust(x, method = "BH")}, as.data.frame(z), SIMPLIFY = T, USE.NAMES = T) %>% `rownames<-`(names)

filter_df=function(x,y){
  x = as.data.frame(x)
  x = x[which(x$cold<0.05 | x$heat<0.05 | x$drought<0.05),]
  sig_names = rownames(x)
  y = PSI[which(PSI$Event %in% sig_names),]
  with(y, 
       by(y, y[,"Event"], 
          function(x) x))
}

test=filter_df(z, filter_df)

violin_plot<-function(x){
  ggplot(x, aes(x=factor(condition)))+
    geom_violin(aes(y=PSI))+
    ylim(0,1)+
    xlab("Conditions")+
    ylab("PSI")+
    guides(fill=FALSE)+
    theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())+
    theme(axis.text.x=element_text(angle=90,hjust=1))}

# make the plots
v_plots=lapply(test, violin_plot) # for each plot

# save the plot
lapply(names(v_plots), function(x){ggsave(filename = paste0("/Users/grantdejong/Desktop/jon-test-plots/",x,".png"),  
                                          plot=v_plots[[x]],
                                          width = 3, height = 3, 
                                          units = "in",
                                          dpi = 800)})
