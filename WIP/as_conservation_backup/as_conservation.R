#!/usr/bin/env Rscript

##########################################################
#### R code for the quantification of AS conservation ####
##########################################################

## Remove redundancy
## Add overlap section to script with venn plots
## Add GO analysis

###################################
### Required packages and setup ###
###################################
rm( list= ls())
args <- commandArgs(trailingOnly = TRUE)

## Libraries
if(length(args) > 2) {
require(devtools, quietly = T, warn.conflicts = FALSE)
require(org.At.tair.db, quietly = T, warn.conflicts = FALSE)
require(plyr, quietly = T, warn.conflicts = FALSE)
require(dplyr, quietly = T, warn.conflicts = FALSE)
require(purrr, quietly = T, warn.conflicts = FALSE)
require(tidyr, quietly = T, warn.conflicts = FALSE)
require(reshape2, quietly = T, warn.conflicts = FALSE)
require(edgeR, quietly = T, warn.conflicts = FALSE)
require(data.table, quietly = T, warn.conflicts = FALSE)
require(grid, quietly = T, warn.conflicts = FALSE)
require(gridExtra, quietly = T, warn.conflicts = FALSE)
require(stats, quietly = T, warn.conflicts = FALSE)
require(forcats, quietly = T, warn.conflicts = FALSE)
require(ggplot2, quietly = T, warn.conflicts = FALSE)
require(car, quietly = T, warn.conflicts = FALSE)
require(emmeans, quietly = T, warn.conflicts = FALSE)
require(clusterProfiler, quietly = T, warn.conflicts = FALSE)
}
  
## Help section
if(length(args) < 1) {
  args <- c("--help")
}
if("--help" %in% args) {
  cat("
Alternative Splicing Analyses of Duplicated Genes
 
Arguments:
  --arg1=junctions.csv        - Table of equivalent homeologous junctions. Must be csv.
  --arg2=junction_counts.tsv  - Table of constitutive vs alternative read counts per 
                              junction. Must be tsv.
  --arg3=sampleinfo.tsv       - Sample info. Must be tsv. Two columns: 
                              (1) sample_data - a unique rep ID (2) samples_names - 
                              condition for each replicate.
  --arg4=groups_for_venn.txt  - OPTIONAL: groups for cross - condition comparisons.
  --help                      - print this text
 
Example:
./as_conservation.R --arg1=junctions.csv --arg2=as_data.tsv --arg3=sample_info.tsv --arg4=[OPTIONAL]groupinfo.txt \n\n\n")
  
  q(save="no")
}
if(length(args) > 4) {
  cat("Too many args. You need (1) junctions and (2) junction read counts (3) sample info.")
}

## Loading the data
print("Loading the data")
junctionids <- read.csv(args[1])
arsenal <- read.table(args[2], sep = '\t', header = TRUE); arsenal$species<-NULL
sampleinfo <- read.table(args[3], sep = '\t', header = T)
if(length(args)==4){
  groupinfo <- read.table(args[4], sep = '\t', header = T)
} else{print("No group-specific info supplied; avoiding cross-condition comparisons.")}

#################
### Functions ###
#################

## Ensure dplyr::select is the default
select <- dplyr::select

## Calculate PSI 
calculate_psi<-function(x){
  l=c()
  for (i in 1:(ncol(x)/2)){
    k=i*2; j=k-1
    l[[i]]<-x[j:k]}
  names(l)=sampleinfo$sample_data
  lapply(l, function(x){x[1]/(x[2]+x[1])})
}

## Filter by PSI
filter_psi<-function(x){lapply(x, function(x){
  x %>% filter_at(vars(contains("alt")), all_vars(.>0)) %>% 
    filter_at(vars(contains("alt")), all_vars(.<1)) %>% na.omit(.)
})}

## Filter by homeologs and ensure each event is equivalent
compare_junctions<-function(x){
  x=unite(x, junction, c("gene_name", "junct"), sep = ".", remove = T)
  a=x[x$junction %in% junctionids$junc1,]
  a=merge(a, junctionids, by.x="junction", by.y="junc1")
  y=merge(a, x, by.x="junc2", by.y="junction", all=F)
  y=y[which(y$class.x==y$class.y),]
  y$event=y$event.x; y$event.x=NULL; y$event.y=NULL
  invisible(y)
}
homologous_only<-function(y){
  if(length(y$unique_code.x=="na:na:+:IR")<1){
    cat("Warning:",deparse(substitute(x)),"has no alternate positions to compare.", sep = "\ ")
    invisible(rm(y))
  }
  else{
    x=y[!y$unique_code.x=="na:na:+:IR",]
    x=separate(x, unique_code.y, c("cy", "cy2"), sep = ":", remove = F) %>% separate(., unique_code.x, c("cx", "cx2"), sep = ":", remove = F)  %>% mutate(leny=as.numeric(.$cy2)-as.numeric(.$cy)) %>%
      mutate(lenx=as.numeric(.$cx2)-as.numeric(.$cx))
    for (i in 1:nrow(x)) {
      if(x$lenx[i] == x$leny[i]){x[i, 'homologous'] <- paste0("YES")}
      else if(x$lenx[i] == x$leny[i]+1){x[i, 'homologous'] <- paste0("YES")}
      else if(x$lenx[i] == x$leny[i]-1){x[i, 'homologous'] <- paste0("YES")}
      else{x[i, 'homologous'] <- paste0("NO")}}
  x=x[which(x$homologous=="YES"),]
  ir=y[y$unique_code.x=="na:na:+:IR",]
  invisible(y[y$unique_code.x %in% x$unique_code.x | y$unique_code.x=="na:na:+:IR",])
  }}

## Summarize qualitative conservation information
event_divergence<-function(x,y){
  d=x[x$junction %in% junctionids$junc1 | x$junction %in% junctionids$junc2,] %>% 
    unite(., event_id, c("junction", "class", "unique_code"), remove = F)
  c=unite(y, event_id_1, c("junction", "class.x", "unique_code.y"), remove = F) %>%
    unite(., event_id_2, c("junc2", "class.x", "unique_code.x"), remove = F)
  d[!(d$event_id %in% c$event_id_1 | d$event_id %in% c$event_id_2),] %>%
    .[!(.$class=="CREX" | .$class=="CPLX"),] %>% # Remove CREX for the time being
    filter_at(vars(contains("alt")), all_vars(.>0)) %>% 
    filter_at(vars(contains("nrm")), all_vars(.>0))
}
summary_conditional <- function(x,y){
  d=lapply(y, function(y){as.data.frame(summary(factor(y$class)))}) # Divergent summary
  c=lapply(x, function(x){as.data.frame(summary(factor(x$class.x)))}) # Conserved summary
  s=mapply(function(x,y){s=merge(x, y, all=T, by="row.names")}, d, c, SIMPLIFY = F)
  s=lapply(s, function(x){x[is.na(x)] = 0; invisible(return(x))})
  s=lapply(s, function(x){x$perc_cons=x[3]/(x[2]+x[3]); invisible(return(x))})
  s=lapply(s, function(x){colnames(x)=c("Event", "Divergent","Q_cons", "perc_cons"); invisible(return(x))})
  return(s)
}

## Quantitative divergence: Basic linear regression models **DEPRECATED**
lm_format<-function(x){unite(x, event_id, c("junction", "class.x", "unique_code.x")) %>% 
    unite(., event_id2, c("junc2", "class.y", "unique_code.y")) %>% 
    unite(., pair, c("event_id", "event_id2"), sep = "|") %>%
    gather(., key = "homeolog", value = "proportion", c(contains("alt.x"), contains("alt.y")), factor_key = T) %>% mutate(homeolog_id=gsub(".*x", 0, .$homeolog)) %>% 
    mutate(homeolog_id=gsub(".*y", 1, .$homeolog_id))}
lm_fit<-function(x){
  fit=with(x, 
           by(x, x[,"pair"], 
              function(x){lm(proportion ~ homeolog_id, data = x)}))
  t=sapply(fit, function(x){summary(x)$coefficients[2,4]})
  t.adj=p.adjust(t,method = "BH")
  length(t.adj[t.adj<0.05])/length(t.adj)
}
LR_tissues_fit<-function(x){
  with(x,
       by(x, x[,"pair"],
          function(x){
            if(length(levels(factor(x$condition)))>1){
              glm(logFC ~ condition, family=gaussian(link = "identity"), data = x)
            } else {
              rm(x)
            }}))
} # if loop makes sure single-tissue junctions are removed

## Quantitative divergence: Beta regression models **WIP** 
#(doesn't work on many genes with low sample size and near zero values; use on deltaPSI data?)
br_fit<-function(x){
  fit=with(x, 
           by(x, x[,"pair"], 
              function(x){fit=betareg(proportion ~ homeolog_id, data = x)}))
}
br_test_eq<-function(x){
  fit=with(x, 
           by(x, x[,"pair"], 
              function(x){fit=betareg(proportion ~ homeolog_id, data = x)
              eq.res=test(contrast(fit, "pairwise"), side = "=", adjust = "fdr", delta = 0.05)
              eq.res[1,6]}))
  p.adjust(fit, method = "BH")
}

## Quantitative divergence: Binomial logit regression models
BLR_format<-function(x){
  set=unite(x, event_id, c("junction", "class.x", "unique_code.x")) %>% 
    unite(., event_id2, c("junc2", "class.y", "unique_code.y")) %>% 
    unite(., pair, c("event_id", "event_id2"), sep = "|")
  # %>% mutate(id=seq(1,nrow(.)))
  alt=set %>% select(-contains("nrm")) %>% gather(key = "homeolog", value = "alt_counts", c(contains("alt.x"), contains("alt.y")), factor_key = T) %>%
    mutate(homeolog_id=gsub(".*x", 0, .$homeolog)) %>% 
    mutate(homeolog_id=gsub(".*y", 1, .$homeolog_id)) %>%
    mutate(sample=gsub("(.*?)_a.*$", "\\1",.$homeolog)) %>%
    mutate(condition=gsub("(.*?)_.*$", "\\1",.$homeolog)) %>%
    select(-homeolog)
  nrm=set %>% select(-contains("alt")) %>% gather(key = "homeolog", value = "nrm_counts", c(contains("nrm.x"), contains("nrm.y")), factor_key = T) %>%
    mutate(homeolog_id=gsub(".*x", 0, .$homeolog)) %>% 
    mutate(homeolog_id=gsub(".*y", 1, .$homeolog_id)) %>%
    mutate(sample=gsub("(.*?)_n.*$", "\\1",.$homeolog)) %>%
    mutate(condition=gsub("(.*?)_.*$", "\\1",.$homeolog)) %>%
    select(-homeolog)
  merge(alt,nrm,by=c("pair", "homeolog_id", "sample", "condition"))}
BLR_fit<-function(x){
  fit=with(x, 
           by(x, x[,"pair"],
              function(x){
                fit=glm(cbind(nrm_counts,alt_counts) ~ homeolog_id, family= binomial(link="logit"), data = x)}))
              }
BLR_test<-function(x){
  test=lapply(x, function(x){
    p=Anova(x, type = 3)
    p$`Pr(>Chisq)`[1]})
  p.adjust(test, method = "BH")
}
BLR_eq_test<-function(x){
  lapply(x, function(x){
    x=ref_grid(x)
    x=regrid(emmeans(x, "homeolog_id", type = "response")) %>%
      pairs(., reverse = T) %>% emmeans::test(., adjust = "fdr", side = "=", delta = 0.05,  df=4) %>% 
      .[1,6]
  })
}

## Muilti-group quantitative divergence: Binomial logit regression models
filter_groups<-function(x,y){
  l=list()
  for(i in names(y)){
    l[[i]]<-with(x, 
                 by(x, x[,"pair"],
                    function(x){
                      if(length(levels(factor(x$condition)))==length(y[[i]])){
                        if(levels(factor(x$condition))==y[[i]]){
                          invisible(return(x))
                        }else{
                          rm(x)} 
                      }else{
                        rm(x)}
                    }))
  }
  invisible(return(l))
}
BLR_fit_multiple<-function(x){
  fit=with(x, 
           by(x, x[,"pair"],
              function(x){
                if(length(levels(factor(x$condition)))>1){
                  fit=glm(cbind(nrm_counts,alt_counts) ~ homeolog_id*condition, family= binomial(link="logit"), data = x)
                } else {
                  rm(x)}
              }))}
BLR_fit_multiple_basic<-function(x){
  if(length(levels(factor(x$condition)))>1){
    fit=glm(cbind(nrm_counts,alt_counts) ~ homeolog_id*condition, family= binomial(link="logit"), data = x)
  } else {
    rm(x)}
} ## Remove?
eq_multi_conditional<-function(x){
  lapply(x, function(x){
    x=ref_grid(x)
    x=regrid(emmeans(x, ~ homeolog_id*condition, type = "response")) %>%
      pairs(., simple="homeolog_id", reverse = T) %>% emmeans::test(., adjust = "fdr", side = "=", delta = 0.05,  df=4) %>% 
      .[1,7] 
  })
}


## Misc purpose
conservation_estimate<-function(x){1-length(x[x<0.05])/length(x)}
set<-function(x){
  set=unite(x, event_id, c("junction", "class.x", "unique_code.x"), remove = F) %>% 
  unite(., event_id2, c("junc2", "class.y", "unique_code.y"), remove = F) %>% 
  unite(., pair, c("event_id", "event_id2"), sep = "|", remove = F)} # Produces pair vector
remove_null<-function(x){
  x[unlist(lapply(x, length) != 0)] 
} # Needed after BLR_fit_multiple

## Quantitative divergence: Mann-Whitney U test **DEPRECATED**
mann_whitney<-function(x){
  x$p_values<-""
  for (i in 1:nrow(x)) {
    data_a = x %>% select(.,contains("alt.x")) %>% .[i,] %>% as.numeric(.) # select alt(?).x?
    data_b = x %>% select(.,contains("alt.y")) %>% .[i,] %>% as.numeric(.) # select alt(?).y?
    x[i, 'p_values']<-wilcox.test(data_a, data_b, alternative = "two.sided" ,exact = F, correct = F)$p.value
  }
  return(x)
} # BH adjustment necessary but problematic with the small replicate number

## Summarize quantitative
conform<-function(x,y){
  df            <- merge(x, y, by= "row.names", all = T)
  rownames(df)  <- df$Row.names
  df$Row.names  <- NULL
  return(df)
}
summary_all<-function(x,d,c,n){
    q=x[sapply(x, function(x) dim(x)[1]) > 0] # Removes samples with missing information
    q=lapply(q, function(x){x %>% `rownames<-`(x$Event) %>% dplyr::select(-Event)}) # For merge
    d=lapply(d, function(d){as.data.frame(summary(factor(d$class.x)))}) # Divergent
    c=lapply(c, function(c){as.data.frame(summary(factor(c$class.x)))}) # Conserved
    n=lapply(n, function(n){as.data.frame(summary(factor(n$class.x)))}) # Non-sig events
    s=mapply(function(q,d,n,c){Reduce(conform ,list(q,d,n,c))}, q,d,n,c, SIMPLIFY = F)
    # Pipe the output to make it more concise?
    s=lapply(s, function(x){x[is.na(x)] = 0; invisible(return(x))})
    s=lapply(s, function(x){x$Perc_Cons=x[6]/(x[1]+x[2]); invisible(return(x))})
    s=lapply(s, function(x){colnames(x)=c("Divergent", "Qual. Cons.","Perc_Q_Cons", "Quant. Div.", "Unsupported", "Conserved", "Perc_Cons"); invisible(return(x))})
    return(s)
  }
plot_event_counts<-function(x){
  x=data.frame(event = rownames(x), div=x$Divergent, quant_div=x$`Quant. Div.`, quant_ns=x$Unsupported, quant_cons=x$Conserved) %>%
    melt(., id.var="event")

limit=sum(x$value)+300 # To ensure the y-axis is never fixed (as that can remove geoms)

# #7FA6B0 - lighter cons colour / #62929E darker cons colour
ggplot(x, aes(x=factor(event), y=as.numeric(value), fill=fct_rev(factor(variable))))+
  geom_bar(stat = "identity", width = .3)+
  scale_y_continuous(limits = c(0, limit))+
  scale_fill_manual(values=c("#7FA6B0", "#5C808E", "#546A7B","#C6C5B9"))+
  xlab("")+
  ylab("")+
  guides(fill=FALSE)+
  theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())+
  theme(axis.text.x=element_text(angle=90,hjust=1))}
generate_pdf<-function(x){
  mapply(function(x,y){
    pdf(file = paste0(y,"_quant_metrics.pdf"), height = 5, width = 12)
    grid.table(x)
    dev.off()}, x, names(x), SIMPLIFY = F)
}
generate_png<-function(x){
  mapply(function(x,y){
    png(file = paste0(y,"_quant_metrics.png"), height = 100, width = 600)
    grid.table(x)
    dev.off()}, x, names(x), SIMPLIFY = F)
}


#######################
### Gene expression ###
#######################
## Incorporate gene expression analysis (count data, DGE)
## Count data? Normalize via FPKM?
## Incoporate cufflinks??

if(any(grepl("GE=", args))){
# Load featureCounts data
count_dir=args[grep("GE=", args)] %>% as.character(.) %>% gsub("GE=(.*?)", "\\1", .)
samples=list.files(count_dir, recursive = T) %>% .[grep("*.txt$", .)]
# Load and organize count dataframe
count_data = data.frame(sample_names = samples %>% gsub("(.*?)/.*$", "\\1", .), samples = samples) %>% 
  mutate(unique_id = makeUnique(as.character(.$sample_names)))
counts = lapply(samples,function(x){
  x = read.table(paste0(count_dir, "/", x), sep = "\t", header = T) %>% .[,c(1,7)]
}) %>% reduce(left_join, by = c("Geneid")) %>% `rownames<-`(.$Geneid) %>% select(-Geneid)
colnames(counts) = count_data$unique_id

# edgeR

DGEList(counts=counts, group = count_data$sample_names, remove.zeros = T)

}

###########
### PSI ###
###########

# Processing dataframe and creating variables of interest
message("Processing dataframe and creating variables of interest...")

# Calculate PSI values
metadata <- subset(arsenal, , c(gene_name, junct, class, unique_code))
counts <- subset(arsenal, , -c(gene_name, junct, class, unique_code))

# Group lists and filter by PSI values
psi <- calculate_psi(counts) %>% as.data.frame(.)
  
psi_conditional <- lapply(levels(factor(sampleinfo$sample_names)), function(z){
  psi[grep(z, names(psi))] %>% as.data.frame(.)}) %>% 
  lapply(., function(x){cbind(x,metadata)}) %>%
  `names<-`(levels(factor(sampleinfo$sample_names))) %>% filter_psi(.)

## PSI homeolog processing
conserved_conditional <- lapply(psi_conditional, function(x){
  compare_junctions(x) %>% homologous_only(.)
})

#conserved_conditional <- conserved_conditional[-which(sapply(conserved_conditional, is.null))]
#?
#summary <- summary_conditional(conserved_conditional, psi_conditional)
#print(summary_conditional(conserved_conditional, psi_conditional))

################################
### Qualitative conservation ###
################################
message("Calculating degree of qualitative conservation...")

## Reads only
qual_conserved <- lapply(levels(factor(sampleinfo$sample_names)), function(z){
  arsenal[grep(z, names(arsenal))] %>% as.data.frame(.) %>% cbind(.,metadata)}) %>%
  `names<-`(levels(factor(sampleinfo$sample_names))) %>%
  lapply(., function(x){compare_junctions(x) %>% 
      homologous_only(.) %>%
      filter_at(vars(contains("alt")), all_vars(.>0)) %>% 
      filter_at(vars(contains("nrm")), all_vars(.>0)) %>%
      na.omit(.)}) # Find a way to filter cases for which nrm > alt

###### What uses qual conserved so you can rewrite the code to add pair names to it for overlap

divergent_events <- lapply(levels(factor(sampleinfo$sample_names)), function(z){
  arsenal[grep(z, names(arsenal))] %>% as.data.frame(.) %>% cbind(.,metadata) %>% 
    unite(., junction, c("gene_name", "junct"), sep = ".", remove = T)}) %>%
  `names<-`(levels(factor(sampleinfo$sample_names))) %>% 
  mapply(event_divergence, ., qual_conserved, SIMPLIFY = F)

summary <- summary_conditional(qual_conserved, divergent_events)
print(summary_conditional(qual_conserved, divergent_events))

#################################
### Quantitative conservation ###
#################################

## Binomial regression with logit link
BLR_log_fit<-lapply(qual_conserved, function(x){BLR_format(x) %>% BLR_fit(.)})
BLR_div<-lapply(BLR_log_fit, function(x){BLR_test(x)}) # Test filtering later
BLR_cons<-lapply(BLR_log_fit, function(x){BLR_eq_test(x)})
rm(BLR_log_fit)

## Print the results (WIP)
CONS_summary=lapply(BLR_cons, function(x){length(x[x<0.05])/length(x)})
DIV_summary=lapply(BLR_div, function(x){length(x[x<0.05])/length(x)})

div=lapply(BLR_div, function(x){x[x<0.05]})
cons=lapply(BLR_cons, function(x){x[x<0.05]}) # reassess the overlap between cons and div

## FIX THIS SUMMARY CODE
# useful debug function:
mapply(function(cons,div){length(which(names(cons) %in% names(div)))/length(cons)}, cons, div)
# FOR NOW filter div > cons:
cons_only=mapply(function(cons,div){cons[!(names(cons) %in% names(div))]}, cons, div)

## List of significantly divergent events 
BLR_div_list<-mapply(function(x,y){
  set=set(x) # adds a key vector to the qual_conserved dataframe (for filtering)
  set[set$pair %in% names(y),]}, qual_conserved, div, SIMPLIFY = F) %>%
  .[sapply(., function(.) dim(.)[1]) > 0]
## List of significantly conserved events
BLR_cons_list<-mapply(function(x,y){
  set=set(x)
  set[set$pair %in% names(y),]}, qual_conserved, cons_only, SIMPLIFY = F) %>% ###EDIT THIS LATER
  .[sapply(., function(.) dim(.)[1]) > 0] # Removes empty dataframes and formats list for tables/barplots
## List of non-significant but qualitatively conserved events
BLR_NS_list<-mapply(function(x,y,z){
  set=set(x)
  set[!(set$pair %in% names(y) | set$pair %in% names(z)),]}, qual_conserved, div, cons_only, SIMPLIFY = F) %>%
  .[sapply(., function(.) dim(.)[1]) > 0]


Quantitative <- summary_all(summary, BLR_div_list, BLR_cons_list, BLR_NS_list)
plots<-lapply(Quantitative, function(x){
  plot_event_counts(x)
})


####################################
### Cross - condition assessment ###
####################################
# Really slow... is there a more optimal way to write this?
# combines all junctions, sorts by condition, runs BLR on specified "group" sortings but only if group info is provided as a file.
if(exists("groupinfo")){
  qual_conserved_grouped<<-lapply(qual_conserved, function(x){BLR_format(x)}) %>% do.call(rbind,.)
  grouped_out<-suppressWarnings(filter_groups(qual_conserved_grouped,groupinfo)) %>% 
    lapply(., function(x){x[lengths(x) != 0] %>% do.call(rbind,.) %>%
        BLR_fit_multiple(.)}) #%>% # heck ya
  cons_group<<-lapply(grouped_out, eq_multi_conditional)
  div_group<<-lapply(grouped_out, BLR_test)
}else{message("You should really include a groupinfo.txt")
}



## Add this first - assess this based on overall tissue overlap (14-wise venn?)
## Can run quant analysis on these then correct for multi-testing afterwards
## filter based on interesting comparisons e.g.: 
## lapply(x, function(x){ x[which(levels(factor(x$condition)))==c("A", "B", "C"))), ] }
## These comparisons will have to be designed a priori

############################
### GEA - AS comparisons ###
############################

## Add this afterwards

#############
### Plots ###
#############
# Place this at the bottom of the doc when all else is completed?
# Do: (1) bar plots (2) violin plots (3) scatterplots?? 

## Condition-specific plots
lapply(names(plots), function(x){ggsave(filename = paste0(x,"CONSPLOT.png"),  
                                        plot=plots[[x]],
                                        width = 2, height = 3, 
                                        units = "in",
                                        dpi = 800)})

## Tissue-specific metrics pdf or png
lapply(Quantitative, function(x){
  x=data.frame(x$Divergent, x$`Qual. Cons.`, signif(x$Perc_Q_Cons, digits=3)*100, x$`Quant. Div.`, x$Unsupported, x$Conserved, signif(x$Perc_Cons, digits=3)*100, row.names = rownames(x))
  colnames(x)=c("Divergent", "Qual. Cons", "% symmetry", "Quant. Div.", "Unsupported", "Conserved", "% Conserved")
  invisible(return(x))}) %>% generate_png(.)

## Violin plots ##
# WIP*
if(any(grepl("plot", args))){
v_cons = BLR_cons_list %>% lapply(., function(x){BLR_format(x) %>% 
    mutate(proportion=alt_counts/(alt_counts+nrm_counts))}) %>%
  do.call(rbind,.) %>%
  with(., 
       by(., .[,"pair"],
          invisible(print)))

v_div = BLR_div_list %>% lapply(., function(x){BLR_format(x) %>% 
    mutate(proportion=alt_counts/(alt_counts+nrm_counts))}) %>%
  do.call(rbind,.) %>%
  with(., 
       by(., .[,"pair"],
          invisible(print)))


violin_plot<-function(x){
  ggplot(x, aes(x=condition))+
    geom_point(aes(y=proportion), colour="#669999ff", data=x[which(factor(x$homeolog_id)=="0"),])+
    geom_point(aes(y=proportion), colour="#00688bff", data=x[which(factor(x$homeolog_id)=="1"),])+
    ylim(0,1)+
    xlab("Root Tissue")+
    ylab("PSI")+
    guides(fill=FALSE)+
    theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())+
    theme(axis.text.x=element_text(angle=90,hjust=1))}

# make the plots
v_plots=lapply(head(v_div, 100), violin_plot) # for each plot
  
# save the plot

if(any(grepl("plot=", args))){
  plot_output=args[grep("plot=", args)] %>% as.character(.) %>% gsub("plot=(.*?)", "\\1", .)
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

###########
### WIP ###
###########

## deltaPSI **DEPRECATED COMPARISON**
#delta<-function(x,y){
#  a=x %>% select(contains(y))
#  a[1]-a[2]
#}
#na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
#delta_list<-function(x){
#  a=list()
#  names = x %>% select(contains("alt.x")) %>% names(.)
#  for(i in levels(sampleinfo$sample_data)){
#    a[[i]]<-ifelse(grepl(i, names), invisible(print(delta(x,i))), NA)
#  }
#  v_names = gsub("alt.x", "_FC", names)
#  a = a %>% na.omit.list(.) %>% lapply(., na.omit.list) %>% as.data.frame %>% 'names<-'(v_names)
#  cbind(x,a)
#}

## Wilcoxon *DEPCRECATED*
#W <- lapply(conserved_conditional, mann_whitney)
#W_summary <- lapply(W, function(x){length(x$p_values[x$p_values<0.05])/length(x$p_values)})
#W_list <- lapply(W, function(x){x[which(x$p_value<0.05),]})

## Beta regression - doesn't work well currently *WORK IN PROGRESS* *only use if no other option*
#BR<-lapply(conserved_conditional, function(x){
#  x= x %>% filter_at(vars(contains("alt")), all_vars(!is.na(.)))
#  x= x %>% filter_at(vars(contains("alt")), all_vars(.<1))
#  x=lm_format(x)
#br_set(x)
#})
#BR<-lapply(BR, function(x){
#  with(x, by(x, x[,"pair"], function(x){betareg(proportion ~ homeolog_id, data = x)}))})