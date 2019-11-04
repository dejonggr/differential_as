###################################
### Required packages and setup ###
###################################
rm( list= ls())
args <- commandArgs(trailingOnly = TRUE)

## Libraries
require(devtools, quietly = T, warn.conflicts = FALSE)
require(plyr, quietly = T, warn.conflicts = FALSE)
require(dplyr, quietly = T, warn.conflicts = FALSE)
require(tidyr, quietly = T, warn.conflicts = FALSE)
require(reshape2, quietly = T, warn.conflicts = FALSE)
require(data.table, quietly = T, warn.conflicts = FALSE)
require(grid, quietly = T, warn.conflicts = FALSE)
require(gridExtra, quietly = T, warn.conflicts = FALSE)
require(stats, quietly = T, warn.conflicts = FALSE)
require(forcats, quietly = T, warn.conflicts = FALSE)
require(ggplot2, quietly = T, warn.conflicts = FALSE)
require(car, quietly = T, warn.conflicts = FALSE)
require(betareg, quietly = T, warn.conflicts = FALSE)
require(lsmeans, quietly = T, warn.conflicts = FALSE)
# Missing venn

## Help section
# ADD OPTIONAL ARGUMENT TO SPECIFY THE STAT TEST MAYBE?
if(length(args) < 1) {
  args <- c("--help")
}
if("--help" %in% args) {
  cat("
      The R Script
      
      Arguments:
      --arg1=junctions.csv        - Table of equivalent homeologous junctions. Must be csv.
      --arg2=junction_counts.tsv  - Table of constitutive vs alternative read counts per junction. Must be tsv.
      --arg3=sampleinfo.tsv       - Sample info. Must be tsv. Two columns: (1) sample_data - a unique rep ID (2) samples_names - condition for each replicate.
      --help                      - print this text
      
      Example:
      ./test.R --arg1=junctions.csv --arg2=as_data.tsv \n\n --arg3=sample_info.tsv")
  
  q(save="no")
}
if(length(args) > 3) {
  cat("Too many args. You need (1) junctions and (2) junction read counts (3) sample info.")
}

## Loading the data
print("Loading the data")
junctionids <- read.csv(args[1])
arsenal <- read.table(args[2], sep = '\t', header = TRUE); arsenal$species<-NULL
sampleinfo <- read.table(args[3], sep = '\t', header = T)

#################
### Functions ###
#################

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

## Rsq for model testing **DEPRECATED**
Rsq <- function( model ){
  fitted.variance <- var(model$fitted)
  total.variance	<- var(model$fitted) + var(model$resid)
  fitted.variance / total.variance
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
BLR_fit_multiple<-function(x){
  fit=with(x, 
           by(x, x[,"pair"],
              function(x){
                fit=glm(cbind(nrm_counts,alt_counts) ~ homeolog_id*condition, family= binomial(link="logit"), data = x)
                p=Anova(fit, type = 3)
                p$`Pr(>Chisq)`[3]}))
  p.adjust(fit,method = "BH")
}

## Misc purpose
conservation_estimate<-function(x){1-length(x[x<0.05])/length(x)}
set<-function(x){
  set=unite(x, event_id, c("junction", "class.x", "unique_code.x"), remove = F) %>% 
    unite(., event_id2, c("junc2", "class.y", "unique_code.y"), remove = F) %>% 
    unite(., pair, c("event_id", "event_id2"), sep = "|", remove = F)}

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
  q=lapply(q, function(x){x %>% `rownames<-`(x$Event) %>% select(-Event)}) # For merge
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
###########
### PSI ###
###########

# Processing dataframe and creating variables of interest
print("Processing dataframe and creating variables of interest...")

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

conserved_conditional <- conserved_conditional[-which(sapply(conserved_conditional, is.null))]

#summary <- summary_conditional(conserved_conditional, psi_conditional)
#print(summary_conditional(conserved_conditional, psi_conditional))

#######################################
### Tissue-specific AS conservation ###
#######################################
## Overlapped datasets

C=lapply(qual_conserved, set)
QC=BLR_cons_list

## Root-tip; QC, HC bad coverage ##

# Qualitiative overlap
v_rt<-venn(list(COL=C$COL$pair, STL=C$STL$pair, LRC=C$LRC$pair, QC=C$QC$pair), ellipse = T, zcolor = "saddlebrown, mediumpurple4, lightseagreen, palegoldenrod")
# venn for printing; v_rt<-venn(list(COL=C$COL$pair, STL=C$STL$pair, LRC=C$LRC$pair, QC=C$QC$pair), ellipse = T, zcolor = "saddlebrown, mediumpurple4, lightseagreen, palegoldenrod", size = 25, cexil = 3, cexsn = 3, borders = F)


# Quantitative overlap
qv_rt<-venn(list(COL=QC$COL$pair, STL=QC$STL$pair, LRC=QC$LRC$pair, QC=QC$QC$pair), ellipse = T, zcolor = "style")
qv_rt.overlap<-calculate.overlap(list("COL"=QC$COL$pair, "STL"=QC$STL$pair, "LRC"=QC$LRC$pair, "QC"=QC$QC$pair))

# Root-tip qualitative overlap list for quantitative analysis
v_rt.overlap<-calculate.overlap(list("COL"=C$COL$pair, "STL"=C$STL$pair, "LRC"=C$LRC$pair, "QC"=C$QC$pair))
v_rt.all<-v_rt.overlap$a6

# Quantitative analysis
q_rt<-list(COL=qual_conserved$COL, STL=qual_conserved$STL, LRC=qual_conserved$LRC, QC=qual_conserved$QC) %>%
  lapply(., function(x){
    x=BLR_format(x)
    x=x[x$pair %in% v_rt.all,]
  }) %>%
  do.call(rbind, .)
q_rt_out<-BLR_fit_multiple(q_rt)
q_rt_cons<-eq_multi_conditional(q_rt_out) %>% unlist(.) %>% p.adjust(., method = "BH")
q_rt_cons_list<-q_rt_cons[which(q_rt_cons<0.05)]
# High conservation %...

# for testing
test=q_rt
test$psi=test$alt_counts/(test$alt_counts+test$nrm_counts)
cons_psi=test[test$pair %in% names(q_rt_cons_list),]


# venn quantitative overlap testing
qtest=test[test$pair %in% qv_rt.overlap$a6,]
ggplot(qtest, aes(x=psi, fill=condition))+
  geom_density(alpha=0.2)

############################################
### Tissue-specific AS rate conservation ###
############################################

FC<-function(x,y){
  a=x %>% select(contains(y))
  log(a[1]/a[2])
}
na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
log_FC_list<-function(x){
  a=list()
  names = x %>% select(contains("alt.x")) %>% names(.)
  for(i in levels(sampleinfo$sample_data)){
    a[[i]]<-ifelse(grepl(i, names), invisible(print(FC(x,i))), NA)
  }
  v_names = gsub("alt.x", "FC", names)
  a = a %>% na.omit.list(.) %>% lapply(., na.omit.list) %>% as.data.frame %>% 'names<-'(v_names)
  cbind(x,a)
}

FC_format<-function(x){unite(x, event_id, c("junction", "class.x", "unique_code.x")) %>% 
    unite(., event_id2, c("junc2", "class.y", "unique_code.y")) %>% 
    unite(., pair, c("event_id", "event_id2"), sep = "|") %>% 
    gather(., key = "sample", value = "logFC", contains("FC"), factor_key = T) %>%
    mutate(condition=as.factor(gsub("(.*?)_.*$", "\\1",.$sample)))}

conserved_conditional_FC<-lapply(conserved_conditional, function(x) log_FC_list(x) %>% select(-contains("alt")))

list_pairs<-lapply(conserved_conditional_FC, function(x){
  x=set(x)
  list(x$pair)})

calculate.overlap(list_pairs)

## Across all conditions ##
# The resulting data frame includes all possible tissue overlaps, difficult to draw conclusions from
# but it is comprehensive 
conserved_conditional_FC_formatted<-lapply(conserved_conditional_FC, function(x) FC_format(x)) %>%
  do.call(rbind, .)
tissue_as_divergence<-LR_tissues_fit(conserved_conditional_FC_formatted) %>%
  .[!sapply(., is.null)] # Removes null (a.k.a single-tissue) junctions from dataset

# do an equivalence test (think of delta)
# compute overlaps for every cell-type then organize junctions by these designations

test=BLR_test(tissue_as_divergence)
test_eq=fceq_test(tissue_as_divergence)
cons_fc=test_eq[which(test_eq<0.05)] %>% unlist(.)


# to inspect the "conserved" FC patterns
psi_formatted=lapply(conserved_conditional, function(x){x=lm_format(x) %>% mutate(homeolog=gsub("(.*?)_.*$", "\\1", .$homeolog))}) %>%
  do.call(rbind,.) %>%
  with(., 
        by(., .[,"pair"],
          invisible(print)))
# you could realistically filter this with anything, it doesn't matter if it contains all 
# conditions as long as they're labelled
psi_FC=psi_formatted[names(psi_formatted) %in% names(cons_fc)]

## violin plot ##

violin_plot<-function(x){
  ggplot(x, aes(x=homeolog))+
    geom_violin(aes(y=proportion), fill="#669999ff", trim=FALSE, data=x[which(factor(x$homeolog_id)=="0"),])+
    geom_violin(aes(y=proportion), fill="#00688bff", trim=FALSE, data=x[which(factor(x$homeolog_id)=="1"),])+
    ylim(0,1)+
    xlab("Root Tissue")+
    ylab("PSI")+
    guides(fill=FALSE)+
    theme_bw() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())+
    theme(axis.text.x=element_text(angle=90,hjust=1))}

# make the plots
v_plots=lapply(psi_FC, violin_plot) # for each plot

# save the plot
lapply(names(v_plots), function(x){ggsave(filename = paste0("~/Desktop/PSI_GRAPHS/",x,".png"),  
                                        plot=v_plots[[x]],
                                        width = 3, height = 3, 
                                        units = "in",
                                        dpi = 800)})


# test to see how many contain mer, elo, mat
lapply(psi_formatted, function(x){homeolog})  



#########

pull=conserved_conditional_FC_formatted[which(conserved_conditional_FC_formatted$pair %in% names(cons_fc)),] %>%
  .[abs(.$logFC)>1 & abs(.$logFC)<1.5,]

t.overlap=lapply(overlapping_qualitatively, function(x){
  x[which(x %in% names(cons_fc))]})



fuck=lapply(conserved_conditional, function(x){
  x=set(x)
  x[x$pair %in% unlist(t.overlap),]})

fuck$COL$COL_1_alt.x

fceq_test<-function(x){
  lapply(x, function(x){
    x=ref_grid(x)
    x=regrid(emmeans(x, "condition", type = "response")) %>%
      pairs(., reverse = T) %>% emmeans::test(., adjust = "fdr", side = "=", delta = 1,  df=4) %>% 
      .[1,6]
  })
}

ggplot(x, aes(x$COL_1_alt.x))+
  geom_density()
# range -5 to 5, good cons measure 5%? so 0.5?











## Across specific conditions groups
foo<-lapply(list("insert conserved_conditional_FC$CONDITION here"), function(x) FC_format(x)) %>%
  do.call(rbind, .)

## FIND A WAY TO MELT EACH ELEMENT OF THE LIST DOWN INTO CORRESPONDING TISSUE TYPES
## You need to work overlaps... which simply can't be done with this data...
# sapply maybe?
