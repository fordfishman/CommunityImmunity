## Ford Fishman

library("ggplot2")
library("tidyr")
library("optparse")
library("scales")
library("MASS")

option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", 
              help="output directory ", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
# the datafile
filename <- opt$file
output_dir <- opt$out

df <- read.csv(file=filename, header = T, row.names = 1)

df$initPopPhageRatio <- df$popinit/df$phageinit

df$hostextinct <- df$pop == 0
df$dominant <- with(df, ifelse(hostextinct, "Extinct Host",ifelse(vulnerable<=immune,"Immune","Vulnerable")) )
df$coexistence <- with(df, 
                        ifelse(hostextinct, "Extinct Host",
                               ifelse(vulnerable>0 & immune>0 & phage >0,"Full Coexistence",
                                      ifelse(vulnerable >0 & immune>0,"Host Coexistence",
                                             ifelse(vulnerable >0 & phage>0,"Vulnerable + Phage",
                                                    ifelse(immune>0 & phage>0,"Immune + Phage",
                                                           ifelse(immune>0,"Immune","Vulnerable")))))) )

# organize histogram stack
df$coexistence <- factor(df$coexistence, levels = c("Extinct Host","Full Coexistence", "Host Coexistence","Immune + Phage","Immune","Vulnerable + Phage","Vulnerable"))

# plotting functions

dominancePlot <- function(df, col, xname=NULL, xlog=F, bins=40){
  if(is.null(xname)){xname = col}
  if(xlog){
    x <- with(df, log10(get(col)))
  } else{
    x <- with(df, get(col))
  }
  p <- ggplot(data=df, aes(x=x, fill=dominant)) +
    geom_histogram(position = "fill",bins=bins,color="black", size = 0.2) +
    scale_y_continuous("Proportion",expand=c(0,0)) +
    scale_fill_manual("", limits = c("Extinct Host","Immune","Vulnerable"),values = c("#666666","palegreen","steelblue2"),drop = F )+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 11),
          axis.ticks = element_line(size = 1),
          axis.ticks.length = unit(5,"pt"),)
  if (xlog){
    p <- p +scale_x_continuous(xname,expand=c(0,0),labels =  math_format(10^.x))
  } else {
    p <- p +scale_x_continuous(xname,expand=c(0,0))
  }
  return(p)
}

coexistencePlot <- function(df, col, xname=NULL, xlog=F, bins=40){
  if(is.null(xname)){xname = col}
  if(xlog){
    x <- with(df, log10(get(col)))
  } else{
    x <- with(df, get(col))
  }
  p <- ggplot(data=df, aes(x=x, fill=coexistence)) +
    geom_histogram(position = "fill",bins=bins,color="black", size = 0.2) +
    scale_y_continuous("Proportion",expand=c(0,0)) +
    scale_fill_manual("", values = c("#666666","gold",  "orange3", "#1B9E77","palegreen", "#7570B3", "steelblue2"), drop = F )+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size = 11),
          axis.ticks = element_line(size = 1),
          axis.ticks.length = unit(5,"pt"),)
  if (xlog){
    p <- p +scale_x_continuous(xname,expand=c(0,0),labels=math_format(10^.x))
  } else {
    p <- p +scale_x_continuous(xname,expand=c(0,0))
  }
  return(p)
}



#histograms

# prob. of forming spacer
pSxName <- "Prob. of spacer formation (per infection)"

pSDom <- dominancePlot(df,"pS",xname=pSxName,xlog=T)
pSDom_dir <- paste(output_dir, "pS_Dom.png", sep="/")
ggsave(plot = pSDom, filename = pSDom_dir, width = 7, height = 5) 

pSCoex <- coexistencePlot(df,"pS",xname=pSxName,xlog=T)
pSCoex_dir <- paste(output_dir, "pS_Coex.png", sep="/")
ggsave(plot = pSCoex, filename = pSCoex_dir, width = 7, height = 5)

# mutation
mxName <- "Prob. of phage mutation (per nucleotide)"

mDom <- dominancePlot(df,"m",mxName,xlog=T)
mDom_dir <- paste(output_dir, "m_Dom.png",sep="/")
ggsave(plot = mDom, filename = mDom_dir, width = 7, height = 5)

mCoex <- coexistencePlot(df,"m",mxName,xlog=T)
mCoex_dir <- paste(output_dir, "m_Coex.png", sep="/")
ggsave(plot = mCoex, filename = mCoex_dir, width = 7, height = 5)

# absorbance
abspxName <- "Absorption Rate (per host-phage interaction)"

abspDom <- dominancePlot(df, "absp", abspxName, xlog=T)
abspDom_dir <- paste(output_dir, "absp_Dom.png",sep="/")
ggsave(plot = abspDom, filename = abspDom_dir, width = 7, height = 5)

abspCoex <- coexistencePlot(df, "absp", abspxName, xlog=T)
abspCoex_dir <- paste(output_dir, "absp_Coex.png", sep="/")
ggsave(plot = abspCoex, filename = abspCoex_dir, width = 7, height = 5)

# burst size
betaxName <- "Burst Size (# of phage)"

betaDom <- dominancePlot(df, "beta", betaxName, xlog=F)
betaDom_dir <- paste(output_dir, "beta_Dom.png", sep="/")
ggsave(plot = betaDom, filename = betaDom_dir, width = 7, height = 5)

betaCoex <- coexistencePlot(df, "beta", betaxName, xlog=F)
betaCoex_dir <- paste(output_dir, "beta_Coex.png", sep="/")
ggsave(plot = betaCoex, filename = betaCoex_dir, width = 7, height = 5)

# decay rate
dxName <- "Phage decay rate (per capita)"

dDom <- dominancePlot(df, "d", dxName, xlog=F)
dDom_dir <- paste(output_dir, "d_Dom.png", sep="/")
ggsave(plot = dDom, filename = dDom_dir, width = 7, height = 5)

dCoex <- coexistencePlot(df, "d", dxName, xlog=F)
dCoex_dir <- paste(output_dir, "d_Coex.png", sep="/")
ggsave(plot = dCoex, filename = dCoex_dir, width = 7, height = 5)

# host birth rate
bxName <- "Host Birth Rate (per capita)"

bDom <- dominancePlot(df,"b",bxName)
bDom_dir <- paste(output_dir, "b_Dom.png", sep="/")
ggsave(plot = bDom, filename = bDom_dir, width = 7, height = 5)

bCoex <- coexistencePlot(df,"b", bxName)
bCoex_dir <- paste(output_dir, "b_Coex.png", sep="/")
ggsave(plot = bCoex, filename = bCoex_dir, width = 7, height = 5)


# competition
axName <- 'Host "Carrying Capacity"'

aDom <- dominancePlot(df,"a",axName, xlog = T)
aDom_dir <- paste(output_dir, "a_Dom.png", sep="/")
ggsave(plot = aDom, filename = aDom_dir, width = 7, height = 5)

aCoex <- coexistencePlot(df,"a",axName, xlog = T)
aCoex_dir <- paste(output_dir, "a_Coex.png", sep="/")
ggsave(plot = aCoex, filename = aCoex_dir, width = 7, height = 5)

# latency
lxName <- "Latency"

lDom <- dominancePlot(df,"l",lxName)
lDom_dir <- paste(output_dir, "l_Dom.png", sep="/")
ggsave(plot = lDom, filename = lDom_dir, width = 7, height = 5)

lCoex <- coexistencePlot(df,"l",lxName)
lCoex_dir <- paste(output_dir, "l_Coex.png", sep="/")
ggsave(plot = lCoex, filename = lCoex_dir, width = 7, height = 5)

# cost of CRISPR
cxName <- "Cost of CRISPR"

cDom <- dominancePlot(df,"c",cxName, xlog = T)
cDom_dir <- paste(output_dir, "c_Dom.png", sep="/")
ggsave(plot = cDom, filename = cDom_dir, width = 7, height = 5)

cCoex <- coexistencePlot(df,"c",cxName, xlog = T)
cCoex_dir <- paste(output_dir, "c_Coex.png", sep="/")
ggsave(plot = cCoex, filename = cCoex_dir, width = 7, height = 5)

