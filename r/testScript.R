## test script
library(ggplot2)
library(tidyr)
library(scales)
library(MASS)
library(here)

wd <- here()
output_dir <- paste0(wd, "/output/testingGrounds/")
filename1 <- paste0(wd, "/output/1000sims/summary.csv")
# filename1 <- "/Users/fordfishman/GitHub/CommunityImmunity/output/RED_test_1000sims.csv"
# filename2 <- "/Users/fordfishman/GitHub/CommunityImmunity/output/RED_test_1000sims_absp1e-08_beta100_c0_d0.1_m1e-07_pS1e-05.csv"
# filename3 <- "/Users/fordfishman/GitHub/CommunityImmunity/output/RED_test_1000sims_a10000000.0_c0_d0.1_l0.9_m1e-07_pS1e-05.csv"
df1 <- read.csv(filename1, header = T, row.names = 1)
# df2 <- read.csv(filename2, header = T, row.names = 1)


df1$initPopPhageRatio <- df1$popinit/df1$phageinit

df1$hostextinct <- df1$pop == 0
df1$dominant <- with(df1, ifelse(hostextinct, "Extinct Host",ifelse(vulnerable<=immune,"Immune","Vulnerable")) )
df1$coexistence <- with(df1, 
                        ifelse(hostextinct, "Extinct Host",
                               ifelse(vulnerable>0 & immune>0 & phage >0,"Full Coexistence",
                                      ifelse(vulnerable >0 & immune>0,"Host Coexistence",
                                             ifelse(vulnerable >0 & phage>0,"Vulnerable + Phage",
                                                    ifelse(immune>0 &phage>0,"Immune + Phage",
                                                           ifelse(immune>0,"Immune",
                                                                  "Vulnerable")))))) )


# organize plot
df1$coexistence <- factor(df1$coexistence, levels = c("Extinct Host","Full Coexistence", "Host Coexistence","Immune + Phage","Immune","Vulnerable + Phage","Vulnerable"))

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
    scale_fill_manual("", values = c("#666666","palegreen","steelblue2") )+
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
    scale_fill_manual("", values = c("#666666","gold",  "orange3", "#1B9E77","palegreen", "#7570B3", "steelblue2") )+
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



# df1.pca <- prcomp(df1[,7:18], scale. = T)
# 
# library(factoextra)
# 
# fviz_eig(df1.pca)
# 
# pca.plot <- fviz_pca_biplot(df1.pca,
#                          col.ind = df1$coexistence,
#                          col.var = "grey",
#                          repel = T)
# pca.plot

#histograms

# prob. of forming spacer
pSxName <- "Prob. of spacer formation (per infection)"

pSDom <- dominancePlot(df1,"pS",xname=pSxName,xlog=T)
pSDom_dir <- paste(output_dir, "pS_Dom.png", sep="/")
ggsave(plot = pSDom, filename = pSDom_dir, width = 7, height = 5) 

pSCoex <- coexistencePlot(df1,"pS",xname=pSxName,xlog=T)
pSCoex_dir <- paste(output_dir, "pS_Coex.png", sep="/")
ggsave(plot = pSCoex, filename = pSCoex_dir, width = 7, height = 5)

# mutation
mxName <- "Prob. of phage mutation (per nucleotide)"

mDom <- dominancePlot(df1,"m",mxName,xlog=T)
mDom_dir <- paste(output_dir, "m_Dom.png",sep="/")
ggsave(plot = mDom, filename = mDom_dir, width = 7, height = 5)

mCoex <- coexistencePlot(df1,"m",mxName,xlog=T)
mCoex_dir <- paste(output_dir, "m_Coex.png", sep="/")
ggsave(plot = mCoex, filename = mCoex_dir, width = 7, height = 5)

# absorbance
abspxName <- "Absorption Rate (per host-phage interaction)"
dominancePlot(df1, "absp", abspxName, xlog=T)
coexistencePlot(df1, "absp", abspxName, xlog=T)

# burst size
betaxName <- "Burst Size (# of phage)"
dominancePlot(df1, "beta", betaxName, xlog=F)
coexistencePlot(df1, "beta", betaxName, xlog=F)

# decay rate
dxName <- "Phage decay rate (per capita)"
dominancePlot(df1, "d", dxName, xlog=F)
coexistencePlot(df1, "d", dxName, xlog=F)

# host birth rate
bxName <- "Host Birth Rate (per capita)"

bDom <- dominancePlot(df1,"b",bxName)
bDom_dir <- paste(output_dir, "b_Dom.png", sep="/")
ggsave(plot = bDom, filename = bDom_dir, width = 7, height = 5)

bCoex <- coexistencePlot(df1,"b", bxName)
bCoex_dir <- paste(output_dir, "b_Coex.png", sep="/")
ggsave(plot = bCoex, filename = bCoex_dir, width = 7, height = 5)


# competition
axName <- 'Host "Carrying Capacity"'

aDom <- dominancePlot(df1,"a",axName, xlog = T)
aDom_dir <- paste(output_dir, "a_Dom.png", sep="/")
ggsave(plot = aDom, filename = aDom_dir, width = 7, height = 5)

aCoex <- coexistencePlot(df1,"a",axName, xlog = T)
aCoex_dir <- paste(output_dir, "a_Coex.png", sep="/")
ggsave(plot = aCoex, filename = aCoex_dir, width = 7, height = 5)

# latency
lxName <- "Latency"
dominancePlot(df1,"l",lxName)
coexistencePlot(df1,"l",lxName)

# cost of CRISPR
cxName <- "Cost of CRISPR"

cDom <- dominancePlot(df1,"c",cxName, xlog = T)
cDom_dir <- paste(output_dir, "c_Dom.png", sep="/")
ggsave(plot = cDom, filename = cDom_dir, width = 7, height = 5)

cCoex <- coexistencePlot(df1,"c",cxName, xlog = T)
cCoex_dir <- paste(output_dir, "c_Coex.png", sep="/")
ggsave(plot = cCoex, filename = cCoex_dir, width = 7, height = 5)


# initial populations
initxName <- "Initial Population Ratio (Host:Phage)"
dominancePlot(df1,"initPopPhageRatio",initxName, xlog = T)
coexistencePlot(df1,"initPopPhageRatio",initxName, xlog = T)

ggplot(df1, aes(x=a,y=c,color = coexistence))+
  geom_point()+
  scale_x_log10() +
  scale_y_log10() 
  # scale_color_manual("", values = c("#666666","palegreen","steelblue2"))

#################################################################################
# df1 <- gather(df1, key = organism, value=density, host:phage)
# df2 <- gather(df2, key = organism, value=richness, HostRichness:PhageRichness)
# 
# ggplot(df1, aes(x=time,y=density))+
#   geom_line(aes(color=organism)) +
#   scale_y_log10() +
#   scale_color_manual("", values = c("steelblue1","firebrick")) +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         axis.text = element_text(size = 11),
#         axis.ticks = element_line(size = 1),
#         axis.ticks.length = unit(5,"pt"),) 
# 
# ggplot(df2, aes(x=time,y=richness))+
#   geom_line(aes(color=organism)) +
#   scale_y_continuous() +
#   scale_color_manual("", values = c("steelblue1","firebrick")) +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         axis.text = element_text(size = 11),
#         axis.ticks = element_line(size = 1),
#         axis.ticks.length = unit(5,"pt"),) 
