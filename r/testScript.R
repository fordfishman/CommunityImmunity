## test script
library(ggplot2)
library(tidyr)
setwd("~/GitHub/CommunityImmunity/r")
filename1 <- "/Users/fordfishman/GitHub/CommunityImmunity/output/moregraphs1_main.csv"
filename2 <- "/Users/fordfishman/GitHub/CommunityImmunity/output/moregraphs1_richness.csv"
df1 <- read.csv(filename1, header = T, row.names = 1)
df2 <- read.csv(filename2, header = T, row.names = 1)


df1 <- gather(df1, key = organism, value=density, host:phage)
df2 <- gather(df2, key = organism, value=richness, HostRichness:PhageRichness)

ggplot(df1, aes(x=time,y=density))+
  geom_line(aes(color=organism)) +
  scale_y_log10() +
  scale_color_manual("", values = c("blue","red")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(5,"pt"),) 

ggplot(df2, aes(x=time,y=richness))+
  geom_line(aes(color=organism)) +
  scale_y_continuous() +
  scale_color_manual("", values = c("blue","red")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(5,"pt"),) 
