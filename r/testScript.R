## test script
library(ggplot2)
library(tidyr)
setwd("~/GitHub/CommunityImmunity/r")
filename1 <- "/Users/fordfishman/GitHub/CommunityImmunity/comim_1000sims_m0_pS0.csv"
# filename2 <- "/Users/fordfishman/GitHub/CommunityImmunity/output/moregraphs1_richness.csv"
df1 <- read.csv(filename1, header = T, row.names = 1)
# df2 <- read.csv(filename2, header = T, row.names = 1)

df1$hostextinct <- df1$pop == 0

lm1 <- lm(data = df1, formula = hostextinct~log10(absp)+beta+l+a)
# 
summary(lm1)
# 
av <- aov(data = df1, formula = hostextinct~log10(absp)+beta+popinit+phageinit+l+a+c+b)
summary(av)

#histograms
# absorbance
ggplot(df1, aes(x=absp, fill=hostextinct)) +
  geom_histogram() +
  scale_x_log10() +
  scale_fill_grey("Host Extinct?", labels = c("No","Yes"))+
  theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 11),
    axis.ticks = element_line(size = 1),
    axis.ticks.length = unit(5,"pt"),)

# burst size
ggplot(df1, aes(x=beta, fill=hostextinct)) +
  geom_histogram() +
  # scale_x_log10() +
  scale_fill_grey("Host Extinct?", labels = c("No","Yes"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(5,"pt"),)

# decay rate
ggplot(df1, aes(x=d, fill=hostextinct)) +
  geom_histogram() +
  # scale_x_log10() +
  scale_fill_grey("Host Extinct?", labels = c("No","Yes"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(5,"pt"),)

# host birth rate
ggplot(df1, aes(x=b, fill=hostextinct)) +
  geom_histogram() +
  # scale_x_log10() +
  scale_fill_grey("Host Extinct?", labels = c("No","Yes"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(5,"pt"),)

# competition
ggplot(df1, aes(x=a, fill=hostextinct)) +
  geom_histogram() +
  scale_x_log10() +
  scale_fill_grey("Host Extinct?", labels = c("No","Yes"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(5,"pt"),)

# latency
ggplot(df1, aes(x=l, fill=hostextinct)) +
  geom_histogram() +
  # scale_x_log10() +
  scale_fill_grey("Host Extinct?", labels = c("No","Yes"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(5,"pt"),)

# initial populations
# absorbance
ggplot(df1, aes(x=popinit, y=phageinit, color=hostextinct)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual("Host Extinct?", labels = c("No","Yes"), values = c("dodgerblue4", "gold"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(5,"pt"),)

#################################################################################
# df1 <- gather(df1, key = organism, value=density, host:phage)
# df2 <- gather(df2, key = organism, value=richness, HostRichness:PhageRichness)
# 
# ggplot(df1, aes(x=time,y=density))+
#   geom_line(aes(color=organism)) +
#   scale_y_log10() +
#   scale_color_manual("", values = c("blue","red")) +
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
#   scale_color_manual("", values = c("blue","red")) +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         axis.text = element_text(size = 11),
#         axis.ticks = element_line(size = 1),
#         axis.ticks.length = unit(5,"pt"),) 
