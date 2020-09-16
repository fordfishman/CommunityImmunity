## test script
library(ggplot2)
library(tidyr)
library(scales)
library(MASS)
library(here)
wd <- here()

# import file
# filename1 <- paste0(wd, "/output/1sims_a1e6_adsp1e-8_b1.2_beta100_d0.1_l0.9_m1e-7_pS1e-6/main.csv")
# filename1 <- paste0(wd, "/output/1sims_t5000_a1e6_adsp1e-8_b1.2_beta100_d0.1_l0.9_m1e-7_pS1e-6/main.csv")
filename1 <- paste0(wd, "/output/1sims_t2000_a1e6_adsp1e-8_b1.2_beta100_c0.01_d0.1_f0_l0.9_m1e-7_pS1e-6_pop1e5_phage1e7/main.csv")
filename2 <- paste0(wd, "/output/1sims_t2000_a1e6_adsp1e-8_b1.2_beta100_c0.01_d0.1_f0_l0.9_m1e-7_pS1e-6_pop1e5_phage1e7/richness.csv")

df1 <- read.csv(filename1, header = T, row.names = 1)
df2 <- read.csv(filename2, header = T, row.names = 1)
df2 <- gather(df2, key = organism, value=richness, HostRichness:PhageRichness)

stepNum <- length( unique( df1$timestep ) ) # total number of timesteps in simulation

InitialHost <- df1$pop[df1$name=="s1"]

InitialHost <- c(InitialHost, rep(0, 2000 - length(InitialHost)))


totalAbundance <- data.frame( # initialize secondary data frame for totals across different categories
  
  Timestep=1:stepNum, 
  Host=rep(0,stepNum), 
  Phage=rep(0,stepNum), 
  InitialHost=InitialHost
)

ind <- 1

# add up total host and phage density at each timestep
for (i in min(df1$timestep):max(df1$timestep)){
  
  df <- subset(df1, timestep==i)
  
  totalAbundance[ind,"Host"] <- with(df, sum(pop[type!="phage"]))
  
  totalAbundance[ind,"Phage"] <- with(df, sum(df$pop[type=="phage"]))
  
  ind <- ind + 1
  
}
totalAbundance$Variant <- totalAbundance$Host - totalAbundance$InitialHost # total of spacer variants at each timestep
totalAbundance2 <- gather(totalAbundance, key = type, value = abundance, -Timestep) # the data in long format

totalAbundance$Host_log <- log10(totalAbundance$Host)
totalAbundance$Phage_log <- log10(totalAbundance$Phage)
totalAbundance2$abundance <- log10(totalAbundance2$abundance)

# df1$dpop_pop <- with(df1, ifelse(dpop==0, NA, dpop/pop))

totalAbundance2$type <- factor(totalAbundance2$type, levels = c("Host","InitialHost","Variant","Phage")) # change legend order

maxPop <- log10(max(df1$pop)) # max host population

## plot of the sums across groups
ggplot(subset(totalAbundance2, type!="Variant"), aes(x=Timestep,y=abundance, group = type)) +
  geom_line(aes(color=type)) +
  scale_y_continuous("Density",labels = math_format(10^.x), expand = c(0,0)) +
  scale_x_continuous("Timestep")+
  scale_color_manual("",values = c("black","blue","green","red"),
                     labels = c("Total Host","Initial Host","Spacer Variant","Phage")) +
  annotation_logticks(sides = "l") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        axis.ticks.length = unit(5,"pt"),)

# all strains: shows coevolution
ggplot(df1, aes(x=timestep,y=log10(pop), group = name))+
  geom_line(aes(color=type)) +
  # scale_y_continuous("Density", limits = c(1, 0.5+maxPop),labels =  math_format(10^.x), ) +
  scale_y_continuous("Density", breaks = seq(0,10,by = 2),labels = math_format(10^.x), expand = c(0,0)) +
  scale_x_continuous("Timestep")+
  scale_color_manual("", labels = c("Initial Strain","Spacer Variant","Phage"),values = c("blue","green","red")) +
  annotation_logticks(sides = "l") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        axis.ticks.length = unit(5,"pt"),)

ggplot(df1, aes(x=timestep,y=(spacers), group = name))+
  geom_line(aes(color=type)) +
  # scale_y_continuous("Density", limits = c(1, 0.5+maxPop),labels =  math_format(10^.x), ) +
  scale_y_continuous("# of spacers", limits=c(0,4)) +
  scale_x_continuous("Timestep")+
  scale_color_manual("", labels = c("Initial Strain","Spacer Variant","Phage"),values = c("blue","green","red")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        axis.ticks.length = unit(5,"pt"),)
ggplot(totalAbundance, aes(x=Host_log,y=Phage_log))+
  geom_line() +
  geom_point(aes(x = totalAbundance[stepNum,"Host_log"], totalAbundance[stepNum,"Phage_log"]), color = "blue") +
  # scale_y_continuous("Density", limits = c(1, 0.5+maxPop),labels =  math_format(10^.x), ) +
  scale_y_continuous("Phage Density",
                     labels = math_format(10^.x), ) +
  scale_x_continuous("Host Density",
                     labels = math_format(10^.x),
                     )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(5,"pt"),)

# ratio of phage to total host over time
ggplot(totalAbundance, aes(x=1:stepNum,y=(Phage/Host)))+
  geom_line() +  # scale_y_continuous("Density", limits = c(1, 0.5+maxPop),labels =  math_format(10^.x), ) +
  scale_y_continuous("Phage/Host", expand = c(0,0)) +
  scale_x_continuous("Timestep")+
  # scale_color_manual("", labels = c("Initial Strain","Spacer Variant","Phage"),values = c("blue","green","red")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(5,"pt"),)
# 
ggplot(subset(df1, timestep>50 & type=="phage"), aes(x=timestep,y=(dpop), group = name))+
  geom_line(aes(color=type)) +
  # scale_y_continuous("Density",labels = math_format(10^.x), ) +
  scale_y_continuous() +
  scale_x_continuous("Timestep")+
  scale_color_manual("", labels = c("Initial Strain","Spacer Variant","Phage"),values = c("blue","green","red")) +
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
  scale_color_manual("", values = c("steelblue1","firebrick")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(5,"pt"),)
