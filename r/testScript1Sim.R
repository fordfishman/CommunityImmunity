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
filename1 <- paste0(wd, "/output/1sims_t5000_a1e6_adsp1e-8_b1.2_beta100_c0.01_d0.1_f0_l0.9_m1e-7_pS1e-6_pop1e5_phage1e7/main.csv")
filename2 <- paste0(wd, "/output/1sims_t5000_a1e6_adsp1e-8_b1.2_beta100_c0.01_d0.1_f0_l0.9_m1e-7_pS1e-6_pop1e5_phage1e7/richness.csv")
filename3 <- paste0(wd, "/output/1sims_t5000_a1e6_adsp1e-8_b1.2_beta100_c0.01_d0.1_f0_l0.9_m1e-7_pS1e-6_pop1e5_phage1e7/full.csv")

df1 <- read.csv(filename1, header = T, row.names = 1)
df2 <- read.csv(filename2, header = T, row.names = 1)
df3 <- read.csv(filename3, header = T, row.names = 1)

df1.1 <- gather(df1, key = type, value = density, N:I)
df1.1$density <- log10(df1.1$density)

## plot of the sums across groups
ggplot(df1.1, aes(x=t,y=density, group = type)) +
  geom_line(aes(color=type)) +
  scale_y_continuous("Density",labels = math_format(10^.x), expand = c(0,0)) +
  scale_x_continuous("Timestep")+
  scale_color_manual("",
                     breaks = c("N","S","I","P"),
                     values = c("black","blue","green","red"),
                     labels = c("Total Host","Initial Host","Spacer Variant","Phage")) +
  annotation_logticks(sides = "l") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        axis.ticks.length = unit(5,"pt"),)

# all strains: shows coevolution
ggplot(df3, aes(x=timestep,y=log10(pop), group = name))+
  geom_line(aes(color=type)) +
  # scale_y_continuous("Density", limits = c(1, 0.5+maxPop),labels =  math_format(10^.x), ) +
  scale_y_continuous("Density", limits = c(0, 8),breaks = seq(0,10,by = 2),labels = math_format(10^.x), expand = c(0,0)) +
  scale_x_continuous("Timestep",)+
  scale_color_manual("", labels = c("Initial Strain","Spacer Variant","Phage"),values = c("blue","green","red")) +
  annotation_logticks(sides = "l") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 11),
        axis.ticks.length = unit(5,"pt"),)

