## Ford Fishman

library("ggplot2")
library("tidyr")
library("optparse")
library("scales")
library("MASS")

option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file root name", metavar="character")
  # make_option(c("-o", "--out"), type="character", default="out.png", 
  #             help="output file name [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# the datafile
dfile1 <- paste0(opt$file, "main.csv")
dfile2 <- paste0(opt$file, "richness.csv")
dfile3 <- paste0(opt$file, "full.csv")

# dfile1 <- paste(path1, "csv", sep=".")
# dfile2 <- paste(path2, "csv", sep=".")
# dfile3 <- paste(path3, "csv", sep=".")

df1 <- read.csv(file=dfile1, header = T, row.names = 1)
df2 <- read.csv(file=dfile2, header = T, row.names = 1)
df3 <- read.csv(file=dfile3, header = T, row.names = 1)

df1.1 <- gather(df1, key = type, value = density, N:I)
df1.1$density <- log10(df1.1$density)

## plots
## plot of the sums across groups
p1 <- ggplot(df1.1, aes(x=t,y=density, group = type)) +
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

output1 <- paste0(opt$file, "total.png")
ggsave(plot = p1, filename = output1, width = 7, height = 5)


# all strains: shows coevolution
p2 <- ggplot(df3, aes(x=timestep,y=log10(pop), group = name))+
  geom_line(aes(color=type)) +
  scale_y_continuous("Density", limits = c(0, 8),breaks = seq(0,10,by = 2),labels = math_format(10^.x), expand = c(0,0)) +
  scale_x_continuous("Timestep",)+
  scale_color_manual("", labels = c("Initial Strain","Spacer Variant","Phage"),values = c("blue","green","red")) +
  annotation_logticks(sides = "l") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(5,"pt"),)

output2 <- paste0(opt$file,"strains.png")
ggsave(plot = p2, filename = output2, width = 7, height = 5)


