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
path1 <- paste(opt$file, "main",sep="")
path2 <- paste(opt$file, "richness",sep="")

dfile1 <- paste(path1, "csv", sep=".")
dfile2 <- paste(path2, "csv", sep=".")

df1 <- read.csv(file=dfile1, header = T, row.names = 1)
df2 <- read.csv(file=dfile2, header = T, row.names = 1)

# df1$vulnerable <- df1$host - df1$resistant

# df1 <- gather(df1, key = organism, value=density, -time)
# df2 <- gather(df2, key = organism, value=richness, HostRichness:PhageRichness)

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

# plots
p1 <- ggplot(totalAbundance2, aes(x=Timestep,y=abundance, group = type)) +
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

output1 <- paste0(opt$file, "total.png")
ggsave(plot = p1, filename = output1, width = 7, height = 5)

p2 <- ggplot(df1, aes(x=timestep,y=log10(pop), group = name))+
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

output2 <- paste0(opt$file,"strains.png")
ggsave(plot = p2, filename = output2, width = 7, height = 5)


