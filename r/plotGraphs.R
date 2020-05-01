## Ford Fishman

library("ggplot2")
library("tidyr")
library("optparse")

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
path1 <- paste(opt$file, "_main",sep="")
path2 <- paste(opt$file, "_richness",sep="")

dfile1 <- paste(path1, "csv", sep=".")
dfile2 <- paste(path2, "csv", sep=".")

df1 <- read.csv(file=dfile1, header = T, row.names = 1)
df2 <- read.csv(file=dfile2, header = T, row.names = 1)

# df1$vulnerable <- df1$host - df1$resistant

df1 <- gather(df1, key = organism, value=density, -time)
df2 <- gather(df2, key = organism, value=richness, HostRichness:PhageRichness)

# plots
p1 <- ggplot(df1, aes(x=time,y=density))+
        geom_line(aes(color=organism)) +
        scale_y_log10() +
        scale_linetype_discrete() +
        scale_color_manual("", values = c("blue","red","green","black")) +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black"),
              axis.text = element_text(size = 11),
              axis.ticks = element_line(size = 1),
              axis.ticks.length = unit(5,"pt"))

output1 <- paste(path1, "png", sep = ".")
ggsave(plot = p1, filename = output1, width = 7, height = 5)

p2 <- ggplot(df2, aes(x=time,y=richness))+
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

output2 <- paste(path2,"png",sep=".")
ggsave(plot = p2, filename = output2, width = 7, height = 5)


