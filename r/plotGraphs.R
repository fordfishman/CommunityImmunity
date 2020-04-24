## Ford Fishman

library("ggplot2")
library("tidyr")
library("optparse")

option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.png", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# the datafile
df <- read.csv(file=opt$file, header = T, row.names = 1)

df1 <- gather(df, key = organism, value=density, host:phage)

p1 <- ggplot(df1, aes(x=time,y=density))+
        geom_line(aes(color=organism)) +
        scale_y_log10() +
        scale_color_manual("", values = c("blue","red")) +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black"),
              axis.text = element_text(size = 11),
              axis.ticks = element_line(size = 1),
              axis.ticks.length = unit(5,"pt"))
# fig3_dir <- paste(figure_dir, "Pres12_6_LambdaVSEpsilon_ConstantRatesHeatmap.png", sep = "")
ggsave(plot = p1, filename = opt$out, width = 7, height = 5)
