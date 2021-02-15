#!/usr/bin/env Rscript
suppressMessages(library(optparse))
suppressMessages(library(lpbrim))
suppressMessages(library(vegan))

option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file root name", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

links <- read.csv(opt$file, header=T, row.names=1)

# nestedness
nodf <- nestednodf(links)


# modularity
B <- as.matrix(links)
B <- B[rowSums(B)>0, colSums(B)>0]

Q <- ''
if (!is.null(dim(B)) ){
  Q <- findModules(B, sparse=FALSE)$Q

}

cat(nodf$statistic['NODF'],Q)
