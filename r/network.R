library(optparse)
library(lpbrim)
library(vegan)

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

path <- paste0(opt$file, "adjacency.csv")
log_file <- paste0(opt$file, "log_file.txt")
# path <- 'C:/Users/fordf/OneDrive/Documents/GitHub/CommunityImmunity/output/1sims_t5000_a1e6_adsp1e-8_b1.2_beta100_c0.01_d0.1_f0_l0.9_m1e-7_pS1e-6_pop1e5_phage1e7/adjacency.csv'

links <- read.csv(path, header=T, row.names=1)
# net <- graph_from_incidence_matrix(links)
# net.bp <- bipartite.projection(net)

# V(net)$color <- c("steel blue", "orange")[V(net)$type+1]
# V(net)$shape <- c("square", "circle")[V(net)$type+1]
# 
# plot(net, vertex.label=NA, vertex.size=7, layout=layout_as_bipartite)

# nestedness
nodf <- nestednodf(links)
write("NODF:",file=log_file, append=T)
write(nodf$statistic['NODF'],file=log_file, append=T)


# modularity
B <- as.matrix(links)
B <- B[rowSums(B)>0, colSums(B)>0]
Q <- findModules(B, sparse=FALSE)$Q
write("\nQ:",file=log_file, append=T)
write(Q,file=log_file, append=T)






