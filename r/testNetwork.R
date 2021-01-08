library(igraph)
library(bipartiteD3)
library(lpbrim)
library(vegan)
library(RInSp)
library(here)
library(tidyr)

path <- 'C:/Users/fordf/OneDrive/Documents/GitHub/CommunityImmunity/output/1sims_t5000_a1e6_adsp1e-8_b1.2_beta100_c0.01_d0.1_f0_l0.9_m1e-7_pS1e-6_pop1e5_phage1e7/adjacency.csv'
# path <- '/mnt/c/Users/fordf/OneDrive/Documents/GitHub/CommunityImmunity/output/1sims_t5000_a1e6_adsp1e-8_b1.2_beta100_c0.01_d0.1_f0_l0.9_m1e-7_pS1e-6_pop1e5_phage1e7/adjacency.csv'
links <- read.csv(path, header=T, row.names=1)
net <- graph_from_incidence_matrix(links)

# nestedness
nodf <- nestednodf(links)
nodf$statistic['NODF']

# modularity
B <- as.matrix(links)
A <- B[rowSums(B)>0, colSums(B)>0]
mod <- findModules(A, sparse=FALSE)


# plots
bipartite_D3(B, PrimaryLab = 'Host Strains', SecondaryLab = 'Phages', SiteNames = 'Infection')
plotMatrixModules(mod, 'frames')

heatmap(B, Colv = NA, Rowv = NA, scale="column")
visweb(A, type = 'nested',xlabel='Phage', ylabel='Host Strain')
mtext("Another \nlabel", side=4, line=1, at=-1, col="green2", las=1)# as_adj_edge_list(net2)
# m <- as_edgelist(net2)
# m <- cbind(m, rep(1, nrow(m)))
# A <- netToAdjacency(m)
# 
# M <- matrix(rbinom(100, 1, 0.3), ncol=10)
# M <- M[rowSums(M)>0, colSums(M)>0]
# mod <- findModules(M, iter=2, sparse=FALSE)