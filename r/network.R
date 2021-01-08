library(igraph)
library(bipartite)
library(lpbrim)
library(vegan)
library(RInSp)
library(here)

path <- 'C:/Users/fordf/OneDrive/Documents/GitHub/CommunityImmunity/output/1sims_t1000_a1e6_adsp1e-8_b1.2_beta100_c0.01_d0.1_f0_l0.9_m1e-7_pS1e-6_pop1e5_phage1e7/adjacency.csv'

links <- read.csv(path, header=T, row.names=1)
net <- graph_from_incidence_matrix(links)
# net.bp <- bipartite.projection(net)

V(net)$color <- c("steel blue", "orange")[V(net)$type+1]
V(net)$shape <- c("square", "circle")[V(net)$type+1]

plot(net, vertex.label=NA, vertex.size=7, layout=layout_as_bipartite)

# nestedness
nodf <- nestednodf(links)
nodf$statistic['NODF']

# modularity
B <- as.matrix(links)
B <- B[rowSums(B)>0, colSums(B)>0]
findModules(B, sparse=FALSE)





