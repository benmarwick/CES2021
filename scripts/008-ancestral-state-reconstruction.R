# from http://www.phytools.org/eqg2015/asr.html
library(phytools)

## read tree from file
anc.tree <- ape::read.nexus(here("output/map_tree.nex"))

## read data
anc.data <- read.csv(here("data/data_for_revbayes.csv"), 
                     row.names=1)
## change this into a vector
anc.data.1    <- as.matrix(anc.data)[,1]
anc.data.1.2  <- as.matrix(anc.data)[,1:2]
anc.data.1.4  <- as.matrix(anc.data)[,1:4]

# Now we can estimate ancestral states. We will also compute variances & 95% confidence intervals for each node:
fit <- fastAnc(anc.tree, 
               anc.data.1, # focus on PC1 shape dimension 
               vars=TRUE, 
               CI=TRUE)

## projection of the reconstruction onto the edges of the tree
obj<-contMap(anc.tree, 
             anc.data.1, 
             plot=FALSE, 
             lwd = 0.5,
             fsize = c(0.1, 0.1))
n <- length(obj$cols)
obj$cols[1:n] <- viridis::viridis(n)

png(here("figures/anc_recon_example.png"),
    height = 15,
    width = 10,
    units = "in",
    res = 400)
plot(obj)
dev.off()

ggsave(here("figures/anc_recon_example_single.png"),
       h = 10, w = 8)


## projection of the tree into phenotype space, take a long time to draw
phenogram(anc.tree, anc.data.1, fsize=0.6, spread.costs=c(1,0))

phylomorphospace(anc.tree,
                 anc.data.1.2[,c(1,2)],
                 xlab="trait 1",
                 ylab="trait 2")

# continuous character mapping and 2D phylomorphospaces using a 
# type of phylogenetic scatterplot. 

x <- 
fancyTree(anc.tree,
          type="scattergram",
          X=anc.data.1.4,
          fsize = 0.1)

# update the colour scheme
for(i in 1:length(x$contMaps)){
n <- length(x$contMaps[[i]]$cols)
x$contMaps[[i]]$cols[1:n] <- viridis::viridis(n)
}

png(here("figures/anc_recon_example_multiple.png"),
    height = 10,
    width = 10,
    units = "in",
    res = 400)
plot(x, fsize = 0.1)
dev.off()
