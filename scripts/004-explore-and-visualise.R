library(tidyverse)
library(ggplot2)
library(ggtree)

# inspect MCMC trace

library(RevGadgets)
mcmc_trace <- readTrace("output/mcmc_log.txt")

ggplot(mcmc_trace[[1]]) +
  aes(Iteration, 
      sigma) +
  geom_line()

plotTrace(mcmc_trace, vars = "sigma")

# explore some plotting methods

outsumfile <- ape::read.nexus("output/map_tree.nex")

# this is only possible if we are looking at the FR subset, the other 
# regions don't have so much chrono data
tip_labels <- 
  names_artefacts_ID_and_period %>% 
  mutate(artefact_ID = str_remove(artefact_ID, "_pseudo_no")) %>% 
  mutate(artefact_ID = str_remove_all(artefact_ID, "_XX")) %>% 
  unite('country_and_period',  
        ID_country, Period, 
        sep = "_",
        remove = FALSE) %>% 
  # taxa names must be first col for ggtree to work with %<+%
  mutate(taxa = row.names(data_for_revbayes) ) %>% 
  relocate(taxa) %>% 
  as_tibble

# if working with all countries
# tip_labels <- data.frame(tip.label = outsumfile$tip.label,
#                         ID_country = str_extract(outsumfile$tip.label, "[A-Z]{2}"))


ggtree(outsumfile) %<+% tip_labels +  # names_artefacts_ID_and_period_unite + 
  theme_tree() + 
  geom_tiplab(size=4, 
              aes(label = ID_country,
                  color = Period)) +
  geom_treescale()  +
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 120, 6, 6)) +
  theme(legend.position = c(0.2, 0.9))

ggsave("figures/map_tree_fr_period.png",
       h = 15, w = 7)


# https://cran.r-project.org/web/packages/MCMCtreeR/vignettes/MCMCtree_plot.html
library(MCMCtreeR, quietly = TRUE, warn.conflicts = FALSE)

outsumfile.rb <-"output/map_tree.nex"

MCMC.tree.plot(analysis.type='revbayes',
               directory.files=outsumfile.rb, 
               cex.tips=0.33,
               plot.type='phylogram', 
               lwd.bar=2,
               add.time.scale=FALSE,
               node.method='bar', 
               col.age='navy')

MCMC.tree.plot(
  directory.files=outsumfile.rb,
  analysis.type = "revbayes",
  cex.tips = 0.2,
  plot.type = "phylogram",
  lwd.bar = 2,
  node.method = "node.length",
  cex.labels = 0.5,
  add.time.scale=FALSE)

# read in the tree space
phy <- ape::read.tree( "output/tree_trace.trees" )

# write it out to a format that crsl4/PhyloNetworks.jl can read
# cf. https://crsl4.github.io/PhyloNetworks.jl/latest/man/inputdata/
# ape::write.tree(phy, "output/phy.trees")
# ape::write.nexus(phy, "output/phy.nexus")

# It takes a moment to show
# phangorn::densiTree(phy, type = "phylogram", cex = 0.5)

# ggtree densitree method
ggdensitree(phy, 
            alpha= .01, 
            colour='steelblue') %<+% tip_labels + 
  theme_tree() + 
  geom_tiplab(size=2, 
              aes(label = ID_country,
                  color = Period)) +
  geom_treescale() 

ggsave("figures/densitree_tree_fr_period.png",
       h = 11, w = 7)


plot(phangorn::maxCladeCred(phy), cex = 0.5)
# plot(consensus(phy))

# Looks promising:
# https://revbayes.github.io/tutorials/intro/revgadgets.html
library(RevGadgets)

outtree.rb <- readTrees(outsumfile.rb)
ape::write.tree(outtree.rb[[1]][[1]]@phylo, "output/out.trees")

outtree.rb1 <- readTrees("output/out.trees")
outtree.rb2 <- readTrees("output/tree_trace.trees")

plot(outtree.rb1[[1]][[1]]@phylo)

# create the plot of the rooted tree
plot <- plotTree(tree = outtree.rb1,
                 # label nodes the with posterior probabilities
                 # node_labels = "posterior", 
                 # node_pp = TRUE,
                 # offset the node labels from the nodes
                 node_labels_offset = 0.005,
                 # make tree lines more narrow
                 line_width = 0.5,
                 # italicize tip labels 
                 tip_labels_italics = FALSE)

# add scale bar to the tree and plot with ggtree
library(ggtree)
plot + geom_treescale(x = -0.35, y = -1)

# https://thibautjombart.github.io/treespace/articles/introduction.html
library("treespace")
library("adegenet")
library("adegraphics")
library("rgl")
library(cowplot)

# use treespace
res <- treespace::treespace(phy, nf=3)

png("figures/treespace_tree_fr_period.png")
plotGroves(res$pco, lab.show=F, lab.cex=1.5)
dev.off()

png("figures/treespace_cluster_fr_period.png")
wm.groves <- findGroves(res, nclust = 3)
plotGroves(wm.groves)
dev.off()

plotGrovesD3(wm.groves)

# https://thibautjombart.github.io/treespace/articles/introduction.html
# get first median tree
median_trees <- medTree(phy, wm.groves$groups)
med.trees <- lapply(median_trees, function(e) ladderize(e$trees[[1]]))

## plot trees
par(mfrow=c(2,3))
for(i in 1:length(med.trees)) plot(med.trees[[i]], main=paste("cluster",i))

## highlight the differences between a pair of median trees
tr1 <- med.trees[[1]]
tr2 <- med.trees[[2]]
tr3 <- med.trees[[3]]

# find the tip differences in advance, to avoid recalculating with each plot
wmTipDiff <- tipDiff(tr1,tr2, sizeOfDifferences=TRUE)

plotTreeDiff(tr1,
             tr2,
             wmTipDiff,
             tipMatch=TRUE,
             treesFacing=TRUE,
             baseCol="black") # the colour used for tips with identical ancestry in the two trees.


ggtree(tr1) %<+% tip_labels +  # names_artefacts_ID_and_period_unite + 
  theme_tree() + 
  geom_tiplab(size=2, 
              aes(label = ID_country,
                  color = Period)) +
  geom_treescale()  +
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 120, 6, 6)) +
  theme(legend.position = c(0.2, 0.9))

ggsave("figures/median_1_tree_fr_period.png",
       h = 7, w = 7)

ggtree(tr2) %<+% tip_labels +  # names_artefacts_ID_and_period_unite + 
  theme_tree() + 
  geom_tiplab(size=2, 
              aes(label = ID_country,
                  color = Period)) +
  geom_treescale()  +
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 120, 6, 6)) +
  theme(legend.position = c(0.2, 0.9))

ggsave("figures/median_2_tree_fr_period.png",
       h = 7, w = 7)

ggtree(tr3) %<+% tip_labels +  # names_artefacts_ID_and_period_unite + 
  theme_tree() + 
  geom_tiplab(size=2, 
              aes(label = ID_country,
                  color = Period)) +
  geom_treescale()  +
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 120, 6, 6)) +
  theme(legend.position = c(0.2, 0.9))

ggsave("figures/median_3_tree_fr_period.png",
       h = 7, w = 7)


