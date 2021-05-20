library(RevGadgets)

#  visualize the branch-specific rates by plotting them using 
# R package RevGadgets.
  
# read in the tree annotated with the branch rates:
log_file <- here("output/relaxed_multivariate_BM_run_1.log")
tree <- readTrees(here("output/relaxed_multivariate_BM_MAP.tre"))


RevGadgets::plotTree(tree, 
         branch_color = viridis::viridis(n = 2, 
                                         begin = 0 , 
                                         end = 1, 
                                         direction = -1),
         color_branch_by="branch_rates") +
  theme(legend.position = c(0.2, 0.8))

ggsave(here("figures/branch_rates_tree_fr_period_mvBM.png"),
       h = 10, w = 10)

phy_mvBM <- ape::read.tree( here("output/relaxed_multivariate_BM.trees" ))


# phangorn::densiTree(phy_mvBM, type = "phylogram", cex = 0.5)


# ggtree densitree method
ggdensitree(phy_mvBM, 
            alpha= .01, 
            colour='steelblue') %<+% tip_labels + 
  theme_tree() + 
  geom_tiplab(size=2, 
              aes(label = ID_country,
                  color = Period)) +
  geom_treescale() 

ggsave(here("figures/densitree_tree_fr_period_mvBM.png"),
       h = 15, w = 7)


##---------------------------------------

# read the MCMC output:

# import trace with our modified function
samples <- readTrace(log_file)
# plot the posterior distribution of the first n correlation 
# parameters rate parameters (the top row of the correlation matrix, 
# which represents how characters two through n are correlated with the first var):
  
RevGadgets::plotTrace(samples, vars=paste0("correlations[",1:5,"]"))



         

