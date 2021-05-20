library(tidyverse)

# attach country and chrono data, only possible for FR data

data_for_revbayes <- nicolas_2016_without_outliers_PCA_as_df_subset_typochron_FR %>% 
  # only want columns that have PC values
  dplyr::select(starts_with("PC"))

# shorten the sample names so the trees are easier to read
row.names(data_for_revbayes) <- 
  str_remove(row.names(data_for_revbayes), 
             "_pseudo_no")

row.names(data_for_revbayes) <- 
  str_remove_all(row.names(data_for_revbayes), 
             "_XX")

# write a nex file for input into RevBayes, just get a few PCs?
write.nexus.data(as.matrix(data_for_revbayes),
                 file.path(here("data/data_for_revbayes.nex")),
                 format = "continuous")

# save this for ASR
data_for_revbayes %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  write_csv(here("data/data_for_revbayes.csv"))

# run RevBayes file to set up and run model, output is exported as files 
# to data directory

# Start a command with results displayed in a terminal buffer
# do tree inference
termId1 <- rstudioapi::terminalExecute('/Applications/rb "003-RevBayes.Rev"')

# Start a command with results displayed in a terminal buffer
# do Marginal Likelihood estimation so we can do Bayes Factors later
termId2 <- rstudioapi::terminalExecute('/Applications/rb "003.1-RevBayes.Rev"')

# quit the terminal
rstudioapi::terminalKill(termId1)
rstudioapi::terminalKill(termId2)

#  go to file 004