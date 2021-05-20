# run RevBayes file to set up and run model, output is exported as files 
# to data directory

# Start a command with results displayed in a terminal buffer
termId1 <- rstudioapi::terminalExecute('/Applications/rb "006-RevBayes-mvBM.Rev"')

termId2 <- rstudioapi::terminalExecute('/Applications/rb "006.1-RevBayes-mvBM.Rev"')

# quit the terminal
rstudioapi::terminalKill(termId1)
rstudioapi::terminalKill(termId2)

#  go to file 007