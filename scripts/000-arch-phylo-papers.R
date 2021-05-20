

library(tidyverse)
library(lubridate)

arch_phylo_papers <- 
  readxl::read_excel(here::here("CES2021/data/Archaeology phylogenetics case studies.xlsx"))

# papers per year
ggplot(arch_phylo_papers) +
  aes(`Publication Year`)  +
  geom_bar() +
  theme_minimal(base_size = 14) +
  labs(x = "")

ggsave(here::here("CES2021/figures/arch-phylo-papers-per-year-plot.png"),
       w = 10,
       h = 3)

arch_phylo_papers_rate <- 
arch_phylo_papers %>%
  mutate(`time period` = (year(floor_date(strptime(`Publication Year`, "%Y"), "2 years")))) %>%
  group_by(`time period`) %>% 
  tally() %>% 
  mutate(`time period` = paste0((`time period`), " - ", lead(`time period`-1) ))

ggplot(arch_phylo_papers_rate) +
  aes(`time period`, n)  +
  geom_col()



# basic method details
arch_phylo_papers_long <- 
arch_phylo_papers %>% 
  select(
         `Inference method`,
         Characters,
         Software ) %>% 
  pivot_longer(c( `Inference method`,
                  Characters,
                  Software))
  
ggplot(arch_phylo_papers_long) +
  aes(value) +
  geom_bar() +
  facet_wrap( ~ name, 
              ncol = 1,
              scales = "free") +
  theme_minimal(base_size = 14) +
  labs(x = "")

ggsave(here::here("CES2021/figures/arch-phylo-papers-panel-plot.png"),
       w = 5,
       h = 10)

# author network

first_names <- c("Michael", "John", "Ree" , "Peter", 
                 "Stephen", "Jelmer" ,  "Robert", "Richard" ,
                 "Briggs" , "Mark" , "Felix", "Ethan" ,  "Carl",
                 "Sean" , "Luke"  , "Jamie", "Fiona" ,  "Charles" ,
                 "Annaarie" , "Randall" , "Niles" , "Colin"   , 
                 "Matthew" , "Ben", "Christopher", "John",  "Thomas"  ,
                 "Danielarcía" , "James" , "Kristen" ,  "Marcelo" , 
                 "Judith"  , "Matt" , "Allison", "Christopher", 
                 "Eric"   , "Erik" , "Ashley" , 
                 "Charlotte", "Sébastien", "David", "JeanLoïce",
                 "Jimena", "Metin")
arch_phylo_papers_coauthors <- 
arch_phylo_papers %>% 
  select(Author) %>% 
  mutate(coauthors = map(Author, ~strsplit(.x, ",|;") %>% unlist)) %>% 
  mutate(coauthors = map(coauthors, str_squish)) %>% 
  mutate(coauthors = map(coauthors, ~str_remove_all(.x, "[[:punct:]]| [[:alpha:]]{1}"))) %>% 
  mutate(coauthors = map(coauthors, ~.x[!.x %in% first_names])) %>% 
  mutate(coauthors = map(coauthors, ~.x[!nchar(.x) == 1]))

coauthors <- arch_phylo_papers_coauthors$coauthors

# how many papers per author and gender of each author
author_details <- 
  table(unlist(coauthors)) %>% 
  enframe %>% 
  mutate(gender = case_when(
    name == "Tripp" ~ "female",
    name == "Prentiss" ~ "female",
    name == "Smallwood" ~ "female",
    name == "JordanFiona" ~ "female",
    name == "Charlin" ~ "female",
    name == "Pevny" ~ "female",
    TRUE ~ "male"
    
  ))


# transform into network structure
coauthors.unique <- unique(unlist(coauthors))[order(unique(unlist(coauthors)))]
bipartite.edges = lapply(coauthors, function(x) {coauthors.unique %in% x})
bipartite.edges = do.call("cbind", bipartite.edges) # dimension is number of authors x number of papers
rownames(bipartite.edges) = coauthors.unique
mat = bipartite.edges %*% t(bipartite.edges) #bipartite to unimode
mat = mat[order(rownames(mat)), order(rownames(mat))]

library(igraph)
library(statnet)
library(ggnetwork)
library(visNetwork)

statnet = as.network(mat, directed = FALSE, names.eval = "edge.lwd", ignore.eval = FALSE)

set.vertex.attribute(statnet, "papers", list(author_details$value))
set.vertex.attribute(statnet, "gender", list(author_details$gender))

ggplot(statnet, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "grey50") +
  geom_nodetext_repel(aes(label = vertex.names),
                      force = 20,
                      size = 8,
                      segment.colour = "white") +
  geom_nodes(aes(size = papers,
                 colour = gender)) +
  theme_blank()

ggsave(here::here("CES2021/figures/arch-phylo-papers-network.png"),
       w = 10,
       h = 10)
  
