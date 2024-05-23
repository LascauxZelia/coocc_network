################################################################################
# R script to globally set up the project
################################################################################


# ==============================================================================
# cleaning memory
# ==============================================================================

rm(list = ls()) 

# ==============================================================================
# Define the local directories
# create objects corresponding to the directories and create directories locally
# ==============================================================================

dir_project  <- getwd()
dir_scripts  <- paste0(dir_project, "/scripts/")
dir_functions <- paste0(dir_scripts, "/functions/")
dir_analyses <- paste0(dir_project, "/analyses/")
dir_data     <- paste0(dir_project, "/data/")

dir_build_networks   <- paste0(dir_analyses, "01_build_networks/")
dir_compare_networks <- paste0(dir_analyses, "02_compare_networks/")
dir_analyse_results  <- paste0(dir_analyses, "03_analyse_results/")
dir_res_global <- paste0(dir_analyse_results, "global_properties/")
dir_res_nodes  <- paste0(dir_analyse_results, "nodes_properties/")
dir_res_edges  <- paste0(dir_analyse_results, "edges_properties/")
dir_res_mods   <- paste0(dir_analyse_results, "modules_properties/")

dir.create(dir_analyses)
dir.create(dir_build_networks)
dir.create(dir_compare_networks)
dir.create(dir_analyse_results)
dir.create(dir_res_global)
dir.create(dir_res_nodes)
dir.create(dir_res_edges)
dir.create(dir_res_mods)


# ==============================================================================
# loading libraries
# ==============================================================================

# source utility functions
for (f in list.files(dir_functions, full.names = T)) { source(f) }

# source package list
source(paste0(dir_scripts, "packages.R"))
install_if_not_there(cran_packages, type = "CRAN")
install_if_not_there(bioc_packages, type = "bioconductor")
sapply(c(cran_packages, bioc_packages), require, character.only = TRUE)

# If you have an error, try to install manually SpiecEasi
library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)

# ==============================================================================
# Some objects for analyses
# ==============================================================================

nms_ranks <- c("Domain","Phylum","Class","Order","Family","Genus","Species", "OTU")

nms_idx_centrality <- c("degree", "betweenness", "closeness", "hub_score")

nms_clust_algo <- c("edge_betweenness", "fast_greedy", "infomap",
                    "leading_eigen", "louvain", "walktrap")
