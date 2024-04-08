################################################################################
# R script to characterize the compositnio of the modules
################################################################################


# LOAD NETWORK RESULTS FROM PREVIOUS STEPS #####################################


#load the results ----

folds <- c("spieceasi_2021-01-20/glasso/")

data <-  load_network_results(paste0(dir_build_networks, folds), sple = FALSE)

data_spl <- load_network_results(paste0(dir_build_networks, folds), sple = TRUE)


# Define the save folder ----

today <- "2021-01-26"# Sys.Date()
dir_save <- paste0(dir_res_mods, today, "/")
dir.create(dir_save)

# categorical metadata ----

metadata <- read.csv(paste0(dir_data, "table_metadata.csv"), 
                     row.names = 1)
row.names(metadata) <- metadata$sample


# phyloseq object ----

dir_load <- paste0(dir_build_networks, "spieceasi_2021-01-20/")
ls_ps_compo <- readRDS(paste0(dir_load, "list_phyloseq_compositional.rds"))

taxo_table <- do.call(rbind, lapply(ls_ps_compo, function(X) { 
  X@tax_table@.Data %>% data.frame()}))
row.names(taxo_table) <- taxo_table$OTU

ps_merged <- merge_phyloseq(ls_ps_compo$Bacteria, ls_ps_compo$Archaea, ls_ps_compo$Eukaryota)


# colors ----

col_domain <- c("#CD2626", "#00B2EE", "#EEAD0E")
names(col_domain) <- c("Archaea", "Bacteria", "Eukaryota")


# the global taxonomy table
taxo_table <- ps_merged@tax_table@.Data %>% data.frame()

# the global OTU table
otu_table <- ps_merged@otu_table@.Data %>% data.frame()


#### COMPARE MODULES  ##########################################################


# load the results ----

folds <- c("spieceasi_2021-01-20/glasso/")

mods <- readRDS(paste0(dir_build_networks, folds, "results_detected_clusters.rds"))
mods_spl <- readRDS(paste0(dir_build_networks, folds, "results_detected_clusters_sple.rds"))

# Compare the number of modules  ----
# We will analyse the composition of the modules identified on the global networks 
#  from each niches there is no need to go at the sample level here.
# However, we'll use the sample level information to determine which algorithm 
# describes the best the differences in terms of modularity between networks.

# number of modules  for each algorithm

tmp <- do.call(cbind, lapply(mods, function(X) do.call(c, lapply(X, length))))

# average number of module accross algorithms 

apply(tmp, 2, mean)
apply(tmp, 2, sd)




# ======================  SIZE OF THE MODULES =====================

# Set of functions to characterize the clusters object
#  that contains "communities" type objects ------------------------------------
# We'll use the walktrap algorithm as it return numbers of modules that are close
# to the average

nm_clust_algo <- "walktrap"

# groups(clust)  # return a list containing the nodes ids for each module
# sizes(clust)   # return a table with the numbre of node in each module
# length(clust)  # return the number of detected modules
# membership(clust) # return a vector of node membership
# modularity(clust) # return the modularity score of the network partitioning

# Get the taxonomy of modules ----

taxonomy_modules <- data.frame(taxo_table, 
                               module = paste0("mod_", 
                                               as.vector(membership(mods[[nm_clust_algo]]))))



# Get the size of the modules  ----

mods_size <- taxonomy_modules %>% group_by(module) %>% dplyr::tally()

mods_size %>% dplyr::summarise(avg = mean(n), sd = sd(n))



# Get the composition for each rank

nms_ranks <- colnames(taxo_table)[1:6]

mods_taxonomy <- lapply(nms_ranks, function(rk) {
  out <- do.call(cbind, lapply(split(taxonomy_modules, taxonomy_modules$module), function(x) {
    table(x[, rk])
  }))
  out[, order(colSums(out))]
}) %>% setNames(nms_ranks)

mods_compo <- lapply(mods_taxonomy, function(x) {
  apply(x, 2, function(xx) xx / sum(xx))
}) 



# ================== TAXONOMIC COMPOSITION OF MODULES ==========================

# create some usefull objects

nms_ranks <- colnames(taxo_table)[1:6]
nms_modules <- mods_size$module

# get the taxonomy table for each module ----

taxonomy_table_per_mods <- split(taxonomy_modules, taxonomy_modules$module)

# list the otus from each module ----

nms_otu_per_mods <- lapply(taxonomy_table_per_mods, 
                              function(x) { as.character(x$OTU) })


# number of OTU per taxaXmodule at each taxonomic rank ----

num_otu_per_mods <- lapply(nms_ranks, function(rk) {
  do.call(cbind, lapply(taxonomy_table_per_mods, function(X) {
     table(X[, rk])
  }))
}) %>% setNames(nms_ranks)

# proportion of OTU in each taxa per module at each taxonomic rank ----

composition_taxa_per_mods <- lapply(num_otu_per_mods, function(X) {
    apply(X, 2, function(xx) xx / sum(xx))
})
