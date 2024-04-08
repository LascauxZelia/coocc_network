################################################################################
# R script to characterize the nodes characteristics
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



#### ANALYSE NODES ATTRIBUTES ##################################################

# Make a table with nodes attributes data ----

nodes <- data$vertices_attr
rownames(nodes) <- nodes$OTU

# > names(nodes)
# [1] "Domain"              "Phylum"              "Class"              
# [4] "Order"               "Family"              "Genus"              
# [7] "Species"             "OTU"                 "taxonomy"           
# [10] "otu"                 "degree"              "betweenness"        
# [13] "closeness"           "eigenvector"         "degree_std"         
# [16] "betweenness_std"     "closeness_std"       "eigenvector_std"    
# [19] "hub_score"           "authority_score"     "edge_betweenness"   
# [22] "fast_greedy"         "infomap"             "label_prop"         
# [25] "leading_eigen"       "louvain"             "walktrap"           
# [28] "pc_edge_betweenness" "pc_fast_greedy"      "pc_infomap"         
# [31] "pc_label_prop"       "pc_leading_eigen"    "pc_louvain"         
# [34] "pc_walktrap"         "z_edge_betweenness"  "z_fast_greedy"      
# [37] "z_infomap"           "z_label_prop"        "z_leading_eigen"    
# [40] "z_louvain"           "z_walktrap"                        


##  Analyze nodes centrality attributes ========================================

# plot the node centrality metrics ----

png(paste0(dir_save, "plot_nodes_centrality_metrics.png"), height = 15, 
    width = 20, unit = "cm", res = 200)

  par(mfrow = c(2,2), mar = c(4,4,1,1), mgp = c(2.75,0.5,0), oma = c(1,1,1,1), las = 1)
  lapply(nms_idx_centrality, function(idx) {
    d <- density(nodes[,idx])
    plot(d, type = "n", main ="", xlab = idx)
    polygon(d)
  })
dev.off()


nodes %>% summarize_at(nms_idx_centrality, funs(mean, sd)) %>% data.frame()



##  Identify OTU with the most extremes attributes values ======================

nms_ranks <- c("Domain","Phylum","Class","Order","Family","Genus","Species", "OTU")

# Get a table for each niche with the OTU names orderd for each index ----

ranked_nodes <- do.call(cbind, lapply(nms_idx_centrality, function(idx) {
    nodes[order(nodes[,idx]), "OTU"]
  }) %>% setNames(nms_idx_centrality))


# Get the rank of OTU for each index, estimate the average and SD ----
# and add the taxonomic information

out <- do.call(cbind, lapply(nms_idx_centrality, function(idx) {
  rank(nodes[,idx], ties.method = "average")
}) %>% setNames(nms_idx_centrality)) %>% data.frame()
out$avg_rank <- apply(out, 1, mean)
out$sd_rank <- apply(out, 1, sd)
nodes_ranks <- out %>% mutate(OTU = nodes$OTU) %>% 
  inner_join(taxo_table, "OTU")# %>% arrange(avg_rank)


df_nodes_rank <-  nodes_ranks %>% 
  arrange(avg_rank) %>% 
  rename_at(nms_idx_centrality, funs(paste0("rk_", nms_idx_centrality))) %>% 
  dplyr::select(-nms_ranks[-8]) %>% 
  left_join(nodes, "OTU") %>% 
  dplyr::select(nms_ranks, nms_idx_centrality, everything())


write.csv(df_nodes_rank, file = paste0(dir_save, "table_all_nodes_with_centrality_ranks.csv"),
          row.names = FALSE)
