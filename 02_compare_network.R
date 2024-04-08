################################################################################
# R script to visualize the networks
################################################################################


# LOAD NETWORK RESULTS FROM PREVIOUS STEPS #####################################


# here for analysis of a singl set of eults

folds <- c("spieceasi_2021-01-20/glasso/")

data <-  load_network_results(paste0(dir_build_networks, folds), sple = FALSE)

data_spl <- load_network_results(paste0(dir_build_networks, folds), sple = TRUE)

# Define the save folder ----

dir_save <- paste0(dir_compare_networks, "res_2021-04-01/")
dir.create(dir_save)

# categorical metadata ----

metadata <- read.csv(paste0(dir_data, "table_metadata.csv"), 
                     row.names = 1)
row.names(metadata) <- metadata$sample


#### COMPARE METRICS AT THE NETWORK LEVEL #########################################


metrics_sple <- data_spl$graph_metrics 

# plot of global network metrics

nms_idx <- c("num_edges", "diameter", "density", "connectance", "avg_distance",
             "centralized_degree", "centralized_betweenness", 
             "hub_score_graph")

png(paste0(dir_save, "plot_network_level_metrics.png"), height = 25, 
    width = 25, unit = "cm", res = 600)

par(mfrow = c(4,2), las = 1, mar = c(3,4,1,1), oma = c(1,1,1,1))
par(las = 1, mar = c(3,10,1,1), oma = c(1,1,1,1))
lapply(nms_idx, function(x) {
  boxplot(metrics_sple[,x] ,
          main = x, xlab = "", ylab = "", horizontal = TRUE) 
})

dev.off()



#### VISUALIZATION OF THE NETWORKS ###############################################

ntwks <- data_spl$networks

coords   <- lapply(ntwks, function(X) layout.fruchterman.reingold(X) )
vec_cols <- lapply(ntwks, function(X) {
  vertex_attr(X)[["Domain"]] %>% as.factor() %>% as.numeric()
})

png(paste0(dir_save, "plot_networks_per_sample.png"), height = 25, 
    width = 25, unit = "cm", res = 600)

par(mfrow = c(5,6), mar = c(1,1,1,1))

lapply(names(ntwks), function(X) {
  plot(ntwks[[X]], layout = coords[[X]], vertex.label = NA,
       main = X, vertex.size = 3, vertex.color = vec_cols[[X]])
})
dev.off()
