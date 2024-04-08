################################################################################
# R script to characterize the networks in the two niches on the basis od
# global networks properties
################################################################################


# LOAD NETWORK RESULTS FROM PREVIOUS STEPS #####################################


#load the results ----

folds <- c("spieceasi_2021-01-20/glasso/")

data <-  load_network_results(paste0(dir_build_networks, folds), sple = FALSE)

data_spl <- load_network_results(paste0(dir_build_networks, folds), sple = TRUE)


# Define the save folder ----

today <- "2021-01-26"# Sys.Date()
dir_save <- paste0(dir_res_global, today, "/")
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


# colors ----

col_domain <- c("#CD2626", "#00B2EE", "#EEAD0E")
names(col_domain) <- c("Archaea", "Bacteria", "Eukaryota")




#### COMPARE METRICS AT THE NETWORK LEVEL #########################################


# sumarize the metrics =========================================================


# for aggregated networks  ----

data$graph_metrics


# for sample level networks ----

metrics_sple <- data_spl$graph_metrics 

write.csv(metrics_sple, file = paste0(dir_save, 
                                      "table_network_metrics_per_sample.csv"))
