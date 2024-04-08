################################################################################
# R script to characterize the networks edges
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



#### ANALYSE EDGES ATTRIBUTES ##################################################


##  Make a table with edges attributes ----

edges <- data$edges_attr

edges_spl <- data_spl$edges_attr


## plot edges weight ----

summ_edg_global <- edges %>% 
  dplyr::summarize(num_links = n(),
                   prop_positive = sum(edge_binary == 1) / num_links,
                   avg_weight = mean(edge_weight))

summ_edg_spl <- edges_spl %>% group_by(sample_name) %>% 
  dplyr::summarize(num_links = n(),
                   prop_positive = sum(edge_binary == 1) / num_links,
                   avg_weight = mean(abs(edge_weight)))



########################## Look at all linkages ################################


# OTU with the strongest linkages

edges %>% arrange(dplyr::desc(edge_weight))


df_avg_edge_weight <- edges %>% group_by(from) %>% 
  summarise(avg_edge_weight = round(mean(edge_weight),3),
            num_edge = n()) %>% 
  arrange(dplyr::desc(avg_edge_weight)) %>% data.frame()
