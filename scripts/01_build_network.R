################################################################################
# R script to build networks 
################################################################################


# Here I'll take a case similar to Lascaux, that is where you have several 
#  datasets obtained through barcoding of different domains of life, which makes
#  things a little more complex than if you have a single OTU table


# LOAD ENVIRONMENTAL AND OTU DATA ##############################################


# Categorical metadata ----

tab_metadata <- read.csv(paste0(dir_data, "table_metadata.csv"))
row.names(tab_metadata) <- tab_metadata$sample


# microbial composition ----

tab_otu <- read.csv(paste0(dir_data, "database_otus_cleaned.csv"), row.names = 1)


# separate the taxonomy and the abundance data ----

tab_taxonomy <- tab_otu[, 29:36] %>% rownames_to_column(var = "OTU") %>% 
  select(Domain, Phylum, Class, Order, Family, Genus, Species, OTU, taxonomy) %>% 
  as.matrix()
row.names(tab_taxonomy) <- tab_taxonomy[, "OTU"]

tab_otu <- tab_otu[, 1:28]


# Create a folder to save the results from today' analysis ----

dir_save <- paste0(dir_build_networks, "spieceasi_2021-01-20/")
dir.create(dir_save)


####  PREPARE DATA FOR ANALYSIS ################################################


### Prepare the data ===========================================================

# store all data into a phyloseq object ----

ps_full <- phyloseq(otu_table(tab_otu, taxa_are_rows = TRUE),
                    tax_table(tab_taxonomy),
                    sample_data(tab_metadata))

# Remove OTU observed in 3 samples or less  ----

ps_filt <- filter_taxa(ps_full, function(x) { 
  sum(x != 0) >= (0.3 * length(x)) 
  }, TRUE)


# further split pĥyloseq object by domains into 3 distinct objects ----

nms_domains <- unique(tab_taxonomy[, "Domain"])

masks <- lapply(nms_domains, function(x) {
  tmp <- ps_filt@tax_table[, "Domain"] == x 
  tmp %>% as.vector()
}) %>% 
  setNames(nms_domains)
  
ps_domain <- lapply(masks, function(XX) { prune_taxa(XX, ps_filt) })



### Transform the data  ========================================================
# This is to add abundance data to the results of networks (i.e. nodes abundances)


# centered log ratio ----

ls_ps_clr <- lapply(ps_domain, function(X) {
  microbiome::transform(X, "clr")
})

# compositional ----

ls_ps_compo <- lapply(ps_domain, function(X) {
  microbiome::transform(X, "compositional")
})

# save the data ----

saveRDS(ls_ps_compo, file = paste0(dir_save, "list_phyloseq_compositional.rds"))
saveRDS(ls_ps_clr, file = paste0(dir_save, "list_phyloseq_clr.rds"))




#### RUN THE SPIECEASI PIPELINE ################################################

### SpiecEasi with multiple data sets - Cross-domain interactions ==============

# SpiecEasi now includes a convenience wrapper for dealing with multiple taxa 
# sequenced on the same samples, such as 16S and ITS, 
# as seen in [Tipton, Müller, et. al. (2018)](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0393-0). 
# It assumes that each taxa is in it's own data matrix and that all samples are 
#  in all data matrices in the same order.
# Quote from this paper:
# In cross-domain studies, different marker genes are amplified and sequenced 
# separately and, hence, do not compete for reads. This implies that technically 
# independent compositions are generated in these studies. Therefore, a naive 
# application of Eq. (4) directly to the combined dataset matrix generated from a 
# simple concatenation of two compositional datasets, would be inappropriate.

# Choose the method and define the directory to save the results ----
# Two methods are proposed to generate the networks, for the details take a look 
#  at the package help. I tend to do both and compare

for ( mthd in c("mb", "glasso")) {
  
  # directory to save results from one method
  # mthd <- "glasso"
  dir_save_res <- paste0(dir_save, mthd, "/")
  dir.create(dir_save_res)

  # Define parameters for bstars method ----
  # Not optimized here but you can try to play with these
  
  if (mthd == "glasso") {
    pars_pulsar <- list(rep.num = 50, thresh = 0.05, ncores = 4)
  } else {
    pars_pulsar <- list(rep.num = 50, thresh = 0.05, ncores = 4)
  }
  
  # Run the pipeline ----
  
  res_se_multdom <-  spiec.easi(ps_domain, method = mthd, lambda.min.ratio = 1e-3,
                                nlambda = 50, sel.criterion = 'bstars', 
                                pulsar.select = TRUE, 
                                pulsar.params = pars_pulsar)
  
  
  # Save the outcome ----
  
  saveRDS(res_se_multdom, 
          file = paste0(dir_save_res, "res_spieceasi.rds"))
  
  # to reload
  # res_se_multdom <- readRDS(paste0(dir_save_res, "res_spieceasi.rds"))
  

  ### Clean up the resulting objects  ==========================================
  # Most of these "cleaning" functions are homemade
  
  ## Extract edge weights ----
  
  e_wght <- extract_edge_wght(res_se = res_se_multdom, ps = ps_filt) 
  
  
  ## Create igraph objects and add edges weights and node names ----
  
  ig  <- make_it_a_graph(res_se_multdom, ps = ps_filt) 

  
  ## Extract subgraphs for each samples ----
  
  ig_sple <- split_graph(ig, ps_filt)
 
  
  
  ####   CHARACTERIZE THE NETWORKS ###############################################

  
  # Estimate overall graph properties ============================================
  
  # at the niche level 
  
  tab_graph_metrics <- get_netwk_metrics(ig, centrality_mode =  c("out", "in", "all")[3])

  
  # for each sample network
  
  tab_graph_metrics_spl <- lapply(ig_sple, function(g) {
      get_netwk_metrics(g, centrality_mode =  c("out", "in", "all")[3])
  }) %>% reformat_as_df(new_var_name = "sample_name")

  
  # Estimate vertices attributes =================================================
  
  vertices_attr <- get_vertices_attr(ig, centrality_mode =  "all", num_iter_max_arpack = 10^5)

  vertices_attr_sple <- lapply(ig_sple, function(g) {
      get_vertices_attr(g, centrality_mode =  "all", num_iter_max_arpack = 10^5)
  })
  
  
  # Estimate edges attributes ====================================================
  
  edges_attr <- data.frame(edge_attr(ig))
  
  edges_attr_sple <- lapply(ig_sple, function(g) {
      data.frame(edge_attr(g))
  })
  
  # tansform into a table
  
  tab_edges_attr_sple <- reformat_as_df(edges_attr_sple, new_var_name = "sample_name")
  
  
  # Identify modules =============================================================
  
  res_modules <- get_modules(ig, num_iter_max_arpack = 10^6, 
                             update_rule_spin = "config",
                             spinglass_implementation = "orig",
                             num_step_wt = 4, num_trials_imap = 10)
  
  res_modules_sple <- lapply(ig_sple, function(g) {
      get_modules(g, num_iter_max_arpack = 10^6, update_rule_spin = "config",
                  spinglass_implementation = "orig",
                  num_step_wt = 4, num_trials_imap = 10)
  })
  
  # Extract module information ---------------------------------------------------
  
  clusters <- res_modules$clusters
  
  clusters_sple <- lapply(res_modules_sple, function(X) { X$clusters })

  
  # format as table
  
  clusters_length <-  clusters$clusters_length
  
  tab_clusters_length_sple <- lapply(res_modules_sple, function(X) {
      X$clusters_length %>% as.data.frame
  }) %>% reformat_as_df(new_var_name = "sample_name")

  
  # Estimate nodes connections within and outside modules ------------------------
  
  module_linkage <- get_inout_module_linkage(ig, clusters, naming_var = "OTU")
  
  module_linkage_sple <- lapply(names(ig_sple), function(x) {
      out <- get_inout_module_linkage(ig_sple[[x]],
                                      clusters_sple[[x]],
                                      naming_var = "OTU")
      out
    }) %>% setNames(names(ig_sple))

  
  
  # Format all the results for export ============================================
  
  
  # Add module membership, nodes connections and taxonomy to vertices attributes -
  
  vertices_attr <- cbind(vertices_attr,
                         res_modules$clusters_tab,
                         module_linkage)
 
  
  vertices_attr_sple <- lapply(names(vertices_attr_sple), function(x) {
      cbind(vertices_attr_sple[[x]], 
            res_modules_sple[[x]]$clusters_tab,
            module_linkage_sple[[x]]
      )
    }) %>% setNames(names(vertices_attr_sple))
  
  
  
  # add nodes attributes to the igraphs objects ----------------------------------
  
  nms_vars <- names(vertices_attr)
  
  for (i in nms_vars) {
      ig <- set_vertex_attr(ig, name = i, value = vertices_attr[, i])
  }

  
  nms_vars <- names(vertices_attr_sple[[1]][[1]])
  
    nms_sple <- names(vertices_attr_sple)
    for (spl in nms_sple) {
      for (i in nms_vars) {
        ig_sple[[spl]] <- set_vertex_attr(ig_sple[[spl]], 
                                                     name = i,
                                                     value =vertices_attr_sple[[spl]][, i])
      }
    }
  
  
  # Make tables with ALL the vertices attributes ---------------------------------
  
  
  tab_vertices_attr <- data.frame(vertex_attr(ig))  
  
  tab_vertices_attr_sple <- lapply(ig_sple, function(X) {
     data.frame(vertex_attr(X)) 
  }) %>% reformat_as_df(new_var_name = "sample_name")
  
  
  
  
  #### Save the results ==========================================================
  
  
  # graph level metrics ----
  
  write.csv(tab_graph_metrics, 
            file = paste0(dir_save_res, "results_table_graph_level_metrics.csv"))
  write.csv(tab_graph_metrics_spl, 
            file = paste0(dir_save_res, "results_table_graph_level_metrics_sple.csv"))
  
  # clusters/modules results ----
  
  write.csv(tab_clusters_length_sple,
            file = paste0(dir_save_res, "results_table_clusters_length_sple.csv"))
  saveRDS(clusters, 
          file = paste0(dir_save_res, "results_detected_clusters.rds"))
  saveRDS(clusters_sple,
          file = paste0(dir_save_res, "results_detected_clusters_sple.rds"))
  
  # full vertices attributes ----
  
  saveRDS(tab_vertices_attr, 
          file = paste0(dir_save_res, "results_vertices_attributes.rds"))
  saveRDS(tab_vertices_attr_sple, 
          file = paste0(dir_save_res, "results_vertices_attributes_sple.rds"))
  
  # full edges attributes  ----
  
  write.csv(edges_attr, 
            file = paste0(dir_save_res, "results_edges_attributes.csv"))
  write.csv(tab_edges_attr_sple, 
            file = paste0(dir_save_res, "results_edges_attributes_sple.csv"))
  
  # the igraph objects ----
  
  saveRDS(ig, 
          file = paste0(dir_save_res, "list_networks.rds"))
  saveRDS(ig_sple, 
          file = paste0(dir_save_res, "list_networks_sple.rds"))
  
} ##eo for loop on mthd
