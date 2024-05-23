################################################################################
# Function to manipulate igraph objects
################################################################################

# extract edge weight and rename them according to the taxa name ----

# res_se: result of a spieceasi analysis
# ps: the phyloseq object used to make the spieceasi analysis

extract_edge_wght <- function(res_se, ps) {
  
  # OTU/ASV names
  nms <- row.names(ps@otu_table)
  
  # extract edges and rename them according to taxa names
  
  if (res_se$est$method == "mb") {
    
    out <- summary(symBeta(getOptBeta(res_se), mode = 'maxabs')) %>% 
      dplyr::rename("from" = i, "to" = j, "edge_weight" = x) %>%
      dplyr::mutate(from = nms[from], 
                    to = nms[to],
                    edge_nm = paste0(from, "_|_", to)) %>% 
      dplyr::select(from, to, edge_nm, edge_weight)
  }
  
  if (res_se$est$method == "glasso") {
    
    secor  <- cov2cor(getOptCov(res_se))
    out <- summary(triu(secor*getRefit(res_se), k=1)) %>% 
      dplyr::rename("from" = i, "to" = j, "edge_weight" = x) %>%
      dplyr::mutate(from = nms[from], 
                    to = nms[to],
                    edge_nm = paste0(from, "_|_", to)) %>% 
      dplyr::select(from, to, edge_nm, edge_weight)
  }
  
  out$edge_binary <- sapply(out$edge_weight, function(x) {
    if (x > 0) { 1 } else { -1 }
  })
  
  return(out)
}

# transform the result of a spieceasi analysis into an igraph object and add
# the edge weights and node names ----

make_it_a_graph <- function(input_spieceasi, ps) {
  
  ew <- extract_edge_wght(res_se = input_spieceasi, ps = ps)
  va <- data.frame(ps@tax_table@.Data) %>% as.list() 
  va <- lapply(va, as.character)
  out <- adjmat_to_igraph(getRefit(input_spieceasi), edge.attr = ew,
                          vertex.att = va, wght = NULL)
  return(out)
}

# Fonction to replace adj2igraph as the original one depend on the function
#  graph.adjacency for which we cannot specify arguments

adjmat_to_igraph <- function (Adj, rmEmptyNodes = FALSE, diag = FALSE, 
                              edge.attr = list(), wght = NULL, 
                              vertex.attr = list(name = 1:ncol(Adj))){
  g <- igraph::graph.adjacency(Adj, mode = "undirected", weighted = wght, 
                               diag = diag)
  if (length(vertex.attr) > 0) {
    for (i in 1:length(vertex.attr)) {
      attr <- names(vertex.attr)[i]
      g <- igraph::set.vertex.attribute(g, attr, index = igraph::V(g), 
                                        vertex.attr[[i]])
    }
  }
  if (length(edge.attr) > 0) {
    for (i in 1:length(edge.attr)) {
      attr <- names(edge.attr)[i]
      g <- igraph::set.edge.attribute(g, attr, index = igraph::E(g), 
                                      edge.attr[[i]])
    }
  }
  if (rmEmptyNodes) {
    ind <- igraph::V(g)$name[which(igraph::degree(g) < 1)]
    g <- igraph::delete.vertices(g, ind)
  }
  g
}

# split an igraph object into several subgraphs 
#   that correspond to each separated sample ----

split_graph <- function(g, ps, group = NULL) {
  
  if (is.null(group)) {
    nms_sples <- sample_names(ps)
    
    ls_g <- lapply(nms_sples, function(spl) {
      vec <- ps@otu_table[,spl]
      g_out <- induced_subgraph(g, which(vec != 0))
      set_vertex_attr(g_out, name = "original_abundance", 
                      index = V(g_out), value = vec[which(vec != 0)])
    }) %>% setNames(nms_sples)
  } else {
    meta <- sample_data(ps)[, group] %>% unlist() %>% as.vector()
    nms_gps <- meta %>% unique()
    otu <- ps@otu_table@.Data
    
    ls_g <- lapply(nms_gps, function(gp) {
      mask <- meta == gp
      vec <- rowMeans(otu[,mask])
      g_out <- induced_subgraph(g, which(vec != 0))
      set_vertex_attr(g_out, name = "original_abundance", 
                      index = V(g_out), value = vec[which(vec != 0)])
    }) %>% setNames(nms_gps)
  } 
  
  return(ls_g)
}

# remove the edges with weight between quantiles (q1, q4) values 
#  to keep only strong associations
trim_edges_by_quantiles <- function(g) { 
  ew <- edge_attr(g)$edge_weight
  thds <- quantile(ew)[c(2,4)]
  mask <- data.table::between(ew, thds[1], thds[2])
  ew[! mask]
  delete_edges(g, which(mask))
}



# load results of network analysis

load_network_results <- function(dir_load, sple = TRUE) {
  
  if (sple) {
    ext_rds <- "_sple.rds"
    ext_csv <- "_sple.csv"
  } else {
    ext_rds <- ".rds"
    ext_csv <- ".csv"
  }
  
  ls_data <- list()
  
  # graph level metrics ----
  path <- paste0(dir_load, "results_table_graph_level_metrics", ext_csv)
  ls_data$graph_metrics <- read.csv(path, row.names = 1)
  
  # full vertices attributes ----
  
  path <- paste0(dir_load, "results_vertices_attributes", ext_rds)
  ls_data$vertices_attr <- readRDS(path)
  
  # full edges attributes ----
  
  path <-  paste0(dir_load, "results_edges_attributes", ext_csv)
  ls_data$edges_attr <- read.csv(path, row.names = 1)
  
  # the igraph objects ----
  
  path <- paste0(dir_load, "list_networks", ext_rds)
  ls_data$networks <- readRDS(path)
  
  return(ls_data)
}


# Functions to get networks metrics

get_netwk_metrics <- function(g, centrality_mode =  c("out", "in", "all")[3]) {
  
  # number of nodes
  num_nodes <- vcount(g)
  # number of links
  num_edges <- ecount(g)
  # is the graph directed ?
  is_g_directed <- is_directed(g)
  # is the graph connected ?
  is_g_connected <- is_connected(g)
  # is the graph weighted ?
  is_g_weighted <- is_weighted(g)
  if (is_g_weighted) { g_weights <- E(g)$weights } else { g_weights <- NULL }
  
  # degree
  centr_deg     <- centr_degree(graph = g, mode = centrality_mode, loops = FALSE, 
                                normalized = FALSE)$centralization
  centr_deg_std <- centr_degree(graph = g, mode = centrality_mode, loops = FALSE, 
                                normalized = TRUE)$centralization 
  # betweenness
  centr_betw     <- centr_betw(graph = g, directed = is_g_directed, nobigint = TRUE, 
                               normalized = FALSE)$centralization
  centr_betw_std <- centr_betw(graph = g, directed = is_g_directed, nobigint = TRUE, 
                               normalized = TRUE)$centralization
  # closeness
  centr_clo     <- centr_clo(g, mode = centrality_mode, 
                             normalized = FALSE)$centralization
  centr_clo_std <- centr_clo(graph = g, mode = centrality_mode, 
                             normalized = TRUE)$centralization  
  
  # eigen vector
  centr_eig     <- centr_eigen(g, directed = is_g_directed, scale = TRUE,
                               options = arpack_defaults, normalized = FALSE)$centralization
  centr_eig_std <- centr_eigen(g, directed = is_g_directed, scale = TRUE,
                               options = arpack_defaults, normalized = TRUE)$centralization
  
  # Graph diameter
  diameter_length <- diameter(g, directed = is_g_directed, 
                              unconnected = !is_g_connected, weights = g_weights)
  
  # estimate hub scores
  hub_score_graph <- hub_score(g, scale = TRUE, weights = g_weights)
  
  # Edge density
  edg_density <- num_edges / num_nodes
  
  # Average distance between nodes
  avg_dist <- mean_distance(g, directed = is_g_directed)
  
  # Edge connectance: num_edges / (num_nodes * ((num_nodes - 1)/2))
  # the function edge_density(g, loops = FALSE) measures connectance !!
  edg_connectance <- num_edges / (num_nodes * ((num_nodes - 1)/2))
  
  # make a table of network-level metrics
  ntwk_metrics <- data.frame(num_nodes = num_nodes,
                             num_edges = num_edges,
                             is_directed  = is_g_directed %>% as.logical(),
                             is_connected = is_g_connected %>% as.logical(),
                             is_weighted  = is_g_weighted %>% as.logical(),
                             diameter     = diameter_length,
                             density      = edg_density,
                             connectance  = edg_connectance,
                             avg_distance = avg_dist,
                             centralized_degree      = centr_deg,
                             centralized_betweenness = centr_betw,
                             centralized_closeness   = centr_clo,
                             centralized_eigenvector = centr_eig,
                             centralized_degree_std      = centr_deg_std,
                             centralized_betweenness_std = centr_betw_std,
                             centralized_closeness_std   = centr_clo_std,
                             centralized_eigenvector_std = centr_eig_std,
                             hub_score_graph = hub_score_graph$value
  )
  return(ntwk_metrics)
}

# Define function to estimate vertices attributes -----------------------------

get_vertices_attr <- function(g, centrality_mode =  c("out", "in", "all")[3],
                              num_iter_max_arpack = 10^5) {
  
  # is the graph directed ?
  is_g_directed <- is_directed(g)
  # is the graph connected ?
  is_g_connected <- is_connected(g)
  # is the graph weighted ?
  is_g_weighted <- is_weighted(g)
  if (is_g_weighted) { g_weights <- E(g)$weights } else { g_weights <- NULL }
  
  # degree
  nodes_deg     <- degree(g, mode = centrality_mode, normalized = FALSE)
  nodes_deg_std <- degree(g, mode = centrality_mode, normalized = TRUE)
  
  # betweenness
  nodes_betw     <- betweenness(g, directed = is_g_directed, weights = g_weights, 
                                nobigint = TRUE, normalized = FALSE)
  nodes_betw_std <- betweenness(g, directed = is_g_directed, weights = g_weights, 
                                nobigint = TRUE, normalized = TRUE)
  
  # closeness
  nodes_clo     <- closeness(g, normalized = FALSE)
  nodes_clo_std <- closeness(g, normalized = TRUE)
  
  # eigen vector
  nodes_eig     <- eigen_centrality(g, directed = is_g_directed, scale = FALSE,
                                    weights = g_weights, options = arpack_defaults)$vector
  nodes_eig_std <- eigen_centrality(g, directed = is_g_directed, scale = TRUE,
                                    weights = g_weights, options = arpack_defaults)$vector
  
  
  # estimate hub scores
  nodes_hub_scor <- hub_score(g, scale = TRUE, weights = g_weights)
  
  # estimate authority score
  arpack_defaults$maxiter <- num_iter_max_arpack
  nodes_auth_scor <- authority_score(g, scale = TRUE, weights = g_weights,
                                     options = arpack_defaults)
  
  
  # make a table with centrality measures at the node level
  centrality_tab <- data.frame(degree      = nodes_deg,
                               betweenness = nodes_betw,
                               closeness   = nodes_clo,
                               eigenvector = nodes_eig,
                               degree_std  = nodes_deg_std,
                               betweenness_std = nodes_betw_std,
                               closeness_std   = nodes_clo_std,
                               eigenvector_std = nodes_eig_std,
                               hub_score       = nodes_hub_scor$vector,
                               authority_score = nodes_auth_scor$vector
  )
  return(centrality_tab %>% rownames_to_column(var = "otu")) 
}


# Define function to estimate edges attributes -----------------------------

get_edges_attr <- function(g, centrality_mode =  c("out", "in", "all")[3]) {
  
  # is the graph directed ?
  is_g_directed <- is_directed(g)
  # is the graph connected ?
  is_g_connected <- is_connected(g)
  # is the graph weighted ?
  is_g_weighted <- is_weighted(g)
  if (is_g_weighted) { g_weights <- E(g)$weights } else { g_weights <- NULL }
  
  # extract edges attributes, here correlation coeff and pval
  edg_attr <- edge_attr(g) %>% bind_cols()
  
  # # Path of the graph diameter
  # diameter_path <- get_diameter(g, directed = is_g_directed, 
  #                               unconnected = !is_g_connected, weights = g_weights)
  
  return(edg_attr)
}

# Define function to identifies modules ----------------------------------------

get_modules <- function(g, num_iter_max_arpack = 10^6, 
                        update_rule_spin = c("config", "random", "simple")[1],
                        spinglass_implementation = c("orig", "neg")[1],
                        num_step_wt = 4, num_trials_imap = 10) {
  
  # is the graph directed ?
  is_g_directed <- is_directed(g)
  # is the graph connected ?
  is_g_connected <- is_connected(g)
  # is the graph weighted ?
  is_g_weighted <- is_weighted(g)
  if (is_g_weighted) { g_weights <- E(g)$weights } else { g_weights <- NULL }
  
  # Run various clustering algorithms
  
  clust_eb <- cluster_edge_betweenness(g, weights = g_weights, directed = is_g_directed, 
                                       edge.betweenness = TRUE, merges = TRUE,
                                       bridges = TRUE, modularity = TRUE, 
                                       membership = TRUE)
  
  clust_fg <- cluster_fast_greedy(g, merges = TRUE, modularity = TRUE,
                                  membership = TRUE, weights = g_weights)
  
  clust_imap <- cluster_infomap(g, e.weights = g_weights, v.weights = NULL,
                                nb.trials = num_trials_imap, modularity = TRUE)
  
  clust_lp <- cluster_label_prop(g, weights = g_weights, initial = NULL,
                                 fixed = NULL)
  
  arpack_defaults$maxiter <- num_iter_max_arpack
  clust_le <- cluster_leading_eigen(g, steps = -1, weights = g_weights, start = NULL, 
                                    options = arpack_defaults, callback = NULL,
                                    extra = NULL, env = parent.frame())
  
  clust_lou <- cluster_louvain(g, weights = g_weights)
  
  # clust_opt <- cluster_optimal(g, weights = g_weights)
  
  if (is_g_connected) {
    clust_spin <- cluster_spinglass(g, weights = g_weights, vertex = NULL, spins = 25,
                                    parupdate = FALSE, start.temp = 1, stop.temp = 0.01,
                                    cool.fact = 0.99, gamma = 1, gamma.minus = 1,
                                    update.rule = update_rule_spin,
                                    implementation = spinglass_implementation)
  }
  
  clust_wt <- cluster_walktrap(g, weights = g_weights, steps = num_step_wt,
                               merges = TRUE, modularity = TRUE, membership = TRUE)
  
  # format the results
  
  clusters <- list(edge_betweenness = clust_eb, fast_greedy = clust_fg, 
                   infomap = clust_imap, label_prop = clust_lp, 
                   leading_eigen = clust_le, louvain = clust_lou,
                   # optimal = clust_opt, spinglass = clust_spin, 
                   walktrap = clust_wt)
  
  clusters_tab <- lapply(clusters, membership) %>% bind_cols()
  
  clusters_length <- map_df(clusters, length)
  
  return(list(clusters        = clusters,
              clusters_tab     = clusters_tab,
              clusters_length = clusters_length
  ))
}

# Define function to estimate nodes connections within and outside modules,
#  i.e.determine participation coefficient and z-score -------------------------

get_inout_module_linkage <- function(g, clusters, naming_var = "OTU") { 
  
  # Participation coefficient : describes the profile of a node interactions 
  #  with species found outside of the module it belongs to
  nodes_partcoeff <- lapply(clusters, function(x) {
    part_coeff(g, x %>% membership)
  }) %>% bind_cols() %>% as.data.frame()
  names(nodes_partcoeff) <- paste0("pc_", names(nodes_partcoeff))
  rownames(nodes_partcoeff) <- vertex_attr(g)[["OTU"]]
  
  # Estimate z-score, a measure of the connectivity from a given vertex to 
  #  other vertices in its module/community 
  nodes_z_score <- lapply(clusters, function(x) {
    within_module_deg_z_score(g, x %>% membership())   
  }) %>% bind_cols() %>% as.data.frame()
  names(nodes_z_score) <- paste0("z_", names(nodes_z_score))
  rownames(nodes_z_score) <- vertex_attr(g)[["OTU"]]
  
  out <- cbind(nodes_partcoeff, nodes_z_score)
  return(out) 
}
