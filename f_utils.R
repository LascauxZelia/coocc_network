################################################################################
# Utility functions used in other scripts
################################################################################


reformat_as_df <- function(input_list, new_var_name = NULL) {
  out <- lapply(names(input_list), function(x) {
    mutate(input_list[[x]], new_var = as.character(rep(x, nrow(input_list[[x]]))))
  }) %>% bind_rows() 
  names(out)[names(out) == "new_var"] <- new_var_name
  out
}


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


# cap_var_props returns the variance explained by the CPCoA transformation 
# for a given CCA object

cap_var_props <- function(cca){
  
  eig_tot <- sum(cca$CCA$eig)
  var_propdf <- cca$CCA$eig / eig_tot
  return(round(var_propdf,3))
}

# cca_ci returns the confidence interval for the variance explained by a 
# CPCoA transformation for a given CCA object using a ANOVA-based permutation 
# test with 5000 permutations by default

cca_ci <- function(cca, permutations = 5000) {
  
  var_tbl <- variability_table(cca)
  p <- permutest(cca, permutations=permutations)
  ci <- quantile(p$F.perm, 
                 c(.05, .95)) * p$chi[1] / var_tbl["total", "inertia"]
  return(ci)
  
}




setrge<-function(x,p=5) {

  # additional range 
  plus<-(p/100)*(max(x,na.rm=T)-min(x,na.rm=T))

  # lower limit
  low<-min(x,na.rm=T)- plus

  # upper limit
  up<-max(x,na.rm=T)+ plus

  return(c(low,up))
} # end of function
