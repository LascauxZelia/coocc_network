

# this function check whether the list of packages provided as input is installed
#  and if not then it installs the packages

# 3 types of packages sources are allowed: CRAN, github and bioconductor

install_if_not_there <- function(pkg_list, type = "CRAN") {
  
  inst <- pkg_list %in% installed.packages()
  
  if (any(! inst)) {
    
    if (type == "CRAN") {
      install.packages(cran_packages[!inst], repos = "http://cran.rstudio.com/")
    }
    
    if (type == "github") {
      devtools::install_github(github_packages[!inst])
    }
  
    if (type == "bioconductor") {
      source("http://bioconductor.org/biocLite.R")
      BiocManager::install(bioc_packages[!inst])
    }
  
  }
}