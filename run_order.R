# Packages

library(fda) 
library(quantreg)
library(matrixStats)
library(Matrix)
library(refund)
library(nloptr)
library(expm)
library(here)

# Load the source files
run_order <- function(){
  source(here("dgp1.R"))
  source(here("dgp2.R"))
  source(here("dgp3.R"))
  source(here("auxilary_functions_FPLS_FPCA.R"))
  source(here("auxilary_functions_fpca.R"))
  source(here("auxilary_functions_fpqr.R"))
  source(here("auxilary_functions_pfqr.R"))
}
