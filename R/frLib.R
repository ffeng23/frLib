###frLib, an accessory R software package for Feng's research work.
#	By Feng @ BU 05/19/2019 Boston University (ffeng@bu.edu)
#
#  Version 0.1
#  For now, 1)it will only try to do some IgSeq data analysis
#
#### Developed By Feng @ BU 2019. All right reserved###

#Below are roxygen2 comments

##########Dependency###########
##  library(MASS)
####################
#this is where the files to be include if necessary
#### # '  @include ELISAtools_IO.R
#'
#'
#' @title Accessory R library
#'
#' @description An R package containing some useful accessory R functions
#'		
#'
#' @details This package is currentlydeveloped to run analysis of IgSeq data. 
#'		
#' 		Please refer to the vignettes to see details.
#' @references Feng, et al 2019 \url{https://doi.org/10.1101/483800}		
"_PACKAGE"


.onUnload <- function (libpath) {
	cat("calling to unload package lib...........\n")
  library.dynam.unload("frLib", libpath)
}