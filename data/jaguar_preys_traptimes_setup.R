
library(overlap)
library(boot) 
library(circular)

##==========================================================
## This program sets up the jaguar-preys camera trap 
## data and loads functions that are needed for subsequent 
## analysis
##
## Note: You need to have 2 additional R packages installed,
##               boot   and   circular
##
##       If you get an error message saying that these are
##       not installed, use the R command
##               install.packages( c("boot", "circular))
##       to install them.
##
##       Alternatively, you can use Load packages ...
##       on the Packages menu.
##==========================================================


    ##--------------------------------------------------
    ## Read the contents of a file that contains various 
    ## functions that are used in the analysis
    ##
    ## Note that this requires two additional R packages
    ## boot and circular  
    ##--------------------------------------------------


source("jaguar_preys_traptimes_ovlcode.r")


   ##----------------------------------------------------------
   ## Read the raw data and create a variable for each species, 
   ## in circular format and with missing values removed
   ##
   ## In the data file there is a column for each species
   ##----------------------------------------------------------


rawtimes <- read.table("jaguar_preys_traptimes_circular.txt")
times <- list(1)
for (j in c(1,5,6,7,8)) {
    tmp <- circular(2*pi*na.omit(rawtimes[,j]))
    times[[j]] <- tmp
}

ow <- options("warn")
options(warn = -1)
ponca = as.circular(na.omit(as.vector(times[[1]])))
tay = as.circular(na.omit(as.vector(times[[5]])))
pec = as.circular(na.omit(as.vector(times[[6]])))
tapir = as.circular(na.omit(as.vector(times[[7]])))
maz.a = as.circular(na.omit(as.vector(times[[8]])))
options(ow) # reset
