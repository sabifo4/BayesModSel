#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir
path_to_file <- getActiveDocumentContext()$path
wd <- paste( dirname( path_to_file ), "/", sep = "" )
wd <- gsub( pattern = "scripts/", replacement = "HPC/MCMCtree/", x = wd )
setwd( wd )

#----------------#
# LOAD LIBRARIES #
#----------------#
# Install if needed by uncommenting line below
# devtools::install_github("dosreislab/mcmc3r", , build_vignettes = TRUE )
library( mcmc3r )

#------------------------------------------------------------------#
# ANALYSE MODELS BASED ON DIFFERENT POINTS OF VIEW ON INTERPRETING #
#                     THE FOSSIL RECORD                            #
#------------------------------------------------------------------#
# Come back to wd
setwd( wd )
# Add seed number for the posterior probabilities
set.seed( 12345 )

## ------------------------------------------------------------------------- ##
## H1: Metazoa is not allowed to be as old as 833 Ma. dRetal15 | pU = 1e-300 ##
## H2: Metazoa is not allowed to be older than 580 Ma. B&M24   | pU = 1e-300 ##
## ------------------------------------------------------------------------- ##

# 1. Get data

##\\ [[ dRetal15pUhb_ILN ]]
#>> Samples collected under hypothesis "dRetal15pUhb"
#>> Maximum age of the root with a hard bound
#>> Clock model: ILN
setwd( "dRetal15pUhb_ILN" )
dRetal15pUhb_ILN <- mcmc3r::stepping.stones( mcmcf = "mcmc.txt",
                                             betaf = "betaweights.txt" )
save( dRetal15pUhb_ILN, file = paste( wd, "dRetal15pUhb_ILN_b128.RData",
                                      sep = "" ) )
# load( paste( wd, "dRetal15pUhb_ILN_b128.RData", sep = "" ) )

##\\ [[ BM23pUhb_ILN ]]
#>> Samples collected under hypothesis "BM23pUhb", same as "BM23" but
#>> Maximum age of the root with a hard bound
#>> Clock model: ILN
setwd( "../BM23pUhb_ILN" )
BM23pUhb_ILN <- mcmc3r::stepping.stones( mcmcf = "mcmc.txt",
                                         betaf = "betaweights.txt" )
save( BM23pUhb_ILN, file = paste( wd, "BM23pUhb_ILN_b128.RData", sep = "" ) )
# load( paste( wd, "BM23pUhb_ILN_b128.RData", sep = "" ) )

##\\ [[ dRetal15pUhb_GBM ]]
#>> Samples collected under hypothesis "dRetal15pUhb"
#>> Maximum age of the root with a hard bound
#>  Clock model: GBM
setwd( "../dRetal15pUhb_GBM" )
dRetal15pUhb_GBM <- mcmc3r::stepping.stones( mcmcf = "mcmc.txt",
                                             betaf = "betaweights.txt" )
save( dRetal15pUhb_GBM, file = paste( wd, "dRetal15pUhb_GBM_b128.RData",
                                      sep = "" ) )
# load( paste( wd, "dRetal15pUhb_GBM_b128.RData", sep = "" ) )

##\\ [[ BM23pUhb_GBM ]]
#>> Samples collected under hypothesis "BM23pUhb"
#>> Maximum age of the root with a hard bound
#>> Clock model: GBM
setwd( "../BM23pUhb_GBM" )
BM23pUhb_GBM <- mcmc3r::stepping.stones( mcmcf = "mcmc.txt",
                                         betaf = "betaweights.txt" )
save( BM23pUhb_GBM, file = paste( wd, "BM23pUhb_GBM_b128.RData", sep = "" ) )
# load( paste( wd, "BM23pUhb_GBM_b128.RData", sep = "" ) )

# 2. Compare models
# Calculate mnln, BFs, and P
mlnl          <- c( dRetal15pUhb_GBM$logml, dRetal15pUhb_ILN$logml,
                    BM23pUhb_GBM$logml, BM23pUhb_ILN$logml )
names( mlnl ) <- c( "dRetal15pUhb-GBM", "dRetal15pUhb-ILN",
                    "BM23pUhb-GBM", "BM23pUhb-ILN" )
mlnl
# dRetal15pUhb-GBM dRetal15pUhb-ILN     BM23pUhb-GBM     BM23pUhb-ILN 
# -674.7174        -697.7024        -720.2386        -764.1428 
se          <- c( dRetal15pUhb_GBM$se, dRetal15pUhb_ILN$se,
                  BM23pUhb_GBM$se, BM23pUhb_ILN$se )
names( se ) <- c( "dRetal15pUhb-GBM", "dRetal15pUhb-ILN",
                  "BM23pUhb-GBM", "BM23pUhb-ILN" )
se
# dRetal15pUhb-GBM dRetal15pUhb-ILN     BM23pUhb-GBM     BM23pUhb-ILN 
# 0.08513776       0.05714799       0.08081300       0.05983390 
BF          <- exp( mlnl - max( mlnl ) )
names( BF ) <- c( "dRetal15pUhb-GBM", "dRetal15pUhb-ILN",
                  "BM23pUhb-GBM", "BM23pUhb-ILN" )
BF
# dRetal15pUhb-GBM dRetal15pUhb-ILN     BM23pUhb-GBM     BM23pUhb-ILN 
# 1.000000e+00     1.041701e-10     1.699836e-20     1.455657e-39 
Pr   <- BF / sum( BF )
Pr
# dRetal15pUhb-GBM dRetal15pUhb-ILN     BM23pUhb-GBM     BM23pUhb-ILN 
# 1.000000e+00     1.041701e-10     1.699836e-20     1.455657e-39
# Use two best models
# 2lnBF_01 = 2( lnml.M0 - lnml.M1 )
BFs  <- 2*( dRetal15pUhb_GBM$logml - dRetal15pUhb_ILN$logml )
BFs
# [1] 45.96999
#
# Ideally, you want the S.E. to be much smaller than
# the log-marginal likelihood difference between the models being tested.
dRetal15pUhb_GBM$logml - dRetal15pUhb_ILN$logml #OK ?
# [1] 22.985
#
# Recently, a function exists to do the steps followed above:
sumBF <- mcmc3r::bayes.factors( dRetal15pUhb_GBM, dRetal15pUhb_ILN,
                                BM23pUhb_GBM, BM23pUhb_ILN )
sumBF
# $bf
# [1] 1.000000e+00 1.041701e-10 1.699836e-20 1.455657e-39
# 
# $logbf
# [1] 0.00000 -22.98500 -45.52117 -89.42536
# 
# $pr
# [1] 1.000000e+00 1.041701e-10 1.699836e-20 1.455657e-39
# 
# $prior
# [1] 0.25 0.25 0.25 0.25
# 
# $pr.ci
# 2.5%        97.5%
#   [1,] 1.000000e+00 1.000000e+00
# [2,] 8.535895e-11 1.269208e-10
# [3,] 1.351798e-20 2.147858e-20
# [4,] 1.190433e-39 1.787259e-39

#----------------#
# OUTPUT RESULTS #
#----------------#
num_mod <- 4
out_mat <- matrix( 0, nrow = num_mod, ncol = 6 )
colnames( out_mat ) <- c( "Diversification model", "logL", "S.E. (delta)",
                          "Pr(M|D)", "2.5% CI", "97.5% CI" )
out_mat[1:num_mod,1] <- c( "Unconstrained | pU = 1e-300 | dRetal15pUhb-GBM",
                           "Unconstrained | pU = 1e-300 | dRetal15pUhb-ILN",
                           "Late Ediacaran | pU = 1e-300 | BM23pUhb-GBM",
                           "Late Ediacaran | pU = 1e-300 | BM23pUhb-ILN" )
count <- 0
for( i in 1:num_mod ){
  out_mat[i,2:6] <- c( round( mlnl[i], 2), round( se[i], 3 ),
                       sumBF$pr[i], sumBF$pr.ci[i,] )
}
# Come back to wd and save out file 
setwd( wd )
write.table( x = out_mat, file = "out_animals_BFs_b128.tsv", sep = "\t",
             quote = FALSE, row.names = FALSE )

#--------------#
# PLOT RESULTS #
#--------------#
# Plot beta points against mean logL
pdf( file = "BpointsVSmlogL_DivModelAnimals.pdf", paper = "a4r",
     width = 0, height = 0 )
par( mfrow = c(1,2) )
##\\ dRetal15pUhb_GBM
plot( dRetal15pUhb_GBM$b, dRetal15pUhb_GBM$mean.logl, pch = 19,
      col = rgb(0,0,0,alpha=0.3),
      xaxs = "i", xlim = c(0,1), xlab = "b", ylab = "mean logL",
      main = "Unconstrained | pU = 1e-300 | dRetal15pUhb-GBM",
      cex = 0.8, ylim = c(-6e5,0),
      cex.main = 1 )
lines( dRetal15pUhb_GBM$b, dRetal15pUhb_GBM$mean.logl )
##\\ dRetal15pUhb_ILN
plot( dRetal15pUhb_ILN$b, dRetal15pUhb_ILN$mean.logl, pch=19, 
      col = rgb(0,0,0,alpha=0.3),
      xaxs = "i", xlim = c(0,1), xlab = "b", ylab = "mean logL",
      main = "Unconstrained | pU = 1e-300 | dRetal15pUhb-ILN",
      cex = 0.8, ylim = c(-6e5,0),
      cex.main = 1 )
lines( dRetal15pUhb_ILN$b, dRetal15pUhb_ILN$mean.logl )
dev.off()
