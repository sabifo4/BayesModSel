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
## H1: the calibration on the age of Placentalia lies between 162.5 and      ##
##     61.66 Ma with probability of violation of the bounds equal to 5%      ##
##     (2.5% each bound). ACetal22 | pL = pU = 0.025                         ##
##                                                                           ##
## H2: the age of Placentalia is constrained to postdate the K-Pg, between   ##
##     66.09-61.66 Ma. B&M24 | pL = 0.025 | pU = 1e-300                      ##
## ------------------------------------------------------------------------- ##

# Get data
##\\ [[ ACetal22_GBM ]]
#>> Samples collected under "Unconstrained" diversification model | Clock: GBM
setwd( paste( wd, "ACetal22_GBM", sep = "" ) )
ACetal22_GBM <- mcmc3r::stepping.stones( mcmcf = "mcmc.txt",
                                         betaf = "betaweights.txt" )
save( ACetal22_GBM, file = paste( wd, "ACetal22_GBM_b256.RData", sep = "" ) )
# load( paste( wd, "ACetal22_GBM_b256.RData", sep = "" ) )

##\\ [[ BM23pUhb_GBM ]]
#>> Samples collected under "Post-K-Pg" diversification model | Clock: GBM
setwd( "../BM23pUhb_GBM" )
BM23pUhb_GBM <- mcmc3r::stepping.stones( mcmcf = "mcmc.txt",
                                         betaf = "betaweights.txt" )
# Warning messages:
#   1: In mcmc3r::stepping.stones(mcmcf = "mcmc.txt", betaf = "betaweights.txt") :
#   unreliable se: var(r_k)/r_k^2 = 1.05184777590924 > 0.1 for b = 0.202472
# 2: In mcmc3r::stepping.stones(mcmcf = "mcmc.txt", betaf = "betaweights.txt") :
#   unreliable se: var(r_k)/r_k^2 = 1.33320198378073 > 0.1 for b = 0.8193618
save( BM23pUhb_GBM, file = paste( wd, "BM23pUhb_GBM_b256.RData", sep = "" ) )
# load( paste( wd, "BM23pUhb_GBM_b256.RData", sep = "" ) )

##\\ [[ ACetal22_ILN ]]
#>> Samples collected under "Unconstrained" diversification model | Clock: ILN
setwd( "../ACetal22_ILN" )
ACetal22_ILN <- mcmc3r::stepping.stones( mcmcf = "mcmc.txt",
                                         betaf = "betaweights.txt" )
save( ACetal22_ILN, file = paste( wd, "ACetal22_ILN_b256.RData", sep = "" ) )
# load( paste( wd, "ACetal22_ILN_b256.RData", sep = "" ) )

##\\ [[ BM23pUhb_ILN ]]
#>> Samples collected under "Post-K-Pg" diversification model | Clock: ILN
setwd( "../BM23pUhb_ILN" )
BM23pUhb_ILN <- mcmc3r::stepping.stones( mcmcf = "mcmc.txt",
                                         betaf = "betaweights.txt" )
# Warning messages:
#   1: In mcmc3r::stepping.stones(mcmcf = "mcmc.txt", betaf = "betaweights.txt") :
#   unreliable se: var(r_k)/r_k^2 = 0.796495840140179 > 0.1 for b = 0.1718552
# 2: In mcmc3r::stepping.stones(mcmcf = "mcmc.txt", betaf = "betaweights.txt") :
#   unreliable se: var(r_k)/r_k^2 = 0.672560656485706 > 0.1 for b = 0.3714468
# 3: In mcmc3r::stepping.stones(mcmcf = "mcmc.txt", betaf = "betaweights.txt") :
#   unreliable se: var(r_k)/r_k^2 = 1.98110156975154 > 0.1 for b = 0.7706045
# 4: In mcmc3r::stepping.stones(mcmcf = "mcmc.txt", betaf = "betaweights.txt") :
#   unreliable se: var(r_k)/r_k^2 = 0.181605132759656 > 0.1 for b = 0.8532152
save( BM23pUhb_ILN, file = paste( wd, "BM23pUhb_ILN_b256.RData", sep = "" ) )
# load( paste( wd, "BM23pUhb_ILN_b256.RData", sep = "" ) )

# Calculate mnln, BFs, and P
mlnl <- c( ACetal22_GBM$logml, BM23pUhb_GBM$logml,
           ACetal22_ILN$logml, BM23pUhb_ILN$logml )
names( mlnl ) <- c( "ACetal22-GBM", "BM23pUhb-GBM",
                    "ACetal22-ILN", "BM23pUhb-ILN" )
mlnl
# ACetal22-GBM BM23pUhb-GBM ACetal22-ILN BM23pUhb-ILN 
# -2509.800    -2609.059    -2612.743    -2774.480
se   <- c( ACetal22_GBM$se, BM23pUhb_GBM$se,
           ACetal22_ILN$se, BM23pUhb_ILN$se )
names( se ) <- c( "ACetal22-GBM", "BM23pUhb-GBM",
                  "ACetal22-ILN", "BM23pUhb-ILN" )
se
# ACetal22-GBM BM23pUhb-GBM ACetal22-ILN BM23pUhb-ILN 
# 0.1571221    1.5565983    0.1235524    1.9217498 
BF <- exp( mlnl - max( mlnl ) )
names( BF ) <- c( "ACetal22-GBM", "BM23pUhb-GBM",
                  "ACetal22-ILN", "BM23pUhb-ILN" )
BF
# ACetal22-GBM  BM23pUhb-GBM  ACetal22-ILN  BM23pUhb-ILN 
# 1.000000e+00  7.802931e-44  1.960515e-45 1.123953e-115
Pr <- BF / sum( BF )
Pr
# ACetal22-GBM  BM23pUhb-GBM  ACetal22-ILN  BM23pUhb-ILN 
# 1.000000e+00  7.802931e-44  1.960515e-45 1.123953e-115
# Use best two models
# 2lnBF_01 = 2( lnml.M0 - lnml.M1 )
BFs  <- 2*( ACetal22_GBM$logml - BM23pUhb_GBM$logml )
BFs
# 198.5185
# Ideally, you want the S.E. to be much smaller than
# the log-marginal likelihood difference between the models being tested.
ACetal22_GBM$logml - BM23pUhb_GBM$logml #OK ? 
# 99.25924 # Yeah!

# Recently, a function exists to do the steps followed above:
sumBF <- mcmc3r::bayes.factors( ACetal22_GBM, BM23pUhb_GBM,
                                ACetal22_ILN, BM23pUhb_ILN )
sumBF
# $bf
# [1]  1.000000e+00  7.802931e-44  1.960515e-45 1.123953e-115
# 
# $logbf
# [1]    0.00000  -99.25924 -102.94312 -264.68043
# 
# $pr
# [1]  1.000000e+00  7.802931e-44  1.960515e-45 1.123953e-115
# 
# $prior
# [1] 0.25 0.25 0.25 0.25
# 
# $pr.ci
# 2.5%         97.5%
#   [1,]  1.000000e+00  1.000000e+00
# [2,]  4.040188e-45  1.608341e-42
# [3,]  1.328052e-45  2.926474e-45
# [4,] 2.805676e-117 4.856903e-114

#----------------#
# OUTPUT RESULTS #
#----------------#
num_mod <- 4
out_mat <- matrix( 0, nrow = num_mod, ncol = 6 )
colnames( out_mat ) <- c( "Diversification model", "logL", "S.E. (delta)",
                          "Pr(M|D)", "2.5% CI", "97.5% CI" )
out_mat[,1] <- c( "Unconstrained | ACetal22-GBM",
                  "Post-K-Pg | BM23pUhb-GBM",
                  "Unconstrained | ACetal22-ILN",
                  "Post-K-Pg | BM23pUhb-ILN" )
for( i in 1:num_mod ){
  out_mat[i,2:6] <- c( round( mlnl[i], 2 ), round( se[i], 2 ),
                       sumBF$pr[i], sumBF$pr.ci[i,] )
}

# Come back to wd and save out file 
setwd( wd )
write.table( x = out_mat, file = "out_mammals_BFs_b256.tsv", sep = "\t",
             quote = FALSE, row.names = FALSE )

#--------------#
# PLOT RESULTS #
#--------------#
# Plot beta points against mean logL
pdf( file = "BpointsVSmlogL_DivModelMammals.pdf", paper = "a4r",
     width = 0, height = 0 )
par( mfrow = c(1,2) )
##\\ ACetal22-GBM
plot( ACetal22_GBM$b, ACetal22_GBM$mean.logl, pch = 19,
      col = rgb(0,0,0,alpha=0.3),
      xaxs = "i", xlim = c(0,1), xlab = "b", ylab = "mean logL",
      main = "Unconstrained | GBM",
      cex = 0.8, ylim = c(-2.5e07,0 ) )
lines( ACetal22_GBM$b, ACetal22_GBM$mean.logl )
##\\ BM23pUhb-GBP
plot( BM23pUhb_GBM$b, BM23pUhb_GBM$mean.logl, pch=19, 
      col = rgb(0,0,0,alpha=0.3),
      xaxs = "i", xlim = c(0,1), xlab = "b", ylab = "mean logL",
      main = "Post-K-Pg | GBM",
      cex = 0.8, ylim = c(-2.5e07,0 ) )
lines( BM23pUhb_GBM$b, BM23pUhb_GBM$mean.logl )
dev.off()
