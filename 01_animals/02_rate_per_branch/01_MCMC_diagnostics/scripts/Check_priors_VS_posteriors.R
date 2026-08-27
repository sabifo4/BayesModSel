#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------------------------------#
# LOAD PACKAGES, FUNCTIONS, AND SET ENVIRONMENT #
#-----------------------------------------------#
# This package lets you find automatically the path to a specific location
# in your file structure
# If you have not installed this package, you will need to install it. 
# You can uncomment the following line to do this:
#install.packages( "rstudioapi" )
library( rstudioapi )
scripts_dir   <- gsub( pattern = "scripts..*", replacement = "scripts/",
                       x = getActiveDocumentContext()$path )
setwd( scripts_dir )
# Load the file with all the functions used throughout this script
source( file = "../../../../src/Functions.R" )
# Run in-house function to set home directory and output directory for ESS
# and convergence tests
# NOTE: This function will create a directory called `plots` and another called
# `ESS_and_chains_convergence` inside the `analyses` directory if you have
# not created them yet
home_dir      <- set_homedir()$home
outchecks_dir <- set_homedir()$ESS
# By now, set the working directory to `home_dir`
setwd( home_dir )


#-------------------------------------------------------------#
# DEFINE GLOBAL VARIABLES -- modify according to your dataset #
#-------------------------------------------------------------#
# First, we will define global variables that we will keep using throughout this
# script.

# 1. Number of chains
num_chains <- 6

# 2. Number of divergence times that have been estimated. One trick to find
# this out quickly is to subtract 1 to the number of species. In this case,
# there are 54 taxa (54), so the number of internal nodes
# is `n_taxa-=54-1=53`.
# Another way to verify this is by opening the `mcmc.txt` file and check the
# header. The first element after `Gen` will have the format of `t_nX`, where
# X will be an integer (i.e., 55). Subtract two to this number 
# (i.e., 55-2=53) and this will be your number of divergence times that are 
# parameters of the MCMC. Please modify the number below so it fits to the 
# dataset you are using. 
num_divt <- 53

# 3. If you set `print= 2`, then you must specify the number of branches
# in the tree. You can calculate this by `2*num_sp-1`. In our case:
# 2*72-1=143
# You also need to specify the number of partitions, as there will be
# a set of "n" branch rates for each partition
num_brate <- 107
num_part  <- 2

# 4. Total number of samples that you collected after generating the
# final `mcmc.txt` files with those from the chains that passed the filters. 
# You can check these numbers in scripts `MCMC_diagnostics_posterior.R` and
# `MCMC_diagnostics_prior.R`. E.g., `sum_post_QC$<name_dataset>$total_samples`
# or `sum_prior_QC$<name_dataset>$total_samples`
#
# dRetal15: The number of lines is 80004 (NODAT-GBM), 100005 (NODAT-ILN),
#           100005 (GBM), and 120006 (ILN); you need to specify one less
# dRetal15pUhb: The number of lines is 80004 (NODAT-GBM), 100005 (NODAT-ILN),
#               60003 (GBM), and 120006 (ILN); you need to specify one less
# BM23: The number of lines is 80004 for both NODAT-GBM NODAT-ILN, 60003 (GBM),
#       and 100005 (ILN); you need to specify one less
# BM23pUhb: The number of lines is 40002 (NODAT-GBM), 120006 (NODAT-ILN),
#           80004 (GBM), and 100005 (ILN); you need to specify one less
#
# NOTE: If you had more than one dataset, you would add another vector of three
# values with the samples for CLK, GBM, and ILN to create `def_samples`
# E.g. two datasts: c( c( 120005, 120005, 120005), c( 120005, 120005, 120005) )
def_samples <- c( c( 80003, 100004, 60002, 120005 ),  # dRetal15pUhb | NDGBM, NDILN, GBM, ILN
                  c( 40001, 120005, 80003, 100004 )    # BM23pUhb | NDGBM, NDILN, GBM, ILN
)

# 5. Quantile percentage that you want to set By default, the variable below is 
# set to 0.975 so the 97.5% and 2.5% quantiles (i.e., 95%CI). If you want to
# change this, however, just modify the value.
perc <- 0.975

# 6. Number of columns in the `mcmc.txt` that are to be deleted as they do not 
# correspond to sample values for divergence times (i.e., the entries are not 
# names following the format `t_nX`). To figure out this number quickly, you 
# can open the `mcmc.txt` file, read the header, and count the number of `mu*`
# elements. Do not count the `lnL` value when looking at 
# `mcmc.txt` files generated when sampling from the posterior -- this is 
# automatically accounted for in the in-house R functions that you will 
# subsequently use. E.g., you expect to see as many `mu[0-9]` as alignment
# blocks you have in your sequence file! E.g., if you had two alignment blocks,
# you would speciy `delcol_prior <- 2`. Please modify the value/s below 
# (depending on having one or more datasets) according to the `mcmc.txt` file
# generated when sampling from the prior (`delcol_prior`)
##> NOTE: If you ran `MCMCtree` with `clock = 2` or `clock = 3` when
##> sampling from the prior, you will also need to count the `sigma2*`
##> columns! We ran `clock = 1` so that the analyses ran quicker, and thus
##> we only have `mu*` columns.
delcol_obj <- 4 # There are 2 partitions, and hence 4 columns (we ran
# analyses with `clock = 2` and `clock = 3`!
# mu[1-2] and sigma2[1-2]

# 7. Path to the directory where the concatenated `mcmc.txt` file has been 
# generated. Note that, if you have run more than one chain in `MCMCtree` for
# each hypothesis tested, you are expected to have generated a concatenated 
# `mcmc.txt` file with the bash script `Combine_MCMC_prior.sh` or any similar 
# approaches.
num_dirs     <- 8
num_datasets <- 2
paths_dat <- all_paths <- c( # dRetal15pUhb
                             paste( home_dir,
                                    "../00_MCMCtree/sum_analyses/00_prior/mcmc_files_dRetal15pUhb_NODAT_GBM/",
                                    sep = "" ),
                             paste( home_dir,
                                    "../00_MCMCtree/sum_analyses/00_prior/mcmc_files_dRetal15pUhb_NODAT_ILN/",
                                    sep = "" ),
                             paste( home_dir,
                                    "../00_MCMCtree/sum_analyses/01_posterior/mcmc_files_dRetal15pUhb_GBM/",
                                    sep = "" ),
                             paste( home_dir,
                                    "../00_MCMCtree/sum_analyses/01_posterior/mcmc_files_dRetal15pUhb_ILN/",
                                    sep = "" ),
                             # BM23pUhb
                             paste( home_dir,
                                    "../00_MCMCtree/sum_analyses/00_prior/mcmc_files_BM23pUhb_NODAT_GBM/",
                                    sep = "" ),
                             paste( home_dir,
                                    "../00_MCMCtree/sum_analyses/00_prior/mcmc_files_BM23pUhb_NODAT_ILN/",
                                    sep = "" ),
                             paste( home_dir,
                                    "../00_MCMCtree/sum_analyses/01_posterior/mcmc_files_BM23pUhb_GBM/",
                                    sep = "" ),
                             paste( home_dir,
                                    "../00_MCMCtree/sum_analyses/01_posterior/mcmc_files_BM23pUhb_ILN/",
                                    sep = "" )
)

# 8. Load a semicolon-separated file with info about calibrated nodes. Note that
# this file is output by script `Merge_node_labels.R`. A summary of its content
# in case you are to generate your own input files:
#
# Each column needs to be separated with semicolons and an extra blank line
# after the last row with calibration information needs to be added. If the
# extra blank is not added, R will complain and will not load the file!
# If you add a header, please make sure you name the column elements as 
# `Calib;node;Prior`. If not, the R function below will deal with the header,
# but make sure you set `head_avail = FALSE` when running `read_calib_f` 
# function below. An example of the content of this file is given below:
#
# ```
# Calib;node;Prior
# ex_n5;5;ST(5.8300,0.0590,0.1120,109.1240)
# ex_n7;7;B(4.1200,4.5200,0.0250,0.0250)
#
# ```
#
# The first column will have the name of the calibration/s that can help you
# identify which node belongs to which calibration. The second column is the
# number given to this node by`MCMCtree` (this information is automatically
# found when you run the script `Merge_node_labels.R`, otherwise you will need
# to keep checking the output file `node_trees.tree` to figure out which node
# is which). The third column is the calibration used for that node in
# `MCMCtree` format.
# 
# [[ NOTES ABOUT ALLOWED CALIBRATION FORMATS]]
#
# Soft-bound calibrations: 
#  E.g.1: A calibration with a minimum of 0.6 and a maximum of 0.8 would with  
#         the default tail probabilities would have the following equivalent 
#         formats:
#         >> B(0.6,0.8) | B(0.6,0.8,0.025,0.025)
#  E.g.2: A calibration with a minimum of 0.6 and a maximum of 0.8 would with  
#         the pL=0.001 and pU=0.025 would have the following format. Note that, 
#         whenever you want to modify either pL or pU, you need to write down 
#         the four  parameters in the format of "B(min,max,pL,pU)":
#         >> B(0.6,0.8,0.001,0.025)
#
# Lower-bound calibrations: 
#  E.g.1: A calibration with a minimum of 0.6 and the default parameters for
#         p = 0.1, c = 1, pL = 0.025:
#         >> L(0.6) | L(0.6,0.1,1,0.025)
#  E.g.2: A calibration with a hard minimum at 0.6, and so pL = 1e-300. 
#         Note that, whenever you want to modify either pL or pU, you need to  
#         write down the four parameters in the format of "L(min,p,c,pL)":
#         >> L(0.6,0.1,1,1e-300)
#
# Upper-bound calibrations: 
#  E.g.1: A calibration with a maximum of 0.8 and the default parameters for
#         pU = 0.025:
#         >> U(0.8) | U(0.8,0.025)
#  E.g.2: A calibration with a hard maximum at 0.8, and so pU = 1e-300. 
#         Note that, if you want to modify pU, you need to write down the two
#         parameters in the format of "U(max,pU)":
#         >> U(0.8,1e-300)
#
# ST distributions: 
#  The format accepted has four parameters: xi (location, mean root age), 
#  omega (scale), alpha (shape), nu (df). Accepted format: 
#  >> ST(5.8300,0.0590,0.1120,109.1240)
#
# SN distributions: 
#  The format accepted has three parameters: xi (location, mean root age), 
#  omega (scale), alpha (shape). Accepted format: 
#  >> SN(5.8300,0.0590,0.1120)  
#
#
# The next command executes the `read_calib_f` in-house function, which reads
# your input files (semicolon-separated files). The path to this directory is 
# what the argument `main_dir` needs. The argument `f_names` requires the name 
# of the file/s that you have used. Argument `dat` requires the same global 
# object that you have created at the beginning of the script.
dat <- c( # dRetal15pUhb
          "Unconstrained-NODAT-GBM", "Unconstrained-NODAT-ILN",
          "Unconstrained-GBM", "Unconstrained-ILN",
          # BM23pUhb
          "UpperEdiacaran-NODAT-GBM", "UpperEdiacaran-NODAT-ILN",
          "UpperEdiacaran-GBM", "UpperEdiacaran-ILN" )
calib_nodes <- read_calib_f( main_dir = paste( home_dir, "calibs/inp_calibs/",
                                               sep = "" ),
                             f_names = c( "Calibnodes_Unconstrained.csv",
                                          "Calibnodes_Unconstrained.csv",
                                          "Calibnodes_Unconstrained.csv",
                                          "Calibnodes_Unconstrained.csv",
                                          "Calibnodes_UpperEdiacaran.csv",
                                          "Calibnodes_UpperEdiacaran.csv",
                                          "Calibnodes_UpperEdiacaran.csv",
                                          "Calibnodes_UpperEdiacaran.csv" ),
                             dat = dat, head_avail = TRUE )


#-----------#
# LOAD DATA #
#-----------#
# Load mcmc files for each dataset
mcmc_dRetal15pUhb_obj <- mcmc_BM23pUhb_obj <- 
  vector( "list", num_dirs/num_datasets )
# dRetal15 - Unconstrained
prior <- c( TRUE, TRUE, FALSE, FALSE )
count <- 1
for( j in 1:4 ) {
  # dRetal15pUhb
  names( mcmc_dRetal15pUhb_obj )[j] <- dat[count]
  cat( "[[ Parsing file for dataset", names( mcmc_dRetal15pUhb_obj )[j],
       " ]]\n" )
  mcmc_dRetal15pUhb_obj[[j]] <- load_dat( mcmc = paste( paths_dat[count],
                                                        "/mcmc.txt", sep = "" ),
                                          delcol = delcol_obj, perc = perc,
                                          def_samples = def_samples[count],
                                          prior = prior[j] )
  count <- count + 4
  # BM23pUhb
  names( mcmc_BM23pUhb_obj )[j] <- dat[count]
  cat( "[[ Parsing file for dataset", names( mcmc_BM23pUhb_obj )[j], " ]]\n" )
  mcmc_BM23pUhb_obj[[j]] <- load_dat( mcmc = paste( paths_dat[count], 
                                                    "/mcmc.txt", sep = "" ),
                                      delcol = delcol_obj, perc = perc,
                                      def_samples = def_samples[count],
                                      prior = prior[j] )
  # Reset counter for the next round
  cat( "\n\n" )
  count <- j + 1
}

# Save datasets -- they are large files (~660Mb and ~618Mb each), so only run
# the commands below if you have enough space to save them
# Some of the figures you shall find in our main text and supplementary material
# are generated with scripts that load these objects so, if you want to 
# reproduce them, you will need to save these objects
if( ! dir.exists( paste( home_dir, "out_RData", sep = "" ) ) ){
  dir.create( paste( home_dir, "out_RData", sep = "" ) ) 
}
# save( file = paste( home_dir, "out_RData/mcmc_dRetal15pUhb_obj.RData",
#       sep = "" ), mcmc_dRetal15pUhb_obj )
# save( file = paste( home_dir, "out_RData/mcmc_BM23pUhb_obj.RData", sep = "" ),
#       mcmc_BM23pUhb_obj )
load( file = paste( home_dir, "out_RData/mcmc_dRetal15pUhb_obj.RData",
                    sep = "" ) )
load( file = paste( home_dir, "out_RData/mcmc_BM23pUhb_obj.RData", sep = "" ) )

#---------------------------#
# PLOTS: prior VS posterior #
#---------------------------#
##\\ UNCONSTRAINED
# Plot calibrated nodes
calib_nodes_dRetal15pUhb <- vector( mode = "list", length = 4 )
count <- 0
for( i in 1:4 ){
  count <- count + 1
  calib_nodes_dRetal15pUhb[[ count ]] <- calib_nodes[[ i ]]
  names( calib_nodes_dRetal15pUhb )[count] <- names( calib_nodes )[i]
}
plot_priorVSpost_v2( dat = "Unconstrained",
                     calib_nodes = calib_nodes_dRetal15pUhb,
                     mcmc_obj = mcmc_dRetal15pUhb_obj,
                     home_dir = home_dir,
                     multi_calib = TRUE,
                     # Colour order:
                     # Calib density, NODAT-GBM, GBM, NODAT-ILN, ILN
                     cols_vec = c( "black", "lightblue", "purple",
                                   "darkolivegreen3", "brown" ) )

##\\ UPPER EDIACARAN
# Plot calibrated nodes
calib_nodes_BM23pUhb <- vector( mode = "list", length = 4 )
count <- 0
for( i in 5:8 ){
  count <- count + 1
  calib_nodes_BM23pUhb[[ count ]] <- calib_nodes[[ i ]]
  names( calib_nodes_BM23pUhb )[count] <- names( calib_nodes )[i]
}
plot_priorVSpost_v2( dat = "UpperEdiacaran",
                     calib_nodes = calib_nodes_BM23pUhb,
                     mcmc_obj = mcmc_BM23pUhb_obj,
                     home_dir = home_dir,
                     multi_calib = TRUE,
                     # Colour order:
                     # Calib density, NODAT-GBM, GBM, NODAT-ILN, ILN
                     cols_vec = c( "black", "lightblue", "purple",
                                   "darkolivegreen3", "brown" ) )


#--------------------------------------------------#
# PLOTS: compare clocks under different strategies #
#--------------------------------------------------#
##\\ GBM | Unconstrained VS Late Ediacaran
# Define vector with user-specified calibrations
calibs_GBM <- matrix( 0, nrow = dim( calib_nodes_dRetal15pUhb$`Unconstrained-GBM` )[1],
                      ncol = 4 )
calibs_GBM[,1:3] <- as.matrix( calib_nodes_dRetal15pUhb$`Unconstrained-GBM` )
calibs_GBM[,4]   <- calib_nodes_BM23pUhb$`UpperEdiacaran-GBM`[,3]
# Resulting matrix, not list!
colnames( calibs_GBM ) <- c( "Name", "Node",
                             "dRetal15pUhb-Prior",
                             "BM23pUhb-Prior" )
# Define vector with all MCMC objects
mcmc_GBM        <- vector( mode = "list", 2 )
mcmc_GBM[[ 1 ]] <- mcmc_dRetal15pUhb_obj$`Unconstrained-GBM`
mcmc_GBM[[ 2 ]] <- mcmc_BM23pUhb_obj$`UpperEdiacaran-GBM`
names( mcmc_GBM ) <- c( "Unconstrained", "UpperEdiacaran" )
# Define colours
cols_vec <- c( "blue", "brown" )
# Plot
plot_compare_clocks( calib_nodes = calibs_GBM,
                     mcmc_obj = mcmc_GBM , home_dir = home_dir,
                     # Colour order:
                     # Calib density: Unconstrained, Late Ediacaran
                     cols_vec = cols_vec,
                     out_suffix = "all_GBM" )

##\\ ILN | Unconstrained VS Late Ediacaran
# Define vector with user-specified calibrations
calibs_ILN <- matrix( 0, nrow = dim( calib_nodes_dRetal15pUhb$`Unconstrained-ILN` )[1],
                      ncol = 4 )
calibs_ILN[,1:3] <- as.matrix( calib_nodes_dRetal15pUhb$`Unconstrained-ILN` )
calibs_ILN[,4]   <- calib_nodes_BM23pUhb$`UpperEdiacaran-ILN`[,3]
# Resulting matrix, not list!
colnames( calibs_ILN ) <- c( "Name", "Node",
                             "dRetal15pUhb-Prior",
                             "BM23pUhb-Prior" )
# Define vector with all MCMC objects
mcmc_ILN        <- vector( mode = "list", 2 )
mcmc_ILN[[ 1 ]] <- mcmc_dRetal15pUhb_obj$`Unconstrained-ILN`
mcmc_ILN[[ 2 ]] <- mcmc_BM23pUhb_obj$`UpperEdiacaran-ILN`
names( mcmc_ILN ) <- c( "Unconstrained", "UpperEdiacaran" )
# Define colours
cols_vec <- c( "blue", "brown" )
# Plot
plot_compare_clocks( calib_nodes = calibs_ILN,
                     mcmc_obj = mcmc_ILN , home_dir = home_dir,
                     # Colour order:
                     # Calib density: Unconstrained, Late Ediacaran
                     cols_vec = cols_vec,
                     out_suffix = "all_ILN" )

