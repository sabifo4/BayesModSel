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
wd <- gsub( pattern = "scripts/", replacement = "", x = wd )
setwd( wd )

#----------------#
# LOAD LIBRARIES #
#----------------#
# Install if needed by uncommenting line below
# devtools::install_github("dosreislab/mcmc3r", , build_vignettes = TRUE )
library( mcmc3r )

#-------------------------------------------------#
# RUN mcmc3r TO CREATE ARCHITECTURE AND CTL FILES #
#-------------------------------------------------#
# Check file architecture is OK 
if( ! dir.exists( "MCMCtree" ) ){
  dir.create( "MCMCtree" )
}
# Run functions to get `beta` vals and generate `MCMCtree`
# control files 
models <- c( "ACetal22_GBM", "BM23pUhb_GBM",
             "ACetal22_ILN", "BM23pUhb_ILN"
            
             )
trees  <- c( "72sp_ACetal22.tree", "72sp_BM23pUhb.tree",
             "72sp_ACetal22.tree", "72sp_BM23pUhb.tree"
             )
b_points <- 256
# Set main wd
setwd( wd )
for( j in 1:length(models) ){
  if( ! dir.exists( paste( wd, "MCMCtree/", models[j], sep = "" ) ) ){
    dir.create( paste( wd, "MCMCtree/", models[j], sep = "" ) )
  }
  # Set current wd where the file will be created
  setwd( paste( wd, "MCMCtree/", models[j], sep = "" ) )
  # Get `beta` values using the SS approach
  b <- mcmc3r::make.beta( n = b_points, method = c( "step-stones" ), a = 5 )
  # Use `beta` values saved in object `b` to generate `MCMCtree` control file
  mcmc3r::make.bfctlf( b, ctlf = paste( "mcmctree_", models[j], ".ctl",
                                        sep = "" ),
                       betaf = "betaweights.txt" )
  # Reset main wd for next iteration
  setwd( wd )
}

# Update the control files generated in the previous loop so that they contain
# the rest of the settings to run `MCMCtree` so that it samples from the power
# posterior
clock_mods <- c( 3, 3, 2, 2 # ACetal22-GBM, BM23pUhb-GBM,
                            # ACetal22-ILN, BM23pUhb-ILN
                )
names( clock_mods ) <- models
for( j in 1:length(models) ){
  cat( "Updating files for model ", models[j], "...\n" )
  for( k in 1:b_points ){
    cat( "Beta ", k, "\n" )
    # Copy this in the control file for beta "k" 
    beta_val <- readLines( con = paste( "MCMCtree/", models[j], "/", k,
                                        "/mcmctree_", 
                                        models[j], ".ctl", sep = "" ) )
    tmp_ctl  <- readLines( con = "../../01_modsel_workflow/00_mcmc3r/control_template/mcmctree_tmp.ctl" )
    tmp_ctl  <- gsub( pattern = "treefile = TREE",
                      replacement = paste( "treefile = ", trees[j], sep ="" ),
                      x = tmp_ctl )
    tmp_ctl  <- gsub( pattern = "clock = RCLOCK",
                      replacement = paste( "clock = ", clock_mods[j], sep ="" ),
                      x = tmp_ctl )
    # Generate an empty character vector to fill out with the same content
    # and then add the BayesFactorBeta at the end
    out_ctl  <- vector( mode = "character", length = length( tmp_ctl )+1 )
    out_ctl[1:length(tmp_ctl)] <- tmp_ctl 
    out_ctl[length(out_ctl)]   <- beta_val
    # Now, write everything and delete the old file
    writeLines( text = out_ctl, 
                con = paste( "MCMCtree/", models[j], "/", k, "/mcmctree_",
                             models[j], "_b", k, ".ctl", sep = "" ) )
    unlink( x = paste( "MCMCtree/", models[j], "/", k, "/mcmctree_", 
                       models[j], ".ctl", sep = "" ) )
  }
}


#-------------------#
# RUN MCMCTREE NOW! #
#-------------------#

