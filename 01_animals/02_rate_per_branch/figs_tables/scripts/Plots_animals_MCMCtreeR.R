#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#--------------#
# LOAD PACKAGE #
#--------------#
library( MCMCtreeR )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir
path_to_file <- getActiveDocumentContext()$path
wd <- paste( dirname( path_to_file ), "/", sep = "" )
setwd( wd )
main_wd            <- gsub( pattern = "figs_tables/scripts/",
                            replacement = "", x = wd )
wd_dRetal15pUhbGBM <- gsub( pattern = "figs_tables/scripts/",
                            replacement = "00_MCMCtree/sum_analyses/01_posterior/mcmc_files_dRetal15pUhb_GBM", x = wd )
wd_dRetal15pUhbILN <- gsub( pattern = "figs_tables/scripts/",
                            replacement = "00_MCMCtree/sum_analyses/01_posterior/mcmc_files_dRetal15pUhb_ILN", x = wd )
wd_BM23pUhbGBM     <- gsub( pattern = "figs_tables/scripts/",
                            replacement = "00_MCMCtree/sum_analyses/01_posterior/mcmc_files_BM23pUhb_GBM", x = wd )
wd_BM23pUhbILN     <- gsub( pattern = "figs_tables/scripts/",
                            replacement = "00_MCMCtree/sum_analyses/01_posterior/mcmc_files_BM23pUhb_ILN", x = wd )
all_wds <- c( wd_dRetal15pUhbGBM, wd_dRetal15pUhbILN,
              wd_BM23pUhbGBM, wd_BM23pUhbILN )
name_dirs <- c( "GBM_Unconstrained", "ILN_Unconstrained",
                "GBM_UpperEdiacaran", "ILN_UpperEdiacaran" )
source( "../../../../src/Functions_plots_MCMCtreeR.R" )
source( "../../../../src/Functions_plots.R" )

#### LOAD DATA ----
#-------------#
# START TASKS #
#-------------#
# Create objects for each dataset that will be required for plotting
plot_obj <- create_plot_obj( abs_path = all_wds, name_entries = name_dirs,
                             tree_pattern = "95HPD\\.tree" )
# Create dir to save RData object
if( ! dir.exists( paste( main_wd, "figs_tables/out_RData", sep = "" ) ) ){
  dir.create( paste( main_wd, "figs_tables/out_RData", sep = "" ) )
}
# Save object, then comment the next line and just load the RData file!
# save( plot_obj, file = paste( main_wd,
#                               "figs_tables/out_RData/plot_obj.RData",
#                               sep = "" ) )
load( file = paste( main_wd, "figs_tables/out_RData/plot_obj.RData",
                    sep = "" ) )

#### PLOTS ----
#----------------#
# START PLOTTING #
#----------------#

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
# 0. GBM | Unconstrained VS Late Ediacaran  #
#----------------------------------------------------------#
# Saved in PDF: plots/Animals_CompareGBMmodels.pdf
nodes_2plot <- c( 55, # Metazoa
                  63, # Bilateria
                  64, # Deuterostomia
                  82  # Protostomia
)
##> GBM | Unconstrained
ind_mat <- which( names( plot_obj$node_ages$GBM_Unconstrained ) %in% nodes_2plot )
transp.col <- adjustcolor( col = "blue", alpha.f = 0.2 )
last.plot  <- mcmc.tree.plot.RETPLOT( phy = plot_obj$phy$GBM_Unconstrained, 
                                      xlim.scale = c(-20,1000),  # wide
                                      #xlim.scale = c(-20,2000),  # less wide
                                      #xlim.scale = c(-20,3000), # narrow
                                      #xlim.scale = c(-20,5000),  # very narrow
                                      node.ages = plot_obj$node_ages$GBM_Unconstrained[ind_mat], 
                                      show.tip.label = TRUE,
                                      analysis.type = "user", cex.tips = 1,
                                      time.correction = 100,
                                      scale.res = c( "Eon", "Period" ),
                                      plot.type = "distributions",
                                      cex.age = 0.6, cex.labels = 0.8,
                                      relative.height = 0.08, 
                                      col.tree = "grey40", no.margin = TRUE,
                                      add.time.scale = TRUE,
                                      grey.bars = FALSE,
                                      density.col = transp.col,
                                      density.border.col = "blue"
)

# Add now the other distribution with this helper function
col.mod = "grey"
transparency = 0.3 
col.plot <- c( "white", "black", "darkgreen", "red", "darkblue", "orange",
               "pink", "lightblue", "cyan", "brown", "purple", "lightgreen" )
##> GBM | Late Ediacaran
transp.col2  <- adjustcolor( col = col.plot[ 7 ], alpha.f = 0.2 )
coords_2plot <- add.extra.dists( phy = plot_obj$phy$GBM_UpperEdiacaran, num.models = 1,
                                 last.plot = last.plot,
                                 node.ages = plot_obj$node_ages$GBM_UpperEdiacaran[ind_mat],
                                 plot.type = "distributions",
                                 time.correction = 100,
                                 density.col = transp.col2,
                                 density.border.col = col.plot[ 10 ],
                                 distribution.height = 0.8,
                                 transparency = transparency,
                                 return_coords = TRUE )

# Add legend
legend( locator(1), legend = c( "GBM, Unconstrained",
                                "GBM, Late Ediacaran"),
        col = c( "blue", col.plot[c(10)] ),
        lty = 1,
        #bty = "n"
        bty = "l"
)

# Get calib names underneath the dists
names( nodes_2plot ) <- c( "Metazoa", "Bilateria",
                           "Deuterostomia", "Protostomia" )
sort_nodes2plot <- sort( nodes_2plot )
for( i in 1:length(nodes_2plot) ){
  text( mean(coords_2plot$coords_k[[ i ]]$x) -30,
        # Bring y axis a bit down
        min(coords_2plot$coords_k[[ i ]]$y) -1,
        names( sort_nodes2plot )[i], cex = 0.8 )    
}


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
# 0. ILN | Unconstrained VS Late Ediacaran               #
#--------------------------------------------------------#
# Saved in PDF: plots/Animals_CompareILNmodels.pdf
nodes_2plot <- c( 55, # Metazoa
                  63, # Bilateria
                  64, # Deuterostomia
                  82  # Protostomia
)
##> ILN | Unconstrained
ind_mat <- which( names( plot_obj$node_ages$ILN_Unconstrained ) %in% nodes_2plot )
transp.col <- adjustcolor( col = "blue", alpha.f = 0.2 )
last.plot  <- mcmc.tree.plot.RETPLOT( phy = plot_obj$phy$ILN_Unconstrained, 
                                      xlim.scale = c(-20,1000),  # wide
                                      #xlim.scale = c(-20,2000),  # less wide
                                      #xlim.scale = c(-20,3000), # narrow
                                      #xlim.scale = c(-20,5000),  # very narrow
                                      node.ages = plot_obj$node_ages$ILN_Unconstrained[ind_mat], 
                                      show.tip.label = TRUE,
                                      analysis.type = "user", cex.tips = 1,
                                      time.correction = 100,
                                      scale.res = c( "Eon", "Period" ),
                                      plot.type = "distributions",
                                      cex.age = 0.6, cex.labels = 0.8,
                                      relative.height = 0.08, 
                                      col.tree = "grey40", no.margin = TRUE,
                                      add.time.scale = TRUE,
                                      grey.bars = FALSE,
                                      density.col = transp.col,
                                      density.border.col = "blue"
)

# Add now the other distribution with this helper function
col.mod = "grey"
transparency = 0.3 
col.plot <- c( "white", "black", "darkgreen", "red", "darkblue", "orange",
               "pink", "lightblue", "cyan", "brown", "purple", "lightgreen" )
##> ILN | Late Ediacaran
transp.col2  <- adjustcolor( col = col.plot[ 7 ], alpha.f = 0.2 )
coords_2plot <- add.extra.dists( phy = plot_obj$phy$ILN_UpperEdiacaran, num.models = 1,
                                 last.plot = last.plot,
                                 node.ages = plot_obj$node_ages$ILN_UpperEdiacaran[ind_mat],
                                 plot.type = "distributions",
                                 time.correction = 100,
                                 density.col = transp.col2,
                                 density.border.col = col.plot[ 10 ],
                                 distribution.height = 0.8,
                                 transparency = transparency,
                                 return_coords = TRUE )

# Add legend
legend( locator(1), legend = c( "ILN, Unconstrained",
                                "ILN, Late Ediacaran"),
        col = c( "blue", col.plot[c(10)] ),
        lty = 1,
        #bty = "n"
        bty = "l"
)

# Get calib names underneath the dists
names( nodes_2plot ) <- c( "Metazoa", "Bilateria",
                           "Deuterostomia", "Protostomia" )
sort_nodes2plot <- sort( nodes_2plot )
for( i in 1:length(nodes_2plot) ){
  text( mean(coords_2plot$coords_k[[ i ]]$x) -30,
        # Bring y axis a bit down
        min(coords_2plot$coords_k[[ i ]]$y) -1,
        names( sort_nodes2plot )[i], cex = 0.8 )    
}
