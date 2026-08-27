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
main_wd        <- gsub( pattern = "figs_tables/scripts/",
                        replacement = "", x = wd )
wd_ACetal22GBM <- gsub( pattern = "figs_tables/scripts/",
                        replacement = "00_MCMCtree/sum_analyses/01_posterior/mcmc_files_ACetal22_GBM", x = wd )
wd_ACetal22ILN <- gsub( pattern = "figs_tables/scripts/",
                        replacement = "00_MCMCtree/sum_analyses/01_posterior/mcmc_files_ACetal22_ILN", x = wd )
wd_BM23pUhbGBM <- gsub( pattern = "figs_tables/scripts/",
                        replacement = "00_MCMCtree/sum_analyses/01_posterior/mcmc_files_BM23pUhb_GBM", x = wd )
wd_BM23pUhbILN <- gsub( pattern = "figs_tables/scripts/",
                        replacement = "00_MCMCtree/sum_analyses/01_posterior/mcmc_files_BM23pUhb_ILN", x = wd )
all_wds <- c( wd_ACetal22GBM, wd_ACetal22ILN,
              wd_BM23pUhbGBM, wd_BM23pUhbILN )
name_dirs <- c( "GBM_Unconstrained", "ILN_Unconstrained",
                "GBM_PostKPg", "ILN_PostKPg" )
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
# Save object, then comment the next line and just load the RData file
# the next time you want to re-run this script!
# save( plot_obj, file = paste( main_wd,
#                               "figs_tables/out_RData/plot_obj.RData",
#                               sep = "" ) )
load( file = paste( main_wd, "figs_tables/out_RData/plot_obj.RData",
                    sep = "" ) )

#### PLOTS ----
#----------------#
# START PLOTTING #
#----------------#

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
# 0. GBM | Unconstrained VS Post-K-Pg #
#--------------------------------------------#
# Saved in PDF: plots/Mammals_CompareGBMmodels.pdf
nodes_2plot <- c( 73, # Mammalia
                  74, # Theria
                  75, # Marsupialia
                  77, # Placentalia
                  78, # Atlantogenata
                  82 # Boreoeutheria
)
##> GBM | Unconstrained
ind_mat <- which( names( plot_obj$node_ages$GBM_Unconstrained ) %in% nodes_2plot )
transp.col <- adjustcolor( col = "blue", alpha.f = 0.2 )
last.plot  <- mcmc.tree.plot.RETPLOT( phy = plot_obj$phy$GBM_Unconstrained, 
                                      xlim.scale = c(-80,300),  # wide
                                      #xlim.scale = c(-50,600),  # less wide
                                      #xlim.scale = c(-50,700), # narrow
                                      #xlim.scale = c(-60,700),  # very narrow
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
##> GBM | Post-K-Pg
transp.col2  <- adjustcolor( col = col.plot[ 10 ], alpha.f = 0.2 )
coords_2plot <- add.extra.dists( phy = plot_obj$phy$GBM_PostKPg, num.models = 1,
                                 last.plot = last.plot,
                                 node.ages = plot_obj$node_ages$GBM_PostKPg[ind_mat],
                                 plot.type = "distributions",
                                 time.correction = 100,
                                 density.col = transp.col2,
                                 density.border.col = col.plot[ 10 ],
                                 distribution.height = 0.8,
                                 transparency = transparency,
                                 return_coords = TRUE )

# Add legend
legend( locator(1), legend = c( "GBM, Unconstrained",
                                "GBM, Post-K-Pg"),
        col = c( "blue", col.plot[c(10)] ),
        lty = 1,
        #bty = "n"
        bty = "l"
)

# Get calib names underneath the dists
names( nodes_2plot ) <- c( "Mammalia", "Theria", "Marsupialia", "Placentalia",
                           "Atlantogenata", "Boreoeutheria" )
sort_nodes2plot <- sort( nodes_2plot )
for( i in 1:length(nodes_2plot) ){
  text( mean(coords_2plot$coords_k[[ i ]]$x) - 20,
        # Bring y axis a bit down
        min(coords_2plot$coords_k[[ i ]]$y) - 1 ,
        names( sort_nodes2plot )[i], cex = 0.8 )    
}

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
# 0. ILN | Unconstrained VS Post-K-Pg #
#--------------------------------------------#
# Saved in PDF: plots/Mammals_CompareILNmodels.pdf
nodes_2plot <- c( 73, # Mammalia
                  74, # Theria
                  75, # Marsupialia
                  77, # Placentalia
                  78, # Atlantogenata
                  82 # Boreoeutheria
)
##> ILN | Unconstrained
ind_mat <- which( names( plot_obj$node_ages$ILN_Unconstrained ) %in% nodes_2plot )
transp.col <- adjustcolor( col = "blue", alpha.f = 0.2 )
last.plot  <- mcmc.tree.plot.RETPLOT( phy = plot_obj$phy$ILN_Unconstrained, 
                                      xlim.scale = c(-80,300),  # wide
                                      #xlim.scale = c(-50,600),  # less wide
                                      #xlim.scale = c(-50,700), # narrow
                                      #xlim.scale = c(-60,700),  # very narrow
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
##> ILN | Post-K-Pg
transp.col2  <- adjustcolor( col = col.plot[ 10 ], alpha.f = 0.2 )
coords_2plot <- add.extra.dists( phy = plot_obj$phy$ILN_PostKPg, num.models = 1,
                                 last.plot = last.plot,
                                 node.ages = plot_obj$node_ages$ILN_PostKPg[ind_mat],
                                 plot.type = "distributions",
                                 time.correction = 100,
                                 density.col = transp.col2,
                                 density.border.col = col.plot[ 10 ],
                                 distribution.height = 0.8,
                                 transparency = transparency,
                                 return_coords = TRUE )

# Add legend
legend( locator(1), legend = c( "ILN, Unconstrained",
                                "ILN, Post-K-Pg"),
        col = c( "blue", col.plot[c(10)] ),
        lty = 1,
        #bty = "n"
        bty = "l"
)

# Get calib names underneath the dists
names( nodes_2plot ) <- c( "Mammalia", "Theria", "Marsupialia", "Placentalia",
                           "Atlantogenata", "Boreoeutheria" )
sort_nodes2plot <- sort( nodes_2plot )
for( i in 1:length(nodes_2plot) ){
  text( mean(coords_2plot$coords_k[[ i ]]$x) - 25,
        # Bring y axis a bit down
        min(coords_2plot$coords_k[[ i ]]$y) - 1 ,
        names( sort_nodes2plot )[i], cex = 0.8 )    
}
