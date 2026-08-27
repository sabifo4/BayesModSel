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
# NOTE: This function will create a directory called `plots` inside the 
# main working directory if you have not created this dir yet
home_dir <- set_workdir()$home
# By now, set the working directory to `home_dir`
setwd( home_dir )

#-------------------------------------------------------------#
# DEFINE GLOBAL VARIABLES -- modify according to your dataset #
#-------------------------------------------------------------#
# First, we will define the global variables that match the settings in our 
# analysis.

# Specify number of partitions
num_parts  <- 4

# Specify objects for dataset names
dat_prior  <- c( "Unconstrained-NODAT-GBM", "Unconstrained-NODAT-ILN",
                 "PostKPg-NODAT-GBM", "PostKPg-NODAT-ILN" )
dirs_prior <- c( "ACetal22_NODAT_GBM", "ACetal22_NODAT_ILN",
                 "BM23pUhb_NODAT_GBM", "BM23pUhb_NODAT_ILN" )
dat_post   <- c( "Unconstrained-GBM", "Unconstrained-ILN",
                 "PostKPg-GBM", "PostKPg-ILN" )
dirs_post  <- c( "ACetal22_GBM", "ACetal22_ILN",
                 "BM23pUhb_GBM", "BM23pUhb_ILN" )

#--------------#
# ANALYSE DATA #
#--------------#
# Create output dir for rategrams
plotdir_rg <-  paste( home_dir, "/plots/rategrams/", sep = "" ) 
if( ! dir.exists( plotdir_rg ) ){
  dir.create( plotdir_rg )
}

# Extract rategrams for each dataset, partition, and chain
##\\ [[ PRIOR ]]
tmp_edge_vec <- vector( mode = "list", length( num_parts ) )
for( i in 1:length( dirs_prior ) ){
  for( j in 1:num_parts ){ # for each partition, get all rategrams
    tmp_path_tt  <- paste( "00_prior/", dirs_prior[i], "/", j, "/",
                           sep = "" )
    tmp_brf      <- list.files( path = tmp_path_tt )
    cat( paste( "Plotting rategram for dataset ", dirs_prior[i],
                " | Partition: ", j, " | Filtered chains: ",
                length( tmp_brf ), "\n", sep = "" ) )
    tmp_tt_edge  <- vector( mode = "list", length( tmp_brf ) )
    # For each chain that passed the filters, there are k files with rategrams
    for( k in 1:length( tmp_brf ) ){
      if( k == 1 ){
        # Save one tree object phylo
        tmp_tt_phylo <- ape::read.tree( file = paste( tmp_path_tt,
                                                      tmp_brf[k], sep = "" ) )
      }
      # Read k rategram
      tmp_list_tt <- ape::read.tree( file = paste( tmp_path_tt,
                                                   tmp_brf[k], sep = "" ) )
      # Save only edge length from k rategram (i.e., est branch rates)
      tmp_tt_edge[[ k ]] <- tmp_list_tt$edge.length
    }
    # Convert list into a matrix/array; then get the mean for rows (1)
    tmp_edge_mean <- apply( simplify2array( tmp_tt_edge ), 1, mean )
    # Replace branch lengths with the mean calculated above and save in tmp
    # vector to then plot the rategram with the average rate across all 
    # partitions
    tmp_tt_phylo$edge.length <- tmp_edge_vec[[ j ]] <- tmp_edge_mean
    # Output rategram with estimated mean branch rates
    ape::write.tree( phy = tmp_tt_phylo, file = paste( "00_prior/",
                                                       dirs_prior[i], "/",
                                                       "brates_tree_",
                                                       dirs_prior[i], "_p", j,
                                                       ".tree", sep = "" ) )
    # Find MRCA for Metazoa, Bilateria, and Deuterostomia
    # Do this only for the first partition, then reuse
    if( j == 1 ){
      tips_Theria        <- c( "homo_sapiens", "monodelphis_domestica" )
      tips_Placentalia   <- c( "homo_sapiens", "choloepus_hoffmanni" )
      tips_Marsupialia   <- c( "sarcophilus_harrisii", "monodelphis_domestica" )
      tips_Boreoeutheria <- c( "homo_sapiens", "sorex_araneus" )
      tips_Atlantogenata <- c( "loxodonta_africana", "choloepus_hoffmanni" )
      mrca_node_Theria        <- phytools::findMRCA( tmp_tt_phylo, tips_Theria )
      mrca_node_Placentalia   <- phytools::findMRCA( tmp_tt_phylo, tips_Placentalia )
      mrca_node_Marsupialia   <- phytools::findMRCA( tmp_tt_phylo, tips_Marsupialia )
      mrca_node_Boreoeutheria <- phytools::findMRCA( tmp_tt_phylo, tips_Boreoeutheria )
      mrca_node_Atlantogenata <- phytools::findMRCA( tmp_tt_phylo, tips_Atlantogenata )
      all_mrca_nodes <- c( mrca_node_Theria, mrca_node_Placentalia, 
                           mrca_node_Marsupialia, mrca_node_Boreoeutheria,
                           mrca_node_Atlantogenata )
    }
    # Find estimated branch rates
    br_Theria        <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Theria ) ]
    br_Placentalia   <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Placentalia ) ]
    br_Marsupialia   <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Marsupialia ) ]
    br_Boreoeutheria <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Boreoeutheria ) ]
    br_Atlantogenata <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Atlantogenata ) ]
    # Plot rategram
    if( ! dir.exists( paste( plotdir_rg, dirs_prior[i], sep = "" ) ) ){
      dir.create( paste( plotdir_rg, dirs_prior[i], sep = "" ) )
    }
    pdf( paste( plotdir_rg, dirs_prior[i], "/", dirs_prior[i], "_rategram_p",
                j, ".pdf", sep = "" ), paper = "a4r" )
    phytools::plotBranchbyTrait( tree = tmp_tt_phylo,
                                 x = tmp_tt_phylo$edge.length,
                                 mode = "edges",
                                 palette = "rainbow", cex = 0.5,
                                 edge.width = 0.5 )
    # Add node labels to the tree
    ape::nodelabels( text = paste( c( "Theria\n", "Placentalia\n",
                                      "Marsupialia\n", "Boreoeutheria\n",
                                      "Atlantogenata\n" ),
                                   c( round( br_Theria, 4 ),
                                      round( br_Placentalia, 4 ),
                                      round( br_Marsupialia, 4 ),
                                      round( br_Boreoeutheria, 4 ),
                                      round( br_Atlantogenata, 4 ) ),
                                   sep = "" ),
                     node = all_mrca_nodes, adj = c(1.1, -0.3), frame = "none",
                     cex = 0.6 )
    # Add the rest of the nodes in a smaller font size
    # rm_nodes       <- which( tmp_tt_phylo$edge[,2] %in% all_mrca_nodes )
    # small_br_2plot <- tmp_tt_phylo$edge.length[ tmp_tt_phylo$edge[,2][-rm_nodes] ]
    # ape::nodelabels( text = round( small_br_2plot, 2 ),
    #                  node = tmp_tt_phylo$edge[,2][-rm_nodes],
    #                  adj = c(1.1, -0.3), frame = "none",
    #                  cex = 0.35 )
    # Close file
    dev.off()
    
  }
  
  # Plot rategram with avg rate across all partitions
  tmp_edge_parts_mean      <- apply( simplify2array( tmp_edge_vec ), 1, mean )
  tmp_tt_phylo$edge.length <- tmp_edge_parts_mean
  cat( "Calculating average branch rates across all partitions...\n\n")
  # Output rategram with estimated mean branch rates
  ape::write.tree( phy = tmp_tt_phylo, file = paste( "00_prior/",
                                                     dirs_prior[i], "/",
                                                     "meanbrates_tree_",
                                                     dirs_prior[i],
                                                     ".tree", sep = "" ) )
  # Extract mean brates
  br_Theria        <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Theria ) ]
  br_Placentalia   <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Placentalia ) ]
  br_Marsupialia   <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Marsupialia ) ]
  br_Boreoeutheria <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Boreoeutheria ) ]
  br_Atlantogenata <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Atlantogenata ) ]
  # Plot rategram
  pdf( paste( plotdir_rg, dirs_prior[i], "/", dirs_prior[i], "_rategram_mean.pdf",
              sep = "" ), paper = "a4r" )
  phytools::plotBranchbyTrait( tree = tmp_tt_phylo,
                               x = tmp_tt_phylo$edge.length,
                               mode = "edges",
                               palette = "rainbow", cex = 0.5,
                               edge.width = 0.5 )
  # Add node labels to the tree
  ape::nodelabels( text = paste( c( "Theria\n", "Placentalia\n",
                                    "Marsupialia\n", "Boreoeutheria\n",
                                    "Atlantogenata\n" ),
                                 c( round( br_Theria, 4 ),
                                    round( br_Placentalia, 4 ),
                                    round( br_Marsupialia, 4 ),
                                    round( br_Boreoeutheria, 4 ),
                                    round( br_Atlantogenata, 4 ) ),
                                 sep = "" ),
                   node = all_mrca_nodes, adj = c(1.1, -0.3), frame = "none",
                   cex = 0.6 )
  # Add the rest of the nodes in a smaller font size
  # rm_nodes       <- which( tmp_tt_phylo$edge[,2] %in% all_mrca_nodes )
  # small_br_2plot <- tmp_tt_phylo$edge.length[ tmp_tt_phylo$edge[,2][-rm_nodes] ]
  # ape::nodelabels( text = round( small_br_2plot, 2 ),
  #                  node = tmp_tt_phylo$edge[,2][-rm_nodes],
  #                  adj = c(1.1, -0.3), frame = "none",
  #                  cex = 0.35 )
  # Close file
  dev.off()
  
}


##\\ [[ POSTERIOR ]]
tmp_edge_vec <- vector( mode = "list", length( num_parts ) )
for( i in 1:length( dirs_post ) ){
  for( j in 1:num_parts ){ # for each partition, get all rategrams
    tmp_path_tt  <- paste( "01_posterior/", dirs_post[i], "/", j, "/",
                           sep = "" )
    tmp_brf      <- list.files( path = tmp_path_tt )
    cat( paste( "Plotting rategram for dataset ", dirs_post[i],
                " | Partition: ", j, " | Filtered chains: ",
                length( tmp_brf ), "\n", sep = "" ) )
    tmp_tt_edge  <- vector( mode = "list", length( tmp_brf ) )
    # For each chain that passed the filters, there are k files with rategrams
    for( k in 1:length( tmp_brf ) ){
      if( k == 1 ){
        # Save one tree object phylo
        tmp_tt_phylo <- ape::read.tree( file = paste( tmp_path_tt,
                                                      tmp_brf[k], sep = "" ) )
      }
      # Read k rategram
      tmp_list_tt <- ape::read.tree( file = paste( tmp_path_tt,
                                                   tmp_brf[k], sep = "" ) )
      # Save only edge length from k rategram (i.e., est branch rates)
      tmp_tt_edge[[ k ]] <- tmp_list_tt$edge.length
    }
    # Convert list into a matrix/array; then get the mean for rows (1)
    tmp_edge_mean <- apply( simplify2array( tmp_tt_edge ), 1, mean )
    # Replace branch lengths with the mean calculated above and save in tmp
    # vector to then plot the rategram with the average rate across all 
    # partitions
    tmp_tt_phylo$edge.length <- tmp_edge_vec[[ j ]] <- tmp_edge_mean
    # Output rategram with estimated mean branch rates
    ape::write.tree( phy = tmp_tt_phylo, file = paste( "01_posterior/",
                                                       dirs_post[i], "/",
                                                       "brates_tree_",
                                                       dirs_post[i], "_p", j,
                                                       ".tree", sep = "" ) )
    # Find MRCA for Metazoa, Bilateria, and Deuterostomia
    # Do this only for the first partition, then reuse
    if( j == 1 ){
      tips_Theria        <- c( "homo_sapiens", "monodelphis_domestica" )
      tips_Placentalia   <- c( "homo_sapiens", "choloepus_hoffmanni" )
      tips_Marsupialia   <- c( "sarcophilus_harrisii", "monodelphis_domestica" )
      tips_Boreoeutheria <- c( "homo_sapiens", "sorex_araneus" )
      tips_Atlantogenata <- c( "loxodonta_africana", "choloepus_hoffmanni" )
      mrca_node_Theria        <- phytools::findMRCA( tmp_tt_phylo, tips_Theria )
      mrca_node_Placentalia   <- phytools::findMRCA( tmp_tt_phylo, tips_Placentalia )
      mrca_node_Marsupialia   <- phytools::findMRCA( tmp_tt_phylo, tips_Marsupialia )
      mrca_node_Boreoeutheria <- phytools::findMRCA( tmp_tt_phylo, tips_Boreoeutheria )
      mrca_node_Atlantogenata <- phytools::findMRCA( tmp_tt_phylo, tips_Atlantogenata )
      all_mrca_nodes <- c( mrca_node_Theria, mrca_node_Placentalia, 
                           mrca_node_Marsupialia, mrca_node_Boreoeutheria,
                           mrca_node_Atlantogenata )
    }
    # Find estimated branch rates
    br_Theria        <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Theria ) ]
    br_Placentalia   <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Placentalia ) ]
    br_Marsupialia   <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Marsupialia ) ]
    br_Boreoeutheria <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Boreoeutheria ) ]
    br_Atlantogenata <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Atlantogenata ) ]
    # Plot rategram
    if( ! dir.exists( paste( plotdir_rg, dirs_post[i], sep = "" ) ) ){
      dir.create( paste( plotdir_rg, dirs_post[i], sep = "" ) )
    }
    pdf( paste( plotdir_rg, dirs_post[i], "/", dirs_post[i], "_rategram_p",
                j, ".pdf", sep = "" ), paper = "a4r" )
    phytools::plotBranchbyTrait( tree = tmp_tt_phylo,
                                 x = tmp_tt_phylo$edge.length,
                                 mode = "edges",
                                 palette = "rainbow", cex = 0.5,
                                 edge.width = 0.5 )
    # Add node labels to the tree
    ape::nodelabels( text = paste( c( "Theria\n", "Placentalia\n",
                                      "Marsupialia\n", "Boreoeutheria\n",
                                      "Atlantogenata\n" ),
                                   c( round( br_Theria, 4 ),
                                      round( br_Placentalia, 4 ),
                                      round( br_Marsupialia, 4 ),
                                      round( br_Boreoeutheria, 4 ),
                                      round( br_Atlantogenata, 4 ) ),
                                   sep = "" ),
                     node = all_mrca_nodes, adj = c(1.1, -0.3), frame = "none",
                     cex = 0.6 )
    # Add the rest of the nodes in a smaller font size
    # rm_nodes       <- which( tmp_tt_phylo$edge[,2] %in% all_mrca_nodes )
    # small_br_2plot <- tmp_tt_phylo$edge.length[ tmp_tt_phylo$edge[,2][-rm_nodes] ]
    # ape::nodelabels( text = round( small_br_2plot, 2 ),
    #                  node = tmp_tt_phylo$edge[,2][-rm_nodes],
    #                  adj = c(1.1, -0.3), frame = "none",
    #                  cex = 0.35 )
    # Close file
    dev.off()
    
  }
  
  # Plot rategram with avg rate across all partitions
  tmp_edge_parts_mean      <- apply( simplify2array( tmp_edge_vec ), 1, mean )
  tmp_tt_phylo$edge.length <- tmp_edge_parts_mean
  cat( "Calculating average branch rates across all partitions...\n\n")
  # Output rategram with estimated mean branch rates
  ape::write.tree( phy = tmp_tt_phylo, file = paste( "01_posterior/",
                                                     dirs_post[i], "/",
                                                     "meanbrates_tree_",
                                                     dirs_post[i],
                                                     ".tree", sep = "" ) )
  # Extract mean brates
  br_Theria        <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Theria ) ]
  br_Placentalia   <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Placentalia ) ]
  br_Marsupialia   <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Marsupialia ) ]
  br_Boreoeutheria <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Boreoeutheria ) ]
  br_Atlantogenata <- tmp_tt_phylo$edge.length[ which( tmp_tt_phylo$edge[,2] == mrca_node_Atlantogenata ) ]
  # Plot rategram
  pdf( paste( plotdir_rg, dirs_post[i], "/", dirs_post[i], "_rategram_mean.pdf",
              sep = "" ), paper = "a4r" )
  phytools::plotBranchbyTrait( tree = tmp_tt_phylo,
                               x = tmp_tt_phylo$edge.length,
                               mode = "edges",
                               palette = "rainbow", cex = 0.5,
                               edge.width = 0.5 )
  # Add node labels to the tree
  ape::nodelabels( text = paste( c( "Theria\n", "Placentalia\n",
                                    "Marsupialia\n", "Boreoeutheria\n",
                                    "Atlantogenata\n" ),
                                 c( round( br_Theria, 4 ),
                                    round( br_Placentalia, 4 ),
                                    round( br_Marsupialia, 4 ),
                                    round( br_Boreoeutheria, 4 ),
                                    round( br_Atlantogenata, 4 ) ),
                                 sep = "" ),
                   node = all_mrca_nodes, adj = c(1.1, -0.3), frame = "none",
                   cex = 0.6 )
  # Add the rest of the nodes in a smaller font size
  # rm_nodes       <- which( tmp_tt_phylo$edge[,2] %in% all_mrca_nodes )
  # small_br_2plot <- tmp_tt_phylo$edge.length[ tmp_tt_phylo$edge[,2][-rm_nodes] ]
  # ape::nodelabels( text = round( small_br_2plot, 2 ),
  #                  node = tmp_tt_phylo$edge[,2][-rm_nodes],
  #                  adj = c(1.1, -0.3), frame = "none",
  #                  cex = 0.35 )
  # Close file
  dev.off()
  
}
