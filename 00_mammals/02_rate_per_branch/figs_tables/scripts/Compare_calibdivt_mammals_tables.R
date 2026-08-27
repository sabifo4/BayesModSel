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
setwd( wd )
main_wd   <- gsub( pattern = "figs_tables/scripts/",
                   replacement = "", x = wd )
divcsv_wd <- gsub( pattern = "figs_tables/scripts/",
                   replacement = "01_MCMC_diagnostics/sum_files_post/sum_divtimes_tsv/", x = wd )
source( "../../../../src/Functions_plots.R" )

#### LOAD DATA ----
#-------------#
# START TASKS #
#-------------#
# Load all output csv tables previously generated
dir_names <- sort( c( "Post-K-Pg", "Unconstrained" ) )
all_divt  <- divt_csv( fpath = divcsv_wd, name_dirs = dir_names )
# Select nodes that are to be plotted:
nodes_calibs <- grep( x = rownames( all_divt$GBM$`PostKPg-GBM_FILT`),
                      pattern = "t_n[0-9]*_" )
nodes_2plot <- gsub( x = rownames( all_divt$GBM$`PostKPg-GBM_FILT`)[ nodes_calibs ],
                     pattern = "_[A-Z].*", replacement = "" )
plots_per_fig <- rep( 1, length( nodes_2plot ) )
names( nodes_2plot ) <- gsub( x = rownames( all_divt$GBM$`PostKPg-GBM_FILT`)[ nodes_calibs ],
                              pattern = "t_n", replacement = "tn" )
# Get indexes
only_nums <- gsub( x = rownames( all_divt$GBM$`PostKPg-GBM_FILT` ),
                   pattern = "_[A-Z]..*", replacement = "" )


#### TABLES ---
#------------------#
# COMPARISON TABLE #
#------------------#
# Get the total number of columns, subtract one as we are not going to keep
# the fourth column with the calibrations (only a check when running 
# MCMC diagnostics)
cols_per_dataset <- length( colnames( all_divt$GBM$`PostKPg-GBM_FILT` ) )-1
short_GBM <- short_ILN <- matrix( 0, nrow = length( nodes_2plot ),
                                  ncol = cols_per_dataset*length(dir_names) )
rownames( short_GBM ) <- rownames( short_ILN ) <- 
  gsub( x = names( nodes_2plot ), pattern = "_", replacement = "|" )
# Prepare a vector to create colnames based on number of datasets
conc_names <- vector( mode = "character", length = length(dir_names)*3 )
start <- end <- 0
for( i in 1:length( dir_names ) ){
  if( i == 1 ){
    start <- 1 
    end   <- start + 2
  }else{
    start <- end + 1
    end <- start + 2
  }
  conc_names[start:end] <- rep( dir_names[i], 3 )
}
# Get colnames ready
colnames( short_GBM ) <- colnames( short_ILN ) <- 
  paste( rep( c( "mean_t", "2.5%q", "97.5%q" ), length( dir_names ) ),
         conc_names, sep = "-" )
# Start saving the output results in the empty vectors for which memory has
# already been allocated in the corresponding cells
start <- end <- 0
for( i in 1:length( dir_names ) ){
  if( i == 1 ){
    start <- 1 
    end   <- cols_per_dataset
  }else{
    start <- end + 1
    end   <- start + cols_per_dataset - 1
  }
  cat( "round ", i, "start = ", start, " | end = ", end, "\n" )
  count <- 0
  for( j in nodes_2plot ){
    count <- count + 1
    tmp_ind <- which( only_nums %in% j )
    short_GBM[count,start:end] <- as.numeric( all_divt$GBM[[ i ]][tmp_ind,1:3] )
    short_ILN[count,start:end] <- as.numeric( all_divt$ILN[[ i ]][tmp_ind,1:3] )
  }
}
# Output the summary tsv files
if( ! dir.exists( paste( main_wd, "figs_tables/tables", sep = "" ) ) ){
  dir.create( paste( main_wd, "figs_tables/tables", sep = "" ) ) 
}
write.table( x = short_GBM, file = paste( main_wd, "figs_tables/tables/",
                                          "compare_calibnodes_divtimes_mammals_GBM.tsv",
                                          sep = "" ),
             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE )
write.table( x = short_ILN, file = paste( main_wd, "figs_tables/tables/",
                                          "compare_calibnodes_divtimes_mammals_ILN.tsv",
                                          sep = "" ),
             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE )

#### PLOTS ----
#----------------#
# START PLOTTING #
#----------------#
# Generate one plot per node indicated above
# If you have run this function inside R script
# `Compare_alldivt_mammals_tables.R`, then you do not need to run it again
# in this script
comparison_plots( sum_obj = all_divt,
                  # Directory `plots` is created inside this directory
                  out_dir = paste( main_wd, "figs_tables/", sep = "" ),
                  plots_per_doc = plots_per_fig,
                  all_nodes = nodes_2plot, lab_nodes = only_nums,
                  data_perclock = 2,
                  # Same order as datasets inside object `all_divt`
                  x_labs = c( "GBM, Post-K-Pg", "GBM, Unconstrained",
                              "ILN, Post-K-Pg", "ILN, Unconstrained" ),
                  # Colours for each dataset
                  points_col =  c( "cyan", "orange",
                                   "cyan", "orange" ),
                  sep_space = rep( c( 0.2, 0.2 ), 2 ),
                  suffix = "",
                  pch_vals = rep( c( 16, 16 ), 2 ),
                  cex_vals = rep( c( 1.3, 1.3 ), 2 ) )


