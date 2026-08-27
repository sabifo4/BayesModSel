#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir "00_Gene_filtering"
path_to_file <- getActiveDocumentContext()$path
script_wd <- paste( dirname( path_to_file ), "/", sep = "" )
wd <- gsub( pattern = "/scripts", replacement = "", x = script_wd )
setwd( wd )

#--------------#
# LOAD OBJECTS #
#--------------#
tt <- ape::read.tree( file = "../../trees/uncalib/uncalib_animal.tree" )

#-------#
# TASKS #
#-------#
# 1. Find number of taxa 
num_sp        <- length( tt$tip.label )
spnames       <- tt$tip.label
phylip_header <- paste( num_sp, "  1", sep = "" )

phylip_header_aln <- paste( num_sp, "  2\n", sep = "" )
spnames_2nuc      <- paste( spnames, "     PR", sep = "" )

# 2. Generate dummy aln
write( x = phylip_header_aln, file = "dummy.txt" )
write( x = spnames_2nuc, file = "dummy.txt", append = TRUE )
write( x = paste( "\n", phylip_header_aln, sep = "" ), file = "dummy.txt",
       append = TRUE )
write( x = spnames_2nuc, file = "dummy.txt", append = TRUE )

