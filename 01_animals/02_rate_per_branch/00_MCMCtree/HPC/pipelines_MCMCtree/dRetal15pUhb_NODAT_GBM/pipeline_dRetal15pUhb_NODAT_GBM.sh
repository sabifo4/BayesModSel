#!/bin/bash
#$ -cwd                    # Run the code from the current directory
#$ -V                      # Export environment to job
#$ -j y                    # Merge the standard output and standard error
#$ -l h_rt=240:00:00       # Limit each task to 10 days
#$ -l h_vmem=1G            # Request 1GB RAM
#$ -t 1-6

#==========================================================================================#
# Contact Sandra Alvarez-Carretero for any doubts about this script: sandra.ac93@gmail.com #
#==========================================================================================#

# --------------------------------------- #
# Creating file structure to run MCMCtree #
# --------------------------------------- # 

# 1. Find global dirs for paths
pipeline_dir=$( pwd )
main_dir=$( echo $pipeline_dir | sed 's/\/brate_animals\/..*/\/brate_animals\//' )
cd $main_dir/Hessian
baseml_dir=$( pwd )
name_inBV=$( echo in.BV )
path_inBV=$( echo $baseml_dir/$name_inBV )
cd $main_dir/MCMCtree/dRetal15pUhb_NODAT_GBM/$SGE_TASK_ID/
home_dir=$( pwd )
cd $main_dir/aln/ 
aln_dir=$( pwd )
name_aln=`ls *txt`
path_aln=$( echo $aln_dir/$name_aln | sed 's/\_/\\\_/g' | sed 's/\//\\\//g' | sed 's/\./\\\./g' )
cd $main_dir/trees/dRetal15pUhb_NODAT_GBM
tree_dir=$( pwd )
name_tree=`ls *tree`
path_tree=$( echo $tree_dir/$name_tree | sed 's/\_/\\\_/g' | sed 's/\//\\\//g' | sed 's/\./\\\./g' )

# 3 Create specific log file
exec 3>&1> >(while read line; do echo "$line" >> $pipeline_dir/log.MCMCtree_r$SGE_TASK_ID".txt"; done;) 2>&1
start=`date`
echo Job starts":" $start

# 4. Start analysis
echo The analyses will take place in directory $home_dir
printf "\n"
# Move to analysis dir
cd $home_dir
# Soft link input files
ln -s $aln_dir/$name_aln $name_aln
ln -s $tree_dir/$name_tree $name_tree
ln -s $path_inBV $home_dir/in.BV

# 5. Run MCMCtree
printf "\nRunning MCMCtree for divergence times estimation ... ...\n"
cd $home_dir
mcmctree *.ctl

# 6. Close
printf "\n"
echo MCMCtree FINISHED"!"
printf "\n"
end=`date`
echo Job ends":" $end

