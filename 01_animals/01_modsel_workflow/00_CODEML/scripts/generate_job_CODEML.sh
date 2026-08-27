#!/bin/bash

# Get args
pipedir=$1 # Path to pipeline dir
name_wd=$2 # Name of the working directory, e.g., `modsel_animals`

# Replace vars in template bash script for job array
cp pipeline_Hessian_CODEML_template.sh $pipedir/pipeline_Hessian.sh
# Replace name of working directory
upd_wd=$( echo $name_wd | sed 's/\//\\\//g' | sed 's/_/\\_/g' )
sed -i 's/WDNAME/'${upd_wd}'/g' $pipedir/pipeline_Hessian.sh