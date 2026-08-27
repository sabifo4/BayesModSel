#!/bin/bash

# Get args
pipeloc=$1     # Path to MCMCtree pipeline dir
runmcmc=$2     # Command to execute MCMCtree
nbeta=$3       # Number of beta values to run the power posteriors
name_wd=$4     # Name of the working directory, e.g., `main`
hypothesis=$5  # Name of tree hypothesis being tested
addCLOCK=$6    # Boolean, if YES, then enable option 7 and modify
change=$7      # Either GBM or ILN. It is required to input something,
               # but only if `addCLOCK` is equal to YES it will be used
               # to replace the var in the template control file!

# First, replace the flags if need be given that some directories
# have either just "HYPOTHESIS" or "HYPOTHESIS_CLOCK"
if [[ $addCLOCK =~ "YES" ]]
then
cp pipeline_MCMCtree_template.sh $pipeloc/$hypothesis"_"$change/pipeline_$hypothesis"_"$change".sh"
out_file=$( echo $pipeloc/$hypothesis"_"$change/pipeline_$hypothesis"_"$change".sh" )
echo Your new pipeline is saved as $out_file
# Replace vars in template bash script for job array
sed -i 's/MCMCtree\/HYPOTHESIS/MCMCtree\/'${hypothesis}'\_'${change}'/g' $out_file
else
echo Your new pipeline is saved as $pipeloc/$hypothesis/pipeline_$hypothesis".sh"
out_file=$( echo $pipeloc/$hypothesis/pipeline_$hypothesis".sh" )
cp pipeline_MCMCtree_template.sh $out_file
# Replace vars in template bash script for job array
sed -i 's/MCMCtree\/HYPOTHESIS/MCMCtree\/'${hypothesis}'/g' $out_file
fi

# Replace rest of vars in template bash script for job array
sed -i 's/calibrated\/HYPOTHESIS/calibrated\/'${hypothesis}'/g' $out_file
sed -i 's/CMDRUN/'${runmcmc}'/' $out_file
sed -i 's/NUMCH/'${nbeta}'/' $out_file

# Replace name of working directory
upd_wd=$( echo $name_wd | sed 's/\//\\\//g' | sed 's/_/\\_/g' )
sed -i 's/WDNAME/'${upd_wd}'/g' $out_file
