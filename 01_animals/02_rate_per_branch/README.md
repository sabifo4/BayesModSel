# Inferring evolutionary rate per branches

## 1. Run `MCMCtree`

Firstly, we need to create the file structure to run our analyses in a HPC. We have provided this file structure here too so that you can also transfer it to your HPC in case you want to reproduce our analyses. If you want to generate this file structure by yourself, make sure that you delete or move the directories we provide inside `02_rate_per_branch` before running the following code (please run the code snippet that suits your OS as specified in the header!):

* If you **are not** a Mac user (or you are a Mac user and have no problems with using the Linux command `sed`), please run the code snippet below:

    ```sh
    # Run from `02_rate_per_branch`
    ## NOTE: Directory `scripts` will be already inside this directory, which contains
    ## the bash scripts that we use to then generate the bash scripts to run job arrays
    for i in `seq 1 6`
    do
    mkdir -p HPC/MCMCtree/{BM23pUhb_GBM/$i,BM23pUhb_ILN/$i,dRetal15pUhb_GBM/$i,dRetal15pUhb_ILN/$i}
    done
    cd HPC
    # Copy alignment
    mkdir aln
    cp ../../../00_data_formatting/aln/dummy/dummy.txt aln
    # Move `in.BV` to directory Hessian
    mkdir Hessian
    cp ../../../01_modsel_workflow/00_CODEML/in.BV Hessian/
    # Copy control file that was already created
    # when carrying out the Bayesian model selection 
    # analysis and remove the `BayesFactor` line
    for i in BM23pUhb_GBM BM23pUhb_ILN dRetal15pUhb_GBM dRetal15pUhb_ILN
    do
    printf "Copying control file for analysis "$i" ...\n"
    for j in `seq 1 6`
    do
    cp ../../../01_modsel_workflow/01_mcmc3r/MCMCtree/$i/1/*ctl MCMCtree/$i/$j/mcmctree_$i.ctl
    sed -i 's/BayesFactorBeta..*//' MCMCtree/$i/$j/*ctl
    sed -i 's/print..*/print\ \=\ 2/' MCMCtree/$i/$j/*ctl
    done
    done
    # Create directories to sample from the prior
    cp -R MCMCtree/BM23pUhb_GBM MCMCtree/BM23pUhb_NODAT_GBM
    cp -R MCMCtree/BM23pUhb_ILN MCMCtree/BM23pUhb_NODAT_ILN
    cp -R MCMCtree/dRetal15pUhb_GBM MCMCtree/dRetal15pUhb_NODAT_GBM
    cp -R MCMCtree/dRetal15pUhb_ILN MCMCtree/dRetal15pUhb_NODAT_ILN
    for i in BM23pUhb_NODAT_GBM BM23pUhb_NODAT_ILN dRetal15pUhb_NODAT_GBM dRetal15pUhb_NODAT_ILN
    do
    printf "Updating control file for analyses that will sample from the prior "$i" ...\n"
    for j in `seq 1 6`
    do
    name_ctl=`ls MCMCtree/$i/$j/`
    new_name_ctl=$( echo mcmctree_$i.ctl )
    mv MCMCtree/$i/$j/$name_ctl MCMCtree/$i/$j/$new_name_ctl
    sed -i 's/usedata..*/usedata\ \=\ 0/' MCMCtree/$i/$j/$new_name_ctl
    done
    done
    # Copy trees
    name_hypotheses=`ls MCMCtree`
    for i in $name_hypotheses
    do
    mkdir -p trees/$i
    done
    # Empirical dataset
    for i in trees/BM23pUhb_*
    do
    cp ../../../00_data_formatting/trees/BM23pUhb/* $i/
    done
    for i in trees/dRetal15pUhb_*
    do
    cp ../../../00_data_formatting/trees/dRetal15pUhb/* $i/
    done
    # Copy scripts
    cp -R ../scripts .
    ```

* If you are a Mac user and have issues with the Linux command `sed`, please run the code snippet below:

    ```sh
    # Run from `02_rate_per_branch`
    ## NOTE: Directory `scripts` will be already inside this directory, which contains
    ## the bash scripts that we use to then generate the bash scripts to run job arrays
    for i in `seq 1 6`
    do
    mkdir -p HPC/MCMCtree/{BM23pUhb_GBM/$i,BM23pUhb_ILN/$i,dRetal15pUhb_GBM/$i,dRetal15pUhb_ILN/$i}
    done
    cd HPC
    # Copy alignment
    mkdir aln
    cp ../../../00_data_formatting/aln/dummy/dummy.txt aln
    # Move `in.BV` to directory Hessian
    mkdir Hessian
    cp ../../../01_modsel_workflow/00_CODEML/in.BV Hessian/
    # Copy control file that was already created
    # when carrying out the Bayesian model selection 
    # analysis and remove the `BayesFactor` line
    for i in BM23pUhb_GBM BM23pUhb_ILN dRetal15pUhb_GBM dRetal15pUhb_ILN
    do
    printf "Copying control file for analysis "$i" ...\n"
    for j in `seq 1 6`
    do
    cp ../../../01_modsel_workflow/01_mcmc3r/MCMCtree/$i/1/*ctl MCMCtree/$i/$j/mcmctree_$i.ctl
    awk '{gsub(/BayesFactorBeta..*/,"")}1' MCMCtree/$i/$j/mcmctree_$i.ctl > MCMCtree/$i/$j/mcmctree_cp_$i.ctl
    awk '{gsub(/print..*/,"print = 2")}1' MCMCtree/$i/$j/mcmctree_cp_$i.ctl > MCMCtree/$i/$j/mcmctree_$i.ctl
    rm MCMCtree/$i/$j/mcmctree_cp_$i.ctl
    done
    done
    # Create directories to sample from the prior
    cp -R MCMCtree/BM23pUhb_GBM MCMCtree/BM23pUhb_NODAT_GBM
    cp -R MCMCtree/BM23pUhb_ILN MCMCtree/BM23pUhb_NODAT_ILN
    cp -R MCMCtree/dRetal15pUhb_GBM MCMCtree/dRetal15pUhb_NODAT_GBM
    cp -R MCMCtree/dRetal15pUhb_ILN MCMCtree/dRetal15pUhb_NODAT_ILN
    for i in BM23pUhb_NODAT_GBM BM23pUhb_NODAT_ILN dRetal15pUhb_NODAT_GBM dRetal15pUhb_NODAT_ILN
    do
    printf "Updating control file for analyses that will sample from the prior "$i" ...\n"
    for j in `seq 1 6`
    do
    name_ctl=`ls MCMCtree/$i/$j/`
    new_name_ctl=$( echo mcmctree_$i.ctl )
    mv MCMCtree/$i/$j/$name_ctl MCMCtree/$i/$j/$new_name_ctl
    awk '{gsub(/usedata..*/,"usedata = 0")}1' MCMCtree/$i/$j/$new_name_ctl > MCMCtree/$i/$j/tmp_$new_name_ctl
    mv MCMCtree/$i/$j/tmp_$new_name_ctl MCMCtree/$i/$j/$new_name_ctl
    done
    done
    # Copy trees
    name_hypotheses=`ls MCMCtree`
    for i in $name_hypotheses
    do
    mkdir -p trees/$i
    done
    # Empirical dataset
    for i in trees/BM23pUhb_*
    do
    cp ../../../00_data_formatting/trees/BM23pUhb/* $i/
    done
    for i in trees/dRetal15pUhb_*
    do
    cp ../../../00_data_formatting/trees/dRetal15pUhb/* $i/
    done
    # Copy scripts
    cp -R ../scripts .
    ```

Now, everything is ready to run `MCMCtree` in a job array manner. To create the file structure to save the output log files when running the job arrays, you can run the following commands:

```sh
# Run from `00_MCMCtree/HPC`
for i in BM23pUhb_GBM BM23pUhb_ILN dRetal15pUhb_GBM dRetal15pUhb_ILN BM23pUhb_NODAT_GBM BM23pUhb_NODAT_ILN dRetal15pUhb_NODAT_GBM dRetal15pUhb_NODAT_ILN
do
mkdir -p pipelines_MCMCtree/$i/
done
```

We have generated a template bash script, [`pipeline_MCMCtree_template.sh`](00_MCMCtree/scripts/pipeline_MCMCtree_template.sh) with flags that will be replaced when running the bash script [`generate_job_MCMCtree.sh`](00_MCMCtree/scripts/generate_job_MCMCtree.sh). These can then be transferred to the HPC (e.g., `rsync`), where they will be submitted. E.g.:

```sh
# Run from `00_MCMCtree`
# Make sure that you have created a directory called `brate_animals`
# in your HPC scratch directory or similar so that these files can
# be successfully copied!
rsync -avz --copy-links HPC/* <uname>@<server>:<path>/brate_animals/
```

Once everything has been transferred to the HPC, just do the following to generate the bash scripts with the job arrays to be later submitted to start the Bayesian model selection analysis:

```sh
# Run from `brate_animals`, the main directory, once you are in the HPC so that the 
# paths are already formatted
home_dir=$( pwd )
cd scripts
chmod 775 *sh
for i in BM23pUhb_GBM BM23pUhb_ILN dRetal15pUhb_GBM dRetal15pUhb_ILN BM23pUhb_NODAT_GBM BM23pUhb_NODAT_ILN dRetal15pUhb_NODAT_GBM dRetal15pUhb_NODAT_ILN
do
# arg1   Path to MCMCtree pipeline dir.
# arg2   Command to execute MCMCtree. E.g. "mcmctree", "mcmctree_4.10.7", etc.
# arg3   Number of chains that will be run (i.e., 6 in our case)
# arg4   Name of the working directory, e.g., `brate_animals` in the HPC
# arg5   Type of hypothesis being tested (i.e., [ACetal22|BM23|BM23pUhb]_[GBM|ILN|NODAT_GBM|NODAT_ILN] )
./generate_job_MCMCtree.sh $home_dir/pipelines_MCMCtree mcmctree 6 brate_animals $i
done
##>NOTE: If you want, you can now transfer these pipelines to your HPC
##>E.g.: Run the following command from your PC
##>      once you are in directory `00_MCMCtree/HPC`
##>      rsync -avz --copy-links <uname>@<server>:<path>/brate_animals/pipelines_MCMCtree .
```

Now, we have everything ready to run `MCMCtree`!

## 2. Analyses with `MCMCtree` when sampling from the prior

### Run `MCMCtree` in the HPC - prior

Firstly, we need to run `MCMCtree` when sampling from the prior (i.e., no data used!) to make sure that the calibration densities that we specified (user-specified prior) are reflected by the marginal density (effective prior) used by `MCMCtree`. You can navigate to each of the directories you created and transferred to the HPC and submit the jobs. E.g.:

```sh
# Run from directories
# `pipelines_MCMCtree/*_NODAT_[GBM|ILN]` 
# in your HPC 
# Please change directories until you are there,
# then run the following commands
# E.g.: if you were inside `BM23_GBM`...
cd BM23_GBM
chmod 775 *sh
qsub pipeline_BM23pUhb_GBM.sh
# Now, just change directories until you have visited
# all *NODAT* directories and submitted all job arrays!
```

### Setting the file structure to analyse `MCMCtree` output - prior

We will now create a `sum_analyses` directory to analyse the `MCMCtree` output. Nevertheless, we first need to transfer the data from the HPC to the corresponding directory on our local PC for further analyses:

```sh
# Run everything from `brate_animals` in your HPC
num_chains=6
mkdir -p tmp_to_transfer/00_prior
cd tmp_to_transfer
for j in BM23pUhb_NODAT_GBM BM23pUhb_NODAT_ILN dRetal15pUhb_NODAT_GBM dRetal15pUhb_NODAT_ILN
do
for i in `seq 1 $num_chains`
do
mkdir -p 00_prior/$j/$i/
# Now, copy the files that are required for sum stats
# We have run 6 chains for analyses sampling from the prior
printf "\n[[ Copying run "$i" for analyses sampling from the prior -- directory "$j" ]]\n\n"
cp -R ../MCMCtree/$j/$i/* 00_prior/$j/$i/
done
grep 'Species tree for FigTree' -A1 ../MCMCtree/$j/1/out.txt | awk 'NR==2' > 00_prior/node_tree_$j.tree
done
```

Now, you can transfer the temporary directory to the local PC, e.g., using `rsync`:

```sh
# Run from `02_rate_per_branch/00_MCMCtree` dir on your local PC
# Please change directories until you are there
# Then run the following commands
mkdir sum_analyses
cd sum_analyses
# Now, trasnfer the data from the HPC
rsync -avz --copy-links <uname>@<logdetails>:<path>/brate_animals/tmp_to_transfer/00_prior .
# You can transfer directory `pipelines_MCMCtree` now or wait until
# you run `MCMCtree` when sampling from the posterior!
```

> [!NOTE]
> If you want to download the output files we generated, please download `02_animals_divtimes_prior.zip` from [our FigShare repository](https://doi.org/10.6084/m9.figshare.32033958) and save directory `pipelines_MCMCtree` inside `01_animals/02_rate_per_branch/00_MCMCtree/sum_analyses`. Then, decompress, `00_prior.zip` and save `00_prior` inside the same directory. Please note that this directory also includes the output files we generate during some of the MCMC diagnostics (e.g., directories `[BM23pUhb|dRetal15pUhb]_NODAT_[GBM|ILN]_*NODAT*/[UpperEdiacaran|Unconstrained]*` and `mcmc_files_*NODAT*`). If you want to run the subsequent steps, you can move these output files to another location and only keep the four directories with the 6 independent chains (i.e., `[BM23pUhb|dRetal15pUhb]_NODAT_[GBM|ILN]/[1-6]`).

### MCMC diagnostics - prior

Now that we have the output files from the different MCMC runs in an organised file structure, we are ready to run MCMC diagnostics!

Firstly, we need to generate a file with calibration information that is compatible with the subsequent scripts. We will use two in-house `R` scripts [`Include_calibrations_MCMCtree.R`](01_MCMC_diagnostics/calibs/scripts/) and [`Merge_node_labels.R`](01_MCMC_diagnostics/scripts/Merge_node_labels.R). The former requires the input csv files that we have already prepared in directory [`calibs/raw`](01_MCMC_diagnostics/calibs/raw/), which follow the format detailed below (this information has also been included in the `R` script):

* Header.
* One row per calibration.
* No spaces at all, semi-colon separated.
* There are 4 columns:
  * Name you want to give to the calibrated node (no spaces!).
  * Name of one of the tips (e.g., tip 1) that leads to MRCA (no spaces!).
  * Name of the other tip (e.g., tip 2) that leads to MRCA (no spaces!).
  * Calibration in `MCMCtree` notation (no spaces!). More details on the `MCMCtree` notation you need to use in the fourth column in the [`PAML` documentation](https://github.com/abacus-gene/paml/blob/master/doc/pamlDOC.pdf).

The output files include trees that can be visualised in graphical interfaces such as `FigTree` or `TreeViewer` (i.e., files will be output in a new directory, `raw/cal_data`). Newly created directory `raw/inp_data` contains formatted tree files that can be used with `MCMCtree`, which we are not using now (we already had them at the beginning of this analysis!). The second script, [`Merge_node_labels.R`](01_MCMC_diagnostics/scripts/Merge_node_labels.R), will then use some of these output files and the uncalibrated tree file (i.e., [`uncalib_animal.tree`](01_MCMC_diagnostics/calibs/raw/uncalib_animal.tree)) to generate one calibration file for each dataset analysed in case there are differences in the tree topologies and the node labels which age is being constrained. They will be saved in a newly created directory called `01_MCMC_diagnostics/calibs/inp_calibs`.

---

> [!NOTE]
> Please note that the script [`Merge_node_labels.R`](01_MCMC_diagnostics/scripts/Merge_node_labels.R) will generate csv output files with calibration information formatted as follows:

```text
Calib;node;Prior
<calibration_tag>;<MCMCtree_node_number>;'<MCMCtree_calibration_format>'
```

The output files generated by [`Merge_node_labels.R`](01_MCMC_diagnostics/scripts/Merge_node_labels.R) will be read by other in-house `R` scripts to extract the calibration information and append it to summary tables and plots.

---

Once the calibration files are ready, we can run the R script [`MCMC_diagnostics_prior.R`](01_MCMC_diagnostics/scripts/MCMC_diagnostics_prior.R) and follow the detailed step-by-step instructions written in the scripts. In a nutshell, the protocol to follow will be the following:

1. Load the `mcmc.txt` files generated after each run.
2. Generate convergence plots with unfiltered chains.
3. Find whether there are major differences between the time estimates sampled across the chains for the same nodes in the 97.5% and the 2.5% quantiles. If so, flag and delete said chains.
4. If some chains have not passed the filters mentioned above, create an object with the chains that have passed the filters.
5. Generate new convergence plots with those chains that passed the filters.
6. Calculate Rhat, tail-ESS, and bulk-ESS to assess chain convergence and mixing with the chains that have passed filters.

> [!IMPORTANT]
> Once you have run the aforementioned `R` script, you have to visit the content of directories `sum_analyses/00_prior/*NODAT_[GBM|ILN]/` to check whether the chains ran under each hypothesis and clock model passed the filters.
>
> An example of how the directory should look like after the MCMC diagnostics ran **if all chains have passed the filters** can be found inside `00_prior/BM23pUhb_NODA_ILN/`. If you visit this directory, you will see that, apart from directories `1` to `6` (i.e., one directory per chain; we ran 6 chains), there is another directory called `UpperEdiacaran-NODAT-ILN`. There are no additional text files or directories with the name of the model, and so all 6 chains have passed the filters. The `tsv` files you should find inside `UpperEdiacaran-NODAT-ILN` have the summarised posterior mean divergence times and corresponding CIs for each node (if the `MCMCtree` control file set `print = 2`, then the branch rates and corresponding CIs should have also been summarised in additional tsv files). Please note that the file with suffix `*all_mean_est.tsv` has both mean divergence times and CIs, and so we recommended you use this file to generate a final table with the timetree inference results (a file with suffix `*all_mean_est_brate.tsv` will also be there if the `MCMCtree` control file set `print = 2`). In addition, if you navigate to directory `out_RData` (i.e., `01_MCMC_diagnostics/out_RData`), you will find individual tsv files for each model analysed with the MCMC summary stats. You can also check the convergence plots inside directory `01_MCMC_diagnostics/plots/ESS_and_chains_convergence`.
>
> **If some chains had NOT passed the filters**, we would have found an additional directory. E.g.: following the directory name used in the example above, we would have seen `BM23pUhb_NODAT_GBM`, `dRetal15pUhb_NODAT_GBM` and `dRetal15pUhb_NODAT_ILN`. A directory with suffix `FILT` will only be created if there are chains that have not passed the filters. In addition, two files will also be created: `chains_kept.txt` and `check_chains.txt`. The former can be used to know which chains have been used to run the final MCMC diagnostics (i.e., results under `*FILT` directory) and the latter to know which chains are the problematic ones, and thus discarded from the final analyses. If there were nodes which Rhat was higher than 1.05, then such chains would be labelled as "problematic", and they would be listed in another text file starting with `problem_nodes_conv_<mod>`, being `<mod>` the flag given to the analysed dataset (e.g., `BM23pUhb_NODAT_GBM`). You would also find the convergence plots with the unfiltered chains inside directory `01_MCMC_diagnostics/plots/ESS_and_chains_convergence`.

The MCMC diagnostics did not find any of the chains problematic after running [our in-house `R` script `MCMC_diagnostics_prior.R`](01_MCMC_diagnostics/scripts/MCMC_diagnostics_prior.R) (i.e., no `chains_kept.txt` or `check_chains.txt` files saved under `*NODAT-[GBM|ILN]*` directories to check which ones passed the filters). Therefore, we used [our in-house bash script `Combine_MCMC.sh`](01_MCMC_diagnostics/scripts/Combine_MCMC.sh) to concatenate those `mcmc.txt` files for the chains that passed the filters in a unique file:

```sh
# Run from `01_MCMC_diagnostics/scripts`
cp Combine_MCMC.sh ../../00_MCMCtree/sum_analyses/00_prior/
# One argument taken: number of chains
cd ../../00_MCMCtree/sum_analyses/00_prior/
## Variables needed
## arg1 --> name of directory where analyses have taken place (e.g., NODAT, GBM, ILN)
## arg2 --> output dir: mcmc_files_NODAT, mcmc_files_GBM, mcmc_files_ILN, etc.
## arg3 --> "`seq 1 36`", "1 2 5", etc. | depends on whether some chains were filtered out or not
## arg4 --> clock model used: ILN, GBM, NODAT
## arg5 --> number of samples specified to collect in the control file to run `MCMCtree`
## arg6 --> 'Y' to generate a directory with files compatible with programs such as `Tracer` to visually
##          inspect traceplots, convergence plots, etc. 'N' otherwise
## arg7 --> if arg6 is 'Y', arg7 needs to have a name for the `mcmcf4traces_<name>` that will be
##          generated. If `arg6` is equal to 'N', then please write `N` too
./Combine_MCMC.sh BM23pUhb_NODAT_GBM mcmc_files_BM23pUhb_NODAT_GBM "2 4" NODAT 20000 Y BM23pUhb_NODAT_GBM
./Combine_MCMC.sh BM23pUhb_NODAT_ILN mcmc_files_BM23pUhb_NODAT_ILN "`seq 1 6`" NODAT 20000 Y BM23pUhb_NODAT_ILN
./Combine_MCMC.sh dRetal15pUhb_NODAT_GBM mcmc_files_dRetal15pUhb_NODAT_GBM "2 4 5 6" NODAT 20000 Y dRetal15pUhb_NODAT_GBM
./Combine_MCMC.sh dRetal15pUhb_NODAT_ILN mcmc_files_dRetal15pUhb_NODAT_ILN "1 2 3 4 6" NODAT 20000 Y dRetal15pUhb_NODAT_ILN
```

The script above will generate directories called `mcmc_files*_NODAT_[GBM|ILN]` inside the `00_prior` directory, where the `mcmc.txt` with the concatenated samples will be saved. In addition, directories with individual `mcmc.txt` files of those chains that passed the filters will be created (i.e., see `mcmcf4traces*_NODAT_[GBM|ILN]` directories); you can open such files with programs like `Tracer` to assess the traces and run other visual MCMC diagnostics.

> [!NOTE]
> When you download `02_animals_divtimes_prior` from [our FigShare repository](https://doi.org/10.6084/m9.figshare.32033958), you will not find the `mcmcf4traces*_NODAT_[GBM|ILN]` directories as they are somewhat redundant. If you want to generate them, please make sure you run the code above. Please note that directories `mcmc_files*_NODAT_[GBM|ILN]` will also be generated, so you should move the ones you downloaded from our FigShare elsewhere so that they do not overwrite (in case you want to compare them).

We will now reuse the dummy alignment with only 2 nucleotides that we had already generated (i.e., check file [`dummy.txt`](00_MCMCtree/HPC/aln/dummy.txt)) so that we can obtain the `FigTree` files with the mean time estimates obtained when using the concatenated `mcmc.txt` files (i.e., all the samples that passed the filters and were collected when running the chains checked in previous steps):

```sh
# Run from `01_MCMC_diagnostics`
# Please change directories until you are 
# in the directory aforementioned
mkdir dummy_aln
cp ../00_MCMCtree/HPC/aln/dummy.txt dummy_aln/
```

We have also generated a dummy control file to read the dummy alignment. Additionally, we have enabled option `print = -1`. This print setting lets `MCMCtree` know that an MCMC is not to be run. Instead, `MCMCtree` is told to read the input files (file with the dummy alignment, the calibrated tree file, and the concatenated `mcmc.txt` file) and summarise the samples in the `mcmc.txt` (those that were collected from those chains that passed the filters!). The final mean estimated divergence times (and branch rates, if `print = 2` was set when running `MCMCtree`) and the corresponding CIs will be written in the output FigTree.tre file

```sh
##> [IMPORTANT] Before running the `for` loop below,
##> please check the commented sections and ammend the
##> commands accordingly depending on whether you are
##> using a program that is exported to the system's path
##> or a binary that needs to be execute with a relative
##> path

# Run from `00_MCMCtree/sum_analyses/00_prior`
# Please change directories until you are 
# in the directory aforementioned
name_dat=( 'BM23pUhb_NODAT_GBM' 'BM23pUhb_NODAT_ILN' 'dRetal15pUhb_NODAT_GBM' 'dRetal15pUhb_NODAT_ILN' )
num_dat=4
count=-1
for i in `seq 1 $num_dat`
do
count=$(( count + 1 ))
printf "\n[[ Analysing dataset "${name_dat[count]}" ]]\n"
base_dir=$( pwd )
cd ../../../01_MCMC_diagnostics/dummy_ctl_files/
ctl_dir=$( pwd )
cd ../../00_MCMCtree/HPC/trees/${name_dat[count]}/
tt_dir=$( pwd )
name_tt=`ls *tree`
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *txt`
cd $base_dir
cd mcmc_files_${name_dat[count]}
printf "[[ Generating tree file for concatenated \"mcmc.txt\"  ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
#
##> [IMPORTANT] Change the command below if you are using a different alias to run
##> `MCMCtree` that is not `mcmctree` and/or add the relative paths if you have not
##> exported `MCMCtree` to the system's path!
##> Please remember we are using PAML v4.9h, but you may have a different version
##> installed on your PC. If you have a version that was released after
##> PAML v4.10.7, please note that you will need to add a "C" as a fourth option
##> in variable "BDparas" in the dummy control. In addition, you will have
##> 95%CIs instead of 95%HPDs -- you may want to also change the suffic of the
##> tree file if that is the case!
mcmctree *ctl
printf "\n"
mv FigTree.tre FigTree_${name_dat[count]}"_95HPD.tree"
cd $base_dir
done
```

We now have our timetree inferred with all the samples collected by all the chains that passed the filters during MCMC diagnostics (when sampling from the prior)! The next step is to plot the calibration densities VS the marginal densities to verify whether there are any serious clashes that may arise because of truncation or problems with the fossil calibrations used. We will use the [in-house `R` script `Check_priors_margVScalib.R`](01_MCMC_diagnostics/scripts//Check_priors_margVScalib.R) to generate these plots.

Once this script has finished, you will see that a new directory `plots/margVScalib` will have been created. Inside this directory, you will find one directory for each individual dataset with individual plots for each node. In addition, all these plots have been merged into a unique document as well (note: some plots may be too small to see for each node, hence why we have generated individual plots). They have been plotted in `JPG`, `PDF`, and `TIFF` format.

> [!NOTE]
> Please note that we have kept directory `sum_files_prior` in this repository (see code snippet below). Nevertheless, if you want to also download the `plots` directory (this is generated and populated during MCMC diagnostics), please download `03_animals_MCMCdiagnostics.zip` from [our FigShare repository](https://doi.org/10.6084/m9.figshare.32033958). Then, decompress this directory and `plots.zip` and save directory `plots` inside `01_animals/02_rate_per_branch/01_MCMC_diagnostics`. Without the `plots` directory, the commands below will not work.

Now, once the MCMC diagnostics have finished, you can extract the relevant output that we used to write our manuscript:

```sh
# Run from `01_MCMC_diagnostics`
mkdir -p sum_files_prior/{MCMC_diagnostics,sum_brates_tsv,sum_divtimes_tsv,sum_divtimes_trees}
cp -R ../00_MCMCtree/sum_analyses/00_prior/mcmc_files*NODAT*/*NODAT*tree sum_files_prior/sum_divtimes_trees
cp -R ../00_MCMCtree/sum_analyses/00_prior/*NODAT*/*/*all_mean_est.tsv sum_files_prior/sum_divtimes_tsv
cp -R ../00_MCMCtree/sum_analyses/00_prior/*NODAT*/*/*all_mean_est_brate.tsv sum_files_prior/sum_brates_tsv
# If there are any problems with nodes that may have not reached convergence, these files will
# be copied; otherwise a warning will be raised saying that there is no such file or directory
cp -R ../00_MCMCtree/sum_analyses/00_prior/*NODAT*/problem_nodes*txt sum_files_prior/MCMC_diagnostics
cp -R plots/ESS_and_chains_convergence/*prior*pdf sum_files_prior/MCMC_diagnostics
cp -R plots/margVScalib sum_files_prior/MCMC_diagnostics
# Keep only PDF files to save space
rm sum_files_prior/MCMC_diagnostics/margVScalib/*/*tif sum_files_prior/MCMC_diagnostics/margVScalib/*/*jpg sum_files_prior/MCMC_diagnostics/margVScalib/*/ind/*tif sum_files_prior/MCMC_diagnostics/margVScalib/*/ind/*jpg
cp -R out_RData/prior*tsv sum_files_prior/MCMC_diagnostics
```

## 3. Analyses with `MCMCtree` when sampling from the posterior

### Run `MCMCtree` in the HPC - posterior

Now that we have verified that there are no issues between the calibration and marginal densities, we can run `MCMCtree` when sampling from the posterior. We will do these analyses under the GBM and ILN relaxed-clock models using the code snippet below:

```sh
# Run from directories
# `pipelines_MCMCtree/*_[GBM|ILN]` in your HPC 
# Please change directories until you are there,
# then run the following commands
# E.g.: if you were inside `dRetal15_GBM`...
cd dRetal15_GBM
chmod 775 *sh
qsub pipeline_dRetal15_GBM.sh
# Now, just change directories until you have visited
# all directories in which `MCMCtree` was run
# to sample from the posterior:
#   'BM23pUhb_GBM' 'BM23pUhb_ILN'
#   'dRetal15pUhb_GBM' 'dRetal15pUhb_ILN' 
```

> [!IMPORTANT]
> When sampling from the posterior, the likelihood is being calculated or approximated, depending on the `userdata` option you set in the control file to run `MCMCtree`. In other words, the larger the dataset, the more time it will take for `MCMCtree` to finish.

### Setting the file structure to analyse `MCMCtree` output - posterior

Firstly, we will summarise the data collected when sampling from the posterior available in the HPC:

```sh
# Run everything from `brate_animals` in your HPC
num_chains=6
mkdir -p tmp_to_transfer/01_posterior
cd tmp_to_transfer
for j in BM23pUhb_GBM BM23pUhb_ILN dRetal15pUhb_GBM dRetal15pUhb_ILN
do
for i in `seq 1 $num_chains`
do
mkdir -p 01_posterior/$j/$i/
# Now, copy the files that are required for sum stats
# We have run 6 chains for analyses sampling from the posterior
printf "\n[[ Copying run "$i" for analyses sampling from the posterior -- directory "$j" ]]\n\n"
cp -R ../MCMCtree/$j/$i/* 01_posterior/$j/$i/
done
done
```

Now, you can transfer the temporary directory to the local PC, e.g., using `rsync`:

```sh
# Run from `02_rate_per_branch/00_MCMCtree/sum_analyses` dir on your local PC
# Please change directories until you are there
# Now, trasnfer the data from the HPC
rsync -avz --copy-links <uname>@<logdetails>:<path>/brate_animals/tmp_to_transfer/01_posterior .
rsync -avz --copy-link <uname>@<logdetails>:<path>/brate_animals/pipelines_MCMCtree .
# Remove blank output files
rm pipelines_MCMCtree/*/*sh.o*
```

> [!NOTE]
> If you want to download the output files we generated, please download `02_animals_divtimes_posterior.zip` from [our FigShare repository](https://doi.org/10.6084/m9.figshare.32033958) and save directory `pipelines_MCMCtree` inside `01_animals/02_rate_per_branch/00_MCMCtree/sum_analyses`. Then, decompress, `01_posterior.zip` and save `01_posterior` inside the same directory. Please note that this directory also includes the output files we generate during some of the MCMC diagnostics (e.g., directories `[BM23pUhb|dRetal15pUhb]_[GBM|ILN]/[UpperEdiacaran|Unconstrained]*` and `mcmc_files_*`). If you want to run the subsequent steps, you can move these output files to another location and only keep the four directories with the 6 independent chains (i.e., `[BM23pUhb|dRetal15pUhb]_[GBM|ILN]/[1-6]`).

### MCMC diagnostics - posterior

Now that we have the output files from the different MCMC runs in an organised file structure, we are ready to check the chains for convergence!

We are going to run the `R` script [`MCMC_diagnostics_posterior.R`](01_MCMC_diagnostics/scripts/MCMC_diagnostics_posterior.R) and follow the detailed step-by-step instructions detailed in the script, which are essentially the same ones we used when analysing the chains when sampling from the prior.

Some chains did not pass the filters (i.e., see `chains_kept.txt` files saved under `01_posterior/*[GBM|ILN]/` directories to check which ones passed the filters; see `check_chains.txt` files saved under `01_posterior/*[GBM|ILN]/` directories to check which chains were flagged due to larger differences than the threshold around the tails), and so we will summarise the samples collected by those chains that did indeed pass the filters:

```sh
# Run from `01_MCMC_diagnostics/scripts`
cp Combine_MCMC.sh ../../00_MCMCtree/sum_analyses/01_posterior/
# One argument taken: number of chains
cd ../../00_MCMCtree/sum_analyses/01_posterior/
## Variables needed
## arg1 --> name of directory where analyses have taken place (e.g., NODAT, GBM, ILN)
## arg2 --> output dir: mcmc_files_NODAT, mcmc_files_GBM, mcmc_files_ILN, etc.
## arg3 --> "`seq 1 36`", "1 2 5", etc. | depends on whether some chains were filtered out or not
## arg4 --> clock model used: ILN, GBM, NODAT
## arg5 --> number of samples specified to collect in the control file to run `MCMCtree`
## arg6 --> 'Y' to generate a directory with files compatible with programs such as `Tracer` to visually
##          inspect traceplots, convergence plots, etc. 'N' otherwise
## arg7 --> if arg6 is 'Y', arg7 needs to have a name for the `mcmcf4traces_<name>` that will be
##          generated. If `arg6` is equal to 'N', then please write `N` too
./Combine_MCMC.sh BM23pUhb_GBM mcmc_files_BM23pUhb_GBM "2 3 4 6" NODAT 20000 Y BM23pUhb_GBM
./Combine_MCMC.sh BM23pUhb_ILN mcmc_files_BM23pUhb_ILN "1 2 3 4 6" NODAT 20000 Y BM23pUhb_ILN
./Combine_MCMC.sh dRetal15pUhb_GBM mcmc_files_dRetal15pUhb_GBM "3 5 6" NODAT 20000 Y dRetal15pUhb_GBM
./Combine_MCMC.sh dRetal15pUhb_ILN mcmc_files_dRetal15pUhb_ILN "`seq 1 6`" NODAT 20000 Y dRetal15pUhb_ILN
```

Once the scripts above have finished, new directories called `mcmc_files*_[GBM|ILN]` and `mcmcf4traces*_[GBM|ILN]` will be created inside `01_posterior/`, respectively.

> [!NOTE]
> When you download `02_animals_divtimes_prior` from [our FigShare repository](https://doi.org/10.6084/m9.figshare.32033958), you will not find the `mcmcf4traces*_[GBM|ILN]` directories as they are somewhat redundant. If you want to generate them, please make sure you run the code above. Please note that directories `mcmc_files*_[GBM|ILN]` will also be generated, so you should move the ones you downloaded from our FigShare elsewhere so that they do not overwrite (in case you want to compare them).

To map the mean time estimates with the filtered chains (mean rate estimates, if you ran `print = 2`!), we need to copy a control file, the calibrated Newick tree, and the previously created dummy alignment:

```sh
##> [IMPORTANT] Before running the `for` loop below,
##> please check the commented sections and ammend the
##> commands accordingly depending on whether you are
##> using a program that is exported to the system's path
##> or a binary that needs to be execute with a relative
##> path

# Run from `00_MCMCtree/sum_analyses/00_prior`
# Please change directories until you are 
# in the directory aforementioned
name_dat=( 'BM23pUhb_GBM' 'BM23pUhb_ILN' 'dRetal15pUhb_GBM' 'dRetal15pUhb_ILN' )
num_dat=4
count=-1
for i in `seq 1 $num_dat`
do
count=$(( count + 1 ))
printf "\n[[ Analysing dataset "${name_dat[count]}" ]]\n"
base_dir=$( pwd )
cd ../../../01_MCMC_diagnostics/dummy_ctl_files/
ctl_dir=$( pwd )
cd ../../00_MCMCtree/HPC/trees/${name_dat[count]}/
tt_dir=$( pwd )
name_tt=`ls *tree`
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *txt`
cd $base_dir
cd mcmc_files_${name_dat[count]}
printf "[[ Generating tree file for concatenated \"mcmc.txt\"  ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
#
##> [IMPORTANT] Change the command below if you are using a different alias to run
##> `MCMCtree` that is not `mcmctree` and/or add the relative paths if you have not
##> exported `MCMCtree` to the system's path!
##> Please remember we are using PAML v4.9h, but you may have a different version
##> installed on your PC. If you have a version that was released after
##> PAML v4.10.7, please note that you will need to add a "C" as a fourth option
##> in variable "BDparas" in the dummy control. In addition, you will have
##> 95%CIs instead of 95%HPDs -- you may want to also change the suffic of the
##> tree file if that is the case!
mcmctree *ctl
printf "\n"
mv FigTree.tre FigTree_${name_dat[count]}"_95HPD.tree"
cd $base_dir
done
```

Now, once the MCMC diagnostics have finished, we can run our [in-house `R` script](01_MCMC_diagnostics/scripts/Check_priors_VS_posteriors.R) to plot the posterior time distributions against the those distributions inferred when the data were not used as well as to compare the posterior time densities inferred under different diversification models.

> [!NOTE]
> Please note that we have kept directory `sum_files_posterior` in this repository (see code snippet below). Nevertheless, if you want to also download the `plots` directory (this is generated and populated during MCMC diagnostics), please download `03_animals_MCMCdiagnostics.zip` from [our FigShare repository](https://doi.org/10.6084/m9.figshare.32033958). Then, decompress this directory and `plots.zip` and save directory `plots` inside `01_animals/02_rate_per_branch/01_MCMC_diagnostics`. Without the `plots` directory, the commands below will not work. If you have followed this step above, you do not need to repeat this (the `plots` directory included output files for both diagnostics with samples collected when the target distributions was the posterior and the prior).

Lastly, you can extract the relevant output that we used to write our manuscript as it follows:

```sh
# Run from `01_MCMC_diagnostics`
# Please change directories until you are 
# in the directory aforementioned
mkdir -p sum_files_post/{MCMC_diagnostics,sum_brates_tsv,sum_divtimes_tsv,sum_divtimes_trees,priorVSpost_plots,compare_plots}
cp -R ../00_MCMCtree/sum_analyses/01_posterior/mcmc_files_*/FigTree*tree sum_files_post/sum_divtimes_trees
cp -R ../00_MCMCtree/sum_analyses/01_posterior/*/*/*all_mean_est.tsv sum_files_post/sum_divtimes_tsv
cp -R ../00_MCMCtree/sum_analyses/01_posterior/*/*/*all_mean_est_brate.tsv sum_files_post/sum_brates_tsv
# If there are any problems with nodes that may have not reached convergence, these files will
# be copied; otherwise a warning will be raised saying that there is no such file or directory
cp -R ../00_MCMCtree/sum_analyses/01_posterior/*/problem_nodes*txt sum_files_post/MCMC_diagnostics
cp -R out_RData/post*tsv sum_files_post/MCMC_diagnostics
cp -R plots/ESS_and_chains_convergence/*post*pdf sum_files_post/MCMC_diagnostics
cp -R plots/Compare_post*pdf sum_files_post/compare_plots
cp -R plots/priorVSpost*pdf sum_files_post/priorVSpost_plots
```

## 4. Summarising branch rate estimates in `R`

Now, we are going to extract the relevant information about the estimated branch rates. In order to do so, we are going to extract the "rategrams" that appear in the `out.txt` files output by `MCMCtree`:

```sh
# Run from `02_Summary_Rates/`
## [[ PRIOR ]]
for i in BM23pUhb_NODAT_GBM BM23pUhb_NODAT_ILN dRetal15pUhb_NODAT_GBM dRetal15pUhb_NODAT_ILN
do
printf "[[ PARSING DATASET "$i" ]] \n"
# Reset counter for every partition
count=-1
# Check if `chains_kept` exists
if [[ -f ../00_MCMCtree/sum_analyses/00_prior/$i/chains_kept.txt ]]
then
ch_passed=$( cat ../00_MCMCtree/sum_analyses/00_prior/$i/chains_kept.txt | tr -d '\r' | tr -d '\n' )
else
ch_passed=$( echo 1 2 3 4 5 6 )
fi
# There are 2 partitions
for j in `seq 1 2`
do
# Update counter according to each partition
count=$(( count + 3 ))
mkdir -p 00_prior/$i/$j
for k in $ch_passed
do
# Grep corresponding rategram for each dataset and partition
# To save it in a unique file, use 'NR==2,NR==5,NR==8,NR==11'
# Expected values:
# NR==2 (rategram p1) | NR==5 (rategram p2)
# NR==8 (rategram p3) | NR==11 (rategram p4)
printf "Extracting rategram from dataset "$i", chain "$k", and partition "$j" ... ...\n"
grep 'rategram' -A1 ../00_MCMCtree/sum_analyses/00_prior/$i/$k/out.txt | awk -v count="$count" 'NR==count' > 00_prior/$i/$j/brates_tree_$i"_p"$j"_ch"$k".tree"
done
printf "  ~> Partition "$j" for all chains extracted! \n\n"
done
done

## [[ POSTERIOR ]]
for i in BM23pUhb_GBM BM23pUhb_ILN dRetal15pUhb_GBM dRetal15pUhb_ILN
do
printf "[[ PARSING DATASET "$i" ]] \n"
# Reset counter for every partition
count=-1
# Check if `chains_kept` exists
if [[ -f ../00_MCMCtree/sum_analyses/01_posterior/$i/chains_kept.txt ]]
then
# Read indexes for chains that passed the filters, remove carriage and new line
ch_passed=$( cat ../00_MCMCtree/sum_analyses/01_posterior/$i/chains_kept.txt | tr -d '\r' | tr -d '\n' )
else
ch_passed=$( echo 1 2 3 4 5 6 )
fi
# There are 2 partitions
for j in `seq 1 2`
do
# Update counter according to each partition
count=$(( count + 3 ))
mkdir -p 01_posterior/$i/$j
for k in $ch_passed
do
# Grep corresponding rategram for each dataset and partition
# To save it in a unique file, use 'NR==2,NR==5,NR==8,NR==11'
# Expected values:
# NR==2 (rategram p1) | NR==5 (rategram p2)
# NR==8 (rategram p3) | NR==11 (rategram p4)
printf "Extracting rategram from dataset "$i", chain "$k", and partition "$j" ... ...\n"
grep 'rategram' -A1 ../00_MCMCtree/sum_analyses/01_posterior/$i/$k/out.txt | awk -v count="$count" 'NR==count' > 01_posterior/$i/$j/brates_tree_$i"_p"$j"_ch"$k".tree"
done
printf "  ~> Partition "$j" for all chains that passed the filters extracted! \n\n"
done
done
```

Once the rategrams have been saved for each dataset, partition, and filtered chain, we can run our in-house `R` script [`Summarise_brates_animals.R`](02_Summary_Rates/scripts/Summarise_brates_animals.R), which will generate one rategram per partition and dataset with the mean estimated branch rates for each partition. The output plots will be saved under directory `02_Summary_Rates/plots/rategrams`, in subdirectories with the dataset name.

## 5. Generating summary tables and plots

If you navigate to directory [figs_table](figs_tables/), you shall find three different scripts that you can use to generate informative summary tables and plots that are ready to use for manuscript submission. Please make sure that you have already generated directory `01_MCMC_diagnostics/sum_files_post/sum_divtimes_tsv` (see sections regarding MCMC diagnostics above) and directories `mcmc_files_*` inside `00_MCMCtree/sum_analyses/0*/mcmc_files*`. In this repository, you shall find the [R scripts](figs_tables/scripts), the [output plots](figs_tables/plots/), and the [output summary tables](figs_tables/tables/). If you want to use our `RData` objects so that you do not need to wait for the scripts to read, parse, and process input files, you can download `04_animals_timetrees_suppmat.zip` [our FigShare repository](https://doi.org/10.6084/m9.figshare.32033958), decompress it, decompress `out_RData.zip`, and save directory `out_RData` inside directory `figs_tables`.

The scripts are self-explanatory so... Just open them, follow the instructions and guidelines, and generate the summary plots and tables!
