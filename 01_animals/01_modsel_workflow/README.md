# Bayesian inference of species divergences

We already have the calibrated tree files and the alignment file but, unfortunately, we do not have the `in.BV` file with the gradient and the Hessian. Consequently, we first need to run `CODEML` so that we can approximate the likelihood calculation during timetree inference when running `MCMCtree` ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)) to save computational time.

## 1. Hessian and gradient inference with `CODEML`

Before running `MCMCtree`, we need to calculate the branch lengths, the gradient, and the Hessian required to approximate the likelihood calculation by `MCMCtree` ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)).

A [template control file](00_CODEML/control_files/prepcodeml.ctl) with the settings used by [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066) and [the file with the rate matrix](00_CODEML/control_files/lg.dat) have been saved in the [`control_files`](00_CODEML/control_files) directory. Note that some of the options in the control file have flags that we will replace when running our in-house pipelines with the correct paths to the input files.

### Set up the file structure

The file structure we will use is the following:

```text
modsel_animals_b128/
  |
  |- alignments/
  |    |- *phy/ # Alignment file in PHYLIP format
  |       
  |- control_files/ # Pre-defined control file with flags to be later replaced with specific settings
  |
  |- Hessian/
  |    |- 1 # Directory to estimate the gradient and the Hessian under LG+G4+F
  |          
  |- pipelines_Hessian # Directory where the pipeline to run `CODEML` will be executed
  |
  |- scripts # Scripts used to prepare control files to run `CODEML
  |
  |- trees
      |- calibrated   # Directory with the calibrated tree for `MCMCtree` -- to be used later!
      |- uncalibrated # Directory with the uncalibrated tree for `CODEML`
```

To create the file structure for the HPC, we will run the following commands from a PC before transferring it to a cluster:

```sh
# Run the following commands from `01_animals/01_modsel_workflow/00_CODEML`
mkdir -p HPC/modsel_animals_b128
cd HPC/modsel_animals_b128 # The main wd will be `modsel_animals_b128`
# Create directory for alignments and create two subdirectories for each
# individual alignment file for the partitions
mkdir -p alignments/{1,2}
# We will be running `CODEML` under one AA model, `model = 3` (empirical + F),
# and for the two partitions in the main alignment file
mkdir -p Hessian/{1,2}/prepare_codeml
# Create dirs for the pipelines, control files, trees, and scripts
mkdir pipelines_Hessian/
mkdir control_files
mkdir -p trees/{uncalibrated,calibrated}
# Create subdirectories for each tree hypothesis under `calibrated`
mkdir -p trees/calibrated/{dRetal15pUhb,BM23pUhb}
```

Once the file structure is created, we can now populate it with the input files: alignment, tree, and control files. We will also add the [`lg.dat` file](00_CODEML/control_files/lg.dat), which has the matrix to enable the LG protein substitution model ([Le and Gascuel, 2008](https://academic.oup.com/mbe/article/25/7/1307/1041491)):

```sh
# Run from `00_CODEML/HPC/modsel_animals_b128`
# First, copy the alignment
cp ../../../../00_data_formatting/aln/aln-2P.phy alignments/
cp ../../../../00_data_formatting/aln/aln_part1.phy alignments/1
cp ../../../../00_data_formatting/aln/aln_part2.phy alignments/2
# Now, copy the the uncalibrated tree
cp ../../../../00_data_formatting/trees/uncalib/uncalib_animal.tree trees/uncalibrated/
# Next, copy the control file and the `lg.dat` file
cp ../../control_files/* control_files/
# Copy scripts
cp -R ../../scripts .
# Now, you can transfer this `modsel_animals_b128` directory to the cluster!
# We show below how to do this with `rsync`, but you may choose a different approach
# Move one directory back so that you are inside `01_PAML/00_CODEML/HPC`
cd ../
rsync -avz --copy-links modsel_animals_b128 <uname>@<server>:<path_to_main_wd>
```

Now, we need to generate other input files to estimate the branch lengths, the gradient, and the Hessian: the input control files for `CODEML`.

To do this in a reproducible manner, you can use the [script `generate_prepcodeml.sh`](00_CODEML/scripts/generate_prepcodeml.sh), which you can find in the [`scripts`](00_CODEML/scripts) directory and which you should have already transferred to the HPC. The [`generate_prepcodeml.sh` script](00_CODEML/scripts/generate_prepcodeml.sh) needs two arguments:

* Argument 1: integer. If `1`, you will enable `model = 2` (empirical model) and, if `2`, you will enable `model = 3` (empirical+F model). Given that we want the latter, we will then type `1`.
* Argument 2: integer. Number of the directory that will be created to host all the input and output files when running `CODEML`, choose from `1` to `n`, being `n` the maximum number of analyses you are running. In this case, we will have two alignments (two alignment block or partitions) and we are running `CODEML` under one AA substitution model, and so we will run this script in a loop during which the second argument will be `1` and then `2`.

```sh
# Run from `modsel_animals_b128/scripts` on the HPC, not on your PC!!
# This is important because relative/absolute paths will be used
# Please change directories until
# you are there. Then, run the following
# commands.
chmod 775 *sh
for i in `seq 1 2`
do
./generate_prepcodeml.sh 1 $i
done
```

To make sure that all the paths have been properly extracted, you can run the following code snippet:

```sh
# Run from `modsel_animals_b128/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
grep 'seqfile' */prepare_codeml/*ctl
grep 'treefile' */prepare_codeml/*ctl
grep 'aaRatefile' */prepare_codeml/*ctl
grep ' ndata' */prepare_codeml/*ctl # We should see 1 as there is only 1 partition per file!
```

### Run `CODEML`

#### Preparing input files

Now that we have the input files (alignment and tree files) and the instructions to run `CODEML` (control file), we will be manually running `MCMCtree` inside each `prepare_codeml` directory (see file structure above) in a special mode that launches `CODEML` for the purpose want: calculating the branch lengths, the gradient, and the Hessian.

```sh
# Run `MCMCtree` from
# `modsel_animals_b128/Hessian/[1,2]/prepare_codeml`
# dir on the HPC. Once you have finished in dir 1,
# ten do the same for dir 2!
# Please change directories until
# you are in there.
# The first command to change directories 
# will work if you are still in 
# `main/Hessian`, otherwise ignore and 
# move to such directory with the command
# that best suits your current directory.
cd 1/prepare_codeml
mcmctree prepcodeml.ctl
```

First, you will see that `MCMCtree` starts parsing the first locus. Then, you will see something like the following printed on your screen (some sections may change depending on the PAML version you have installed on your PC!):

```text
*** Locus 1 ***
running codeml tmp0001.ctl

AAML in paml version 4.9h, March 2018
ns = 54         ls = 20220
Reading sequences, sequential format..
Reading seq #54: Xenoturbel       
Sequences read..

20220 site patterns read, 22811 sites
Counting frequencies..

    11448 bytes for distance
  6470400 bytes for conP
   647040 bytes for fhK
  5000000 bytes for space
```

As soon as you see the last line, you will see that the `tmp000X*` files will have been created, and hence you can stop this run by typing `ctrl+C` on the terminal where you have run such command. Please do not stop the run until the `tmp000X*` files are created. Only then, you will see that the control file you will later use to run `CODEML` has been created (i.e., `tmp0001.ctl`):

```sh
# Run from the `modsel_animals_b128/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
grep 'seqfile' */*/tmp0001.ctl | wc -l # You should get as many datasets as you have, in this case 2
```

Note that, when we ran the commands above, we were not interested in running `CODEML` or `MCMCtree`. To this end, any options that are used by `MCMCtree` and not `CODEML` (i.e., `ndata`, `usedata`, `clock`, `cleandata`, `BDparas`, `rgene_gamma`, `sigam2_gamma`, `print`, `burnin`, `sampfreq`, `nsample`, ) will be ignored, and so they do not appear in the newly generated `tmp0001.ctl` file!

---

**NOTE**: we just executed `MCMCtree` with option `usedata = 3` so that it generates the `tmp000*` files that `CODEML` will later need to calculate the branch lengths, the gradient, and the Hessian. We do this analysis in two steps given that there are restrictions in the HPC we are using that do not allow us to run `CODEML` + `MCMCtree` in one unique job within a reasonable amount of time. In addition, this analysis helps to clearly see the steps required to generate the branch lengths, the gradient, and the Hessian that `MCMCtree` will then use to approximate the likelihood calculation. In short, this is what you will be doing:

1. Run `MCMCtree` to generate the `tmp000*` files.
2. Modify the `tmp0001.ctl` file according to the dataset and the AA substitution model to be used.
3. Run `CODEML` using the `tmp000*` files so that it calculates the branch lengths, the gradient, and the Hessian; which will then be saved in a file called `rst2`.
4. Rename the `rst2` file as `in.BV` file.

---

Once all `tmp000*` files are generated for all alignments, we need to make sure that the following options have been enabled:

1. The (absolute or relative) path to the `lg.dat` file with the relevant LG matrix should be the argument of variable `aaRatefile` in the `tmp000*.ctl`.
2. Following the `Tutorial 4: Approximate likleihood with protein data`, one of the sections in the [`MCMCtree Tutorials` document](https://github.com/abacus-gene/paml/blob/master/doc/MCMCtree.Tutorials.pdf), we set the template control file to use gamma rates among sites instead of the default model, which uses no gamma rates. As we had already defined these options in the template control file, these have been already included when we previously ran `MCMCtree` to obtain the `tmp000*.ctl` files. In other words, the tmp control file should already have (i) `fix_alpha = 0` and `alpha = 0.5` (options that will enable the search to estimate the value of alpha parameter of the Gamma distribution that is used to model rate heterogeneity across sites; the starting value is $\alpha = 0.5$ as specified in the settings) and (ii) `ncatG = 4` (option that sets the number of categories required to discretise the afprementined Gamma distribution; a total of 4). We can double check everything to make sure we have made no mistake.
3. In addition, we need to make sure that option `method = 1` is enabled, which will speed up the computation of the branch lengths, gradient, and Hessian:

```sh
# Run from the `modsel_animals_b128/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
sed -i 's/method\ \=\ 0/method\ \=\ 1/' */*/tmp0001.ctl
grep 'method = 1' */*/tmp0001.ctl | wc -l # You should get as many as datasets you have
grep 'ncatG' */*/tmp0001.ctl  # You should see `ncatG = 4`
grep 'alpha' */*/tmp0001.ctl # You should see `alpha = 0.5` snd `fix_alpha = 0`
grep 'lg.dat' */*/tmp0001.ctl # You should see the absolute path in your file structure to this file
```

#### Executing `CODEML`

Now that we have the input files and control files ready to run `CODEML` for each partition, we can run `CODEML`!

We are using a template bash script to run `CODEML` (i.e., see script `pipeline_Hessian_CODEML_template.sh` in the [`scripts` directory](00_CODEML/scripts/pipeline_Hessian_CODEML_template.sh)), which flags will be replaced with relevant values when calling another bash script (`generate_job_CODEML.sh`, also saved in the [`scripts` directory](00_CODEML/scripts/generate_job_CODEML.sh)). Please note that the second bash script will edit the template bash script according to the AA substitution model chosen:

```sh
# Run from `modsel_animals_b128` dir on your HPC. Please change directories until
# you are there. Then, run the following
# commands:
home_dir=$( pwd )
cd scripts
chmod 775 *sh
# Arg1: path to where the pipeline will be executed: `pipelines_Hessian`
# Arg2: name of the main working directory, e.g., `modsel_animals_b128`. Note that only "_" or "/" are 
#       allowed to have as the names of the wd!
./generate_job_CODEML.sh $home_dir/pipelines_Hessian modsel_animals_b128
```

Now, we just need to go to the `pipelines_Hessian` directory and run the script that will have been generated using the commands above:

```sh
# Run from `modsel_animals_b128/pipelines_Hessian` dir on your HPC.
# Please change directories until
# you are there. Then, run the following
# commands.
#
# If you list the content of this directory,
# you will see the pipeline you will need 
# to execute in a bash script called
# `pipeline_Hessian.sh`
ls *
# Now, execute this bash script
chmod 775 *sh
qsub pipeline_Hessian.sh
```

Once `CODEML` finishes, we are ready to generate the `in.BV` files that we will later use when running `MCMCtree` with the approximate likelihood calculation:

```sh
# Run from dir `modsel_animals_b128/Hessian/` dir on your HPC
# Please change directories until
# you are there. Then, run the following
# commands.
home_dir=$( pwd )
for i in `seq 1 2`
do
cd $home_dir/$i
printf "[[ Accesing dir"$home_dir/$i" ]]\n"
printf "Adding Hessian and gradient for dir"$i" in \"in.BV\" file ... ...\n\n"
if [[ $i -eq 1 ]]
then
cat rst2 > $home_dir/in.BV
else
printf "\n" >> $home_dir/in.BV
cat rst2 >> $home_dir/in.BV
fi
done
```

We can also retrieve our results from the HPC for future reference:

```sh
# Run from `00_CODEML`
mkdir results_HnG
cd results_HnG
rsync -avz --copy-links <uname>@<server>:<path_to_main_wd>/modsel_animals_b128 .
# Now, extract in.BV -- we will save it inside
# `00_CODEML` for easy access in case `resilts_HnG` cannot
# be stored within this repository
cp modsel_animals_b128/Hessian/in.BV ../
```

Now, we can now proceed to run our Bayesian model selection analysis with `mcmc3r` and `MCMCtree`! If you need to refresh the concepts of Bayes factors, please read the [**Introduction to Bayes factor that we wrote for the first example of this study**](../../00_mammals/01_modsel_workflow/README.md#introduction-to-bayes-factors).

## 2. Run `mcmc3r` to select $\beta$ values

We are going to use the R package `mcmc3r` ([dos Reis et al. 2018](https://academic.oup.com/sysbio/article/67/4/594/4802240?login=false)) to (i) generate the file architecture needed to run `MCMCtree` so it can collect samples from the various power posteriors, (ii) obtain the $\beta$ values according to the method described by [Xie et al. (2011)](https://pubmed.ncbi.nlm.nih.gov/21187451/), and (iii) estimate the marginal likelihood and the Bayes factors once `MCMCtree` has finished. In this section, we will focus on the first two analyses.

We will use the `mcmc3r` function `mcmc3r::make.beta` to generate the $\beta$ values required to enable `MCMCtree` to sample from the corresponding power posterior. We have decided to summarise the samples collected across $k=128$ power posteriors to approximate the marginal likelihood (i.e., $\beta=0$ will collect samples from the prior, $\beta=1$ from the posterior, and then we will have 126 additional power posteriors with $\beta$ larger than 0 and smaller than 1 collecting samples from target distributions so that a path is traced from the prior to the posterior). Once the $\beta$ values are calculated, we will update the `MCMCtree` control files that will execute the various power posteriors so that these beta values are included. Note that we have decided to run the analyses both under the ILN and the GBM relaxed-clock models.

Once you run lines 1-55 (at the time of writing) in [our in-house R script](01_mcmc3r/scripts/Generate_BFs_inp.R), you should see the following new directories created under [`01_mcmc3r`](01_mcmc3r):

```text
01_mcmc3r
    |- control_template
    |- MCMCtree         # <-- new directory!
    |   |- BM23pUhb_GBM
    |   |   |- [1-128]
    |   |        |- mcmctree_s1_BM23pUhb_GBM.ctl    
    |   |- BM23pUhb_ILN
    |   |   |- [1-128]
    |   |        |- mcmctree_s1_BM23pUhb_ILN.ctl   
    |   |- dRetal15pUhb_GBM
    |   |   |- [1-128]
    |   |         |- mcmctree_s1_dRetal15pUhb_GBM.ctl
    |   |- dRetal15pUhb_ILN
    |       |- [1-128]
    |             |- mcmctree_s1_dRetal15pUhb_ILN.ctl
    |- scripts
```

If you open one of these control files, you will see that only one line with the corresponding beta value can be found (e.g., `BayesFactorBeta = 1e-300`). Once you run lines 56-95 (at the time of writing), you will then update the content of the control file so that the rest of options are specified (i.e., all the flags will be replaced with file names or specific integers that enable the evolutionary model under which we will run these analyses). At the same time, this part of the script will update the file names to the following:

```text
01_mcmc3r
    |- control_template
    |- MCMCtree         # <-- new directory!
    |   |- [dRetal15pUhb|BM23pUhb]
    |       |- [1-128]
    |           |- mcmctree*b[1-128].ctl  # <-- new file name including the number of the beta value being used, matching the directory
    |   
    |- scripts
```

Note that our in-house R script is using [a template control file](01_mcmc3r/control_template/mcmctree_tmp.ctl) with the same settings that [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066) specified to run their timetree inference analyses. You will see that the `treefile` setting has a flag `TREE`, which is updated when running the R script to the file name that has the diversification hypothesis being tested.

## 3. Run `MCMCtree`

We have already transferred the partitioned alignment file to the HPC and we have just generated the `in.BV` file, so we are almost ready to run `MCMCtree`!

Firstly, we need to transfer the `MCMCtree` directory we previously generated with the control files and the scripts that we will use to prepare the job arrays we will submit to the HPC and create the directory with calibrated trees:

```sh
# Run from `02_MCMCtree` on your PC
# Make sure that you have directories `scripts` and `alignments` inside `modsel_animals_b128`,
# which you should if you have followed this tutorial and ran `CODEML` before
rsync -avz --copy-links scripts/*sh <uname>@<server>:<path>/modsel_animals_b128/scripts
cd ../01_mcmc3r
rsync -avz --copy-links MCMCtree <uname>@<server>:<path>/modsel_animals_b128
# Transfer a dummy alignment that we have generated to save computations resources
##> NOTE: for details on how this alignment was generated, please read the comments
##> in the R script "Generate_dummy_aln.R", which can be found inside
##> "00_data_formatting/aln/dummy"
cd ../../00_data_formatting/aln/dummy/
rsync -avz --copy-links dummy.txt <uname>@<server>:<path>/modsel_animals_b128/alignments
```

Now, we need to create a directory where the job arrays will be running from:

```sh
# Run from `modsel_animals_b128` on the HPC
for i in dRetal15pUhb_GBM dRetal15pUhb_ILN BM23pUhb_GBM BM23pUhb_ILN
do
mkdir -p pipelines_MCMCtree/$i/
done
```

Once everything has been transferred to the cluster and the directory for the pipelines has been created, just run the following code snippet to generate the bash scripts with the job arrays to be later submitted to start the Bayesian model selection analysis:

```sh
# Run from `modsel_animals_b128`, the main directory, once you are in the cluster so that the 
# paths are already formatted
home_dir=$( pwd )
cd scripts
chmod 775 *sh
for i in dRetal15pUhb_GBM dRetal15pUhb_ILN BM23pUhb_GBM BM23pUhb_ILN
do
# arg1   Path to MCMCtree pipeline dir.
# arg2   Command to execute MCMCtree. E.g. "mcmctree", "mcmctree_4.10.
#        7", etc.
# arg3   Number of power posteriors that will be run (i.e., k=128 in our case)
# arg4   Name of the working directory, e.g., `modsel_animals_b128` in the cluster
# arg5   Type of hypothesis being tested
./generate_job_MCMCtree.sh $home_dir/pipelines_MCMCtree mcmctree 128 modsel_animals_b128 $i
done
##>NOTE: If you want, you can now transfer these pipelines to your HPC,
##>      which is why we have directory `pipelines_MCMCtree` inside `02_MCMCtree`
##>E.g.: Run the following command from your PC
##>      once you are in directory `02_MCMCtree`
##>      rsync -avz --copy-links <uname>@<server>:<path>/modsel_mammals_b256/pipelines_MCMCtree .
```

Now, we just need to go to `pipelines_MCMCtree` in the cluster and submit the job arrays for each hypothesis to run `MCMCtree` that the previous set of command will have generated!

> [!NOTE]
> If you want to download the output files we generated to continue the analysis, please download `01_animals_modsel.zip` from [FigShare repository](https://doi.org/10.6084/m9.figshare.32033958). Then, decompress this file and `03_mcmc3r.zip` and save directory `03_mcmc3r` inside `01_animals/01_modsel_workflow/03_mcmc3r`. This directory includes files that you may have already created and that you will create if you run the R script inside directory `03_mcmc3r/scripts`. You can move some of these output files elsewhere if you want to try to reproduce our results and compare output files.

## 4. Run `mcmc3r` to estimate marginal likelihoods values

Once `MCMCtree` has finished, we will retrieve the likelihood samples collected by `MCMCtree` from the power posteriors for each of the competing models and load them in an R session. We will then use `mcmc3r` for the following tasks:

* Estimate the marginal likelihood of each model
* Use the estimated marginal likelihood to calculate Bayes factors (BFs) and the posterior probabilities.
* Use the obtained results to select the best-fitting model according to this Bayesian model selection analysis.

To retrieve the `mcmc.txt` files, we will use `rsync`:

```sh
# Run from `03_mcmc3r`
mkdir HPC
cd HPC
rsync -avz --copy-links <uname>@<server>:<path>/modsel_animals_b128/MCMCtree .
rsync -avz --copy-links <uname>@<server>:<path>/modsel_animals_b128/pipelines_MCMCtree .
# Remove *sh.o* files as they are empty
rm pipelines_MCMCtree/*/*sh.o*
```

Now, we just need to run [our R in-house script](03_mcmc3r/scripts/Calculate_mlnL_BFs.R) to estimate the marginal likelihoods, compute the BFs, and the posterior probabilities!

The results show that the the **_Unconstrained_** model under the GBM relaxed-clock model is the best-fitting model. The output file [`out_animals_BFs_b128.tsv`](03_mcmc3r/HPC/MCMCtree/out_animals_BFs_b128.tsv) summarises our observations, which you can also find in the table below:

<table>
<!-- HEADER -->
<tr>
<th>Diversification model</th>
<th>logL ± S.E. (delta approximation)</th>
<th>Pr(M|D)	2.5% CI	97.5% CI</th>
</tr>
<!-- FIRST ROW -->
<tr>
<td>Unconstrained (dRetal15pUhb) GBM rates</td><td>-674.72 ± 0.09</td><td>1.00 (1.00,1.00)</td>
</tr>
<!-- SECOND ROW -->
<td>Unconstrained (dRetal15pUhb) ILN rates</td><td>-697.70 ± 0.06</td><td>0 (0,0)</td>
</tr>
<!-- THIRD ROW -->
<tr>
<td>Late Ediacaran (BM23pUhb) GBM rates</td><td>-720.24 ± 0.08</td><td>0 (0,0)</td>
</tr>
<!-- FOURTH ROW -->
<tr>
<td>Late Ediacaran (BM23pUhb) ILN rates</td><td>-764.14 ± 0.06</td><td>0 (0,0)</td>
</tr>

</table>
