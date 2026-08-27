# Bayesian model selection

Given that we already have the calibrated tree files for each hypothesis, a dummy alignment to check that the same taxa names are used, and the `in.BV` file with the branch lengths, the gradient, and the Hessian required to approximate the likelihood calculation... We have everything we need to run `mcmc3r` and `MCMCtree`!

Before proceeding to Bayesian model selection, however, please take some time to familiarise with concepts such as "marginal likelihood", "power posterior", and "Bayes factors" in the next section [`Introduction to Bayes factors`](README.md#introduction-to-bayes-factors), which are key to understand how Bayesian model selection can help us test competing hypotheses and find the best-fitting one for the dataset being analysed.

## Introduction to Bayes factors

When carrying out a **Bayesian model selection** analysis, we need to calculate **Bayes factors** so that we can assess how fitting the contested hypotheses are given our dataset. Consequently, we will first need to estimate the **marginal likelihood**, also known as the normalising constant. In Bayesian timetree inference, we always avoid calculating this term during MCMC -- the marginal likelihood is a multidimensional integral the calculation of which tends to be computationally unfeasible! For instance, when using the Metropolis-Hastings algorithm ([Metropolis et al., 1953](https://aip.scitation.org/doi/10.1063/1.1699114) and [Hastings 1970](https://academic.oup.com/biomet/article-abstract/57/1/97/284580)), the marginal likelihood cancels out when the newly proposed value is being evaluated (i.e., should this value be kept or discarded?). Without using mathematical notation, this step can be written as follows:

$$\frac{\textrm{unnormalised posterior'}}{\textrm{unnormalised posterior}}=\frac{\frac{\textrm{prior'}\times \textrm{ likelihood'}}{\textrm{marginal likelihood'}}}{\frac{\textrm{prior}\times \textrm{ likelihood}}{\textrm{marginal likelihood}}}=\frac{\textrm{prior'}\times \textrm{ likelihood'}}{\textrm{prior}\times \textrm{ likelihood}}$$

In the notation above, you can see that the terms that end with `'` correspond to those calculated with the newly proposed value(s) for the parameter(s) of interest in a given chain iteration. Those without this symbol had already been calculated in the previous chain iteration with the value(s) that had been accepted for the parameter(s) of interest. Both marginal likelihood terms must integrate to 1, and thus they cancel out. Now, the main question is: how do we **estimate** this normalising constant that we have been avoiding to deal with during Bayesian timetree inference?

There are various **approximation methods** that can be used to estimate the marginal likelihood, being thermodynamic integration (TI, [Gelmand and Meng 1998](https://www.jstor.org/stable/2676756); [Lartillot and Philippe 2006](https://doi.org/10.1080/10635150500433722); [Lepage et al. 2008](https://doi.org/10.1093/molbev/msm193); [Friel and Pettitt 2008](https://academic.oup.com/jrsssb/article/70/3/589/7109555)) and the stepping-stone (SS, [Xie et al. 2011](https://doi.org/10.1093/sysbio/syq085)) approach those most widely used in the field of phylogenetics. In this study, we will be using the **SS approach to estimate the marginal likelihood**. You can think of this method as trying to build a path from the prior to the posterior by running $k$ MCMCs that will collect samples for our parameter(s) of interest from a target distribution: the so-called "**power posterior**". When sampling from the power posteriors, we will focus on the sampled **likelihood** values: we need those to estimate our marginal likelihood! Without using mathematical notation, the power posterior can be expressed as follows:

$$\textrm{power posterior}_{\beta}\propto \textrm{prior}\times \textrm{likelihood}^{\beta}$$

We can see that the likelihood term is raised to the power of $\beta$. Without knowing much  about what $\beta$ values one should use, we can guess the outcome of using $\beta=0$ and $\beta=1$: anything raised to the power of 0 becomes 1 and anything raised to the power of 1 does not change. Consequently, when $\beta=0$, the previous equation can be rewritten as follows:

$$\textrm{power posterior}_{\beta = 0}\propto \textrm{prior}\times \textrm{likelihood}^{\beta = 0}\propto \textrm{prior}$$

In the equation above, we can see that the **target distribution** (i.e., the power posterior when $\beta=0$) from which samples will be collected during the MCMC is proportional to the **prior distribution**. When $\beta=1$, however, we obtain the equation that we are used to seeing in Bayesian timetree inference (i.e., prior times the likelihood):

$$\textrm{power posterior}_{\beta = 1}\propto \textrm{prior}\times \textrm{likelihood}^{\beta = 1}\propto \textrm{prior}\times \textrm{likelihood}$$

If we were now to select several values of $\beta$ that range between 0 and 1, we would be collecting samples from different target distributions (i.e., power posteriors, one for each selected $\beta$ value) that bridge the prior ($\beta = 0$) and the posterior ($\beta = 1$). In other words, we are using "stepping stones" to help build a path between the prior and the posterior. It is noteworthy that the selection of $\beta$ values is not trivial. Following [Xie et al. 2011](https://doi.org/10.1093/sysbio/syq085), most of the $\beta$ values are sampled from evenly spaced quantiles of a Beta($\alpha$,1.0). The authors suggest a value of $\alpha=0.3$, which yields a positively skewed distribution so that half of the sampled $\beta$ values are less than 0.1 (and thus boost sampling in the unstable area of the parameter space).

Now that the concept of Bayesian model selection has been introduced, let's focus on the protocol we have developed as part of our study!

## 1. Run `mcmc3r` to select $\beta$ values

We are going to use the R package `mcmc3r` ([dos Reis et al. (2018)](https://academic.oup.com/sysbio/article/67/4/594/4802240?login=false)) to (i) generate the file architecture needed to run `MCMCtree` so it can collect samples from the various power posteriors, (ii) obtain the $\beta$ values according to the method described by [Xie et al. (2011)](https://pubmed.ncbi.nlm.nih.gov/21187451/), and (iii) estimate the marginal likelihood and the Bayes factors once `MCMCtree` has finished. In this section, we will focus on the first two analyses.

We will use the `mcmc3r` function `mcmc3r::make.beta` to generate the $\beta$ values required to enable `MCMCtree` to sample from the corresponding power posterior. We have decided to summarise the samples collected across $k=256$ power posteriors to approximate the marginal likelihood (i.e., $\beta=0$ will collect samples from the prior, $\beta=1$ from the posterior, and then we will have 126 additional power posteriors with $\beta$ larger than 0 and smaller than 1 from which samples will be collected with the aim to trace a path from the prior to the posterior). Once the $\beta$ values are calculated, we will update the `MCMCtree` control files that will execute the various power posteriors so that these beta values are included.

If you want to generate the input files inside directory [`MCMCtree`](00_mcmc3r/MCMCtree/), make sure that you copy the current directory `MCMCtree` with our results and/or rename it differently (e.g., `MCMCtree_orig`). Once you have done that, **if you run lines 1-56** (at the time of writing) in [our in-house R script](00_mcmc3r/scripts/Generate_BFs_inp.R), you should see that a new directory `MCMCtree` is created inside directory [`00_mcmc3r`](00_mcmc3r) with the following subdirectories:

```text
00_mcmc3r
    |- control_template
    |- MCMCtree         # <-- new directory!
    |   |- ACetal22_GBM
    |   |   |- [1-256]
    |   |        |- mcmctree_ACetal22_GBM_b.ctl
    |   |- ACetal22_ILN
    |   |   |- [1-256]
    |   |        |- mcmctree_ACetal22_ILN.ctl
    |   |- BM23pUhb_GBM
    |   |   |- [1-256]
    |   |        |- mcmctree_BM23_GBM.ctl
    |   |- BM23pUhb_ILN
    |      |- [1-256]
    |            |- mcmctree_BM23_ILN.ctl 
    |   
    |- scripts
```

If you open one of the control files that you shall find inside the subdirectories shown above, you will see that only one line with the corresponding beta value can be found (e.g., `BayesFactorBeta = 1e-300` for directory `1` inside all the directories aforementioned). **Once you run lines 61-92** (at the time of writing), you will then update the content of the control file so that the rest of options are included (i.e., the same options we used in our first study). At the same time, this part of the script will update the file names to the following:

```text
00_mcmc3r
    |- control_template
    |- MCMCtree         # <-- new directory!
    |   |- <name_directory>
    |       |- [1-256]
    |            |- mcmctree_<name_directory>_b[1-256].ctl  # <-- new file name including the number of the beta value being used, matching the directory
    |- scripts
```

Note that our in-house R script is using [a template control file](00_mcmc3r/control_template/mcmctree_tmp.ctl) with the same settings we specified to run the analysis with this backbone tree in [Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1) (i.e., any of the `mcmctree.ctl` control files inside `run[1-4]/mcmctree_GBM/` directories you shall find when you decompress [the following file available in our Dropbox archive](https://www.dropbox.com/s/53mdfyc47hukkrh/SeqBayesS1_MCMCtree_mainT2_posterior_newchrono.zip?dl=0)). You will see that the `treefile` and the `clock` settings have two flags, `TREE` and `RCLOCK`, respectively; which are replaced with the file name that has the tree hypothesis that will be evaluated (i.e., `72sp_ACetal22.tree`, `72sp_BM23.tree`, or `72sp_BM23pUhb.tree`) when running the R script. A summary of the settings, while given in our study, can be found below:

* Random seed number (`seed`): by setting this variable to `-1`, a random seed number will be generated for an analysis based on the time stamp.
* Path to the sequence file (`seqfile`): this variable specifies the path to the alignment file. Given that we will have a soft link to our alignment file (`dummy.txt`) in the same directory where this control file is, we just need to type the file name. Nevertheless, note that absolute and relative paths can also be used.
* Path to the calibrated tree file (`treefile`): this variable specifies the path to the calibrated file. You can see that the name of the tree file corresponding to the hypothesis being tested will be here (i.e., `72sp_ACetal22.tree` and `72sp_BM23.treepUhb`). As mentioned above, there will be a soft link to this tree file in the same directory where this control file will be, and thus only the file name is required. Nevertheless, note that absolute and relative paths can also be used.
* Path to the output file (`outfile`): this variable specifies the path to the alignment file. The output file is to be saved in the same directory where this control file is, and thus no relative paths are needed.
* Number of partitions (`ndata`): the partitioned input file has four alignment blocks, and thus this variable is set to `4`.
* Sequence type (`seqtype`): the alignment file has nucleotide sequences, and thus this variable is set to `0`.
* Type of method to calculate the likelihood (`usedata`): this variable is set to `0` when sampling from the prior (i.e., the alignment will not be used, and thus the likelihood will equal to 1) and to `2` when sampling from the posterior and using the approximate likelihood calculation (i.e., the `in.BV` file will be read and used). Note that there are two additional options that can be enabled when `usedata = 2`: (i) the path to the `in.BV` file (useful if the file name is not `in.BV`) and (ii) the type of approximation being used. In our case, we have not changed the name of the `in.BV` file and we will have its soft link saved in the same directory where the control file shall be. Consequently, there is no need to include the path to this file after option `2` (same as in [Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1)), although it would have been possible to include `./in.BV` as in `usedata = 2 ./in.BV` if users want to include this option. We will be using the arcsin approximation, the default option, and thus there is no need to include a third argument. If the approximation approach is to be changed, then an integer corresponding to the chosen approximation needs to be included (e.g., 0: NT, 1: SQRT, 2: LOG, 3: ARCSIN (default); `usedata = 2 1` would be using the square root approximation).
* Type of clock model (`clock`): when sampling from the prior, this option is not really used, and thus any of the three possible options (i.e., `1`: strict clock model, `2`: independent-rates log-normal model, `3`: geometric Brownian motion [i.e., autocorrelated-rates] model) will be ignored. We have run the analyses under the autocorrelated-rates model or Geometric Brownian Motion (GBM, `clock = 3`, ([Thorne et al. (1998)](https://pubmed.ncbi.nlm.nih.gov/9866200/), [Yang and Rannala (2006)](https://pubmed.ncbi.nlm.nih.gov/16177230/))) and the independent-rates log-normal model (ILN, `clock = 2`). Note that the former was the model which we showed in [Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1) to better fit our data (see Fig. 2a and 2b in [Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1)).
* Dealing with ambiguity data (`cleandata`): we set this variable to `0` as our sequences had already been processed and there was no ambiguity data to be removed. All the filtering and processing steps are described in the Methods section in [Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1) and detailed in our [GitHub repository](https://github.com/sabifo4/mammals_dating/).
* Prior on divergence times (`BDparas`): the prior on divergence times that `MCMCtree` uses is based on the birth-death process with species sampling of [Yang and Rannala (1997)](https://pubmed.ncbi.nlm.nih.gov/9214744/). We chose the values for the birth ($\lambda$) and death ($\mu$) parameters of the model to be equal to 1, while the sampling frequency was $\rho=0.1$. These parameters yield a uniform kernel density that has an approximate uniform distribution. This kernel density will specify the prior distribution for the uncalibrated nodes, and thus will be used together with the user-specified priors for the calibrated nodes to build the joint prior. The resulting joint prior is the effective prior that `MCMCtree` will use when sampling from the posterior. Users must examine this effective prior to make sure it is sensible and that it represents their knowledge of the species and the relevant fossil record. If need be, they may have to change the fossil calibrations that the effective prior looks reasonable. One can retrieve such effective priors when running `usedata = 0`, which can then be plotted against the user-specified priors to check whether major discrepancies are to be found due to possible truncation issues. For more information on the matter, please see [dos Reis et al. (2015)](https://pubmed.ncbi.nlm.nih.gov/26603774/).
* Rate prior (`rgene_gamma`): we assumed an approximate evolutionary rate of 0.05 substitutions per site per 100Myr (mean rate =~ 1e-10 substitutions per site per year). We used a diffuse gamma prior with $\alpha=2$ and $\beta=40$ (i.e., $\text{mean rate}=\frac{\alpha}{\beta}=\frac{2}{40}=0.05$) to represent our uncertainty regarding the average substitution rate of nuclear genes in mammals.
* Prior on the variation of the rate (`sigma2_gamma`): given that our phylogeny has deep divergences, the clock may have been violated, and thus we have fixed a mean for the `sigma2` parameter (i.e., clock variation) as 0.1 using a gamma prior with $\alpha=1$ and $\beta=10$ (i.e., $\text{mean sigma2}=\frac{\alpha}{\beta}=\frac{1}{10}=0.1$).
* MCMC burn-in (`burnin`): we chose to discard the samples collected during the first **150,000** iterations as burn-in following the same settings as in [Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1).
* Sampling frequency during the MCMC (`sampfreq`): we will sample every **400** iterations as in [Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1).
* Number of samples to be collected during the MCMC (`nsample`): we chose to collect a total of 20,000 samples following the same settings as in [Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1). Therefore, we ran a total of $150000+5000\times 20000=100150000$ iterations (burnin + sampling_frequency*num_samples).
* Type of data that will be written out (`print`): we specified value `1` so that the samples collected for all model parameters are written out in file `mcmc.txt`.

## 2. Run `MCMCtree`

Now, we need to create the file structure to run our analyses in a HPC. We have provided this file structure here too so that you can also transfer it to your HPC in case you want to reproduce our analyses. If you want to generate this file structure by yourself, make sure that you delete or move the directories we provide inside `01_MCMCtree/HPC` before running the following code:

```sh
# Run from `01_MCMCtree`
## NOTE: Directory `scripts` will be already inside this directory, which contains
## the bash scripts that we use to then generate the bash scripts to run job arrays
mkdir HPC
cd HPC
# Copy alignment
mkdir aln
cp ../../../00_data_formatting/aln/dummy.txt aln
# Move `in.BV` to directory Hessian
mkdir Hessian
cp ../../../00_data_formatting/aln/in.BV Hessian
# Copy `MCMCtree` dir previously generated with `mcmc3r` output
cp -R ../../00_mcmc3r/MCMCtree/ .
# Copy trees
name_hypotheses=`ls MCMCtree`
for i in $name_hypotheses
do
mkdir -p trees/$i
done
for i in trees/ACetal22_*
do
cp ../../../00_data_formatting/trees/ACetal22/* $i/
done
for i in trees/BM23pUhb_*
do
cp ../../../00_data_formatting/trees/BM23pUhb/* $i/
done
# Copy scripts
cp -R ../scripts .
```

Now, everything is ready to run `MCMCtree` in a job array manner. To create the file structure to save the output log files when running the job arrays, you can run the following commands:

```sh
# Run from `01_MCMCtree/HPC`
for i in ACetal22_GBM ACetal22_ILN BM23pUhb_ILN BM23pUhb_GBM
do
mkdir -p pipelines_MCMCtree/$i/
done
```

We have generated a template bash script, [`pipeline_MCMCtree_template.sh`](01_MCMCtree/scripts/pipeline_MCMCtree_template.sh) with flags that will be replaced when running the bash script [`generate_job_MCMCtree.sh`](01_MCMCtree/scripts/generate_job_MCMCtree.sh). These can then be transferred to the cluster (e.g., `rsync`), where they will be submitted. E.g.:

```sh
# Run from `01_MCMCtree`
# Make sure that you have created a directory called `modsel_mammals_b256`
# in your HPC scratch directory or similar so that these files can
# be successfully copied!
rsync -avz --copy-links HPC/* <uname>@<server>:<path>/modsel_mammals_b256/
```

Once everything has been transferred to the cluster, just do the following to generate the bash scripts with the job arrays to be later submitted to start the Bayesian model selection analysis:

```sh
# Run from `modsel_mammals_b256`, the main directory, once you are in the cluster so that the 
# paths are already formatted
home_dir=$( pwd )
cd scripts
chmod 775 *sh
for i in ACetal22_GBM ACetal22_ILN BM23pUhb_ILN BM23pUhb_GBM
do
# arg1   Path to MCMCtree pipeline dir.
# arg2   Command to execute MCMCtree. E.g. "mcmctree", "mcmctree_4.10.7", etc.
# arg3   Number of power posteriors that will be run (i.e., k=256 in our case)
# arg4   Name of the working directory, e.g., `modsel_mammals_b256` in the cluster
# arg5   Type of hypothesis being tested (i.e., ACetal22_GBM, BM23_GBM, ACetal22_ILN, BM23_ILN )
./generate_job_MCMCtree.sh $home_dir/pipelines_MCMCtree mcmctree 256 modsel_mammals_b256 $i
done
##>NOTE: If you want, you can now transfer these pipelines to your HPC
##>E.g.: Run the following command from your PC
##>      once you are in directory `01_MCMCtree/HPC`
##>      rsync -avz --copy-links <uname>@<server>:<path>/modsel_mammals_b256/pipelines_MCMCtree .
```

Now, we just need to load the job arrays for each hypothesis to run `MCMCtree`!

> [!NOTE]
> If you want to download the output files we generated to continue the analysis, please download `00_mammals_modsel.zip` from [FigShare repository](https://doi.org/10.6084/m9.figshare.32033958). Then, decompress this file and `02_mcmc3r.zip` and save directory `02_mcmc3r` inside `00_mammals/01_modsel_workflow/02_mcmc3r`. This directory includes files that you may have already created and that you will create if you run the R script inside directory `02_mcmc3r/scripts`. You can move some of these output files elsewhere if you want to try to reproduce our results and compare output files.

## 3. Run `mcmc3r` to estimate marginal likelihood values

Once `MCMCtree` has finished, we will retrieve the likelihood samples collected by `MCMCtree` from the power posteriors for each of the analyses and load them in an R session. We will then use `mcmc3r` for the following tasks:

* Estimate the marginal likelihood of each model
* Use the estimated marginal likelihood to calculate Bayes factors (BFs) and the posterior probabilities.
* Use the obtained results to select the best-fitting model according to this Bayesian model selection analysis.

To retrieve the `mcmc.txt` files, we will use `rsync`:

```sh
# Run from `02_mcmc3r`
mkdir HPC
cd HPC
rsync -avz --copy-links <uname>@<server>:<path>/modsel_mammals_b256/MCMCtree .
rsync -avz --copy-links <uname>@<server>:<path>/modsel_mammals_b256/pipelines_MCMCtree .
# Once transferred, remove the `s.oh*` files as they are empty
rm pipelines_MCMCtree/*/*sh.o*
```

Now, we just need to run [our R in-house script](02_mcmc3r/scripts/Calculate_mlnL_BFs.R) to estimate the marginal likelihoods, compute the BFs, and the posterior probabilities!

The results show that the **_Unconstrained_** model under GBM fits the data better than any of the other competing diversification models tested, which are summarised in the output file [`out_mammals_BFs_b256.R`](02_mcmc3r/HPC/MCMCtree/out_mammals_BFs_b256.tsv) and displayed in the table below:

<table>
<!-- HEADER -->
<tr>
<th>Diversification model</th>
<th>logL ± S.E. (delta approximation)</th>
<th>Pr(M|D)	2.5% CI	97.5% CI</th>
</tr>
<!-- FIRST ROW -->
<tr>
<td>Unconstrained (ACetal22) GBM rates</td><td>-2509.80 ± 0.16</td><td>1 (1, 1)</td>
</tr>
<!-- SECOND ROW -->
<tr>
<td>Post-K-Pg (same as BM24 except<br>for hard maximum age to constrain node Placentalia)<br>GBM rates</td><td>-2609.06 ± 1.56</td><td>7.80 × 10<sup>–44</sup> (4.04 × 10<sup>–45</sup>, 1.61 × 10<sup>–42</sup>)</td>
</tr>
<!-- THIRD ROW -->
<tr>
<td>Unconstrained (ACetal22) ILN rates</td><td>-2612.74 ± 0.12</td><td>1.96 × 10<sup>–45</sup> (1.33 × 10<sup>–45</sup>, 2.93 × 10<sup>–45</sup>)</td>
</tr>
<!-- FOURTH ROW -->
<tr>
<td>Post-K-Pg (same as BM24 except<br>for hard maximum age to constrain node Placentalia)<br>ILN rates</td><td>-2774.48 ± 1.92</td><td>1.12 × 10<sup>–115</sup> (2.81 × 10<sup>–117</sup>, 4.86 × 10<sup>–114</sup>)</td>
</tr>

</table>
