# Testing contentious diversification models using a Bayesian model selection approach

## What will you find in this repository?

In this repository, you will find step-by-step guidelines to follow a protocol we have established to test different evolutionary hypotheses regarding the interpretation of the fossil record under Bayesian model selection. We have applied our protocol to assess competing hypotheses regarding two contentious nodes of the Tree of Life: the origin of animals and the origin of placental mammals.

## How is this repository organised?

The contents of this repository are categorised under the following file structure:

* [**Test case 1: origin of placental mammals**](00_mammals/README.md): we use our protocol for Bayesian model selection to test three different hypotheses concerning the origin of placental mammals when using a 72-taxa phylogeny and the molecular alignment consisting of 15,268 whole genomes assembled by [Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1):
  * **_Unconstrained_**: this diversification model does not constrain the age of node Placentalia to any specific geological event. Instead, the maximum age of the calibration density used to constrain the age of node Placentalia is informed by fossil <i>†Juramaia</i>, at ca. 162.5 Ma ([Benton et al. (2011)](https://palaeo-electronica.org/content/pdfs/424.pdf)). This interpretation of the fossil record was also followed by a study of ours, [Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1), in which, when combining genomic data and fossil evidence across the mammal tree using a Bayesian sequential subtree approach, we estimated the age of node Placentalia at ca. 80 Ma. [Our `Fossil_calibrations_justification.pdf` file](https://github.com/sabifo4/mammals_dating/blob/main/calibrations/Fossil_calibrations_justifications.pdf) contains all the details justifying the calibration densities we used across the 72-taxa phylogeny we used in our study. In this study, we have used [the same calibrated tree file](00_mammals/00_data_formatting/trees/ACetal22/72sp_ACetal22.tree) as in [Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1), which includes all the node age constraints under the _Unconstrained_ diversification model in `MCMCtree` notation. We tested this model under both **GBM** (Geometric Brownian motion, or autocorrelated) and **ILN** (independent log-normal) relaxed-clock models.
  * **_Post-K-Pg_**: based on a literal interpretation of the fossil record, [Budd and Mann (2024)](https://doi.org/10.1093/sysbio/syad057) suggest that the age of node Placentalia **cannot** pre-date the K-Pg mass extinction event (66.09 Ma), which disagrees with the mammal evolutionary timeline we inferred in [Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1) and the _Unconstrained_ diversification model. Consequently, the maximum age of node Placentalia under this model has been constrained using a hard bound. The [calibrated tree file](00_mammals/00_data_formatting/trees/BM23pUhb/72sp_BM23pUhb.tree) we have used in this study includes all the node age constraints under hypothesis _Post-K-Pg_ in `MCMCtree` notation; this is the same file that was used by [Budd and Mann (2024)](https://doi.org/10.1093/sysbio/syad057) but includes the tail probability "1e-300" instead of "0.025" to enforce a hard maximum on node Placentalia. We tested this model under both **GBM** and **ILN** relaxed-clock models.
* [**Test case 2: origin of animals**](01_animals/README.md): following the same protocol we used for test case 1, we then proceeded to test two different hypotheses regarding the origin of animals given the 2-partition alignment and corresponding phylogeny with 54 animal taxa assembled by [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066):
  * **_Unconstrained_**: this diversification model sets the maximum age of the calibration density used to constrain the age of node Metazoa  based on the Bitter springs preserved biota at ca. 833 Ma ([Benton et al. (2011)](https://palaeo-electronica.org/content/pdfs/424.pdf)). This interpretation of the fossil record was followed by [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066); their inferred evolutionary timeline placed the age of node Metazoa within the Cryogenian period. Based on this model, we now derive a stricter model in which the calibration density used to constrain the age of node Metazoa has a hard maximum bound. The [calibrated tree file](01_animals/00_data_formatting/trees/dRetal15pUhb/54sp_s1_dRetal15pUhb.tree) is the same file that was used by [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066) when they were testing their calibration strategy 1 (i.e., original file name: `strategy1.tree`; this file and others can be downloaded from their Supplementary Material [available on FigShare](https://figshare.com/articles/dataset/Uncertainty_in_the_timing_of_origin_of_animals_and_the_limits_of_precision_in_molecular_timescales/1525089?file=2230069) except from the tail probability set to "1e-300" instead of "0.001" to enforce a hard maximum on node Metazoa). We tested this model under both **GBM** and **ILN** relaxed-clock models.
  * **_Late Ediacaran_**: based on a literal interpretation of the fossil record, [Budd and Mann (2024)](https://doi.org/10.1093/sysbio/syad057) suggest that the age of node Metazoa **cannot** be older than 580 Ma based on the oldest Ediacaran assemblages ([Μatthews et al. 2021](https://doi.org/10.1130/B35646.1)) that set the boundary between the Late Ediacaran and the Cambrian period, which disagrees with the animal evolutionary timeline inferred by [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066). Based on this model, we derived a stricter one in which the calibration density used to constrain the age of node Metazoa has a hard maximum bound. The [calibrated tree file](01_animals/00_data_formatting/trees/BM23pUhb/54sp_s1_BM23pUhb.tree) we have used in this study includes all the node age constraints under hypothesis _Late Ediacaran_ in `MCMCtree` notation; this is the same file that was used by [Budd and Mann (2024)](https://doi.org/10.1093/sysbio/syad057) but includes the tail probability "1e-300" instead of "0.001" to enforce a hard maximum on node Placentalia. We tested this model under both **GBM** and **ILN** relaxed-clock models.

In order to easily navigate our analysis workflow and reproduce our results, please start with the `README.md` file under directory `00_data_formatting` that you shall find [for the mammals dataset](01_mammals/00_data_formatting/README.md) and for the [animals dataset](01_animals/00_data_formatting/README.md) -- you can start with the one you prefer! Then, you just need to follow the links you shall find throughout the document, which will eventually lead you to the next `README.md` files within the other subdirectories -- the protocol consists of step-by-step guidelines that you can follow to reproduce our results from start to end!

## What can you do with the content of this repository?

You can use the data and instructions provided in this repository to...

* ... understand how you can parse and format your input data when using `PAML` and the `mcmc3r` `R` package.
* ... understand how you can run `PAML` software and the `mcmc3r` `R` package for Bayesian model selection analysis when testing competing evolutionary hypotheses regarding node age constraints for timetree inference.
* ... understand how you can run `MCMCtree` when enabling `print = 2` so that both species divergence times and branch rates are printed on the output file (e.g., `mcmc.txt`); then summarise both results.
* ... reproduce the results we report and discuss in our study.

## Software

Before you go through the content of this repository, please make sure you have the following software installed on your PCs/HPCs:

* **`PAML`**: we used `PAML` v4.9h as this is the same version we used to analyse the data in [Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1). We could not use the older `PAML` version used by [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066), and thus we also analysed the dataset with `PAML` v4.9h. If you have installed [the latest PAML version](https://github.com/abacus-gene/paml/releases), which is available from the [`PAML` GitHub repository](https://github.com/abacus-gene/paml), please note that some values may change despite using the same seed numbers as some algorithms and part of the code may have changed from version to version. If you need to install `PAML`, please follow [the instructions on the PAML GitHub repository according to your OS](https://github.com/abacus-gene/paml/wiki/Installation).

* **`R`** and **`RStudio`**: please download [`R`](https://cran.r-project.org/) and [`RStudio`](https://posit.co/download/rstudio-desktop/) as you will need them to run most of our in-house `R` scripts. If you are a Windows user, please make sure that you have the correct version of `RTools` installed, which will allow you to install packages from the source code if required. For instance, if you have R v4.1.2, then installing `RTools4.0` shall be fine. If you have another `R` version installed on your PC, please check whether you need to install a later version (e.g., `RTools 4.4`, `RTools 4.5`, etc.). For more information on which version you should download, [please go to the CRAN website by following this link and download the version you need](https://cran.r-project.org/bin/windows/Rtools/).

    Before you proceed, however, please make sure that you install the following packages too:

    ```R
    # Run from the R console in RStudio
    # Check that you have at least R v4.1.2 (scripts were first written when working
    # with this version)
    # We have run our scripts with R v4.4.2
    version$version.string
    # Now, install the packages we will be using
    # Note that it may take a while if you have not 
    # installed all these software before
    install.packages( c( 'rstudioapi', 'ape', 'phytools', 'sn', 'stringr', 'rstan' ), dep = TRUE )
    ## NOTE: If you are a Windows user and see the message
    ## "Do you want to install from sources the 
    ## packages which need compilarion?", please make sure that
    ## you have installed `RTools`
    ## as aforementioned before proceeding.
    ## The versions we used for each of the packages aforementioned
    ## are the following:
    ##   rstudioapi: v0.17.1 (at least use v0.14)
    ##   ape: v5.8.1 (at least use v5.7.1)
    ##   phytools: v2.4.4 (at least use v1.5.1)
    ##   sn: v2.1.1 (at least use v2.1.1)
    ##   stringr: v1.5.1 (at least use v1.5.0)
    ##   rstan: v2.32.6 (at least use v2.21.7)
    ## Now, instal the `mcmc3r` R package with vignettes. If you do not have `devtools` installed,
    ## please install this package before running the next command
    devtools::install_github( "dosreislab/mcmc3r", build_vignettes = TRUE )
    ```

* Your preferred graphical interface for phylogeny visualisation and manipulation (e.g., **`FigTree`**, **`TreeViewer`**, etc.): you can use your preferred graphical interface to display tree topologies with/without branch lengths and with/without additional labels. You can then decide what you want to be displayed by selecting the options that you require for that to happen. For instance, you can download [`TreeViewer` from their website](https://treeviewer.org/) ((Bianchini & Sánchez-Baracaldo 2024)[https://doi.org/10.1002/ece3.10873]) and [the latest stable pre-compiled binaries from the `FigTree` GitHub repository](https://github.com/rambaut/figtree/releases).

* **`Tracer`**: you can use this graphical interface to visually assess the MCMCs you have run during your analyses (e.g., chain efficiency, chain convergence, autocorrelation, etc.). You can [download the latest stable pre-compiled binaries from the `Tracer` GitHub repository](https://github.com/beast-dev/tracer/releases/).

## Data analysis

If you have gone through the previous sections and have a clear understanding of the dataset we used, the workflow we followed (which you shall follow if you want to reproduce our analyses!), and have installed the required software... Then you are ready to go!

If you want to start reproducing our analyses for our test case 1 (testing diversification models for the origin of placental mammals), please [follow this link to get started](00_mammals/00_data_formatting/README.md). Once you go through all the details surrounding how our input data were formatted, you will find additional links in the `README.md` file that will guide you through the next steps! If you want to start reproducing our analyses for test case 2 (testing diversification models for the origin of animals), then please [follow this link and get started with our second set of analyses](01_animals/00_data_formatting/README.md)!

Happy reproducible timetree inference! :smile:

## Data reproducibility

The step-by-step guidelines you shall find in this repository will help you reproduce our analyses! Nevertheless, there are some large files that we have been unable to commit to this repository due to size limit. Instead, you shall find them in our [FigShare archive](https://figshare.com/s/badcc59a192fa4fd0722).

If you are re-running our analyses, please download the relevant files so that you can compare them to your output files. In addition, some of the large files may be required to generate some of the figures and tables included in our manuscript. If you also want to reproduce, please make sure that you save them in the directory indicated in the R script or modify such path to where you have now saved them.

## Citations and Contact

If you use our protocol to run a Bayesian model selection analysis, please **cite**...

* ... our study: [Álvarez-Carretero et al. (2026). Mol Biol Evol msag211](https://doi.org/10.1093/molbev/msag211).
* ... the `mcmc3r` R package: [dos Reis et al. (2018). Syst Biol 67:594-614](https://doi.org/10.1093/sysbio/syy001).
* ... the `PAML` package: [Yang (2007). Mol Biol Evol 24:1586–1591](https://doi.org/10.1093/molbev/msm088).

If you find any **technical issues** or any **formatting problems** in this repository...

* ... please raise an issue so we can fix it as soon as possible! :muscle:

If you have any **other questions** not related with technical or formatting problems or you would like to **discuss our study**...

* ... please contact [@sabifo4](https://github.com/sabifo4/), the developer and manager of this repository, <a href="mailto://sandra.ac93@gmail.com">via e-mail</a>.
