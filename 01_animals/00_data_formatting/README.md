# Input files

The input files used for this analysis are the following:

* [**Tree files following "dRetal15" calibrations**](trees/dRetal15/): tree files that were used in [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066) to constrain the node ages of their 54-taxa animal tree. These files can be retrieved from the [FigShare archive](https://dx.doi.org/10.6084/m9.figshare.1525089) that was uploaded as a supplementary material for their study. Specifically, we are using the tree files used to evaluate two different strategies:
  * [**`54sp_s1_dRetal15.tree`**](trees/dRetal15/54sp_s1_dRetal15.tree): we downloaded the file called `strategy1.tree` from their [FigShare archive](https://dx.doi.org/10.6084/m9.figshare.1525089) and renamed it to `54sp_s1_dRetal15.tree`. This input tree file is used to evaluate the fitness of calibration strategy 1 following [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066).
  * [**`54sp_s2_dRetal15.tree`**](trees/dRetal15/54sp_s2_dRetal15.tree): we downloaded the file called `strategy2.tree` from their [FigShare archive](https://dx.doi.org/10.6084/m9.figshare.1525089) and renamed it to `54sp_s2_dRetal15.tree`. This input tree file is used to evaluate the fitness of calibration strategy 2 following [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066).

* [**Tree files following "BM23" calibrations**](trees/BM23/): the authors did not make their input files available as part of [their published point of view](https://doi.org/10.1093/sysbio/syad057), but they kindly shared them with us when we requested them. Following their analysis, the files we are using as input tree files are the following:
  * [**`54sp_s1_BM23.tree`**](trees/BM23/54sp_s1_BM23.tree): the node age constraints used in this input tree file are the same as in [`54sp_s1_dRetal15.tree`](trees/dRetal15/54sp_s1_dRetal15.tree) except for the maximum age of node Metazoa. [Budd and Mann (2023)](https://doi.org/10.1093/sysbio/syad057) used this input calibrated tree file in their inference analyses. We have renamed this file as [`54sp_s1_BM23.tree`](trees/dRetal15/54sp_s1_dRetal15.tree), and used it to compare the fitness of the different constraint for Metazoa against the one used by [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066).
  * [**`54sp_s2_BM23.tree`**](trees/BM23/54sp_s2_BM23.tree): in this input calibrated tree file, the node age constraints are the same as in [`54sp_s2_dRetal15.tree`](trees/dRetal15/54sp_s2_dRetal15.tree) except for the maximum age to constrain node Metazoa, which follows [Budd and Mann (2023)](https://doi.org/10.1093/sysbio/syad057). We named this file [`54sp_s2_BM23.tree`](trees/BM23/54sp_s2_BM23.tree), and used it to compare the fitness of the different constraint for Metazoa against the one used by [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066).

* [**Alignment files**](aln): we are using the same partitioned alignment assembled in [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066) and used in [Budd and Mann (2023)](https://doi.org/10.1093/sysbio/syad057) for timetree inference. In addition, we have generated a dummy alignment file to use when sampling from the multiple power posteriors to save computational resources. These files are the following:
  * [**Partitioned alignment file, `aln-2P.phy`**](aln/aln-2P.phy): the same partitioned alignment file used in [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066), which can be retrieved from the [FigShare archive](https://dx.doi.org/10.6084/m9.figshare.1525089) that was uploaded as a supplementary material for their study.
  * [**Alignment files with each partition separately, `aln_P1.phy` and `aln_P2.phy`**](aln/): we need the branch lengths, the gradient, and the Hessian for each partition to run `MCMCtree` with the approximate likelihood calculation ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)). Therefore, we extracted the individual partitions from [their partitioned alignment file](aln/aln-2P.phy) using the following code snippet:

  ```sh
  # Run from `00_data_formatting/aln`
  pattern="^[0-9]* "
  count=0
  while IFS= read -r line
  do
  if [[ $line =~ $pattern ]]
  then
  count=$(( count + 1 ))
  echo Header $count
  echo "$line" > aln_part$count".phy"
  else
  echo "$line" >> aln_part$count".phy"
  fi
  done < aln-2P.phy
  ```

  * [**Dummy alignment file, `dummy.txt`**](aln/dummy/dummy.txt): this alignment file can be used to save computational resources when the `in.BV` file (i.e., branch lengths, gradient, and Hessian) is available. In that way, less memory is allocated to calculate site patterns from the alignment -- all the information concerning the alignment that will be used when approximating the likelihood calculation during the MCMC is saved in the `in.BV` file!

## Reproducibility

In order for users to reproduce the input calibrated tree files we have used in this study, we have written an in-house R script: [`Include_calibrations.R`](trees/uncalib/scripts/Include_calibrations.R). The input files required to run this script are the following:

* [`uncalib_animal.tree`, tree topology without labels](trees/uncalib/uncalib_animal.tree): this uncalibrated tree is used to run `CODEML` (i.e., this program is used to analyse amino acid data, `BASEML` is used for nucleotide data). In order to generate this file, we are going to pick one of the input calibrated tree files from [dos Reis et al. (2015)](i.e., [strategy 1](trees/dRetal15/54sp_s1_dRetal15.tree)) and run the following code snippet:

  ```sh
  # Run from `00_data_formatting/trees/uncalib`
  # Run only if you want to reproduce our results as we already
  # provide users with this uncalibrated tree
  mkdir original
  mv uncalib_animal.tree original/
  cp ../dRetal15/54sp_s1_dRetal15.tree uncalib_animal.tree
  sed -i "s/'//g" uncalib_animal.tree
  sed -i 's/B([0-9]*\.[0-9]*,[0-9]*\.[0-9]*,[0-9]*\.[0-9]*,[0-9]*\.[0-9]*)//g' uncalib_animal.
  ```

* [`uncalib_animal_flags.tree`, tree topology with labels](trees/uncalib/uncalib_animal_flags.tree): we generated a copy of one of the calibrated tree files and manually incorporated the flags to label the 34 nodes which ages are being constrained by fossil calibrations. These flags are in the format of `[CALIBRATION]`, where `CALIBRATION` can be replaced with the name of the node being calibrated (e.g., METAZOA, PROTOSTOMIA, etc.). These labels are replaced with the corresponding `MCMCtree` notation to enable the constraints by [our in-house R script](trees/uncalib/scripts/Include_calibrations.R) (see below). It is important to note that the labels need to be written within square brackets and that the last line of this file must be a blank line!
* [Matching text files](trees/uncalib/add_calibs/): the following matching text files have been generated:
  * [`Calib_converter_s1_dRetal15.txt`](trees/uncalib/add_calibs/Calib_converter_s1_dRetal15.txt)
  * [`Calib_converter_s2_dRetal15.txt`](trees/uncalib/add_calibs/Calib_converter_s2_dRetal15.txt)
  * [`Calib_converter_s1_BM23.txt`](trees/uncalib/add_calibs/Calib_converter_s1_BM23.txt)
  * [`Calib_converter_s2_BM23.txt`](trees/uncalib/add_calibs/Calib_converter_s2_BM23.txt)
 These files include the labels used to flag the nodes to be calibrated in the tree topologies above (e.g., MAMMALIA, PLACENTALIA, etc.), followed by a pipe delimiter |, and then the calibration in `MCMCtree` notation within single quotation marks (e.g., soft-bound calibrations:`'B(min_age,max_age)'`; skew-normal densities:`SN(xi,omega,alpha)`). The [`Include_calibrations.R`](trees/uncalib/scripts/Include_calibrations.R) R script will use these text files to replace the labels in the uncalibrated tree topology (i.e., [`uncalib_animal_flags.tree`](trees/uncalib/uncalib_animal_flags.tree)) with the corresponding `MCMCtree` notation. It is important that the last line of these files is a blank line or, otherwise, the R script will not work!

Once you run [the `Include_calibrations.R` R script](trees/uncalib/scripts/Include_calibrations.R), you will see a newly created directory inside `scripts` (i.e., `scripts/reproduce`) where the four calibrated tree files will be output:

* [**`54sp_s1_dRetal15.tree`**](trees/dRetal15/54sp_s1_dRetal15.tree): the calibrated tree file following strategy 1 in [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066).
* [**`54sp_s1_BM23.tree`**](trees/BM23/54sp_s1_BM23.tree): the calibrated tree file used in [Budd and Mann (2023)](https://doi.org/10.1093/sysbio/syad057), same constraints as in [`54sp_s1_dRetal15.tree`](trees/dRetal15/54sp_s1_dRetal15.tree) except for node Metazoa.
* [**`54sp_s2_dRetal15.tree`**](trees/dRetal15/54sp_s2_dRetal15.tree): the calibrated tree file following strategy 2 in [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066).
* [**`54sp_s2_BM23.tree`**](trees/BM23/54sp_s2_dRetal15.tree): the calibrated tree file used in [Budd and Mann (2023)](https://doi.org/10.1093/sysbio/syad057), same constraints as in [`54sp_s2_dRetal15.tree`](trees/dRetal15/54sp_s2_dRetal15.tree) except for node Metazoa.

> **NOTE**: I have not yet run this part as I am to confirm that the flags I used in [`uncalib_animal_flags.tree`](trees/uncalib/uncalib_animal_flags.tree) are indeed the names I should be using. I will need to verify this as well with the matching text files in [the `add_calibs` directory](trees/uncalib/add_calibs/).
