# Input files

The input files used for this analysis are the following:

* [**Tree file following "dRetal15" calibrations**](trees/dRetal15pUhb/): tree file that was used in [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066) to constrain the node ages of their 54-taxa animal tree. This file (and other tree files they used) can be retrieved from the [FigShare archive](https://dx.doi.org/10.6084/m9.figshare.1525089) that was uploaded as a supplementary material for their study. Specifically, we are using the following tree file:
  * [**`54sp_s1_dRetal15pUhb.tree`**](trees/dRetal15pUhb/54sp_s1_dRetal15pUhb.tree): we downloaded the file called `strategy1.tree` from their [FigShare archive](https://dx.doi.org/10.6084/m9.figshare.1525089) and renamed it to `54sp_s1_dRetal15pUhb.tree`. This input tree file is used to evaluate the fitness of calibration strategy 1 following [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066).

* [**Tree file following "BM23" calibrations**](trees/BM23pUhb/): the authors did not make their input files available as part of [their published point of view](https://doi.org/10.1093/sysbio/syad057), but they kindly shared them with us when we requested them. Following their analysis, the file we are using as input tree files is the following:
  * [**`54sp_s1_BM23pUhb.tree`**](trees/BM23pUhb/54sp_s1_BM23pUhb.tree): the node age constraints used in this input tree file are the same as in [`54sp_s1_dRetal15pUhb.tree`](trees/dRetal15pUhb/54sp_s1_dRetal15pUhb.tree) except for the maximum age of node Metazoa. [Budd and Mann (2023)](https://doi.org/10.1093/sysbio/syad057) used this input calibrated tree file in their inference analyses. We have renamed this file as [`54sp_s1_BM23pUhb.tree`](trees/BM23pUhb/54sp_s1_BM23pUhb.tree), and used it to compare the fitness of the different constraint for Metazoa against the one used by [dos Reis et al. (2015)](https://doi.org/10.1016/j.cub.2015.09.066).

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

