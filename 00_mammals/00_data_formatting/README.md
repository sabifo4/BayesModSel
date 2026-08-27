# Input files

## Input tree files

* **Tree file with node age constraints following the _Unconstrained_ model**: we are using the same tree file that we used in [Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1) to constrain the node ages of the mammal backbone tree with 72 taxa. We constrained the age of node Placentalia to be between 61.66 Ma and 162.5 Ma (soft bounds) (see [our `Fossil_calibrations_justification.pdf` file with all the details justifying all node age constraints](https://github.com/sabifo4/mammals_dating/blob/main/calibrations/Fossil_calibrations_justifications.pdf)). Regarding node Placentalia, our justifications when we published our article were the following:

  ```text
  The human-tenrec split is equivalent to the origin of the clade comprising Boreoeutheria (Laurasiatheria and Euarchontoglires) and Atlantogenata (Xenarthra and Afrotheria). Following Benton et al. 11 but with the geochronology revised following [23] and [6].

  [6] Hesselbo, S. P., Ogg, J. G., Ruhl, M., Hinnov, L. A. & Huang, C. J. in Geologic Time Scale 2020
  955-1021 (2020).
  [23]  Speijer, R. P., Pälike, H., Hollis, C. J., Hooker, J. J. & Ogg, J. G. in Geologic Time Scale 2020
  1087-1140 (2020).
  ```

  In this repository, we have saved this input tree file as [`72sp_ACetal22.tree`](trees/ACetal22/72sp_ACetal22.tree); the tag `72sp` refers to the number of taxa in the tree topology and `ACetal22` to the initials of the first author in the study followed by the year it was published ([Álvarez-Carretero et al. 2022](https://www.nature.com/articles/s41586-021-04341-1)).

> [!IMPORTANT]
>
> If you are keen on obtaining the input tree file we generated and used as part of our study ([Álvarez-Carretero et al. 2022](https://www.nature.com/articles/s41586-021-04341-1)), you can download it from the following sites:
>
>* [**Our Dropbox archive**](https://www.dropbox.com/s/53mdfyc47hukkrh/SeqBayesS1_MCMCtree_mainT2_posterior_newchrono.zip?dl=0) as specified in our [GitHub repository ("01_SeqBayes_S1/02_MCMCtree/README.md")](https://github.com/sabifo4/mammals_dating/tree/main/01_SeqBayes_S1/02_MCMCtree). Note that the calibrations used to constrain the age of various nodes in this tree are those that we used when we repeated our analyses following the updated geochronology as of September 2021. You can download any of tree files inside `02_MCMCtree_posterior_newchrono/run[1-4]/mcmctree_GBM/tree.txt` -- same tree file used for the 4 independent chains we ran.
>* [**Our FigShare archive**](https://figshare.com/articles/dataset/Data_for_A_Species-Level_Timeline_of_Mammal_Evolution_Integrating_Phylogenomic_Data_/14885691), as specified both in our article ([Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1)) and [our GitHub repository (main page)](https://github.com/sabifo4/mammals_dating/tree/main). You can access the `00_main_tree_T2_72sp-updated-geochronology.tree` file in directory `trees/00_step1` (i.e., this is how we renamed our calibrated tree file when we archived our results after publishing our analyses).

* **Tree file with node age constraints following the _Post-K-Pg_ model**: while the authors did not make their input files available as part of [their supplementary material](https://doi.org/10.1093/sysbio/syad057), they kindly shared them with us after requesting them. As stated in their article, their **_Post-K-Pg_** model included the following constraints:

  ```text
  Quote: "we re-ran the analysis using a uniform calibration for the placentals and other deep nodes within them from the age of their oldest fossil down to the K-Pg boundary, with soft bounds on both ages. We note that similar calibration ranges are used for 11 other nodes in their analysis (e.g., Caniformia, Chiroptera, Carnivora, etc)".
  ```

  We have renamed their input tree file as [`72sp_BM23.tree`](trees/BM23/72sp_BM23.tree); the tag `72sp` refers to the number of taxa in the tree topology and `BM23` to the initials of the authors in the study followed by the year it was published ([Budd and Mann 2023](https://doi.org/10.1093/sysbio/syad057)). Please note that all node age constraints are the same as in [Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1) (i.e., [`72sp_ACetal22.tree`](trees/ACetal22/72sp_ACetal22.tree)) except for the maximum age of node Placentalia and that of the nodes within this clade, which are fixed to the K-Pg boundary: 66.09Ma. The nodes which constraints change in comparison to `ACetal22` after following the **_Post K-Pg_** diversification model are the following: Primates, Anthropoidea, Strepsirhini, Rodentia, Lagomorpha, Euungulata, Artiodactyla, "Whip-Rum" (i.e., Cetruminantia), Carnivora, Caniformia, Chiroptera, Euarchontoglires, Lipotyphla, Afrotheria, Paenungulata, and Xenarthra.

* **Tree file with node age constraints following the _Post-K-Pg_ model**: we updated the input tree file [`72sp_BM23.tree`](trees/BM23/72sp_BM23.tree) so that the calibration on node Placentalia has a hard bound on its maximum age: instead of using a tail percentage of 2.5%, we used "1e-300" as the closest value to 0 that does not cause numerical issues but still is used as a strict hard bound. To this end, all calibration densities are the same except for the one constraining the age of node Placentalia. We have saved the updated tree file as [`72sp_BM23pUhb.tree`](trees/BM23pUhb/72sp_BM23pUhb.tree); the suffix `pUhb` is used to inform that the tail **p**ercentage of the **U**pper bound in the calibration density to constrain the age of node Placentalia has a **h**ard **b**ound.

## Input alignment files

* **[`in.BV`](aln/in.BV)**: the branch lengths, the gradient, and the Hessian we inferred with `BASEML` ([Yang 2007](https://academic.oup.com/mbe/article/24/8/1586/1103731)) as part of the analyses in our study ([Álvarez-Carretero et al. 2022)](https://www.nature.com/articles/s41586-021-04341-1)) are saved in the [`in.BV` file](aln/in.BV). The vector of branch lengths, the gradient, and the Hessian matrix are subsequently used by `MCMCtree` to approximate the likelihood calculation following the approach described by [dos Reis and Yang (2011)](https://academic.oup.com/mbe/article/28/7/2161/1051613). Please note that the [`in.BV` file](aln/in.BV) we obtained as part of our study ([Álvarez-Carretero et al. 2022)](https://www.nature.com/articles/s41586-021-04341-1)) is the same one you can find in this repository. You can also retrieve this file from the following sites:
  * [**Our Dropbox archive**](https://www.dropbox.com/s/w6xnoleo4ssjalv/SeqBayesS1_BASEML_method1_T2hyp.zip?dl=0): as specified in our [GitHub repository ("01_SeqBayes_S1/01_BASEML_02_Hessian/README.md")](https://github.com/sabifo4/mammals_dating/tree/main/01_SeqBayes_S1/01_BASEML/02_Hessian), you can find the `in.BV.02p1-p4` file (i.e., this is how we named our `in.BV` file while carrying out our analyses) in directory `000_Hessian/02_atlantogenata_tarver2016`.
  * [**Our FigShare archive**](https://figshare.com/articles/dataset/Data_for_A_Species-Level_Timeline_of_Mammal_Evolution_Integrating_Phylogenomic_Data_/14885691): as specified both in our study ([Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1)) and [our GitHub repository (main page)](https://github.com/sabifo4/mammals_dating/tree/main), you can access the `00_main_tree_T2_in.BV` file (i.e., this is how we renamed our `in.BV` file when we archived our results after publishing our analyses so that it was clearer for users to know the tree topology that had been fixed to calculate the MLEs of branch lengths, the gradient, and the Hessian) in directory `inBV/00_step1/`.
* [**Dummy alignment file, `dummy.txt`**](aln/dummy.txt): this alignment file can be used when the `in.BV` file has already been generated (i.e., when `MCMCtree` is used with option `usedata = 2`; everything this program needs to approximate the likelihood calculation is saved in the `in.BV` file) as `MCMCtree` will read a smaller file to carry out one of the first checks: `MCMCtree` first reads the input files (sequence alignment file and tree file) and checks that the tags used for taxa present in the alignment and the phylogeny are the same. To this end, `MCMCtree` will need to allocate less memory when reading a smaller "dummy" alignment instead of the large partitioned alignment file, and will also be faster! In this study, we used the same `dummy.txt` file that we used in ([Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1)): you can retrieve the original file from [our Dropbox archive](https://www.dropbox.com/s/53mdfyc47hukkrh/SeqBayesS1_MCMCtree_mainT2_posterior_newchrono.zip?dl=0) inside directory `02_MCMCtree_posterior_newchrono/run[1-4]/mcmctree_GBM`.

> [!IMPORTANT]
>
> If you are keen on obtaining the partitioned alignment file we generated and used in [Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1), you can retrieve it from the following sites:
>
>* [**Our Dropbox archive**](https://www.dropbox.com/s/mrvzzvd4o6qqyqk/000_alignments.zip?dl=0): as detailed in [our GitHub repository ("01_SeqBayes_S1/00_Gene_filtering/README.md")](https://github.com/sabifo4/mammals_dating/tree/main/01_SeqBayes_S1/00_Gene_filtering), you can download the compressed file `000_alignments.zip` and, once decompressed, you can locate the partitioned alignment file `alignment_4part.aln`. The file structure you shall see is the following:
>
> ```text
> 000_alignments/
>   |- part1/                    # <-- genes included in the first partition + analyses using the `fasta-phylip-partition` pipeline
>   |- part2/                    # <-- genes included in the first partition + analyses using the `fasta-phylip-partition` pipeline
>   |- part3/                    # <-- genes included in the first partition + analyses using the `fasta-phylip-partition` pipeline
>   |- part4/                    # <-- genes included in the first partition + analyses using the `fasta-phylip-partition` pipeline
>   |- alignment_4parts.aln      # <-- this is the partitioned alignment you are looking for!
>   |- count_missingdat_72sp.pl  # <-- in-house bash script to count missing data in our alignment blocks
>   |- README.md                 # <-- further information about our analyses with this dataset -- they are also detailed in our GitHub repository
> ```
>
>* [**Our FigShare archive**](https://figshare.com/articles/dataset/Data_for_A_Species-Level_Timeline_of_Mammal_Evolution_Integrating_Phylogenomic_Data_/14885691): as specified both in our study ([Álvarez-Carretero et al. (2022)](https://www.nature.com/articles/s41586-021-04341-1)) and [our GitHub repository (main page)](https://github.com/sabifo4/mammals_dating/tree/main), you can access the `alignment_4parts.aln` file in directory `aln/00_step1`.
>
> **[NOTE]**: please choose one repository or another to retrieve the partitioned alignment. We are unable to host this file on GitHub as it is larger than 100Mb; the limit allowed to store files on this site. Note that we do not need to use the partitioned alignment for our Bayesian model selection analyses as we already have the `in.BV` file with the branch lengths, the gradient, and the Hessian required to approximate the likelihood calculation. The only thing that `MCMCtree` needs to find is a sequence alignment with sequence names that match those in the input tree file: this condition needs to be meet so that the program can run. We can reduce the memory required to read the alignment file if we use a dummy alignment file instead of the large phylogenomic alignment with four partitions, thus being the main reason for using the `dummy.txt` file.

----

Once you have gotten familiar with how our input data have been formatted and renamed, you are ready to continue with Bayesian model selection. You can continue your journey by following the steps described in [this `README.md` file](../01_modsel_workflow/README.md)!
