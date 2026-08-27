          seed = -1
       seqfile = dummy.txt
      treefile = 54sp_s1_dRetal15pUhb.tree
      mcmcfile = mcmc.txt
       outfile = out.txt

         ndata = 2  * 2 partitions
       seqtype = 2  * 0: nucleotides; 1:codons; 2:AAs
       usedata = 2  * 0: no data (prior); 1:exact likelihood; 0: NT, 1: SQRT, 2: LOG, 3: ARCSIN
                    * 2:approximate likelihood; 3:out.BV (in.BV)
                    * NOTE: option `aaRatefile` not needed as the `in.BV` is being used now!
         clock = 3   * 1: global clock; 2: independent rates; 3: correlated rates

         model = 3        * models for AAs or codon-translated AAs:
                          *     0:poisson, 1:proportional,2:Empirical,3:Empirical+F
                          *     6:FromCodon, 8:REVaa_0, 9:REVaa(nr=189)
         alpha = 0.5      * alpha for gamma rates at sites
         ncatG = 4        * No. categories in discrete gamma

     cleandata = 0        * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0    * birth, death, sampling

   rgene_gamma = 2 40     * gammaDir prior for rate for genes
  sigma2_gamma = 1 10     * gammaDir prior for sigma^2     (for clock=2 or 3)

         print = 1        * 0: no mcmc sample; 1: everything except branch rates 2: everything
        burnin = 100000
      sampfreq = 1000 
       nsample = 20000
BayesFactorBeta = 0.0042785088007804
