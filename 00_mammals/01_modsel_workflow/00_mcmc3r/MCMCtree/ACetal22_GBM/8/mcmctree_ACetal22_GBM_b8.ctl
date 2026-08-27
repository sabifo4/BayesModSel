          seed = -1
       seqfile = dummy.txt
      treefile = 72sp_ACetal22.tree
      mcmcfile = mcmc.txt
       outfile = out.txt

         ndata = 4
       seqtype = 0    * 0: nucleotides; 1:codons; 2:AAs
       usedata = 2    * 0: no data (prior); 1:exact likelihood;
                      * 2:approximate likelihood | 0: NT, 1: SQRT, 2: LOG, 3: ARCSIN (default: 3);
                      * 3:out.BV (in.BV)
         clock = 3    * 1: global clock; 2: independent rates; 3: correlated rates

         model = 4    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0.5  * alpha for gamma rates at sites
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0.1    * birth, death, sampling

   rgene_gamma = 2 40   * gammaDir prior for rate for genes
  sigma2_gamma = 1 10   * gammaDir prior for sigma^2     (for clock=2 or 3)

         print = 1   * 0: no mcmc sample; 1: everything except branch rates 2: everything
        burnin = 150000
      sampfreq = 400 
       nsample = 20000
BayesFactorBeta = 1.52858774526976e-08
