# fly_XQTL
## GitHub for the Drosophila X-QTL paper

### There are folders for notes, scripts, helperfiles:  
-notes contain text files to reproduce results in the paper  
-scripts are scripts that are required (typically bash/R/python)  
-helpfiles are various files used by the scripts   
### notes:  
* raw.data.txt  
  * code to massage the raw data, including fastq file to treatments, and SRA info
* make.bams.txt
  * code to make bams for each file, including code to sub-sample to lower coverage, or combine samples to make pseudo-high-coverage data
* call.SNPs.txt
  * code to call SNPs by sample.  We call at a fixed set of HQ SNPs based on founders. In the pools for each SNP we calculate the frequency of the ALT allele
* call.haplotypes.txt
 * call haplotypes separately by pooled sample using "limSolve" (see paper and Linder et al 2020).  Requires a pooled sample, and founder bams
* figure1.txt
  * code to generate representative haplotypes for Figure 1
* figure2.txt
  * code to generate Figure 2 = distribution of LOD scores at causative site
    * this is a simple model of equally frequent founders and Gaussian effects, only a marker at the QTL is simulated
    * we also simulate a QTL of zero effect to get a null distribution
* resolution.coarse.txt
  * code to look at resolution under coarse markers (markers every 200kb) including code for
    * Figure3, Supp Figure 2A, and Supp NRILs_looselinkage
* resolution.fine.txt
  * code to look at resolution under fine markers (markers every 10kb) for a subset of conditions (smaller region, n_gens_random_mate = 4), including code for
    * Supp Figure 2B
* manhattan.txt
  * code to reproduce the Manhattan plots (Fig 4 and Supp Fig 4), including code to make manhattan plots more generally (qqman wasn't flexible enough!)
* QQplot.txt
  * code to reproduce the QQ plot of Supp Figure 3
* haplotype.change.at.QTL.txt
  * code to reproduce Figure 5, showing changes in haplotypes at QTL (and local LOD scores)

### if analyzing real data
* raw.data.txt -> make.bams.txt -> call.SNPs.txt -> call.haplotypes.txt -> manhattan.txt -> haplotype.change.at.QTL.txt
