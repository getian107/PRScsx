# PRS-CSx

**PRS-CSx** is a Python based command line tool that integrates GWAS summary statistics and external LD reference panels from multiple populations to improve cross-population polygenic prediction. Posterior SNP effect sizes are inferred under coupled continuous shrinkage (CS) priors across populations. PRS-CSx is an extension of the Bayesian polygenic prediction method PRS-CS (https://github.com/getian107/PRScs), described in the article:

T Ge, CY Chen, Y Ni, YCA Feng, JW Smoller. Polygenic Prediction via Bayesian Regression and Continuous Shrinkage Priors. *Nature Communications*, 10:1776, 2019.

Preprint of PRS-CSx is available on medrxiv: https://www.medrxiv.org/content/10.1101/2020.12.27.20248738v2


## Recent Version History

**July 29, 2021**: Changed default MCMC parameters.

**Jun 4, 2021**: Expanded reference panels to five populations.

**May 26, 2021**: Added suggestions for limiting the number of threads in scipy when running PRS-CS (see Computational Efficiency section below).

**Apr 6, 2021**: Added projection of the LD matrix to its nearest non-negative definite matrix.

**Mar 4, 2021**: LD reference panels constructed using the UK Biobank data are now available. 

**Jan 4, 2021**: Improved the accuracy and robustness of random sampling from the generalized inverse Gaussian distribution. Prediction accuracy will probably slightly improve over the initial release.

**Dec 20, 2020**: Repository made public.


## Getting Started

- Clone this repository using the following git command:

    `git clone https://github.com/getian107/PRScsx.git`

    Alternatively, download the source files from the github website (`https://github.com/getian107/PRScsx`)
    
- Download the LD reference panels and extract files:

    LD reference panels constructed using the 1000 Genomes Project phase 3 samples:
    
     [AFR reference](https://www.dropbox.com/s/mq94h1q9uuhun1h/ldblk_1kg_afr.tar.gz?dl=0 "AFR reference") (~4.44G);
     `tar -zxvf ldblk_1kg_afr.tar.gz`
     
     [AMR reference](https://www.dropbox.com/s/uv5ydr4uv528lca/ldblk_1kg_amr.tar.gz?dl=0 "AMR reference") (~3.84G);
     `tar -zxvf ldblk_1kg_amr.tar.gz`
        
     [EAS reference](https://www.dropbox.com/s/7ek4lwwf2b7f749/ldblk_1kg_eas.tar.gz?dl=0 "EAS reference") (~4.33G);
     `tar -zxvf ldblk_1kg_eas.tar.gz`
        
     [EUR reference](https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=0 "EUR reference") (~4.56G);
     `tar -zxvf ldblk_1kg_eur.tar.gz`
     
     [SAS reference](https://www.dropbox.com/s/hsm0qwgyixswdcv/ldblk_1kg_sas.tar.gz?dl=0 "SAS reference") (~5.60G);
     `tar -zxvf ldblk_1kg_sas.tar.gz`
    
    LD reference panels constructed using the UK Biobank data ([Notes](https://www.dropbox.com/s/y3hsc15kwjxwjtd/UKBB_ref.txt?dl=0 "Notes")):
    
     [AFR reference](https://www.dropbox.com/s/dtccsidwlb6pbtv/ldblk_ukbb_afr.tar.gz?dl=0 "AFR reference") (~4.93G);
     `tar -zxvf ldblk_ukbb_afr.tar.gz`
     
     [AMR reference](https://www.dropbox.com/s/y7ruj364buprkl6/ldblk_ukbb_amr.tar.gz?dl=0 "AMR reference") (~4.10G);
     `tar -zxvf ldblk_ukbb_amr.tar.gz`
    
     [EAS reference](https://www.dropbox.com/s/fz0y3tb9kayw8oq/ldblk_ukbb_eas.tar.gz?dl=0 "EAS reference") (~5.80G);
     `tar -zxvf ldblk_ukbb_eas.tar.gz`
    
     [EUR reference](https://www.dropbox.com/s/t9opx2ty6ucrpib/ldblk_ukbb_eur.tar.gz?dl=0 "EUR reference") (~6.25G);
     `tar -zxvf ldblk_ukbb_eur.tar.gz`
    
     [SAS reference](https://www.dropbox.com/s/nto6gdajq8qfhh0/ldblk_ukbb_sas.tar.gz?dl=0 "SAS reference") (~7.37G);
     `tar -zxvf ldblk_ukbb_sas.tar.gz`
    
    Note that these files are identical to the reference panels used in **PRS-CS**. Therefore, there is no need to download again if you are already using **PRS-CS**.

- Download the SNP information file and put it in the same folder containing the reference panels:

    1000 Genomes reference: [SNP info](https://www.dropbox.com/s/rhi806sstvppzzz/snpinfo_mult_1kg_hm3?dl=0 "SNP info") (~106M)
    
    UK Biobank reference: [SNP info](https://www.dropbox.com/s/oyn5trwtuei27qj/snpinfo_mult_ukbb_hm3?dl=0 "SNP info") (~108M)
    
- PRScsx requires Python packages **scipy** (https://www.scipy.org/) and **h5py** (https://www.h5py.org/) installed.
 
- Once Python and its dependencies have been installed, running

    `./PRScsx.py --help` or `./PRScsx.py -h`

    will print a list of command-line options.
    

## Using PRS-CSx

`
python PRScsx.py --ref_dir=PATH_TO_REFERENCE --bim_prefix=VALIDATION_BIM_PREFIX --sst_file=SUM_STATS_FILE --n_gwas=GWAS_SAMPLE_SIZE --pop=POPULATION --out_dir=OUTPUT_DIR --out_name=OUTPUT_FILE_PREFIX [--a=PARAM_A --b=PARAM_B --phi=PARAM_PHI --n_iter=MCMC_ITERATIONS --n_burnin=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR --chrom=CHROM --meta=META_FLAG --seed=SEED]
`
 - PATH_TO_REFERENCE (required): Full path to the directory that contains the SNP information file and LD reference panels. If the 1000 Genomes reference is used, the folder would contain the SNP information file `snpinfo_mult_1kg_hm3` and one or more of the LD reference files: `ldblk_1kg_afr`, `ldblk_1kg_amr`, `ldblk_1kg_eas`, `ldblk_1kg_eur`, `ldblk_1kg_sas`; if the UK Biobank reference is used, the folder would contain the SNP information file `snpinfo_mult_ukbb_hm3` and one or more of the LD reference files: `ldblk_ukbb_afr`, `ldblk_ukbb_amr`, `ldblk_ukbb_eas`, `ldblk_ukbb_eur`, `ldblk_ukbb_sas`.

 - VALIDATION_BIM_PREFIX (required): Full path and the prefix of the bim file for the target (validation/testing) dataset. This file is used to provide a list of SNPs that are available in the target dataset.

 - SUM_STATS_FILE (required): Full path and the file name of the GWAS summary statistics. Multiple GWAS summary statistics files are allowed and should be separated by comma. Summary statistics files must have the following format (including the header line and order of the columns):


```
    SNP          A1   A2   BETA      P
    rs4970383    C    A    -0.0064   4.7780e-01
    rs4475691    C    T    -0.0145   1.2450e-01
    rs13302982   A    G    -0.0232   2.4290e-01
    ...
```
Or:
```
    SNP          A1   A2   OR        P
    rs4970383    A    C    0.9825    0.5737                 
    rs4475691    T    C    0.9436    0.0691
    rs13302982   A    G    1.1337    0.0209
    ...
```
where SNP is the rs ID, A1 is the effect allele, A2 is the alternative allele, BETA/OR is the effect/odds ratio of the A1 allele, P is the p-value of the effect. In fact, BETA/OR is only used to determine the direction of an association. Therefore if z-scores or even +1/-1 indicating effect directions are presented in the BETA column, the algorithm should still work properly.

 - GWAS_SAMPLE_SIZE (required): Sample sizes of the GWAS, in the same order of the GWAS summary statistics files, separated by comma.

 - POPULATION (required): Population of the GWAS sample, in the same order of the GWAS summary statistics files, separated by comma. For both the 1000 Genomes reference and the UK Biobank reference, AFR, AMR, EAS, EUR and SAS are allowed.

 - OUTPUT_DIR (required): Output directory of the posterior effect size estimates.

 - OUTPUT_FILE_PREFIX (required): Output filename prefix of the posterior effect size estimates.

 - PARAM_A (optional): Parameter a in the gamma-gamma prior. Default is 1.

 - PARAM_B (optional): Parameter b in the gamma-gamma prior. Default is 0.5.

 - PARAM_PHI (optional): Global shrinkage parameter phi. If phi is not specified, it will be learnt from the data using a fully Bayesian approach. This usually works well for polygenic traits with very large GWAS sample sizes (hundreds of thousands of subjects). For GWAS with limited sample sizes (including most of the current disease GWAS), fixing phi to 1e-2 (for highly polygenic traits) or 1e-4 (for less polygenic traits), or doing a small-scale grid search (e.g., phi=1e-6, 1e-4, 1e-2, 1) to find the optimal phi value in the validation dataset often improves perdictive performance.

 - MCMC_ITERATIONS (optional): Total number of MCMC iterations. Default is 1,000 * the number of discovery populations.

 - MCMC_BURNIN (optional): Number of burnin iterations. Default is 500 * the number of discovery populations. Both `--n_iter` and `--n_burnin` need to be specified to overwrite their default values.

 - MCMC_THINNING_FACTOR (optional): Thinning factor of the Markov chain. Default is 5.

 - CHROM (optional): The chromosome on which the model is fitted, separated by comma, e.g., --chrom=1,3,5. Parallel computation for the 22 autosomes is recommended. Default is iterating through 22 autosomes (can be time-consuming).

 - META_FLAG (optional): If True, return combined SNP effect sizes across populations using an inverse-variance-weighted meta-analysis of the population-specific posterior effect size estimates. Default is False.

 - SEED (optional): Non-negative integer which seeds the random number generator.


## Output

For each input GWAS, PRS-CSx writes posterior SNP effect size estimates for each chromosome to the user-specified directory. The output file contains chromosome, rs ID, base position, A1, A2 and posterior effect size estimate for each SNP. If `--meta=True`, meta-analyzed posterior effect sizes will also be written to the output directory. An individual-level polygenic score can be produced by concatenating output files from all chromosomes and then using `PLINK`'s `--score` command (https://www.cog-genomics.org/plink/1.9/score). If polygenic scores are generated by chromosome, use the 'sum' modifier so that they can be combined into a genome-wide score.

We recommend calculating one polygenic score for each discovery population using population-specific posterior SNP effect size estimates and learn a linear combination of the polygenic scores that most accurately predicts the trait in the validation dataset. The predictive performance of the method can be assessed in an independent dataset, using the optimal global shrinkage parameter phi (if grid search is used) and weights for the linear combination learnt in the validation dataset. We reccommend standardizing the polygenic scores (i.e., converting the scores to zero mean and unit variance) in both validation and testing datasets before learning/applying the linear combination.


## Computational Efficiency

PRS-CSx relies on scipy packages, which automatically use all available cores on a compute node. This can be problematic when running PRS-CSx on a compute cluster; PRS-CSx jobs may interfere with other jobs running on the same node, reducing computational efficiency. To resolve this issue, including the following code in the script to specify the number of threads in scipy:

```
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS
```
For example, to use a single thread for the computation, set `N_THREADS=1`.


## Test Data

The test data contains EUR and EAS GWAS summary statistics and a bim file for 1,000 SNPs on chromosome 22.
An example to use the test data:

`
python PRScsx.py --ref_dir=path_to_ref --bim_prefix=path_to_bim/test --sst_file=path_to_sumstats/EUR_sumstats.txt,path_to_sumstats/EAS_sumstats.txt --n_gwas=200000,100000 --pop=EUR,EAS --chrom=22 --phi=1e-2 --out_dir=path_to_output --out_name=test
`

The test data analysis would be finished in approximately 1 min when using 8Gb of RAM. 


## Data release

Posterior SNP effect size estimates for traits examined in the [PRS-CSx publication](https://www.nature.com/articles/s41588-022-01054-7 "PRS-CSx publication") are available [here](https://www.dropbox.com/sh/5v1bzlukxoor9fi/AACA580wl_gNKapqWvx3siOza?dl=0 "here").

Posterior SNP effect size estimates for the trans-ancestry T2D PRS developed and validated in the Genome Medicine publication are available [here](https://www.dropbox.com/s/c3r880fdthj0vhh/psteff_t2d_mahajan_media_bbj?dl=0 "here").


## Support

Please direct any problems or questions to Tian Ge (tge1@mgh.harvard.edu).


