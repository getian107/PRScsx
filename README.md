# PRS-CSx

**PRS-CSx** is a Python based command line tool that integrates GWAS summary statistics and external LD reference panels from multiple populations to improve cross-population polygenic prediction. Posterior SNP effect sizes are inferred under coupled continuous shrinkage (CS) priors across populations. PRS-CSx is an extension of the Bayesian polygenic prediction method PRS-CS (https://github.com/getian107/PRScs), described in the article:

T Ge, CY Chen, Y Ni, YCA Feng, JW Smoller. Polygenic Prediction via Bayesian Regression and Continuous Shrinkage Priors. *Nature Communications*, 10:1776, 2019.


## Getting Started

- Clone this repository using the following git command:

    `git clone https://github.com/getian107/PRScsx.git`

    Alternatively, download the source files from the github website (`https://github.com/getian107/PRScsx`)
    
- Download the LD reference computed using the 1000 Genomes phase 3 samples, and extract files:

    [EUR reference](https://www.dropbox.com/s/p9aqanhxvxaqv8k/ldblk_1kg_eur.tar.gz?dl=0 "EUR reference") (~4.56G);
    `tar -zxvf ldblk_1kg_eur.tar.gz`

    [EAS reference](https://www.dropbox.com/s/o2yo2x7icu1xtpn/ldblk_1kg_eas.tar.gz?dl=0 "EAS reference") (~4.33G);
    `tar -zxvf ldblk_1kg_eas.tar.gz`
    
    [AFR reference](https://www.dropbox.com/s/mq94h1q9uuhun1h/ldblk_1kg_afr.tar.gz?dl=0 "AFR reference") (~4.44G);
    `tar -zxvf ldblk_1kg_afr.tar.gz`

    Note that these files are identical to the reference panels used in **PRS-CS**, and therefore there is no need to download again if you are already using **PRS-CS**.

- Download the SNP information file and put it in the same folder containing the reference panels:
    [SNP info](https://www.dropbox.com/s/rhi806sstvppzzz/snpinfo_mult_hm3?dl=0 "SNP info") (~77M)
    
- PRScsx requires Python packages **scipy** (https://www.scipy.org/) and **h5py** (https://www.h5py.org/) installed.
 
- Once Python and its dependencies have been installed, running

    `./PRScsx.py --help` or `./PRScsx.py -h`

    will print a list of command-line options.
    

## Using PRS-CSx

`
python PRScsx.py --ref_dir=PATH_TO_REFERENCE --bim_prefix=VALIDATION_BIM_PREFIX --sst_file=SUM_STATS_FILE --n_gwas=GWAS_SAMPLE_SIZE --pop=POPULATION --out_dir=OUTPUT_DIR --out_name=OUTPUT_FILE_PREFIX [--a=PARAM_A --b=PARAM_B --phi=PARAM_PHI --n_iter=MCMC_ITERATIONS --n_burnin=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR --chrom=CHROM --meta=META_FLAG --seed=SEED]
`
- PATH_TO_REFERENCE (required): Full path to the directory that contains the SNP information file `snpinfo_mult_hm3` and `ldblk_1kg_eur` and/or `ldblk_1kg_eas` and/or `ldblk_1kg_afr`.

 - VALIDATION_BIM_PREFIX (required): Full path and the prefix of the bim file for the validation/testing set.

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
where SNP is the rs ID, A1 is the reference/effect allele, A2 is the alternative allele, BETA/OR is the effect/odds ratio of the reference allele, P is the p-value of the effect. In fact, BETA/OR is only used to determine the direction of an association, and therefore if z-scores or even +1/-1 indicating effect directions are presented in the BETA column, the algorithm should still work properly.

 - GWAS_SAMPLE_SIZE (required): Sample sizes of the GWAS, in the same order of the GWAS summary statistics files, separated by comma.

 - POPULATION (required): Population of the GWAS sample (EUR, EAS or AFR), in the same order of the GWAS summary statistics files, separated by comma.

 - OUTPUT_DIR (required): Output directory of the posterior effect size estimates.

 - OUTPUT_FILE_PREFIX (required): Output filename prefix of the posterior effect size estimates.

 - PARAM_A (optional): Parameter a in the gamma-gamma prior. Default is 1.

 - PARAM_B (optional): Parameter b in the gamma-gamma prior. Default is 0.5.

 - PARAM_PHI (optional): Global shrinkage parameter phi. If phi is not specified, it will be learnt from the data using a fully Bayesian approach. This usually works well for polygenic traits with very large GWAS sample sizes (hundreds of thousands of subjects). For GWAS with limited sample sizes (including most of the current disease GWAS), fixing phi to 1e-2 (for highly polygenic traits) or 1e-4 (for less polygenic traits), or doing a small-scale grid search (e.g., phi=1e-6, 1e-4, 1e-2, 1) to find the optimal phi value often improves perdictive performance.

 - MCMC_ITERATIONS (optional): Total number of MCMC iterations. Default is 1,000.

 - MCMC_BURNIN (optional): Number of burnin iterations. Default is 500.

 - MCMC_THINNING_FACTOR (optional): Thinning of the Markov chain. Default is 5.

 - CHROM (optional): The chromosome on which the model is fitted, separated by comma, e.g., --chrom=1,3,5. Parallel computation for the 22 autosomes is recommended. Default is iterating through 22 autosomes (can be time-consuming).

 - META_FLAG (optional): If True, return combined effect size estimates across populations using an inverse-variance-weighted meta-analysis, in addition to the population-specific posterior effect size estimates. Default is False.

 - SEED (optional): Non-negative integer which seeds the random number generator.


## Output

For each input GWAS, PRS-CSx writes posterior SNP effect size estimates for each chromosome to the user-specified directory. The output file contains chromosome, rs ID, base position, A1, A2 and posterior effect size estimate for each SNP. If `--meta=True`, meta-analyzed posterior effect sizes will also be written to the output directory. An individual-level polygenic score can be produced by concatenating output files from all chromosomes and then using `PLINK`'s `--score` command (https://www.cog-genomics.org/plink/1.9/score). If polygenic scores are generated by chromosome, use the 'sum' modifier so that they can be combined into a genome-wide score.


## Test Data

The test data contains EUR and EAS GWAS summary statistics and a bim file for 1,000 SNPs on chromosome 22.
An example to use the test data:

`
python PRScsx.py --ref_dir=path_to_ref --bim_prefix=path_to_bim/test --sst_file=path_to_sumstats/EUR_sumstats.txt,path_to_sumstats/EAS_sumstats.txt --n_gwas=200000,100000 --pop=EUR,EAS --chrom=22 --phi=1e-2 --out_dir=path_to_output --out_name=test
`


## Support

Please direct any problems or questions to Tian Ge (tge1@mgh.harvard.edu).


