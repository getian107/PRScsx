#!/usr/bin/env python

"""
PRS-CSx: a method that integrates GWAS summary statistics and external LD reference panels from multiple populations to improve cross-population polygenic prediction.
Posterior SNP effect sizes are inferred under coupled continuous shrinkage (CS) priors across populations.
PRS-CSx is an extension of the Bayesian polygenic prediction method PRS-CS.

References: T Ge, CY Chen, Y Ni, YCA Feng, JW Smoller. Polygenic Prediction via Bayesian Regression and Continuous Shrinkage Priors.
           Nature Communications, 10:1776, 2019.

           Ruan Y, Feng YCA, Chen CY, Lam M, Stanley Global Asia Initiatives, Sawa A, Martin AR, Qin S, Huang H, Ge T.
           Improving Polygenic Prediction in Ancestrally Diverse Populations.
           MedRxiv preprint, 2021.

Usage:
python PRScsx.py --ref_dir=PATH_TO_REFERENCE --bim_prefix=VALIDATION_BIM_PREFIX --sst_file=SUM_STATS_FILE --n_gwas=GWAS_SAMPLE_SIZE
                 --pop=POPULATION --out_dir=OUTPUT_DIR --out_name=OUTPUT_FILE_PREFIX
                 [--a=PARAM_A --b=PARAM_B --phi=PARAM_PHI --n_iter=MCMC_ITERATIONS --n_burnin=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR
                  --chrom=CHROM --meta=META_FLAG --seed=SEED]

 - PATH_TO_REFERENCE: Full path to the directory that contains the SNP information file and LD reference panels.
                      If the 1000 Genomes reference is used, the folder would contain the SNP information file snpinfo_mult_1kg_hm3
                      and one or more of the LD reference files: ldblk_1kg_afr, ldblk_1kg_amr, ldblk_1kg_eas, ldblk_1kg_eur, ldblk_1kg_sas;
                      if the UK Biobank reference is used, the folder would contain the SNP information file snpinfo_mult_ukbb_hm3 
                      and one or more of the LD reference files: ldblk_ukbb_afr, ldblk_ukbb_amr, ldblk_ukbb_eas, ldblk_ukbb_eur, ldblk_ukbb_sas.

 - VALIDATION_BIM_PREFIX: Full path and the prefix of the bim file for the target (validation/testing) dataset.
                          This file is used to provide a list of SNPs that are available in the target dataset.

 - SUM_STATS_FILE: Full path and the file name of the GWAS summary statistics.
                   Multiple GWAS summary statistics files are allowed and should be separated by comma.
                   Summary statistics files must have the following format (including the header line and order of the columns):
              
                   SNP          A1   A2   BETA      P
                   rs4970383    C    A    -0.0064   4.7780e-01
                   rs4475691    C    T    -0.0145   1.2450e-01
                   rs13302982   A    G    -0.0232   2.4290e-01
                   ...

                Or:

                   SNP          A1   A2   OR        P
                   rs4970383    A    C    0.9825    0.5737                 
                   rs4475691    T    C    0.9436    0.0691
                   rs13302982   A    G    1.1337    0.0209
                   ...

 - GWAS_SAMPLE_SIZE: Sample sizes of the GWAS, in the same order of the GWAS summary statistics files, separated by comma.

 - POPULATION: Population of the GWAS sample, in the same order of the GWAS summary statistics files, separated by comma.
               For both the 1000 Genomes reference and the UK Biobank reference, AFR, AMR, EAS, EUR and SAS are allowed.

 - OUTPUT_DIR: Output directory of the posterior effect size estimates.

 - OUTPUT_FILE_PREFIX: Output filename prefix of the posterior effect size estimates.

 - PARAM_A (optional): Parameter a in the gamma-gamma prior. Default is 1.

 - PARAM_B (optional): Parameter b in the gamma-gamma prior. Default is 0.5.

 - PARAM_PHI (optional): Global shrinkage parameter phi. If phi is not specified, it will be learnt from the data using a fully Bayesian approach.
                         This usually works well for polygenic traits with very large GWAS sample sizes (hundreds of thousands of subjects).
                         For GWAS with limited sample sizes (including most of the current disease GWAS),
                         fixing phi to 1e-2 (for highly polygenic traits) or 1e-4 (for less polygenic traits),
                         or doing a small-scale grid search (e.g., phi=1e-6, 1e-4, 1e-2, 1) to find the optimal phi value 
                         in the validation dataset often improves perdictive performance.

 - MCMC_ITERATIONS (optional): Total number of MCMC iterations. Default is the number of discovery populations * 1,000.

 - MCMC_BURNIN (optional): Number of burnin iterations. Default is the number of discovery populations * 500.

 - MCMC_THINNING_FACTOR (optional): Thinning factor of the Markov chain. Default is 5.

 - CHROM (optional): The chromosome on which the model is fitted, separated by comma, e.g., --chrom=1,3,5.
                     Parallel computation for the 22 autosomes is recommended. Default is iterating through 22 autosomes (can be time-consuming).

 - META_FLAG (optional): If True, return combined SNP effect sizes across populations using 
                         an inverse-variance-weighted meta-analysis of the population-specific posterior effect size estimates. Default is False.

 - SEED (optional): Non-negative integer which seeds the random number generator.

"""


import os
import sys
import getopt

import parse_genet
import mcmc_gtb
import gigrnd


def parse_param():
    long_opts_list = ['ref_dir=', 'bim_prefix=', 'sst_file=', 'a=', 'b=', 'phi=', 'n_gwas=', 'pop=',
                      'n_iter=', 'n_burnin=', 'thin=', 'out_dir=', 'out_name=', 'chrom=', 'meta=', 'seed=', 'help']

    param_dict = {'ref_dir': None, 'bim_prefix': None, 'sst_file': None, 'a': 1, 'b': 0.5, 'phi': None, 'n_gwas': None, 'pop': None,
                  'n_iter': None, 'n_burnin': None, 'thin': 5, 'out_dir': None, 'out_name': None, 'chrom': range(1,23), 'meta': 'FALSE', 'seed': None}

    print('\n')

    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_opts_list)          
        except:
            print('* Option not recognized.')
            print('* Use --help for usage information.\n')
            sys.exit(2)

        for opt, arg in opts:
            if opt == "-h" or opt == "--help":
                print(__doc__)
                sys.exit(0)
            elif opt == "--ref_dir": param_dict['ref_dir'] = arg
            elif opt == "--bim_prefix": param_dict['bim_prefix'] = arg
            elif opt == "--sst_file": param_dict['sst_file'] = arg.split(',')
            elif opt == "--a": param_dict['a'] = float(arg)
            elif opt == "--b": param_dict['b'] = float(arg)
            elif opt == "--phi": param_dict['phi'] = float(arg)
            elif opt == "--n_gwas": param_dict['n_gwas'] = list(map(int,arg.split(',')))
            elif opt == "--pop": param_dict['pop'] = arg.split(',')
            elif opt == "--n_iter": param_dict['n_iter'] = int(arg)
            elif opt == "--n_burnin": param_dict['n_burnin'] = int(arg)
            elif opt == "--thin": param_dict['thin'] = int(arg)
            elif opt == "--out_dir": param_dict['out_dir'] = arg
            elif opt == "--out_name": param_dict['out_name'] = arg
            elif opt == "--chrom": param_dict['chrom'] = arg.split(',')
            elif opt == "--meta": param_dict['meta'] = arg.upper()
            elif opt == "--seed": param_dict['seed'] = int(arg)
    else:
        print(__doc__)
        sys.exit(0)


    if param_dict['ref_dir'] == None:
        print('* Please specify the directory to the reference panel using --ref_dir\n')
        sys.exit(2)
    elif param_dict['bim_prefix'] == None:
        print('* Please specify the directory and prefix of the bim file for the target (validation/testing) dataset using --bim_prefix\n')
        sys.exit(2)
    elif param_dict['sst_file'] == None:
        print('* Please provide at least one summary statistics file using --sst_file\n')
        sys.exit(2)
    elif param_dict['n_gwas'] == None:
        print('* Please provide the sample size of the GWAS using --n_gwas\n')
        sys.exit(2)
    elif param_dict['pop'] == None:
        print('* Please specify the population of the GWAS sample using --pop\n')
        sys.exit(2)
    elif param_dict['out_dir'] == None:
        print('* Please specify the output directory using --out_dir\n')
        sys.exit(2)
    elif param_dict['out_name'] == None:
        print('* Please specify the prefix of the output file using --out_name\n')
        sys.exit(2)
    elif (len(param_dict['sst_file']) != len(param_dict['n_gwas']) or 
          len(param_dict['sst_file']) != len(param_dict['pop'])):
        print('* Length of sst_file, n_gwas and pop does not match\n')
        sys.exit(2)

    n_pop = len(param_dict['pop'])
    if param_dict['n_iter'] == None or param_dict['n_burnin'] == None:
        param_dict['n_iter'] = n_pop*1000
        param_dict['n_burnin'] = n_pop*500

    for key in param_dict:
        print('--%s=%s' % (key, param_dict[key]))

    print('\n')
    return param_dict


def main():
    param_dict = parse_param()
    n_pop = len(param_dict['pop'])

    for chrom in param_dict['chrom']:
        print('##### process chromosome %d #####' % int(chrom))

        if os.path.isfile(param_dict['ref_dir'] + '/snpinfo_mult_1kg_hm3'):
            ref = '1kg'
            ref_dict = parse_genet.parse_ref(param_dict['ref_dir'] + '/snpinfo_mult_1kg_hm3', int(chrom), ref)
        elif os.path.isfile(param_dict['ref_dir'] + '/snpinfo_mult_ukbb_hm3'):
            ref = 'ukbb'
            ref_dict = parse_genet.parse_ref(param_dict['ref_dir'] + '/snpinfo_mult_ukbb_hm3', int(chrom), ref)

        vld_dict = parse_genet.parse_bim(param_dict['bim_prefix'], int(chrom))

        sst_dict = {}
        for pp in range(n_pop):
            sst_dict[pp] = parse_genet.parse_sumstats(ref_dict, vld_dict, param_dict['sst_file'][pp], param_dict['pop'][pp], param_dict['n_gwas'][pp])

        ld_blk = {}
        blk_size = {}
        for pp in range(n_pop):
            ld_blk[pp], blk_size[pp] = parse_genet.parse_ldblk(param_dict['ref_dir'], sst_dict[pp], param_dict['pop'][pp], int(chrom), ref)

        snp_dict, beta_dict, frq_dict, idx_dict = parse_genet.align_ldblk(ref_dict, vld_dict, sst_dict, n_pop, int(chrom))

        mcmc_gtb.mcmc(param_dict['a'], param_dict['b'], param_dict['phi'], snp_dict, beta_dict, frq_dict, idx_dict, param_dict['n_gwas'], ld_blk, blk_size,
            param_dict['n_iter'], param_dict['n_burnin'], param_dict['thin'], param_dict['pop'], int(chrom),
            param_dict['out_dir'], param_dict['out_name'], param_dict['meta'], param_dict['seed'])

        print('\n')


if __name__ == '__main__':
    main()


