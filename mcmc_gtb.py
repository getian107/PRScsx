#!/usr/bin/env python

"""
Markov Chain Monte Carlo (MCMC) sampler for cross-ethnic polygenic prediction with continuous shrinkage (CS) priors - PRS-CSx.

"""


import scipy as sp
from scipy import linalg 
from numpy import random
import gigrnd


def mcmc(a, b, phi, snp_dict, beta_mrg, frq_dict, idx_dict, n, ld_blk, blk_size, n_iter, n_burnin, thin, pop, chrom, out_dir, out_name, meta, seed):
    print('... MCMC ...')

    # seed
    if seed != None:
        random.seed(seed)

    # derived stats
    n_pst = (n_iter-n_burnin)/thin
    n_pop = len(pop)
    p_tot = len(snp_dict['SNP'])

    p = {}
    n_blk = {}
    het = {}
    for pp in range(n_pop):
        p[pp] = len(beta_mrg[pp])
        n_blk[pp] = len(ld_blk[pp])
        het[pp] = sp.sqrt(2.0*frq_dict[pp]*(1.0-frq_dict[pp]))

    n_grp = sp.zeros((p_tot,1))
    for jj in range(p_tot):
        for pp in range(n_pop):
            if jj in idx_dict[pp]:
                n_grp[jj] += 1

    # initialization
    beta = {}
    sigma = {}
    for pp in range(n_pop):
        beta[pp] = sp.zeros((p[pp],1))
        sigma[pp] = 1.0

    psi = sp.ones((p_tot,1))

    if phi == None:
        phi = 1.0; phi_updt = True
    else:
        phi_updt = False

    # space allocation
    beta_est = {}
    beta_sq_est = {}
    sigma_est = {}
    for pp in range(n_pop):
        beta_est[pp] = sp.zeros((p[pp],1))
        beta_sq_est[pp] = sp.zeros((p[pp],1))
        sigma_est[pp] = 0.0

    psi_est = sp.zeros((p_tot,1))
    phi_est = 0.0

    # MCMC
    for itr in range(1,n_iter+1):
        if itr % 100 == 0:
            print('--- iter-' + str(itr) + ' ---')

        for pp in range(n_pop):
            mm = 0; quad = 0.0
            psi_pp = psi[idx_dict[pp]]
            for kk in range(n_blk[pp]):
                if blk_size[pp][kk] == 0:
                    continue
                else:
                    idx_blk = range(mm,mm+blk_size[pp][kk])
                    dinvt = ld_blk[pp][kk]+sp.diag(1.0/psi_pp[idx_blk].T[0])
                    dinvt_chol = linalg.cholesky(dinvt)
                    beta_tmp = linalg.solve_triangular(dinvt_chol, beta_mrg[pp][idx_blk], trans='T') \
                               + sp.sqrt(sigma[pp]/n[pp])*random.randn(len(idx_blk),1)
                    beta[pp][idx_blk] = linalg.solve_triangular(dinvt_chol, beta_tmp, trans='N')
                    quad += sp.dot(sp.dot(beta[pp][idx_blk].T, dinvt), beta[pp][idx_blk])
                    mm += blk_size[pp][kk]

            err = max(n[pp]/2.0*(1.0-2.0*sum(beta[pp]*beta_mrg[pp])+quad), n[pp]/2.0*sum(beta[pp]**2/psi_pp))
            sigma[pp] = 1.0/random.gamma((n[pp]+p[pp])/2.0, 1.0/err)

        delta = random.gamma(a+b, 1.0/(psi+phi))

        xx = sp.zeros((p_tot,1))
        for pp in range(n_pop):
            xx[idx_dict[pp]] += n[pp]*beta[pp]**2/sigma[pp]

        for jj in range(p_tot):
            while True:
                try:
                    psi[jj] = gigrnd.gigrnd(a-0.5*n_grp[jj], 2.0*delta[jj], xx[jj])
                except:
                    continue
                else:
                    break
        psi[psi>1] = 1.0

        if phi_updt == True:
            w = random.gamma(1.0, 1.0/(phi+1.0))
            phi = random.gamma(p_tot*b+0.5, 1.0/(sum(delta)+w))

        # posterior
        if (itr > n_burnin) and (itr % thin == 0):
            for pp in range(n_pop):
                beta_est[pp] = beta_est[pp] + beta[pp]/n_pst
                beta_sq_est[pp] = beta_sq_est[pp] + beta[pp]**2/n_pst
                sigma_est[pp] = sigma_est[pp] + sigma[pp]/n_pst

            psi_est = psi_est + psi/n_pst
            phi_est = phi_est + phi/n_pst

    # convert standardized beta to per-allele beta
    for pp in range(n_pop):
        beta_est[pp] /= het[pp]
        beta_sq_est[pp] /= het[pp]**2

    # meta
    if meta == 'TRUE':
        vv = sp.zeros((p_tot,1))
        zz = sp.zeros((p_tot,1))
        for pp in range(n_pop):
            vv[idx_dict[pp]] += 1.0/(beta_sq_est[pp]-beta_est[pp]**2)
            zz[idx_dict[pp]] += 1.0/(beta_sq_est[pp]-beta_est[pp]**2)*beta_est[pp]
        mu = zz/vv

    # write posterior effect sizes
    for pp in range(n_pop):
        if phi_updt == True:
            eff_file = out_dir + '/' + '%s_%s_pst_eff_a%d_b%.1f_phiauto_chr%d.txt' % (out_name, pop[pp], a, b, chrom)
        else:
            eff_file = out_dir + '/' + '%s_%s_pst_eff_a%d_b%.1f_phi%1.0e_chr%d.txt' % (out_name, pop[pp], a, b, phi, chrom)

        snp_pp = [snp_dict['SNP'][ii] for ii in idx_dict[pp]]
        bp_pp = [snp_dict['BP'][ii] for ii in idx_dict[pp]]
        a1_pp = [snp_dict['A1'][ii] for ii in idx_dict[pp]]
        a2_pp = [snp_dict['A2'][ii] for ii in idx_dict[pp]]

        with open(eff_file, 'w') as ff:
            for snp, bp, a1, a2, beta in zip(snp_pp, bp_pp, a1_pp, a2_pp, beta_est[pp]):
                ff.write('%d\t%s\t%d\t%s\t%s\t%.6e\n' % (chrom, snp, bp, a1, a2, beta))

    if meta == 'TRUE':
        if phi_updt == True:
            eff_file = out_dir + '/' + '%s_META_pst_eff_a%d_b%.1f_phiauto_chr%d.txt' % (out_name, a, b, chrom)
        else:
            eff_file = out_dir + '/' + '%s_META_pst_eff_a%d_b%.1f_phi%1.0e_chr%d.txt' % (out_name, a, b, phi, chrom)

        with open(eff_file, 'w') as ff:
            for snp, bp, a1, a2, beta in zip(snp_dict['SNP'], snp_dict['BP'], snp_dict['A1'], snp_dict['A2'], mu):
                ff.write('%d\t%s\t%d\t%s\t%s\t%.6e\n' % (chrom, snp, bp, a1, a2, beta))

    # print estimated phi
    if phi_updt == True:
        print('... Estimated global shrinkage parameter: %1.2e ...' % phi_est )

    print('... Done ...')


