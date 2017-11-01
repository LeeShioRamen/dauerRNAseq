"""
Compute Bayesian probabilities of Dauer formation
using Markov chain Monte Carlo (MCMC)
"""
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import pymc3 as pm

sns.set(style='white', palette='Set2', font_scale=2)

# get the data
data_dir = './data/sbt1dauercounts.txt'
dauer = pd.read_csv(data_dir)
dauer = dauer[['strain','dauer','non-dauer']]
# pool data by strain
pooled = dauer.groupby('strain', as_index=False).sum()
# get n number of trials
pooled['n'] = pooled['dauer'] + pooled['non-dauer']
# rename columns, n->trial and r->dauers
df = pooled[['strain','n', 'dauer']].rename(columns={'dauer':'r'})
# rename wild type for convenience
df = df.replace('wild type', 'WT')

# dataframe to store MCMC traces
df_mcmc = pd.DataFrame()
for strain in df.strain:
    # construct pymc3 bayesian model
    with pm.Model() as model:
        # uniform prior for p, bounds 0 and 1 since it is a probability
        p = pm.Uniform('p', lower=0, upper=1)
        # number of trials and dauer formation events
        n, r = df.loc[df['strain']==strain, ['n', 'r']].values.flatten()
        # binomial likelihood
        y = pm.Binomial('y', n=n, p=p, observed=r)
        # get the samples
        trace = pm.sample(draws=2000, tune=2000, init='advi+adapt_diag', njobs=4)
    # save to dataframe
    df_mcmc[strain] = trace.get_values('p')

# Make a summary dataframe
# Because this is pretty gaussian, median is pretty much the same
# as the mean, although we might as well compute it since we have the distribution
# one standard deviation covers 68% of the distribution, might be better to
# keep 95% HPD

# get strain names first
strains = [*df_mcmc.columns]
# dataframe to store summary statistics
df_summary = pd.DataFrame(index=['median', '_hpd', 'hpd_', 'mean','std'],
        columns=strains)
for strain in strains:
    # median
    df_summary.loc['median', strain] = np.median(df_mcmc[strain])
    # 95% highest posterior density
    df_summary.loc[['_hpd', 'hpd_'], strain] = pm.hpd(df_mcmc[strain], alpha=0.05)
    # mean
    df_summary.loc['mean', strain] = np.mean(df_mcmc[strain])
    # standard deviation
    df_summary.loc['std', strain] = np.std(df_mcmc[strain])

# plot sample histograms
for strain in strains:
     plt.hist(df_mcmc[strain], bins=100, normed=True, histtype='step', linewidth=2)
plt.xlabel('prob. of dauer, $p$')
plt.ylabel(r'$P(p\mid d, n)$')
plt.legend(strains, loc='upper center');
sns.despine()
plt.tight_layout()
#plt.savefig('./output/probdistrib_dauer.pdf', transparent=True, bbox_inches='tight')
plt.close('all')


# Subtract samples to get difference in dauer prob. and make summary dataframe
df_diff = pd.DataFrame()
df_diff_summary = pd.DataFrame()
strains_woWT = [s for s in strains if 'WT' not in s]
for spair in [('WT', s) for s in strains_woWT]:
    p_name = '{0} - {1}'.format(*spair)
    df_diff[p_name] = df_mcmc[spair[0]] - df_mcmc[spair[1]]
    df_diff_summary[p_name] = df_summary[spair[0]] - df_summary[spair[1]]
    # correct stdev
    df_diff_summary.loc['std', p_name] = df_diff[p_name].std()


# Plot histograms of differences with mode +- 95% HPD
for i, col in enumerate(df_diff):
    _ = plt.hist(df_diff[col], bins=100, normed=True, 
                 histtype='step', linewidth=2)

# Label axes
plt.xlabel('prob. of dauer, $p$')
plt.ylabel(r'$P(p\mid d, n)$')
plt.legend(df_diff.columns, loc=0);
sns.despine()
plt.tight_layout()
#plt.savefig('./output/probdiff_dauer.pdf', transparent=True, bbox_inches='tight')
# Save summary dataframes
df_summary.to_csv('./output/probdauer_summary.csv')
df_diff_summary.to_csv('./output/probdiffdauer_summary.csv')
