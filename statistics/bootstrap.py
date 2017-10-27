import matplotlib
from matplotlib import pyplot as plt

import pandas as pd
import numpy as np
import seaborn as sns
import itertools

# ########################################################################## #
#                        Bootstrap tests                                     #
# ########################################################################## #

def draw_bs_sample(data):
    """
    Draw a bootstrap sample from a 1D data set.
    """
    return np.random.choice(data, size=len(data))

def draw_bs_reps(data, stat_fun, size=10000):
    """
    Draw boostrap replicates computed with stat_fun from 1D data set.
    """
    return np.array([stat_fun(draw_bs_sample(data)) for _ in range(size)])

def frac2boolean(data, yes_label='dauer', no_label='non-dauer'):
    """
    Convert number of YES/NO dauer data into a boolean 1D dataset for bootstrap
    """
    yes, no = int(data[yes_label]), int(data[no_label])
    boolean = np.concatenate((np.ones(yes, dtype=int), np.zeros(no, dtype=int)))
    # Check that I'm not messing anything up
    assert sum(boolean) == yes and len(boolean) == yes+no
    # returning a tuple because for some reason pandas apply does not work with arrays
    return np.random.permutation(boolean)

def add_ci(data, labels_col, values_col, summ_col, c_interval, test_func, 
        bs_size=10000):
    """
    Compute confidence interval with bootstrap

    Arguments
    ---------
    data: dataframe
        tidy dataframe containing labeled sets to compare
    labels_col: string
        name of column containing set labels
    value_col: string
        name of column containing values
    summ_col: string
        name of column for summary values (e.g. mean)
    test_func: function
        function to summarize sets (e.g. np.mean)
    c_interval: tuple of ints
        confidence interval to compute (e.g. for 99%, c_interval=(1,99))
    bs_size: int
        number of bootstrap replicates

    Returns
    ---------
    pooled: dataframe
        dataframe with pooled data, confidence interval and bs samples
    """

    # pool data by strain
    pooled = data.groupby(labels_col)

    # summary dataframe
    summ_df = pd.DataFrame()
    summ_df[summ_col] = pooled[values_col].apply(test_func)

    # dataframe to store bootstraps
    bs_samples = pd.DataFrame()
    # lists for confidence intervals
    CI_, _CI = [], []

    # compute bootstrap samples and confidence intervals
    for l, group in pooled:
        bs_samples[l] = draw_bs_reps(group[values_col], test_func, bs_size)
        CI_.append(summ_df.loc[l].values[0] \
            - np.percentile(bs_samples[l], c_interval[0]))
        _CI.append(- summ_df.loc[l].values[0] \
            + np.percentile(bs_samples[l], c_interval[1]))
    summ_df['CI_'] = CI_
    summ_df['_CI'] = _CI
    return summ_df.reset_index(), bs_samples

def plot_ci(pooled, labels_col, summ_col, ylabel, title=''):
    """
    Plot data with confidence interval (output from add_ci func)
    """
    # plot with errorbars
    fig, ax = plt.subplots(1)
    ax.errorbar(x=np.arange(len(pooled)), y=pooled[summ_col], fmt='o', xerr=None,
            yerr=(pooled['CI_'], pooled['_CI']), capsize=10,capthick=1, color='#BA6102', alpha=0.8)

    labels = pooled[labels_col].tolist()
    ax.set_xticklabels(labels, rotation=60)
    tick_locator = matplotlib.ticker.LinearLocator(numticks=4)
    ax.xaxis.set_major_locator(tick_locator)
    plt.ylabel(ylabel)
    plt.title(summ_col)
    plt.tight_layout()
    plt.margins(x=0.05, y=0.05)
    sns.despine()

# ########################################################################## #
#                        Permutation test                                    #
# ########################################################################## #

def draw_perm_sample(x, y):
    """Generate a permutation sample."""
    concat_data = np.concatenate((x, y))
    np.random.shuffle(concat_data)
    return concat_data[:len(x)], concat_data[len(x):]

def draw_perm_reps(x, y, test_func, size=10000):
    """
    Generate array of permutation replicates.
    """
    out = np.empty(size)
    for i in range(size):
        x_perm, y_perm = draw_perm_sample(x, y)
        out[i] = test_func(x_perm, y_perm)
    return out

def diff_mean(x, y):
    """
    Compute difference of means
    """
    return np.abs(np.mean(x) - np.mean(y))

def diff_median(x, y):
    """
    Compute difference of means
    """
    return np.abs(np.median(x) - np.median(y))

def diff_passthresh(x, y):
    """
    Compute difference of means
    """
    return np.abs(np.sum(x >= 60) - np.sum(y >= 60)) 

def diff_var(x, y):
    """
    Compute difference of variance
    """
    return np.abs(np.var(x) - np.var(y))

def compute_pval(s1, s2, test_func=diff_mean, stat_func=draw_perm_reps, size=10000):
    """
    Compute p-value from comparing datasets s1 and s2 based on test_func
    using stat_func
    """
    # Compute test statistic for original data set
    diff = test_func(s1, s2)
    # Draw replicates and compute test statistic
    perm_reps = stat_func(s1, s2, test_func, size=size)
    # Compute p-value
    pval = np.sum(perm_reps >= diff) / len(perm_reps)
    return pval

def pairwise_pvals(data, labels_col, values_col, test_func=diff_mean, stat_func=draw_perm_reps, size=10000):
    """
    Compute pairwise p-values from statistical comparison using test_func and
    stat_func of entry sets in data

    Arguments
    ---------
    data: dataframe
        tidy dataframe containing labeled sets to compare
    labels_col: string
        name of column containing set labels
    value_col: string
        name of column containing values
    stat_func: function
        function to compute p-values (e.g. permutation)
    test_func: function
        function to compare sets (e.g. difference of means)

    Returns
    ---------
    p_vals: dataframe
        df with p_values from statistical comparison
    """

    # lists for p-value dataframe
    set1, set2, pval = [], [], []

    for (s1_label, s1), (s2_label, s2) in \
        itertools.combinations_with_replacement(data.groupby(labels_col), 2):
        pair = s1_label +'_'+ s2_label
        # Compute test statistic for original data set
        # Draw replicates and compute test statistic
        # Compute p-value
        _pval = perm_pval(s1[values_col].values, s2[values_col].values,
                        test_func, stat_func, size)
        # save it
        set1.append(s1_label)
        set2.append(s2_label)
        pval.append(_pval)

    # Store p-values in dataframe
    p_vals = pd.DataFrame(columns = ['set1','set2','p-value'])
    p_vals['set1'] = set1
    p_vals['set2'] = set2
    p_vals['p-value'] = pval
    return p_vals

def triang_heatmap(p_vals, title=''):
    """
    Plot triangular heatmap of p-values
    p_vals: dataframe
    """
    # pivot dataframe to make heatmap, transpose to align with labels
    p_vals = p_vals.pivot(index='set1', columns='set2', values='p-value').transpose()

    # plot heatmap of p-values by pairwise difference of percent dauer
    fig, ax = plt.subplots(1)
    sns.heatmap(p_vals, annot=p_vals, robust=True, ax=ax, fmt='0.2f', cmap=plt.cm.viridis, cbar=True, square=True)
    ax.set_xticklabels([name for name in p_vals.columns], rotation=60, fontsize=20)
    ax.set_yticklabels(reversed([name for name in p_vals.columns]), rotation=0, fontsize=20)
    ax.set_ylabel('')
    ax.set_xlabel('')
    plt.title(title)
    plt.tight_layout()
