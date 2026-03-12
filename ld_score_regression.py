import pandas as pd
import gzip
import shutil
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns
import statsmodels.api as sm
import scipy

def load_ld_scores(path_to_ld_scores):

    ldscores = []
    for i in range(1,23):
        ldscore_file = pd.read_csv(path_to_ld_scores + f'/LDscore.{i}.l2.ldscore.gz',sep='\s+')
        
        ldscores.append(ldscore_file)
    ldscore_df = pd.concat(ldscores).reset_index(drop=True)

    ldscore_df.index = ldscore_df['CHR'].astype(str) + '.' + ldscore_df['SNP'].astype(str)
    return ldscore_df

def load_p_values(gwas_results_df, chr_col_name, snp_col_name):

    df = gwas_results_df.copy()
    df.index = df[chr_col_name].astype(str) + '.' + df[snp_col_name].astype(str)
    return df

def integrate_data(gwas_results_df, chr_col_name, snp_col_name, pval_col_name,
                   path_to_ld_scores):

    ldscore_df = load_ld_scores(path_to_ld_scores)
    df = load_p_values(gwas_results_df, chr_col_name, snp_col_name)
    common_indices = list(set(ldscore_df.index) & set(df.index))

    common_df = pd.concat([df[[chr_col_name,snp_col_name,pval_col_name]].loc[common_indices],
           ldscore_df[['L2']].loc[common_indices]],axis=1)

    common_df = common_df[common_df[pval_col_name] != 0].copy()

    common_df['chisq'] = scipy.stats.chi2.isf(common_df[pval_col_name], 1)

    common_df = common_df.dropna()
    
    common_df['L2 Score Bin'] = (pd.qcut(common_df['L2'], q=50, labels=False) + 1)

    score_bins = common_df[['L2','chisq','L2 Score Bin']].groupby('L2 Score Bin').mean()

    return score_bins, len(common_indices)

def ld_score_regression(l2_score, chisq, sample_size, num_snp, outpath):

    weights = 1 / l2_score

    X = l2_score
    y = chisq
    X = sm.add_constant(X)
    model = sm.WLS(y, X, weights=weights)
    results = model.fit()

    r2 = results.rsquared
    intercept = results.params['const']
    slope = results.params['L2']
    a = (intercept - 1)  / sample_size
    h2 = slope  * num_snp / sample_size

    ax=sns.scatterplot(x=l2_score, y=chisq,hue=weights,palette='coolwarm')
    ax.axline(xy1=(0, intercept), slope=slope, color='r')
    handles, labels = ax.get_legend_handles_labels()

    extra_label = f"r\u00B2 = {r2:.3e}"
    extra_handle = mpatches.Patch(color='none', label=extra_label)

    handles.append(extra_handle)
    labels.append(extra_label)

    extra_label = f"a = {a:.3e}"
    extra_handle = mpatches.Patch(color='none', label=extra_label)

    handles.append(extra_handle)
    labels.append(extra_label)

    extra_label = f"h\u00B2 = {h2:.3e}"
    extra_handle = mpatches.Patch(color='none', label=extra_label)

    handles.append(extra_handle)
    labels.append(extra_label)

    ax.legend(handles, labels, title='Regression Weight',loc='best')

    ax.set_ylabel('-log\u2081\u2080p-value')
    ax.set_xlabel('LD Score Bin')

    ax.set_title(f'LD Score Regression')

    plt.savefig(outpath,bbox_inches='tight')

