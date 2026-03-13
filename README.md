# LD Score Regression Implementation

CSE 284 Project, Winter 2026

Galen Heuer and Samuel Wright

LD Scores downloaded from 1000 Genomes: https://zenodo.org/records/10515792

GWAS on late-onset Alzheimer's (Kunkle et al., Nature Genetics 2019) to be used with published summary statistics: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007511/

# Example usage on publicly available late-onset Alzheimer's dataset

Download the LD score file "1000G_Phase3_ldscores.tgz" from the /data folder in this repository.

Download "Kunkle_etal_Stage1_results.txt" file from the link above for the GWAS summary statistics from the Kunkle et al. paper and read the data into a python notebook as a Pandas dataframe.

Load in the summary statistics as a Pandas dataframe:

```python
import sys
import pandas as pd

alz = pd.read_csv('/path/to/pvalues/Kunkle_etal_Stage1_results.txt', sep=' ')
```

Load in integrate_data and ld_score_regression function from the "ld_score_regression.py" file available in this repository:

```python
sys.path.insert(1, '/path/to/py/file/ld_score_regression.py')

from ld_score_regression import integrate_data,ld_score_regression
```

Run the integrate_data function, by providing the dataframe you loaded in previously, the names of the chromosome, snp name and p-value columns of this dataframe, respectively, as well as the folder file path for the LD score file you downloaded. This function retrieves a binned LD score dataframe along with the number of SNPs that were common between the LD score file and the summary statistics, to be used to calculate the heritability in the next step.

```python
score_bins, num_snp = integrate_data(alz,'Chromosome','MarkerName','Pvalue',
               '/path/to/ldscore/folder/LDscore')
```

Run the LD score regression function, specifying the 'L2' LD score and the 'chisq' columns of the score_bins dataframe obtained from the last step. Also specify the sample size of the GWAS (63,926 per Kunkle et al.) and number of snps (11,480,632 in the original dataframe from Kunkle et al.), as well as the desired output path for the LD score regression figure, which reports the R<sup>2</sup>, a (contribution of confounding) and heritability.

```python
ld_score_regression(score_bins['L2'],score_bins['chisq'],
                   94437,num_snp,'/out/path/ld_score_graph.pdf')
```

<img width="585" height="455" alt="image" src="https://github.com/user-attachments/assets/32f586e1-cee8-425f-832c-7c67805f3ce4" />


# Test on simulated data

To validate the tool’s ability to distinguish between true polygenic signals and population inflation (as per reviewer feedback), we use two simulated datasets generated from 1000 Genomes (Chromosome 22) genotypes.

### Execution

Run the benchmarks by importing the module functions. This will generate regression plots and calculate the SNP-heritability ($h^2$) and the intercept-based bias ($a$).

```python
from ld_score_regression import integrate_data, ld_score_regression

# Scenario A: Polygenic Trait
# Simulated with h2 = 0.5 using real LD structure
df_poly, n_snps = integrate_data(pd.read_csv('simulated_data/sim_polygenic.tsv', sep='\t'), 
                                 'CHR', 'SNP', 'Pvalue', 'simulated_data/sim_ldscores.txt')
ld_score_regression(df_poly['L2'], df_poly['nlogp'], 2500, n_snps, 'plot_polygenic.png')

# Scenario B: Inflated Trait
# Simulated with h2 = 0 but injected p-value bias
df_infl, n_snps = integrate_data(pd.read_csv('simulated_data/sim_inflated.tsv', sep='\t'), 
                                 'CHR', 'SNP', 'Pvalue', 'simulated_data/sim_ldscores.txt')
ld_score_regression(df_infl['L2'], df_infl['nlogp'], 2500, n_snps, 'plot_inflated.png')
```
### Interpreting results

Compare your generated plots and statistics against the table below. This validation confirms that the tool correctly attributes signal to either genetic architecture or confounding factors.

| Metric | Polygenic Scenario | Inflated Scenario | Significance |
| :--- | :--- | :--- | :--- |
| **Slope ($h^2$)** | **High & Positive** | **Near Zero (Flat)** | A high slope indicates the tool successfully captured the simulated heritability ($h^2$). |
| **Intercept ($a$)** | **$\approx 1/N$** | **Significantly $> 1/N$** | A high intercept relative to $1/N$ proves the tool identifies non-genetic bias (inflation). |
| **Visual Trend** | Upward diagonal | Horizontal / Flat | Confirms whether the $-\log_{10}(P)$ values are correlated with LD scores. |

# Benchmarking against SUMRHE

To benchmark the performance of our tool, we tested the SUMRHE method which calculates partitioned heritability by performing LD score regression on unweighted and weighted LD scores. The Github for this tool (https://github.com/sriramlab/SUMRHE) provides example usage for this tool with toy datasets, including an LD score table with weighted and unweighted LD scores, so to compare performance with our tool, we tested the performance of each method on this dataset. The SUMRHE method calculated a heritability of 0.05908 using the unweighted LD scores and 0.21531 using the weighted LD scores, for a sum of 0.27439 for total heritability. To get this value for heritability we followed the tutorial specified in the README for this tool using the function:

```python
python3 ../src/sumrhe.py --pheno ./sim_50k_h2_0.25_p_0.01.sumstat \
                  --ldscores ./double_uniform_10k_stoc_k100.gw.ldscore.gz \
                  --annot ./double_uniform_0.2.txt \
                  --verbose \
                  --njack 1000
```

To implement this method on our tool, we downloaded the sim_50k_h2_0.25_p_0.01.sumstat and double_uniform_10k_stoc_k100.gw.ldscore.gz and read these into a python notebook. For our method to work properly, we need to use positive weights only for the statsmodels.api implementation of WLS, so we removed rows with negative LD scores from the simulated LD score data as well.

```python
sumstat = pd.read_csv("/Users/samuelwright/Downloads/sim_50k_h2_0.25_p_0.01.sumstat",sep=' ')

ldscore_toy = pd.read_csv('/Users/samuelwright/Downloads/double_uniform_10k_stoc_k100.gw.ldscore',sep='\t')

ldscore_toy = ldscore_toy[(ldscore_toy['L2_0'] > 0) & (ldscore_toy['L2_1'] > 0)]
```
This data was structured somewhat differently from the datasets we had built the tool around, necessitating some manual formatting of the data instead of using the integrate_data function. For example, instead of including a p-value column, this data included a Z score column, which we took the square of this column to obtain a chi-square score.

```python
sumstat.index = sumstat['SNP']

ldscore_toy.index = ldscore_toy['SNP']

combined = pd.concat([ldscore_toy[['SNP','L2_0','L2_1']],sumstat[['N','Z']]],axis=1).dropna()
```

We then make bins as in our implementation for the unweighted LD scores.

```python
combined['Score Bins_0'] = pd.qcut(combined['L2_0'], q=50, labels=False) + 1

combined['L2'] = combined['L2_0']

score_bins = combined[['L2','Z2','Score Bins_0']].groupby('Score Bins_0').mean()
```

We can then use our function for performing and plotting LD score regression, using the specified sample size in the SUMRHE Github.

```python
ld_score_regression(score_bins['L2'], score_bins['Z2'], 10000, len(combined),'/out/path/ld_score_graph.pdf')
```

Which results in a heritability of 0.1873.

We repeat this process for the weighted LD scores.

```python
combined['Score Bins_1'] = pd.qcut(combined['L2_1'], q=50, labels=False) + 1

combined['L2'] = combined['L2_1']

score_bins = combined[['L2','Z2','Score Bins_1']].groupby('Score Bins_1').mean()

ld_score_regression(score_bins['L2'], score_bins['Z2'], 10000, len(combined),'/out/path/ld_score_graph.pdf')
```

Which results in a heritability of 0.1203. While the heritabilities from the individual regressions on the weighted and unweighted LD scores were quite different from SUMRHE's implementation, it's notable that the sum of our heritabilities of 0.3076 is quite close to SUMRHE's sum of 0.27439. While this result may be incidental since our model has to remove negative LD scores to run properly as well as SUMRHE's likely more sophisticated method, it is interesting that the total heritabilities are similar, perhaps indicating that the core method is still capturing the modes variability we are hoping to.
