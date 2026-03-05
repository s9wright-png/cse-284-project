# cse-284-project
LD Score Regression Implementation

CSE 284 Project, Winter 2026

Galen Heuer and Samuel Wright

LD Scores downloaded from 1000 Genomes: https://zenodo.org/records/10515792

GWAS on late-onset Alzheimer's (Kunkle et al., Nature Genetics 2019) to be used with published summary statistics: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007511/

# Example usage

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

Run the LD score regression function, specifying the 'L2' LD score and the 'nlogp' columns of the score_bins dataframe obtained from the last step. Also specify the sample size of the GWAS (94,437 per Kunkle et al.), as well as the desired output path for the LD score regression figure, which reports the R<sup>2</sup>, a (contribution of confounding) and heritability.

```python
ld_score_regression(score_bins['L2'],score_bins['nlogp'],
                   94437,num_snp,'/out/path/ld_score_graph.pdf')
```

<img width="584" height="455" alt="image" src="https://github.com/user-attachments/assets/c14fbe31-634d-49c1-993f-6e97786d3aef" />

# Test on simulated data (TO BE FINALIZED)

Run LD score regression function on simulated data to see if the tool can distinguish traits with inflated p-values vs. true polygenic architectures.

Load the simulated data from /simulated_data in this repository:
```python
inflated = pd.read_csv('/cse-284-project/simulateddata/sim_inflated.tsv', sep='\t')
polygenic = pd.read_csv('/cse-284-project/simulateddata/sim_polygenic.tsv', sep='\t')

```

# Next steps

Additionally, benchmark performance of our function against LDSC tool.
