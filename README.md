# Over-optimism in unsupervised microbiome analysis: Insights from network learning and clustering

## Project description 

In recent years, unsupervised analysis of microbiome data, such as microbial network
analysis and clustering, has increased in popularity. Many new statistical and
computational methods have been proposed for these tasks. This multitude of analysis
strategies poses a challenge for researchers, who are often unsure which method(s) to
use and might be tempted to try different methods on their dataset to look for the
"best" ones. However, if only the best results are selectively reported, this may cause
over-optimism: the "best" method is overly fitted to the specific dataset, and the results
might be non-replicable on validation data. 

In our illustrative study, we aim to quantify
over-optimism effects in the context of unsupervised microbiome analysis. We model the approach of a hypothetical
microbiome researcher who undertakes three unsupervised research tasks: clustering of
bacterial genera, hub detection in microbial networks, and differential microbial network
analysis. While these tasks are unsupervised, the researcher might still have certain
expectations as to what constitutes interesting results. We translate these expectations
into concrete evaluation criteria that the hypothetical researcher might want to
optimise. We then randomly split an exemplary dataset from the American Gut Project
into discovery and validation sets multiple times. For each research task, multiple
method combinations (e.g. methods for data normalization, network generation and/or
clustering) are tried on the discovery data, and the combination that yields the best
result according to the evaluation criterion is chosen. While the hypothetical researcher
might only report this result, we also apply the "best" method combination to the
validation dataset. The results are then compared between discovery and validation
data.

Here are the instructions for reproducing our results. 

## Setup and data

First, clone this repository. We used data from the American Gut Project (AGP). The OTU table and the metadata are stored in the files
`otu_table__BODY_HABITAT_UBERON_feces_json.biom` and `metadata__BODY_HABITAT_UBERON_feces__.txt` which can be downloaded from [Zenodo](https://doi.org/10.5281/zenodo.6652711). Put these files in the `data` folder. 

## Running the code 

The software versions we used are listed at the end of the README. 

All commands should be run in an R console with the base directory of this repository as the working directory.

### Preprocessing 

```
source("code/preprocess/preprocessing.R")
```
Starting from the OTU table, this script performs preprocessing (sample filtering, OTU filtering and agglomeration to genus level). The result is the phyloseq object `ag.genus.rds`. Instead of running the preprocessing script, you can download the file
`ag.genus.rds` directly from [Zenodo](https://doi.org/10.5281/zenodo.6652711) and put it in the `data` folder.


### Optimization on the discovery data and validation on the validation data 

The following scripts perform, for each research task in turn, the following steps:
1. For varying sample sizes, discovery and validation sets are sampled 50 times each.
2. For each sampling, different method combinations are tried on the discovery data, and the "best" combination is chosen.
3. The "best" method combination is applied to the validation data.

The results are stored in the `results` folder. 

Clustering bacterial genera:
```
source("code/50 splits/researchtask1.R")
```

Hub detection:
```
source("code/50 splits/researchtask2.R")
```

Differential network analysis: 
```
source("code/50 splits/researchtask3.R")
```

### Generating descriptive statistics and plots

```
source("code/analyse results/generate_statistics_and_plots.R")
```
This script generates the tables in the `descriptive statistics` folder (Tables 1-3 in the manuscript), as well as the figures in the `plots` folder. 

## Software versions

We used R version 4.0.4 and the following R package versions:

```
ggnewscale_0.4.5      ggbeeswarm_0.7.0.9000 gridExtra_2.3         ggplot2_3.3.3        
orca_1.1-1            reticulate_1.24       dplyr_1.0.2           readxl_1.3.1         
mclust_5.4.7          igraph_1.2.6          mixedCCA_1.4.6        MASS_7.3-53.1        
NetCoMi_1.0.2         SpiecEasi_1.1.0       phyloseq_1.34.0      

```

Additionally, for research task 1 (clustering), we used Python version 3.6.13 und numpy version 1.19.5. 
