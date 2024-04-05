# Pipeline
## Data
- [X] from scPerturb
    - ReplogleWeissman2022: K562 and RPE1
- [X] from VIPER splicing
    - gene expression and splicing of SF perturbations
    - [X] cancer splicing programs
    
## Analysis
### Splicing factor activity analysis using single cell perturbation datasets
1. [X] Evaluate how well using empirical SF-exon and SF-gene networks recapitulates SF activity
    - SF-gene networks recapitulate SF-exon activities
    - FIGURES
        - [X] recall SF-exon and SF-gene vs random
        - [X] distribution of correlations between SF-exon and SF-gene activities in CCLE dataset
        - [X] best and worse examples of correlations in a scatter plot
2. [ ] Evaluate robustness empirical SF-gene networks
    - FIGURES
        - [ ] recall SF-gene pruned randomly vs random 
        - [ ] recall SF-gene pruned with different logFC thresholds
3. [ ] Evaluate how well activities in the same perturbation experiment in bulk and singl-cell RNA seq correlate
    - activities measured with bulk and single cell are equivalent
    - FIGURES
        - [X] distribution of correlations between activities in bulk (ENCORE K562) and single cell (Replogle K562)
        - [X] best and worst examples of correlations in scatter plot connected to the distribution panel
        - [X] evaluate reproducibility between single cell datasets using Replogle K562 genome-wide and essential as independent experiments
        - [ ] Does the reproducibility change with bulk CRISPR KO vs shRNA KD? Replogle single cell was CRISPRi
        - [ ] Correlation between changes in gene expression in bulk vs single cell. Is it comparable to the correlation with activity?
        - [ ] Uncertainty of estimated protein activities: number of expressed/detected genes in single cells vs activity correlations with bulk
        - in single cell, should we compute median across replicates?

### Charting the regulation of the cancer splicing program
4. [ ] Confirm that there is a switch of cancer splicing program during carcinogenesis
    - SC tumor initiation in intestine: https://www.nature.com/articles/s41392-022-00881-8#Sec9

5. [ ] Compute splicing factor activities for all perturbations in Replogle K562 genome-wide and RPE1
    - explore which KOs cause the switch of cancer programs
    - FIGURES
        - [ ] plot fold change between programs in K562 and RPE1 in a scatter plot -> correlated? One is already cancer the other immortalized
        - [ ] GSEA in K562 and RPE1 of the fold changes to SF