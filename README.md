# Pipeline
## Data
- [X] from scPerturb
    - ReplogleWeissman2022: K562 and RPE1
- [X] from VIPER splicing
    - gene expression and splicing of SF perturbations
    - [X] cancer splicing programs
    
## Analysis
### Splicing factor activity analysis using single cell perturbation datasets
1. [ ] Evaluate how well using empirical SF-exon and SF-gene networks recapitulates SF activity
    - SF-gene networks recapitulate SF-exon activities
    - FIGURES
        - recall SF-exon and SF-gene vs random
        - distribution of correlations between SF-exon and SF-gene activities in CCLE dataset
        - best and worse examples of correlations in a scatter plot
2. [ ] Evaluate how well activities in the same perturbation experiment in bulk and singl-cell RNA seq correlate
    - activities measured with bulk and single cell are equivalent
    - FIGURES
        - distribution of correlations between activities in bulk (ENCORE K562) and single cell (Replogle K562)
        - best and worst examples of correlations in scatter plot connected to the distribution panel
        - evaluate reproducibility between single cell datasets using Replogle K562 genome-wide and essential as independent experiments
        - in single cell, should we compute median across replicates?

### Charting the regulation of the cancer splicing program
3. [ ] Compute splicing factor activities for all perturbations in Replogle K562 genome-wide and RPE1
    - explore which KOs cause the switch of cancer programs
    - FIGURES
        - plot fold change between programs in K562 and RPE1 in a scatter plot -> correlated? One is already cancer the other immortalized
        - GSEA in K562 and RPE1 of the fold changes to SF