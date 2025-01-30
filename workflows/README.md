## workflows
### Execute them in the following order
1. `01-obtain_data`: Download data
2. `02-preprocess_data`: Preprocess downloaded data for the project
3. `03-new_empirical_network`: Inference of SF-exon networks
4. `04-sf_programs_in_differentiation`: Validation of SF activity estimation with VIPER and SF-exon networks
5. `06-carcinogenic_switch_regulation`: Identifying a recurrent cancer splicing program
6. `07-prepare_submission`: Prepare supplementary tables for submission

### Stucture of each workflow
Inside each workflow folder you'll find a `README.md` file explaining the corresponding workflow outline as well as Snakefile(s) to run the workflows. In most cases there is also a `<workflow_dir>/scripts` folder containing helper scripts for those specific workflows.

```{shell}
workflows/<workflow_dir>/
├── <workflow_snakefile>.smk
├── README.md
└── scripts
    └── <helper_script>
```