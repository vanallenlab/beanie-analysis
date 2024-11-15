# BEANIE

## Dissecting tumor cell programs through group biology estimation in clinical single-cell transcriptomics


## Paper Title

[Preprint](https://www.biorxiv.org/content/10.1101/2021.10.22.465130v2.full) | [Conference Submission](https://aacrjournals.org/cancerres/article/83/7_Supplement/1120/722439/Abstract-1120-Dissecting-tumor-cell-programs) | [Cite Us](link)

This repository contains data and code for reproducing analysis presented in the paper. To use the BEANIE method for your own analysis, please see installation instructions on [this github repo](https://github.com/vanallenlab/beanie). 

Main Findings:
<details>
    
  <summary>
      <b>With the growth of clinical cancer single-cell RNA sequencing (scRNA-seq) studies, robust differential expression methods for case/control analyses (e.g., treatment responders vs. non-responders) using gene signatures are pivotal to nominate hypotheses for further investigation. However, current commonly used methods for these analyses present multiple challenges for interpretation: (i) they produce a large number of false positives; (ii) they do not adequately represent the patient-specific hierarchical structure of clinical scRNA-seq data—a problem especially prominent for the tumor compartment due to the higher patient specificity exhibited by tumor cells; and (iii) they do not account for sample-driven confounders that arise from the variability in the number of recovered cells per sample and patient biology. Here, we present a novel nonparametric statistical method, BEANIE, to address these issues for investigating differential expression of gene signatures between clinically relevant groups within scRNA-seq data.</b>
  </summary>

1. We benchmark BEANIE's performance relative to conventional methods used in differential gene signature expression analysis scRNA-seq datasets using simulated datasets, and find superior sensitivity-specificity trade-off.
2. We demonstrate BEANIE's use in real-world clinical datasets in breast cancer, lung cancer and melanoma.
3. Overall, BEANIE provides a methodological strategy to inform biological insights into unique and shared differentially expressed gene signatures across different tumor states, with utility in single-study, meta-analysis, and cross-validation across cell types.
</details>

<!-- ### Updates
 -->

### Reproducing Results
The main code files are organized as follows in this repository.

```bash
beanie-analysis/
	├── analysis
    		├── brca
    		├── luad
    		└── mel
	├── data
    		├── brca
    		├── lung
    		├── melanoma
    		├── resource_utilization
    		└── signatures
	├── resource_utilization
    		├── results
    		└── scripts
	├── simulations/scripts
    		├── create_simulations.py
    		├── simulation_baseline.py
    		├── simulation_benchmark.py
    		└── simulation_perturbation.py
	└── process_list_autogen.csv
```

#### Setup Instructions
To reproduce the analysis, re-create the conda environment.

```
conda env create --file=beanie_analysis_environment.yaml
```

#### Dataset
The dataset pre-processing scripts are present in the corresponding directory within the data folder. Raw and processed data is also available for download within the same repository.

#### Reproducing Results

**Simulations** : Run scripts in the following order (replace the location of seed dataset file in create_simulations.py with your own seed data to be used).
```
python create_simulations.py
python simulation_baseline.py
python simulation_perturbation.py
python_benchmark.py
```

**Results** : Jupyter notebooks are available in the `analysis` folder, within each cancer type subfolder, to use the raw data, identify tumor states and run analysis. Each folder is organized as follows-

	├── brca
        ├── inputs --> contains adata.h5ad and metad.csv files for each of the tumor states in the analysis.
        ├── supplementary tables --> consolidated results presented in the supplementary tables of the paper.
        ├── outputs --> individual full result files from beanie analysis and benchmark methods.
            ├── ts1_beanie_de.csv --> BEANIE results for DE analysis on tumor state 1.
            ├── ts1_benchmark_de.csv --> Benchmark methods results for DE analysis on tumor state 1.
            ├── ts1_noise_benchmark_de.csv --> Immune cell surface marker results for DE analysis on tumor state 1 using benchmark methods
            ├── ts1_noise_de.csv --> Immune cell surface marker results for DE analysis on tumor state 1 using BEANIE
            ├── ts1_sig_scores.csv --> signature scores  
            └── ...
        └── ... --> Jupyter notebooks for identifying tumor states from processed data.

### Issues
Please open new issue threads specifying the issue with the codebase or report issues directly to sjohri@g.harvard.edu.

### License and Terms of use
The source code is licensed under the GPL-2.0 license, which can be found under the `LICENSE` file. 

### Citation

If you find this work useful, please cite our paper -

>Johri S., Bi K., Titchen B., Fu J., Conway J., Crowdis J., Vokes N., Fan Z, Fong L., Park J., Liu D., He MX., Van Allen E. (2021) Dissecting tumor cellprograms through group biology estimation in clinical single cell transcriptomics. biorxiv.