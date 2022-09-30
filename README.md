# Clustering and annotation of TCR repertoires
This repository contains several tutorials that are linked to the chapter on *Clustering and annotation of TCR repertoires* of the *Computational Vaccine Design* book from Methods in Molecular Biology.

## Repository overview:

```
.
├── data
│   ├── examples
│   │   ├── P1_0.tsv
│   │   └── P1_15.tsv
│   ├── imgt_genes
│   │   ├── j_genes.tsv
│   │   └── v_genes.tsv
│   ├── tcrex_in
│   │   ├── P1_0
│   │   └── P1_15
│   └── tcrex_out
│       ├── P1_0
│       └── P1_15
├── notebooks
│   ├── part1_clustering.ipynb
│   ├── part2_splitfiles.ipynb
│   ├── part3_concatfiles.ipynb
│   ├── part4_tcrex_stats.ipynb
│   └── part5_cluster_annotation.ipynb
├── README.md
├── results
│   ├── clustcr
│   │   ├── P1_0_cluster_features.tsv
│   │   ├── P1_0_clusters.tsv
│   │   ├── P1_15_cluster_features.tsv
│   │   └── P1_15_clusters.tsv
│   ├── figures
│   ├── P1_0_clusters_tcrex.tsv
│   ├── P1_15_clusters_tcrex.tsv
│   └── tcrex
│       ├── P1_0_tcrex.tsv
│       └── P1_15_tcrex.tsv
├── src
│    └── tools.py
└── environment.yml
```

## How to use

You can follow along with the analysis workflow, either with the example data or your own, by running the code from the notebooks provided in the [./notebooks/](./notebooks/) folder. This requires the installation of a few dependencies. To do this, set up a new conda file, either through the use of the YAML file

```bash
conda env create -f environment.yml
```

or from scratch with the following sequence of commands:

```
conda create --name bookchapter
conda activate bookchapter
conda install python=3.9
conda install clustcr -c svalkiers -c bioconda -c pytorch -c conda-forge
conda install statsmodels
conda install ipykernel
```

If you wish to use your own data, populate the [./data/](./data/) directory with the desired files and adjust the corresponding paths in the notebook. Note that the V and J genes should not include allele information, as this may conflict with the input requirements for the TCRex software. For more information, visit the TCRex instructions page at [https://tcrex.biodatamining.be/instructions/](https://tcrex.biodatamining.be/instructions/).

## Contact

When encountering error, please write an issue at the [issues page](https://github.com/svalkiers/MMB_antigen_receptor_repertoire_tools/issues) with a detailed explanation of the error and provide an example of how it can be recreated.

