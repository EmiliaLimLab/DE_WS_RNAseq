# DE_WS_RNAseq

Here we present the scripts and data used in the RNA-seq analysis of our submitted manuscript entitled "Woodsmoke and diesel exhaust: distinct transcriptomic profiles in the human airway epithelium" in AJRCMB.

## Raw sequence data, processing, and gene quantification

RNA-seq data can be found on the European Genome-Phenome Archive (EGA) under accession EGAXXXXXX.

QC and gene quantification was performed using the [`nf-core/rnaseq`](https://github.com/nf-core/rnaseq) nextflow pipeline. Relevant files can be found in `data/nfcore-rnaseq`.

## Analysis scripts and corresponding data

Code for all RNA-seq analyses presented in the paper can be found in their respective subdirectories under `scripts/`. Specifically, code is provided for:

* Sample PCA: `scripts/01_PCA.R`
* Preparing CIBERSORTx inputs: `scripts/02_CIBERSORTx`
* Differential expression analysis: `scripts/03_DEA.R`
* Gene set enrichment/pathway analysis: `scripts/04_GSEA.R`

For convenience, outputs for majority of these scripts are also available in `data/`. In particular:

* `data/QC/GSC-2581_mergeFA_PCA.tsv`: contains sample PCA data ouptut by `scripts/01_PCA.R`
* `data/DEA/GSC-2581_mergeFA_deseq-inclBasal-*.tsv`: contains DEA results output by `scripts/03_DEA.R` which compares DE vs. FA, WS vs. FA, and IC vs. FA
* `data/GSEA/GSC-2581_mergeFA_deseq-inclBasal-*_GSEA-KEGG.tsv`: contains GSEA results output by `scripts/04_GSEA.R` from DE vs. FA and WS vs. FA ranked gene lists

If you are interested in running these analysis scripts, relevant input files will be automatically read in and results will be written to `output/`.

### Resource requirements for cell deconvolution with CIBERSORTx

The core scRNA-seq dataset from the [HCLA](https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293) was downloaded and randomly downsampled. This subset was then used to generate a signature matrix for ALI cultured epithelial cells. A script is provided to perform the downsampling (`scripts/02_CIBERSORTx/02A_SubsampleHCLA.R`); however, this requires about 35 GB of RAM. 

Alternatively, all initial and intermediary files needed to run CIBERSORTx can be found at  https://www.bcgsc.ca/downloads/lim_lab/DE_WS_RNAseq/. Refer to the `README` for more information on each of the files.

## Figure generation
The script used to generate plots Figure 1A-D: `scripts/plot_figure1A-D.R`.
