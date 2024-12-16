# Analyzing Adult Mouse Brain scATAC-seq
Compiled on August 20, 2024

## Reference
The analysis is based on the workflow outlined by the Stuart Lab. For more details, refer to the original vignette:  
[Mouse Brain scATAC-seq Vignette](https://stuartlab.org/signac/articles/mouse_brain_vignette)

## Data Sources
The following data files are used in this analysis:
- **Filtered peak-barcode matrix (.h5)**
- **Single-cell metadata (.csv)**
- **Fragments file (.tsv.gz)**
- **Fragments index file (.tsv.gz.tbi)**

These files can be downloaded using the following `wget` commands:

```bash
wget http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_filtered_peak_bc_matrix.h5
wget http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_singlecell.csv
wget http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_fragments.tsv.gz
wget http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_fragments.tsv.gz.tbi
