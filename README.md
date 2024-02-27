# IP4CI
IP4CI: an interpretable pathway-based computational method for cross-species cell-type identification from single cell data

IP4CI is an R package designed for cell-type annotation of single-cell RNA-seq data with its applicability for cross-species datasets. It aims to enable users to align single-cell datasets for reliable cell type annotation with improved biological interpretation. IP4CI takes as input two scRNA-seq datasets, along with their cell-type labels, and a list of gene-sets or pathways. First, it encompasses the standard pre-processing workflow for scRNA-seq data as in Seurat. Compared to other cell-type annotation tools, IP4CI transforms each gene expression matrix into a pathway activity matrix. To improve interpretability, IP4CI can extract biological information by providing users with dimensional reduction, cell clustering, cell annotation, and ranking of biological discoveries of pathways and its associated gene members.
IP4CI has been successfully installed on Mac OS X, and Linux. It is hosted on GitHub to view and clone the repository. Instructions, documentation, and tutorials can be found at (https://github.com/harvardinformatics/IP4CI/). Improvements and new features will be added on a regular basis, users can post on the GitHub page with any questions or inquiries to contribute.


