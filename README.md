**IP4CI**: an interpretable pathway-based computational method for cross-species cell-type identification from single cell data

**input**: two seurat objects, along with their metadata
Note: All the gene names should be converted to human gene

**installation**:
*From terminal type:*
R CMD INSTALL --no-multiarch --with-keep.source IP4CI
*From R Studio:*
locate IP4CI folder and open project file IP4CI.Rproj

*Demo using pancreas data*
1. create a folder to hold experiment data & results
2. create sub-folder to hold data/ and move data into it
   for single cell data , please make sure they are Two Seurat Objects, one for each dataset, where counts & metadata are available in each object  


