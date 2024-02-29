################################  main run #####################################
#```{r clean_env}

rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
#```
#R CMD INSTALL --preclean --no-multiarch --with-keep.source IP4CI
library(IP4CI)
# run all required pkgs
#```{r run_packages}
run_packages()
#```
################################################################################

# Amy dataset

# variables needed by the main scripts
db_opt = 1
if(db_opt == 1)
{
  expr_dir = '~brain'
  obj1_fname='m.RDS'
  obj2_fname='h.RDS'
  typeColumnToUse.list=c('xl','yl') # later replaced by 'type'
  assaySlotToUse.list=list('xa','ya')
}
if(db_opt == 2)
{
  expr_dir = '~pancreas'
  obj1_fname='m.RDS'
  obj2_fname='h.RDS'
  typeColumnToUse.list=c('xl','yl') # later replaced by 'type'
  assaySlotToUse.list=list('xa','ya')
}
id1='mouse'
id2='human'
id.list=c(id1,id2)
obj.fname.list =c(obj1_fname,obj2_fname)
processObjOpt = T
genes2keep='commonVarG'
convertG.list = c(FALSE,FALSE)
dataType = 'norm' # use normalized obj data slot
geneType = 'var' # use variable genes
ct = as.numeric(args[1])
message('ct:  ',ct)
pathwaydatabase = 'reactome.db'
filter_opt = 'pvalue'
filter_cutoff = 0.05
sample_cutoff = 0.3
topRanked = 100

runIP4CI_ct(expr_dir, id.list,obj.fname.list,
         convertG.list,
         typeColumnToUse.list,assaySlotToUse.list,
         genes2keep,
         processObjOpt,
         pathwaydatabase,dataType,geneType,
         filter_opt,filter_cutoff,sample_cutoff,
         topRanked,
         ct)
