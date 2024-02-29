#' @title  createDir
#' @description create directory folder if it does not already exists
#' @param dirName directory folder name
#' @examples
#' dirName='mainDir'
#' createDir(dirName)
# then run roxygenize()

#' @export
createDir <- function(dirName)
{
  if (!dir.exists(dirName)){
    dir.create(dirName)
  } else {
    print(paste0("Dir ", dirName, " already exists!"))
  }
}
