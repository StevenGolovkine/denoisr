################################################################################
#                     Utils functions for the package                          #
################################################################################

#' Convert \code{funData} objects into right format for this package
#' 
#' @param data An object of the class \code{funData::funData}
#' @param norm Boolean, if TRUE, the sampling are normalized on [0, 1]
#' 
#' @return A list, where each element represents a curve. Each curve is defined
#'  as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#'  @export
transform_funData <- function(data, norm = TRUE){
  t <- funData::argvals(data)[[1]]
  x <- funData::X(data)
  
  data_list <- list()
  for(i in 1:funData::nObs(data)){
    if(norm) t <- (t - min(t)) / (max(t) - min(t))
    data_list[[i]] <- list(t = t, x = x[i,])
  }
  
  data_list
}

#' Convert \code{irregFunData} objects into right format for this package
#' 
#' @param data An object of the class \code{funData::irregFunData}
#' @param norm Boolean, if TRUE, the sampling are normalized on [0, 1]
#' 
#' @return A list, where each element represents a curve. Each curve is defined
#'  as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#'  @export
transform_irregFunData <- function(data, norm = TRUE){
  t <- funData()
  x <- data@X
  
  data_list <- list()
  if(norm) {
    data_list <- purrr::map2(t, x, 
                             ~ list(t = (.x - min(.x)) / (max(.x) - min(.x)), x = .y))
  } else {
    data_list <- purrr::map2(t, x, ~ list(t = .x, x = .y))
  }
  
  data_list
}

#' Convert \code{multiFunData} objects into right format for this package
#' 
#' @param data An object of the class \code{funData::multiFunData}
#' @param norm Boolean, if TRUE, the sampling are normalized on [0, 1]
#'
#' @return A list, where each element represents a curve. Each curve is defined
#'  as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#'  @export
transform_multiFunData <- function(data, norm = TRUE){
  
  data_list <- list()
  cpt <- 1
  for(fun_data in data){
    if(inherits(fun_data, 'funData')) {
      data_list[[cpt]] <- transform_funData(fun_data, norm)
    } else if(inherits(fun_data, 'irregFunData')) {
      data_list[[cpt]] <- transform_irregFunData(fun_data, norm)
    } else if(inherits(fun_data, 'multiFunData')){
      data_list[[cpt]] <- transform_multiFunData(fun_data, norm)
    } else{
      stop('Something went wrong with one of the functional data object!')
    }
    cpt <- cpt + 1
  }
  data_list
}
