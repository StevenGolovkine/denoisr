################################################################################
#                     Utils functions for the package                          #
################################################################################

#' Convert \code{funData} objects into right format for this package
#' 
#' @param data An object of the class \code{funData::funData}
#' @param norm Boolean, if TRUE, the sampling are normalized on \eqn{[0, 1]}
#' 
#' @return A list, where each element represents a curve. Each curve is defined
#'  as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#'  @export
funData2list <- function(data, norm = TRUE){
  t <- funData::argvals(data)[[1]]
  x <- funData::X(data)
  
  data_list <- list()
  for(i in 1:funData::nObs(data)){
    if(norm) t <- (t - min(t)) / (max(t) - min(t))
    data_list[[i]] <- list(t = t, x = x[i,])
  }
  
  data_list
}

#' Convert comprehensive lists into \code{funData::funData} objects.
#' 
#' We assume that we \strong{know} that the curves are on the same interval.
#' 
#' @importFrom magrittr %>%
#'  
#' @param data_list A list, where each element represents a curve. Each curve is
#'  defined as list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'   }
#'   
#' @return  An object of the class \code{funData::funData}
#' @export
list2funData <- function(data_list){
  argvals <- data_list[[1]]$t
  obs <- data_list %>% 
    purrr::map_dfc("x") %>% 
    as.matrix() %>% 
    t()
  funData::funData(argvals = argvals, X = obs)
}

#' Convert \code{irregFunData} objects into right format for this package
#' 
#' @param data An object of the class \code{funData::irregFunData}
#' @param norm Boolean, if TRUE, the sampling are normalized on \eqn{[0, 1]}
#' 
#' @return A list, where each element represents a curve. Each curve is defined
#'  as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#'  @export
irregFunData2list <- function(data, norm = TRUE){
  t <- data@argvals
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

#' Convert comprehensive lists into \code{funData::irregFunData} objects.
#' 
#' @importFrom magrittr %>%
#'  
#' @param data_list A list, where each element represents a curve. Each curve is
#'  defined as list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'   }
#'   
#' @return  An object of the class \code{funData::irregFunData}
#' @export
list2irregFunData <- function(data_list){
  argvalsList <- data_list %>% purrr::map("t")
  obsList <- data_list %>% purrr::map("x")
  funData::irregFunData(argvals = argvalsList, X = obsList)
}

#' Convert \code{multiFunData} objects into right format for this package
#' 
#' @param data An object of the class \code{funData::multiFunData}
#' @param norm Boolean, if TRUE, the sampling are normalized on \eqn{[0, 1]}
#'
#' @return A list, where each element represents a curve. Each curve is defined
#'  as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#'  @export
multiFunData2list <- function(data, norm = TRUE){
  
  data_list <- list()
  cpt <- 1
  for(fun_data in data){
    if(inherits(fun_data, 'funData')) {
      data_list[[cpt]] <- funData2list(fun_data, norm)
    } else if(inherits(fun_data, 'irregFunData')) {
      data_list[[cpt]] <- irregFunData2list(fun_data, norm)
    } else if(inherits(fun_data, 'multiFunData')){
      data_list[[cpt]] <- multiFunData2list(fun_data, norm)
    } else{
      stop('Something went wrong with one of the functional data object!')
    }
    cpt <- cpt + 1
  }
  data_list
}
