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
#' @export
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
#' @export
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
#' @export
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

#' Check the input data and return a list in the right list format
#' 
#' @param data An object from the package \code{funData}
#' 
#' @return A list, where each element represents a curve. Each curve is defined 
#'  as a list with two entries:
#'  \itemize{
#'   \item \strong{$t} The sampling points
#'   \item \strong{$x} The observed points
#'  }
#' @export
checkData <- function(data){
  if (inherits(data, 'funData')){
    data_ <- funData2list(data)
  } else if (inherits(data, 'irregFunData')) {
    data_ <- irregFunData2list(data)
  } else {
    stop('Wrong data type!')
  }
  data_
}


#' Generate logarithmic spaced sequence
#' 
#' @param from Numeric
#' @param to Numeric
#' @param length.out Numeric
#' 
#' @return Vector
#' @export
lseq <- function(from = 1, to = 100, length.out = 51) {
  exp(seq(log(from), log(to), length.out = length.out))
}


#' Compute different kernels
#' 
#' @param u Vector of numeric, points to estimate the kernel.
#' @param type Integer, used kernel. If 1, uniform kernel, If 2, Epanechnikov
#'  kernel. If 3, biweight kernel.
#'  
#' @return Vector of numeric
#' @export
kernel <- function(u, type = 1){
  indicator <- function(u) 2 * stats::dunif(u, -1, 1)
  switch(type,
         indicator(u) / 2,
         0.75 * (1 - u**2) * indicator(u),
         0.9375 * (1 - u**2)**2 * indicator(u))
}
