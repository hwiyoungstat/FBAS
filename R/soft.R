#' soft
#'
#' @param x value
#' @param t threshold
#'
#' @return
#' @export

soft <- function(x,t){
  temp <- (x-t)*((x-t)>0) - (-x-t)*((-x-t)>0)
  return(temp)
}

