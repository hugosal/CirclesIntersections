#' @title Simple validation test of the output of Librino_N()

#' @description Checks if the sum of the areas of intersection of each
#' circle adds to the total area of the circle, as it should.

#' @param librino A named numeric vector with the  from resulting from [Librino_N()]
#' @param radii  Numeric vectors of length N with the radius of each circle.
#'
#' @return  TRUE if all the partitions of the circles add to their total area,
#' else a numeric vector with the number of the circles that failed this test.


#' @author Hugo Salinas \email{hugosal@comunidad.unam.mx}.

#' @examples
#' # Example of intersection areas including a Reuleaux triangle
#' x <- c(0, 1, 0.5)
#' y <-c(0, 0, sqrt(1-0.5**2))
#' radii <- c(1, 1, 1)
#' intersections <- Librino_N(centers_x = x, centers_y = y, radii = radii)
#' validate_Librino(librino = unlist(intersections, use.names = TRUE), radii = radii)

#' @export

validate_Librino <- function(librino, radii){
  passed_tests <- logical(length(radii))

  for(n in seq_along(radii)){
    summ <- 0
    for (i in seq_along(librino)){
      nombres <- as.numeric(strsplit(names(librino)[i], ":")[[1]])
      if (n %in% nombres){
        summ <- summ + abs(librino[i])
      }
      difference <- (radii[n]**2*pi) - summ}
    passed_tests[n] <- abs(difference) < 1e-5 # arbitrary small number
  }
  if (all(passed_tests)){TRUE
  }else{
      which(!passed_tests)}
}

