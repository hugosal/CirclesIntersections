#' @title Basic plot of circles

#' @description Generates a very simple diagram of circles, useful
#' for testing the functions of this package. The text in the plot indicates the
#' area of each circle. The number of each circle is written on top of itself.

#' @param centers_x,centers_y,radii  Numeric vectors of length N with the
#' coordinates of the center of the circles in x and y, and the radius
#' respectively.

#' @param npoints The number of points to plot each circle with. Defaults to 500

#' @author Hugo Salinas \email{hugosal@comunidad.unam.mx}.

#' @examples

#' # Example of intersection areas including a Reuleaux triangle
#' x <- c(0, 1, 0.5)
#' y <-c(0, 0, sqrt(1-0.5**2))
#' radii <- c(1, 1, 1)
#' plot_circles_simple(centers_x = x, centers_y = y, radii = radii)
#'
#' # Example with more circles
#' x2 <- c(0, 4, 2, 4, 5)
#' y2 <- c(1, 5, 4, 2, 1)
#' radii2 <- c(1, 4 ,2, 2, 1)
#' plot_circles_simple(centers_x = x2, centers_y = y2, radii = radii2)


#' @importFrom grDevices col2rgb rainbow rgb
#' @importFrom graphics grid polygon text

#' @export
#'
plot_circles_simple <- function(centers_x, centers_y, radii, npoints=500){

  plot(1, xlab="x", ylab="y", type="n",
       xlim = c(min(centers_x) - max(radii), max(centers_x) + max(radii)*2),
       ylim = c(min(centers_y) - max(radii), max(centers_y) + max(radii)*2))
  grid()
  colors <- col2rgb(rainbow(length(centers_y)), alpha = 0.2)

  for (i in 1:length(centers_y)){
    xes <- seq(from = centers_x[i] - radii[i],
               to = centers_x[i] + radii[i],
               length.out = npoints )

    yes1 <- centers_y[i] + sqrt(radii[i]**2 - (xes-centers_x[i])**2)
    yes2 <- centers_y[i] - sqrt(radii[i]**2 - (xes-centers_x[i])**2)

    xes <- c(xes[!is.nan(yes1)], rev(xes[!is.nan(yes2)]))
    yes <- c(yes1[!is.nan(yes1)], rev(yes2[!is.nan(yes2)]))
    polygon( xes, yes,
             col=  rgb(colors[1, i], colors[2, i],
                       colors[3, i], 255/2, maxColorValue=255))
    text(centers_x[i], centers_y[i] + radii[i], i)
    text(x =  max(centers_x) + max(radii ) , y= max(centers_y)+max(radii)*2-i+1,
         bquote("C"~ .(i)~  A == .(round(radii[i]**2*pi ,3)) ))

  }
}
