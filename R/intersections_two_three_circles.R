#' @title Auxiliary functions for computing circle intersection areas

#' @description
#'  `intersection_two_circles()` Returns the area of intersection of two circles.
#'
#'  `intersection_three_circles()` Returns the area of intersection of three circles.
#'
#' @inheritParams Librino_N
#'
#' @export
intersection_two_circles <- function(centers_x, centers_y, radii){
  r1 <- radii[1]
  r2 <- radii[2]
  D <- sqrt( (centers_x[1] - centers_x[2])**2 + (centers_y[1] - centers_y[2])**2 )

  if (D <= abs(r1 - r2)){ # one circle contains the other
    return (pi * min(r1, r2)**2)

  }else if (D > r1 + r2){ # no intersection
    return (0)
  }

  aa <- (r1**2 - r2**2 + D**2) / (2 * D)
  bb <- (r2**2 - r1**2 + D**2) / (2 * D)
  th1 <- 2 * acos(aa / r1)
  th2 <- 2 * acos(bb / r2)

  if (th1 > pi){
    a1 <- (r2**2 * (th2 + sin((2 * pi) - th2)))/2
    } else {
      a1 <- (r2**2 * (th2 - sin(th2)))/2}

  if (th2 > pi){
    a2 <- (r1**2 * (th1 + sin((2 * pi) - th1)))/2
    } else {
      a2 <- (r1**2 * (th1 - sin(th1)))/2
    }
  a1 + a2
  }

#' @rdname intersection_two_circles
intersection_three_circles <- function(centers_x, centers_y, radii){

  # circles must be sorted according to their radius
  order_rad <- order(radii, decreasing = TRUE)
  cx <- centers_x[order_rad]
  cy <- centers_y[order_rad]
  r <- radii[order_rad]
  Dist <- as.matrix(dist(matrix(c(cx, cy), ncol = 2)))

  if (any(combn(1:3, 2,
                FUN = function(x){Dist[x[1], x[2] ] > r[x[1]] + r[x[2]] } ))){return(0)}

  # count which circles are inside of another circle
  circles_inside <-  combn(1:3, 2,
                     FUN = function(x){
                       r[x[1]] >=  Dist[x[1], x[2]] + r[x[2]]} )

  # count which circles intersect
  circles_intersecting <- combn(1:3, 2,
                 FUN = function(x){
                   if ((r[x[1]] - r[x[2]]) < Dist[x[1], x[2]] &
                       Dist[x[1], x[2]] < (r[x[1]] + r[x[2]])){
                     x
                     }else{
                       c(-1,-1)}
                   })

  counts_intersect <- circles_intersecting[, circles_intersecting[1,] != -1]

  if (any(circles_inside)){

    container_name <- c(1, 1, 2)
    contained_name <- c(2, 3, 3)

    maximum_container <- min(container_name[circles_inside])
    circle_most_inside <- max(contained_name[circles_inside])
    other_circle <- (1:3)[-c(maximum_container, circle_most_inside)]
    return(intersection_two_circles(centers_x = cx[c(circle_most_inside, other_circle)],
                                    centers_y = cy[c(circle_most_inside, other_circle)],
                                    radii =  r[c(circle_most_inside, other_circle)]))

    }

  intersec_points_contained <- c(FALSE, FALSE , FALSE)
  # compute points of overlap between  a pair of circles, and see if
  # the third other circle (g) contains both of them
  for (g in 1:3){
    c1 <- (1:3)[-g][1]
    c2 <- (1:3)[-g][2]
    D_c1_c2 <- Dist[c1, c2]

    intsrs_points <- two_circles_inters_points(centers_x = cx[c(c1, c2)],
                              centers_y =  cy[c(c1, c2)],
                              radii = r[c(c1, c2)],
                              distance = D_c1_c2)

    dist_matrix <- as.matrix(dist(matrix(c(cx[g], intsrs_points[, 1],
                                           cy[g], intsrs_points[, 2]),
                                         ncol = 2  )))
    if (all(dist_matrix[1, 2:3] < r[g])){
      intersec_points_contained[g] <- TRUE
      }
    }

  if (sum(intersec_points_contained) == 2){

    containing_circle <- which(!intersec_points_contained)
    containing_circle_area <- r[containing_circle]**2 * pi

    area_remove <- 0

    for (f in (1:3)[-containing_circle]){
      area_inters <- intersection_two_circles(cx[c(f, containing_circle )],
                               cy[c(f, containing_circle )],
                                r[c(f, containing_circle )])

      area_remove <- area_remove + (containing_circle_area - area_inters) }

    return(containing_circle_area - area_remove)

  }else if (sum(intersec_points_contained) == 1){

    circs <- which(!intersec_points_contained)

    return(intersection_two_circles(centers_x = cx[circs],
                                    centers_y = cy[circs],
                                    radii =  r[circs]))
  }

    if( !(r[1] - r[2]) < Dist[1, 2] &  Dist[1, 2]  < (r[1] + r[2])){return(0)}

  rs_1 <- r[1]**2
  rs_2 <- r[2]**2
  rs_3 <- r[3]**2

  d1_2 <- Dist[1, 2]
  d2_3 <- Dist[2, 3]
  d1_3 <- Dist[1, 3]

  x12 <- (rs_1 - rs_2 + d1_2**2) / (2 * d1_2)
  y12 <- sqrt(((2 * d1_2**2) * (rs_1 + rs_2)) - (rs_1 - rs_2)**2 - d1_2**4 )/(2 * d1_2)

  costheta <- (d1_2**2 + d1_3**2 - d2_3**2)/(2 * d1_2* d1_3)
  sintheta <- sqrt(1 - costheta**2)

  costhetap <- -(d1_2**2 + d2_3**2 - d1_3**2)/(2 * d1_2 * d2_3)
  sinthetap <- sqrt(1 - costhetap**2)

  if (! (x12 - (d1_3 * costheta))**2 + (y12 - (d1_3 * sintheta))**2 < rs_3 ){return(0)}

  if (! (x12 - (d1_3 * costheta))**2 + (y12 + (d1_3 * sintheta))**2 > rs_3 ){return(0)}

  x13p <- (rs_1 - rs_3 + d1_3**2) / (2 * d1_3)
  y13p <- -sqrt((2 * d1_3**2 * (rs_1 + rs_3)) - (rs_1 - rs_3)**2 - d1_3**4)/(2 * d1_3)

  x13 <- (x13p * costheta) - (y13p * sintheta)
  y13 <- (x13p * sintheta) + (y13p * costheta)

  x23pp <- (rs_2 - rs_3 + d2_3**2)/(2 * d2_3)
  y23pp <- sqrt((2 * d2_3**2 * (rs_2 + rs_3)) - (rs_2 - rs_3)**2 - d2_3**4)/(2 * d2_3)

  x23 <- (x23pp * costhetap) - (y23pp * sinthetap) + d1_2
  y23 <- (x23pp * sinthetap) + (y23pp * costhetap)

  c1 <- sqrt((x13 - x12)**2 + (y13 - y12)**2) # i3 j2 k1
  c2 <- sqrt((x12 - x23)**2 + (y12 - y23)**2) # i1 j3 k2
  c3 <- sqrt((x23 - x13)**2 + (y23 - y13)**2) # i2 j1 k3

  A <- (sqrt((c1 + c2 + c3) * (c2 + c3 - c1) * (c1 + c3 - c2) * (c1 + c2 - c3))/4) +
    ((rs_1 * asin(c1 / (2 * r[1]))) - ( (c1/ 4) * sqrt((4 * rs_1) - c1**2))) +
    ((rs_2 * asin(c2 / (2 * r[2]))) - ( (c2/ 4) * sqrt((4 * rs_2) - c2**2)))

  if ((d1_3 * sintheta) < (y13 + ((y23 - y13) / (x23 - x13)) * ((d1_3 * costheta) - x13) )){
      arc_3 <- (r[3]**2 * pi) - ((rs_3 * asin(c3 / (2 * r[3])))  - ( ((c3/ 4) * sqrt((4 * rs_3) - c3**2) )))
      }else{
    arc_3 <-((rs_3 * asin(c3 / (2 * r[3])))  - ( ((c3/ 4) * sqrt((4 * rs_3) - c3**2) ) ))
    }
  A + arc_3
  }


two_circles_inters_points <- function(centers_x, centers_y, radii,
                                      distance,
                                      circle_numbers = NULL){

  a <- (radii[1]**2 - radii[2]**2 + distance**2) / (2 * distance)
  h <- sqrt(radii[1]**2 - a**2)
  x2 <- centers_x[1] + a * (centers_x[2] - centers_x[1])/distance
  y2 <- centers_y[1] + a * (centers_y[2] - centers_y[1])/distance

  dx <- h*(centers_y[2]-centers_y[1])/distance
  dy <- h*(centers_x[2]-centers_x[1])/distance

  x3_1 <- x2+dx
  x3_2 <- x2-dx

  y3_1 <- y2-dy
  y3_2 <- y2+dy
  if (! is.null(circle_numbers)){
    names1 <- paste(c(circle_numbers, "A"), collapse =  ":" )
    names2 <- paste(c(circle_numbers, "B"), collapse =  ":" )
  }else{
    names1 <- NULL
    names2 <- NULL}
  matrix(c(x3_1, x3_2, y3_1, y3_2), nrow = 2,
         dimnames = list(c(names1, names2),
                         c(NULL)))
}
