#' @title Implementation of Librino's algorithm for computing circle intersection areas

#' @description This function computes the exclusive areas of intersection of N
#' circles.

#' @details This is an implementation of Librino, Levorato, and Zorzi (2014)
#' algorithm for computation of the intersection areas of an arbitrary
#' number of circles. Given a set of circles this function computes all
#' the possible intersection areas of the possible subsets of the circles.

#' @param centers_x,centers_y,radii  Numeric vectors of length N with the
#' coordinates of the center of the circles in x and y, and the radius
#' respectively.

#' @return A list of length N containing the areas of exclusive intersection.
#' The position of each element in the list indicates the number of intersecting
#' circles. The first element of the list corresponds to the
#' area of non-overlap of every circle, the second element is the pairwise
#' area of intersection. Up to the last element of the list which corresponds to
#' the area of intersection of all circles.
#' Each of the elements of the list is a named numeric vector corresponding to the
#' area of intersection between a set of circles. The names of the vector
#' indicate the number of the circles in the intersection.

#' @author Hugo Salinas \email{hugosal@comunidad.unam.mx}.

#' @references

#' Librino, F., Levorato, M., & Zorzi, M. (2014). An algorithmic solution for
#' computing circle intersection areas and its applications to wireless
#' communications. Wireless Communications and Mobile Computing, 14, 1672–1690.

#' @examples

#' # Example of intersection areas including a Reuleaux triangle
#' x <- c(0, 1, 0.5)
#' y <- c(0, 0, sqrt(1-0.5**2))
#' radii <- c(1, 1, 1)
#' intersections <- Librino_N(centers_x = x, centers_y = y, radii = radii)

#' # Example with more circles
#' x2 <- c(0, 4, 2, 4, 5)
#' y2 <- c(1, 5, 4, 2, 1)
#' radii2 <- c(1, 4 ,2, 2, 1)
#' intersections2 <- Librino_N(centers_x = x2, centers_y = y2, radii = radii2)

#' @importFrom stats dist
#' @importFrom utils combn

#' @export

Librino_N <- function(centers_x, centers_y, radii){

  N <- length(centers_x)

  if (N != length(centers_y) | N != length(radii) ){
    stop("arguments must have the same length")
  }

  if (any(radii <= 0)){
    stop("radii must be greater than 0")
  }

  # This function returns the decimal position of non zero element in a binary code
  get_circle_number_from_binary <- function(binary){
    which(strsplit(binary, "")[[1]]==1)
  }

  # This function return the transition matrix from n to n + k
  get_transition_M_n_to_k <- function(n, k, transitions){
    name_transition <- paste(n + k - 1, ":", n + k, sep = "")
    a_bar <- transitions[[name_transition]]
    if (k > 1){
      for (i in 1:(k - 1)){
        name_transition <- paste(n + k - i - 1, ":", n + k - i, sep = "")
        a_bar <- a_bar %*% transitions[[name_transition]]}
    }
    a_bar/factorial(k)
  }

  # This function returns a subtrellis given the nelement terminal node
  # in the nvert subset
  get_sub_trellis <- function(original, nvert, nelement, start){
    final_node <- names(original[[nvert]][nelement])
    reduced_trelis <- original[1:(nvert - 1)]
    n_final_node <- strsplit(final_node, "")[[1]]
    for (back in (nvert - 1):start){
      included_nodes <- sapply(names(reduced_trelis[[back]]), function(x){
        sum(!(strsplit(x, "")[[1]] == n_final_node)) == (nvert - back)
      })
      reduced_trelis[[back]] <- reduced_trelis[[back]][included_nodes]
    }
    reduced_trelis
  }

  # This function returns the transition matrix of a corresponding
  # list of areas of  overlap a_r
  transition_from_a_i <- function(a, start){
    transition_matrices_a <- list()
    # 1 is 0:1
    # 2 is 1:2 ...
    for (r in start:(length(a)- 1 )){
      next_lab <- paste(r, ":", r + 1, sep = "")
      names_rows <- lapply(names(a[[r + 1]]), function(x) strsplit(x, split = "")[[1]])
      names_cols <- lapply(names(a[[r]]),  function(x) strsplit(x, split = "")[[1]])

      this_matrix <- outer(X = 1:length(names_rows), Y =  1:length(names_cols),
                  FUN = Vectorize(function(x, y){
                    sum(!(names_rows[[x]] ==  names_cols[[y]])) == 1}))

      rownames(this_matrix) <- names(a[[r+1]])
      colnames(this_matrix) <- names(a[[r]])

      transition_matrices_a[[next_lab]] <- this_matrix
    }
    transition_matrices_a
  }

  areas <- pi * radii**2

  if (N==2){
    intersection <- intersection_two_circles(centers_x, centers_y, radii)
    return (list("1:2" = intersection,
                 "1" = areas[1] - intersection,
                 "2" = areas[2] - intersection ))
  }else{

    d <- as.matrix(stats::dist(matrix(c(centers_x, centers_y), ncol = 2)))
    vector_dec_from_bin <- numeric(2**N)

    names(vector_dec_from_bin) <- sapply(X = 1:length(vector_dec_from_bin),
                                         FUN = function(x){
                              bin <- paste(rev(as.integer(intToBits(x))), collapse="")
                              substr(bin, start = nchar(bin) - N + 1 , stop = nchar(bin))
                                         } )

    sum_binary <- sapply(names(vector_dec_from_bin),
                         FUN = function(x){sum(as.numeric(strsplit(x, "")[[1]]))})

    a_i <- lapply(1:N, function(x) sort(which(sum_binary == x), decreasing = TRUE) )

    # # Transition matrices of the complete trellis

    transition_matrices <- transition_from_a_i(a_i, start = 1)

    # Areas of a_1 are the areas of each circle

    a_i[[1]] <- sapply(names(a_i[[1]]), function(x){
      areas[get_circle_number_from_binary(x)]})

    # Areas of a_2 are the areas of each pairs of circles intersections

    a_i[[2]] <- sapply(names(a_i[[2]]), function(x){
      this_pair <- get_circle_number_from_binary(x)
      intersection_two_circles(centers_x = centers_x[this_pair],
                               centers_y = centers_y[this_pair],
                               radii = radii[this_pair])    })

    # Areas of a_3 are the areas of each threecircle
    a_i[[3]] <- sapply(names(a_i[[3]]), function(x){
      this_three <- get_circle_number_from_binary(x)
      intersection_three_circles(centers_x = centers_x[this_three],
                                 centers_y = centers_y[this_three],
                                 radii = radii[this_three])
    })

    # Areas of intersection of a_n n>= 4 are computed using a_bar vectors

    if (N >= 4){
      for (u in 4:length(a_i)){
        for (v in seq_along(a_i[[u]])){
          another_one <- u != 4 # this is becasue if  u=4, a_1 may be necessary
          minimum_depth <- ifelse (u == 4, 1, u - 2) # 1 - another_one
          Abar <- vector(mode = "list", length = u - 1 - another_one)
          subtrelis <- get_sub_trellis(nvert = u, nelement =  v, original = a_i, start = minimum_depth )
          transicion_sub <- transition_from_a_i(subtrelis, start = minimum_depth)

          # dont compute if there is no intersection
          if(any(subtrelis[[length(subtrelis)]]==0)){
            a_i[[u]][v] <-0
            next}

          for (e in minimum_depth:(u - 1 - another_one )){
            product_sum <- 0
            if (e < u - 1){
              for (j in (e + 1):(u - 1)){
                tnk <- get_transition_M_n_to_k(n = e,
                                         k = j - e,
                                         transitions = transicion_sub)
                this_rep <- (-1)**(j - e + 1) * (t(tnk) %*% as.matrix(subtrelis[[j]]))
                product_sum <- product_sum + this_rep
              }
            }
            Abar[[e]] <- subtrelis[[e]] - product_sum
          }

          if (u == 4){
            agam <- min(Abar[[1]])
            bgam <- max(-Abar[[2]])
            cgam <- min(Abar[[3]])

            if (any(subtrelis[[3]] <1e-6) ){
              a_i[[u]][v] <-  min(Abar[[u-1]])  # m = 0

            }else{
              if (cgam > bgam & # c > b
                  abs(cgam - cgam) < 1e-5){ # a = c

                four_circles <- get_circle_number_from_binary(names(a_i[[u]][v]))
                combinations <- utils::combn(four_circles, m = 2)

                # get all circles intersection points
                intersc_pts <- lapply(1:ncol(combinations), FUN = function(x){
                  c1 <- combinations[1, x]
                  c2 <- combinations[2, x]
                  D_c1_c2 <- d[c1, c2]

                  two_circles_inters_points(centers_x = centers_x[c(c1, c2)],
                                            centers_y = centers_y[c(c1, c2)],
                                            radii = radii[c(c1, c2)],
                                            distance = D_c1_c2,
                                            circle_numbers = c(c1, c2))
                })

                intersections_mat <- do.call(rbind, intersc_pts)

                # Get which intersection points are inside each circle,
                # not considering the intersections with that circle

                inters_pts_inside_circle <- lapply(four_circles, function(x){
                  distances_to_inters <- as.matrix(stats::dist(rbind(
                    matrix(c(centers_x[x], centers_y[x]), nrow = 1),
                    intersections_mat)))[, 1]

                  not_including <-sapply(rownames(intersections_mat), function(y){
                    ! as.character(x) %in% strsplit(y, split = ":")[[1]]
                  })
                  unname(which(distances_to_inters[-1] < radii[x] & not_including))
                })

                the_one <- logical(4)

                # Find if the intersection of the points
                # inside one of the circles with the other circles is
                # just one element

                for (j in seq_along(inters_pts_inside_circle)){
                  the_one[j] <- all(sapply(seq_along(inters_pts_inside_circle)[-j],
                                           function(x){
                                             length(intersect(inters_pts_inside_circle[[j]],
                                                              inters_pts_inside_circle[[x]]))==1
                                           }))
                }

                if (all(unlist(lapply(inters_pts_inside_circle,
                                      function(x) length(x) == 3)))){
                  if (sum(the_one)==1){

                    a_i[[u]][v] <-   min(Abar[[u-1]])
                  }else{

                    a_i[[u]][v] <- max(-Abar[[u-2]])
                  }
                }else{
                  a_i[[u]][v] <-  max(-Abar[[u-2]])
                }

              }else{
                a_i[[u]][v] <-  max(-Abar[[u-2]])
              }}
          }else{
            a_i[[u]][v] <- max(-Abar[[u-2]])
          }}
      }
    }

    # Compute exclusive intersection areas using a_i vectors

    Intersections_final <- list()

    for (i in 1:N){
      This_intersection_area <- a_i[[i]]
      product_sum <- 0
      if (i < N){
        product_sum <- product_sum + apply( sapply((i+1):N, function(j){

          tnk <- get_transition_M_n_to_k(n = i,
                                   k = j - i,
                                   transitions = transition_matrices)

          ((-1)**(j - i + 1)) * (t(tnk) %*% as.matrix(a_i[[j]]))}), MARGIN = 1, sum)
      }
      Intersections_final[[i]] <- as.matrix( This_intersection_area - product_sum )
    }

    for (l in length(Intersections_final):1 ){
      this_vect <- Intersections_final[[l]]
      rownames(this_vect)<- sapply(rownames(this_vect), function(x){
        paste(get_circle_number_from_binary(x), collapse = ":")})
      if (length(this_vect)==0){
        Intersections_final[[l]] <- NULL
      }else{
        Intersections_final[[l]] <- this_vect[, 1]
      }
    }
    Intersections_final
  }
}
