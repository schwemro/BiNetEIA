#' @title Double-checking the amount of vertices
#'
#' @description Double-checks the amount of vertices in the network if they are equal to the amount in the interaction matrix. Otherwise the missing vertices will be added to the graph.
#' @param g Network as igraph object
#' @param no_v number of vertices

check_vertices <- function(g, no_v){
  if (length(V(g)) < no_v){
    x <- no_v - length(V(g))
    g <- add_vertices(g,x)
  }
  return(g)
}
