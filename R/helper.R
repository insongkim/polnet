#' Convert edges to igraph object
#' @param edges Matrix or data.frame or igraph of connection strengths as counts
#' @param group1.id A character string indicating the name of group1 identifier
#' variable in the \code{edges} data.frame. It is required in the case of data.frame.
#' @param group2.id A character string indicating the name of group2 identifier
#' variable in the \code{edges} data.frame. It is required in the case of data.frame.
#' @param count.id A character string indicating the name of count identifier
#' variable in the \code{edges} data.frame. The variable must be numeric.
#' @param vertices A data frame with vertex metadata, or NULL.

#' @import igraph
#' @import tidyr
#' @return igraph object
#' @useDynLib polnet, .registration = TRUE
#' @export

edges_to_igraph <- function(edges, 
                            group1.id = NULL,
                            group2.id = NULL,
                            count.id = NULL,
                            group1.cluster = NULL,
                            group2.cluster = NULL) {
  ## Warning for missing parameter
  if (class(edges)=="data.frame"&is.null(group1.id))
    stop("'group1.id' should be provided")
  if (class(edges)=="data.frame"&is.null(group2.id))
    stop("'group2.id' should be provided")
  if (class(edges)%in%c("data.frame", "igraph")&is.null(count.id))
    stop("'count.id' should be provided")
  
  if (class(edges)=="data.frame") {
    edges <- edges[,c(group1.id, group2.id, count.id)]
    edges <- tidyr::spread(edges, group2.id, count.id)
    rownames(edges) <- edges[,group1.id]
    edges <- edges[,-1]
    edges <- as.matrix(edges)
  }
  g <- igraph::graph_from_incidence_matrix(edges, directed=F)

  if (!is.null(group1.cluster)) {
    V(g)$cluster[V(g)$type==F] <- group1.cluster
  }
  if (!is.null(group2.cluster)) {
    V(g)$cluster[V(g)$type==T] <- group2.cluster
  }
    
  return(g)
}
