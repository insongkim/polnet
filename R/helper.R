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
                            vertices = NULL) {
  if (class(edges)=="matrix") {
    edges_mat <- edges
    
    if (is.null(rownames(edges_mat))) rownames(edges_mat) <- paste0("G1_",1:nrow(edges_mat))
    if (is.null(colnames(edges_mat))) colnames(edges_mat) <- paste0("G2_",1:ncol(edges_mat))
    edges <- as.data.frame(edges_mat)
    
    edges_gathered <- tidyr::gather(edges, "group2", "count")
    edges_gathered$group1 <- rep(rownames(edges),ncol(edges_mat))
    edges_gathered <- edges_gathered[,c("group1","group2","count")]
    
    out <- igraph::graph.data.frame(edges_gathered, directed=F, vertices=vertices)
  } else {
    edges <- edges[,c(group1.id, group2.id, count.id)]
    out <- igraph::graph.data.frame(edges_gathered, directed=F, vertices=vertices)
  }
  return(out)
}
