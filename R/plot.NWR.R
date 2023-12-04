#' plot.NWR -  A function that plots the estimated coefficient vector(s) of an NWR object. 
#' @param NWR An NWR object.
#' @param seed The random seed to set if you want the plot to change orientation every time you run the function.
#' @param boxplot_args Additional arguments to go into the boxplot() function. 
#' @param igraph_args Additional arguments to go into the igraph() funciton. 
#' @return Boxplot(s) of the distribution of the coefficient matrix and network plots for each coefficient vector. 
#' @export
plot.NWR <- function(NWR, seed = NULL,boxplot_args,igraph_args){
  
  if(!(is.null(seed))){
    set.seed(seed)
  }
  opar <- par()
  on.exit(suppressWarnings(par(opar)))
  
  if(missing(boxplot_args)){
    boxplot_args = list()
  }
  boxplot_args$x = data.frame(NWR$beta)
  boxplot_args$main = "Boxplot of Beta Estimates"
  boxplot_args$names = colnames(NWR$beta)
  do.call(boxplot,boxplot_args)
  # boxplot(data.frame(NWR$beta), main = "Boxplot of Beta Estimates",
  #         names = colnames(NWR$beta),...)
  
  readline("Hit Enter/Return for next plot")
  
  
  A1 <- diag(rowSums(NWR$A > 0)) %*% NWR$A
  
  Agraph <- igraph::graph_from_adjacency_matrix(A1, mode = "undirected")
  Agraph <- tidygraph::as_tbl_graph(Agraph)
  A_lo <- igraph::layout_nicely(Agraph)
  
  
  for(i in 1:NCOL(NWR$beta)){
    layout(matrix(c(1,1,1,1,0,
                    1,1,1,1,2,
                    1,1,1,1,2,
                    1,1,1,1,0),
                  4,5,byrow=T))
    
    par(mar = c(0,0,5,0))
    
    plot(Agraph,
         layout = A_lo,
         vertex.label = NA,
         vertex.frame = NA,
         vertex.size = 7,
         vertex.color = scales::dscale(NWR$beta[,i] %>% cut(50),viridis::viridis))
    mtext(colnames(NWR$beta)[i], side = 3, line = 2)
    par(mar = c(0,1,0,0))
    plot3D::colkey(col=scales::dscale(sort(NWR$beta[,i]) %>% cut(50),viridis::viridis),
                   clim = range(NWR$beta[,i]),
                   width = 5)
    
    readline("Hit Enter/Return for next plot")
  }
}