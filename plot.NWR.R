plot.NWR <- function(NWR_dat, seed = NULL,boxplot_args,igraph_args){
  
  if(!(is.null(seed))){
    set.seed(seed)
  }
  opar <- par()
  on.exit(suppressWarnings(par(opar)))
  
  if(missing(boxplot_args)){
    boxplot_args = list()
  }
  boxplot_args$x = data.frame(NWR_dat$beta)
  boxplot_args$main = "Boxplot of Beta Estimates"
  boxplot_args$names = colnames(NWR_dat$beta)
  do.call(boxplot,boxplot_args)
  # boxplot(data.frame(NWR_dat$beta), main = "Boxplot of Beta Estimates",
  #         names = colnames(NWR_dat$beta),...)
  
  readline("Hit Enter/Return for next plot")
  
  
  A1 <- diag(rowSums(NWR_dat$A > 0)) %*% NWR_dat$A
  
  Agraph <- igraph::graph_from_adjacency_matrix(A1, mode = "undirected")
  Agraph <- tidygraph::as_tbl_graph(Agraph)
  A_lo <- igraph::layout_nicely(Agraph)
  
  
  for(i in 1:NCOL(NWR_dat$beta)){
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
         vertex.color = scales::dscale(NWR_dat$beta[,i] %>% cut(50),viridis::viridis))
    mtext(colnames(NWR_dat$beta)[i], side = 3, line = 2)
    par(mar = c(0,1,0,0))
    plot3D::colkey(col=scales::dscale(sort(NWR_dat$beta[,i]) %>% cut(50),viridis::viridis),
                   clim = range(NWR_dat$beta[,i]),
                   width = 5)
    
    readline("Hit Enter/Return for next plot")
  }
}