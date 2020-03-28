library (igraph)

g = read_graph ("~/ankit/projects/social network analysis/karate.gml", format = "gml")
cp_set (g, 0.9)
#####################################################################################################################
#                                               Main Function                                                       #
#####################################################################################################################

cp_set <- function (G, beta) {
  
  # Rank the nodes
  U$nodes = rerank_nodes (G)
  k = (2 * length (E(G))) / length (V(G))
  
  for (i in 1:length (U))
    U$rd = compute_RD (U$nodes, G)
  
  Cset = FindCoreSet (G, U, beta)
  CPset = FindCPSet (G, Cset, NumC)
}


rerank_nodes <- function (lg) {
  
  # Rank the given nodes
  # Use any centrality measure for finding the first node; we used: betweeness
  U = c (which.max (igraph::betweeness (lg)))
  count = 1
  # Populate neighbours group
  ng = c(neighbors (lg, U [1]))
  
  while (count != length (V(lg))) {
    
    for (i in 1:length (ng)) {
      
    }
    
    # repopulate neighbours group
    for (i in 1:length (U)) {
      n = neighbors (lg, U [i])
      for (j in 1:length (n)) {
        if (!(n [i] %in% U)) {
          ng = c (ng, n)
        } 
      }
    }
    
    
  }
  
}