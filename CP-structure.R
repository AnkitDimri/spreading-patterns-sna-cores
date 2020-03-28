library (igraph)

g = read_graph ("~/ankit/projects/social network analysis/karate.gml", format = "gml")
cp_set (g, 0.9)
#####################################################################################################################
#                                               Main Function                                                       #
#####################################################################################################################

cp_set <- function (G, beta) {
  
  # Rank the nodes
  U = c()
  U$nodes = rerank_nodes (G)
  #k = (2 * length (E(G))) / length (V(G))
  print (U$nodes)
  #for (i in 1:length (U))
   # U$rd = compute_RD (U$nodes, G)
  
  #Cset = FindCoreSet (G, U, beta)
  #CPset = FindCPSet (G, Cset, NumC)
}


rerank_nodes <- function (lg) {
  
  # Rank the given nodes
  # Use any centrality measure for finding the first node; we used: betweeness
  U = c()
  U$n = c (which.max (igraph::betweenness (lg)))
  count = 1
  maxD = max (igraph::degree (lg))
  # Populate neighbours group
  ng = c()
  ng$p = c()
  ng$n = c(neighbors (lg, U [1]))
  
  while (count != length (V(lg))) {
    
    ng$p = c()
    for (i in 1:length (ng$n))
      ng$p = c (ng$p, length (intersect (neighbors (lg, ng$n [i]), U$n)) + (igraph::degree (lg, ng$n [i]) / maxD))
    
    add_node = ng$n [which.max (ng$p)]
    # add to list
    U$n = c (U$n, add_node)
    # repopulate neighbours group
    # remove from neighbours group
    ng$n = ng$n [ng$n != add_node]
    # add its neighbours to the group
    n = neighbors (lg, add_node)
    add_n = n [! n %in% intersect (n, U$n)]
    ng$n = c (ng$n, add_n)
    count = count + 1
    
  }
  
  return (U$n)
}

