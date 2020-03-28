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
