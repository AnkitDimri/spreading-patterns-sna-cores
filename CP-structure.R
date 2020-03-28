library (igraph)

g = read_graph ("~/ankit/projects/social network analysis/karate.gml", format = "gml")
cp_set (g, 0.9, 4)
#####################################################################################################################
#                                               Main Function                                                       #
#####################################################################################################################

cp_set <- function (G, beta, alpha) {
  
  # Rank the nodes
  U = c()
  U$nodes = rerank_nodes (G)
  k = (2 * length (E(G))) / length (V(G))
  
  U$rd = compute_RD (U$nodes, G, alpha)
  
  plot (x = c(1:length (V(G))) , y = U$rd, type = "l")
  
  Cset = FindCoreSet (U, beta, alpha)
  CPset = FindCPSet (G, Cset)
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


compute_RD <- function (U, lg, alpha) {
  
  rd = c (0)
  subg = c (U [1])
  
  for (i in 2:alpha) {
    if (i <= length (U)) {
      subg = c (subg, U [i])
      rd = c (rd, CD (subgraph (lg, subg)))
    }
  }

  for (i in (alpha+1):length (U)) {
    subg = subg [2:alpha]
    subg = c (subg, U [i])
    rd = c (rd, CD (subgraph (lg, subg)))
  }

  return (rd)
}

CD <- function (lg) {
  return ((2* length (E(lg))) / (length (V(lg)) * (length (V(lg)) - 1)))
}


FindCoreSet <- function (U, beta, alpha) {
  
  CoreSet = list ()
  numC = 1
  core = c ()
  
  for (i in alpha:length (U$nodes)) {
    
    if (U$rd [i] >= beta)
      core = union (core, c (U$nodes [(i-alpha+1):i]))
    else
      if (U$rd [i - 1] >= beta && i > alpha) {
        CoreSet [[numC]] = core
        numC = numC + 1
        core = c ()
      }
  }
  
  if (U$rd [length (U$nodes)] < beta)
    numC = numC - 1
  
  return (CoreSet)
}


FindCPSet <- function (lg, Cset) {
  
}