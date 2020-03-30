library (igraph)

g = read_graph ("~/ankit/projects/social network analysis/karate.gml", format = "gml")

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
  CPSet = FindCPSet (G, Cset)
  
  # Print the CP structures
  for (i in 1:length (Cset)) {
    cat ("\n\n Core ", i, ": ", Cset [[i]])
    for (j in 1:length (CPSet [[i]]))
      cat ("\n Level ", j, " :", CPSet [[i]] [[j]])
  }
  cat ("\n\n Overlapping nodes: ", CPSet$overlapping, "\n\n")
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
  
  
  numC = length (Cset)
  U = c ()
  ng = c ()
  active = c ()
  level = 1
  CPSet = list ()
  
  # create CPSet
  for (i in 1:numC) {
    CPSet [[i]] = list ()
  }
  
  # fill core vertices in universal CSet
  u_cset = Cset
  
  # get all the vertex in the core sets
  for (i in 1:numC)
    U = c (U, Cset [[i]])
  
  # get neighbours group
  for (i in 1:length (U)) {
    n = neighbors (lg, U [i])
    add_n = n [! n %in% union (intersect (n, U), intersect (n, ng))]
    ng = c (ng, add_n)
  }
   
  while (length (ng) != 0) {
   
    # create level
    for (j in 1:numC) {
      CPSet [[j]] [[level]] = c (-1)
    }
    
    new_ng = c ()
    u_cset_add = list ()
    for (j in 1:numC)
      u_cset_add [[j]] = c (-1)
    
    # for complete set of neighbouring vertices
    while (length (ng) != 0) {
      
      n =  neighbors (lg, ng [1])
      n = n [n %in% intersect (n, U)]
      
      count_c = c () # count connections
      for (j in 1:numC) {
        count_c = c (count_c, length (intersect (u_cset [[j]], n)))
      }
      cp = which (count_c %in% max (count_c))
      if (length (cp) > 1)
        active = c (active, ng [1])

      # fill level
      for (j in 1:length (cp)) {
        if (CPSet [[cp [j]]] [[level]] [1] == -1) {
          CPSet [[cp [j]]] [[level]] = c (CPSet [[cp [j]]] [[level]], ng [1])
          u_cset_add[[cp [j]]] = c (u_cset_add [[cp [j]]], ng [1])
          CPSet [[cp [j]]] [[level]] = CPSet [[cp [j]]] [[level]] [2:2]
        }
        else {
          CPSet [[cp [j]]] [[level]] = c (CPSet [[cp [j]]] [[level]], ng [1])
          u_cset [[cp [j]]] = c (u_cset [[cp [j]]], ng [1]) 
        }
      }
      #print (u_cset)
      
      # add new neighbours
      n = neighbors (lg, ng [1])
      new_ng = c (new_ng, n [! n %in% union (intersect (n, U), union (intersect (n, ng), intersect (n, new_ng)))])
      # add it in U
      U = c (U, ng [1])
      # remove the neighbour
      ng = ng [ng != ng [1]]
    } # for every neighbour
  
    # create new neighboouring group
    ng = new_ng
    level = level + 1
    for (j in 1:length (u_cset)) {
      u_cset_add [[j]] == u_cset_add [[j]] [u_cset_add [[j]] != -1]
      u_cset [[j]] = c (u_cset [[j]], u_cset_add [[j]])
    } 
    
  }
  
  # redistribute active vertices
  overlapping = c ()
  for (i in 1:length (active)) {
    n =  neighbors (lg, active [i])
    n = n [! n %in% intersect (n, active)]
    count_c = c ()
    remove_cp = c ()
    
    count_c = c () # count connections
    for (j in 1:numC) {
      count_c = c (count_c, length (intersect (u_cset [[j]], n)))
    }
    cp = which (count_c %in% max (count_c))
    if (length (cp) > 1)
      overlapping = c (overlapping, active [i])
    else {
      remove_cp =  c(1:length (count_c))
      remove_cp = remove_cp [! remove_cp %in% cp]
      # remove the active nodes from the given cp sets
      for (j in 1:length (remove_cp)) {
        for (k in 1:length (CPSet [[remove_cp [j]]])) {
          CPSet [[remove_cp [j]]] [[k]] = CPSet [[remove_cp [j]]] [[k]] [CPSet [[remove_cp [j]]] [[k]] != active [i]]
        }
      }
    }
  }
  
  CPSet$overlapping = overlapping
  
  return (CPSet)
}

cp_set (g, 0.9, 4)

