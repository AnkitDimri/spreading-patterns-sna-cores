library (igraph)

#####################################################################################################################
#                                               Main Function                                                       #
#####################################################################################################################

cp_set <- function (G, beta, alpha) {
  
  # Rank the nodes
  U = c()
  U$nodes = rerank_nodes (G)
  k = (2 * length (E(G))) / length (V(G))
  
  U$rd = compute_RD (U$nodes, G, alpha)
  
  x = c(1:length (V(G)))
  plot (x, y = U$rd, type = "l",xlab = "Rank",ylab = "RD")
  lines (x, y = beta + x * 0, col = "red")
  text (x, y = U$rd, labels = U$nodes, cex=0.9, font=2)
  
  Cset = FindCoreSet (U, beta, alpha)
  CPSet = FindCPSet (G, Cset)

  # Remove empty levels
  for (i in 1:length (Cset))
    for (j in 1:length (CPSet [[i]]))
      if (CPSet [[i]] [[j]] == -1 || length (CPSet [[i]] [[j]])  == 0) {
        CPSet [[i]] = CPSet [[i]] [[1:(j-1)]]
        break
      }
  
  # Make adjacency matrix representation
  mat = matrix (0, nrow = length (V (G)), ncol = length (V (G)))
  for (i in 1:length (Cset)) {
    for (j in 1:length (Cset [[i]])) {
      for (k in 1:length (Cset [[i]])) {
        if (are.connected (G, Cset [[i]] [j], Cset [[i]] [k]))
          mat [Cset [[i]] [j], Cset [[i]] [k]] = 1
      }
    }
  }
  
  for (i in 1:length (Cset)) {
    level_abv = Cset [[i]]
    for (j in 1:length (CPSet [[i]])) {
      for (k in 1:length (CPSet [[i]] [[j]])) {
        for (l in 1:length (level_abv)) {
          if (are.connected (G, CPSet [[i]] [[j]] [k], level_abv [l])) {
            mat [CPSet [[i]] [[j]] [k], level_abv [l]] = 1
            mat [level_abv [l], CPSet [[i]] [[j]] [k]] = 1
          }
        }
      }
      level_abv = CPSet [[i]] [[j]]
    }
  }
  
  heatmap (mat, Rowv = NA, Colv = NA, scale="none")
  
  # Print the CP structures
  for (i in 1:length (Cset)) {
    cat ("\n\n Core ", i, ": ", Cset [[i]])
    for (j in 1:length (CPSet [[i]]))
      cat ("\n Level ", j, " :", CPSet [[i]] [[j]])
  }
  cat ("\n\n Overlapping nodes: ", CPSet$overlapping, "\n\n")

  
  newlist <- list("cset" = Cset,"cpset" = CPSet)  
  return(newlist)
  
}


rerank_nodes <- function (lg) {
  
  # Rank the given nodes
  # Use any centrality measure for finding the first node; we used: betweeness
  U = c()
  U$n = c (V (lg) [which.max (igraph::betweenness (lg)) [[1]]])
  count = 1
  maxD = max (igraph::degree (lg))
  # Populate neighbours group
  ng = c()
  ng$p = c()
  ng$n = c(neighbors (lg, U$n [1]))
  
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
      rd = c (rd,CD (subgraph (lg, subg)))
      
    }
  }
  
  for (i in (alpha+1):length (U)) {
    subg = subg [2:alpha]
    subg = c (subg, U [i])
    rd = c (rd,CD (subgraph (lg, subg)))
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






# KARATE CLUB (alpha = 4, beta = 0.9)
g = read_graph ("datasets/karate.gml", format = "gml")
op<-cp_set (g, 0.9, 4)

get_col <- c("red","green","yellow","blue","purple","orange","cyan","pink","lightblue","lightgreen","deepskyblue4","darkseagreen3")
col <- c()
index = 1
for(i in 1:length(op$cset)){
  col[op$cset[[i]]] <- c(get_col[index])
  index = index + 1
}
for(i in 1:(length(op$cpset))){
  for(j in 1:length(op$cpset[[i]])){
    col[op$cpset[[i]][[j]]] <- c(get_col[index])
    index = index + 1
  }
  
}

regions = list()
ri = 1
for (i in 1:length (op$cset)) {
  regions [[ri]] = op$cset [[i]]
  ri = ri+1
  unionset = op$cset [[i]]
  
  for (j in 1:length (op$cpset [[i]])) {
    unionset = union (unionset, op$cpset [[i]] [[j]])
    regions [[ri]] = unionset
    ri = ri+1
  }
}

layoutt <-layout.fruchterman.reingold(g)
plot (as.directed(g, mode = c ("mutual")), mark.groups = regions, layout = layoutt, vertex.color = col, edge.arrow.size = 0, xlim = c(ceiling (min (layoutt [,1])), ceiling (max (layoutt [,1]))), ylim = c(ceiling (min (layoutt [,2])), ceiling (max (layoutt [,2]))), rescale = F, vertex.size = 50)










# DOLPHIN NETWORK (alpha = 5, beta = 0.9)
D = read.csv("datasets/dolphin.csv", header = F)
D = data.frame (D)
gd = make_empty_graph (n = 62)
# To make graph sorted and traversable for function
for (x in 1:nrow (D))
  gd = gd + edge (D [x, "V1"], D [x, "V2"])
gd = as.undirected (gd, mode = "collapse")
op = cp_set (gd, 0.9, 5)


get_col <- c("red","green","yellow","blue","purple","orange","cyan","pink","lightblue","lightgreen","deepskyblue4","darkseagreen3")
col <- c()
index = 1
for(i in 1:length(op$cset)){
  col[op$cset[[i]]] <- c(get_col[index])
  index = index + 1
}
for(i in 1:(length(op$cpset))){
  for(j in 1:length(op$cpset[[i]])){
    col[op$cpset[[i]][[j]]] <- c(get_col[index])
    index = index + 1
  }
  
}

regions = list()
ri = 1
for (i in 1:length (op$cset)) {
  regions [[ri]] = op$cset [[i]]
  ri = ri+1
  unionset = op$cset [[i]]
  
  for (j in 1:length (op$cpset [[i]])) {
    unionset = union (unionset, op$cpset [[i]] [[j]])
    regions [[ri]] = unionset
    ri = ri+1
  }
}

layoutt <-layout.fruchterman.reingold(gd)
plot (as.directed(gd, mode = c ("mutual")), mark.groups = regions, layout = layoutt, vertex.color = col, edge.arrow.size = 0, xlim = c(ceiling (min (layoutt [,1])), ceiling (max (layoutt [,1]))), ylim = c(ceiling (min (layoutt [,2])), ceiling (max (layoutt [,2]))), rescale = F, vertex.size = 50)









#Football (alpha = 10, beta = 0.8)
g = read_graph ("datasets/football.gml", format = "gml")
op<-cp_set (g, 0.6, 12)

get_col <- c("red","green","yellow","blue","purple","orange","cyan","pink","lightblue","lightgreen","deepskyblue4","darkseagreen3")
col <- c()
index = 1
for(i in 1:length(op$cset)){
  col[op$cset[[i]]] <- c(get_col[index])
  index = index + 1
}
for(i in 1:(length(op$cpset))){
  for(j in 1:length(op$cpset[[i]])){
    col[op$cpset[[i]][[j]]] <- c(get_col[index])
    index = index + 1
  }
  
}

regions = list()
ri = 1
for (i in 1:length (op$cset)) {
  regions [[ri]] = op$cset [[i]]
  ri = ri+1
  unionset = op$cset [[i]]
  
  for (j in 1:length (op$cpset [[i]])) {
    unionset = union (unionset, op$cpset [[i]] [[j]])
    regions [[ri]] = unionset
    ri = ri+1
  }
}

layoutt <-layout.fruchterman.reingold(g)
plot (as.directed(g, mode = c ("mutual")), mark.groups = regions,vertex.label=V(g), layout = layoutt, vertex.color = col, edge.arrow.size = 0, xlim = c(ceiling (min (layoutt [,1])), ceiling (max (layoutt [,1]))), ylim = c(ceiling (min (layoutt [,2])), ceiling (max (layoutt [,2]))), rescale = F, vertex.size = 50)















###################
# US Airport Network

D = read.csv("datasets/USairport500.csv", header = F, sep = "\t")
D = data.frame (D)
D = D [, 1:2]
D
gd = make_empty_graph (n = 332)
# To make graph sorted and traversable for function
for (x in 1:nrow (D))
  gd = gd + edge (D [x, "V1"], D [x, "V2"])
gd = as.undirected (gd, mode = "collapse")
op = cp_set (gd, 0.9, 5)


get_col <- c("red","green","yellow","blue","purple","orange","cyan","pink","lightblue","lightgreen","deepskyblue4","darkseagreen3")
col <- c()
index = 1
for(i in 1:length(op$cset)){
  col[op$cset[[i]]] <- c(get_col[index])
  index = index + 1
}
for(i in 1:(length(op$cpset))){
  for(j in 1:length(op$cpset[[i]])){
    col[op$cpset[[i]][[j]]] <- c(get_col[index])
    index = index + 1
  }
  
}

regions = list()
ri = 1
for (i in 1:length (op$cset)) {
  regions [[ri]] = op$cset [[i]]
  ri = ri+1
  unionset = op$cset [[i]]
  
  for (j in 1:length (op$cpset [[i]])) {
    unionset = union (unionset, op$cpset [[i]] [[j]])
    regions [[ri]] = unionset
    ri = ri+1
  }
}

layoutt <-layout.fruchterman.reingold(gd)
plot (as.directed(gd, mode = c ("mutual")), mark.groups = regions, layout = layoutt, vertex.color = col, edge.arrow.size = 0, xlim = c(ceiling (min (layoutt [,1])), ceiling (max (layoutt [,1]))), ylim = c(ceiling (min (layoutt [,2])), ceiling (max (layoutt [,2]))), rescale = F, vertex.size = 120)
