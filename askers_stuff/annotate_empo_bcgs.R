library(igraph)

empo = read.table("empo_v3.csv", sep = ";", stringsAsFactors = F, header = T)

adj_names = unique(unname(unlist(empo)))
adj = matrix(0, nrow = length(adj_names), ncol = length(adj_names), dimnames = list(adj_names,adj_names))
diag(adj) = 1

empolevel = rep(0, length(adj_names))
for (row in 1:nrow(empo)){
  for (col in 2:ncol(empo)){
    row_cur = which(adj_names == empo[row, col-1])
    col_cur = which(adj_names == empo[row, col])
    print(paste(empo[row, col-1], empo[row, col], row_cur, col_cur))
    
    empolevel[col_cur] =  col
    adj[row_cur, col_cur] = 1
  }
} 
empolevel[1] = 1 #Fix "EMP sample"
g = graph_from_adjacency_matrix(adj, mode = "undirected", diag = F)
V(g)$size = (5-empolevel)*10
plot(g)


bigscape =  read.table("visualize_bcgs/Network_Annotations_Full.tsv", sep = "\t", stringsAsFactors = F, header = T)

