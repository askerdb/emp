library(igraph); library(reshape2); library(plyr)

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


#Generating the mapping
#(bigscape) askerbrejnrod@Askers-MacBook-Pro-2 data 2 % for i in *tar.gz;do zgrep  -a 'Original ID' ${i}|cut -f 1,3 -d ":"| awk -v i="${i}" '{print i $0}' ; done > contig_environment_mapping.txt
contigmap = read.table("contig_environment_mapping.txt", sep = " ")
contigmap = contigmap[,c("V1", "V17")]
colnames(contigmap) = c("Envfile", "contig")
contigmap$Env = sub( "\\.tar\\.gz", "", contigmap$Envfile)
bigscape =  read.table("Network_Annotations_Full.tsv", sep = "\t", stringsAsFactors = F, header = T)
bigcontig = merge(contigmap, bigscape, by.x = "contig", by.y = "Description")
# A few is missing in this merge
empo_map = read.table("emp500_metadata_shotgun_assembly_scaffolds.txt", sep = "\t", quote = "", header = T)
bigempo = merge(bigcontig, empo_map, by.x = "Env", by.y = "sample_name")

barp = ddply(bigempo[,c("BiG.SCAPE.class", "Product.Prediction", "sample_type", "empo_3", "empo_2", "empo_1")], 
      .(BiG.SCAPE.class, sample_type, empo_3), summarise, Count = length(sample_type))
ggplot(barp, aes(x = empo_3, y = Count, fill = BiG.SCAPE.class))+geom_col(stat = "identity", position="stack")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
