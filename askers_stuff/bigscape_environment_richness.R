empo = read.table("empo_v3.csv", sep = ";", stringsAsFactors = F, header = T)
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
             .(BiG.SCAPE.class, sample_type), summarise, Count = length(sample_type))

barpagg = dcast(barp, sample_type ~ BiG.SCAPE.class, value.var = "Count")
barpagg[is.na(barpagg)] = 0
row.names(barpagg) = barpagg$sample_type
barpagg$sample_type = NULL
pheatmap::pheatmap(barpagg)

H <- diversity(barpagg)
S <- specnumber(barpagg) ## rowSums(BCI > 0) does the same...
J <- H/log(S)

