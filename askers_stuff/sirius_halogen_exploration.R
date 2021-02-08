library(vegan); library(cowplot); library(ape); library(ztable)
siriusmfidabund = read.table("1907_EMPv2_INN_GNPS_quant.csv", sep = ",", header = T)
#Explore MF identifications
siriurmf = read.table("data/formula_identifications.tsv", sep = "\t", header = T)
siriurmf = subset(siriurmf, ZodiacScore > 0.9)
siriurmf$shared_id = unname(sapply(siriurmf$id, function(x) strsplit(x, "_")[[1]][5]) )

siriurmf$Halogen = grepl("Br|I|F|Cl", siriurmf$molecularFormula)
siriurmf$HalogenType = unlist(sapply(siriurmf$molecularFormula, function(x) c("Br", "I", "Cl", "F", "Not halogen", "Multiple halogen substitutions")[c(grepl("Br", x), grepl("I", x), grepl("Cl", x), grepl("F", x), !grepl("Br|I|F|Cl", x), grep("Br|I|F|Cl", x) > 1)][1]))
siriurmf_halo = siriurmf[grep("Br|I|F|Cl", siriurmf$molecularFormula),]
ggplot(siriurmf, aes(x = retentionTimeInSeconds, y = ionMass, color = HalogenType)) + 
  facet_wrap(~HalogenType) + geom_point() + theme_bw()
ggplot(siriurmf, aes(  ionMass, fill = HalogenType)) + geom_density(alpha = 0.5) + theme_bw()
ztable(data.frame(table(siriurmf$HalogenType)))
siriurmfjoint = merge(siriurmf, siriusmfidabund, by.x = "shared_id", by.y = "row.ID", sort = F )
siriurmfjoint$IntensitySum = rowSums(siriurmfjoint[,33:834])
ggplot(siriurmfjoint, aes(log10(IntensitySum), fill =  HalogenType)) + geom_density(alpha = 0.7) + facet_grid(HalogenType ~ .) + theme_bw()

#For all halogens
siriurmfjointhal = subset(siriurmfjoint, Halogen == T)

plot_halogen_richness = function(halogen_joint_table, siriusmeta, title = ""){
  siriurmfjointhal = halogen_joint_table
siriurmfjointhalonly = siriurmfjointhal[,33:833]
colnames(siriurmfjointhalonly) = sub("X", "", gsub("\\.", "_", sub(".mzML.cropped.Peak.area","", colnames(siriurmfjointhalonly))))
siriusmeta = read.table("qiime2_metadata.tsv", sep = "\t", header =T, comment.char = "", quote = "")

siriusmeta$new_filename = gsub("\\.", "_", sub(".mzML", "", gsub("-", "_", siriusmeta$metabo_v2_filename_2019)))
siriushalojoint = merge(t(siriurmfjointhalonly), siriusmeta, by.x = 0, by.y = "new_filename")
row.names(siriushalojoint) = siriushalojoint$Row.names; siriushalojoint$Row.names = NULL
siriushalojoint01 =  siriushalojoint[, 1:(ncol(siriushalojoint) - ncol(siriusmeta))]
siriushalojoint01[siriushalojoint01 > 0] = 1
siriushalojoint$HalogenRichness = rowSums(siriushalojoint01)
ggplot(siriushalojoint, aes(x = env_biome, y = HalogenRichness, color = empo_3))  + theme_bw() + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + ggtitle(title)
}

p1 = plot_halogen_richness(subset(siriurmfjoint, Halogen == T), siriusmeta, title = "All halogen compounds")
pbr = plot_halogen_richness(subset(siriurmfjoint, HalogenType == "Br"), siriusmeta, title = "Br compounds")
pi = plot_halogen_richness(subset(siriurmfjoint, HalogenType == "I"), siriusmeta, title = "I compounds")
pf = plot_halogen_richness(subset(siriurmfjoint, HalogenType == "F"), siriusmeta, title = "F compounds")
pcl = plot_halogen_richness(subset(siriurmfjoint, HalogenType == "Cl"), siriusmeta, title = "Cl compounds")
library(cowplot)

plot_grid(p1, pbr, pi, pf, pcl)


#Plot PCoA 

siriurmfjointhal = subset(siriurmfjoint, Halogen == T)

siriurmfjointhalonly = siriurmfjointhal[,33:833]
colnames(siriurmfjointhalonly) = sub("X", "", gsub("\\.", "_", sub(".mzML.cropped.Peak.area","", colnames(siriurmfjointhalonly))))
siriusmeta = read.table("qiime2_metadata.tsv", sep = "\t", header =T, comment.char = "", quote = "")
  
siriusmeta$new_filename = gsub("\\.", "_", sub(".mzML", "", gsub("-", "_", siriusmeta$metabo_v2_filename_2019)))
siriushalojoint = merge(t(siriurmfjointhalonly), siriusmeta, by.x = 0, by.y = "new_filename")
row.names(siriushalojoint) = siriushalojoint$Row.names; siriushalojoint$Row.names = NULL
siriushalojointabund =  siriushalojoint[, 1:(ncol(siriushalojoint) - ncol(siriusmeta)+1)]
siriurmfjointhaldist = vegdist(siriushalojointabund)
siriurmfjointhaldistfit = pcoa(siriurmfjointhaldist)
plot(siriurmfjointhaldistfit$values$Relative_eig[1:10]*100, ylab = "Percent variance explained", xlab = "Principal coordinate", main = "Scree plot")
siriurmfjointhaldistfitjoint = data.frame(siriurmfjointhaldistfit$vectors[,1:2], siriushalojoint[, ncol(siriurmfjointhalonly):ncol(siriushalojoint)])
ggplot(siriurmfjointhaldistfitjoint, aes(x = Axis.1, y = Axis.2, color = empo_2)) + geom_point() + theme_bw()

ggplot(siriurmfjointhaldistfitjoint, aes(x = Axis.1, y = Axis.2, color = env_biome)) + geom_point() + theme_bw() + facet_wrap(~empo_3)

#Project halogen type on pcoa

siriurmfhalotype = merge(t(siriushalojointabund), siriurmfjointhal,  by = 0, sort = F, all.x = T)

siriurmfjointhalagg = aggregate(t(siriushalojointabund), by = list(siriurmfhalotype$HalogenType), FUN = sum )
row.names(siriurmfjointhalagg) = siriurmfjointhalagg$Group.1;siriurmfjointhalagg$Group.1 = NULL
siriurmfjointhalaggn = sweep(siriurmfjointhalagg, 2,  colSums(siriurmfjointhalagg),  FUN = '/')
x = siriurmfjointhaldistfit
Y = t(siriurmfjointhalaggn)
plot.axes = c(1, 2)
pr.coo <- x$vectors
n <- nrow(Y)
points.stand <- scale(pr.coo[, plot.axes])
S <- cov(Y, points.stand)
U <- S %*% diag((x$values$Eigenvalues[plot.axes]/(n - 
                                                    1))^(-0.5))
colnames(U) <- colnames(pr.coo[, plot.axes])
U = data.frame(U)
U$Halogen = row.names(U)

hull <- siriurmfjointhaldistfitjoint %>%
  group_by(lcms_batch) %>%
  slice(chull(Axis.1, Axis.2))

ggplot(siriurmfjointhaldistfitjoint, aes(x = Axis.1, y = Axis.2, color = empo_2)) + geom_point() + theme_bw() + 
  geom_segment(data = data.frame(U), aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), arrow = arrow(length = unit(0.5, "cm")), inherit.aes = F) + 
  geom_label(data = data.frame(U), aes(x = Axis.1, y = Axis.2, label = Halogen), inherit.aes = F) + stat_ellipse(aes(group = lcms_samplingmethod, linetype = lcms_samplingmethod))
  #geom_polygon(data = hull, alpha = 0.5)


ggplot(siriurmfjointhaldistfitjoint, aes(x = Axis.1, y = Axis.2, color = empo_2)) + geom_point() + theme_bw() + 
  geom_segment(data = data.frame(U), aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), arrow = arrow(length = unit(0.2, "cm")), color = "grey", linetype = "dashed", inherit.aes = F) + 
  geom_text(data = data.frame(U), aes(x = Axis.1, y = Axis.2, label = Halogen),color = "grey", inherit.aes = F) + facet_wrap(~emp500_title)


#Explore identified compounds
siriusmfid = read.table("data/compound_identifications.tsv", sep = "\t", header = T, comment.char = "")
siriusmfid$shared_id = unname(sapply(siriusmfid$id, function(x) strsplit(x, "_")[[1]][5]) )
siriusmfid$Halogen = grepl("Br|I|F|Cl", siriusmfid$molecularFormula)
siriusmfid$HalogenType = unlist(sapply(siriusmfid$molecularFormula, function(x) c("Br", "I", "Cl", "F", "Not halogen", "Multiple halogen substitutions")[c(grepl("Br", x), grepl("I", x), grepl("Cl", x), grepl("F", x), !grepl("Br|I|F|Cl", x), grep("Br|I|F|Cl", x) > 1)][1]))
siriusmfid = siriusmfid#subset(siriusmfid, Halogen == T)
ggplot(siriurmf, aes(x = retentionTimeInSeconds, y = ionMass, color = HalogenType)) + 
  facet_wrap(~HalogenType) + geom_point() + theme_bw()
siriusmfidjoint = merge(siriusmfid, siriusmfidabund, by.x = "shared_id", by.y = "row.ID", sort = F )
siriusmfidjoint$IntensitySum = rowSums(siriusmfidjoint[,33:837])
#siriusmfidjoint$IntensityPercentile =  sapply(siriusmfidjoint$IntensitySum, function(x) sum( x >siriusmfidjoint$IntensitySum)/ *100) #This is too complicated

#table of halogen type counts
ggplot(siriusmfidjoint, aes(log10(IntensitySum), fill =  HalogenType)) + geom_density(alpha = 0.7) + facet_grid(HalogenType ~ .) + theme_bw()
library(ztable);options(ztable.type="viewer")
ztable(data.frame(table(siriusmfidjoint$HalogenType)))

#Top abundant compounds of each type
library(tidyr); library(dplyr)
siriusmfidjointtop = subset(siriusmfidjoint, name != "null" & HalogenType != "Not halogen") %>% group_by(HalogenType) %>% top_n(n=10)
siriusmfidjointtoppp = siriusmfidjointtop[, c("HalogenType", "name", "molecularFormula", "IntensitySum")][order(siriusmfidjointtop$HalogenType,siriusmfidjointtop$IntensitySum),]
ztable(data.frame(siriusmfidjointtoppp))
