Regen_Nemve1_normalized_cpm.txt

rm(list = ls())
gc()
library(WGCNA)

regen_WGNCA_data <- read.table(file='data/Regen_Nemve1_normalized_cpm.txt', header = T, quote = "") %>%
  mutate_all(function(x) log(x+1,2)) %>%
  t()

#Need to pick soft powers. This takes a while
powers = c(c(1:10), seq(from = 12, to=40, by=2))
disableWGCNAThreads()
sft = pickSoftThreshold(regen_WGNCA_data, powerVector = powers, verbose = 5, networkType = "signed")
sft$powerEstimate

{cex1 = 0.9;
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.80,col="red")
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}

softPower = 18 # (Read WGCNA tutorial to learn how to pick your power)
rankExpr = rank(colMeans(regen_WGNCA_data))
rankConn = rank(softConnectivity((regen_WGNCA_data),type="signed",power=softPower)) 

enableWGCNAThreads()
library(flashClust)

adjacency = adjacency((regen_WGNCA_data),power=softPower,type="signed"); 
diag(adjacency)=0
Degree <- rowSums(adjacency)


Degree=Degree/max(Degree); gc()
dissTOM = 1-TOMsimilarity(adjacency, TOMType="signed") 
geneTree = flashClust(as.dist(dissTOM), method="average")

mColor=NULL
for (ds in 0:3){
  tree = cutreeHybrid(dendro = geneTree, pamStage=F, minClusterSize = (30-3*ds), deepSplit = ds, distM = dissTOM)
  mColor=cbind(mColor,labels2colors(tree$labels)); }


plotDendroAndColors(geneTree, mColor, paste("dpSplt =",0:3), main = "",dendroLabels=FALSE); 

modules = mColor[,1] # (Chosen based on plot below)
dynamicColors= labels2colors(modules)

#print the module sizes
table(modules)

PCs = moduleEigengenes((regen_WGNCA_data), colors=modules)
ME  = PCs$eigengenes
MEDiss= 1-cor(ME)
METree= flashClust(as.dist(MEDiss), method= "average")
MEs= orderMEs(ME)

plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

#optional merge
MEDissThres= 0.15
merge= mergeCloseModules(fischer_WGNCA_data, modules, cutHeight= MEDissThres, verbose =3)
mergedColors= merge$colors

library(reshape)
ME_plot <- as.data.frame(t(MEs))
colnames(ME_plot) = c(-1,-1,-1,0,0,0,2,2,2,4,4,4,8,8,8,12,12,12,16,16,16,20,20,20,24,24,24,36,36,36,48,48,48,60,60,60,72,72,72,96,96,96,120,120,120,144,144,144)

library(ggplot2)
row.names(ME_plot) <- rownames(ME_plot)
ME_plotM <- melt(as.matrix(ME_plot), id.vars="row.names", value.name="value")

co = c("MEred")

ME_plot_color <- ME_plotM[ME_plotM$X1 == paste(co),]
group <- c(1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2)

p <-  ggplot(ME_plot_color, aes(x=as.numeric(as.character(ME_plot_color$X2)), y=value, colour=paste(co), group = group)) + 
  geom_line(size =1.5) +
  scale_x_continuous(minor_breaks = NULL, breaks=c(0,2,4,8,12,16,20,24,36,48,60,72,96,120,144)) +
  ylab("Eigengene expression") +
  xlab("Hours Post Fertilization")

print(p)

# for(co in rownames(ME_plot)){
#   ME_plot_color <- ME_plotM[ME_plotM$X1 == paste(co),]
#   group <- c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3)
#   
#   p <-  ggplot(ME_plot_color, aes(x=as.numeric(as.character(ME_plot_color$X2)), y=value, colour=paste(co), group = group)) + 
#     geom_line(size =1.5) +
#     scale_x_continuous(minor_breaks = NULL, breaks=c(0,2,4,8,12,16,20,24,36,48,60,72,96,120,144)) +
#     ylab("Eigengene expression") +
#     xlab("Hours Post Amputation")
#   
#   file = paste("plots/WGCNA/ME_lines_R/ME_R_", co, ".pdf", sep="")
#   pdf(file = file, width = 8, height = 6)
#   print(p)
#   dev.off()
# }


traits <- read.table(file="data/regen_traits.txt", sep='\t', header=T, quote="",row.names=1)

nGenes = ncol(regen_WGNCA_data)
nSamples = nrow(regen_WGNCA_data)
moduleTraitCor = cor(ME, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
library('viridis')
library(RColorBrewer)
colors = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(200)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = c("UC","0hpf","2hpf","4hpf","8hpf","12hpf","16hpf","20hpf","24hpf","36hpf","48hpf","60hpf","72hpf","96hpf","120hpf","144hpf"),
               yLabels = names(ME),
               colorLabels = FALSE,
               colors = colors,
               textMatrix = textMatrix,
               setStdMargins = F,
               cex.text = 0.5,
               zlim = c(-1,1),
               verticalSeparator.x = 1:dim(moduleTraitCor)[2],
               verticalSeparator.col="grey80",
               horizontalSeparator.y = 1:dim(moduleTraitCor)[1],
               horizontalSeparator.col="grey80",
               main = paste("Module-hpf relationships")
)

annotations = read.table(url("https://raw.githubusercontent.com/ScientistJake/Spiralbase_VLM_tutorial/main/data/Nemve1_annotated_uprotnames.txt"), sep='\t', header=T, quote="")
regen_data <- read.table(file='data/Regen_Nemve1_normalized_cpm.txt', header = T, quote = "")

master = merge(regen_data,annotations, by.x='row.names',by.y='ContigName', all.x=T)
master$module = modules
write.table(master,file="data/master_regen_annotations.txt", sep='\t',row.names = F,quote=F)

# go for visant
TOM = TOMsimilarity(adjacency, TOMType="signed") 

module = "cyan";
probes = colnames(regen_WGNCA_data)
inModule = (modules==module);
modProbes = probes[inModule];

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,file = paste("data/regen_VisANTInput-", module, "_0.05.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0.05,
                            #below is helper function where you can replace the gene ID with a common name. Just pass it two vectors, one of the IDs and one of the common names
                            probeToGene = data.frame(annotations$ContigName, annotations$Description)
)
write.table(vis, file= paste("data/regen_VisANTInput-", module, "_0.05.txt", sep=""),sep='\t',row.names = F,quote=F)

#now export just the top 250 connections
nTop = 250;
IMConn = softConnectivity(regen_WGNCA_data[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],file = paste("data/regen_VisANTInput-", module, "-top",nTop,".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annotations$ContigName, annotations$Description)
)
write.table(vis, file=paste("data/regen_VisANTInput-", module, "-top",nTop,".txt", sep=""),sep='\t',row.names = F,quote=F)
