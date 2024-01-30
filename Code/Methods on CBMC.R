library(MOFA2)
library(Seurat)
library(SeuratData)
library(dplyr)
library(ggplot2)
library(mclust)
library(pROC)
library(caret)


# Load CBMC data
cbmc=LoadData(ds='cbmc')

# Normalization, feature seleciton
DefaultAssay(cbmc) <- 'RNA'
cbmc <- NormalizeData(cbmc) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(cbmc) <- 'ADT'
VariableFeatures(cbmc) <- rownames(cbmc[["ADT"]])
cbmc <- NormalizeData(cbmc, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca', verbose=FALSE)

# WNN method
cbmc <- FindMultiModalNeighbors(
  cbmc, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:9), modality.weight.name = "RNA.weight"
)

cbmc <- RunUMAP(cbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
cbmc <- FindClusters(cbmc, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)

p1 <- DimPlot(cbmc, reduction = 'wnn.umap', group.by = 'rna_annotations', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(cbmc, reduction = 'wnn.umap', group.by = 'protein_annotations', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2

# OMC method
rna.emb = t(cbmc@assays$RNA@scale.data)
adt.emb = t(cbmc@assays$ADT@scale.data)

model=lm(adt.emb~rna.emb)

adt.reduce=model$residuals
rna.emb = cbmc@reductions$pca@cell.embeddings[,1:30]

n.rna.pc = ncol(cbmc[['pca']])-20
n.adt.pc = ncol(cbmc[['apca']])

colnames(adt.reduce)=c(1:10)
colnames=colnames(adt.reduce)
colnames=paste("apca_",colnames,sep="")
colnames(adt.reduce)=colnames

# Visualization
integrate.data = t(scale(data.frame(cbind(rna.emb,adt.reduce))))

cbmc.integrate = cbmc
cbmc.integrate@assays$integrate=CreateAssayObject(data = integrate.data)

DefaultAssay(cbmc.integrate) <- 'integrate'
cbmc.integrate@assays$integrate@key = "integrate_"
VariableFeatures(cbmc.integrate) <- rownames(cbmc.integrate[["integrate"]])

cbmc.integrate@reductions$integrateRed = CreateDimReducObject(embeddings = t(as.matrix(integrate.data)),
                                                              assay = "integrate",key = 'integrateRed_')

cbmc.integrate <- RunUMAP(cbmc.integrate, reduction = 'integrateRed',dims = 1:(n.adt.pc+n.rna.pc),
                          reduction.name = "integrate.umap", reduction.key = "integrateUMAP_")
cbmc.integrate <- FindNeighbors(cbmc.integrate, reduction = 'integrateRed',dims = 1:(n.adt.pc+n.rna.pc)
                                , assay = 'integrate',graph.name = "intg_graph")

cbmc.integrate <- FindClusters(cbmc.integrate , graph.name = "intg_graph"
                               ,resolution = 0.878,method=3)

DimPlot(cbmc.integrate, reduction = 'integrate.umap', group.by = 'rna_annotations',repel = TRUE,label = TRUE, label.size = 4.5)+ggtitle(NULL)+xlab("UMAP1")+ylab("UMAP2")

adjustedRandIndex(cbmc.integrate$seurat_clusters,cbmc$rna_annotations)

# MOFA+ method
mofa <- create_mofa(cbmc, assays = c("RNA","ADT"))
mofa
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 15  
mofa <- prepare_mofa(mofa,
                     model_options = model_opts
)
mofa <- run_mofa(mofa)  

factors=get_factors(mofa)

# Visualization
cbmc.integrate = cbmc
cbmc.integrate@assays$integrate=CreateAssayObject(data = t(factors$group1))

DefaultAssay(cbmc.integrate) <- 'integrate'
cbmc.integrate@assays$integrate@key = "integrate_"
VariableFeatures(cbmc.integrate) <- rownames(cbmc.integrate[["integrate"]])

cbmc.integrate@reductions$integrateRed = CreateDimReducObject(embeddings = t(as.matrix(t(factors$group1))),
                                                              assay = "integrate",key = 'integrateRed_')

cbmc.integrate <- RunUMAP(cbmc.integrate, reduction = 'integrateRed',dims = c(1:length(factors$group1[1,])),
                          reduction.name = "integrate.umap", reduction.key = "integrateUMAP_")
cbmc.integrate <- FindNeighbors(cbmc.integrate, reduction = 'integrateRed',dims = c(1:length(factors$group1[1,]))
                                , assay = 'integrate',graph.name = "intg_graph")

cbmc.integrate <- FindClusters(cbmc.integrate , graph.name = "intg_graph"
                               ,resolution = 0.45,method=3)

p4 <- DimPlot(cbmc.integrate, reduction = 'integrate.umap', label = TRUE,
              repel = TRUE, label.size = 5) + NoLegend() +ggtitle(NULL)
p5 <- DimPlot(cbmc.integrate, reduction = 'integrate.umap', group.by = 'rna_annotations', label = TRUE,
              repel = TRUE, label.size = 5) + NoLegend() +ggtitle(NULL)

p4 <- DimPlot(cbmc.integrate, reduction = 'integrate.umap', label = TRUE,
              repel = TRUE, label.size = 5)+ NoLegend() +ggtitle(NULL)+xlab("UMAP1(OMIC)")+ylab("UMAP2(OMIC)") 
p5 <- DimPlot(cbmc.integrate, reduction = 'integrate.umap', group.by = 'rna_annotations', label = TRUE,
              repel = TRUE, label.size = 4.5)+ggtitle(NULL)+xlab("UMAP1")+ylab("UMAP2") 

p4+p5
