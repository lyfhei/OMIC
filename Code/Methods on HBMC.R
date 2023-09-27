library(MOFA2)
library(Seurat)
library(SeuratData)
library(dplyr)
library(ggplot2)
library(mclust)
library(pROC)
library(caret)

# Load HBMC dataset
bm=LoadData('bmcite')

# Normalization, feature selection
DefaultAssay(bm) <- 'RNA'
bm <- NormalizeData(bm) %>% FindVariableFeatures(nfeatures=1000) %>% ScaleData() %>% RunPCA()

DefaultAssay(bm) <- 'ADT'
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>%
  ScaleData() %>% RunPCA(reduction.name = 'apca')

# WNN method
bm <- FindMultiModalNeighbors(
  bm, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)

bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
bm <- FindClusters(bm, graph.name = "wsnn", algorithm = 3, resolution = 1, verbose = FALSE)

p1 <- DimPlot(bm, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) +ggtitle(NULL)+NoLegend()+xlab("UMAP1(WNN)")+ylab("UMAP2(WNN)") 
p2 <- DimPlot(bm, reduction = 'wnn.umap', group.by = 'celltype.l2', label = TRUE, repel = TRUE, label.size = 4.5) +ggtitle(NULL)+NoLegend()+xlab("UMAP1")+ylab("UMAP2") 
p1 + p2

adjustedRandIndex(bm$celltype.l2,bm$seurat_clusters)


# OMIC method
n.adt.pc = 25

rna.emb = t(bm@assays$RNA@scale.data)
adt.emb = t(bm@assays$ADT@scale.data)

model=lm(adt.emb~rna.emb)

adt.reduce=model$residuals
rna.emb = bm@reductions$pca@cell.embeddings[,1:30]

colnames(adt.reduce)=c(1:25)
colnames=colnames(adt.reduce)
colnames=paste("apca_",colnames,sep="")
colnames(adt.reduce)=colnames

#Visualization
integrate.data = t(scale(data.frame(cbind(rna.emb,adt.reduce))))
bm.integrate = bm
bm.integrate@assays$integrate=CreateAssayObject(data = integrate.data)

DefaultAssay(bm.integrate) <- 'integrate'
bm.integrate@assays$integrate@key = "integrate_"
VariableFeatures(bm.integrate) <- rownames(bm.integrate[["integrate"]])

bm.integrate@reductions$integrateRed = CreateDimReducObject(embeddings = t(as.matrix(integrate.data)),
                                                              assay = "integrate",key = 'integrateRed_')

bm.integrate <- RunUMAP(bm.integrate, reduction = 'integrateRed',dims = 1:50,
                          reduction.name = "integrate.umap", reduction.key = "integrateUMAP_")
bm.integrate <- FindNeighbors(bm.integrate, reduction = 'integrateRed',dims = 1:50
                                , assay = 'integrate',graph.name = "intg_graph",k.param =20)

bm.integrate <- FindClusters(bm.integrate , graph.name = "intg_graph"
                               ,resolution = 1.1,method=3)

p4 <- DimPlot(bm.integrate, reduction = 'integrate.umap', label = TRUE,
              repel = TRUE, label.size = 5)+ NoLegend() +ggtitle(NULL)+xlab("UMAP1(OMIC)")+ylab("UMAP2(OMIC)") 
p5 <- DimPlot(bm.integrate, reduction = 'integrate.umap', group.by = 'celltype.l2', label = TRUE,
              repel = TRUE, label.size = 4.5)+ NoLegend() +ggtitle(NULL)+xlab("UMAP1")+ylab("UMAP2") 

p4+p5


# MOFA+ method
bm=LoadData('bmcite')

DefaultAssay(bm) <- 'RNA'
bm <- NormalizeData(bm) %>% FindVariableFeatures(nfeatures=2000) %>% ScaleData() %>% RunPCA()

DefaultAssay(bm) <- 'ADT'
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>%
  ScaleData() %>% RunPCA(reduction.name = 'apca')


mofa <- create_mofa(bm, assays = c("RNA","ADT"))
mofa
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 15  
mofa <- prepare_mofa(mofa,
                     model_options = model_opts
)
mofa <- run_mofa(mofa)  

factors=get_factors(mofa)

# Visualization
bm.integrate = bm
bm.integrate@assays$integrate=CreateAssayObject(data = t(factors$group1))

DefaultAssay(bm.integrate) <- 'integrate'
bm.integrate@assays$integrate@key = "integrate_"
VariableFeatures(bm.integrate) <- rownames(bm.integrate[["integrate"]])

bm.integrate@reductions$integrateRed = CreateDimReducObject(embeddings = t(as.matrix(t(factors$group1))),
                                                              assay = "integrate",key = 'integrateRed_')

bm.integrate <- RunUMAP(bm.integrate, reduction = 'integrateRed',dims = c(1:length(factors$group1[1,])),
                          reduction.name = "integrate.umap", reduction.key = "integrateUMAP_")
bm.integrate <- FindNeighbors(bm.integrate, reduction = 'integrateRed',dims = c(1:length(factors$group1[1,]))
                                , assay = 'integrate',graph.name = "intg_graph")

bm.integrate <- FindClusters(bm.integrate , graph.name = "intg_graph"
                               ,resolution = 0.23,method=3)

adjustedRandIndex(bm.integrate$intg_graph_res.0.23,bm.integrate$celltype.l2)