library(Seurat)
library(SeuratData)
library(dplyr)
library(ggplot2)
library(mclust)
library(pROC)
library(caret)

#Load CBMC data
cbmc=LoadData(ds='cbmc')

#Nomarlization and feature selection
DefaultAssay(cbmc) <- 'RNA'
cbmc <- NormalizeData(cbmc) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(cbmc) <- 'ADT'
VariableFeatures(cbmc) <- rownames(cbmc[["ADT"]])
cbmc <- NormalizeData(cbmc, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca', verbose=FALSE)

#OMIC method
rna.emb = t(cbmc@assays$RNA@scale.data)
adt.emb = t(cbmc@assays$ADT@scale.data)

model=lm(adt.emb~rna.emb)

adt.reduce=model$residuals
rna.emb = cbmc@reductions$pca@cell.embeddings[,1:30]

colnames(adt.reduce)=c(1:10)
colnames=colnames(adt.reduce)
colnames=paste("apca_",colnames,sep="")
colnames(adt.reduce)=colnames


#Clustering
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

#Find cell marker for each resulted cluster
DefaultAssay(cbmc.integrate) <- 'RNA'
cbmc.markers <- FindAllMarkers(cbmc.integrate, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.75)
cbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

marker_index=list()

for(i in 1:14){
  marker_index[[i]]=which(cbmc.markers$cluster==i-1)
}


ADT=t(cbmc.integrate$ADT@scale.data)
colnames(ADT)=paste('adt_',colnames(ADT))

#Divide into training set and testing set  
set.seed(123)
separate <- sample(c(rep(0,6032), rep(1,2585)))
train=ADT[separate==0,]
test=ADT[separate==1,]

cluster=rep(0,8617)
cluster_list=c()
index_list=c()
model_list=c()
roc_train=c()
roc_test=c()

#Perform logistic regression on each cluster using only ADT information  
for(i in 1:14){
    index_list[[i]]=which(cbmc.integrate$seurat_clusters==(i-1))
    cluster[index_list[[i]]]=1
    cluster_list[[i]]=cluster
    cluster=rep(0,8617)
    
    
    cluster_train=cluster_list[[i]][separate==0]
    cluster_test=cluster_list[[i]][separate==1]
    
    train=cbind(train,cluster_train)
    train=as.data.frame(train)
    test=cbind(test,cluster_test)
    test=as.data.frame(test)
    
    temp <- (0.5*train$cluster_train*(1/sum(train$cluster_train== 1) - 1/sum(train$cluster_train == 0)) + 0.5/sum(train$cluster_train == 0))
    dat_weights = temp*length(train$cluster_train)
    
    g1=glm(cluster_train~.,data=train, weights=dat_weights, family = 'binomial')
    
    model_list[[i]]=g1
    
    predict=ifelse(g1$fitted.values>0.5,1,0)
    c1=confusionMatrix(data=factor(predict, levels = c(1,0)), reference = factor(cluster_train, levels = c(1,0)))

    predict2=predict(g1,newdata=test,'response')
    predict3=ifelse(predict2>0.5,1,0)
    c2=confusionMatrix(data=factor(predict3, levels = c(1,0)), reference = factor(cluster_test, levels = c(1,0)))
    
    roc_train=c(roc_train, roc(cluster_train, predict)$auc[1])
    roc_test=c(roc_test, roc(cluster_test, predict3)$auc[1])

    
    train=ADT[separate==0,]
    test=ADT[separate==1,]
}
  
cluster=rep(0,8617)
cluster_list=c()
index_list=c()
model_list2=c()
roc_train2=c()
roc_test2=c()

#Perform logistic regression on each cluster using only RNA information  
for(i in 1:14){
    
    pbmc <- CreateSeuratObject(counts = cbmc.integrate@assays$RNA@counts[cbmc.markers$gene[marker_index[[i]]],])
    pbmc <- NormalizeData(pbmc)%>% ScaleData()
    
    rna.emb=t(pbmc@assays$RNA@scale.data)
    separate <- sample(c(rep(0,6032), rep(1,2585)))
    train=rna.emb[separate==0,]
    test=rna.emb[separate==1,]
    
    index_list[[i]]=which(cbmc.integrate$seurat_clusters==(i-1))
    cluster[index_list[[i]]]=1
    cluster_list[[i]]=cluster
    cluster=rep(0,8617)
    
    
    cluster_train=cluster_list[[i]][separate==0]
    cluster_test=cluster_list[[i]][separate==1]
    
    train=cbind(train,cluster_train)
    train=as.data.frame(train)
    test=cbind(test,cluster_test)
    test=as.data.frame(test)
    
    temp <- (0.5*train$cluster_train*(1/sum(train$cluster_train== 1) - 1/sum(train$cluster_train == 0)) + 0.5/sum(train$cluster_train == 0))
    dat_weights = temp*length(train$cluster_train)
    
    g1=glm(cluster_train~.,data=train, weights=dat_weights, family = 'binomial')
    
    model_list2[[i]]=g1
    
    
    
    predict=ifelse(g1$fitted.values>0.5,1,0)
    c1=confusionMatrix(data=factor(predict, levels = c(1,0)), reference = factor(cluster_train, levels = c(1,0)))
    
    
    predict2=predict(g1,newdata=test,'response')
    predict3=ifelse(predict2>0.5,1,0)
    c2=confusionMatrix(data=factor(predict3, levels = c(1,0)), reference = factor(cluster_test, levels = c(1,0)))
    
    roc_train2=c(roc_train2, roc(cluster_train, predict)$auc[1])
    roc_test2=c(roc_test2, roc(cluster_test, predict3)$auc[1])

    
    train=rna.emb[separate==0,]
    test=rna.emb[separate==1,]
    
  }
  
  cluster=rep(0,8617)
  cluster_list=c()
  index_list=c()
  model_list3=c()
  roc_train3=c()
  roc_test3=c()

  #Perform logistic regression on each cluster using only RNA information  
  for(i in 1:14){
    pbmc <- CreateSeuratObject(counts = cbmc.integrate@assays$RNA@counts[cbmc.markers$gene[marker_index[[i]]],])
    pbmc <- NormalizeData(pbmc) %>% ScaleData()
    
    rna.emb0=t(pbmc@assays$RNA@scale.data)
    rna.emb = cbind(rna.emb0,ADT)
    separate <- sample(c(rep(0,6032), rep(1,2585)))
    train=rna.emb[separate==0,]
    test=rna.emb[separate==1,]
    
    
    
    
    index_list[[i]]=which(cbmc.integrate$seurat_clusters==(i-1))
    cluster[index_list[[i]]]=1
    cluster_list[[i]]=cluster
    cluster=rep(0,8617)
    
    
    cluster_train=cluster_list[[i]][separate==0]
    cluster_test=cluster_list[[i]][separate==1]
    
    train=cbind(train,cluster_train)
    train=as.data.frame(train)
    test=cbind(test,cluster_test)
    test=as.data.frame(test)
    
    temp <- (0.5*train$cluster_train*(1/sum(train$cluster_train== 1) - 1/sum(train$cluster_train == 0)) + 0.5/sum(train$cluster_train == 0))
    dat_weights = temp*length(train$cluster_train)
    
    g1=glm(cluster_train~.,data=train, weights=dat_weights, family = 'binomial')
    
    model_list3[[i]]=g1
    
    predict=ifelse(g1$fitted.values>0.5,1,0)
    c1=confusionMatrix(data=factor(predict, levels = c(1,0)), reference = factor(cluster_train, levels = c(1,0)))
    
    
    
    predict2=predict(g1,newdata=test,'response')
    predict3=ifelse(predict2>0.5,1,0)
    c2=confusionMatrix(data=factor(predict3, levels = c(1,0)), reference = factor(cluster_test, levels = c(1,0)))
    
    roc_train3=c(roc_train3, roc(cluster_train, predict)$auc[1])
    roc_test3=c(roc_test3, roc(cluster_test, predict3)$auc[1])
    
    train=rna.emb[separate==0,]
    test=rna.emb[separate==1,]
  }


#Find coefficients and p-value for CD8 T cluster    
model_cd8t_both=model_list3[[6]]
summary=summary(model_cd8t_both)
pvalue=summary$coefficients[,4]
pvalue[-1]
  
coefficients=data.frame(summary$coefficients[,1])
predictors=rownames(coefficients)[-1]
coefficients=as.data.frame(coefficients[-1,])
coefficients <- as.data.frame(cbind(predictors, coefficients))
coefficients <- coefficients %>% mutate(sign=factor(ifelse(.$`coefficients[-1, ]`>0,'positive','negative'), levels=c('positive','negative')))
coefficients=data.frame(cbind(coefficients, pvalue[-1]))
coefficients <- coefficients %>% mutate(abs.value=ifelse(.$pvalue..1.<0.01, abs(.$coefficients..1...), NA))
coefficients <- subset(coefficients, coefficients$pvalue..1.<0.05)
coefficients <- cbind(coefficients, 'CD8 T')
coefficients <- subset(coefficients,pvalue..1.<0.05)
colnames(coefficients)[6]='Cluster'
coefficients1 <- coefficients
  

#Find coefficients and p-value for CD16+ Mono cluster  
model_cd16pmono_both=model_list3[[7]]
summary=summary(model_cd16pmono_both)
pvalue=summary$coefficients[,4]
pvalue[-1]
  
coefficients=data.frame(summary$coefficients[,1])
predictors=rownames(coefficients)[-1]
coefficients=as.data.frame(coefficients[-1,])
coefficients <- as.data.frame(cbind(predictors, coefficients))
coefficients <- coefficients %>% mutate(sign=factor(ifelse(.$`coefficients[-1, ]`>0,'positive','negative'), levels=c('positive','negative')))
coefficients=data.frame(cbind(coefficients, pvalue[-1]))
coefficients <- coefficients %>% mutate(abs.value=ifelse(.$pvalue..1.<0.05, abs(.$coefficients..1...), NA))
coefficients <- cbind(coefficients, 'CD16+ Mono')
coefficients <- subset(coefficients,pvalue..1.<0.05)
colnames(coefficients)[6]='Cluster'
coefficients2 <- coefficients
  
  
#Find coefficients and p-value for CD14+ Mono cluster   
model_cd14pmono_both=model_list3[[1]]
summary=summary(model_cd14pmono_both)
pvalue=summary$coefficients[,4]
pvalue[-1]
  
coefficients=data.frame(summary$coefficients[,1])
predictors=rownames(coefficients)[-1]
coefficients=as.data.frame(coefficients[-1,])
coefficients <- as.data.frame(cbind(predictors, coefficients))
coefficients <- coefficients %>% mutate(sign=factor(ifelse(.$`coefficients[-1, ]`>0,'positive','negative'), levels=c('positive','negative')))
coefficients=data.frame(cbind(coefficients, pvalue[-1]))
coefficients <- coefficients %>% mutate(abs.value=ifelse(.$pvalue..1.<0.05, abs(.$coefficients..1...), NA))
coefficients <- cbind(coefficients, 'CD14+ Mono')
coefficients <- subset(coefficients,pvalue..1.<0.05)
colnames(coefficients)[6]='Cluster'
coefficients3 <- coefficients
  
coefficients=as.data.frame(rbind(coefficients1, coefficients2, coefficients3))
  
#Plot significant coefficients values  
ggplot(coefficients, aes(x=Cluster,y=predictors, colour=sign, size=abs.value)) + geom_point()+labs(x=NULL)
  
