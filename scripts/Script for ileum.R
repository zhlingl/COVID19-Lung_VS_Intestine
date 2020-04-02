countdata = function(gsmname){
   dirname = file.path(gsmdir,gsmname)
   t = Read10X(data.dir = dirname)
   t = CreateSeuratObject(counts = t)
   t[["percent.mt"]] <- PercentageFeatureSet(t, pattern = "^MT-")
   t <- subset(t, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25)

    ################Normalizing the data
	t <- NormalizeData(t)

	#############Scaling the data
	all.genes <- rownames(t)
	t <- ScaleData(t,features = all.genes)

	############Identification of highly variable features (feature selection)
	t <- FindVariableFeatures(t, selection.method = "vst", nfeatures = 2000)

	##############Perform linear dimensional reduction
	t <- RunPCA(t, features = VariableFeatures(object = t))

  t <- FindNeighbors(t, dims = 1:15)
	t <- FindClusters(t, resolution = 0.5)

	#############Run non-linear dimensional reduction (UMAP/tSNE)
	t <- RunTSNE(t,dims = 1:15)

	t.test <- as.SingleCellExperiment(t)
    t.main <- SingleR(test=t.test, ref=human, labels=human$label.main)
    #t.fine <- SingleR(test=t.test, ref=human, labels=human$label.fine)
    t[['main']]=t.main$labels
    
    celldefine = c()
	for ( i in c(1:length(levels(Idents(t))))){
	  i=i-1
	  cellsname = rownames(t[[]][t[[]][,6]==i,])
	  cellmain = t[[]][cellsname,7]
	  out= rownames(as.matrix(sort(table(cellmain),decreasing = T)))[1]
	  celldefine=rbind(celldefine,out)
	}
    
    new.cluster.ids = celldefine[,1]
    names(new.cluster.ids) = levels(t)
    t = RenameIdents(t,new.cluster.ids)
    
    dataname = paste0(gsmname,'.RData')
    dataname = file.path(dirname,dataname)
    save(t,file = dataname)

    
   
}

readtable = function(gsmname){
  dirname = file.path(gsmdir,gsmname)
  t = Read10X(data.dir = dirname)
  t = CreateSeuratObject(count = t)

  t[["percent.mt"]] <- PercentageFeatureSet(t, pattern = "^MT-")
  t <- subset(t, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
  t = NormalizeData(t)
  return(t)  
}

setwd("/data/tangzh/2019_nCov/data/ileal/result")

gsmdir = "/data/tangzh/2019_nCov/data/ileal"
gsmname = list.files(gsmdir)[grep('GSM',list.files(gsmdir))]
all.list =list()
for (i in gsmname){
	all.list[[i]] = readtable(i)
}

all.merge = merge(all.list[[1]],all.list[2:22],add.cell.ids=names(all.list) ,merge.data = T)
all.merge <- FindVariableFeatures(all.merge,nfeatures=2000)
all.merge <- ScaleData(all.merge,vars.to.regress=c('percent.mt','nCount_RNA'))
all.merge <- RunPCA(all.merge, features = VariableFeatures(object = all.merge))

pdfname = "/home/tangzh/tangzh/elbow.pdf"
pdf(file=pdfname,width=16,height=9)
print( ElbowPlot(all.merge))
dev.off()

if(F){
# all.merge <- FindNeighbors(all.merge,dims=1:16,verbose=T)
# all.merge <- FindClusters(all.merge,resolution=0.5,verbose=T)
# all.merge <- RunTSNE(all.merge,dims=1:16)
# 
# cellid = matrix(nrow = length(rownames(all.merge[[]])), ncol =1)
# rownames(cellid) = rownames(all.merge[[]])
# 
# for (i in rownames(cellid)){
#   cellid[i,1]=substr(i,1,10)
# }
# all.merge$cellid = cellid[,1]
# 
# pdfname="sample-tsne.pdf"
# pdf(file=pdfname,width=16,height=9)
# print(DimPlot(all.merge,reduction = 'tsne',group.by = 'cellid',))
# dev.off()
# 
# library(SingleR)
# load(file = '/data/tangzh/RData/whole-RNA/human2.RData')
# 
# t = as.SingleCellExperiment(all.merge)
# t.main = SingleR(test = t, ref= human, labels = human$label.main)
# all.merge[['main']] = t.main$labels
# 
# celldefine = c()
#   for ( i in c(1:length(levels(Idents(all.merge))))){
#     i=i-1
#     cellsname = rownames(all.merge[[]][all.merge[[]][,7]==i,])
#     cellmain = all.merge[[]][cellsname,8]
#     out= rownames(as.matrix(sort(table(cellmain),decreasing = T)))[1]
#     celldefine=rbind(celldefine,out)
#   }
# 
# new.cluster.ids = celldefine[,1]
# names(new.cluster.ids) = levels(all.merge)
# all.merge = RenameIdents(all.merge,new.cluster.ids)
#    
# library(ggplot2)
# plot.celltype <- DimPlot(all.merge,label=T,repel=T) + labs(title='Cell Type') + 
#   theme(plot.title = element_text(hjust = 0.5),legend.position='none')
# plot.dot <- DotPlot(all.merge,features=c('SLC6A19','CTSL','ADAM17','TMPRSS2','ACE2')) + RotatedAxis()
# 
# vln.list <- lapply(c('ACE2','TMPRSS2','ADAM17','CTSL','SLC6A19'),function(x) VlnPlot(all.merge,features = x,pt.size = 0,y.max = 4) + geom_jitter(alpha=0.2,size=0.5))
# plot.vln <- CombinePlots(plots = vln.list,ncol = 2)
# 
# pdf('/data/tangzh/2019_nCov/data/ileal/result/celltype.pdf',width = 5,height = 5)
# plot.celltype
# dev.off()
# 
# pdf('/data/tangzh/2019_nCov/data/ileal/result/dotplot.pdf',width = 6,height = 4)
# plot.dot
# dev.off()
# 
# pdf('/data/tangzh/2019_nCov/data/ileal/result/vlnplot.pdf',width = 12,height = 10)
# plot.vln
# dev.off()
}#泽华师兄的

all.merge <- RunUMAP(all.merge,dims=1:16)
all.merge <- FindNeighbors(all.merge,dims=1:16,verbose=T)
all.merge <- FindClusters(all.merge,resolution=0.8,verbose=T)
DimPlot(all.merge)

#亚型注释
library(SingleR)
hpca <- HumanPrimaryCellAtlas()
pred <- SingleR(test=all.merge@assays$RNA@data,
                ref=hpca,method='cluster',fine.tune=T,labels=hpca$label.main,cluster=all.merge@meta.data$seurat_clusters)
id.list <- pred$pruned.labels
names(id.list) <- rownames(pred)
all.merge.annot <- RenameIdents(all.merge,id.list)

id.list['17'] <- 'Mast Cell' # TPSAB1
id.list['6'] <- 'Paneth Cell' #C1QB C1QA
id.list['18'] <- 'Fibroblast' #DCN LUM
# 注释用的marker来源于
# CellMarker
# panglaodb
# wikipedia

all.merge.reannot <- RenameIdents(all.merge,id.list)

library(ggplot2)
plot.celltype <- DimPlot(all.merge.reannot,label=T,repel=T) + labs(title='Cell Type') + 
  theme(plot.title = element_text(hjust = 0.5),legend.position='none')
plot.dot <- DotPlot(all.merge.reannot,features=c('SLC6A19','CTSL','ADAM17','TMPRSS2','ACE2')) + RotatedAxis()

pdf('../plot/ileum_celltype.pdf',width = 5,height = 5)
plot.celltype
dev.off()

pdf('../plot/ileum_dotpot.pdf',width = 6,height = 4)
plot.dot
dev.off()

