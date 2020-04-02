if(F){
  require(GEOquery)
  small.intestine.list <- as.vector(unlist(read.table('./small_intestine.txt')))
  data.list <- sapply(small.intestine.list,function(x) getGEOSuppFiles(x,baseDir='.'))
}

require(ggplot2)
require(Seurat)
require(patchwork)

setwd('project/chenwx/scRNA_analysis/data/data_small_intestine')
intestine_list <- as.vector(unlist(read.table('../small_intestine.txt',stringsAsFactors=F)))
data.list <- lapply(intestine_list,function(x) read.table(paste(x,list.files(x),sep='/'),
                                                      stringsAsFactors=F,header=T,row.names=1))
names(data.list) <- intestine_list

obj.list <- lapply(1:3,function(x) CreateSeuratObject(data.list[[x]],project=names(data.list[x]),min.cells=3))
obj.list <- lapply(obj.list,function(x) PercentageFeatureSet(x,pattern='^MT-',col.name= 'percent.mt'))
obj.list <- lapply(obj.list, function(x) NormalizeData(x))

#QC
if(T){
for(i in 1:length(obj.list)){
  assign(paste('plot_',as.character(i),'_1',sep = ''),FeatureScatter(obj.list[[i]],feature1='nCount_RNA',feature2='nFeature_RNA',group.by='orig.ident',pt.size=0.1))
  assign(paste('plot_',as.character(i),'_2',sep = ''),FeatureScatter(obj.list[[i]],feature1='nCount_RNA',feature2='percent.mt',group.by='orig.ident',pt.size=0.1))
}
plot_1_1+plot_1_2+plot_2_1+plot_2_2+plot_3_1+plot_3_2 + plot_layout(ncol=2)
lapply(obj.list,function(x) summary(x[['percent.mt']]))#看mtDNA
}
lapply(obj.list.annot, function(x) VlnPlot(x,features = c('nCount_RNA','nFeature_RNA','percent.mt')))
#过滤细胞
obj.list[[1]] <- subset(obj.list[[1]],subset=nCount_RNA>200 & nCount_RNA<2000 & nFeature_RNA<1000 & percent.mt<30)
obj.list[[2]] <- subset(obj.list[[2]],subset=nCount_RNA>200 & nCount_RNA<5000 & nFeature_RNA<2000 & percent.mt<30)
obj.list[[3]] <- subset(obj.list[[3]],subset=nCount_RNA>200 & nCount_RNA<4000 & nFeature_RNA<1400 & percent.mt<30)


#后续常规流程
obj.list <- lapply(obj.list, function(x) FindVariableFeatures(x,nfeatures=2000))
obj.list <- lapply(obj.list, function(x) ScaleData(x,vars.to.regress=c('percent.mt','nCount_RNA')))
obj.list <- lapply(obj.list, function(x) RunPCA(x,verbose = F))
CombinePlots(plots = lapply(obj.list,function(x) ElbowPlot(x,ndims = 30)))#看看用多少个PC
pc.vector <- c(10,20,20) #1用10个PC，2，3用20个PC
obj.list <- lapply(1:3, function(x) RunUMAP(obj.list[[x]],dims=1:pc.vector[x]))
obj.list <- lapply(1:3, function(x) FindNeighbors(obj.list[[x]],dims=1:pc.vector[x],verbose=T))
obj.list <- lapply(obj.list, function(x) FindClusters(x,resolution=0.5,verbose=T))

#亚型注释
require(SingleR)
hpca <- HumanPrimaryCellAtlas()
pred.list <- lapply(obj.list, function(x) SingleR(test=x@assays$RNA@data,
                                                  ref=hpca,method='cluster',fine.tune=T,labels=hpca$label.main,cluster=x@meta.data$seurat_clusters))
id.list <- lapply(pred.list, function(x) x$pruned.labels)
for(i in 1:length(id.list)){
  names(id.list[[i]]) <- rownames(pred.list[[i]])
}
CombinePlots(plots=lapply(obj.list,function(x) DimPlot(x,label=T,repel=T)+NoLegend()),ncol=3)
id.list[[2]]['9'] <- 'Mast cell' #TPSAB1
id.list[[3]]['3'] <- 'Mast cell' #TPSAB1
id.list[[1]]['5'] <- 'Paneth cell' #C1QB C1QA
id.list[[2]]['5'] <- 'Paneth cell' #C1QB C1QA
id.list[[3]][c('14','4')] <- 'Paneth cell' #C1QB C1QA
id.list[[1]]['3'] <- 'Enterocyte progenitor cell' #ALDOB
id.list[[2]][c('0','6','12')] <- 'Enterocyte progenitor cell' #ALDOB
id.list[[3]]['11'] <- 'Enterocyte progenitor cell' #ALDOB
id.list[[1]]['8'] <- 'Fibroblast' #DCN LUM
id.list[[2]]['7'] <- 'Fibroblast' #DCN LUM
id.list[[3]][c('0','9','13')] <- 'Fibroblast' #DCN LUM
id.list[[2]]['10'] <- 'Smooth muscle cell' #TAGLN
id.list[[3]]['6'] <- 'Smooth muscle cell' #TAGLN

# 注释用的marker来源于
# CellMarker
# panglaodb
# wikipedia

obj.list.annot <- lapply(1:3, function(x) RenameIdents(obj.list[[x]],id.list[[x]]))#注释


#画图看结果
intestine.section <- c('Duodenum','Jejunum','Ileum')
dimplots <- lapply(obj.list.annot,function(x) DimPlot(x,label=T,repel=T)+NoLegend())
dimplots <- lapply(1:3,function(x) dimplots[[x]]+labs(title=intestine.section[x])+theme(plot.title = element_text(hjust = 0.5)))
CombinePlots(plots=dimplots,ncol=3)

# vlnplots <- lapply(obj.list.annot,function(x) VlnPlot(x,features = 'SLC6A19'))
# vlnplots <- lapply(1:3,function(x) vlnplots[[x]]+labs(title=intestine.section[x])+theme(plot.title = element_text(hjust = 0.5)))
# CombinePlots(plots=vlnplots,ncol=3)

dotplots<- lapply(obj.list.annot,function(x) DotPlot(x,features = c('SLC6A19','CTSL','ADAM17','TMPRSS2','ACE2'))+RotatedAxis())
dotplots <- lapply(1:3,function(x) dotplots[[x]]+labs(title=intestine.section[x])+theme(plot.title = element_text(hjust = 0.5)))
CombinePlots(plots=dotplots,ncol=3)

pdf('../plot/celltype.pdf',width=15,height=5)
CombinePlots(plots=dimplots,ncol=3)
dev.off()

pdf('../plot/dotplot.pdf',width=18,height=4)
CombinePlots(plots=dotplots,ncol=3)
dev.off()

#找marker
marker.list.celltype <- lapply(obj.list.annot,function(x) FindAllMarkers(x))
marker.list.seurat.clusters <- lapply(obj.list,function(x) FindAllMarkers(x))

#检查marker会用到的代码
#lapply(marker.list.celltype,function(x) head(x[x$cluster=='B_cell',],10))
#Reduce(intersect,lapply(marker.list.celltype,function(x) head(x[x$cluster=='B_cell','gene'],10)))
#CombinePlots(plots=lapply(obj.list.annot,function(x) FeaturePlot(x,features='SSR4')),ncol=3)
#marker.temp <- FindMarkers(obj.list[[3]],ident.1='13') 检查特定一类的marker
