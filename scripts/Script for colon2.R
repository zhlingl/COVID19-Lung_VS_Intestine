#上一次尝试将结肠三个部分一起分析，后来发现批次效应太过严重，根本不能一起注释
#因此现在重写脚本，分别对三个样本进行分析
require(ggplot2)
require(Seurat)
require(patchwork)
setwd('project/chenwx/scRNA_analysis/data/data_colon')
colon_list <- as.vector(unlist(read.table('../colon_list.txt',stringsAsFactors=F)))
data.list <- lapply(colon_list,function(x) read.table(paste(x,list.files(x),sep='/'),
                                                      stringsAsFactors=F,header=T,row.names=1))
names(data.list) <- colon_list

obj.list <- lapply(1:3,function(x) CreateSeuratObject(data.list[[x]],project=names(data.list[x]),min.cells=3))
obj.list <- lapply(obj.list,function(x) PercentageFeatureSet(x,pattern='^MT-',col.name= 'percent.mt'))
obj.list <- lapply(obj.list, function(x) NormalizeData(x))

#QC
for(i in 1:length(obj.list)){
  assign(paste('plot_',as.character(i),'_1',sep = ''),FeatureScatter(obj.list[[i]],feature1='nCount_RNA',feature2='nFeature_RNA',group.by='orig.ident',pt.size=0.1))
  assign(paste('plot_',as.character(i),'_2',sep = ''),FeatureScatter(obj.list[[i]],feature1='nCount_RNA',feature2='percent.mt',group.by='orig.ident',pt.size=0.1))
}
plot_1_1+plot_1_2+plot_2_1+plot_2_2+plot_3_1+plot_3_2 + plot_layout(ncol=2)
#过滤细胞
obj.list[[1]] <- subset(obj.list[[1]],subset=nCount_RNA>200 & nCount_RNA<2000 & nFeature_RNA<800 & percent.mt<40)
obj.list[[2]] <- subset(obj.list[[2]],subset=nCount_RNA>200 & nCount_RNA<2000 & nFeature_RNA<1000 & percent.mt<40)
obj.list[[3]] <- subset(obj.list[[3]],subset=nCount_RNA>200 & nCount_RNA<4000 & nFeature_RNA<1800 & percent.mt<40)


#后续常规流程
obj.list <- lapply(obj.list, function(x) FindVariableFeatures(x,nfeatures=2000))
obj.list <- lapply(obj.list, function(x) ScaleData(x,vars.to.regress=c('percent.mt','nCount_RNA')))
obj.list <- lapply(obj.list, function(x) RunPCA(x,verbose = F))
CombinePlots(plots = lapply(obj.list,function(x) ElbowPlot(x,ndims = 30)))#看看用多少个PC：20个

pc.vector <- c(15,15,20) #1,2用15个PC，2，3用20个PC
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
id.list[[1]]['6'] <- 'Fibroblast' #DCN LUM
id.list[[2]][c('2','8')] <- 'Fibroblast' #DCN LUM
id.list[[3]][c('9','11')] <- 'Fibroblast' #DCN LUM
id.list[[1]]['8'] <- 'Smooth muscle cell' #TAGLN
id.list[[2]]['9'] <- 'Smooth muscle cell' #TAGLN
id.list[[3]]['12'] <- 'Smooth muscle cell' #TAGLN
id.list[[1]]['9'] <- 'Mast cell' #TPSAB1
id.list[[2]]['6'] <- 'Mast cell' #TPSAB1
id.list[[2]]['5'] <- 'Neurons' #S100B
id.list[[1]]['10'] <- 'Neurons' #S100B
id.list[[1]]['7'] <- 'Epithelial cell' #OLFM4 AGR2 与另外两批数据对比得到 

# 注释用的marker来源于
# CellMarker
# panglaodb
# wikipedia


obj.list.annot <- lapply(1:3, function(x) RenameIdents(obj.list[[x]],id.list[[x]]))


#画图看结果
colon.section <- c('Ascending Colon','Sigmoid Colon','Transverse Colon')
dimplots <- lapply(obj.list.annot,function(x) DimPlot(x,label=T,repel=T)+NoLegend())
dimplots <- lapply(1:3,function(x) dimplots[[x]]+labs(title=colon.section[x])+theme(plot.title = element_text(hjust = 0.5)))
CombinePlots(plots=dimplots,ncol=3)

# vlnplots <- lapply(obj.list.annot,function(x) VlnPlot(x,features = 'SLC6A19'))
# vlnplots <- lapply(1:3,function(x) vlnplots[[x]]+labs(title=colon.section[x])+theme(plot.title = element_text(hjust = 0.5)))
# CombinePlots(plots=vlnplots,ncol=3)

dotplots<- lapply(obj.list.annot,function(x) DotPlot(x,features = c('SLC6A19','CTSL','ADAM17','TMPRSS2','ACE2'))+RotatedAxis())
dotplots <- lapply(1:3,function(x) dotplots[[x]]+labs(title=colon.section[x])+theme(plot.title = element_text(hjust = 0.5)))
CombinePlots(plots=dotplots,ncol=3)

pdf('../plot_colon/celltype.pdf',width=15,height=5)
CombinePlots(plots=dimplots,ncol=3)
dev.off()

pdf('../plot_colon/dotplot.pdf',width=18,height=4)
CombinePlots(plots=dotplots,ncol=3)
dev.off()

#找marker
marker.list.celltype <- lapply(obj.list.annot,function(x) FindAllMarkers(x))
marker.list.seurat.clusters <- lapply(obj.list,function(x) FindAllMarkers(x))

#检查marker会用到的代码
if(F){
lapply(marker.list.celltype,function(x) head(x[x$cluster=='B_cell',],10))
Reduce(intersect,lapply(marker.list.celltype,function(x) head(x[x$cluster=='B_cell','gene'],10)))
CombinePlots(plots=lapply(obj.list.annot,function(x) FeaturePlot(x,features='SSR4')),ncol=3)
marker.temp <- FindMarkers(obj.list[[3]],ident.1='13') #检查特定一类的marker
}
