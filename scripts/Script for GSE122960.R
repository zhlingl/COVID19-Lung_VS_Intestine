#因为服务器上不能用Read10X_h5函数（hdf5r包要用sudo安装），所以现在本地读取了之后再上传服务器
#读取数据
if(F){
setwd('project/chenwx/scRNA_analysis/data/data_lung')
library(Seurat)
data.list <- lapply(list.files(), function(x) Read10X_h5(paste(x,list.files(x)[1],sep='/')))
names(data.list) <- list.files()
saveRDS(data.list,'../GSE122960_mat_list_filtered.rds')
}

#上传之后
all.list <- readRDS('../GSE122960_mat_list_filtered.rds')
obj.list <- lapply(1:8,function(x) CreateSeuratObject(all.list[[x]],project=names(all.list[x]),min.cells=3,min.features=200))
obj.list <- lapply(obj.list,function(x) PercentageFeatureSet(x,pattern='^MT-',col.name= 'percent.mt'))
#QC
obj.list <- lapply(obj.list,function(x) subset(x, subset = nFeature_RNA < 6000 & percent.mt < 10))
#merge
obj.list <- lapply(obj.list, function(x) NormalizeData(x))
all.merge <- merge(obj.list[[1]],obj.list[2:8],merge.data=T)


all.merge <- FindVariableFeatures(all.merge,nfeatures=2000)
all.merge <- ScaleData(all.merge,vars.to.regress=c('percent.mt','nCount_RNA'))
all.merge <- RunPCA(all.merge,verbose = F)
ElbowPlot(all.merge)#看看用多少个PC：16个
all.merge <- RunUMAP(all.merge,dims=1:16)
all.merge <- FindNeighbors(all.merge,dims=1:16,verbose=T)
all.merge <- FindClusters(all.merge,resolution=0.5,verbose=T)
DimPlot(all.merge)

#亚型注释
library(SingleR)
hpca <- HumanPrimaryCellAtlas()
pred <- SingleR(test=all.merge@assays$RNA@data,
                ref=hpca,method='cluster',fine.tune=T,labels=hpca$label.main,cluster=all.merge@meta.data$seurat_clusters)
id.list <- pred$pruned.labels
names(id.list) <- rownames(pred)
all.merge.annot <- RenameIdents(all.merge,id.list)


#修改上皮细胞的注释
id.list[c('0','2','3','4','11','19')] <- 'AT2' #SFTPC
id.list['9'] <- 'AT1' #AGER
id.list['17'] <- 'Ciliated Cells' #TPPP3
id.list['5'] <- 'Club Cells' #SLGB3A2
id.list[c('1','6','7','8','10','12','21','22')] <- 'Dust Cells' #lung macrophage
# 注释用的marker来源于
# CellMarker
# panglaodb
# wikipedia

all.merge.reannot <- RenameIdents(all.merge,id.list)

#补充一些分组信息
#性别
male.list <- c('GSM3489185','GSM3489197')
all.merge.reannot[['sex']] <- ifelse(all.merge.reannot@meta.data$orig.ident %in% male.list,'male','female')
#sample
recode.sample <- unique(all.merge.reannot@meta.data$orig.ident)
all.merge.reannot[['sample']] <- unlist(lapply(all.merge.reannot@meta.data$orig.ident,
                                                    function(x) which(recode.sample==x)))

#画图
library(ggplot2)
plot.celltype <- DimPlot(all.merge.reannot,label=T,repel=T) + labs(title='Cell Type') + 
  theme(plot.title = element_text(hjust = 0.5),legend.position='none')
plot.dot <- DotPlot(all.merge.reannot,features=c('CTSL','ADAM17','TMPRSS2','ACE2')) + RotatedAxis()

# vln.list <- lapply(c('ACE2','TMPRSS2','ADAM17','CTSL'),function(x) VlnPlot(all.merge.reannot,features = x,pt.size = 0,split.by = 'sex',y.max = 4) + geom_jitter(alpha=0.2,size=0.5))
# plot.vln <- CombinePlots(plots = vln.list,ncol = 2)

pdf('D:/2019n-CoV/出图/lung_fig1_fixed_GSE122960.pdf',width = 5,height = 5)
plot.celltype
dev.off()

pdf('D:/2019n-CoV/出图/lung_fig3_fixed_GSE122960.pdf',width = 6,height = 4)
plot.dot
dev.off()

#pdf('D:/2019n-CoV/出图/lung_fig2_fixed_GSE122960.pdf',width = 12,height = 10)
#plot.vln
#dev.off()


#统计阳性细胞比例
type.list <- as.vector(unique(all.merge.reannot@active.ident))
lapply(type.list,function(x) sum(all.merge.reannot@assays$RNA@counts['ACE2',all.merge.reannot@active.ident==x]>0)/table(all.merge.reannot@active.ident)[x])

