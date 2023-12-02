library("qs")
library("EnhancedVolcano")
# the reason of large number of neutrophil--doublets

#------------DEG&GO for cell types Oct27 
#---------input
int_flash<-qread("/fs/ess/PCON0022/guoqi/microbe_sci/results/integrate/mit20/int_annotation_0614.qs")
Idents(int_flash)<-int_flash$celltype
DefaultAssay(int_flash)<-'RNA'

#Astrocytes
umap_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Astrocytes/")
celltype<-c("Astrocytes")
deg_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Astrocytes/DEG")
go_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Astrocytes/GO")


#B cells
umap_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Bcells/")
celltype<-c("B cells")
deg_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Bcells/DEG")
go_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Bcells/GO")


#Monocytes 
umap_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Monocytes/")
celltype<-c("Monocytes")
deg_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Monocytes/DEG")
go_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Monocytes/GO")
#-------highlight cell types
Idents(int_flash)<-int_flash$celltype
highlight_ct_umap(umap_path,int_flash,celltype)
#--------deg volcano plot
deg_plot(deg_path,int_flash,celltype)
#-------GO_function
goenrih_function(celltype,deg_path,go_path)
#------go plot
library(tidyverse)
go_plot(go_path)


#Ependymal cells
#---------input
umap_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Ependymalcells/")
celltype<-c("Ependymal cells")
deg_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Ependymalcells/DEG")
go_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Ependymalcells/GO")
#-------highlight cell types
Idents(int_flash)<-int_flash$celltype
highlight_ct_umap(umap_path,int_flash,celltype)
#--------deg volcano plot
deg_plot(deg_path,int_flash,celltype)
#-------GO_function
goenrih_function(celltype,deg_path,go_path)
#------go plot
library(tidyverse)
go_plot(go_path)



#Endothelial cells
#---------input
umap_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Endothelialcells/")
celltype<-c("Endothelial cells")
deg_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Endothelialcells/DEG")
go_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Endothelialcells/GO")
#-------highlight cell types
Idents(int_flash)<-int_flash$celltype
highlight_ct_umap(umap_path,int_flash,celltype)
#--------deg volcano plot
deg_plot(deg_path,int_flash,celltype)
#-------GO_function
goenrih_function(celltype,deg_path,go_path)
#------go plot
library(tidyverse)
go_plot(go_path)


#Erythroid cells
#---------input
umap_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Erythroidcells/")
celltype<-c("Erythroid cells")
deg_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Erythroidcells/DEG")
go_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Erythroidcells/GO")
#-------highlight cell types
Idents(int_flash)<-int_flash$celltype
highlight_ct_umap(umap_path,int_flash,celltype)
#--------deg volcano plot
deg_plot(deg_path,int_flash,celltype)
#-------GO_function
goenrih_function(celltype,deg_path,go_path)
#------go plot
library(tidyverse)
go_plot(go_path)



#Intermediate progenitors
#---------input
umap_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Intermediateprogenitors/")
celltype<-c("Intermediate progenitors")
deg_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Intermediateprogenitors/DEG")
go_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Intermediateprogenitors/GO")
#-------highlight cell types
Idents(int_flash)<-int_flash$celltype
highlight_ct_umap(umap_path,int_flash,celltype)
#--------deg volcano plot
deg_plot(deg_path,int_flash,celltype)
#-------GO_function
goenrih_function(celltype,deg_path,go_path)
#------go plot
library(tidyverse)
go_plot(go_path)

#Microglia
#---------input
umap_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Microglia/")
celltype<-c("Microglia")
deg_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Microglia/DEG")
go_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Microglia/GO")
#-------highlight cell types
Idents(int_flash)<-int_flash$celltype
highlight_ct_umap(umap_path,int_flash,celltype)
#--------deg volcano plot
deg_plot(deg_path,int_flash,celltype)
#-------GO_function
goenrih_function(celltype,deg_path,go_path)
#------go plot
library(tidyverse)
go_plot(go_path)

#MDMs (error 0 cell)
#---------input
umap_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/MDMs/")
celltype<-c("MDMs")
deg_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/MDMs/DEG")
go_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/MDMs/GO")
#-------highlight cell types
Idents(int_flash)<-int_flash$celltype
highlight_ct_umap(umap_path,int_flash,celltype)
#--------deg volcano plot
deg_plot_injury(deg_path,int_flash,celltype)
#-------GO_function
goenrih_function(celltype,deg_path,go_path)
#------go plot
library(tidyverse)
go_plot(go_path)

#Oligodendrocyte lineage (error 0 cells)
#---------input
umap_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Oligodendrocyte lineage/")
celltype<-c("Oligodendrocyte lineage")
deg_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Oligodendrocyte lineage/DEG")
go_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Oligodendrocyte lineage/GO")
#-------highlight cell types
Idents(int_flash)<-int_flash$celltype
highlight_ct_umap(umap_path,int_flash,celltype)
#--------deg volcano plot
deg_plot_injury(deg_path,int_flash,celltype)
#-------GO_function
goenrih_function(celltype,deg_path,go_path)
#------go plot
library(tidyverse)
go_plot(go_path)


#Pericytes error (Cell group 2 has fewer than 3 cells)
#---------input
umap_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Pericytes/")
celltype<-c("Pericytes")
deg_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Pericytes/DEG")
go_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Pericytes/GO")
#-------highlight cell types
Idents(int_flash)<-int_flash$celltype
highlight_ct_umap(umap_path,int_flash,celltype)
#--------deg volcano plot
deg_plot(deg_path,int_flash,celltype)
#-------GO_function
goenrih_function(celltype,deg_path,go_path)
#------go plot
library(tidyverse)
go_plot(go_path)

#Platelets
umap_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Platelets/")
celltype<-c("Platelets")
deg_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Platelets/DEG")
go_path<-c("/fs/ess/PCON0022/guoqi/microbe_sci/results/subset/Platelets/GO")



#-----------------------function
#---umap
highlight_ct_umap<-function(umap_path,object,celltype){
  setwd(umap_path)
  cell_names <-  WhichCells(object, idents = celltype)
  ct_umap<-DimPlot(object = int_flash, cells.highlight = cell_names, order = TRUE,pt.size=0.5)+ 
    scale_color_manual(labels = c("Others", celltype), values = c("grey", "black")) +
    labs(color = "legend title")
  ggsave(
    filename = paste0("UMAP_",celltype,"_highlighted.tiff"),
    plot = ct_umap,
    device = "tiff",
    width = 12,
    height = 10,
    units = "in",
    dpi = 150
  )
}

#------deg and deg plot
deg_plot<-function(deg_path,object,cell_type){
  setwd(deg_path)
  Idents(object)<-object$celltype
  ct_object<-subset(object,idents=cell_type)
  Idents(ct_object)<-ct_object$spinalcord
  #sc
  ct_object_sc<-subset(ct_object,idents="SC")
  Idents(ct_object_sc)<-ct_object_sc$microbiota
  deg_ct_flash_sc<-FindMarkers(ct_object_sc,ident.1 = "Colonized germ",ident.2 = "Germ free")
  deg_ct_flash_sc$gene<-rownames(deg_ct_flash_sc)
  #voline plot
  library(EnhancedVolcano)
  volcanoplot<-EnhancedVolcano(deg_ct_flash_sc, lab = deg_ct_flash_sc$gene,
                               x = 'avg_log2FC',
                               y = 'p_val_adj',
                               #selectLab = c('Fyb','Ngp'),
                               subtitle=NULL,
                               pCutoff =0.05,
                               FCcutoff = 1,
                               title="Colonized_SC vs Germfree_SC",
                               xlim = c(-2.5,2.5),
                               boxedLabels = TRUE,
                               colAlpha = 4/5,
                               legendLabSize = 14,
                               legendIconSize = 4.0,
                               drawConnectors = TRUE,
                               #widthConnectors = 1.0,
                               colConnectors = 'black')
  ggsave(
    plot = volcanoplot,
    filename = paste0("DEG_",cell_type,"_SC_volcano.tiff"),
    device = "tiff",
    dpi = 150,
    width = 10,
    height = 10,
    units = "in"
  )
  
  #-------sci
  ct_object_sci<-subset(ct_object,idents="SCI")
  Idents(ct_object_sci)<-ct_object_sci$microbiota
  deg_ct_flash_sci<-FindMarkers(ct_object_sci,ident.1 = "Colonized germ",ident.2 = "Germ free")
  deg_ct_flash_sci$gene<-rownames(deg_ct_flash_sci)
  volcanoplot2<-EnhancedVolcano(deg_ct_flash_sci, lab = deg_ct_flash_sci$gene,
                                x = 'avg_log2FC',
                                y = 'p_val_adj',
                                #selectLab = c('Fyb','Ngp'),
                                subtitle=NULL,
                                pCutoff =0.05,
                                FCcutoff = 1,
                                title="Colonized_SCI vs Germfree_SCI",
                                xlim = c(-2.5,4),
                                boxedLabels = TRUE,
                                colAlpha = 4/5,
                                legendLabSize = 14,
                                legendIconSize = 4.0,
                                drawConnectors = TRUE,
                                #widthConnectors = 1.0,
                                colConnectors = 'black')
  ggsave(
    plot = volcanoplot2,
    filename = paste0("DEG_",cell_type,"_SCI_volcano.tiff"),
    device = "tiff",
    dpi = 150,
    width = 10,
    height = 10,
    units = "in"
  )
  
  #separate up and down significant genes
  deg_ct_flash_sc_0.05<-deg_ct_flash_sc[deg_ct_flash_sc$p_val_adj<0.05,]
  deg_ct_flash_sci_0.05<-deg_ct_flash_sci[deg_ct_flash_sci$p_val_adj<0.05,]
  deg_ct_sc_up<-deg_ct_flash_sc_0.05[deg_ct_flash_sc_0.05$avg_log2FC>0,]
  deg_ct_sc_down<-deg_ct_flash_sc_0.05[deg_ct_flash_sc_0.05$avg_log2FC<0,]
  deg_ct_sci_up<-deg_ct_flash_sci_0.05[deg_ct_flash_sci_0.05$avg_log2FC>0,]
  deg_ct_sci_down<-deg_ct_flash_sci_0.05[deg_ct_flash_sci_0.05$avg_log2FC<0,]
  #------save deg
  write.csv(deg_ct_flash_sci,paste0("deg_",cell_type,"_sci_all.csv"))
  write.csv(deg_ct_flash_sc,paste0("deg_",cell_type,"_sc_all.csv"))
  
  write.csv(deg_ct_sci_up,paste0("deg_",cell_type,"_sci_up_0.05.csv"))
  write.csv(deg_ct_sci_down,paste0("deg_",cell_type,"_sci_down_0.05.csv"))
  write.csv(deg_ct_sc_up,paste0("deg_",cell_type,"_sc_up_0.05.csv"))
  write.csv(deg_ct_sc_down,paste0("deg_",cell_type,"_sc_down_0.05.csv"))
}
deg_plot_injury<-function(deg_path,object,cell_type){
  setwd(deg_path)
  Idents(object)<-object$celltype
  ct_object<-subset(object,idents=cell_type)
  Idents(ct_object)<-ct_object$spinalcord
  #-------sci
  ct_object_sci<-subset(ct_object,idents="SCI")
  Idents(ct_object_sci)<-ct_object_sci$microbiota
  deg_ct_flash_sci<-FindMarkers(ct_object_sci,ident.1 = "Colonized germ",ident.2 = "Germ free")
  deg_ct_flash_sci$gene<-rownames(deg_ct_flash_sci)
  volcanoplot2<-EnhancedVolcano(deg_ct_flash_sci, lab = deg_ct_flash_sci$gene,
                                x = 'avg_log2FC',
                                y = 'p_val_adj',
                                #selectLab = c('Fyb','Ngp'),
                                subtitle=NULL,
                                pCutoff =0.05,
                                FCcutoff = 1,
                                title="Colonized_SCI vs Germfree_SCI",
                                xlim = c(-2.5,2.5),
                                boxedLabels = TRUE,
                                colAlpha = 4/5,
                                legendLabSize = 14,
                                legendIconSize = 4.0,
                                drawConnectors = TRUE,
                                #widthConnectors = 1.0,
                                colConnectors = 'black')
  ggsave(
    plot = volcanoplot2,
    filename = paste0("DEG_",cell_type,"_SCI_volcano.tiff"),
    device = "tiff",
    dpi = 150,
    width = 10,
    height = 10,
    units = "in"
  )
  
  #separate up and down significant genes
  deg_ct_flash_sci_0.05<-deg_ct_flash_sci[deg_ct_flash_sci$p_val_adj<0.05,]
  deg_ct_sci_up<-deg_ct_flash_sci_0.05[deg_ct_flash_sci_0.05$avg_log2FC>0,]
  deg_ct_sci_down<-deg_ct_flash_sci_0.05[deg_ct_flash_sci_0.05$avg_log2FC<0,]
  #------save deg
  write.csv(deg_ct_flash_sci,paste0("deg_",cell_type,"_sci_all.csv"))
  write.csv(deg_ct_sci_up,paste0("deg_",cell_type,"_sci_up_0.05.csv"))
  write.csv(deg_ct_sci_down,paste0("deg_",cell_type,"_sci_down_0.05.csv"))
}
#-------GO_function
#GO (all enrichment are less than 0.05) (remember to upload to github)
GOenrichment_mouse <- function(Seurat.DEGs = NULL) {
  library(org.Mm.eg.db)
  library(clusterProfiler)
  genes.use <- Seurat.DEGs
  #Obtain a GO object
  GO_pathway <-
    enrichGO(
      gene = genes.use,
      OrgDb = org.Mm.eg.db,
      ont = "BP",
      keyType = "SYMBOL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.4
    )
  if (is.null(GO_pathway)) {
    GO_simplied_res <-
      as.data.frame(cbind(NO_results = "NO_result_for_GO"))
  } else{
    #GO_res <- simplify(GO_pathway)
    GO_simplied_res <-
      GO_pathway@result[GO_pathway@result$p.adjust < 0.05, ]
    dim(GO_simplied_res)
    if (dim(GO_simplied_res)[1] == 0) {
      GO_simplied_res <-
        as.data.frame(cbind(NO_results = "NO_result_for_GO"))
    }
  }
  return(GO_simplied_res)
}
goenrih_function<-function(cell_type,deg_path,go_path){
  down_degs<-list.files(path=deg_path,pattern = "down")
  up_degs<-list.files(path=deg_path,pattern = "up")
  deg_files<-c(up_degs,down_degs)
  #read deg files
  for(i in seq_along(deg_files)){
    setwd(deg_path)
    deg_file<-read.csv(deg_files[i])
    go_file<-GOenrichment_mouse(deg_file$gene)
    name<-strsplit(deg_files[i],split = "_")
    setwd(go_path)
    write.csv(go_file,paste("GO",name[[1]][2],name[[1]][3],name[[1]][4],name[[1]][5],sep = "_"))
  }
}
#------go plot
go_plot<-function(go_path){
  setwd(go_path)
  go_filename<-list.files(go_path,pattern = ".csv")
  for(i in seq_along(go_filename)){
    data<-read.csv(go_filename[i])
    if(ncol(data)!=2){
      temp_sort<-data[with(data, order(Count,decreasing = T)), ]
      temp_sort$Description<-factor(temp_sort$Description,levels=(temp_sort$Description))
      count_sum<-as.numeric(unlist(strsplit(temp_sort$GeneRatio[1],split="/"))[2])
      temp_sort$GeneRatio_num<-c(temp_sort$Count)/count_sum
      library(ggplot2)
      p1=ggplot(temp_sort[1:10,], # you can replace the numbers to the row number of pathway of your interest
                aes(x = GeneRatio_num, y = Description)) + 
        geom_point(aes(size = Count, color = p.adjust)) +
        theme_bw(base_size = 14) +
        theme(axis.text = element_text(size = 14, face = "bold"),
        )+
        scale_colour_gradient(limits=NULL, low="red",high="blue") +
        ylab(NULL) +
        ggtitle("GO enrichment")
      #define name
      plot_filename<-strsplit(go_filename[i],split = "[.]")[[1]][1]
      ggsave(
        plot = p1,
        filename = paste0(plot_filename,".05","_enrichment.tiff"),
        device = "tiff",
        dpi = 150,
        width = 13,
        height = 10,
        units = "in"
      )
    }
  }
}
