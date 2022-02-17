  #### Packages ####
  library(tidyverse)
  library(dplyr)
  library(gdata)
  library(FactoMineR)
  library(factoextra)
  library(devtools)
  library(ggpubr)
  library(logisticPCA)
  library(glmpca)
  library(rstatix)
  library(corrr)
  library(cluster)
  library(readxl)
  library(plotly)
  library(writex)
  library('FeatureImpCluster')
  library(flexclust)
  library(logisticPCA)
  library(NbClust)
  library(clValid) 
  library(patchwork)
  library(grid)
  library(gridExtra)
  
  
  #### Dataset & Matrix #####
  Genes_Hallmark <- readRDS('Oncofusion_proteins_project/Data_sets/Genes_in_Hallmark_GOs.rds')
  unique(names(Genes_Hallmark))
  onco_genes <- read.table('Oncofusion_proteins_project/Oncogene_fusions/Data_sets/oncofusion_genes.txt', sep="",header= T)
  onco_genes <- unique(colnames(onco_genes))
  matrix_generation <- function(file, gene){
    dataSet<- file
    bin_matrix <-matrix(data=NA, nrow=length(gene), 
                        ncol = length(unique(dataSet)), dimnames = c(list(gene),list(unique(names(dataSet)))))
    for (i in 1:length(gene)){
      for (j in 1:length(unique(dataSet))){
        if(is.element(gene[i],unlist(unique(dataSet[j])))) {
          bin_matrix[i,j] <- 1
        }
        else {
          bin_matrix[i,j] <- 0
        }
      }
    }
    return(bin_matrix)
  }

  hallmark_GO_terms_matrix <- matrix_generation(Genes_Hallmark, onco_genes)
  write.table(hallmark_GO_terms_matrix,'hallmark_GO_terms_matrix.txt')
  hallmark_GO_terms_matrix <- read.delim('Oncofusion_proteins_project/Oncogene_fusions/hallmark_GO_terms_matrix.txt', header = T, sep = "")
  
  #### PCA & Clustering ####
  
  get_cluster_terms <- function(data_f,cluster){
    n_clust <- unique(cluster)
    data_f <- as.data.frame(hallmark_GO_terms_matrix_tot)
    GO_terms <- list()
    for (j in n_clust){
      dd <- data_f%>% filter(data_f[1] ==j)
      GO_terms[[paste0('cluster',j)]] <- colSums(dd[,-1]%>% select_if(colSums(dd[,-1])!=0))
      
    }
    return(GO_terms)
    
  }
  ### PCA
  set.seed(1)
  
  
  logpca_cv <- cv.lpca(as.matrix(hallmark_GO_terms_matrix), ks = 8, ms= 1:10) 
  log_pca <- logisticPCA(hallmark_GO_terms_matrix, k=8,m= which.min(logpca_cv))
  
  ### Clustering k-means
  pc_hallmark <- data.frame(log_pca$PCs)
  gap_statistics <- fviz_nbclust(pc_hallmark, kmeans,  method= "gap_stat") # 6 clusters

  km_hallmark <- eclust(pc_hallmark, "kmeans",k = 6, nstart = 25, graph = F)

  silhouette_plot <- fviz_silhouette(km_hallmark,print.summary = F)

  clust_hallmark <- as.factor(km_hallmark$cluster)
  

  pca_plot <- ggplot(pc)+ ggtitle('Parent proteins by PC1 and PC2')+geom_point(aes(x=V1, y= V2,colour = as.factor(clust_hallmark)))+ 
    theme_linedraw()+ labs(color = "Clusters")+ xlab('PC1')+ylab('PC2')
  
  
  ### Loadings
  hallmark_GO_terms_matrix_tot <- cbind(clust_hallmark,hallmark_GO_terms_matrix)
  
  hallmark_GO_terms <-get_cluster_terms(hallmark_GO_terms_matrix_tot,sort(clust_hallmark))
  
  ### Getting feature importance
  Hallmark_kcca <- as.kcca(km_hallmark,pc_hallmark)
  kcca_plot <- barplot(Hallmark_kcca) 

  barplot(Hallmark_kcca)
  
  ### To make data sets comparable
  loadings_Hallmark_PC1 <- as.data.frame(log_pca$U[,1]) 
  colnames(loadings_Hallmark_PC1) <- "PC1"
  rownames(loadings_Hallmark_PC1) <- unique(names(Genes_Hallmark))
  loadings_Hallmark_PC1 <- head(loadings_Hallmark_PC1%>% arrange(desc(abs(PC1))),200)
  loadings_Hallmark_PC1  <- rownames_to_column(loadings_Hallmark_PC1)
  loadings_Hallmark_PC2 <- as.data.frame(log_pca$U[,2]) 
  colnames(loadings_Hallmark_PC2) <- "PC2"
  rownames(loadings_Hallmark_PC2) <- unique(names(Genes_Hallmark))
  loadings_Hallmark_PC2 <- head(loadings_Hallmark_PC2%>% arrange(desc(abs(PC2))),200)
  loadings_Hallmark_PC2  <- rownames_to_column(loadings_Hallmark_PC2)
  loadings_Hallmark_PC3 <- as.data.frame(log_pca$U[,3]) 
  colnames(loadings_Hallmark_PC3) <- "PC3"
  rownames(loadings_Hallmark_PC3) <- unique(names(Genes_Hallmark))
  loadings_Hallmark_PC3 <- head(loadings_Hallmark_PC3%>% arrange(desc(abs(PC3))),200)
  loadings_Hallmark_PC3  <- rownames_to_column(loadings_Hallmark_PC3)
  

  ### Most important loadings plot
   PC1_loadings<- ggplot(loadings_Hallmark_PC1[1:9,],aes( x= PC1,y= rowname))+ geom_bar(stat = 'identity')+
    ggtitle("PC 1 Loadings") + xlab('r') + ylab('')
   PC2_loadings<- ggplot(loadings_Hallmark_PC2[1:6,],aes( x= PC2,y= rowname))+ geom_bar(stat = 'identity')+
     ggtitle("PC 2 Loadings") + xlab('r') + ylab('')
   PC3_loadings<- ggplot(loadings_Hallmark_PC3[1:6,],aes( x= PC3,y= rowname))+ geom_bar(stat = 'identity')+
     ggtitle("PC 3 Loadings") + xlab('r') + ylab('')
   PC_Loaddings <- PC1_loadings+PC2_loadings+PC3_loadings+ plot_layout(design=layout) +
     plot_annotation(title= 'Principal component loadings') 
   
   layout <- "
   AB
   AC
   "

  
  loading_cluster1 <-head(loadings_Hallmark_PC3[is.element(loadings_Hallmark_PC3$rowname,str_replace(rownames(as.data.frame(hallmark_GO_terms$cluster1)),'O.','O:'))==T,],15)
  loading_cluster2 <- head(loadings_Hallmark_PC1[is.element(loadings_Hallmark_PC1$rowname,str_replace(rownames(as.data.frame(hallmark_GO_terms$cluster2)),'O.','O:'))==T,],15)
  loading_cluster3 <- head(loadings_Hallmark_PC1[is.element(loadings_Hallmark_PC1$rowname,str_replace(rownames(as.data.frame(hallmark_GO_terms$cluster3)),'O.','O:'))==T,],15)
  loading_cluster4 <-head(loadings_Hallmark_PC2[is.element(loadings_Hallmark_PC2$rowname,str_replace(rownames(as.data.frame(hallmark_GO_terms$cluster4)),'O.','O:'))==T,],15)
  loading_cluster5 <- head(loadings_Hallmark_PC1[is.element(loadings_Hallmark_PC1$rowname,str_replace(rownames(as.data.frame(hallmark_GO_terms$cluster5)),'O.','O:'))==T,],10)
  loading_cluster6 <- head(loadings_Hallmark_PC1[is.element(loadings_Hallmark_PC1$rowname,str_replace(rownames(as.data.frame(hallmark_GO_terms$cluster6)),'O.','O:'))==T,],10)
  
  write_xlsx(list(cluster1 =loading_cluster1,cluster2=loading_cluster2,cluster3 =loading_cluster3,cluster4 =loading_cluster4), 'hallmarks_4_clusters_highest_loadings_top10.xlsx')
  

  
  ### Gene list for each cluster
  
  df_genes_cluster <- as.data.frame(hallmark_GO_terms_matrix_tot)
  for (i in sort(unique(df_genes_cluster$clust_hallmark))){
    write.table(t(rownames(filter(df_genes_cluster,clust_hallmark==i))),paste0('hallmark_genes_cluster',i,'.txt'),quote = F, sep = " ", row.names = F, col.names = F)
  }
  
  dim(filter(df_genes_cluster,clust_hallmark==1))
  
  ### Plots
  gap_statistics  
  silhouette_plot
  pca_plot
  PC_Loaddings 
  barplot(Hallmark_kcca)
  

  ggarrange(PC_Loaddings,ggarrange(silhouette_plot,pca_plot, nrow=2, 
                                   labels = c("b", "c"),hjust=0.05,font.label = list(size = 13)), nrow=1, 
            labels = c("a"),hjust=0.05,font.label = list(size = 13))
  
  
  

  
  
  