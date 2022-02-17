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
library(writexl)
library('FeatureImpCluster')
library(flexclust)
library(logisticPCA)
library(NbClust)
library(openxlsx)
library(patchwork)
library(grid)
library(gridExtra)


 
  #### Dataset preparation 6 clusters ####
  
  for (i in 1:6){
    assign(paste0('cluster',i),colnames(read.table(paste0('hallmark_genes_cluster',i,'.txt'),
                                                         sep="", header = T)))
  }

  genes_clusters <- list(cluster1=cluster1,cluster2=cluster2,
                            cluster3=cluster3,cluster4=cluster4,cluster5=cluster5,cluster6=cluster6)
  parent_genes <- colnames(read.csv('Oncofusion_proteins_project/Oncogene_fusions/Data_sets/oncofusion_genes.txt',
                           header = T,sep=""))
  
  df <- read_excel('Oncofusion_proteins_project/Oncogene_fusions/Data_sets/oncogene_fusions.xlsx')
  Fusion_Genes<- df %>% select(FusionGene, Hgene, Tgene)
  length(genes_clusters)


  Cluster_func <- function(dd, genes){
    dd$Cluster_Hgene <- NA
    dd$Cluster_Tgene <- NA
    for (i in 1: nrow(dd)){
      for (j in 1: length(genes)){
        if(is.element(dd$Hgene[i],unlist(genes[j]))){
          dd$Cluster_Hgene[i] <- j
        }
      }
    }
    for (k in 1: nrow(dd)){
      for (l in 1: length(genes)){
        if(is.element(dd$Tgene[k],unlist(genes[l]))){
          dd$Cluster_Tgene[k] <- l
        }
      }
    }
    return(dd)
  }

  Fusion_Genes_6 <- Cluster_func(Fusion_Genes,genes_clusters)
  Fusion_Genes_6$Cluster_comb <- as.factor(paste(Fusion_Genes_6$Cluster_Hgene,Fusion_Genes_6$Cluster_Tgene, sep="_"))  
  write_csv(Fusion_Genes,'Fusion_Genes_by_Clsuter.csv')
  Count_by_cluster_comb <- Fusion_Genes_6%>% group_by(Cluster_comb)%>% summarise(count= n())
  write.csv(Count_by_cluster_comb,'Count_by_6_cluster_comb.csv')
  View(Count_by_cluster_comb)
  
  Barplot_fusions <- ggplot(Count_by_cluster_comb)+geom_boxplot(mapping=aes(x=Cluster_comb,y=count))+ 
    ggtitle('Cluster Combinations in Fusion Genes') + xlab('Cluster fusions') + ylab('Count')+
    theme_linedraw()
  
  
  
  #### Random samples from the two parent genes list
  
  
  combinations <- crossing(var1=1:6,var2= 1:6)
  cluster_combinations <- paste(combinations$var1,combinations$var2, sep= '_')
  random_data <- data.frame(Cluster_Comb = as.factor(cluster_combinations))
  random_data$Cluster_Comb

  
  
  
  #### Fraction of between cluster fusions ####
  
  #### We take for further analysis 6 clusters:
  ## for real fusion genes: 
  Fraction_between_cluster_6 <-subset(Fusion_Genes_6, Cluster_Hgene != Cluster_Tgene)
  Fraction_between_cluster_6 <- length(Fraction_between_cluster_6$FusionGene)/length(Fusion_Genes_6$FusionGene)
  
  ## for random pairs: 
  
  N <- 1000
  Fraction_between_cluster_random_6 <- data.frame(matrix(nrow= N, ncol=1, dimnames = c(list(c(1:N)),list(c('Between_Cluster_Fraction')))))
  
  for (i in 1:N) {
    Hgene <- sample(parent_genes,length(Fusion_Genes_6$Hgene), replace = T)
    Tgene <- sample(parent_genes,length(Fusion_Genes_6$Hgene), replace = T) 
    random_data_set <- data.frame('Hgene'= Hgene,'Tgene'=Tgene)
    random_data_set <- Cluster_func(random_data_set,genes_clusters)
    Fraction <-subset(random_data_set, Cluster_Hgene != Cluster_Tgene)
    Fraction  <- length(Fraction$Hgene)/length(random_data_set$Hgene)
    Fraction_between_cluster_random_6$Between_Cluster_Fraction[i] <- Fraction
  }
  
  write.csv(Fraction_between_cluster_random, 'Fraction_between_6_cluster_random.csv')
  
  View(Fraction_between_cluster_random_6)
  
  
  
  ### standardizing the dataset
  
  stand_Fraction_between_cluster_random_6 <- Fraction_between_cluster_random%>% 
    mutate_at('Between_Cluster_Fraction',~(scale(.)%>% as.vector))
  Z_score_Fraction_between_cluster_6 <- (Fraction_between_cluster_6-mean(Fraction_between_cluster_random_6$Between_Cluster_Fraction))/sd(Fraction_between_cluster_random_6$Between_Cluster_Fraction)
  
  ## getting p-value of two sided test
  p_value <- 2*pnorm(abs(Z_score_Fraction_between_cluster_6),lower.tail = F)
  
  
  Between_cluster_hist <- ggplot(Fraction_between_cluster_random_6)+ geom_histogram(aes(x=Between_Cluster_Fraction),bins = 60, color = 'lightgrey', fill = 'grey') +
    theme_linedraw() +
    geom_vline(aes(xintercept=Fraction_between_cluster_6), color= 'red', size = 1)+
    geom_text(aes(x=Fraction_between_cluster_6, label=4e-13, y = 0), 
              colour ='grey', angle= 90, vjust=1.5, hjust=-2, size=5)+
    xlab('Fraction between cluster fusions')+
    ylab('Count')
  
  combinations <- crossing(var1=1:6,var2= 1:6)
  cluster_combinations <- paste(combinations$var1,combinations$var2, sep= '_')
  
  Fraction_data_set_6 <- data.frame(matrix(nrow= N+1, ncol=(length(cluster_combinations)+2), dimnames = c(list(c('True_pairs',c(1:N))),list(c('Between_Cluster','Within_Cluster',cluster_combinations)))))
  Standadized_data_set_6 <- data.frame(matrix(nrow= N+1, ncol=(length(cluster_combinations)+2), dimnames = c(list(c('True_pairs',c(1:N))),list(c('Between_Cluster','Within_Cluster',cluster_combinations)))))
  Final_data_set_6 <- data.frame(matrix(nrow=(length(cluster_combinations)+2),ncol= 2,dimnames = c(list(c('Between_Cluster','Within_Cluster',cluster_combinations)),list(c('Z_score','pvalue')))))
  
  
  ## data frames
  Fraction_data_set_6$Between_Cluster <- c(Fraction_between_cluster_6,Fraction_between_cluster_random_6$Between_Cluster_Fraction)
  Standadized_data_set_6$Between_Cluster <- c(Z_score_Fraction_between_cluster_6, stand_Fraction_between_cluster_random_6$Between_Cluster_Fraction)
  Final_data_set_6$Z_score[1] <- Z_score_Fraction_between_cluster_6
  Final_data_set_6$pvalue[1] <- p_value
  
  
  
  ## within clusters
  Subset <-subset(Fusion_Genes_6, Cluster_Hgene == Cluster_Tgene)
  Subset <- length(Subset$FusionGene)/length(Fusion_Genes_6$FusionGene)
  Fraction_data_set_6[1,]$Within_Cluster <- Subset
  
  for (i in 1:N) {
    Hgene <- sample(parent_genes,length(Fusion_Genes_6$Hgene), replace = T)
    Tgene <- sample(parent_genes,length(Fusion_Genes_6$Hgene), replace = T) 
    random_data_set <- data.frame('Hgene'= Hgene,'Tgene'=Tgene)
    random_data_set <- Cluster_func(random_data_set,genes_clusters)
    Fraction <-subset(random_data_set, Cluster_Hgene == Cluster_Tgene)
    Fraction  <- length(Fraction$Hgene)/length(random_data_set$Hgene)
    Fraction_data_set_6[-1,]$Within_Cluster[i] <- Fraction
  }
  
  stand<- Fraction_data_set_6%>% 
    mutate_at('Within_Cluster',~(scale(.)%>% as.vector))
  Standadized_data_set_6$Within_Cluster <- stand$Within_Cluster
  p_value <- 2*pnorm(abs(Standadized_data_set_6$Within_Cluster[1]), lower.tail = F)
  Final_data_set_6$Z_score[2] <- Standadized_data_set_6$Within_Cluster[1]
  Final_data_set_6$pvalue[2] <- p_value
  
  Within_cluster_hist <- ggplot(Fraction_data_set_6)+ geom_histogram(aes(x=Within_Cluster),bins = 60, color = 'lightgrey', fill = 'grey') +
    theme_linedraw() +
    geom_vline(aes(xintercept=Fraction_data_set_6$Within_Cluster[1]), color= 'red', size = 1)+
    geom_text(aes(x=Within_Cluster[1], label=4e-12, y = 0), 
              colour ='grey', angle= 90, vjust=1.5, hjust=-2, size=5)+
    xlab('Fraction within cluster fusions')+
    ylab('Count')
  
  ### for all cluster combinations
  for (j in 1:length(combinations$var1)){
    S <-subset(Fusion_Genes_6, Cluster_Hgene == combinations$var1[j] & Cluster_Tgene==combinations$var2[j])
    S <- length(S$FusionGene)/length(Fusion_Genes_6$FusionGene)
    Fraction_data_set_6[1,j+2] <- S
  }
  for (j in 1:length(combinations$var1)){
    print(j)
    for (i in 1:N) {
      Hgene <- sample(parent_genes,length(Fusion_Genes_6$Hgene), replace = T)
      Tgene <- sample(parent_genes,length(Fusion_Genes_6$Hgene), replace = T) 
      random_data_set <- data.frame('Hgene'= Hgene,'Tgene'=Tgene)
      random_data_set <- Cluster_func(random_data_set,genes_clusters)
      S <-subset(random_data_set, Cluster_Hgene == combinations$var1[j] & Cluster_Tgene==combinations$var2[j])
      S <- length(S$Hgene)/length(Fusion_Genes_6$Hgene)
      Fraction_data_set_6[-1,][i,j+2] <- S
    }
  }
  for (N in 3:length(colnames(Fraction_data_set_6))){
    stand<- Fraction_data_set_6%>% 
      mutate_at(colnames(Fraction_data_set_6)[N],~(scale(.)%>% as.vector))
    Standadized_data_set_6[,N] <- stand[,N]
    p_value <- 2*pnorm(abs(Standadized_data_set_6[1,N]), lower.tail = F)
    Final_data_set_6$Z_score[N] <- Standadized_data_set_6[1,N]
    Final_data_set_6$pvalue[N] <- p_value 
  }
  Mean <- rep(NA,length(colnames(Fraction_data_set_6)))
  pvalue_right <- rep(NA,length(colnames(Fraction_data_set)))
  
  for (N in 1:length(colnames(Fraction_data_set_6))){
    Mean[N] <- mean(Fraction_data_set_6[,N])
    
  }
  
  Final_data_set_6<- cbind(Mean,Final_data_set_6)
  
  
  all_data_sets_6 <- list('Fraction_per_combination'= Fraction_data_set_6, 'Standardized_data'=Standadized_data_set_6, 'Pvalues'=Final_data_set_6)
  write.xlsx(all_data_sets_6, 'Fractions_per_cluster_combination_6_cluster.xlsx', rowNames= TRUE, colNames = TRUE)
  
  ### count per H and T gene
  random_data_set$Cluster_Hgene <- as.factor(random_data_set$Cluster_Hgene)
  random_data_set$Cluster_Tgene <- as.factor(random_data_set$Cluster_Tgene)
  random_data_set %>%group_by(Cluster_Hgene,Cluster_Tgene )%>% summarise(count = n())
  
  summarise(random_data_set, count_H = n_distinct(Cluster_Hgene),count_T = n_distinct(Cluster_Tgene))
  random_data_set %>%group_by(Cluster_Hgene)%>% summarise(count = n())

   d <- count(random_data_set, Cluster_Hgene)
   t <- count(random_data_set, Cluster_Tgene)
  merge_obs$Hgene <- merge_obs$Hgene/840
  merge_obs$Tgene <- merge_obs$Tgene/840
  
  
  d <- count(Fusion_Genes_6, Cluster_Hgene)
  t <- count(Fusion_Genes_6, Cluster_Tgene)
  
  merge_true <- merge(d,t,by.x = "Cluster_Hgene", by.y = "Cluster_Tgene" )
  colnames(merge_true) <- c("Cluster", "Hgene", "Tgene")
  
  merge_true$Hgene <- merge_true$Hgene/840
  merge_true$Tgene <- merge_true$Tgene/840
  
  
  ### have a look at fusions beteween cluster 1 and 3
  
  Cluster_fusions_6_5_5 <- Fusion_Genes_6 %>% filter(Cluster_Hgene==5 & Cluster_Tgene==5)
  
  Clusterr_fusions_6_5_6 <- Fusion_Genes_6 %>% filter(Cluster_Hgene==5 & Cluster_Tgene==6)
  
  Clusterr_fusions_6_6_5 <- Fusion_Genes_6 %>% filter(Cluster_Hgene==6 & Cluster_Tgene==5)
  Clusterr_fusions_6_6_6 <- Fusion_Genes_6 %>% filter(Cluster_Hgene==6 & Cluster_Tgene==6)
  Clusterr_fusions_6_5_1 <- Fusion_Genes_6 %>% filter(Cluster_Hgene==5 & Cluster_Tgene==1)
  Clusterr_fusions_6_5_4 <- Fusion_Genes_6 %>% filter(Cluster_Hgene==5 & Cluster_Tgene==4)
  
  
  all_data_sets_6 <- list('Cluster_5_5'= Cluster_fusions_6_5_5, 
                          'Cluster_5_6'=Clusterr_fusions_6_5_6, 
                          'Cluster_5_4'=Clusterr_fusions_6_5_4,
                          'Cluster_5_1'=Clusterr_fusions_6_5_1,
                          'Cluster_6_5'=Clusterr_fusions_6_6_5,
                          'Cluster_6_6'=Clusterr_fusions_6_6_1)
  
  write.xlsx(all_data_sets_6, 'Gene_fusions_significant_clusters.xlsx', rowNames= F, colNames = TRUE)
  
  fusion_5_5 <- ggplot(Fraction_data_set_6)+ geom_histogram(aes(x=X5_5),bins = 60, color = 'lightgrey', fill = 'grey') +
    theme_linedraw() +
    geom_vline(aes(xintercept=X5_5[1]), color= 'red', size = 1)+
    geom_text(aes(x=Fraction_data_set_6$X5_5[1], label=5e-13, y = 0), 
              colour ='grey', angle= 90, vjust=1, hjust=0, size=4)+
    xlab('Fraction 5/5 cluster fusion')+
    ylab('Count')
  
  fusion_5_6 <- ggplot(Fraction_data_set_6)+ geom_histogram(aes(x=X5_6),bins = 50, color = 'lightgrey', fill = 'grey') +
    theme_linedraw() +
    geom_vline(aes(xintercept=X5_6[1]), color= 'red', size = 1)+
    geom_text(aes(x=Fraction_data_set_6$X5_6[1], label=3.7e-07, y = 0), 
              colour ='grey', angle= 90, vjust=1, hjust=0, size=4)+
    xlab('Fraction 5/6 cluster fusion')+
    ylab('Count')

  fusion_5_4 <- ggplot(Fraction_data_set_6)+ geom_histogram(aes(x=X5_4),bins = 40, color = 'lightgrey', fill = 'grey') +
    theme_linedraw() +
    geom_vline(aes(xintercept=X5_4[1]), color= 'red', size = 1)+
    geom_text(aes(x=Fraction_data_set_6$X5_4[1], label=0.00953, y = 0), 
              colour ='grey', angle= 90, vjust=2, hjust=-2, size=5)+
    xlab('Fraction 5/4 cluster fusion')+
    ylab('Count')
  
  fusion_5_1 <- ggplot(Fraction_data_set_6)+ geom_histogram(aes(x=X5_1),bins = 40, color = 'lightgrey', fill = 'grey') +
     theme_linedraw() +
    geom_vline(aes(xintercept=X5_1[1]), color= 'red', size = 1)+
    geom_text(aes(x=Fraction_data_set_6$X5_1[1], label=0.00464, y = 0), 
              colour ='grey', angle= 90, vjust=2, hjust=, size=5)+
    xlab('Fraction 5/1 cluster fusion')+
    ylab('Count')
  
  fusion_6_5 <- ggplot(Fraction_data_set_6)+ geom_histogram(aes(x=X6_5),bins = 40, color = 'lightgrey', fill = 'grey') +
    theme_linedraw() +
    geom_vline(aes(xintercept=X6_5[1]), color= 'red', size = 1)+
    geom_text(aes(x=Fraction_data_set_6$X6_5[1], label=8.8e-09, y = 0), 
              colour ='grey', angle= 90, vjust=1, hjust=-0, size=4)+
    xlab('Fraction 6/5 cluster fusion')+
    ylab('Count')
  
  
  
  
  ### Plots
  
  Barplot_fusions 
  Between_cluster_hist
  Within_cluster_hist
  fusion_5_5
  fusion_5_6 
  fusion_5_4 ### not
  fusion_5_1  ### not
  fusion_6_5 
  
  
  ggarrange(fusion_5_4,fusion_5_1,nrow=1, 
            labels = c("a", "b"),hjust=0.05,font.label = list(size = 13))
  
  fusions_5 <- fusion_5_6 |(fusion_5_5/fusion_6_5)
  fusions_5 <- fusions_5 +plot_annotation(title = 'Expeced fraction from randomization analysis')
  
  with_between <- Between_cluster_hist / Within_cluster_hist
  with_between <- with_between+
    plot_annotation(title = 'Expeced fraction from randomization analysis \nof within and between cluster fusions')
  
  ggarrange(fusion_5_4,fusion_5_1,nrow=1, 
            labels = c("a", "b"),hjust=0.05,font.label = list(size = 13))
  
  l <- ggarrange(Barplot_fusions,fusions_5,
            nrow= 2,labels = c("a","b"),
            hjust=0.05, vjust= 1, font.label = list(size = 13))
  
  n <- ggarrange(with_between,nrow = 1,labels = c("c"),
            hjust=0.05,font.label = list(size = 13))

  
  l+n + plot_layout(widths = c(2, 1.5))
 