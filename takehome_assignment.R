#CS4115-Computational biology
#Clustering and trees_Take-home assignment
#S12764
#2015s15420

#installing packages using bioconductor
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('multtest')
BiocManager::install('pheatmap')
install.packages('ClassDiscovery')
install.packages('factoextra')
#installing package which include hierarchical clustering libraries
BiocManager::install('fastcluster')
#calling the installed library
library('multtest')
#getting the pathway
getwd()
#setting the pathway
setwd("C:/Users/thula/Documents/second semester-4th/Computational biology")

#loading in-built dataset
?golub
data(golub)
golub.gnames=as.matrix(golub.gnames)
row.names(golub)=golub.gnames[,3]
#checking the dimensions of the dataset
dim(golub)



#########################
#2.
#a)
##Hierarchical clustering on distance matrix built up using pearson correlation method
#Creating distance matrix using pearson correlation (using distanceMarix function in classDiscovery package libraries)
library('ClassDiscovery')
d=distanceMatrix(golub, 'pearson')
#applying hierarchical clustering on distance matrix cretaed using pearson correlation (hclust function in fastcluster package libraries)
library("fastcluster")
d_t=distanceMatrix(t(golub),'pearson')
hc_pearson <- hclust(d, method = "ward.D2", members=NULL)
names(hc_pearson)
hr_pearson <- hclust(d_t, method = "ward.D2", members=NULL)
#plotting the tree build using hierarchical clustering method

plot(hc_pearson,labels = NULL, hang = 0.1, check = TRUE,
     axes = TRUE, frame.plot = FALSE, ann = TRUE,
     main = "Cluster Dendrogram Pearson",
     sub = NULL, xlab = NULL, ylab = "Height")
plot(hr_pearson,labels = NULL, hang = 0.1, check = TRUE,
     axes = TRUE, frame.plot = FALSE, ann = TRUE,
     main = "Cluster Dendrogram Pearson",
     sub = NULL, xlab = NULL, ylab = "Height")

#b)
##Hierarchical clustering on distance matrix built up using eucledian distance based method
#Creating distance matrix using eucledian distance based method (using distanceMarix function in classDiscovery package libraries)
d_t=distanceMatrix(t(golub),'euclid')
d=distanceMatrix(golub, 'euclid')
#applying hierarchical clustering on distance matrix cretaed using eucledian distance based method (hclust function in fastcluster package libraries)
hc_euclid <- hclust(d, method = "ward.D2", members=NULL)
names(hc_euclid)
hr_euclid=hclust(d_t, method = "ward.D2", members=NULL)
#plotting the tree build using hierarchical clustering method
plot(hc_euclid,labels = NULL, hang = 0.1, check = TRUE,
     axes = TRUE, frame.plot = FALSE, ann = TRUE,
     main = "Cluster Dendrogram Euclidean",
     sub = NULL, xlab = NULL, ylab = "Height")
plot(hr_euclid,labels = NULL, hang = 0.1, check = TRUE,
     axes = TRUE, frame.plot = FALSE, ann = TRUE,
     main = "Cluster Dendrogram Euclidean",
     sub = NULL, xlab = NULL, ylab = "Height")
?hclust
############################
#3.
#Application of k-means clustering on dataset
#to select number of clusters elbow plots need to be drawn
#for that first principal component analysis was done on the dataset. Then elbow plot and scree plot was drawn
PCA=prcomp(golub,center = T, scale=T)
summary(PCA)
#Screeplot
screeplot(PCA, main = "Screeplot of first 10 PC's")
#abline of eigenvalue=2 - to select PC's cut by the abline
abline(h = 2, col="red")
#Elbow plot
plot(PCA, type="l", main = "Elbow plot of first 10 PC's")
#abline of eigenvalue=2 - to select PC's above the abline
abline(h = 2, col="red")

#K-means clustering was done using k=2 (sing kmeans function in cluster package libraries)
library(cluster)
clustering=kmeans(t(golub), centers = 2, nstart = 38)

#comparison between already given clusters in golub.cl vector and the clusters made from kmeans clustering
table(clustering$cluster, golub.cl)
#drawing of slihouette plot to visualize clusters after doing kmeans clustering
#calculation of all the pairwise distances between observations in the dataset (daisy function in cluster package libraries)
dissE <- daisy(t(golub))
#Convert to squared eucledian distance
dE2   <- dissE^2
#plotting of silhouette plot
sil1 <- silhouette(clustering$cluster, dE2)
plot(sil1)
#plotting of kmeans clusters
library('factoextra')
fviz_cluster(clustering, data = t(golub))

#Plotting already available class data
dissE <- daisy(t(golub))
#Convert to squared eucledian distance
dE2   <- dissE^2
#plotting of silhouette plot
sil2 <- silhouette(golub.cl, dE2)
plot(sil2)
#####################
#5.
#constructing heat map to the dataset (using heatmap/2 function in gplots package libraries )
library(gplots)
#loading the gene names as row names
golub.gnames=as.matrix(golub.gnames)
row.names(golub)=golub.gnames[,3]
#Heatmap 
heatmap.2(golub, Rowv=as.dendrogram(hr_pearson), Colv=as.dendrogram(hc_pearson), col=wheel(38, sat = 1),
          scale="column", density.info="none", trace="none")
##############


                    