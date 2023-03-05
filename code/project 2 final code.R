install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("pheatmap")
#importing data from ncbi
library(GEOquery)
gse <- getGEO("GSE1739", GSEMatrix = TRUE ,AnnotGPL=TRUE)
gse <- getGEO(filename="GSE1739_series_matrix.txt.gz")
show(gse)
gse<-gse[[1]]
#print 20 gene expression matrix
exprs(gse[1:20,]) 
#grouping data into two samples
normal<-exprs(gse[1:20,1:4])
patient<-exprs(gse[1:20,5:14])
print(normal) 
print(patient)
#applying T-Test
t.test(normal,patient)
#boxplot
boxplot(normal ,notch=F,xlab="samples" , ylab="value", boxwex=0.6,main="GSE1739-normal samples"
        ,outline=FALSE,las=2,col=c("blue"))
boxplot(patient ,notch=F,xlab="samples" ,ylab="value", boxwex=0.6,main="GSE1739-patient samples",
        outline=FALSE,las=2,col=c("red"))
#boxplot for final results
data<-data.frame(normal,patient)
names(data)<-c("control61","control62","control63","control64","treat65","treat66"
  ,"treat67","treat68","treat69","treat70","treat71","treat72","treat73","treat74")
data
a<-c("blue","blue","blue","blue")
b<-c("red","red","red","red","red","red","red","red","red","red")
boxplot(data,col=c(a,b)) 
#scatterplot
plot(exprs(gse[1:20,]))
#heatmap
heatmap(exprs(gse[1:20,]))



