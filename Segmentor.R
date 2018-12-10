library(Segmentor3IsBack)
library(neuroblastoma)
library(PeakSegDP)

#Getting the Data
data(chr11ChIPseq,package="PeakSegDP")
data(neuroblastoma,package="neuroblastoma")
x=chr11ChIPseq$regions$chromStart
y=neuroblastoma$profiles$logratio[1:20]

#Performing Segmentation
x_seg=Segmentor(x,model=1,Kmax=15)
y_seg=Segmentor(y,model=2,Kmax=5)
z_seg=Segmentor(c(1,2,2,1),model=1,Kmax=4)

#Getting Change Points
x_seg@breaks
y_seg@breaks
z_seg@breaks

#Plotting for case of 6 segments in X
plot(x)
abline(v=2,col="red")
abline(v=3,col="red")
abline(v=4,col="red")
abline(v=7,col="red")
abline(v=8,col="red")
abline(v=16,col="red")
