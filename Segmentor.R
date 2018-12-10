library(Segmentor3IsBack)
library(neuroblastoma)
library(PeakSegDP)

#Getting the Data
data(chr11ChIPseq,package="PeakSegDP")
data(neuroblastoma,package="neuroblastoma")
x=chr11ChIPseq$regions$chromStart
y=neuroblastoma$profiles$logratio[1:20]
z=c(1,2,2,1)

#Performing Segmentation
x_seg=Segmentor(x,model=1,Kmax=15)
y_seg=Segmentor(y,model=2,Kmax=5)
z_seg=Segmentor(z,model=1,Kmax=4)

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

#Plotting for case of 5 segments in Y
plot(y)
abline(v=1,col="blue")
abline(v=7,col="blue")
abline(v=8,col="blue")
abline(v=11,col="blue")
abline(v=20,col="blue")

#Plotting for case of 3 segments in Z (c(1,2,2,1))
plot(z)
abline(v=1,col="green")
abline(v=3,col="green")
abline(v=4,col="green")
