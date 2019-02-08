## 1. Easy Test

## Required Packages
The following packages are required to perform the segmentation.

```{r, message=FALSE}
#Load packages
if(!require(Segmentor3IsBack)){
  install.packages(Segmentor3IsBack)
  library(Segmentor3IsBack)
}

if(!require(neuroblastoma)){
  install.packages("neuroblastoma")
  library(neuroblastoma)
}

if(!require(PeakSegDP)){
  install.packages("PeakSegDP")
  library(PeakSegDP)
}
```

## Extracting Data from desired packages and generating artificial datasets
```{r, message=FALSE}
#Getting the Data
data(chr11ChIPseq,package="PeakSegDP")
data(neuroblastoma,package="neuroblastoma")
appended_gaussian=c(rpois(200,4),rpois(200,1),rpois(200,2.2))  #Artificial Dataset with 3 different normal distributions of length 200 appended together
PeakSegDP_dataset=chr11ChIPseq$regions$chromStart
neuroblastoma_dataset=neuroblastoma$profiles$logratio[1:20]
edge_case=c(1,2,2,1)
```
![alt tag](https://user-images.githubusercontent.com/37847118/51668016-41998f80-1fe7-11e9-9e6f-c242201aebf0.png)
```{r, message=FALSE}
#Generating Jumping Average Artificial Dataset
mu=function(n){
  if(n==0){return(0)}
  else{
    return(mu(n-1)+(n/16)) }
}
muvect=vector()
for(i in 1:5){
  muvect[i]=mu(i) }
jumping_average=vector()
jumping_average[1]=0
jumping_average[2]=0
for(i in 3:23){
  jumping_average[i]=0.6*jumping_average[i-1]-0.5*jumping_average[i-2]+rnorm(200,muvect[1],1.5)[i] }
for(i in 24:44){
  jumping_average[i]=0.6*jumping_average[i-1]-0.5*jumping_average[i-2]+rnorm(200,muvect[2],1.5)[i] }
for(i in 45:65){
  jumping_average[i]=0.6*jumping_average[i-1]-0.5*jumping_average[i-2]+rnorm(200,muvect[3],1.5)[i] }
for(i in 66:86){
  jumping_average[i]=0.6*jumping_average[i-1]-0.5*jumping_average[i-2]+rnorm(200,muvect[4],1.5)[i] }
for(i in 87:107){
  jumping_average[i]=0.6*jumping_average[i-1]-0.5*jumping_average[i-2]+rnorm(200,muvect[5],1.5)[i] }
jumping_average=jumping_average[3:107]
```
![alt tag](https://user-images.githubusercontent.com/37847118/51668283-de5c2d00-1fe7-11e9-86ae-8600b13e75ad.png)
```{r, message=FALSE}
#Generating Moving Frequency Dataset
o=vector()
o[1]=1
for(i in 2:10){
  o[i]=o[i-1]*log(exp(1)+(i/2))
}
jumping_frequency=vector()
for(i in 1:50){
  jumping_frequency[i]=sin(o[1]*i)+rnorm(5000,0,0.8)[i]
}
for(i in 51:100){
  jumping_frequency[i]=sin(o[2]*i)+rnorm(5000,0,0.8)[i]
}
for(i in 101:150){
  jumping_frequency[i]=sin(o[3]*i)+rnorm(5000,0,0.8)[i]
}
for(i in 151:200){
  jumping_frequency[i]=sin(o[4]*i)+rnorm(5000,0,0.8)[i]
}
for(i in 201:250){
  jumping_frequency[i]=sin(o[5]*i)+rnorm(5000,0,0.8)[i]
}
```
## Applying the Segmentor
```{r, message=FALSE}
appended_gaussian_seg=Segmentor(appended_gaussian,model=1,Kmax=5)
PeakSegDP_dataset_seg=Segmentor(PeakSegDP_dataset,model=1,Kmax=15)
neuroblastoma_dataset_seg=Segmentor(neuroblastoma_dataset,model=2,Kmax=5)
edge_case_seg=Segmentor(edge_case,model=1,Kmax=4)
jumping_average_seg=Segmentor(jumping_average,model = 2,Kmax = 6)
jumping_frequency_seg=Segmentor(jumping_frequency,Kmax = 6,model = 2)
```
## Plotting
The vertical lines in the plots correspond to change-points.

```{r, message=FALSE}
#appended_gaussian
AGchoose=SelectModel(appended_gaussian_seg,penalty = "oracle")
plot(appended_gaussian,col='red')
abline(v=getBreaks(appended_gaussian_seg)[AGchoose, 1:AGchoose],col='blue')
```
![alt tag](https://user-images.githubusercontent.com/37847118/52477751-0ee7bd80-2bc9-11e9-8515-158116cd0725.png)


The algorithm gives the predicted set of change-points, i.e c(200,400,600).

```{r, message=FALSE}
#PeakSegDP_dataset
PSDchoose=SelectModel(PeakSegDP_dataset_seg,penalty = "oracle")
plot(PeakSegDP_dataset,col='dark red')
abline(v=getBreaks(PeakSegDP_dataset_seg)[PSDchoose, 1:PSDchoose],col='blue')
```
![alt tag](https://user-images.githubusercontent.com/37847118/52477745-0e4f2700-2bc9-11e9-8391-cf3a1e512586.png)
```{r, message=FALSE}
#Plotting for case of 5 segments in neuroblastoma_dataset
plot(neuroblastoma_dataset)
abline(v=1,col="blue")
abline(v=7,col="blue")
abline(v=8,col="blue")
abline(v=11,col="blue")
abline(v=20,col="blue")
```
![alt tag](https://user-images.githubusercontent.com/37847118/49724993-b3648680-fc90-11e8-9729-41ccb76dc9fa.png)
```{r, message=FALSE}
#edge_case
ECchoose=SelectModel(edge_case_seg,penalty = "oracle")
plot(edge_case,col="red")
abline(v=getBreaks(edge_case_seg)[ECchoose, 1:ECchoose],col='blue')
```
![alt tag](https://user-images.githubusercontent.com/37847118/52477749-0ee7bd80-2bc9-11e9-874f-6072d3572183.png)

```{r, message=FALSE}
#Jumping_average
JAchoose=SelectModel(jumping_average_seg,penalty = "oracle")
plot(jumping_average,col="red")
abline(v=getBreaks(jumping_average_seg)[JAchoose, 1:JAchoose],col='blue')
```
![alt tag](https://user-images.githubusercontent.com/37847118/52477748-0e4f2700-2bc9-11e9-9293-f6e989a77445.png)


It is clearly observed that the algorithm fails to produce a valid set of change-points for the case of Jumping average artificial dataset.

```{r, message=FALSE}
#Changing Frequency
JFchoose=SelectModel(jumping_frequency_seg,penalty = "oracle")
plot(jumping_frequency,col="red")
abline(v=getBreaks(jumping_frequency_seg)[JFchoose, 1:JFchoose],col='blue')
```
![alt tag](https://user-images.githubusercontent.com/37847118/52477747-0e4f2700-2bc9-11e9-8789-0c26e0b5f928.png)


The algorithm does not bring expected change-points.
