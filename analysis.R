
rm(list=ls())
DSJ40<-read.delim("DSJ40.txt",header=T)
DataSHQ = DSJ40[,18:29]
library(Scale)
names(DataSHQ)<-c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10","Q11","Q12")
DataSHQ<-as.vector(DataSHQ)

#KMO检验
  library(MASS)
  X<- cor(as.matrix(DataSHQ))
  iX<- ginv(X)
  S2<- diag(diag((iX^-1)))
  AIS<- S2%*%iX%*%S2              # anti-image covariance matrix                               
  Dai<- sqrt(diag(diag(AIS)))
  AIR <- ginv(Dai)%*%AIS%*%ginv(Dai)        # anti-image correlation matrix
  a <- apply((AIR - diag(diag(AIR)))^2, 2, sum)
  AA <- sum(a)
  b <- apply((X - diag(nrow(X)))^2, 2, sum)
  BB <- sum(b)
  kmo <- BB/(AA+BB)         # overall KMO statistic
  kmo

#bartlett球形检验
library(psych)
cortest.bartlett(cor(DataSHQ), n=106)

#因子分析，方差最大，因子载荷矩阵，
covariances<-cov(DataSHQ)
correlations<-cov2cor(covariances)
correlations
library(psych)
fa.parallel(correlations,n.obs=106,fa="both",n.iter=100,main="scree plots with parallel analysis")
fa.varimax<-fa(correlations,nfactors=2,rotate="varimax",fm="pa")
fa.varimax$loadings

factor.plot(fa.varimax,labels=rownames(fa.varimax$loadings))  ##正交结果图

fa<-factanal(~.,factors=2,data=DataSHQ,scores="regression")  ##回归法计算因子得分
scs<-fa$scores

##将S表格导出到excel
#setwd("D:\研\李正帮老师\数据文章Code")
#write.table(scs,"scs.csv",sep=",")

#正常人与非正常得分比较  t检验
#SHQ总分比较
SDataSHQ<-apply(DataSHQ,1,sum)
norm<-SDataSHQ[58:106]
dnorm<-SDataSHQ[-(58:106)]
outcome<-data.frame(
x=c(norm,dnorm),
a=factor(rep(1:2,c(49,57)))
)
t.test(x~a,alternative = "two.sided", var.equal = FALSE,data=outcome)

#第一因子比较
F1DataSHQ<-scs[,1]
norm1<-F1DataSHQ[58:106]
dnorm1<-F1DataSHQ[-(58:106)]
outcome1<-data.frame(
x1=c(norm1,dnorm1),
a1=factor(rep(1:2,c(49,57)))
)
t.test(x1~a1,alternative = "two.sided", var.equal = FALSE,data=outcome1)

#第二因子比较
F2DataSHQ<-scs[,2]
norm2<-F2DataSHQ[58:106]
dnorm2<-F2DataSHQ[-(58:106)]
outcome2<-data.frame(
x2=c(norm2,dnorm2),
a2=factor(rep(1:2,c(49,57)))
)
t.test(x2~a2,alternative = "two.sided", var.equal = FALSE,data=outcome2)


#SHQ得分与声源定位的相关性(耳聋55例)
Data<-DSJ40[1:57,c(18:29,31)]
dataSHQ<-Data[c(-7,-36),] 
SHQ<-apply(dataSHQ[,-13],1,sum)
RMS<-dataSHQ[,13]
F1<-scs[c(1:6,8:35,37:57),1]
F2<-scs[c(1:6,8:35,37:57),2]
newdata<-data.frame(F1,F2,RMS,SHQ)
##将S表格导出到excel
#setwd("D:\研\李正帮老师\数据文章Code")
#write.table(newdata,"newdata.csv",sep=",")

newdata1<-read.csv("newdata1.csv",header=T)
covariances<-cov(newdata1)
correlations<-cov2cor(covariances)   
correlations
cor.test(newdata1$RMS,newdata1$SHQ,alternative="two.side",method="pearson")

 
#人工耳蜗植入前后的SHQ得分比较
#SHQ总分比较
before<-SDataSHQ[19:57]
after<-SDataSHQ[1:18]
compare<-data.frame(
x3=c(before,after),
a3=factor(rep(1:2,c(39,18)))
)
t.test(x3~a3,alternative = "two.sided", var.equal = FALSE,data=compare)

#第一因子比较(与总分比较结果一样)
#before1<-scs[19:57,1]
#after1<-scs[1:18,1]
#compare1<-data.frame(
#x4=c(before,after),
#a4=factor(rep(1:2,c(39,18)))
#)
#t.test(x4~a4,alternative = "two.sided", var.equal = FALSE,data=compare1)





