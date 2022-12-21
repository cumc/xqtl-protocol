#based on MRAID, science advances,2022 https://xzlab.org/papers/2022_Yuan_etal_SA.pdf 

setwd("/disks/DATATMP/PhD_AlexandrePelletier/postdoc/xqtl-pipeline/")
# install.packages("devtools")
library(devtools)
install_github("yuanzhongshang/MRAID")

#install.packages("data.table")

out<-"outputs/01-MRAID_test"
dir.create(out,recursive = T)


library(MRAID)
library(data.table)

#Summary data with 100 candidate correlated SNPs, n1=30000, n2=30000, 
data_path<-"code/mendelian_randomization/MRAID/example/" #clone the github repository of MRAID (https://github.com/yuanzhongshang/MRAID)

#load the Zscore vector for the exposure
x<-fread(file.path(data_path,"Zscore_1.txt"))
Zscore_1<-x$x

#load the Zscore vector for the outcome
y<-fread(file.path(data_path,"Zscore_2.txt"))
Zscore_2<-y$x


#load the LD matrix in exposure GWAS data 
zx<-fread(file.path(data_path,"Sigma1sin.txt"))
Sigma1sin<-as.matrix(data.frame(zx,row.names = "V1"))

#load the LD matrix in outcome GWAS data 
zy<-fread(file.path(data_path,"Sigma2sin.txt"))
Sigma2sin<-as.matrix(data.frame(zy,row.names = "V1"))

#load the sample size
samplen1=30000
samplen2=30000

#run MRAID 
?MRAID
result<-MRAID(Zscore_1, Zscore_2, Sigma1sin, Sigma2sin,
              samplen1, samplen2, 
              Gibbsnumber=1000,burninproportion=0.2,
              pi_beta_shape=0.5,
pi_beta_scale=4.5,pi_c_shape=0.5,pi_c_scale=9.5,pi_1_shape=0.5,pi_1_scale=1.5,pi_0_shape=0.05,pi_0_scale=9.95,option=1)

result
# $causal_effect
# [1] 0.07167195
# 
# $causal_pvalue
# [1] 0.2958619
# 
# $correlated_pleiotropy_effect
# [1] -0.00347539
# 
# $sigmabeta
# [1] 0.01766178
# 
# $sigmaeta
# [1] 0.01002248
# 
# $sigma_error_1
# [1] 0.9914672
# 
# $sigma_error_2
# [1] 1.000612



