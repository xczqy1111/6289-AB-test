set.seed(11235)
library(carat)
##1.data pre-processing
df<- read.csv("WA_Marketing-Campaign+2.csv")
df1<-df[,c("MarketSize","AgeOfStore","Promotion","week","SalesInThousands")]

table(df1$MarketSize)
table(df1$Promotion)
table(df1$week)

levels(df1$MarketSize)<-list("0"="Small","1"="Medium","2"="Large")
levels(df1$MarketSize)

#summary statistics of ages of stores
summary(df1$AgeOfStore)

#discretize age of stores into categorical variables, by 1st quantile, median, and 3rd quantile
for (i in 1:nrow(df1)) {
  if (df1$AgeOfStore[i]<=4) {
    df1$age.store.ind[i]<-0
  }
  else if (df1$AgeOfStore[i]<=7 & df1$AgeOfStore[i]>4) {
    df1$age.store.ind[i]<-1
  }
  else if (df1$AgeOfStore[i]<=12 & df1$AgeOfStore[i]>7) {
    df1$age.store.ind[i]<-2
  }
  else if (df1$AgeOfStore[i]>12) {
    df1$age.store.ind[i]<-3
  }
}

#discretize sales into binary outcome
summary(df1$SalesInThousands)
df1$sales.ind<-ifelse(df1$SalesInThousands>=50.2,1,0)

#mu for two treatments
mu1<-mean(df1[which(df1$Promotion==1),5])
mu2<-mean(df1[which(df1$Promotion==2),5])

#success rate for two treatments
p1<-mean(df1[which(df1$Promotion==1),7])
p2<-mean(df1[which(df1$Promotion==2),7])

#sample size for two treatments
n1<-nrow(df1[which(df1$Promotion==1),])
n2<-nrow(df1[which(df1$Promotion==2),])

df2<-df1[which(df1$Promotion==1|df1$Promotion==2),]

##2.redesign treatment assignment
##2.1response-adaptive design

rep<-1000
diff<-0.2
#power function for response-adaptive design
power1<-function(xa,xb,na,nb){
  pa<-xa/na
  pb<-xb/nb
  p<-(xa+xb)/(na+nb)
  vpool<-p*(1-p)
  t<-(pa-pb-diff)/sqrt(vpool*(1/na+1/nb))
  if(t>qnorm(0.95)){
    i=1
  }
  else if (t<=qnorm(0.95)){
    i=0
  }
  return(i)
}

n<-n1+n2

#(i) complete randomization 
na_r_matrix<-NULL
i_r_matrix<-NULL
for(i in 1:rep){
  ya_r_matrix<-NULL
  yb_r_matrix<-NULL
  sa_r_matrix<-NULL
  for(j in 1:n){
    sa_r<-rbinom(1,1,0.5)
    if(sa_r==1){
      ya_r<-rbinom(1,1,p1)
      ya_r_matrix<-c(ya_r,ya_r_matrix)
    }
    else if(sa_r==0){
      yb_r<-rbinom(1,1,p2)
      yb_r_matrix<-c(yb_r,yb_r_matrix)
    }
    sa_r_matrix<-c(sa_r,sa_r_matrix)
  }
  na_r<-sum(sa_r_matrix)
  nb_r<-n-sum(sa_r_matrix)
  xa_r<-sum(ya_r_matrix)
  xb_r<-sum(yb_r_matrix)
  i_r<-power1(xa_r,xb_r,na_r,nb_r)
  na_r_matrix<-c(na_r,na_r_matrix)
  i_r_matrix<-c(i_r,i_r_matrix)
}
#power for complete randomization
power_r<-mean(i_r_matrix)
#mean of na for complete randomization
mean_r<-mean(na_r_matrix)
#sd of na for complete randomization
sd_r<-sd(na_r_matrix)

#(ii) play-the-winner rule starts with (1,1)
na_rpw_matrix<-NULL
i_rpw_matrix<-NULL
for(i in 1:rep){
  ya_rpw_matrix<-NULL
  yb_rpw_matrix<-NULL
  sa_rpw_matrix<-NULL
  for(j in 1:n){
    sa_rpw<-rbinom(1,1,(1+sum(ya_rpw_matrix)+(length(yb_rpw_matrix)-sum(yb_rpw_matrix)))/(2+length(ya_rpw_matrix)+length(yb_rpw_matrix)))
    if(sa_rpw==1){
      ya_rpw<-rbinom(1,1,p1)
      ya_rpw_matrix<-c(ya_rpw,ya_rpw_matrix)
    }
    else if(sa_rpw==0){
      yb_rpw<-rbinom(1,1,p2)
      yb_rpw_matrix<-c(yb_rpw,yb_rpw_matrix)
    }
    sa_rpw_matrix<-c(sa_rpw,sa_rpw_matrix)
  }
  na_rpw<-sum(sa_rpw_matrix)
  nb_rpw<-n-sum(sa_rpw_matrix)
  xa_rpw<-sum(ya_rpw_matrix)
  xb_rpw<-sum(yb_rpw_matrix)
  i_rpw<-power1(xa_rpw,xb_rpw,na_rpw,nb_rpw)
  na_rpw_matrix<-c(na_rpw,na_rpw_matrix)
  i_rpw_matrix<-c(i_rpw,i_rpw_matrix)
}
#power for RPW(1,1)
power_rpw<-mean(i_rpw_matrix)
#mean of na for RPW(1,1)
mean_rpw<-mean(na_rpw_matrix)
#sd of na for RPW(1,1)
sd_rpw<-sd(na_rpw_matrix)

#(iii) DBCD gamma=0
n0=5

g<-function(rho,x,gamma){
  (rho*(rho/x)^gamma)/((rho*(rho/x)^gamma)+((1-rho)*((1-rho)/(1-x))^gamma))
}

na_dbcd_r0_matrix<-NULL
i_dbcd_r0_matrix<-NULL
for(i in 1:rep){
  ya_dbcd_r0_matrix<-NULL
  yb_dbcd_r0_matrix<-NULL
  sa_dbcd_r0_matrix<-NULL
  xa_ini<-sum(rbinom(n0,1,p1))
  xb_ini<-sum(rbinom(n0,1,p2))
  for(j in 1:(n-2*n0)){
    pa_hat<-(xa_ini+sum(ya_dbcd_r0_matrix)+0.5)/(n0+length(ya_dbcd_r0_matrix)+1)
    pb_hat<-(xb_ini+sum(yb_dbcd_r0_matrix)+0.5)/(n0+length(yb_dbcd_r0_matrix)+1)
    qa_hat<-1-pa_hat
    qb_hat<-1-pb_hat
    rho<-qb_hat/(qa_hat+qb_hat)
    x<-(n0+length(ya_dbcd_r0_matrix))/(n0+length(ya_dbcd_r0_matrix)+n0+length(yb_dbcd_r0_matrix))
    sa_dbcd_r0<-rbinom(1,1,g(rho,x,gamma=0))
    if(sa_dbcd_r0==1){
      ya_dbcd_r0<-rbinom(1,1,p1)
      ya_dbcd_r0_matrix<-c(ya_dbcd_r0,ya_dbcd_r0_matrix)
    }
    else if(sa_dbcd_r0==0){
      yb_dbcd_r0<-rbinom(1,1,p2)
      yb_dbcd_r0_matrix<-c(yb_dbcd_r0,yb_dbcd_r0_matrix)
    }
    sa_dbcd_r0_matrix<-c(sa_dbcd_r0,sa_dbcd_r0_matrix)
  }
  na_dbcd_r0<-sum(sa_dbcd_r0_matrix)+n0
  nb_dbcd_r0<-n-(sum(sa_dbcd_r0_matrix)+n0)
  xa_dbcd_r0<-sum(ya_dbcd_r0_matrix)+xa_ini
  xb_dbcd_r0<-sum(yb_dbcd_r0_matrix)+xb_ini
  i_dbcd_r0<-power1(xa_dbcd_r0,xb_dbcd_r0,na_dbcd_r0,nb_dbcd_r0)
  na_dbcd_r0_matrix<-c(na_dbcd_r0,na_dbcd_r0_matrix)
  i_dbcd_r0_matrix<-c(i_dbcd_r0,i_dbcd_r0_matrix)
}
#power for DBCD(gamma=0)
power_dbcd_r0<-mean(i_dbcd_r0_matrix)
#mean of na for DBCD(gamma=0)
mean_dbcd_r0<-mean(na_dbcd_r0_matrix)
#sd of na for DBCD(gamma=0)
sd_dbcd_r0<-sd(na_dbcd_r0_matrix)

#(iv) DBCD gamma=2
n0=5

g<-function(rho,x,gamma){
  (rho*(rho/x)^gamma)/((rho*(rho/x)^gamma)+((1-rho)*((1-rho)/(1-x))^gamma))
}

na_dbcd_r2_matrix<-NULL
i_dbcd_r2_matrix<-NULL
for(i in 1:rep){
  ya_dbcd_r2_matrix<-NULL
  yb_dbcd_r2_matrix<-NULL
  sa_dbcd_r2_matrix<-NULL
  xa_ini<-sum(rbinom(n0,1,p1))
  xb_ini<-sum(rbinom(n0,1,p2))
  for(j in 1:(n-2*n0)){
    pa_hat<-(xa_ini+sum(ya_dbcd_r2_matrix)+0.5)/(n0+length(ya_dbcd_r2_matrix)+1)
    pb_hat<-(xb_ini+sum(yb_dbcd_r2_matrix)+0.5)/(n0+length(yb_dbcd_r2_matrix)+1)
    qa_hat<-1-pa_hat
    qb_hat<-1-pb_hat
    rho<-qb_hat/(qa_hat+qb_hat)
    x<-(n0+length(ya_dbcd_r2_matrix))/(n0+length(ya_dbcd_r2_matrix)+n0+length(yb_dbcd_r2_matrix))
    sa_dbcd_r2<-rbinom(1,1,g(rho,x,gamma=2))
    if(sa_dbcd_r2==1){
      ya_dbcd_r2<-rbinom(1,1,p1)
      ya_dbcd_r2_matrix<-c(ya_dbcd_r2,ya_dbcd_r2_matrix)
    }
    else if(sa_dbcd_r2==0){
      yb_dbcd_r2<-rbinom(1,1,p2)
      yb_dbcd_r2_matrix<-c(yb_dbcd_r2,yb_dbcd_r2_matrix)
    }
    sa_dbcd_r2_matrix<-c(sa_dbcd_r2,sa_dbcd_r2_matrix)
  }
  na_dbcd_r2<-sum(sa_dbcd_r2_matrix)+n0
  nb_dbcd_r2<-n-(sum(sa_dbcd_r2_matrix)+n0)
  xa_dbcd_r2<-sum(ya_dbcd_r2_matrix)+xa_ini
  xb_dbcd_r2<-sum(yb_dbcd_r2_matrix)+xb_ini
  i_dbcd_r2<-power1(xa_dbcd_r2,xb_dbcd_r2,na_dbcd_r2,nb_dbcd_r2)
  na_dbcd_r2_matrix<-c(na_dbcd_r2,na_dbcd_r2_matrix)
  i_dbcd_r2_matrix<-c(i_dbcd_r2,i_dbcd_r2_matrix)
}
#power for DBCD(gamma=2)
power_dbcd_r2<-mean(i_dbcd_r2_matrix)
#mean of na for DBCD(gamma=2)
mean_dbcd_r2<-mean(na_dbcd_r2_matrix)
#sd of na for DBCD(gamma=2)
sd_dbcd_r2<-sd(na_dbcd_r2_matrix)

#(v) erade alpha=0.5
n0=5

na_erade_matrix<-NULL
i_erade_matrix<-NULL
for(i in 1:rep){
  ya_erade_matrix<-NULL
  yb_erade_matrix<-NULL
  sa_erade_matrix<-NULL
  xa_ini<-sum(rbinom(n0,1,p1))
  xb_ini<-sum(rbinom(n0,1,p2))
  for(j in 1:(n-2*n0)){
    pa_hat<-(xa_ini+sum(ya_erade_matrix)+0.5)/(n0+length(ya_erade_matrix)+1)
    pb_hat<-(xb_ini+sum(yb_erade_matrix)+0.5)/(n0+length(yb_erade_matrix)+1)
    qa_hat<-1-pa_hat
    qb_hat<-1-pb_hat
    rho<-qb_hat/(qa_hat+qb_hat)
    prop_a<-(length(ya_erade_matrix)+n0)/(length(ya_erade_matrix)+n0+length(yb_erade_matrix)+n0)
    if (prop_a==rho){
      sa_erade<-rbinom(1,1,0.5)
    }
    else if (prop_a>rho){
      sa_erade<-rbinom(1,1,0.5*rho)
    }    
    else if (prop_a<rho){
      sa_erade<-rbinom(1,1,(1-0.5*(1-rho)))
    }     
    if(sa_erade==1){
      ya_erade<-rbinom(1,1,p1)
      ya_erade_matrix<-c(ya_erade,ya_erade_matrix)
    }
    else if(sa_erade==0){
      yb_erade<-rbinom(1,1,p2)
      yb_erade_matrix<-c(yb_erade,yb_erade_matrix)
    }
    sa_erade_matrix<-c(sa_erade,sa_erade_matrix)
  }
  na_erade<-sum(sa_erade_matrix)+n0
  nb_erade<-n-(sum(sa_erade_matrix)+n0)
  xa_erade<-sum(ya_erade_matrix)+xa_ini
  xb_erade<-sum(yb_erade_matrix)+xb_ini
  i_erade<-power1(xa_erade,xb_erade,na_erade,nb_erade)
  na_erade_matrix<-c(na_erade,na_erade_matrix)
  i_erade_matrix<-c(i_erade,i_erade_matrix)
}
#power for erade
power_erade<-mean(i_erade_matrix)
#mean of na for erade
mean_erade<-mean(na_erade_matrix)
#sd of na for erade
sd_erade<-sd(na_erade_matrix)

#(vi) DBCD gamma=2 with sqrt rho
n0=5

g<-function(rho,x,gamma){
  (rho*(rho/x)^gamma)/((rho*(rho/x)^gamma)+((1-rho)*((1-rho)/(1-x))^gamma))
}

na_dbcd_sqrt_r2_matrix<-NULL
i_dbcd_sqrt_r2_matrix<-NULL
for(i in 1:rep){
  ya_dbcd_sqrt_r2_matrix<-NULL
  yb_dbcd_sqrt_r2_matrix<-NULL
  sa_dbcd_sqrt_r2_matrix<-NULL
  xa_ini<-sum(rbinom(n0,1,p1))
  xb_ini<-sum(rbinom(n0,1,p2))
  for(j in 1:(n-2*n0)){
    pa_hat<-(xa_ini+sum(ya_dbcd_sqrt_r2_matrix)+0.5)/(n0+length(ya_dbcd_sqrt_r2_matrix)+1)
    pb_hat<-(xb_ini+sum(yb_dbcd_sqrt_r2_matrix)+0.5)/(n0+length(yb_dbcd_sqrt_r2_matrix)+1)
    qa_hat<-1-pa_hat
    qb_hat<-1-pb_hat
    rho<-sqrt(pa_hat)/(sqrt(pa_hat)+sqrt(pb_hat))
    x<-(n0+length(ya_dbcd_sqrt_r2_matrix))/(n0+length(ya_dbcd_sqrt_r2_matrix)+n0+length(yb_dbcd_sqrt_r2_matrix))
    sa_dbcd_sqrt_r2<-rbinom(1,1,g(rho,x,gamma=2))
    if(sa_dbcd_sqrt_r2==1){
      ya_dbcd_sqrt_r2<-rbinom(1,1,p1)
      ya_dbcd_sqrt_r2_matrix<-c(ya_dbcd_sqrt_r2,ya_dbcd_sqrt_r2_matrix)
    }
    else if(sa_dbcd_sqrt_r2==0){
      yb_dbcd_sqrt_r2<-rbinom(1,1,p2)
      yb_dbcd_sqrt_r2_matrix<-c(yb_dbcd_sqrt_r2,yb_dbcd_sqrt_r2_matrix)
    }
    sa_dbcd_sqrt_r2_matrix<-c(sa_dbcd_sqrt_r2,sa_dbcd_sqrt_r2_matrix)
  }
  na_dbcd_sqrt_r2<-sum(sa_dbcd_sqrt_r2_matrix)+n0
  nb_dbcd_sqrt_r2<-n-(sum(sa_dbcd_sqrt_r2_matrix)+n0)
  xa_dbcd_sqrt_r2<-sum(ya_dbcd_sqrt_r2_matrix)+xa_ini
  xb_dbcd_sqrt_r2<-sum(yb_dbcd_sqrt_r2_matrix)+xb_ini
  i_dbcd_sqrt_r2<-power1(xa_dbcd_sqrt_r2,xb_dbcd_sqrt_r2,na_dbcd_sqrt_r2,nb_dbcd_sqrt_r2)
  na_dbcd_sqrt_r2_matrix<-c(na_dbcd_sqrt_r2,na_dbcd_sqrt_r2_matrix)
  i_dbcd_sqrt_r2_matrix<-c(i_dbcd_sqrt_r2,i_dbcd_sqrt_r2_matrix)
}
#power for DBCD(sqrt_rho gamma=2)
power_dbcd_sqrt_r2<-mean(i_dbcd_sqrt_r2_matrix)
#mean of na for DBCD(sqrt_rho gamma=2)
mean_dbcd_sqrt_r2<-mean(na_dbcd_sqrt_r2_matrix)
#sd of na for DBCD(sqrt_rho gamma=2)
sd_dbcd_sqrt_r2<-sd(na_dbcd_sqrt_r2_matrix)

power<-c(power_r,power_rpw,power_dbcd_r0,power_dbcd_r2,power_erade,power_dbcd_sqrt_r2)
mean<-c(mean_r,mean_rpw,mean_dbcd_r0,mean_dbcd_r2,mean_erade,mean_dbcd_sqrt_r2)
sd<-c(sd_r,sd_rpw,sd_dbcd_r0,sd_dbcd_r2,sd_erade,sd_dbcd_sqrt_r2)

results<-round(rbind(power,mean,sd),4)
colnames(results)<-c("CR","RPW","DBCD_0","DBCD_2","ERADE","DBCD_sqrt_2")
results

##2.2 covariate-adaptive design
strata_ind=matrix(rep(0,144),ncol=3)
k=1

for (i1 in 0:2) {
  for (i2  in 0:3) {
    for (i3 in 1:4) {
      strata_ind[k,]=c(i1,i2,i3)
      k=k+1
    }
  }
}

strata_find=function(dat){
  res=rep(0,nrow(strata_ind))
  n=dim(dat)[2]
  for (i in 1:nrow(strata_ind)) {
    k=0
    for (j in 1:nrow(dat)){
      temp<-(dat[j,]==strata_ind[i,])
      k=k+all(temp)
    }
    res[i]=k
  }
  return(res)
}

#counting function
count0<-function(str){
  temp<-ifelse(str==0,1,0)
  return(temp)
}
count1<-function(str){
  temp<-ifelse(str==1,1,0)
  return(temp)
}

count2<-function(str){
  temp<-ifelse(str==2,1,0)
  return(temp)
}
count3<-function(str){
  temp<-ifelse(str==3,1,0)
  return(temp)
}
count4<-function(str){
  temp<-ifelse(str>=4,1,0)
  return(temp)
}

#imbalance function
imba<-function(data,t){
  ta<-ifelse(t=="A",1,-1)
  data<-rep(1,length(data))
  temp<-abs(sum(ta*data))
  return(temp)
}

df3<-df2[,c("MarketSize","age.store.ind","week")]
#Table 3
table3<-function(n,repli){
  str0_matrix<-NULL
  str1_matrix<-NULL
  str2_matrix<-NULL
  str3_matrix<-NULL
  str4_matrix<-NULL
  for (i in 1:repli){
    data<-df3
    str0<-sum(count0(strata_find(data)))
    str1<-sum(count1(strata_find(data)))
    str2<-sum(count2(strata_find(data)))
    str3<-sum(count3(strata_find(data)))
    str4<-sum(count4(strata_find(data)))
    
    str0_matrix<-c(str0_matrix,str0)
    str1_matrix<-c(str1_matrix,str1)
    str2_matrix<-c(str2_matrix,str2)
    str3_matrix<-c(str3_matrix,str3)
    str4_matrix<-c(str4_matrix,str4)
  }
  return(c(mean(str0_matrix),mean(str1_matrix),mean(str2_matrix),mean(str3_matrix),mean(str4_matrix)))
}

#Table 4
table4<-function(n,repli,method){
  t4<-NULL
  for (i in 1:repli){
    data<-df3
    if (method=="CR"){
      assignments<-rbinom(n,1,0.5)
      assignments<-ifelse(assignments == 1, "A", "B")
      ovimba<-abs(sum(ifelse(assignments=="A",1,-1)))
    }
    else if (method=="StrPBR"){
      assignments<-StrPBR(data,bsize=4)$assignments
      ovimba<-abs(sum(ifelse(assignments=="A",1,-1)))
    }
    else if (method=="PS"){
      assignments<-PocSimMIN(data,c(0.35,0.35,0.3),p=0.75)$assignments
      ovimba<-abs(sum(ifelse(assignments=="A",1,-1)))
    }
    else if (method=="HuHu"){
      assignments<-HuHuCAR(data,c(0.35,0.35,0.1,0.1,0.1),p=0.75)$assignments
      ovimba<-abs(sum(ifelse(assignments=="A",1,-1)))
    }
    t4<-cbind(ovimba,t4)
  }
  mean<-mean(t4)
  median<-median(t4)
  quantile95<-quantile(t4,0.95)
  result<-c(mean,median,quantile95)
  return(result)
}

#Table 5
table5<-function(n,repli,method){
  imba_market_matrix<-NULL
  imba_age_matrix<-NULL
  imba_week_matrix<-NULL
  for (i in 1:repli){
    data<-df3
    if (method=="CR"){
      assignments<-rbinom(n,1,0.5)
      assignments<-ifelse(assignments == 1, "A", "B")
      data1<-data.frame(data,assignments,stringsAsFactors = TRUE)
      
      imba.market.small<-imba(data1[data1$MarketSize==0,1],data1[data1$MarketSize==0,4])
      imba.market.medium<-imba(data1[data1$MarketSize==1,1],data1[data1$MarketSize==1,4])
      imba.market.large<-imba(data1[data1$MarketSize==2,1],data1[data1$MarketSize==2,4])
      
      imba.age0<-imba(data1[data1$age.store.ind==0,2],data1[data1$age.store.ind==0,4])
      imba.age1<-imba(data1[data1$age.store.ind==1,2],data1[data1$age.store.ind==1,4])
      imba.age2<-imba(data1[data1$age.store.ind==2,2],data1[data1$age.store.ind==2,4])
      imba.age3<-imba(data1[data1$age.store.ind==3,2],data1[data1$age.store.ind==3,4])
      
      imba_week1<-imba(data1[data1$week==1,3],data1[data1$week==1,4])
      imba_week2<-imba(data1[data1$week==2,3],data1[data1$week==2,4])
      imba_week3<-imba(data1[data1$week==3,3],data1[data1$week==3,4])
      imba_week4<-imba(data1[data1$week==4,3],data1[data1$week==4,4])
      
      imba_market<-c(imba.market.small,imba.market.medium,imba.market.large)
      imba_age<-c(imba.age0,imba.age1,imba.age2,imba.age3)
      imba_week<-c(imba_week1,imba_week2,imba_week3,imba_week4)
    }
    else if (method=="StrPBR"){
      assignments<-StrPBR(data,bsize=4)$assignments
      data1<-data.frame(data,assignments,stringsAsFactors = TRUE)
      
      imba.market.small<-imba(data1[data1$MarketSize==0,1],data1[data1$MarketSize==0,4])
      imba.market.medium<-imba(data1[data1$MarketSize==1,1],data1[data1$MarketSize==1,4])
      imba.market.large<-imba(data1[data1$MarketSize==2,1],data1[data1$MarketSize==2,4])
      
      imba.age0<-imba(data1[data1$age.store.ind==0,2],data1[data1$age.store.ind==0,4])
      imba.age1<-imba(data1[data1$age.store.ind==1,2],data1[data1$age.store.ind==1,4])
      imba.age2<-imba(data1[data1$age.store.ind==2,2],data1[data1$age.store.ind==2,4])
      imba.age3<-imba(data1[data1$age.store.ind==3,2],data1[data1$age.store.ind==3,4])
      
      imba_week1<-imba(data1[data1$week==1,3],data1[data1$week==1,4])
      imba_week2<-imba(data1[data1$week==2,3],data1[data1$week==2,4])
      imba_week3<-imba(data1[data1$week==3,3],data1[data1$week==3,4])
      imba_week4<-imba(data1[data1$week==4,3],data1[data1$week==4,4])
      
      imba_market<-c(imba.market.small,imba.market.medium,imba.market.large)
      imba_age<-c(imba.age0,imba.age1,imba.age2,imba.age3)
      imba_week<-c(imba_week1,imba_week2,imba_week3,imba_week4)
    }
    else if (method=="PS"){
      assignments<-PocSimMIN(data,c(0.35,0.35,0.3),p=0.75)$assignments
      data1<-data.frame(data,assignments,stringsAsFactors = TRUE)
      
      imba.market.small<-imba(data1[data1$MarketSize==0,1],data1[data1$MarketSize==0,4])
      imba.market.medium<-imba(data1[data1$MarketSize==1,1],data1[data1$MarketSize==1,4])
      imba.market.large<-imba(data1[data1$MarketSize==2,1],data1[data1$MarketSize==2,4])
      
      imba.age0<-imba(data1[data1$age.store.ind==0,2],data1[data1$age.store.ind==0,4])
      imba.age1<-imba(data1[data1$age.store.ind==1,2],data1[data1$age.store.ind==1,4])
      imba.age2<-imba(data1[data1$age.store.ind==2,2],data1[data1$age.store.ind==2,4])
      imba.age3<-imba(data1[data1$age.store.ind==3,2],data1[data1$age.store.ind==3,4])
      
      imba_week1<-imba(data1[data1$week==1,3],data1[data1$week==1,4])
      imba_week2<-imba(data1[data1$week==2,3],data1[data1$week==2,4])
      imba_week3<-imba(data1[data1$week==3,3],data1[data1$week==3,4])
      imba_week4<-imba(data1[data1$week==4,3],data1[data1$week==4,4])
      
      imba_market<-c(imba.market.small,imba.market.medium,imba.market.large)
      imba_age<-c(imba.age0,imba.age1,imba.age2,imba.age3)
      imba_week<-c(imba_week1,imba_week2,imba_week3,imba_week4)
    }
    else if (method=="HuHu"){
      assignments<-HuHuCAR(data,c(0.35,0.35,0.1,0.1,0.1),p=0.75)$assignments
      data1<-data.frame(data,assignments,stringsAsFactors = TRUE)
      
      imba.market.small<-imba(data1[data1$MarketSize==0,1],data1[data1$MarketSize==0,4])
      imba.market.medium<-imba(data1[data1$MarketSize==1,1],data1[data1$MarketSize==1,4])
      imba.market.large<-imba(data1[data1$MarketSize==2,1],data1[data1$MarketSize==2,4])
      
      imba.age0<-imba(data1[data1$age.store.ind==0,2],data1[data1$age.store.ind==0,4])
      imba.age1<-imba(data1[data1$age.store.ind==1,2],data1[data1$age.store.ind==1,4])
      imba.age2<-imba(data1[data1$age.store.ind==2,2],data1[data1$age.store.ind==2,4])
      imba.age3<-imba(data1[data1$age.store.ind==3,2],data1[data1$age.store.ind==3,4])
      
      imba_week1<-imba(data1[data1$week==1,3],data1[data1$week==1,4])
      imba_week2<-imba(data1[data1$week==2,3],data1[data1$week==2,4])
      imba_week3<-imba(data1[data1$week==3,3],data1[data1$week==3,4])
      imba_week4<-imba(data1[data1$week==4,3],data1[data1$week==4,4])
      
      imba_market<-c(imba.market.small,imba.market.medium,imba.market.large)
      imba_age<-c(imba.age0,imba.age1,imba.age2,imba.age3)
      imba_week<-c(imba_week1,imba_week2,imba_week3,imba_week4)
    }
    imba_market_matrix<-rbind(imba_market,imba_market_matrix)
    imba_age_matrix<-rbind(imba_age,imba_age_matrix)
    imba_week_matrix<-rbind(imba_week,imba_week_matrix)
  }
  result<-list(apply(imba_market_matrix,2,mean),apply(imba_age_matrix,2,mean),apply(imba_week_matrix,2,mean))
  return(result)
}

#Table 6
table6<-function(n,repli,method){
  result2_matrix<-NULL
  result3_matrix<-NULL
  for (i in 1:repli){
    data<-df3
    if (method=="CR"){
      assignments<-rbinom(n,1,0.5)
      assignments<-ifelse(assignments == 1, "A", "B")
      data1<-data.frame(data,assignments,stringsAsFactors = TRUE)
      dataA<-data1[data1$assignments=="A",1:3]
      dataB<-data1[data1$assignments=="B",1:3]
      if (sum(which(strata_find(data1[,1:3])==2))!=0){
        index2<-which(strata_find(data1[,1:3])==2)
        diff2<-abs(strata_find(dataA)[index2]-strata_find(dataB)[index2])
      }
      if (sum(which(strata_find(data1[,1:3])==2))==0){
        diff2<-NULL
      }
      if (sum(which(strata_find(data1[,1:3])==3))!=0){
        index3<-which(strata_find(data1[,1:3])==3)
        diff3<-abs(strata_find(dataA)[index3]-strata_find(dataB)[index3])
      }
      if (sum(which(strata_find(data1[,1:3])==3))==0){
        diff3<-NULL
      }
      result2_matrix<-c(result2_matrix,as.vector(diff2))
      result3_matrix<-c(result3_matrix,as.vector(diff3))
    }
    else if (method=="StrPBR"){
      assignments<-StrPBR(data,bsize=4)$assignments
      data1<-data.frame(data,assignments,stringsAsFactors = TRUE)
      dataA<-data1[data1$assignments=="A",1:3]
      dataB<-data1[data1$assignments=="B",1:3]
      if (sum(which(strata_find(data1[,1:3])==2))!=0){
        index2<-which(strata_find(data1[,1:3])==2)
        diff2<-abs(strata_find(dataA)[index2]-strata_find(dataB)[index2])
      }
      if (sum(which(strata_find(data1[,1:3])==2))==0){
        diff2<-NULL
      }
      if (sum(which(strata_find(data1[,1:3])==3))!=0){
        index3<-which(strata_find(data1[,1:3])==3)
        diff3<-abs(strata_find(dataA)[index3]-strata_find(dataB)[index3])
      }
      if (sum(which(strata_find(data1[,1:3])==3))==0){
        diff3<-NULL
      }
      result2_matrix<-c(result2_matrix,as.vector(diff2))
      result3_matrix<-c(result3_matrix,as.vector(diff3))
    }
    else if (method=="PS"){
      assignments<-PocSimMIN(data,c(0.35,0.35,0.3),p=0.75)$assignments
      data1<-data.frame(data,assignments,stringsAsFactors = TRUE)
      dataA<-data1[data1$assignments=="A",1:3]
      dataB<-data1[data1$assignments=="B",1:3]
      if (sum(which(strata_find(data1[,1:3])==2))!=0){
        index2<-which(strata_find(data1[,1:3])==2)
        diff2<-abs(strata_find(dataA)[index2]-strata_find(dataB)[index2])
      }
      if (sum(which(strata_find(data1[,1:3])==2))==0){
        diff2<-NULL
      }
      if (sum(which(strata_find(data1[,1:3])==3))!=0){
        index3<-which(strata_find(data1[,1:3])==3)
        diff3<-abs(strata_find(dataA)[index3]-strata_find(dataB)[index3])
      }
      if (sum(which(strata_find(data1[,1:3])==3))==0){
        diff3<-NULL
      }
      result2_matrix<-c(result2_matrix,as.vector(diff2))
      result3_matrix<-c(result3_matrix,as.vector(diff3))
    }
    else if (method=="HuHu"){
      assignments<-HuHuCAR(data,c(0.35,0.35,0.1,0.1,0.1),p=0.75)$assignments
      data1<-data.frame(data,assignments,stringsAsFactors = TRUE)
      dataA<-data1[data1$assignments=="A",1:3]
      dataB<-data1[data1$assignments=="B",1:3]
      if (sum(which(strata_find(data1[,1:3])==2))!=0){
        index2<-which(strata_find(data1[,1:3])==2)
        diff2<-abs(strata_find(dataA)[index2]-strata_find(dataB)[index2])
      }
      if (sum(which(strata_find(data1[,1:3])==2))==0){
        diff2<-NULL
      }
      if (sum(which(strata_find(data1[,1:3])==3))!=0){
        index3<-which(strata_find(data1[,1:3])==3)
        diff3<-abs(strata_find(dataA)[index3]-strata_find(dataB)[index3])
      }
      if (sum(which(strata_find(data1[,1:3])==3))==0){
        diff3<-NULL
      }
      result2_matrix<-c(result2_matrix,as.vector(diff2))
      result3_matrix<-c(result3_matrix,as.vector(diff3))
    }
  }
  prob2<-length(which(result2_matrix==2))/length(result2_matrix)
  prob0<-length(which(result2_matrix==0))/length(result2_matrix)
  prob3<-length(which(result3_matrix==3))/length(result3_matrix)
  prob1<-length(which(result3_matrix==1))/length(result3_matrix)
  mean2<-mean(result2_matrix)
  mean3<-mean(result3_matrix)
  return(list(c(prob2,prob0),c(prob3,prob1),c(mean2,mean3)))
}  

diff<-mu1-mu2

y0<-df2$SalesInThousands
t0<--(df2$Promotion-2)
z1<-df3$MarketSize
d = data.frame(z1)
d$z1[d$z1=='Small']=1
d$z1[d$z1=='Medium']=2
d$z1[d$z1=='Large']=3
z1 = as.numeric(d$z1) 
z2<-df3$age.store.ind
z3<-df3$week

fit0<-lm(y0~t0+z1+z2)
sigma<-summary(fit0)$sigma
summary(fit0)
b0<-fit0$coefficients[1]
b1<-fit0$coefficients[2]
b2<-fit0$coefficients[3]
b3<-fit0$coefficients[4]

#power for covariate-adaptive design
repli<-1000
power2<-function(n,method,diff){
  dat<-df3
  #complete randomization 
  if (method=="CR"){
    t=rbinom(n,1,0.5)
  }
  else if(method=="StrPBR"){
    #stratified permuted block randomization
    t<-StrPBR(dat,bsize=4)$assignments
    t<-ifelse(t == "A", 1, 0)
  }
  else if(method=="PS"){
    #Pocock and Simon’s method 
    t<-PocSimMIN(dat,weight=c(0.35,0.35,0.3),p=0.75)$assignments
    t<-ifelse(t == "A", 1, 0)
  }
  else if(method=="HH"){
    #Hu and Hu 
    t<-HuHuCAR(dat,c(0.35,0.35,0.1,0.1,0.1),p=0.75)$assignments
    t<-ifelse(t == "A", 1, 0)
  }
  
  y<-b0+diff*t+b2*z1+b3*z2+rnorm(n,0,sigma)
  
  fit1<-lm(y~t)
  fit2<-lm(y~t+z1)
  fit3<-lm(y~t+z2)
  fit4<-lm(y~t+z1+z2)
  var1<-vcov(fit1)
  var2<-vcov(fit2)
  var3<-vcov(fit3)
  var4<-vcov(fit4)
  test1 = fit1$coefficients[2]/sqrt(var1[2,2])
  test2 = fit2$coefficients[2]/sqrt(var2[2,2])
  test3 = fit3$coefficients[2]/sqrt(var3[2,2])
  test4 = fit4$coefficients[2]/sqrt(var4[2,2])
  i_power2 = rep(0,4)
  i_power2[1]<-ifelse(abs(test1) < qnorm(0.975), 0, 1)
  i_power2[2]<-ifelse(abs(test2) < qnorm(0.975), 0, 1)
  i_power2[3]<-ifelse(abs(test3) < qnorm(0.975), 0, 1)
  i_power2[4]<-ifelse(abs(test4) < qnorm(0.975), 0, 1)
  return(i_power2)
}

power2_cr_matrix<-NULL
power2_StrPBR_matrix<-NULL
power2_PS_matrix<-NULL
power2_HH_matrix<-NULL

#diff=0 (Type I error)
for(i in 1:repli){
  power2_cr<-power2(n,method="CR",0*diff)
  power2_StrPBR<-power2(n,method="StrPBR",0*diff)
  power2_PS<-power2(n,method="PS",0*diff)
  power2_HH<-power2(n,method="HH",0*diff)
  power2_cr_matrix<-rbind(power2_cr,power2_cr_matrix)
  power2_StrPBR_matrix<-rbind(power2_StrPBR,power2_StrPBR_matrix)
  power2_PS_matrix<-rbind(power2_PS,power2_PS_matrix)
  power2_HH_matrix<-rbind(power2_HH,power2_HH_matrix)
}

apply(power2_cr_matrix,2,mean)
apply(power2_StrPBR_matrix,2,mean)
apply(power2_PS_matrix,2,mean)
apply(power2_HH_matrix,2,mean)

#diff=0.2*treatment effect
for(i in 1:repli){
  power2_cr<-power2(n,method="CR",0.2*diff)
  power2_StrPBR<-power2(n,method="StrPBR",0.2*diff)
  power2_PS<-power2(n,method="PS",0.2*diff)
  power2_HH<-power2(n,method="HH",0.2*diff)
  power2_cr_matrix<-rbind(power2_cr,power2_cr_matrix)
  power2_StrPBR_matrix<-rbind(power2_StrPBR,power2_StrPBR_matrix)
  power2_PS_matrix<-rbind(power2_PS,power2_PS_matrix)
  power2_HH_matrix<-rbind(power2_HH,power2_HH_matrix)
}

apply(power2_cr_matrix,2,mean)
apply(power2_StrPBR_matrix,2,mean)
apply(power2_PS_matrix,2,mean)
apply(power2_HH_matrix,2,mean)

#diff=0.5*treatment effect
for(i in 1:repli){
  power2_cr<-power2(n,method="CR",0.5*diff)
  power2_StrPBR<-power2(n,method="StrPBR",0.5*diff)
  power2_PS<-power2(n,method="PS",0.5*diff)
  power2_HH<-power2(n,method="HH",0.5*diff)
  power2_cr_matrix<-rbind(power2_cr,power2_cr_matrix)
  power2_StrPBR_matrix<-rbind(power2_StrPBR,power2_StrPBR_matrix)
  power2_PS_matrix<-rbind(power2_PS,power2_PS_matrix)
  power2_HH_matrix<-rbind(power2_HH,power2_HH_matrix)
}

apply(power2_cr_matrix,2,mean)
apply(power2_StrPBR_matrix,2,mean)
apply(power2_PS_matrix,2,mean)
apply(power2_HH_matrix,2,mean)

#diff=treatment effect
for(i in 1:repli){
  power2_cr<-power2(n,method="CR",diff)
  power2_StrPBR<-power2(n,method="StrPBR",diff)
  power2_PS<-power2(n,method="PS",diff)
  power2_HH<-power2(n,method="HH",diff)
  power2_cr_matrix<-rbind(power2_cr,power2_cr_matrix)
  power2_StrPBR_matrix<-rbind(power2_StrPBR,power2_StrPBR_matrix)
  power2_PS_matrix<-rbind(power2_PS,power2_PS_matrix)
  power2_HH_matrix<-rbind(power2_HH,power2_HH_matrix)
}

apply(power2_cr_matrix,2,mean)
apply(power2_StrPBR_matrix,2,mean)
apply(power2_PS_matrix,2,mean)
apply(power2_HH_matrix,2,mean)

#diff=1.5*treatment effect
for(i in 1:repli){
  power2_cr<-power2(n,method="CR",1.5*diff)
  power2_StrPBR<-power2(n,method="StrPBR",1.5*diff)
  power2_PS<-power2(n,method="PS",1.5*diff)
  power2_HH<-power2(n,method="HH",1.5*diff)
  power2_cr_matrix<-rbind(power2_cr,power2_cr_matrix)
  power2_StrPBR_matrix<-rbind(power2_StrPBR,power2_StrPBR_matrix)
  power2_PS_matrix<-rbind(power2_PS,power2_PS_matrix)
  power2_HH_matrix<-rbind(power2_HH,power2_HH_matrix)
}

apply(power2_cr_matrix,2,mean)
apply(power2_StrPBR_matrix,2,mean)
apply(power2_PS_matrix,2,mean)
apply(power2_HH_matrix,2,mean)

##adjusted t-test
df4<-df3
df4$MarketSize<-as.numeric(df3$MarketSize)
df4$age.store.ind<-df3$age.store.ind+1

#adjusted t-test StrPBR
i_adjusted_matrix<-NULL
for (i in 1:repli) {
  t.str<-StrPBR(df4,bsize=4)$assignments
  t.str<-ifelse(t.str == "A", 2, 1)
  y<-b0+b2*z1+b3*z2+rnorm(n,0,sigma)
  
  df5<-t(cbind(df4,t.str,y))
  corr.test(df5)$p.value
  i_adjusted<-ifelse(corr.test(df5)$p.value>0.05,0,1)
  i_adjusted_matrix<-c(i_adjusted_matrix,i_adjusted)
}

adjust.type1.StrPBR<-mean(i_adjusted_matrix)

#adjusted t-test PS
i_adjusted_matrix<-NULL
for (i in 1:repli) {
  t.ps<-PocSimMIN(df4,weight=c(0.35,0.35,0.3),p=0.75)$assignments
  t.ps<-ifelse(t.ps == "A", 2, 1)
  y<-b0+b2*z1+b3*z2+rnorm(n,0,sigma)
  
  df5<-t(cbind(df4,t.ps,y))
  corr.test(df5)$p.value
  i_adjusted<-ifelse(corr.test(df5)$p.value>0.05,0,1)
  i_adjusted_matrix<-c(i_adjusted_matrix,i_adjusted)
}

adjust.type1.PS<-mean(i_adjusted_matrix)

#adjusted t-test HuHu
i_adjusted_matrix<-NULL
for (i in 1:repli) {
  t.hu<-HuHuCAR(df4,c(0.35,0.35,0.1,0.1,0.1),p=0.75)$assignments
  t.hu<-ifelse(t.hu == "A", 2, 1)
  y<-b0+b2*z1+b3*z2+rnorm(n,0,sigma)
  
  df5<-t(cbind(df4,t.hu,y))
  corr.test(df5)$p.value
  i_adjusted<-ifelse(corr.test(df5)$p.value>0.05,0,1)
  i_adjusted_matrix<-c(i_adjusted_matrix,i_adjusted)
}

adjust.type1.HuHu<-mean(i_adjusted_matrix)
#table 3
table3(n,1)

#table 4
table4(n,1000,"CR")
table4(n,1000,"StrPBR")
table4(n,1000,"PS")
table4(n,1000,"HuHu")

#table 5
table5(n,1000,"CR")
table5(n,1000,"StrPBR")
table5(n,1000,"PS")
table5(n,1000,"HuHu")

#table 6
table6(n,100,"CR")
table6(n,100,"StrPBR")
table6(n,100,"PS")
table6(n,100,"HuHu")
