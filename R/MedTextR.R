#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#




 rm(list=ls(all=TRUE))


library(foreign)
library(boot)
library(gdata)
library(ReadImages)
library(MASS)
library(RGtk2)

dtr = function (dd, nn1,nn2)
{ rr= dd*sqrt(nn1*nn2/(nn1+nn2)/(nn1+nn2-2+dd^2*(nn1*nn2)/(nn1+nn2)))
  return(rr) 
 }

 




effsz = function (eff, dic)
{ 
  esmall=.1
  emed =.3
  elarge=.5
  if (dic == 1) {esmall=.2}
  if (dic == 1) {emed =.5}
  if (dic == 1) {elarge =.8}
  btto =" less than small"
  if (abs(eff)>esmall) btto=" and a small effect size"
  if (abs(eff)>emed) btto=" and a medium effect size"
  if (abs(eff)>elarge) btto=" and a large effect size"
  return(btto) 
 }
effszm = function (dic)
{  ess=' (r = '
  if (dic == 1)ess=' (d = '
   return(ess)
}




cpow  = function (r, numb,df,alph,collin,rsqp)
 {z2T = qnorm(c(1-alph/2),mean=0,sd=1,lower.tail=TRUE)
 tr = r*sqrt((numb-4-df)*(1-collin)/(1-rsqp))
 es = (2*tr/sqrt(numb))**2
 powr= pnorm(c(sqrt((numb/4)*es)/(1+1.21*(z2T -1.06)/(numb- df - 2)) - z2T), mean=0,  sd=1,lower.tail=TRUE) + 1- pnorm(c(sqrt((numb/4)*es)/(1+1.21*(z2T-1.06)/(numb-2-df)) +  z2T), mean=0, sd=1,lower.tail=TRUE)
 return(powr)
 }

rround  = function (num,dig)
 {
if (dig==2) return(sprintf("%.2f", round(num,2)))
if (dig==3) return(sprintf("%.3f", round(num,3)))
if (dig==1) return(sprintf("%.1f", round(num,1)))
}

rrround  = function (num,dig)
 {
if (substring(as.character(rround(num,dig)),1,2)=="0.") return(substring(as.character(rround(num,dig)),2,5))
if (substring(as.character(rround(num,dig)),1,2)=="-0") return(gsub("-0","-", rround(num,dig)))
if (substring(as.character(rround(num,dig)),1,2)=="1.") return(substring(as.character(rround(num,dig)),1,5))
if (substring(as.character(rround(num,dig)),1,2)=="-1") return(substring(as.character(rround(num,dig)),1,6))
}



pval = function(num) 
 { if (num>=.001) xxx= paste(" (p = ",rrround(num,3),")",sep="")
   if (num<.001) xxx= " (p < .001)"
   return(xxx)
}

knt = function(pvv,alph) 
 { xxx=""
   if (pvv >= alph) xxx=" not"
   return(xxx)
}



Medtexty = function(button,user.data){

cat("", sep="\n")

cat("MedTextR has begun.", sep="\n")




ifilename <- filename$getText() 


trt = substr(ifilename, nchar(ifilename)-3+1, nchar(ifilename))
if (trt=='sav') OrDa = read.spss (ifilename,use.value.labels=FALSE,max.value.labels=Inf,to.data.frame=TRUE)
if (trt=='csv') OrDa = read.csv(file=ifilename,head=TRUE,sep=",")







# Input marcro agruments

MaDa <- c(MaDa=1)

MaDa <- as.data.frame(MaDa)

# Create dataset with new names


xvar1 <- xvar$getText() 
mvar1 <- mvar$getText() 
yvar1 <- yvar$getText() 
names(OrDa)[names(OrDa)==xvar1] <- "xvar"
names(OrDa)[names(OrDa)==yvar1] <- "yvar"
names(OrDa)[names(OrDa)==mvar1] <- "mvar"
MaDa$mn  <- mn$getText()
MaDa$xn <- xn$getText()
MaDa$yn <- yn$getText()
MaDa$drt1 <- drt$getText()
MaDa$ofilename <- ofile$getText()
MaDa$alpha <- alpha$getText()
MaDa$clist <- clis$getText()
MaDa$covn <- cnam$getText()
MaDa$trials <- ntrials$getText()
MaDa$rfilename <- rfile$getText()

MaDa$rfilename = paste(MaDa$drt1,MaDa$rfilename,sep="")  


MaDa$ofilename = paste(MaDa$drt1,MaDa$ofilename,sep="")  



MaDa$alpha <- as.numeric(MaDa$alpha)
MaDa$trials <- as.numeric(MaDa$trials)

MaDa$figgy=paste(MaDa$drt1,'MTRfigure.png',sep="")
MaDa$F1=paste(MaDa$drt1,'MTFig1.png',sep="")
MaDa$F2=paste(MaDa$drt1,'MTFig2.png',sep="")
MaDa$labx1 <- tlab1$getText()
MaDa$labx2 <- tlab2$getText()
MaDa$xm1 <- xm$active
MaDa$z2T=qnorm(c(1-MaDa$alpha/2),mean=0,sd=1,lower.tail=TRUE)



cat("MedTextR Output",file=MaDa$rfilename,sep="\n",append=FALSE)
cat("\n", file = MaDa$rfilename, append = TRUE)


 













MaDa$df1 = NA
MaDa$df2 = NA
MaDa$df3 = NA
MaDa$pcov1 = NA
MaDa$pcov2 = NA
MaDa$Fcov1 = NA
MaDa$Fcov2 = NA
MaDa$x2p = NA
MaDa$x2mp = NA


# Listwise deletion

MaDa$ntot = nrow(OrDa)
RaDa <- OrDa[!is.na(OrDa$xvar),]
RaDa <- RaDa[!is.na(RaDa$mvar),]
RaDa <- RaDa[!is.na(RaDa$yvar),]

#Covariate information
MaDa$ncov=length(strsplit(MaDa$clist,',')[[1]])
MaDa$df1 = MaDa$ncov
cvs= unlist(strsplit(MaDa$clist,","))
if (MaDa$ncov==0) cat=""
if (MaDa$ncov==1) cat="+c1"
if (MaDa$ncov==2) cat="+c1+c2"
if (MaDa$ncov==3) cat="+c1+c2+c3"
if (MaDa$ncov==4) cat="+c1+c2+c3+c4"
if (MaDa$ncov==5) cat="+c1+c2+c3+c4+c5"
if (MaDa$ncov==6) cat="+c1+c2+c3+c4+c5+c6"
if (MaDa$ncov==7) cat="+c1+c2+c3+c4+c5+c6+c7"
if (MaDa$ncov==8) cat="+c1+c2+c3+c4+c5+c6+c7+c8"
if (MaDa$ncov==9) cat="+c1+c2+c3+c4+c5+c6+c7+c8+c9"
if (MaDa$ncov==10) cat="+c1+c2+c3+c4+c5+c6+c7+c8+c9+c10"
#Change this to a loop someway
 RaDa$c1 = RaDa[[cvs[1]]]
 RaDa$c2 = RaDa[[cvs[2]]]
 RaDa$c3 = RaDa[[cvs[3]]]
 RaDa$c4 = RaDa[[cvs[4]]]
 RaDa$c5 = RaDa[[cvs[5]]]
 RaDa$c6 = RaDa[[cvs[6]]]
 RaDa$c7 = RaDa[[cvs[7]]]
 RaDa$c8 = RaDa[[cvs[8]]]
 RaDa$c9 = RaDa[[cvs[9]]]
 RaDa$c10 = RaDa[[cvs[10]]]
covs= list(RaDa$c1,RaDa$c2,RaDa$c3,RaDa$c4,RaDa$c5,RaDa$c6,RaDa$c7,RaDa$c8,RaDa$c9,RaDa$c10)
#Above should be a loop.
 x <- 1:10-MaDa$ncov
        for(i in seq(along=x)) { 
        covs <- covs[-(MaDa$ncov+1)]                
        }
#Remove missing covariate cases
#Below should be a loop.
if (MaDa$ncov>0) RaDa <- RaDa[!is.na(RaDa$c1),]
if (MaDa$ncov>1) RaDa <- RaDa[!is.na(RaDa$c2),]
if (MaDa$ncov>2) RaDa <- RaDa[!is.na(RaDa$c3),]
if (MaDa$ncov>3) RaDa <- RaDa[!is.na(RaDa$c4),]
if (MaDa$ncov>4) RaDa <- RaDa[!is.na(RaDa$c5),]
if (MaDa$ncov>5) RaDa <- RaDa[!is.na(RaDa$c6),]
if (MaDa$ncov>6) RaDa <- RaDa[!is.na(RaDa$c7),]
if (MaDa$ncov>7) RaDa <- RaDa[!is.na(RaDa$c8),]
if (MaDa$ncov>8) RaDa <- RaDa[!is.na(RaDa$c9),]
if (MaDa$ncov>9) RaDa <- RaDa[!is.na(RaDa$c10),]

# Add new variables to the dataset
 RaDa$m2=RaDa$mvar*RaDa$mvar
 RaDa$x2=RaDa$xvar*RaDa$xvar
 RaDa$xm=RaDa$mvar*RaDa$xvar

# Descriptives
MaDa$xmin=min(RaDa$xvar)
MaDa$xmax=max(RaDa$xvar)
MaDa$xmean=mean(RaDa$xvar)
MaDa$xsd=sd(RaDa$xvar)
MaDa$mmin=min(RaDa$mvar)
MaDa$mmax=max(RaDa$mvar)
MaDa$mmean=mean(RaDa$mvar)
MaDa$msd=sd(RaDa$mvar)
MaDa$ymin=min(RaDa$yvar)
MaDa$ymax=max(RaDa$yvar)
MaDa$ymean=mean(RaDa$yvar)
MaDa$ysd=sd(RaDa$yvar)
MaDa$nxn=length(RaDa$xvar)
MaDa$nmiss = MaDa$ntot - MaDa$nxn 
MaDa$df2= MaDa$nxn-MaDa$ncov-2
MaDa$df3=MaDa$df2-1

zz="Descriptives"
out<-capture.output(as.name(zz))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)
cat("\n", file = MaDa$rfilename, append = TRUE)
out<-capture.output(summary(RaDa, na.rm = TRUE))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)

 MaDa$xmeann=(MaDa$xmean-MaDa$xmin)/(MaDa$xmax-MaDa$xmin)
 zzz =(MaDa$xmax-MaDa$xmin)*sqrt(MaDa$nxn*MaDa$xmeann*(1-MaDa$xmeann)/(MaDa$nxn-1.0))
 xoox=abs(zzz-MaDa$xsd)/(zzz+MaDa$xsd)
 MaDa$di=0
 if (xoox <.0000001) {MaDa$di=1}




# Step 1

dog = paste("yvar ~ xvar",cat,sep="")

step1 <- lm(dog, data=RaDa)

ss1 <- summary(step1)
ss <- ss1[4]
sss <-  as.data.frame(ss)
MaDa$ccp=sss[2,4]
MaDa$ccc=sss[2,1]
MaDa$int1=sss[1,1]
MaDa$cse=sss[2,2]
step1.ci <- confint(step1, level= 1-MaDa$alpha)
MaDa$ccl=step1.ci[2,1]
MaDa$ccu=step1.ci[2,2]
MaDa$rsq1 <-  summary(step1)$r.squared
MaDa$evar1 <- (summary(step1)$sigma)^2
MaDa$ccb=MaDa$ccc*MaDa$xsd/MaDa$ysd




cat("\n", file = MaDa$rfilename, append = TRUE)
zz2="Step 1 Results"
out<-capture.output(as.name(zz2))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)
out<-capture.output(dog,summary(step1))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)
out<-capture.output(step1.ci)
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)









# Step 2

dog = paste("mvar ~ xvar",cat,sep="")


step2 <- lm(dog , data=RaDa)
ss2 <- summary(step2)
ss <- ss2[4]
sss <-  as.data.frame(ss)
MaDa$aap=sss[2,4]
MaDa$aaa=sss[2,1]
MaDa$ase=sss[2,2]
MaDa$int2=sss[1,1]
step2.ci <- confint(step2, level= 1-MaDa$alpha)
MaDa$aal=step2.ci[2,1]
MaDa$aau=step2.ci[2,2]
MaDa$rsq2 <-  summary(step2)$r.squared
MaDa$evar2 <- (summary(step2)$sigma)^2

MaDa$aab=MaDa$aaa*MaDa$xsd/MaDa$msd
x <- model.matrix(step2)
 lev <- hat(x)
 gs <- summary(step2)
 RaDa$stud <- step2$res/(gs$sig*sqrt(1-lev))
 RaDa$no1 <- ifelse(abs(RaDa$stud)>3.5, 1, 0)
 MaDa$nout1 <- sum(RaDa$no1)






cat("\n", file = MaDa$rfilename, append = TRUE)
zz2="Step 2 Results"
out<-capture.output(as.name(zz2))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)
out<-capture.output(as.name(dog),summary(step2))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)
out<-capture.output(step2.ci)
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)

 




# Steps 3 and 4

dog = paste("yvar ~ xvar + mvar",cat,sep="")


step34 <- lm(dog, data=RaDa)
ss34 <- summary(step34)
ss <- ss34[4]
sss <-  as.data.frame(ss)
MaDa$cpp=sss[2,4]
MaDa$cpc=sss[2,1]
MaDa$bbp=sss[3,4]
MaDa$bbb=sss[3,1]
MaDa$int3=sss[1,1]
MaDa$cpse=sss[2,2]
MaDa$bse=sss[3,2]
step34.ci <- confint(step34, level=1-MaDa$alpha)
MaDa$cpl=step34.ci[2,1]
MaDa$cpu=step34.ci[2,2]
MaDa$bbl=step34.ci[3,1]
MaDa$bbu=step34.ci[3,2]
MaDa$rsq3 <-  summary(step34)$r.squared
MaDa$evar3 <- (summary(step34)$sigma)^2 
MaDa$bbbb=MaDa$bbb*MaDa$msd/MaDa$ysd
MaDa$cpb=MaDa$cpc*MaDa$xsd/MaDa$ysd

 x <- model.matrix(step34)
 lev <- hat(x)
 gs <- summary(step34)
 RaDa$stud <- step34$res/(gs$sig*sqrt(1-lev))
 RaDa$no2 <- ifelse(abs(RaDa$stud)>3.5, 1, 0)
 MaDa$nout2 <- sum(RaDa$no2)



cat("\n", file = MaDa$rfilename, append = TRUE)
zz="Steps 3 and 4 Results"
out<-capture.output(as.name(zz))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)
out<-capture.output(as.name(dog),summary(step34))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)
out<-capture.output(step34.ci)
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)

if (MaDa$ncov>0)

{

# Step 2 (without covariates)

dog = "mvar ~ xvar"


step2w <- lm(dog , data=RaDa)
MaDa$rsq4 <-  summary(step2w)$r.squared




MaDa$Fcov1= (MaDa$rsq2-MaDa$rsq4)*(MaDa$df2)/((1-MaDa$rsq2)*MaDa$df1)
MaDa$Fcov1 <- as.numeric(MaDa$Fcov1)
MaDa$pcov1= 1- pf(MaDa$Fcov1, MaDa$df1,MaDa$df2)

cat("\n", file = MaDa$rfilename, append = TRUE)
zz="Step 2 without covariates"
out<-capture.output(as.name(zz))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)
out<-capture.output(as.name(dog),step2w)
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)



# Step 3 & 4 (without covariates)

dog = "yvar ~ xvar + mvar"


step34w <- lm(dog, data=RaDa)
ss34w <- summary(step34w)
MaDa$rsq5 <-  summary(step34w)$r.squared



MaDa$Fcov2= (MaDa$rsq3-MaDa$rsq5)*(MaDa$df3)/((1-MaDa$rsq3)*MaDa$df1)
MaDa$Fcov2 <- as.numeric(MaDa$Fcov2)
MaDa$pcov2= 1- pf(MaDa$Fcov2, MaDa$df1,MaDa$df3)

cat("\n", file = MaDa$rfilename, append = TRUE)
zz="Steps 3 and 4 without covariates"
out<-capture.output(as.name(zz))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)
out<-capture.output(as.name(dog),step34w)
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)




}


cat1=cat
if (cat=="") cat1="0"
dog = paste("yvar ~ ",cat1,sep="")
stepx <- lm(dog, data=RaDa)
MaDa$rsqY.C <-  summary(stepx)$r.squared

dog = paste("mvar ~ ",cat1,sep="")
stepx <- lm(dog, data=RaDa)
MaDa$rsqM.C <-  summary(stepx)$r.squared

dog = paste("xvar ~ ",cat1,sep="")
stepx <- lm(dog, data=RaDa)
MaDa$rsqX.C <-  summary(stepx)$r.squared

dog = paste("yvar ~ xvar",cat,sep="")
stepx <- lm(dog, data=RaDa)
MaDa$rsqY.CX <- summary(stepx)$r.squared

dog = paste("yvar ~ mvar",cat,sep="")
stepx <- lm(dog, data=RaDa)
MaDa$rsqY.CM <- summary(stepx)$r.squared

dog = paste("mvar ~ xvar",cat,sep="")
stepx <- lm(dog, data=RaDa)
MaDa$rsqM.CX <- summary(stepx)$r.squared

dog = paste("xvar ~ mvar",cat,sep="")
stepx <- lm(dog, data=RaDa)
MaDa$rsqX.CM <- summary(stepx)$r.squared

MaDa$aar=MaDa$aab*(sqrt((1-MaDa$rsqX.C)/(1-MaDa$rsqM.C)))
MaDa$ccr=MaDa$ccb*(sqrt((1-MaDa$rsqX.C)/(1-MaDa$rsqY.C)))
MaDa$cpr=MaDa$cpb*(sqrt((1-MaDa$rsqX.CM)/(1-MaDa$rsqY.CM)))
MaDa$bbr=MaDa$bbbb*(sqrt((1-MaDa$rsqM.CX)/(1-MaDa$rsqY.CX)))






# bootstap using Kelley's MBESS



ieest <- mediation(RaDa$xvar, RaDa$mvar, RaDa$yvar, ncv=MaDa$ncov,catyr=cat, conf.level = 1-MaDa$alpha, B = MaDa$trials,RaDa$c1,RaDa$c2,RaDa$c3)


cat("\n", file = MaDa$rfilename, append = TRUE)
zz="Bootstrap Results"
out<-capture.output(as.name(zz))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)
out<-capture.output(ieest$Bootstrap.Results)
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)


ss1 <- ieest[4]
ss2 = as.data.frame(ss1)
MaDa$blci <- ss2[1,2]
MaDa$buci <- ss2[1,3]
MaDa$bie <- ss2[1,1]



IE <- bsrep$V1

MaDa$bieb = mean(IE)

MaDa$bootsd = sqrt(var(IE))
MaDa$bootp=1-2*abs(.5- mean(IE+MaDa$bie-MaDa$bieb>0)) 




if (MaDa$di < 1) 

# X to M nonlinear

{


dog = paste("mvar ~ xvar + x2",cat,sep="")


xsquare <- lm(dog, data=RaDa)

ssx2 <- summary(xsquare)
ss <- ssx2[4]
sss <-  as.data.frame(ss)
MaDa$x2mp=sss[3,4]
MaDa$x2mb=sss[3,1]


cat("\n", file = MaDa$rfilename, append = TRUE)
zz="Nonlinear effect of X on M Results"
out<-capture.output(as.name(zz))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)
out<-capture.output(as.name(dog),xsquare)
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)


# X to Y nonlinear

dog = paste("yvar ~ xvar + mvar + x2",cat,sep="")


msquare <- lm(dog, data=RaDa)

ssm2 <- summary(msquare)
ss <- ssm2[4]
sss <-  as.data.frame(ss)
MaDa$x2p=sss[4,4]
MaDa$x2b=sss[4,1]


cat("\n", file = MaDa$rfilename, append = TRUE)
zz="Nonlinear effect of X on Y Results"
out<-capture.output(as.name(zz))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)
out<-capture.output(as.name(dog),msquare)
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)

}

# M to Y  nonlinear

dog = paste("yvar ~ xvar + mvar + m2",cat,sep="")


msquare <- lm(dog, data=RaDa)

ssm2 <- summary(msquare)
ss <- ssm2[4]
sss <-  as.data.frame(ss)
MaDa$m2p=sss[4,4]
MaDa$m2b=sss[4,1]

cat("\n", file = MaDa$rfilename, append = TRUE)
zz="Nonlinear effect of M on Y Results"
out<-capture.output(as.name(zz))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)
out<-capture.output(as.name(dog),msquare)
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)


# xm interaction

dog = paste("yvar ~ xvar + mvar+ xm",cat,sep="")


nonlin <- lm(dog, data=RaDa)

ssxm <- summary(nonlin)
ss <- ssxm[4]
sss <-  as.data.frame(ss)
MaDa$xmp=sss[4,4]
MaDa$xmb=sss[4,1]

cat("\n", file = MaDa$rfilename, append = TRUE)
zz="Interaction of X and M Results"
out<-capture.output(as.name(zz))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)
out<-capture.output(as.name(dog),nonlin)
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)



# Robust Regression

 
if (MaDa$ncov==0)  robust1 <- rlm(yvar ~ xvar , data=RaDa)
if (MaDa$ncov==1)  robust1 <- rlm(yvar ~ xvar + c1 , data=RaDa)
if (MaDa$ncov==2)  robust1 <- rlm(yvar ~ xvar + c1 + c2, data=RaDa)
if (MaDa$ncov==3)  robust1 <- rlm(yvar ~ xvar + c1 + c2 + c3, data=RaDa)
if (MaDa$ncov==3)  robust1 <- rlm(yvar ~ xvar + c1 + c2 + c3 + c4, data=RaDa)
if (MaDa$ncov==3)  robust1 <- rlm(yvar ~ xvar + c1 + c2 + c3 + c4 + c5, data=RaDa)
if (MaDa$ncov==3)  robust1 <- rlm(yvar ~ xvar + c1 + c2 + c3 + c4 + c5 + c6, data=RaDa)
if (MaDa$ncov==3)  robust1 <- rlm(yvar ~ xvar + c1 + c2 + c3 + c4 + c5 + c6 + c7, data=RaDa)
if (MaDa$ncov==3)  robust1 <- rlm(yvar ~ xvar + c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8, data=RaDa)
if (MaDa$ncov==3)  robust1 <- rlm(yvar ~ xvar + c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9, data=RaDa)
if (MaDa$ncov==3)  robust1 <- rlm(yvar ~ xvar + c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10, data=RaDa)

 rss1 <- summary(robust1)
 rss <- rss1[4]
 rsss <-  as.data.frame(rss)
 MaDa$rccc=rsss[2,1]
 MaDa$rint1=rsss[1,1]
 MaDa$rcse=rsss[2,2]
 MaDa$rccp= 2*(1-pnorm(c(abs(MaDa$rccc)/MaDa$rcse),mean=0,sd=1,lower.tail=TRUE))

 MaDa$rccl= MaDa$rccc - MaDa$z2T*MaDa$rcse
 MaDa$rccu=  MaDa$rccc + MaDa$z2T*MaDa$rcse
 MaDa$rccb=MaDa$ccc*MaDa$xsd/MaDa$ysd
 rss <- rss1[6]
 MaDa$resd1 <- as.data.frame(rss)

cat("\n", file = MaDa$rfilename, append = TRUE)
zz="Step 1 Robust Results"
out<-capture.output(as.name(zz))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)
out<-capture.output(rss1)
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)





 if (MaDa$ncov==0)  robust2 <- rlm(mvar ~ xvar, data=RaDa)
 if (MaDa$ncov==1)  robust2 <- rlm(mvar ~ xvar + c1, data=RaDa)
 if (MaDa$ncov==2)  robust2 <- rlm(mvar ~ xvar + c1 + c2, data=RaDa)
 if (MaDa$ncov==3)  robust2 <- rlm(mvar ~ xvar+ c1 + c2 + c3, data=RaDa)
 if (MaDa$ncov==3)  robust2 <- rlm(mvar ~ xvar+ c1 + c2 + c3 + c4, data=RaDa)
 if (MaDa$ncov==3)  robust2 <- rlm(mvar ~ xvar+ c1 + c2 + c3 + c4 + c5, data=RaDa)
 if (MaDa$ncov==3)  robust2 <- rlm(mvar ~ xvar+ c1 + c2 + c3 + c4 + c5 + c6, data=RaDa)
 if (MaDa$ncov==3)  robust2 <- rlm(mvar ~ xvar+ c1 + c2 + c3 + c4 + c5 + c6 + c7, data=RaDa)
 if (MaDa$ncov==3)  robust2 <- rlm(mvar ~ xvar+ c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8, data=RaDa)
 if (MaDa$ncov==3)  robust2 <- rlm(mvar ~ xvar+ c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9, data=RaDa)
 if (MaDa$ncov==3)  robust2 <- rlm(mvar ~ xvar+ c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8 + c9 + c10, data=RaDa)

 rss2 <- summary(robust2)
 rss <- rss2[4]
 rsss <-  as.data.frame(rss)
 MaDa$raaa=rsss[2,1]
 MaDa$rint2=rsss[1,1]
 MaDa$rase=rsss[2,2]
 MaDa$raap= 2*(1-pnorm(c(abs(MaDa$rccc)/MaDa$rcse),mean=0,sd=1,lower.tail=TRUE))
 MaDa$raal= MaDa$raaa - MaDa$z2T*MaDa$rase
 MaDa$raau=  MaDa$raaa + MaDa$z2T*MaDa$rase
 MaDa$raab=MaDa$aaa*MaDa$xsd/MaDa$msd
 rss <- rss2[6]
 MaDa$resd2 <- as.data.frame(rss)

cat("\n", file = MaDa$rfilename, append = TRUE)
zz="Step 2 Robust Results"
out<-capture.output(as.name(zz))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)
out<-capture.output(rss2)
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)

 


 if (MaDa$ncov==0)  robust34 <- rlm(yvar ~ xvar + mvar, data=RaDa)
 if (MaDa$ncov==1)  robust34 <- rlm(yvar ~ xvar + mvar + c1, data=RaDa)
 if (MaDa$ncov==2)  robust34 <- rlm(yvar ~ xvar + mvar+ c1 + c2, data=RaDa)
 if (MaDa$ncov==3)  robust34 <- rlm(yvar ~ xvar + mvar + c1 + c2 + c3, data=RaDa)
 if (MaDa$ncov==4)  robust34 <- rlm(yvar ~ xvar + mvar + c1 + c2 + c3+ c4, data=RaDa)
 if (MaDa$ncov==5)  robust34 <- rlm(yvar ~ xvar + mvar + c1 + c2 + c3+ c4+ c5, data=RaDa)
 if (MaDa$ncov==6)  robust34 <- rlm(yvar ~ xvar + mvar + c1 + c2 + c3+ c4+ c5+ c6, data=RaDa)
 if (MaDa$ncov==7)  robust34 <- rlm(yvar ~ xvar + mvar + c1 + c2 + c3+ c4+ c5+ c6+ c7, data=RaDa)
 if (MaDa$ncov==8)  robust34 <- rlm(yvar ~ xvar + mvar + c1 + c2 + c3+ c4+ c5+ c6+ c7+ c8, data=RaDa)
 if (MaDa$ncov==9)  robust34 <- rlm(yvar ~ xvar + mvar + c1 + c2 + c3+ c4+ c5+ c6+ c7+ c8+ c9, data=RaDa)
 if (MaDa$ncov==10)  robust34 <- rlm(yvar ~ xvar + mvar + c1 + c2 + c3+ c4+ c5+ c6+ c7+ c8+ c9+ c10, data=RaDa)


 rss34 <- summary(robust34)
 rss <- rss34[4]
 rsss <-  as.data.frame(rss)
 MaDa$rbbb=rsss[3,1]
 MaDa$rcpc=rsss[2,1]
 MaDa$rint3=rsss[1,1]
 MaDa$rbse=rsss[3,2]
 MaDa$rcpse=rsss[2,2]
 MaDa$rbbp= 2*(1 - pnorm(c(abs(MaDa$rbbb)/MaDa$rbse),mean=0,sd=1,lower.tail=TRUE))
 MaDa$rcpp= 2*(1 - pnorm(c(abs(MaDa$rcpc)/MaDa$rcpse),mean=0,sd=1,lower.tail=TRUE))
 MaDa$rbbl = MaDa$rbbb - MaDa$z2T*MaDa$rbse
 MaDa$rbbu = MaDa$rbbb + MaDa$z2T*MaDa$rbse
 MaDa$rcpl = MaDa$rcpc - MaDa$z2T*MaDa$rcpse
 MaDa$rcpu = MaDa$rcpc + MaDa$z2T*MaDa$rcpse
 MaDa$rbbbb=MaDa$rbbb*MaDa$msd/MaDa$ysd
 MaDa$rcpb=MaDa$rcpc*MaDa$xsd/MaDa$ysd
 rss <- rss34[6]
 MaDa$resd34 <- as.data.frame(rss)

cat("\n", file = MaDa$rfilename, append = TRUE)
zz="Steps 3 and 4 Robust Results"
out<-capture.output(as.name(zz))
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)
out<-capture.output(rss34)
cat(out,file=MaDa$rfilename,sep="\n",append=TRUE)


 
# Save MaDa as a csv file

# write.table(MaDa,file="MaDa.csv",sep=",",row.names=F)




# Preliminary stuff
 if (is.na(MaDa$pcov1) == "TRUE") {MaDa$pcov1 = 1}
 if (is.na(MaDa$pcov2) == "TRUE") {MaDa$pcov2 = 1}
 if (is.na(MaDa$x2p) == "TRUE") {MaDa$x2p = 1}
 if (is.na(MaDa$x2mp) == "TRUE") {MaDa$x2mp = 1}
 if (is.na(MaDa$m2p) == "TRUE") {MaDa$m2p = 1}


 MaDa$wnum = 0
 MaDa$wrn = ""

# Start Program
 
 if (MaDa$di==0) MaDa$txt9 = "        Note that effect sizes are partial correlations (r)."
 if (MaDa$di==1) MaDa$txt9 = paste("        Note that effect sizes are partial correlations (r) unless the predictor is ",MaDa$xn," where it is Cohen's d.",sep="")
 if (MaDa$di==0) MaDa$txt9 = paste(MaDa$txt9,"  Because an indirect effect is the product of two effect sizes, the effect size is the product of partial correlations (r*r).",sep="")
 if (MaDa$di==1) MaDa$txt9 = paste(MaDa$txt9,"  Because an indirect effect is the product of two effect sizes, the effect size is Cohen's d times the partial correlation (d*r).",sep="")
 tto="equals"
 if (MaDa$ncov==1) tto="and the covariate equal"
 if (MaDa$ncov==2) tto="and the covariates equal"
 if (MaDa$di==1) MaDa$txt9 = paste(MaDa$txt9,"  All predicted means presume that the mediator ",tto," zero.",sep="")
 if (MaDa$labx1 == "" & MaDa$di==1) MaDa$labx1="One"
 if (MaDa$labx2 == "" & MaDa$di==1) MaDa$labx2="Two"
 if (MaDa$labx2 == "Two" & MaDa$di==1) {MaDa$wnum = MaDa$wnum+1}
 if (MaDa$labx2 == "Two" & MaDa$di==1) {MaDa$wrn=paste(MaDa$wrn,'  ', tmp,'.  No "values" were given to ',MaDa$xn,', and so "One" and "Two" have been assigned.')}

 

 MaDa$ymeann=(MaDa$ymean-MaDa$ymin)/(MaDa$ymax-MaDa$ymin)
 zzz =(MaDa$ymax-MaDa$ymin)*sqrt(MaDa$nxn*MaDa$ymeann*(1-MaDa$ymeann)/(MaDa$nxn-1))
 xoox=abs(zzz-MaDa$ysd)/(zzz+MaDa$ysd)
 if (xoox < .0000001) MaDa$wnum=MaDa$wnum+1
 if (xoox < .0000001) {MaDa$wrn=paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),".  The variable ",MaDa$yn," is a dichotomy and logistic regression and not ordinary regression should be used.  ",sep="")}
 
 MaDa$mmeann=(MaDa$mmean-MaDa$mmin)/(MaDa$mmax-MaDa$mmin)
 zzz =(MaDa$mmax-MaDa$mmin)*sqrt(MaDa$nxn*MaDa$mmeann*(1-MaDa$mmeann)/(MaDa$nxn-1))
 xoox=abs(zzz-MaDa$msd)/(zzz+MaDa$msd)
 if (xoox < .0000001) {MaDa$dimm=1}
 if (xoox < .0000001) {MaDa$wnum=MaDa$wnum+1}
 if (xoox < .0000001) {MaDa$wrn=paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),".  The variable ",MaDa$mn," is a dichotomy and logistic regression and not ordinary regression should be used.",sep="")}
 
 if (.05 <= MaDa$alpha & MaDa$nxn > 400) {MaDa$wnum = MaDa$wnum + 1}
 if  (.05 <= MaDa$alpha & MaDa$nxn > 400) {MaDa$wrn = paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),".  Given the large sample size, alpha might be lowered.  ",sep="")} 
 if (MaDa$nxn < 41) {MaDa$wnum=MaDa$wnum+1}
 if (MaDa$nxn < 41) {MaDa$wrn=paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),".  The small sample size might preclude a mediational analysis.  ",sep="")} 
 xkx=0
 if (MaDa$mmax*MaDa$mmin >0.00000001) {xkx=1}
 if (MaDa$di == 1 & xkx == 1) {MaDa$wnum=MaDa$wnum+1}
 if (MaDa$di == 1 & xkx == 1) {MaDa$wrn=paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),".  As zero does not appear to be meaningful value for ",MaDa$mn,", grand mean centering the mediator might be considered.",sep="")}
 zzz =(MaDa$mmax-MaDa$mmin)*sqrt(MaDa$nxn*MaDa$mmeann*(1-MaDa$mmeann)/(MaDa$nxn-1))
 xoox=abs(zzz-MaDa$msd)/(zzz+MaDa$msd)
 if (xoox < .0000001) {MaDa$wnum=MaDa$wnum+1}
 MaDa$dimm=0
 if (xoox < .0000001) MaDa$wrn=paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),".  The variable ",MaDa$mn," is a dichotomy and logistic regression and not ordinary regression should be used.",sep="")
 if (xoox < .0000001) MaDa$dimm=1
 MaDa$xkx=0
 if (MaDa$mmax*MaDa$mmin >0.00000001) {MaDa$xkx=1}
 if (MaDa$di == 1 & MaDa$xkx == 1) {MaDa$wnum=MaDa$wnum+1}
 if (MaDa$di == 1 & MaDa$xkx == 1) {MaDa$wrn=paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),".  As zero does not appear to be meaningful value for ",MaDa$mn,", grand mean centering the mediator might be considered.",sep="")}

 if (MaDa$nout1>0) {MaDa$wnum=MaDa$wnum+1}
 if (MaDa$nout1>1) {MaDa$wrn = paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),".  There are ",round(MaDa$nout1,digits=0)," outliers (studentized residual greater than 3.5) for the variable ",MaDa$mn,".  Examine the output to see what observations are considered to be outliers.",sep="")}
 if (MaDa$nout1 == 1) {MaDa$wrn = paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),".  There is one outlier (studentized residual greater than 3.5) for the variable ",MaDa$mn,".  Examine the output to see what observations are considered to be outliers.",sep="")}

 if (MaDa$nout2>0) MaDa$wnum = MaDa$wnum+1
 if (MaDa$nout2>1) {MaDa$wrn = paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),".  There are ",round(MaDa$nout2,digits=0)," outliers (studentized residual greater than 3.5) for the variable ",MaDa$yn,".  Examine the output to see what observations are considered to be outliers.",sep="")}
 if (MaDa$nout2 == 1) {MaDa$wrn = paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),".  There is one outlier (studentized residual greater than 3.5) for the variable ",MaDa$yn,".  Examine the output to see what observations are considered to be outliers.",sep="")}

 MaDa$aarr = MaDa$aar
 MaDa$cprr = MaDa$cpr
 MaDa$ccrr = MaDa$ccr


 if (MaDa$aarr**2>.359999999) {MaDa$wnum=MaDa$wnum+1}
 if (MaDa$aarr**2 > .35999999) {MaDa$wrn=paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),".  There is severe multicollinearity of ",rrround(MaDa$aarr,3)," between ",MaDa$xn," and ",MaDa$mn," which may compromise the Step 3 and 4 tests.",sep="")}


 if (MaDa$rsqX.C>.35999999) {MaDa$wnum=MaDa$wnum+1}
 if (MaDa$rsqX.C>.35999999 & MaDa$ncov>1) {MaDa$wrn=paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),".  There is severe multicollinearity of ",rrround(MaDa$rsqX.C,3)," between ",MaDa$xn," and the covariates which may compromise the Step 1 and 4 tests.",sep="")}
 if (MaDa$rsqX.C>.35999999 & MaDa$ncov==1) {MaDa$wrn=paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),".  There is severe multicollinearity of ",rrround(MaDa$rsqX.C,3)," between ",MaDa$xn," and the covariate which may compromise the Step 1 and 4 tests.",sep="")}

 if (MaDa$rsqM.C>.35999999) {MaDa$wnum=MaDa$wnum+1}
 if (MaDa$rsqM.C>.35999999 & MaDa$ncov>1) {MaDa$wrn=paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),".  There is severe multicollinearity of ",rrround(MaDa$rsqM.C,3)," between ",MaDa$mn," and the covariates which may compromise the Step 3 test.",sep="")}
 if (MaDa$rsqM.C>.35999999 & MaDa$ncov==1) {MaDa$wrn=paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),".  There is severe multicollinearity of ",rrround(MaDa$rsqM.C,3)," between ",MaDa$mn," and the covariate which may compromise the Step 3 test.",sep="")}




 if (MaDa$di == 1) {MaDa$ccr=MaDa$ccc/(sqrt(MaDa$evar1)*(MaDa$xmax-MaDa$xmin))}
 if (MaDa$di == 1) {MaDa$aar=MaDa$aaa/(sqrt(MaDa$evar2)*(MaDa$xmax-MaDa$xmin))}
 if (MaDa$di == 1) {MaDa$cpr=MaDa$cpc/(sqrt(MaDa$evar3)*(MaDa$xmax-MaDa$xmin))}
 


 MaDa$txt1 = paste("        The causal variable or X is ",MaDa$xn,",",sep="") 
 if (MaDa$xm1 == "TRUE")   {MaDa$txt1 = paste(MaDa$txt1," a manipulated variable,", sep="")}
 if (MaDa$di == 1) {MaDa$txt1 = paste(MaDa$txt1," and is a dichotomy", sep="")}
 lg=100*(MaDa$xmean -MaDa$xmin)/(MaDa$xmax-MaDa$xmin)
 ug=100 - lg
 if (MaDa$di == 1) {MaDa$txt1=paste(MaDa$txt1," with ",round(lg,digit=1),"% ",MaDa$labx1," and ",round(ug,digit=1),"% ",MaDa$labx2,",",sep="")}  

 
 MaDa$txt1=paste(MaDa$txt1," the outcome variable or Y variable is ",MaDa$yn,", and ",sep="") 
 MaDa$txt1 = paste(MaDa$txt1,"the mediator or M is ",MaDa$mn,".",sep="") 
 MaDa$txt1=paste(MaDa$txt1,"  The causal mediational model is as follows:  The variable ",MaDa$xn," is presumed to cause",sep="")
 MaDa$txt1=paste(MaDa$txt1," ",MaDa$mn,", which in turn is presumed to cause ",MaDa$yn,".",sep="")
 MaDa$txt1=paste(MaDa$txt1,"  If there were complete mediation, then the causal effect of ",MaDa$xn," on ",MaDa$yn," controlling for ",MaDa$mn," would be zero.",sep="")
ttu="s"
if (MaDa$ncov==1) ttu="" 
ttx=""
if (length(MaDa$covn)>1) ttx =paste(" (",MaDa$covn,")",sep="")
if (MaDa$ncov>0) MaDa$txt1=paste(MaDa$txt1,"  For all analyses, the effects of the covariate",ttu,ttx," are controlled.",sep="")



# Descriptives


 
 MaDa$txt8=paste("        There are a total of ",round(MaDa$ntot,digit=0)," observations",sep="")
 if (MaDa$nmiss==0) MaDa$txt8=paste(MaDa$txt8," with no missing data.",sep="")
 if (MaDa$nmiss>0) MaDa$txt8=paste(MaDa$txt8,".",sep="")

 if (MaDa$nmiss > 1) {MaDa$txt8=paste(MaDa$txt8,"  However, ",round(MaDa$nmiss,digit=0)," observations are missing for one or more of the variables, which makes the sample size for analysis equal to ",round(MaDa$nxn,digit=0),".",sep="")}
 if (MaDa$nmiss == 1) {MaDa$txt8=paste(MaDa$txt8,"  However, one case is missing for one or more of the variables, which makes the sample size for analysis equal to ",round(MaDa$nxn,digit=0),".",sep="")}
 if (MaDa$nmiss > 0) {MaDa$txt8=paste(MaDa$txt8,"  Thus, listwise deletion was used.",sep="")}
 if (MaDa$nmiss > 9) MaDa$wnum=MaDa$wnum+1
 if (MaDa$nmiss > 9) {MaDa$wrn = paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),".  There are sufficient missing data so as to make the use of listwise missing data option less than optimal.  Better strategies for missing data (e.g., multiple imputation or" ,sep="")}
 if (MaDa$nmiss > 9) {MaDa$wrn = paste(MaDa$wrn,"  ",round(MaDa$wnum,digit=0),". full information maximum likelihood) should be considered.",sep="")}
 MaDa$txt8=paste(MaDa$txt8,"  The means and standard deviations are presented in Table 1.",sep="")
 ttu=as.character(round(MaDa$df1,digit=0))
 if (MaDa$df1 == 2) ttu="two" 
 if (MaDa$df1 == 3) ttu="three" 
 ttx=""
 if (length(MaDa$covn)>1) ttx =paste(" (",MaDa$covn,")",sep="")

 if (MaDa$df1 > 1) {MaDa$txt8 = paste(MaDa$txt8,"  There are ",ttu," covariates",ttx,", that are controlled in all analyses.",sep="")}
 if (MaDa$df1 == 1) {MaDa$txt8 = paste(MaDa$txt8,"  There is one covariate",ttx," that is controlled in all analyses.",sep="")}
 if (MaDa$df1>0) ttu=pval( MaDa$pcov1)
 if (MaDa$pcov1>= MaDa$alpha & MaDa$df1==1) {tte="covariate does not explain"}
 if (MaDa$pcov1>= MaDa$alpha & MaDa$df1>1) {tte="covariates do not explain"}
 if (MaDa$pcov1< MaDa$alpha & MaDa$df1==1)  {tte="covariate explains"}
 if (MaDa$pcov1< MaDa$alpha & MaDa$df1>1)  {tte="covariates explain"}
 tti=paste(" F(", round(MaDa$df1, digit = 0),", ", round(MaDa$df2, digit = 0),") = ", rround(MaDa$Fcov1,  3), ttu,".",sep="")
 if (MaDa$df1>0) MaDa$txt8 = paste(MaDa$txt8,"  The ", tte," a statistically significant amount of variance of the mediator, ", MaDa$mn,",", tti,sep="")
 tto= rround(sqrt(MaDa$evar2),  3)
 tti= rround(MaDa$evar2,  3)
 tte=""
 if (MaDa$df1>1) tte=" and the covariates"
 if (MaDa$df1==1) tte=" and the covariate"
  
 MaDa$txt8 = paste(MaDa$txt8,"  The unexplained variance in ",MaDa$mn," is equal to ", tti," (sd = ", tto,") controlling for ",MaDa$xn, tte,", with a multiple correlation for the regression equation of ", rrround(sqrt(MaDa$rsq2),  3),".",sep="")
 if (MaDa$df1>0) ttu=pval(MaDa$pcov2) 
 if (MaDa$pcov2 >= MaDa$alpha & MaDa$df1==1)  {tte="covariate does not explain"}
 if (MaDa$pcov2 >= MaDa$alpha & MaDa$df1>1)  {tte="covariates do not explain"}
 if (MaDa$pcov2 < MaDa$alpha & MaDa$df1==1)  {tte="covariate explains"}
 if (MaDa$pcov2 < MaDa$alpha & MaDa$df1>1) {tte="covariates explain"}
 tto=sqrt(MaDa$evar3)
 tti=paste(", F(",round(MaDa$df1,digit=0),", ",round(MaDa$df3,digit=0),") = ",rround(MaDa$Fcov2,3),ttu,".",sep="")
 if (MaDa$df1>0) MaDa$txt8=paste(MaDa$txt8,"  The ",tte," a statistically significant amount of variance of the outcome, ", MaDa$yn,tti,sep="")
 tte=""
 if (MaDa$df1>1) tte=", and the covariates"
 if (MaDa$df1 == 1) tte=", and the covariate"
 tta=" and"
 if (MaDa$df1>0) tta=","
 MaDa$txt8 = paste(MaDa$txt8,"  The unexplained variance in ",MaDa$yn," is equal to ",rround(MaDa$evar3,3)," (sd = ",rround(tto,3),") controlling for ",MaDa$xn,tta," ",MaDa$mn,tte,sep="") 
 MaDa$txt8 = paste(MaDa$txt8,", with a multiple correlation for the regression equation of ",rrround(sqrt(MaDa$rsq3),3),".",sep="")
 tte= round(MaDa$alpha,digit=4)
 if (MaDa$alpha==.005) tte=".005"
 if (MaDa$alpha==.05) tte=".05"
 if (MaDa$alpha==.01) tte=".01"
 if (MaDa$alpha==.02) tte=".02"
 if (MaDa$alpha==.10) tte=".10"
 
 MaDa$txt8 = paste(MaDa$txt8,"  For all the analyses, alpha is set at ",tte,".",sep="") 

# Baron and Kenny steps
 tta1= round(100-100*MaDa$alpha, digit = 0)
if (MaDa$alpha < .01) tta1= round(100-100*MaDa$alpha, digit = 1)

 tto= rround(MaDa$ccc,  3)
 ttu=pval(MaDa$ccp) 
 tmp="        The results of the four Baron and Kenny (1986) steps, which are summarized in Table 2, are as follows."
 tmp = paste(tmp," The effect of ",MaDa$xn," on ",MaDa$yn," or path c is equal to ",tto, ttu,sep="")
 tta= rround(MaDa$ccl,  3)
 tte= rround(MaDa$ccu,  3)
 tmp = paste(tmp,", with a ",tta1,"% confidence interval of ",tta," to ", tte,sep="")
 ttu=paste(tmp, effsz(MaDa$ccr,MaDa$di),effszm(MaDa$di), rrround(MaDa$ccr,  3),").",sep="")
 lowmean=MaDa$int1+ MaDa$ccc*MaDa$xmin
 himean=MaDa$int1+ MaDa$ccc*MaDa$xmax
 tta= rround(lowmean,  3)
 tmp= rround(himean,  3)
 if (MaDa$di == 1) ttu=paste(ttu," The least squares mean for ",MaDa$labx1," is equal to ",tta," and the mean for ",MaDa$labx2," is equal to ", tmp,".",sep="")
 MaDa$txt7=paste( ttu,"  Step 1 has",knt(MaDa$ccp,MaDa$alpha)," been passed.",sep="")
 tto= rround(MaDa$aaa,  3)
 ttu=pval(MaDa$aap)
 MaDa$txt4 = paste(" The effect of ",MaDa$xn," on ",MaDa$mn," or path a is equal to ", tto , ttu,sep="")
 tta= rround(MaDa$aal,  3)
 tte= rround(MaDa$aau,  3)
 MaDa$txt4= paste(MaDa$txt4,", with a ",tta1,"% confidence interval of ",tta," to ", tte,sep="")
 ttu=paste(MaDa$txt4,effsz(MaDa$aar,MaDa$di),effszm(MaDa$di), rrround(MaDa$aar,  3),").",sep="")
 lowmean=MaDa$int2+MaDa$aaa*MaDa$xmin
 himean=MaDa$int2+MaDa$aaa*MaDa$xmax
 tta= rround(lowmean,  3)
 tmp= rround(himean,  3)
 if (MaDa$di == 1) ttu=paste( ttu," The least squares mean for ",MaDa$labx1," is equal to ",tta," and the mean for ",MaDa$labx2," is equal to ", tmp,".",sep="")
 MaDa$txt3=paste( ttu,"  Step 2 has", knt(MaDa$aap,MaDa$alpha)," been passed.",sep="")

 tto= rround(MaDa$bbb,  3)
 ttu=pval(MaDa$bbp)
 MaDa$txt5 = paste(" The effect of ",MaDa$mn," on ",MaDa$yn," controlling for ",MaDa$xn," or path b is equal to ", tto, ttu,sep="")
 tta= rround(MaDa$bbl,  3)
 tte= rround(MaDa$bbu,  3)
 MaDa$txt5 = paste(MaDa$txt5,", with a ",tta1,"% confidence interval of ",tta," to ", tte,sep="")
 ttu=paste(MaDa$txt5,effsz(MaDa$bbr,0),effszm(0), rrround(MaDa$bbr,  3),").",sep="")
 MaDa$txt4=paste( ttu,"  Step 3 has", knt(MaDa$bbp,MaDa$alpha)," been passed.",sep="")
 tto= rround(MaDa$cpc,  3)
 ttu=pval(MaDa$cpp)
 MaDa$txt6=paste(" The effect of ",MaDa$xn," on ",MaDa$yn," controlling for ",MaDa$mn," or path c' is equal to ", tto, ttu,sep="")
 tta= rround(MaDa$cpl,  3)
 tte= rround(MaDa$cpu,  3)
 MaDa$txt6 = paste(MaDa$txt6,", with a ",tta1,"% confidence interval of ",tta," to ", tte,sep="")
 ttu=paste(MaDa$txt6,effsz(MaDa$cpr,MaDa$di),effszm(MaDa$di), rrround(MaDa$cpr,  3),").",sep="")
 lowmean=MaDa$int1+MaDa$cpc*MaDa$xmin
 himean=MaDa$int1+MaDa$cpc*MaDa$xmax
 tta= rround(lowmean,  3)
 tmp= rround(himean,  3)
 if (MaDa$di == 1) ttu=paste( ttu,"  The least squares mean for ",MaDa$labx1," is equal to ",tta," and the mean for ",MaDa$labx2," is equal to ", tmp,".",sep="")
 MaDa$txt5=paste( ttu,"  Step 4 has", knt(MaDa$alpha,MaDa$cpp)," been passed.",sep="")
 MaDa$txt2 = paste(MaDa$txt7," ",MaDa$txt3," ",MaDa$txt4," ",MaDa$txt5,sep="")





 MaDa$txt2 = paste(MaDa$txt2,"  A mediational diagram for unstandardized estimates is contained in Figure 1 (see ",MaDa$F1,") and for standardized estimates is contained in Figure 2 (see ",MaDa$F2,").",sep="")
 MaDa$txt2 = paste(MaDa$txt2,"  (In contemporary analyses, Baron and Kenny (1986) steps are no longer reported, but rather total, direct, and indirect effects are reported and tested.)",sep="")

# Robust Results


 tto= rround(MaDa$rccc,  3)
 ttu=pval(MaDa$rccp)
 MaDa$txt11="        The results of the four Baron and Kenny (1986) steps using robust regression are as follows:"
 MaDa$txt11 = paste(MaDa$txt11,"  Huber weighting (Huber, 1964) is used and observations with small residuals are given more weight than observations with larger residuals.",sep="")
 MaDa$txt11 = paste(MaDa$txt11," The effect of ",MaDa$xn," on ",MaDa$yn," or path c is equal to ",tto, ttu,sep="")
 tta= rround(MaDa$rccl,  3)
 tte= rround(MaDa$rccu,  3)
 MaDa$txt11 = paste(MaDa$txt11,", with a ",tta1,"% confidence interval of ",tta," to ", tte,".",sep="")
 MaDa$txt11=paste(MaDa$txt11,"  Step 1 has", knt(MaDa$rccp,MaDa$alpha)," been passed.",sep="")
 f1=0
 if (MaDa$rccp >= MaDa$alpha & MaDa$ccp < MaDa$alpha) f1 = 1
 if (MaDa$rccp < MaDa$alpha & MaDa$ccp >= MaDa$alpha) f1 = 1
 if (f1 == 1) MaDa$txt11=paste(MaDa$txt11,"  For Step 1, the robust regression results lead to a different conclusion.",sep="")
 tto= rround(MaDa$raaa,  3)
 ttu=pval(MaDa$raap)
 MaDa$txt11 = paste(MaDa$txt11," The effect of ",MaDa$xn," on ",MaDa$mn," or path a is equal to ", tto , ttu,sep="")
 tta= rround(MaDa$raal,  3)
 tte= rround(MaDa$raau,  3)
 MaDa$txt11= paste(MaDa$txt11,", with a ",tta1,"% confidence interval of ",tta," to ", tte,".",sep="")
 MaDa$txt11=paste(MaDa$txt11,"  Step 2 has", knt(MaDa$raap,MaDa$alpha)," been passed.",sep="")

 f2=0
 if (MaDa$raap >= MaDa$alpha & MaDa$aap < MaDa$alpha) f2 = 1
 if (MaDa$raap < MaDa$alpha & MaDa$aap >= MaDa$alpha) f2 = 1
 if (f2 == 1) MaDa$txt11=paste(MaDa$txt11,"  For Step 2, the robust regression results lead to a different conclusion.",sep="")
 tto= rround(MaDa$rbbb,  3)
 ttu=pval(MaDa$rbbp)
 MaDa$txt11 = paste(MaDa$txt11," The effect of ",MaDa$mn," on ",MaDa$yn," controlling for ",MaDa$xn," or path b is equal to ", tto, ttu,sep="")
 tta= rround(MaDa$rbbl,  3)
 tte= rround(MaDa$rbbu,  3)
 MaDa$txt11 = paste(MaDa$txt11,", with a ",tta1,"% confidence interval of ",tta," to ", tte,".",sep="")
 if (tto == " ")ttu=paste(MaDa$txt5,".",sep="")
 MaDa$txt11=paste(MaDa$txt11,"  Step 3 has", knt(MaDa$rbbp,MaDa$alpha)," been passed.",sep="")


 f3=0
 if (MaDa$rbbp >= MaDa$alpha & MaDa$bbp < MaDa$alpha) f3 = 1
 if (MaDa$rbbp < MaDa$alpha & MaDa$bbp >= MaDa$alpha) f3 = 1
 if (f3 == 1) MaDa$txt11=paste(MaDa$txt11,"  For Step 3, the robust regression results lead to a different conclusion.",sep="")
 tto= rround(MaDa$rcpc,  3)
 ttu=pval(MaDa$rcpp)
 MaDa$txt11=paste(MaDa$txt11," The effect of ",MaDa$xn," on ",MaDa$yn," controlling for ",MaDa$mn," or path c' is equal to ", tto, ttu,sep="")
 tta= rround(MaDa$rcpl,  3)
 tte= rround(MaDa$rcpu,  3)
 MaDa$txt11 = paste(MaDa$txt11,", with a ",tta1,"% confidence interval of ",tta," to ", tte,".",sep="")
 MaDa$txt11=paste(MaDa$txt11,"  Step 4 has", knt(MaDa$alpha,MaDa$rcpp)," been passed.",sep="")


 f4=0
 if (MaDa$rcpp >= MaDa$alpha & MaDa$cpp < MaDa$alpha) f4 = 1
 if (MaDa$rcpp < MaDa$alpha & MaDa$cpp >= MaDa$alpha) f4 = 1
 if (f4 == 1) MaDa$txt11=paste(MaDa$txt11,"  For Step 4, the robust regression results lead to a different conclusion.",sep="")
 if (f1+f2+f3+f4>0) MaDa$wnum=MaDa$wnum+1
 if (f1+f2+f3+f4>0) MaDa$wrn = paste(MaDa$wrn,"  ", round(MaDa$wnum, digit = 0),".  The results using robust regression lead to different conclusions than that using ordinary least squares.",sep="")
 if (f1+f2+f3+f4==0) MaDa$txt11=paste(MaDa$txt11,"  Robust methods yield essentially the same conclusions as ordinary least squares.",sep="")
# Indirect Effects
 MaDa$ide=MaDa$aaa*MaDa$bbb
 MaDa$pcent= 100*(MaDa$ccc-MaDa$cpc)/MaDa$ccc
 MaDa$txt6= rround(MaDa$ide,  3)
 ttu=paste("        The indirect effect of ",MaDa$xn," on ",MaDa$yn," or ab is equal to ",MaDa$txt6,sep="")
 MaDa$ides=MaDa$aar*MaDa$bbr
 tto ="less than small"
 if (abs(MaDa$ides)>.01) tto="small"
 if (abs(MaDa$ides)>.09) tto="medium"
 if (abs(MaDa$ides)>.25) tto="large"
 if (abs(MaDa$ides)>.02 & MaDa$di==1) tto="small"
 if (abs(MaDa$ides)>.15 & MaDa$di==1) tto="medium"
 if (abs(MaDa$ides)>.40 & MaDa$di==1) tto="large"
 if (tto !=" ") tto="smaller than small"
 if (MaDa$di == 1)MaDa$ess='d*r'
 if (MaDa$di == 0)MaDa$ess='r*r'
 tte= rround(MaDa$ides,  3)
 MaDa$txt7=paste( ttu,".")
 if (tto !=" ")MaDa$txt7=paste( ttu,", with a ", tto," effect size (",MaDa$ess," = ", tte,"),",sep="")
 tto= rround(MaDa$cpc,  3)
 MaDa$txt6 = paste(MaDa$txt7," and the direct effect is equal to ", tto,".",sep="")
 tmp= rround(MaDa$pcent,2)
 if (MaDa$pcent<0) tmp="zero"
 tto= rround(MaDa$ides, 3)
 ttu= rround(MaDa$cpc,  3)
 MaDa$txt3=paste(MaDa$txt6,"  The percentage of the total effect of ",rround(MaDa$ccc,  3)," (c' + ab) that is mediated is equal to ", tmp," percent.",sep="")
 aaii=0
 if ( MaDa$bbp < MaDa$alpha |  MaDa$aap < MaDa$alpha) aaii=1
 if (aaii==1 & abs( MaDa$aab)>abs(MaDa$bbbb)) MaDa$txt3=paste(MaDa$txt3,'  The mediator is said to be "proximal" (Hoyle & Kenny, 1999) in that standardized path a is greater than standardized path b.  Thus, ',MaDa$mn,' is "closer" to ',MaDa$xn," than to ",MaDa$yn,".",sep="")
 if (aaii==1 & abs( MaDa$aab)<abs(MaDa$bbbb)) MaDa$txt3=paste(MaDa$txt3,'  The mediator is said to be "distal" (Hoyle & Kenny, 1999) in that standardized path b is greater than standardized path a.  Thus, ',MaDa$mn,' is "closer" to ',MaDa$yn," than to ",MaDa$xn,".",sep="")
 sobel = sqrt((MaDa$aaa*MaDa$bse)**2+(MaDa$bbb*MaDa$ase)**2)
 s_Z= MaDa$ide/sobel
 s_p=2*(1-pnorm(s_Z))
 
 MaDa$txt4=paste(MaDa$txt3,"  The Sobel standard error is equal to ",rround(sobel,  3),",",sep="")
 tto= rround(s_Z,  3)
 ttu=pval(s_p)
 MaDa$txt7=paste(" which makes the Z test of the indirect effect equal to ", tto, ttu,".",sep="")
 MaDa$txt6=paste(MaDa$txt4,MaDa$txt7,sep="")
 if (s_p >= MaDa$alpha) MaDa$txt3=paste(MaDa$txt6,"  Because the Sobel test is not statistically significant, it is concluded that the indirect effect is not significantly different from zero.",sep="")
 if (s_p < MaDa$alpha) MaDa$txt3=paste(MaDa$txt6,"  Because the Sobel test is statistically significant, it is concluded that the indirect effect is significantly different from zero.",sep="")
 ttu=pval(MaDa$bootp)
 MaDa$txt6=paste("        The bootstrap estimated indirect effect (before bias correction) is ", rround(MaDa$bieb, 3), ttu," with a standard error of ", rround(MaDa$bootsd, 3)," (Preacher & Hayes, 2008).",sep="")
 MaDa$alci=100*(1-MaDa$alpha)
 tte= rround(MaDa$blci,  3)
 ttu= rround(MaDa$buci,  3)
 MaDa$txt6=paste(MaDa$txt6,"  The ", MaDa$alci," percent bias-corrected bootstrap confidence interval (",round(MaDa$trials, digit = 0)," trials) is from ", tte," to ", ttu,", and because zero is",sep="")
 tto=" not" 
 if (MaDa$blci*MaDa$buci < 0) tto=""
 tte=" not" 
 if (MaDa$blci*MaDa$buci > 0)tte=""
 MaDa$txt6=paste(MaDa$txt6, tto," in the confidence interval, it is concluded that the indirect effect is", tte," different from zero.",sep="")
 MaDa$txt6=paste(MaDa$txt6,"  Table 3 presents the total, direct, and indirect effects.  (In contemporary analyses, the bootstrapped test, and not the Sobel test, is reported.)",sep="")

# Assumptions 
 MaDa$txt4=paste("        For the estimates below to be valid, it must be assumed that there is no measurement error in ",MaDa$mn,sep="") 
 if (MaDa$xm1 != "TRUE")  MaDa$txt4=paste(MaDa$txt4," and ",MaDa$xn,".",sep="")
 if (MaDa$xm1 == "TRUE")  MaDa$txt4=paste(MaDa$txt4,".",sep="")
 MaDa$txt4=paste(MaDa$txt4,"  Additionally, it must be assumed that there are no unmeasured common causes of ",MaDa$mn," and ",MaDa$yn,sep="")
 if (MaDa$xm1 != "TRUE")MaDa$txt4=paste(MaDa$txt4," or of ",MaDa$xn," and ",MaDa$mn,".",sep="")
 if (MaDa$xm1 == "TRUE")MaDa$txt4=paste(MaDa$txt4,".",sep="")
 MaDa$txt4=paste(MaDa$txt4,"  It must be assumed that ",MaDa$yn," does not cause ",MaDa$mn,sep="")
 if (MaDa$xm1 != "TRUE")MaDa$txt4=paste(MaDa$txt4," and that ",MaDa$mn," does not cause ",MaDa$xn,".",sep="")
 if (MaDa$xm1 == "TRUE") MaDa$txt4=paste(MaDa$txt4,".",sep="")
 MaDa$txt4=paste(MaDa$txt4,"  Finally, it must be assumed that ",MaDa$xn," and ",MaDa$mn," do not interact to cause ",MaDa$yn,".",sep="")


# Nonlinearity and Interaction.
 if (MaDa$di+MaDa$dimm <2) MaDa$txt7 = "        The tests of nonlinearity and interaction are as follows:"
 if (MaDa$di+MaDa$dimm==2) MaDa$txt7 = "        The test of interaction is as follows:"

# XM Interaction
 tto= rround(MaDa$xmb,  3)
 tte=knt(MaDa$xmp,MaDa$alpha)
 MaDa$txt7 = paste(MaDa$txt7,"  The interactive effect of ",MaDa$xn," and ",MaDa$mn, " on ",MaDa$yn," is ", tto," and is", tte," statistically significant", pval(MaDa$xmp),".",sep="")


 if (MaDa$di==1)  MaDa$txt7 = paste(MaDa$txt7,"  Because ",MaDa$xn," is a dichotomy, its quadratic effects cannot be measured.",sep="")

if(MaDa$di<1)
{

#X squared on M 
if (is.na(MaDa$x2mb) == "FALSE") tto= rround(MaDa$x2mb,  3)
 if ( MaDa$di==0 &  tto !=".") ttu = pval(MaDa$x2mp)
 if ( MaDa$di<1 &  tto !=".")MaDa$txt7 = paste(MaDa$txt7,"  The quadratic effect of ",MaDa$xn," squared on ",MaDa$mn," is ", tto," and is", knt(MaDa$x2mp,MaDa$alpha)," statistically significant", ttu,".",sep="")
 if (MaDa$di<1 &  tto==".")  MaDa$txt7 = paste(MaDa$txt7,"  The quadratic effect of ",MaDa$xn," squared on ",MaDa$mn," cannot be estimated because of insufficient range in ",MaDa$xn,".",sep="")


#X squared on Y 
 if (is.na(MaDa$x2b) == "FALSE") tto= rround(MaDa$x2b,  3)
 if (MaDa$di==0 & tto != ".") MaDa$txt7 = paste(MaDa$txt7,"  The quadratic effect of ",MaDa$xn," squared on ",MaDa$mn," is ", tto," and is ", knt(MaDa$x2p,MaDa$alpha),"statistically significant", pval(MaDa$x2p),".",sep="")
 if (MaDa$di==0 & tto==".")  MaDa$txt7 = paste(MaDa$txt7,"  The quadratic effect of ",MaDa$xn," squared on ",MaDa$mn," cannot be estimated because of insufficient range in ",MaDa$xn,".",sep="")

}
#M squared on Y 
 tto= rround(MaDa$m2b, 3)

 if (MaDa$dimm==0 &  tto != ".")MaDa$txt7 = paste(MaDa$txt7,"  The quadratic effect of ",MaDa$mn," squared on ",MaDa$yn," is ", tto," and is", knt(MaDa$m2p,MaDa$alpha)," statistically significant", pval(MaDa$m2p),".",sep="")
 if (MaDa$dimm==0 &  tto==".")  MaDa$txt7 = paste(MaDa$txt7,"  The quadratic effect of ",MaDa$mn," squared on ",MaDa$yn," cannot be estimated because of insufficient range in ",MaDa$mn,".",sep="")
 if (MaDa$di==1) MaDa$x2p = 1

#Summary
 if (MaDa$x2p<MaDa$alpha | MaDa$x2p<MaDa$alpha | MaDa$m2p<MaDa$alpha) MaDa$txt7 = paste(MaDa$txt7,"  There are concerns about nonlinear effects and either a data transformation or a nonlinear term might be advisable.",sep="")
 if (MaDa$x2p>MaDa$alpha & MaDa$m2p>MaDa$alpha & MaDa$x2p>MaDa$alpha & MaDa$di*MaDa$dimm==0) MaDa$txt7 = paste(MaDa$txt7,"  There is then no evidence of nonlinear effects.",sep="")
 if (MaDa$x2p<MaDa$alpha) MaDa$wnum=MaDa$wnum+1
 if (MaDa$x2p<MaDa$alpha) MaDa$wrn = paste(MaDa$wrn,"  ", round(MaDa$wnum, digit = 0),".  There is evidence that the effect of  ",MaDa$xn," on ", MaDa$mn," is nonlinear and either a data transformation or a nonlinear term might be advisable.",sep="")
 MaDa$txt7 = paste(MaDa$txt7,"  The linear interactive effect of ",MaDa$xn," and ",MaDa$mn," is", knt(MaDa$xmp,MaDa$alpha)," statistically significant", pval(MaDa$xmp),sep="")
 if (MaDa$xmp< MaDa$alpha) MaDa$txt7 = paste(MaDa$txt7,"  This interaction should be added to the mediational model.",sep="")
 if (MaDa$m2p<MaDa$alpha) MaDa$wnum=MaDa$wnum+1
 if (MaDa$m2p<MaDa$alpha) MaDa$wrn = paste(MaDa$wrn,"  ", round(MaDa$wnum, digit = 0),".  There is evidence that the effect of ",MaDa$mn," on ", MaDa$yn," is nonlinear and either a data transformation or a nonlinear term might be advisable.",sep="")
 if (MaDa$x2p<MaDa$alpha) MaDa$wnum=MaDa$wnum+1
 if  (MaDa$x2p<MaDa$alpha) MaDa$wrn = paste(MaDa$wrn,"  ", round(MaDa$wnum, digit = 0),".  There is evidence that ",MaDa$xn," on ", MaDa$yn," is nonlinear and some sort of data transformation or allowance for nonlinearity is advisable.",sep="")


 txx= ifelse (MaDa$rsq3<.1,1,0) 
 if (txx==1)MaDa$wnum=MaDa$wnum+1
 tmp= round(MaDa$wnum, digit = 0)
 tte=""
 if (MaDa$df1>1) tte=" and the covariates"
 if (MaDa$df1==1) tte=" and the covariate"
 if (txx==1) MaDa$wrn = paste(MaDa$wrn,"  ", round(MaDa$wnum, digit = 0),".  The variables ",MaDa$xn," and ", MaDa$mn, tte," explain less than 1 percent of the variance of ",MaDa$yn,".",sep="")
 
#Add an additional warning for the mediator as the outcome 1 percent. 

 if (MaDa$xmp<MaDa$alpha) MaDa$wnum=MaDa$wnum+1
 if (MaDa$xmp<MaDa$alpha) MaDa$wrn = paste(MaDa$wrn,"  ", round(MaDa$wnum, digit = 0),".  The variables ",MaDa$xn," and ", MaDa$mn," interact to explain ",MaDa$yn,".",sep="")

 if (MaDa$trials<5000) {MaDa$wnum=MaDa$wnum+1
 MaDa$wrn = paste(MaDa$wrn,"  ", round(MaDa$wnum, digit = 0),".  It is recommended that there be a minimum of 5,000 trials for valid confidence intervals.",sep="")}




# Power Analysis.

 tte= round(MaDa$alpha,digit=4)
 if (MaDa$alpha==.005) tte=".005"
 if (MaDa$alpha==.05) tte=".05"
 if (MaDa$alpha==.01) tte=".01"
 if (MaDa$alpha==.02) tte=".02"
 if (MaDa$alpha==.10) tte=".10"
 ttu= round(MaDa$nxn,digit=0)
 MaDa$txt5=paste("        In this section, theoretical power analyses are computed using the study's sample size of ", ttu," with an alpha of ", tte,".",sep="")
 MaDa$txt5=paste(MaDa$txt5,"  Baron and Kenny (1986) terminology is used.",sep="")
 if (MaDa$ncov>1) MaDa$txt5=paste(MaDa$txt5,"  Power calculations include the effects of the covariates: the loss of degrees of freedom, collinearity, and reduced error variance.",sep="")
 if (MaDa$ncov==1) MaDa$txt5=paste(MaDa$txt5,"  Power calculations include the effects of the covariates: the loss of degrees of freedom, collinearity, and reduced error variance.",sep="")

# Step 1 
 r = .09
 if (MaDa$di==1) r = dtr(.5,lg*MaDa$nxn,lg*MaDa$nxn)*.3
 collin = MaDa$rsqX.C
 rsqp = MaDa$rsqY.C +r^2*(1-MaDa$rsqX.C)
 MaDa$pow= cpow(r,MaDa$nxn,MaDa$df1,MaDa$alpha,collin,rsqp)
 tto= rrround(MaDa$pow,2)
 if (tto =="1.00") tto= "virtually 1"
 if (MaDa$di==0) MaDa$txt5=paste(MaDa$txt5,"  The power of the Step 1 test is ", tto,", assuming that direct effect (path c') is zero and that paths a and b have a moderate effect size (r = .3).",sep="")
 if (MaDa$di==1) MaDa$txt5=paste(MaDa$txt5,"  The power of the Step 1 test is ", tto,", assuming that direct effect (path c') is zero and that paths a and b have a moderate effect size (d = .5 for path a and r = .3 for path b).",sep="")

 r = .3
 if (MaDa$di==1) r = dtr(.5,lg*MaDa$nxn,lg*MaDa$nxn)
 collin = MaDa$rsqX.C
 rsqp = MaDa$rsqY.C +r^2*(1-MaDa$rsqX.C)
 MaDa$pow= cpow(r,MaDa$nxn,MaDa$df1,MaDa$alpha,collin,rsqp)
 tto= rrround(MaDa$pow,2)
 if (tto =="1.00") tto= "virtually 1"
 if (MaDa$di==0)  MaDa$txt5=paste(MaDa$txt5,"  The power of the Step 1 test is ", tto,", assuming that direct effect (path c') is moderate (r = .3) and that the other paths are zero.",sep="")
 if (MaDa$di==1)  MaDa$txt5=paste(MaDa$txt5,"  The power of the Step 1 test is ", tto,", assuming that direct effect (path c') is moderate (d = .5) and that the other paths are zero.",sep="")


# Step 2
 r= .3
 if (MaDa$di==1) r = dtr(.5,lg*MaDa$nxn,lg*MaDa$nxn)
 collin = MaDa$rsqX.C
 rsqp = MaDa$rsqM.C+r^2*(1-MaDa$rsqX.C)
 MaDa$pow= cpow(r,MaDa$nxn,MaDa$df1,MaDa$alpha,collin,rsqp)
 tto= rrround(MaDa$pow,2)
 if (tto=="1.00") tto= "virtually 1"
 ttu= round(MaDa$alpha,  3)
 tte=paste(", given an alpha of ", ttu,",")
 if (MaDa$alpha==.05) tte=""
 if (MaDa$di==0)  MaDa$txt5=paste(MaDa$txt5,"  The power of the Step 2 test or a is ", tto,", assuming that effect size is moderate (r = .3).",sep="")
 if (MaDa$di==1)  MaDa$txt5=paste(MaDa$txt5,"  The power of the Step 2 test or a is ", tto,", assuming that effect size is moderate (r = .5).",sep="")

 
# Step 3 (and 4)
 r= .3
 collin = MaDa$rsqM.CX
 rsqp = MaDa$rsqY.C+r^2*(1-MaDa$rsqM.C)
 MaDa$pow= cpow(r,MaDa$nxn+1,MaDa$df1,MaDa$alpha,collin,rsqp)
 tto= rrround(MaDa$pow,2)
 if (tto=="1.00") tto= "virtually 1"
 ttu= round(MaDa$alpha,  3)
 tte=paste(", given an alpha of ", ttu,",",sep="")
 if (MaDa$alpha==.05) tte=""
 if (MaDa$ncov==0) MaDa$txt5=paste(MaDa$txt5,"  The power of the Step 3 (path b) and Step 4 (path c') tests is ", tto,sep="")
if (MaDa$ncov==0) MaDa$txt5=paste(MaDa$txt5,", assuming that the tested path has a moderate effect size (r = .3) and the other path is zero, and the correlation between ",MaDa$xn," and ",MaDa$mn," is ", rrround(MaDa$aarr,  3)," (the actual correlation between those variables).",sep="")
if (MaDa$ncov>0 & MaDa$di==0) MaDa$txt5=paste(MaDa$txt5,"  The power of the Step 3 (path b) test is ", tto,sep="")
if (MaDa$ncov>0 & MaDa$di==0) MaDa$txt5=paste(MaDa$txt5,", assuming that the path b has a moderate effect size (r = .3) and path c' is zero, and the partial correlation between ",MaDa$xn," and ",MaDa$mn," is ", rrround(MaDa$aarr, 3)," (the actual correlation between those variables).",sep="")

# Step 4
 r= .3
 if (MaDa$di==1) r = dtr(.5,lg*MaDa$nxn,lg*MaDa$nxn)
 collin = MaDa$rsqX.CM
 rsqp = MaDa$rsqY.C+r^2*(1-MaDa$rsqX.C)
 MaDa$pow= cpow(r,MaDa$nxn+1,MaDa$df1,MaDa$alpha,collin,rsqp)
 tto= rrround(MaDa$pow,2)
 if (tto=="1.00") tto= "virtually 1"
 ttu= round(MaDa$alpha,  3)
 tte=paste(", given an alpha of ", ttu,",",sep="")
 if (MaDa$alpha==.05) tte=""
 if (MaDa$ncov>0 & MaDa$di==0) MaDa$txt5=paste(MaDa$txt5,"  The power of the Step 4 (path c') test is ", tto,sep="")
 if (MaDa$ncov>0 & MaDa$di==0) MaDa$txt5=paste(MaDa$txt5,", assuming that the path c' has a moderate effect size (r = .3) and path b is zero, and the partial correlation between ",MaDa$xn," and ",MaDa$mn," is ", rround(MaDa$aarr,  3)," (the actual correlation between those variables).",sep="")

 if (MaDa$di==1) MaDa$txt5=paste(MaDa$txt5,"  The power of the Step 4 (path c') test is ", tto,sep="")
 if (MaDa$di==1) MaDa$txt5=paste(MaDa$txt5,", assuming that the path c' has a moderate effect size (d = .5) and path b is zero, and the partial correlation between ",MaDa$xn," and ",MaDa$mn," is ", rround(MaDa$aarr,  3)," (the actual correlation between those variables).",sep="")



# Indirect Effect
# r2=.09 
# if (MaDa$di==1) r2 = dtr(.5,lg*MaDa$nxn,lg*MaDa$nxn)*.3
# This needs to be compared to the old one.
# es=r2/sqrt(r2*(1-MaDa$rsqX.C)/(1-MaDa$rsqM.C-r2*(1-MaDa$rsqX.C))/(MaDa$nxn- MaDa$df1-3)+r2*(1-MaDa$rsqM.CX)/(1-MaDa$rsqY.C-r2*(1-MaDa$rsqM.C))/(MaDa$nxn- MaDa$df1-2))

ra=.3
rb=.3 
 if (MaDa$di==1) ra = dtr(.5,lg*MaDa$nxn,lg*MaDa$nxn)

 es=ra*rb/sqrt(ra^2*(1-MaDa$rsqX.C)/(1-MaDa$rsqM.C-ra^2*(1-MaDa$rsqX.C))/(MaDa$nxn- MaDa$df1-3)+rb^2*(1-MaDa$rsqM.CX)/(1-MaDa$rsqY.C-rb^2*(1-MaDa$rsqM.C))/(MaDa$nxn- MaDa$df1-2))


 MaDa$pow=1-pnorm(c(MaDa$z2T-es),mean=0,sd=1,lower.tail=TRUE)+ pnorm(c(-MaDa$z2T - es), mean=0, sd=1,lower.tail=TRUE)
 tto= rrround(MaDa$pow,2)
 if (tto=="1.00") tto= "virtually 1"
 MaDa$txt5=paste(MaDa$txt5,"  A conservative estimate of power of the test of the indirect effect",sep="")
 MaDa$txt5=paste(MaDa$txt5," is ", tto," assuming that a and b have moderate effect sizes and that the direct effect is zero.  Again, all of these power calculations are hypothetical based on the assumption of moderate effect sizes.",sep="")

# Summary Statement.

 MaDa$txt0 = "        Here is an attempt to summarize the results, but they need to be carefully verified by the investigator."
 
 tto= rround(MaDa$cpc,  3)
 ttu = pval(MaDa$cpp) 

 tte=""
 if (MaDa$cpp> MaDa$alpha) tte= " not"
 MaDa$txt0=paste(MaDa$txt0,"  The direct effect from ",MaDa$xn," to ",MaDa$yn, " equals ", tto," and is", tte," statistically significant", ttu,".",sep="")
 tta= rround(MaDa$bbb,  3)
 tte="increases"
 if (MaDa$cpc<0) tte="decreases"
 if (MaDa$di != 1) MaDa$txt0=paste(MaDa$txt0,"  As ",MaDa$xn," increases by one unit, ",MaDa$yn," ", tte," by ",tta," units.",sep="")
 if (MaDa$di==1) MaDa$txt0=paste(MaDa$txt0,"  The predicted mean difference between the ",MaDa$labx2," and ",MaDa$labx1," groups on ",MaDa$yn," equals ", tto,".",sep="")
 tto= rround(MaDa$ide,  3)

 ttu = pval(MaDa$bootp)

 MaDa$txt0=paste(MaDa$txt0,"  The indirect effect from ",MaDa$xn," to ",MaDa$yn," equals ", tto," and is", knt(MaDa$bootp,MaDa$alpha)," statistically significant", ttu,".",sep="")
 ttu="increases"
 if (MaDa$ide<0) ttu="decreases"
 tta= rround(abs(MaDa$ide),  3)
 if (MaDa$di != 1) MaDa$txt0=paste(MaDa$txt0,"  For the indirect effect, as ",MaDa$xn," increases by one unit, ",MaDa$yn," ", ttu," indirectly via ",MaDa$mn," by ",tta," units.",sep="")
 if (MaDa$di==1) MaDa$txt0=paste(MaDa$txt0,"  For the indirect effect, the predicted mean difference indirectly via ",MaDa$mn," between the ",MaDa$labx2," and ",MaDa$labx1," groups on ",MaDa$yn," equals ", tto,".",sep="")



 if (MaDa$pcent >= 80 & MaDa$bootp<MaDa$alpha & MaDa$cpp > MaDa$alpha) MaDa$txt0=paste(MaDa$txt0,"  There is evidence of complete mediation.")
 if (MaDa$pcent >= 80 & MaDa$bootp<MaDa$alpha & MaDa$cpp <MaDa$alpha) MaDa$txt0=paste(MaDa$txt0,"  There is evidence of nearly complete mediation, despite the effect from ",MaDa$xn," to ",MaDa$yn," is statistically significant .",sep="")
 if (MaDa$pcent < 80 & MaDa$bootp<MaDa$alpha) MaDa$txt0=paste(MaDa$txt0,"  There is evidence of partial mediation of the effect of ",MaDa$xn," on ",MaDa$yn," given that the indirect effect is statistically significant",sep="") 
 if (MaDa$pcent < 80 & MaDa$bootp<MaDa$alpha) MaDa$txt0=paste(MaDa$txt0," but the percentage of the total effect mediated is less than 80 percent.",sep="")
 if (MaDa$bootp>MaDa$alpha) MaDa$txt0=paste(MaDa$txt0,"  There is not evidence of mediation because the indirect effect is not statistically significant.",sep="")


 



 if ((MaDa$pcent >= 100 | MaDa$pcent <0) & (MaDa$bootp<MaDa$alpha | MaDa$cpp<MaDa$alpha)) MaDa$txt0=paste(MaDa$txt0,"  There is inconsistent mediation (MacKinnon, Fairchild, & Fritz, 2007) in that the indirect effect", ttu,") and the direct effect", tto,") have opposite signs",sep="")
 if ((MaDa$pcent >= 100 | MaDa$pcent <0) & (MaDa$bootp<MaDa$alpha | MaDa$cpp<MaDa$alpha)) MaDa$txt0=paste(MaDa$txt0,", which is akin to a suppressor effect, that makes path c small.",sep="")

 if (MaDa$wnum==0) MaDa$wrn=" "
 if (MaDa$wnum>1) MaDa$wrn=paste("WARNINGS:",MaDa$wrn,sep="")
 if (MaDa$wnum==1) MaDa$wrn=paste("WARNING:",MaDa$wrn,sep="")

tab1 <- matrix(c(
 "Variable", "Mean" ,"Standard Deviation",
 MaDa$xn,rround(MaDa$xmean,3),rround(MaDa$xsd,3),
 MaDa$mn,rround(MaDa$mmean,3),rround(MaDa$msd,3),
 MaDa$yn,rround(MaDa$ymean,3),rround(MaDa$ysd,3)),ncol=3,byrow=TRUE)

xxx1 = rround(MaDa$ccl,3)
yy1 = rround(MaDa$ccu,3)
xxx2 = rround(MaDa$aal,3)
yy2 = rround(MaDa$aau,3)
xxx3 = rround(MaDa$bbl,3)
yy3 = rround(MaDa$bbu,3)
xxx4 = rround(MaDa$cpl,3)
yy4 = rround(MaDa$cpu,3)

xxx5 = rround(MaDa$blci,3)
yy5 = rround(MaDa$buci,3)


yyy1=ifelse(MaDa$ccp<.001,"<.001",rrround(MaDa$ccp,3))
yyy2=ifelse(MaDa$aap<.001,"<.001",rrround(MaDa$aap,3))
yyy3=ifelse(MaDa$bbp<.001,"<.001",rrround(MaDa$bbp,3))
yyy4=ifelse(MaDa$cpp<.001,"<.001",rrround(MaDa$cpp,3))

yyy5=ifelse(MaDa$bootp<.001,"<.001",rrround(MaDa$bootp,3))


tta1=paste(tta1,"% CI",sep="")
tab2 <- matrix(c(
 "Step","Path","Estimate","Lower",tta1,"Upper","Beta","r","p",
" 1  "," c  ",rround(MaDa$ccc,3),xxx1,"to  ",yy1,rround(MaDa$ccb,3),rrround(MaDa$ccrr,3),yyy1,
" 2  "," a  ",rround(MaDa$aaa,3),xxx2,"to  ",yy2,rround(MaDa$aab,3),rrround(MaDa$aarr,3),yyy2,
" 3  "," b  ",rround(MaDa$bbb,3),xxx3,"to  ",yy3,rround(MaDa$bbbb,3),rrround(MaDa$bbr,3),yyy3,
" 4  "," c' ",rround(MaDa$cpc,3),xxx4,"to  ",yy4,rround(MaDa$cpb,3),rrround(MaDa$cprr,3),yyy4),
  ,ncol=9,byrow=TRUE)

tab3 <- matrix(c(
 "Effect",  "Path(s)","Estimate","Lower",tta1,"Upper","p",
 "Total",   "  c   ",rround(MaDa$ccc,3),xxx1,"to  ",yy1,yyy1,
 "Direct",  "  c'  ",rround(MaDa$cpc,3),xxx4,"to  ",yy4,yyy4,
 "Indirect"," ab   ",rround(MaDa$bie,3),xxx5,"to  ",yy5,yyy5),
  ,ncol=7,byrow=TRUE)





 
fileMT<-file(MaDa$ofilename,"w")



writeLines(MaDa$wrn, fileMT)
writeLines(" ", fileMT)
writeLines("MEDIATIONAL MODEL", fileMT)
writeLines(" ", fileMT)
writeLines(MaDa$txt1, fileMT)
writeLines(" ", fileMT)
writeLines(MaDa$txt4, fileMT)
writeLines(" ", fileMT)
writeLines(MaDa$txt9, fileMT)
writeLines(" ", fileMT)
writeLines("RESULTS", fileMT)
writeLines(" ", fileMT)
writeLines("Descriptive Statistics", fileMT)
writeLines(" ", fileMT)
writeLines(MaDa$txt8, fileMT)
writeLines(" ", fileMT)
writeLines("Power", fileMT)
writeLines(" ", fileMT)
writeLines(MaDa$txt5, fileMT)
writeLines(" ", fileMT)
writeLines("The Four Steps", fileMT)
writeLines(" ", fileMT)
writeLines(MaDa$txt2, fileMT)
writeLines(" ", fileMT)
writeLines(MaDa$txt11, fileMT)
writeLines(" ", fileMT)
writeLines("The Indirect Effect", fileMT)
writeLines(" ", fileMT)
writeLines(MaDa$txt3, fileMT)
writeLines(" ", fileMT)
writeLines(MaDa$txt6, fileMT)
writeLines(" ", fileMT)
writeLines("Tests of Nonlinearity and Interaction", fileMT)
writeLines(" ", fileMT)
writeLines(MaDa$txt7, fileMT)
writeLines(" ", fileMT)
writeLines("Overall Summary", fileMT)
writeLines(" ", fileMT)
writeLines(MaDa$txt0, fileMT)
writeLines(" ", fileMT)
writeLines("                                       References", fileMT)
writeLines(" ", fileMT)
writeLines("      Baron, R. M., & Kenny, D. A.  (1986).  The moderator-mediator variable distinction in social psychological research: Conceptual, strategic and statistical considerations. Journal of Personality and Social Psychology, 51, 1173-1182.", fileMT)
writeLines("      Hoyle, R. H., & Kenny, D. A.  (1999).  Sample size, reliability, and tests of statistical mediation.  In R. H. Hoyle (Ed.), Statistical strategies for small sample research (pp. 195-222).  Thousand Oaks, CA:  Sage." , fileMT) 
writeLines("      Huber, P. J. (1964). Robust estimation of a location parameter. Annals of Mathematical Statistics, 35, 73-101.", fileMT)
writeLines("      MacKinnon, D. P., Fairchild, A. J., & Fritz, M. S.  (2007). Mediation analysis. Annual Review of Psychology, 58, 593-614.", fileMT)
writeLines("      Preacher, K. J., & Hayes, A. F. (2008).  Asymptotic and resampling strategies for assessing and comparing indirect effects in multiple mediator models.  Behavior Research Methods, 40, 879-891.", fileMT)
writeLines(" ", fileMT)
writeLines("                                       Table 1", fileMT)
writeLines(" ", fileMT)
write.table(format(tab1, justify="right"),row.names=F, col.names=F, quote=F,fileMT)
writeLines(" ", fileMT)
writeLines("                                       Table 2", fileMT)
writeLines(" ", fileMT)
write.table(format(tab2, justify="right"),row.names=F, col.names=F, quote=F,fileMT)
writeLines(" ", fileMT)
writeLines("                                       Table 3", fileMT)
writeLines(" ", fileMT)
write.table(format(tab3, justify="right"),row.names=F, col.names=F, quote=F,fileMT)



tte=""
if  (MaDa$cpp < MaDa$alpha) tte="*"
tto=""
if  (MaDa$ccp < MaDa$alpha) tto="*"
tte1= paste(rround(MaDa$cpc,3),tte," (",rround(MaDa$ccc,3),tto,")",sep="")
tte2=paste(rround(MaDa$cpb,3),tte," (",rround(MaDa$ccb,3),tto,")",sep="")

tti=""
if  (MaDa$aap < MaDa$alpha) tti="*"
tti1 =paste(rround(MaDa$aaa,3),tti,sep="")
tti2=paste(rround(MaDa$aab,3),tti,sep="")

ttu=""
if  (MaDa$bbp < MaDa$alpha) ttu="*"
ttu1 = paste(rround(MaDa$bbb,3),ttu,sep="")
ttu2=paste(rround(MaDa$bbbb,3),ttu,sep="")

if (MaDa$wnum> 0)
cat(MaDa$wrn, sep="")
cat("", sep="\n")

 
cat("Run is now finished and the following files can now be viewed:", sep="\n")
cat("   The text file that describes and tables the mediation results: ",MaDa$ofilename, sep="")
cat("", sep="\n")
cat("   The two mediation figures: ",MaDa$F1," (unstandardized) ",MaDa$F1," (standardized)", sep="")
cat("", sep="\n")
cat("   The R output of the results: ",MaDa$rfilename,".", sep="")





img <- read.jpeg(MaDa$figgy)
png(filename =MaDa$F1, width = 1200, height = 800,   pointsize = 12, bg = "white")
plot(img)
text(400,95,MaDa$yn,cex=1.8)
text(56,95,MaDa$xn,cex=1.8)
text(253,295,MaDa$mn,cex=1.8)
text(253,100,tte1,cex=1.8) 
text(171,200,tti1,cex=1.8) 
text(371,200,ttu1,cex=1.8) 
title(main = "Mediation Diagram: Unstandardized Estimates",cex.main=2.4)
dev.off()


img1 <- read.jpeg(MaDa$figgy)
png(filename =MaDa$F2, width = 1200, height = 800,   pointsize = 12, bg = "white")

plot(img1)
text(400,95,MaDa$yn,cex=1.8)
text(56,95,MaDa$xn,cex=1.8)
text(253,295,MaDa$mn,cex=1.8)
text(253,100,tte2,cex=1.8) 
text(171,200,tti2,cex=1.8) 
text(371,200,ttu2,cex=1.8) 
title(main = "Mediation Diagram: Standardized Estimates",cex.main=2.4)
dev.off()

close(fileMT)






}




# bootstrap function using a modified version Kelley's MBESS


mediation <- function (x,mediator,dv,ncv=0,catyr="", conf.level = 0.95,  B = 1000,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10) 
{ncv1=ncv
 cog=catyr

# Here is the internal mediation function.
    .mediation <- function(x = x, mediator = mediator, dv = dv,ncv=ncv,catyr=catyr,conf.level = conf.level,c1=c1,c2=c2,c3=c3,c4=c4,c5=c5,c6=c6,c7=c7,c8=c8,c9=c9,c10=c10) 
{

ncv=ncv1


if (ncv1==0)    Data <- na.omit(cbind(x, mediator, dv))
if (ncv1==1)    Data <- na.omit(cbind(x, mediator, dv,c1))
if (ncv1==2)    Data <- na.omit(cbind(x, mediator, dv,c1,c2))
if (ncv1==3)    Data <- na.omit(cbind(x, mediator, dv,c1,c2,c3))
if (ncv1==4)    Data <- na.omit(cbind(x, mediator, dv,c1,c2,c3,c4))
if (ncv1==5)    Data <- na.omit(cbind(x, mediator, dv,c1,c2,c3,c4,c5))
if (ncv1==6)    Data <- na.omit(cbind(x, mediator, dv,c1,c2,c3,c4,c5,c6))
if (ncv1==7)    Data <- na.omit(cbind(x, mediator, dv,c1,c2,c3,c4,c5,c6,c7))
if (ncv1==8)    Data <- na.omit(cbind(x, mediator, dv,c1,c2,c3,c4,c5,c6,c7,c8))
if (ncv1==9)    Data <- na.omit(cbind(x, mediator, dv,c1,c2,c3,c4,c5,c6,c7,c8,c9))
if (ncv1==10)    Data <- na.omit(cbind(x, mediator, dv,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10))


            N <- dim(Data)[1]
Dater <- (as.data.frame(Data))
         

dog = paste("mediator~x",cog,sep="")
            M.on.X <- lm(dog,data=Dater)
M.on.X

ss2 <- summary(M.on.X)
ss <- ss2[4]
sss <-  as.data.frame(ss)
path.a=sss[2,1]

SE.M_X <-  sss[2,2]


dog = paste("dv ~ x + mediator",cog,sep="")

            Y.on.X.and.M <- lm(dog,data=Dater)
ss2 <- summary(Y.on.X.and.M)
ss <- ss2[4]
sss <-  as.data.frame(ss)
path.b=sss[3,1]


SE.Y_XM <-  sss[3,2]

                            
	
        
         
            
        ab <- path.a * path.b
        

        

        
        Indirect.Effect <- c(Estimate = ab)
        
#  Stores results of bootstapping. 
        
               if (sum(x == 0 | x == 1) != N) {
              
              
              Effect.Sizes <- rbind(
                  Indirect.Effect = Indirect.Effect)
              
                Results.mediation <- list(
                  Y.on.X = Regression.of.Y.on.X, 
                  M.on.X = Regression.of.M.on.X, 
                  Y.on.X.and.M = Regression.of.Y.on.X.and.M,
                  Effect.Sizes = Effect.Sizes)
                
                 
            }
            if (sum(x == 0 | x == 1) == N) {
              

           ES <- c(Estimate = ifelse((path.a==0 & path.b==0), 0, (((path.a * path.b)/(sqrt(SE.Y_XM^2 * path.a^2 + SE.M_X^2 * path.b^2))) * (sqrt(1/sum(x == 0) + 1/sum(x == 1))))))

                  Effect.Sizes <- rbind(
                  Indirect.Effect = Indirect.Effect   )
                
                Results.mediation <- list(
                     Effect.Sizes = Effect.Sizes)
                
                
                  
            
        }
        
        return(Results.mediation)
    }

cat("Bootstrap resampling has begun.  This process takes a considerable amount of time.  Please be patient.", sep="\n")

  
# Here is the internal bootstrap function.
.mediation.bs <- function(x,mediator,dv,ncv,catyr,conf.level = conf.level, B = B,c1=c1,c2=c2,c3=c3,c4=c4,c5=c5,c6=c6,c7=c7,c8=c8,c9=c9,c10=c10)
{
 
if (ncv1==0)Data <- na.omit(cbind(x, mediator, dv))
if (ncv1==1)Data <- na.omit(cbind(x, mediator, dv,c1))
if (ncv1==2)Data <- na.omit(cbind(x, mediator, dv,c1,c2))
if (ncv1==3)Data <- na.omit(cbind(x, mediator, dv,c1,c2,c3))
if (ncv1==4)Data <- na.omit(cbind(x, mediator, dv,c1,c2,c3,c4))
if (ncv1==5)Data <- na.omit(cbind(x, mediator, dv,c1,c2,c3,c4,c5))
if (ncv1==6)Data <- na.omit(cbind(x, mediator, dv,c1,c2,c3,c4,c5,c6))
if (ncv1==7)Data <- na.omit(cbind(x, mediator, dv,c1,c2,c3,c4,c5,c6,c7))
if (ncv1==8)Data <- na.omit(cbind(x, mediator, dv,c1,c2,c3,c4,c5,c6,c7,c8))
if (ncv1==9)Data <- na.omit(cbind(x, mediator, dv,c1,c2,c3,c4,c5,c6,c7,c8,c9))
if (ncv1==10)Data <- na.omit(cbind(x, mediator, dv,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10))


N <- dim(Data)[1]

Effect.Size.Names <- rownames(Result$Effect.Sizes)

       Boot.This <- function(Data, g) 
         {
          # Applies the internal mediation function to the data for the rows identified by g.

if (ncv1==0)  Output.mediation <- .mediation(x = Data[g, 1], mediator = Data[g, 2], dv = Data[g, 3], conf.level = conf.level)
if (ncv1==1)  Output.mediation <- .mediation(x = Data[g, 1], mediator = Data[g, 2], dv = Data[g, 3],c1 = Data[g, 4], conf.level = conf.level)
if (ncv1==2)  Output.mediation <- .mediation(x = Data[g, 1], mediator = Data[g, 2], dv = Data[g, 3],c1 = Data[g, 4],c2 = Data[g, 5], conf.level = conf.level)
if (ncv1==3)  Output.mediation <- .mediation(x = Data[g, 1], mediator = Data[g, 2], dv = Data[g, 3],c1 = Data[g, 4],c2 = Data[g, 5],c3 = Data[g, 6], conf.level = conf.level)
if (ncv1==4)  Output.mediation <- .mediation(x = Data[g, 1], mediator = Data[g, 2], dv = Data[g, 3],c1 = Data[g, 4],c2 = Data[g, 5],c3 = Data[g, 6],c4 = Data[g, 7], conf.level = conf.level)
if (ncv1==5)  Output.mediation <- .mediation(x = Data[g, 1], mediator = Data[g, 2], dv = Data[g, 3],c1 = Data[g, 4],c2 = Data[g, 5],c3 = Data[g, 6],c4 = Data[g, 7],c5 = Data[g, 8], conf.level = conf.level)
if (ncv1==6)  Output.mediation <- .mediation(x = Data[g, 1], mediator = Data[g, 2], dv = Data[g, 3],c1 = Data[g, 4],c2 = Data[g, 5],c3 = Data[g, 6],c4 = Data[g, 7],c5 = Data[g, 8],c6 = Data[g, 9], conf.level = conf.level)
if (ncv1==7)  Output.mediation <- .mediation(x = Data[g, 1], mediator = Data[g, 2], dv = Data[g, 3],c1 = Data[g, 4],c2 = Data[g, 5],c3 = Data[g, 6],c4 = Data[g, 7],c5 = Data[g, 8],c6 = Data[g, 9],c7 = Data[g, 10], conf.level = conf.level)
if (ncv1==8)  Output.mediation <- .mediation(x = Data[g, 1], mediator = Data[g, 2], dv = Data[g, 3],c1 = Data[g, 4],c2 = Data[g, 5],c3 = Data[g, 6],c4 = Data[g, 7],c5 = Data[g, 8],c6 = Data[g, 9],c7 = Data[g, 10],c8 = Data[g, 11], conf.level = conf.level)
if (ncv1==9)  Output.mediation <- .mediation(x = Data[g, 1], mediator = Data[g, 2], dv = Data[g, 3],c1 = Data[g, 4],c2 = Data[g, 5],c3 = Data[g, 6],c4 = Data[g, 7],c5 = Data[g, 8],c6 = Data[g, 9],c7 = Data[g, 10],c8 = Data[g, 11],c9 = Data[g, 12], conf.level = conf.level)
if (ncv1==10)  Output.mediation <- .mediation(x = Data[g, 1], mediator = Data[g, 2], dv = Data[g, 3],c1 = Data[g, 4],c2 = Data[g, 5],c3 = Data[g,6],c4 = Data[g, 7],c5 = Data[g, 8],c6 = Data[g, 9],c7 = Data[g, 10],c8 = Data[g, 11],c9 = Data[g, 12],c10 = Data[g, 13], conf.level = conf.level)

          as.numeric(Output.mediation$Effect.Sizes) # Returns only the values of effect size estimates
         }
        

boot.out <- boot(data = Data, statistic=Boot.This, R = B, stype = "i")


Bootstrap.Replicates <- as.data.frame(boot.out$t)
colnames(Bootstrap.Replicates) <- Effect.Size.Names

bsrep <<- as.data.frame(boot.out$t)



# write.table(Bootstrap.Replicates,file="ggMD.csv",sep=",",row.names=F)

    
Number.of.Effect.Sizes <- length(Effect.Size.Names)
Get.CIs.for.These <- c(1:Number.of.Effect.Sizes) # All estimates of effect size (this is the index for each submitted to a loop)

 

BCa.BS.Results <- matrix(NA, Number.of.Effect.Sizes, 3)
colnames(BCa.BS.Results) <- c("Estimate", "CI.Lower_BCa", "CI.Upper_BCa") 
rownames(BCa.BS.Results) <- Effect.Size.Names

for(k in 1:Number.of.Effect.Sizes)
{
  


{
BS.Results <- c(NA, NA)
BS.Results <- try(boot.ci(boot.out = boot.out, index=Get.CIs.for.These[k], conf = conf.level, type = c("bca"))$bca[4:5], silent=TRUE)
if(is.numeric(BS.Results)==FALSE) BS.Results <- c(NA, NA)
BCa.BS.Results[k, 1:3] <- c(Result$Effect.Sizes[k], BS.Results)
}
}

 BS.Output <- BCa.BS.Results

return(BS.Output) # Bootstrap statistics are output here.
return(BCa.BS.Results)
}
    
    


# Overall Results 
Result <-   .mediation(x = x, mediator = mediator, dv = dv,ncv=ncv,  conf.level = conf.level,c1=c1,c2=c2, c3=c3,c4=c4,c5=c5,c6=c6,c7=c7,c8=c8,c9=c9,c10=c10)

# Result if bootstrapping is done.
 

BS.Result <- .mediation.bs(x=x, mediator=mediator, dv=dv,   conf.level=conf.level, B = B,c1=c1,c2=c2, c3=c3,c4=c4,c5=c5,c6=c6,c7=c7,c8=c8,c9=c9,c10=c10)
Result <- list(Y.on.X=Result$Y.on.X, M.on.X=Result$M.on.X, Y.on.X.and.M=Result$Y.on.X.and.M, Bootstrap.Results=BS.Result)

return(Result)
return(BCa.BS.Results)
return(BS.Results)


}



# Create window
window = gtkWindow()
# Add title
window["title"] = "Basic Information for MedTextR"
# Add a frame
 frame = gtkFrameNew("")
window$add(frame)
# Create vertical container for file name entry
vbox = gtkVBoxNew(FALSE, 8)
vbox$setBorderWidth(24)
frame$add(vbox)
# Add horizontal container for every widget line
hbox = gtkHBoxNew(FALSE, 8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("Data file name and location:")
hbox$packStart(label,FALSE,FALSE,0)
# Add entry in the second column; named "filename"
filename = gtkEntryNew()
filename$setWidthChars(50)
label$setMnemonicWidget(filename)
filename$setText("C:/mydata.csv")
label = gtkLabel("(csv or sav file; use forward slash (/) not backslash)")
hbox$packStart(filename,FALSE,FALSE,0)
hbox$packStart(label,FALSE,FALSE,0)




hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("Variable Names in the Dataset")
hbox$packStart(label,FALSE,FALSE,0)

hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("     Causal variable:")
hbox$packStart(label,FALSE,FALSE,0)
xvar = gtkEntryNew()
xvar$setWidthChars(40)
xvar$setText("X")
label$setMnemonicWidget(xvar)
hbox$packStart(xvar,FALSE,FALSE,0)

hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("     Outcome variable:")
hbox$packStart(label,FALSE,FALSE,0)
yvar = gtkEntryNew()
yvar$setWidthChars(40)
yvar$setText("Y")
label$setMnemonicWidget(yvar)
hbox$packStart(yvar,FALSE,FALSE,0)


hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("     Mediator:")
hbox$packStart(label,FALSE,FALSE,0)
mvar = gtkEntryNew()
mvar$setWidthChars(40)
mvar$setText("M")
label$setMnemonicWidget(mvar)
hbox$packStart(mvar,FALSE,FALSE,0)

hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("     Covariates:")
hbox$packStart(label,FALSE,FALSE,0)
clis = gtkEntryNew()
clis$setWidthChars(100)
clis$setText("")
label = gtkLabel("separate by commas")
label$setMnemonicWidget(clis)
hbox$packStart(clis,FALSE,FALSE,0)
hbox$packStart(label,FALSE,FALSE,0)








hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("Variable Names in English for MedTextR")
hbox$packStart(label,FALSE,FALSE,0)


hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("     Causal Variable:")
hbox$packStart(label,FALSE,FALSE,0)
xn = gtkEntryNew()
xn$setWidthChars(40)
xn$setText("Causal Variable")
label$setMnemonicWidget(xn)
hbox$packStart(xn,FALSE,FALSE,0)



hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("     Outcome variable:")
hbox$packStart(label,FALSE,FALSE,0)
yn = gtkEntryNew()
yn$setWidthChars(40)
yn$setText("Outcome")
label$setMnemonicWidget(yn)
hbox$packStart(yn,FALSE,FALSE,0)


hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("     Mediator:")
hbox$packStart(label,FALSE,FALSE,0)
mn = gtkEntryNew()
mn$setWidthChars(40)
mn$setText("Mediator")
label$setMnemonicWidget(mn)
hbox$packStart(mn,FALSE,FALSE,0)

hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("     Covariates:")
hbox$packStart(label,FALSE,FALSE,0)
cnam = gtkEntryNew()
cnam$setWidthChars(100)
cnam$setText("")
label$setMnemonicWidget(cnam)
hbox$packStart(cnam,FALSE,FALSE,0)





hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("Causal variable manipulated (check if yes):")
hbox$packStart(label,FALSE,FALSE,0)
xm = gtkCheckButton()
xm$active <- FALSE
hbox$packStart(xm,FALSE,FALSE,0)
label$setMnemonicWidget(xm)




hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("Label for the Causal Variable if a Dichotomy")
hbox$packStart(label,FALSE,FALSE,0)


hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("     Label for level 1:")
hbox$packStart(label,FALSE,FALSE,0)
tlab1 = gtkEntryNew()
tlab1$setWidthChars(40)
tlab1$setText("Level 1")
label$setMnemonicWidget(tlab1)
hbox$packStart(tlab1,FALSE,FALSE,0)



hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("     Label for level 2:")
hbox$packStart(label,FALSE,FALSE,0)
tlab2 = gtkEntryNew()
tlab2$setWidthChars(40)
tlab2$setText("Level 2")
label$setMnemonicWidget(tlab2)
hbox$packStart(tlab2,FALSE,FALSE,0)








hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("Alpha:")
hbox$packStart(label,FALSE,FALSE,0)
alpha = gtkEntryNew()
alpha$setWidthChars(10)
alpha$setText(".05")
label$setMnemonicWidget(alpha)
hbox$packStart(alpha,FALSE,FALSE,0)

hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("Number of bootstrap trials:")
hbox$packStart(label,FALSE,FALSE,0)
ntrials = gtkEntryNew()
ntrials$setWidthChars(8)
ntrials$setText("5000")
label = gtkLabel("minimum of 500")
label$setMnemonicWidget(ntrials)
hbox$packStart(ntrials,FALSE,FALSE,0)
hbox$packStart(label,FALSE,FALSE,0)

hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("Directory to write files and place MTRfigure.png:")
hbox$packStart(label,FALSE,FALSE,0)
drt = gtkEntryNew()
drt$setWidthChars(30)
drt$setText("C:/")
label = gtkLabel("(use forward slash (/) not backslash)")
label$setMnemonicWidget(drt)
hbox$packStart(drt,FALSE,FALSE,0)
hbox$packStart(label,FALSE,FALSE,0)




hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("MedTextR output file name:")
hbox$packStart(label,FALSE,FALSE,0)
ofile = gtkEntryNew()
ofile$setWidthChars(50)
ofile$setText("MedTextR.txt")
label = gtkLabel("(do not include directory name)")
label$setMnemonicWidget(ofile)
hbox$packStart(ofile,FALSE,FALSE,0)
hbox$packStart(label,FALSE,FALSE,0)




hbox = gtkHBoxNew(FALSE,8)
vbox$packStart(hbox, FALSE, FALSE, 0)
label = gtkLabelNewWithMnemonic("R statistical summaries file name:")
hbox$packStart(label,FALSE,FALSE,0)
rfile = gtkEntryNew()
rfile$setWidthChars(50)
rfile$setText("rfile.txt")
label = gtkLabel("(do not include directory name)")
label$setMnemonicWidget(rfile)
hbox$packStart(rfile,FALSE,FALSE,0)
hbox$packStart(label,FALSE,FALSE,0)


# Add button
the.buttons = gtkHButtonBoxNew()
the.buttons$setBorderWidth(5)
vbox$add(the.buttons)
the.buttons$setLayout("spread")
the.buttons$setSpacing(40)
 buttonOK = gtkButtonNewFromStock("gtk-ok")
 gSignalConnect(buttonOK, "clicked", Medtexty)
 the.buttons$packStart(buttonOK,fill=F)
buttonCancel = gtkButtonNewFromStock("gtk-close")
gSignalConnect(buttonCancel, "clicked", window$destroy)
the.buttons$packStart(buttonCancel,fill=F)







