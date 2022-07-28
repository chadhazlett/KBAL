rm(list=ls())
library(KBAL)
#Run Lalonde example as in paper:
data(lalonde)
lalonde$nodegr=as.numeric(lalonde$educ<=11)
xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
attach(lalonde)

#Raw diff-in-means: way off, -$15205
mean(re78[nsw==1])-mean(re78[nsw==0])

#OLS with covariates:
summary(lm(re78~nsw+., data=lalonde[,xvars]))

#Kbal at new defaults:
kbalout= kbal(allx=lalonde[,xvars],treatment=lalonde$nsw,
                ebal.tol=1e-6, printprogress =TRUE)
summary(lm(re78~nsw,w=kbalout$w, data = lalonde))

# Examine bias bound due to remaining imbalances; 
#note that using KRLS, gamma = c'Kc for the related regression is approximately 55
kbalout$biasbound.orig*sd(re78)*sqrt(55)
kbalout$biasbound.opt*sd(re78)*sqrt(55)
plot(x=kbalout$dist.record[1,],y=sd(re78)*sqrt(55)*kbalout$dist.record[2,],ylab="biasbound", xlab="Num. dims of K balanced", pch=16)
