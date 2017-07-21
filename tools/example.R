#Run Lalonde example as in paper:
data(lalonde)
lalonde$nodegr=as.numeric(lalonde$educ<=11)
xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
attach(lalonde)

#Raw diff-in-means: way off, -$15205
mean(re78[nsw==1])-mean(re78[nsw==0])

#OLS with covariates:
summary(lm(re78~nsw+., data=lalonde[,xvars]))

#Kbal at defaults: $1806
kbalout=kbal(D=nsw,X=lalonde[,xvars])
summary(lm(re78~nsw,w=kbalout$w))
plot(x=seq(2:41),kbalout$dist.record[2:41],
ylab="L1 imbalance", xlab="Num. dims of K balanced")


#Kbal with mean balance ensured first, at defaults: $1806
kbalout_mean=kbal_meanfirst(D=nsw,X=lalonde[,xvars])
summary(lm(re78~nsw,w=kbalout_mean$w))
plot(x=seq(2:41),kbalout_mean$dist.record[2:41],
     ylab="L1 imbalance", xlab="Num. dims of K balanced")
