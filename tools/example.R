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
kbalout=kbal(D=nsw,X=lalonde[,xvars], method="ebal")
summary(lm(re78~nsw,w=kbalout$w))

#Kbal with mean balance ensured first, at defaults
kbalout_mean=kbal_meanfirst(D=nsw,X=lalonde[,xvars], sigma=5)
summary(lm(re78~nsw,w=kbalout_mean$w))

#Plot both:
plot(kbalout$dist.record[1:40], pch=16,
     ylab="L1 imbalance", xlab="Num. dims of K balanced")
points(kbalout_mean$dist.record[1:40], col=2, pch=16)
legend("topright", col=c(1,2), pch=16, legend=c("full kbal","kbal after mean"))
