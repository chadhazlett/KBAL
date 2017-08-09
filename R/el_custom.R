#==========================================
# workspace stuff
#==========================================
#data(lalonde)
#lalonde$nodegr=as.numeric(lalonde$educ<=11)
#xvars=c("age","black","educ","hisp","married","re74","re75","nodegr","u74","u75")
#attach(lalonde)

#Kbal at defaults: $1806
#kbalout=kbal(D=nsw,X=lalonde[,xvars])
#D=nsw
#X=lalonde[,xvars]

#====================================================
# We'd be giving it...
#meanK2tx= apply(rbind(K2[D==1,],K2[D==1,]),2,mean)
#K2_0=K2[D==0,]
#z=t(t(K2_0)-meanK2tx)
#bal.out.pc=try(el_custom(z=z))

#===================================================

#' @export
llog_custom = function(z, eps){
  ans = z
  avoidNA = !is.na(z)
  lo = (z < eps) & avoidNA
  ans[lo] = log(eps) - 1.5 + 2 * z[lo]/eps - 0.5 * (z[lo]/eps)^2
  ans[!lo] = log(z[!lo])
  ans
}

#' @export
log_el_custom=function(par,z,eps){
  ncon=ncol(z)
  N=nrow(z)
  gamma=par
  arg = (N + t(gamma)%*%(-t(z)))
  #used to be:  arg = (1 + t(gamma)%*%(t(z)-eta_mat))

  log_el=-sum(llog_custom(z=arg,eps=eps))
  return(log_el)
}


#' A custsom built EL optimizer.
#' @description  Custom built EL optimizer.
#' @param z Matrix of constraints -- weighted sum must equal 0.
#' @param eps Used to deal with problem of taking a log near 0.
#' @export

el_custom=function(z, sumw.tol=0.05, eps=NULL){
  gam.init=rep(0, ncol(z))
  opt.gamma=optim(par=gam.init, method="BFGS", fn=log_el_custom, z=z, eps=eps, control=list(fnscale=1))
  gam.opt=opt.gamma$par

  N=nrow(z)

  if (is.null(eps)) eps=1/N

  arg_temp = (N + t(gam.opt)%*%(-t(z)))

  #just use 1/x instead instead of the derivative of the pseudo-log
  w=as.numeric(1/arg_temp)
  sum.w=sum(w)

  #scale: should sum to 1 when actually applied:
  w_scaled=w/sum.w

  if (abs(1-sum.w)<=sumw.tol){log_el=-sum(log(w_scaled))}
  if (abs(1-sum.w)>=sumw.tol){log_el=-sum(log(w_scaled))-10^4*(1+abs(1-sum.w))}

  R=list()
  R$w=w
  R$sumw=sum.w
  R$log_el=log_el
  R$el.gamma=gam.opt
  return(R)
}
