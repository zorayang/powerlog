f = function(x){expit(x)*expit(-x)}
expit = function(x){1/(1+exp(-x))}
logit = function(p){log(p/(1-p))}

get_samp_size = function(i0, i1, i2, b1, b2 = NULL, rho = NULL,
                         targetpwr = 0.80, alpha = 0.05, reg = "uni"){
  switch(reg,
         uni = {
           sd = sqrt(i0/(i2*i0-i1^2))
           za = qnorm(1 - alpha/2)
           zb = qnorm(targetpwr)
           n = (sd*(za + zb))^2/(b1^2)
         },
         multi = {
           ed = i1^2-i0*i2
           num = b2^2*(ed) + 2*b1*b2*(ed)*rho+b1^2*(i1^2*rho^2-i0*i2*rho^2 + i0^2*(rho^2-1))
           denom = i0*(-ed)*(b1^2+b2^2+2*b1*b2*rho)*(rho^2-1)
           var_n = num/denom
           ncp <- uniroot(function(ncp){
             targetpwr-pchisq(qchisq(1-alpha, df=1), ncp=ncp, df=1, lower=FALSE)
           }, c(1E-4,10000))$root
           n = ncp*var_n/(b1^2)
         }
  )
  return(n)
}

calc_samp_size = function(b0, b1, targetpwr, alpha, dist = "normal",
                    mu = 0, sigma = 1,
                    lwr = 0, upr = 1,
                    ctr = 0, scl = 1, df = 10,
                    succ = 1, fail = 0, p = 0.5)
{
  switch(dist,
         normal = {
           i0 = integrate(function(x){    f(b0+b1*x)*dnorm(x, mu, sigma)}, -Inf,Inf)$value
           i1 = integrate(function(x){  x*f(b0+b1*x)*dnorm(x, mu, sigma)}, -Inf,Inf)$value
           i2 = integrate(function(x){x^2*f(b0+b1*x)*dnorm(x, mu, sigma)}, (mu-5*sigma), (mu+5*sigma))$value
         },
         unif = {
           i0 = integrate(function(x){    f(b0+b1*x)*dunif(x, lwr, upr)}, -Inf,Inf)$value
           i1 = integrate(function(x){  x*f(b0+b1*x)*dunif(x, lwr, upr)}, -Inf,Inf)$value
           i2 = integrate(function(x){x^2*f(b0+b1*x)*dunif(x, lwr, upr)}, -Inf,Inf)$value
         },
         t = { # replace dt with density of scaled t
           i0 = integrate(function(x){    f(b0+b1*x)*dt.scaled(x, mean=ctr, sd = scl, df = df)}, -Inf,Inf)$value
           i1 = integrate(function(x){  x*f(b0+b1*x)*dt.scaled(x, mean=ctr, sd = scl, df = df)}, -Inf,Inf)$value
           i2 = integrate(function(x){x^2*f(b0+b1*x)*dt.scaled(x, mean=ctr, sd = scl, df = df)}, (ctr-5*scl^2*df/(df-2)), (ctr+5*scl^2*df/(df-2)))$value
         },
         binary = {
           i0 = p*f(b0+b1*(succ)) + (1-p)*f(b0+b1*(fail))
           i1 = p*(succ)*f(b0+b1*(succ)) + (1-p)*(fail)*f(b0+b1*(fail))
           i2 = p*(succ^2)*f(b0+b1*(succ)) + (1-p)*(fail^2)*f(b0+b1*(fail))
         }
  )
  return(get_samp_size(i0=i0, i1=i1, i2=i2, b1=b1, b2 = NULL, rho = NULL, targetpwr=targetpwr, alpha=alpha, reg = "uni"))
}


approx_samp_size = function(b0, b1, b2, rho, targetpwr, alpha, dist = "normal", mu=0, sigma=1){
  s = sqrt(b1^2+b2^2+2*b1*b2*rho)
  i0 = integrate(function(x){    f(b0+s*x)*dnorm(x,0,1)}, (mu-5*sigma), (mu+5*sigma))$value
  i1 = integrate(function(x){  x*f(b0+s*x)*dnorm(x,0,1)}, (mu-5*sigma), (mu+5*sigma))$value
  i2 = integrate(function(x){x^2*f(b0+s*x)*dnorm(x,0,1)}, (mu-5*sigma), (mu+5*sigma))$value
  return(get_samp_size(i0=i0, i1=i1, i2=i2, b1=b1, b2=b2, rho=rho, targetpwr=targetpwr, alpha=alpha, reg = "multi"))
}


calc_samp_size_emp = function(xmat, b0, b1, b2 = NULL, targetpwr = 0.8, alpha = 0.05, reg = "uni"){
  n = nrow(xmat)
  switch (reg,
          uni = {
            p = expit(xmat%*%c(b0, b1))
            rho = NULL
          },
          multi = {
            p = expit(xmat%*%c(b0, b1, b2))
            rho = cor(xmat[,2], xmat[,3])
          }
  )
  i0 = mean(p*(1-p))
  i1 = mean(p*(1-p)*xmat[,2])
  i2 = mean(p*(1-p)*(xmat[,2])^2)
  return(get_samp_size(i0=i0, i1=i1, i2=i2, b1=b1, b2=b2, rho=rho, targetpwr=targetpwr, alpha=alpha, reg=reg))
}


