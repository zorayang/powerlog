f = function(x){expit(x)*expit(-x)}
expit = function(x){1/(1+exp(-x))}
logit = function(p){log(p/(1-p))}

get_pwr = function(i0, i1, i2, b1, b2 = NULL, rho = NULL, n, alpha, reg){
  switch(reg,
         uni = {
           sd = sqrt(i0/(i2*i0-i1^2))
           za = qnorm(1 - alpha/2)
           zb = sqrt(n*b1^2)/sd-za
           pwr = pnorm(zb)
         },
         multi = {
           ed = i1^2-i0*i2
           num = b2^2*(ed) + 2*b1*b2*(ed)*rho+b1^2*(i1^2*rho^2-i0*i2*rho^2 + i0^2*(rho^2-1))
           denom = i0*(-ed)*(b1^2+b2^2+2*b1*b2*rho)*(rho^2-1)
           var_n = num/denom
           pwr = pchisq(qchisq((1-alpha), df=1), ncp= n*b1^2/var_n, df=1, lower=FALSE)
         }
  )
  return(pwr)
}

calc_pwr = function(b0, b1, n, alpha =.05, dist = "normal",
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
  return(get_pwr(i0=i0, i1=i1, i2=i2, b1=b1, b2=NULL, rho=NULL, n=n, alpha=alpha, reg = "uni"))
}



approx_pwr = function(b0, b1, b2, rho, n, alpha =.05, dist = "normal", mu = 0, sigma = 1){
  s = sqrt(b1^2+b2^2+2*b1*b2*rho)
  i0 = integrate(function(x){    f(b0+s*x)*dnorm(x, mu, sigma)}, -Inf,Inf)$value
  i1 = integrate(function(x){  x*f(b0+s*x)*dnorm(x, mu, sigma)}, -Inf,Inf)$value
  i2 = integrate(function(x){x^2*f(b0+s*x)*dnorm(x, mu, sigma)}, (mu-5*sigma), (mu+5*sigma))$value
  return(get_pwr(i0=i0, i1=i1, i2=i2, b1=b1, b2=b2, rho=rho, n=n, alpha=alpha, reg="multi"))
}


calc_pwr_emp = function(b0, b1, b2 = NULL, xmat, alpha = .05, reg = "uni"){
    n = nrow(xmat)
    switch (reg,
            uni = {
              rho = NULL
              p = expit(xmat%*%c(b0, b1))
            },
            multi = {
              rho = cor(xmat[,2], xmat[,3])
              p = expit(xmat%*%c(b0, b1, b2))
            }
    )
    i0 = mean(p*(1-p))
    i1 = mean(p*(1-p)*xmat[,2])
    i2 = mean(p*(1-p)*(xmat[,2])^2)
    return(get_pwr(i0=i0, i1=i1, i2=i2, b1=b1, b2=b2, rho=rho, n=n, alpha=alpha, reg=reg))
}
