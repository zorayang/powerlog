f = function(x){expit(x)*expit(-x)}
expit = function(x){1/(1+exp(-x))}
logit = function(p){log(p/(1-p))}

propn = function(b0, b1, b2=NULL, mu, sigma, multi = F){
  if(multi)
  {
    integrate(function(x){
      expit(b0 + b1*x + b2*x) * dnorm(x, mean = mu, sd = sigma)
      }, mu-5*sigma, mu+5*sigma)$value

  }else
    integrate(function(x){
      expit(b0 + b1*x) * dnorm(x, mean = mu, sd = sigma)
      }, mu-5*sigma, mu+5*sigma)$value
}

propu = function(b0, b1, lwr, upr){
  integrate(function(x){
    expit(b0 + b1*x) * dunif(x, lwr, upr)
    }, lwr, upr)$value
}

propt = function(b0, b1, ctr, scl, df){
  sd = scl^2 * df/(df-2)
  integrate(function(x){
    expit(b0 + b1*x) * dt.scaled(x, mean = ctr, sd = scl, df = df)
    }, (ctr-5*sd), (ctr+5*sd))$value
}

propb = function(b0, b1, succ, fail, p){
  expit(b0+b1*(succ))*p + expit(b0+b1*(fail))*(1-p)
}

propemp = function(betavec, xmat){
  n = nrow(xmat)
  sum(expit(as.vector(xmat %*% betavec)))/n
}


findb0emp = function(b1, b2 = NULL, prop, xmat, reg = "uni", lower, upper){
  switch (reg,
          "uni"= {
            result = uniroot(function(b0){
              propemp(c(b0,b1), xmat)-prop},
              interval = c(lower, upper), extendInt = "yes")
          },

          "multi" = {
            result = uniroot(function(b0){
              propemp(c(b0,b1,b2), xmat)-prop},
              interval = c(lower, upper), extendInt = "yes")
          }
  )
  return(result$root)
}


findb0 = function(b1, prop, dist = "normal", b2=NULL, multi=F,
                  mu = 0, sigma = 1,
                  lwr = 0, upr = 1,
                  ctr = 0, scl = 1, df = 10,
                  succ = 1, fail = 0, p = 0.5){

  switch(dist,
         normal = {
           if(multi){
             result = uniroot(function(b0){
               propn(b0=b0, b1=b1, b2=b2, multi=T, mu=mu, sigma=sigma)-prop},
               interval = c(mu-5*sigma, mu+5*sigma), extendInt = "yes")
           }else
             result = uniroot(function(b0){
               propn(b0=b0, b1=b1, mu=mu, sigma=sigma)-prop},
               interval = c(mu-5*sigma, mu+5*sigma), extendInt = "yes")
         },

         unif = {
           result = uniroot(function(b0){
             propu(b0, b1, lwr, upr)-prop},
             interval = c(lwr, upr), extendInt = "yes")
         },

         t = {
           result = uniroot(function(b0){
             propt(b0, b1, ctr, scl, df)-prop},
             interval = c(ctr-5*scl^2*df/(df-2), ctr+5*scl^2*df/(df-2)), extendInt = "yes")
         },

         binary = {
           result = uniroot(function(b0){
             propb(b0, b1, succ, fail, p)-prop},
             interval = c(fail, succ), extendInt = "yes")
         }
  )
  return(result$root)
}

calc_effect_size = function(n, b2=NULL, rho=NULL, prop, targetpwr=0.80, alpha=0.05,
                            reg = "uni", dist = "normal", lower=-10, upper=10, ...){
  switch(reg,
         uni = {
           powerb1 = function(b1){
             b0 = findb0(b1=b1, b2=NULL, prop=prop, dist=dist, multi=F, ...)
             return(calc_pwr(b0=b0, b1=b1, n=n, alpha=alpha, dist=dist, ...))
           }
           b1 = uniroot(function(b1){powerb1(b1)-targetpwr}, c(lower, upper), extendInt="yes")$root
         },

         multi = {
           powerb1 = function(b1){
             b0 = findb0(b1=b1, b2=b2, prop=prop, dist="normal", multi = T, ...)
             return(approx_pwr(b0=b0, b1=b1, b2=b2, rho=rho, n=n, alpha=alpha, dist=dist))
           }
           b1 = uniroot(function(b1){powerb1(b1)-targetpwr}, c(lower, upper), extendInt="yes")$root
         }
  )
  return(b1)
}


calc_effect_size_emp = function(xmat, prop, targetpwr, alpha = 0.05,
                                reg = "uni", b2 = NULL, lower=-10, upper=10){
  switch(reg,
         uni = {
           powerb1 = function(b1){
             b0 = findb0emp(b1=b1, prop=prop, xmat=xmat, reg = "uni", lower=lower, upper=upper)
             return(calc_pwr_emp(b0=b0, b1=b1, b2=NULL, xmat=xmat, alpha=alpha, reg="uni"))
           }
           b1 = uniroot(function(b1){powerb1(b1)-targetpwr}, c(lower, upper), extendInt = "yes")$root
         },

         multi = {
           powerb1 = function(b1){
             b0 = findb0emp(b1=b1, b2=b2, prop=prop, xmat=xmat, reg = "multi", lower=lower, upper=upper)
             return(calc_pwr_emp(b0=b0, b1=b1, b2=b2, xmat=xmat, alpha=alpha, reg = "multi"))
           }
           b1 = uniroot(function(b1){powerb1(b1)-targetpwr}, c(lower, upper), extendInt = "yes")$root
         }
  )
  return(b1)
}



