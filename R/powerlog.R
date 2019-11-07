#' logit: the logit function
#' @param p a proportion between 0 and 1
#' @return log(p/(1-p))
#' @noRd
logit = function(p){log(p/(1-p))}

#' expit: the expit function
#' @param x a real number
#' @return 1/(1+exp(-x))
#' @noRd
expit = function(x){1/(1+exp(-x))}

#' f: the function defined by the S-B method
#' @param x a real number
#' @return expit(x)*expit(-x)
#' @noRd
f = function(x){expit(x)*expit(-x)}

#' propn: find the proportion in a normal distribution
#' @param b0 the intercept term in a logistic regression model
#' @param b1 the effect size in a logistic regression model
#' @param mu the mean of a normal distribution
#' @param sigma the standard deviation of a normal distribution
#' @param multi univariate if == F, multivariate if == T
#' @return the proportion under a normal distribution
#' @noRd
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

#' propu: find the proportion in a uniform distribution
#' @param b0 the intercept term in a logistic regression model
#' @param b1 the effect size in a logistic regression model
#' @param lwr the lower bound of a uniform distribution
#' @param upr the upper bound of a uniform distribution
#' @return the proportion under a uniform distribution
#' @noRd
propu = function(b0, b1, lwr, upr){
  integrate(function(x){
    expit(b0 + b1*x) * dunif(x, lwr, upr)
  }, lwr, upr)$value
}

#' propt: find the proportion in a scaled and shifted t-distribution
#' @param b0 the intercept term in a logistic regression model
#' @param b1 the effect size in a logistic regression model
#' @param ctr the center of target t-distribution
#' @param scl the scale factor of target t-distribution
#' @param df the degrees of freedom of target t-distribution
#' @return the proportion under a scaled and shifted t-distribution
#' @noRd
propt = function(b0, b1, ctr, scl, df){
  sd = scl^2 * df/(df-2)
  integrate(function(x){
    expit(b0 + b1*x) * dt.scaled(x, mean = ctr, sd = scl, df = df)
  }, (ctr-5*sd), (ctr+5*sd))$value
}

#' propb: find the proportion in a binary distribution
#' @param b0 the intercept term in a logistic regression model
#' @param b1 the effect size in a logistic regression model
#' @param succ a numerical value representing "success" in a binary distribution
#' @param fail a numerical value representing "failure" in a binary distribution
#' @param p the proportion of "successes" in target binary distribution
#' @return the proportion under a binary distribution
#' @noRd
propb = function(b0, b1, succ, fail, p){
  expit(b0+b1*(succ))*p + expit(b0+b1*(fail))*(1-p)
}

#' propemp: find the proportion in an emperical distribution
#' @param betavec
#' @param xmat
#' @noRd
propemp = function(betavec, xmat){
  n = nrow(xmat)
  sum(expit(as.vector(xmat %*% betavec)))/n
}

#' findb0emp: find the intercept in a logistic model from an emperical distribution
#' @param b1 hypothesized effect size
#' @param b2 hypothesized adjustment variable effect size
#' @param prop A number between 0 and 1 that represents case/success/outcome==1 proportion
#' @param xmat design matrix
#' @param reg either "uni" for univaraite calculations or "multi" for bivariate calculations
#' @param lower lower bound of the hypothesized effect size b1
#' @param upper upper bound of the hypothesized effect size b1
#' @return calculated intercept value in the logistic model given
#' the hypothesized effect sizes \code{b1} and possibly \code{b2},
#' proportion \code{prop} and design matrix \code{xmat},
#' with \code{reg} specifying the a univariate or multivariate regression model,
#' seaching from \code{lower} to \code{upper} of the hypothesized effect size \code{b1}.
#' @examples
#' findb0emp(b1 = 0.33, b2 = 0.13, prop = 0.5, xmat = design_matrix, reg = "multi", lower = 0.03, upper = 0.83)
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

#' findb0: find the intercept in a logistic model
#' @param b1 hypothesized effect size
#' @param prop A number between 0 and 1 that represents case/success/outcome==1 proportion
#' @param dist supposed data distribution
#' @param b2 hypothesized adjustment variable effect size, if multi
#' @param multi univariate if == F, multivariate if == T
#' @param mu the mean of supposed normal distribution
#' @param sigma the standard deviation of supposed normal distribution
#' @param lwr the lower bound of supposed uniform distribution
#' @param upr the upper bound of supposed uniform distribution
#' @param ctr the center of supposed t distribution
#' @param scl the scale of supposed t distribution
#' @param df the degree of freedom of supposed t distribution
#' @param succ the value for "success" in a supposed binary distribution, usually 1
#' @param fail the value for "failure" in a supposed binary distribution, usually 0
#' @param p the proportion of "successes" in a supposed binary distribution
#' @return calculated intercept value in the logistic model given
#' the effect size \code{b1}, proportion \code{prop} and distribution parameters \code{mu}, \code{sigma} or so forth.
#' @examples
#' findb0(b1 = 0.33, prop = 0.5, dist = "t", ctr = 3, scl = 6, df = 9)
findb0 = function(b1, prop, dist = "normal", b2 = NULL, multi = F,
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


##############################################################
##############################################################
### -------------- Effect Size Calculations -------------- ###
##############################################################
##############################################################


#' calc_effect_size: find effect size satisfying a target power and sample size in a logistic model
#' @param n sample size
#' @param b2 hypothesized adjustment variable effect size, if reg == "multi"
#' @param rho hypothesized correlations between variable of interest X1 and adjument variable X2
#' @param prop A number between 0 and 1 that represents case/success/outcome==1 proportion
#' @param targetpwr target power, default = 0.80
#' @param alpha significance level
#' @param dist supposed data distribution
#' @param mu the mean of supposed normal distribution
#' @param sigma the standard deviation of supposed normal distribution
#' @param lwr the lower bound of supposed uniform distribution
#' @param upr the upper bound of supposed uniform distribution
#' @param ctr the center of supposed t distribution
#' @param scl the scale of supposed t distribution
#' @param df the degree of freedom of supposed t distribution
#' @param succ the value for "success" in a supposed binary distribution, usually 1
#' @param fail the value for "failure" in a supposed binary distribution, usually 0
#' @param p the proportion of "successes" in a supposed binary distribution
#' @param reg either "uni" for univaraite calculations or "multi" for bivariate calculations
#' @param lower lower bound of the hypothesized effect size b1
#' @param upper upper bound of the hypothesized effect size b1
#' @return calculated effect size from
#' sample size \code{n}, proportion \code{prop}, target power \code{targetpwr} and significance level \code{alpha}
#' in a data distrbution \code{dist} with potential distribution parameters.
#' If multivariate, also use the hypothesized adjustment variable effect size \code{b2},
#' and correlation \code{rho} between variable of interest X1 and adjustment variable X2.
#' Seaching from \code{lower} to \code{upper} of the hypothesized effect size of X1.
#' @examples
#' findb0(b1 = 0.33, prop = 0.5, dist = "t", ctr = 3, scl = 6, df = 9)
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

#' calc_effect_size_emp: find effect size satisfying a target power and sample size in a logistic model given design matrix
#' @param xmat design matrix
#' @param prop A number between 0 and 1 that represents case/success/outcome==1 proportion
#' @param targetpwr target power, default = 0.80
#' @param alpha significance level
#' @param reg either "uni" for univaraite calculations or "multi" for bivariate calculations
#' @param b2 hypothesized adjustment variable effect size, if reg == "multi"
#' @param lower lower bound of the hypothesized effect size b1
#' @param upper upper bound of the hypothesized effect size b1
#' @return calculated effect size from
#' design matrix \code{xmat}, proportion \code{prop}, target power \code{targetpwr} and significance level \code{alpha}.
#' If multivariate, also use the hypothesized adjustment variable effect size \code{b2},
#' and correlation \code{rho} between variable of interest X1 and adjustment variable X2.
#' Seaching from \code{lower} to \code{upper} of the hypothesized effect size of X1.
#' @examples
#' findb0(b1 = 0.33, prop = 0.5, dist = "t", ctr = 3, scl = 6, df = 9)
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


########################################################
########################################################
### -------------- Power Calculations -------------- ###
########################################################
########################################################


#' get_pwr: approximate power calculations given information matrix
#' @param i0 I11 in the information matrix of esimated b1
#' @param i1 I12 and I21 in the information matrix
#' @param i2 I22 in the information matrix
#' @param b1 hypothesized effect size
#' @param b2 hypothesized adjustment variable effect size, if reg == "multi"
#' @param rho hypothesized correlations between variable of interest X1 and adjustment variable X2
#' @param n current sample size
#' @param alpha significance level
#' @param reg either "uni" for univaraite calculations or "multi" for bivariate calculations
#' @return calculated power given information matrix
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


#' calc_pwr: approximate power given data distribution and its parameters
#' @param b0 hypothesized intercept of the logistic regression model
#' @param b1 hypothesized effect size of variable of interest X1
#' @param n current sample size
#' @param alpha significance level
#' @param dist supposed data distribution
#' @param mu the mean of supposed normal distribution
#' @param sigma the standard deviation of supposed normal distribution
#' @param lwr the lower bound of supposed uniform distribution
#' @param upr the upper bound of supposed uniform distribution
#' @param ctr the center of supposed t distribution
#' @param scl the scale of supposed t distribution
#' @param df the degree of freedom of supposed t distribution
#' @param succ the value for "success" in a supposed binary distribution, usually 1
#' @param fail the value for "failure" in a supposed binary distribution, usually 0
#' @param p the proportion of "successes" in a supposed binary distribution
#' @return calculated power
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

#' approx_pwr: approximate power calculations
#' @param b0 hypothesized intercept of the logistic regression model
#' @param b1 hypothesized effect size of variable of interest X1
#' @param b2 hypothesized effect size of adjustment variable X2
#' @param rho hypothesized correlations between variable of interest X1 and adjustment variable X2
#' @param n current sample size
#' @param alpha significance level
#' @param dist supposed data distribution. only supports "normal" at the moment.
#' @param mu mean of the supposed normal distribution of X1
#' @param sigma standard deviation of the supposed normal distribution of X1
#' @return approximate power
approx_pwr = function(b0, b1, b2, rho, n, alpha =.05, dist = "normal", mu = 0, sigma = 1){
  s = sqrt(b1^2+b2^2+2*b1*b2*rho)
  i0 = integrate(function(x){    f(b0+s*x)*dnorm(x, mu, sigma)}, -Inf,Inf)$value
  i1 = integrate(function(x){  x*f(b0+s*x)*dnorm(x, mu, sigma)}, -Inf,Inf)$value
  i2 = integrate(function(x){x^2*f(b0+s*x)*dnorm(x, mu, sigma)}, (mu-5*sigma), (mu+5*sigma))$value
  return(get_pwr(i0=i0, i1=i1, i2=i2, b1=b1, b2=b2, rho=rho, n=n, alpha=alpha, reg="multi"))
}

#' calc_pwr_emp: approximate power from design matrix
#' @param b0 hypothesized intercept of the logistic regression model
#' @param b1 hypothesized effect size of variable of interest X1
#' @param b2 hypothesized effect size of adjustment variable X2
#' @param xmat design matrix
#' @param alpha significance level
#' @param reg either "uni" for univaraite calculations or "multi" for bivariate calculations
#' @return approximate power from design matrix
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

##############################################################
##############################################################
### -------------- Sample Size Calculations -------------- ###
##############################################################
##############################################################

#' get_samp_size: approximate sample size
#' @param i0 I11 in the information matrix of esimated b1
#' @param i1 I12 and I21 in the information matrix
#' @param i2 I22 in the information matrix
#' @param b1 hypothesized effect size between 0 to 10000
#' @param b2 hypothesized adjustment variable effect size, if reg == "multi"
#' @param rho hypothesized correlations between variable of interest and adjument variable.
#' @param targetpwr target power, default to 0.80
#' @param alpha significance level
#' @param reg either "uni" for univaraite calculations or "multi" for bivariate calculations
#' @return calculated power given parameters
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

#' calc_samp_size: approximate sample size given data distribution and its parameters
#' @param b0 hypothesized intercept of the logistic regression model
#' @param b1 hypothesized effect size of variable of interest X1
#' @param targetpwr power target, by default = 0.80
#' @param alpha significance level
#' @param dist supposed data distribution
#' @param mu the mean of supposed normal distribution
#' @param sigma the standard deviation of supposed normal distribution
#' @param lwr the lower bound of supposed uniform distribution
#' @param upr the upper bound of supposed uniform distribution
#' @param ctr the center of supposed t distribution
#' @param scl the scale of supposed t distribution
#' @param df the degree of freedom of supposed t distribution
#' @param succ the value for "success" in a supposed binary distribution, usually 1
#' @param fail the value for "failure" in a supposed binary distribution, usually 0
#' @param p the proportion of "successes" in a supposed binary distribution
#' @return approximate sample size
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

#' approx_samp_size: approximate sample size given multivariate-normal distribution and its parameters
#' @param b0 hypothesized intercept of the logistic regression model
#' @param b1 hypothesized effect size of variable of interest X1
#' @param b2 hypothesized effect size of adjustment variable X2
#' @param rho hypothesized correlations between variable of interest X1 and adjustment variable X2
#' @param targetpwr power target, by default = 0.80
#' @param alpha significance level
#' @param dist supposed data distribution. only supports "normal" at the moment.
#' @param mu mean of the supposed normal distribution of X1
#' @param sigma standard deviation of the supposed normal distribution of X1
#' @return approximate sample size
approx_samp_size = function(b0, b1, b2, rho, targetpwr, alpha, dist = "normal", mu=0, sigma=1){
  s = sqrt(b1^2+b2^2+2*b1*b2*rho)
  i0 = integrate(function(x){    f(b0+s*x)*dnorm(x,0,1)}, (mu-5*sigma), (mu+5*sigma))$value
  i1 = integrate(function(x){  x*f(b0+s*x)*dnorm(x,0,1)}, (mu-5*sigma), (mu+5*sigma))$value
  i2 = integrate(function(x){x^2*f(b0+s*x)*dnorm(x,0,1)}, (mu-5*sigma), (mu+5*sigma))$value
  return(get_samp_size(i0=i0, i1=i1, i2=i2, b1=b1, b2=b2, rho=rho, targetpwr=targetpwr, alpha=alpha, reg = "multi"))
}

#' calc_samp_size_emp: approximate sample size from design matrix
#' @param xmat design matrix
#' @param b0 hypothesized intercept of the logistic regression model
#' @param b1 hypothesized effect size of variable of interest X1
#' @param b2 hypothesized effect size of adjustment variable X2
#' @param targetpwr power target, by default = 0.80
#' @param alpha significance level
#' @param reg either "uni" for univaraite calculations or "multi" for bivariate calculations
#' @return approximate sample size from design matrix
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
