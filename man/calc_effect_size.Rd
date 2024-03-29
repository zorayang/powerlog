% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/powerlog.R
\name{calc_effect_size}
\alias{calc_effect_size}
\title{calc_effect_size: find effect size satisfying a target power and sample size in a logistic model}
\usage{
calc_effect_size(n, b2 = NULL, rho = NULL, prop, targetpwr = 0.8,
  alpha = 0.05, reg = "uni", dist = "normal", lower = -10,
  upper = 10, ...)
}
\arguments{
\item{n}{sample size}

\item{b2}{hypothesized adjustment variable effect size, if reg == "multi"}

\item{rho}{hypothesized correlations between variable of interest X1 and adjument variable X2}

\item{prop}{A number between 0 and 1 that represents case/success/outcome==1 proportion}

\item{targetpwr}{target power, default = 0.80}

\item{alpha}{significance level}

\item{reg}{either "uni" for univaraite calculations or "multi" for bivariate calculations}

\item{dist}{supposed data distribution}

\item{lower}{lower bound of the hypothesized effect size b1}

\item{upper}{upper bound of the hypothesized effect size b1}

\item{mu}{the mean of supposed normal distribution}

\item{sigma}{the standard deviation of supposed normal distribution}

\item{lwr}{the lower bound of supposed uniform distribution}

\item{upr}{the upper bound of supposed uniform distribution}

\item{ctr}{the center of supposed t distribution}

\item{scl}{the scale of supposed t distribution}

\item{df}{the degree of freedom of supposed t distribution}

\item{succ}{the value for "success" in a supposed binary distribution, usually 1}

\item{fail}{the value for "failure" in a supposed binary distribution, usually 0}

\item{p}{the proportion of "successes" in a supposed binary distribution}
}
\value{
calculated effect size from
sample size \code{n}, proportion \code{prop}, target power \code{targetpwr} and significance level \code{alpha}
in a data distrbution \code{dist} with potential distribution parameters.
If multivariate, also use the hypothesized adjustment variable effect size \code{b2},
and correlation \code{rho} between variable of interest X1 and adjustment variable X2.
Seaching from \code{lower} to \code{upper} of the hypothesized effect size of X1.
}
\description{
calc_effect_size: find effect size satisfying a target power and sample size in a logistic model
}
\examples{
findb0(b1 = 0.33, prop = 0.5, dist = "t", ctr = 3, scl = 6, df = 9)
}
