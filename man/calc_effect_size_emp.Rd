% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/powerlog.R
\name{calc_effect_size_emp}
\alias{calc_effect_size_emp}
\title{calc_effect_size_emp: find effect size satisfying a target power and sample size in a logistic model given design matrix}
\usage{
calc_effect_size_emp(xmat, prop, targetpwr, alpha = 0.05, reg = "uni",
  b2 = NULL, lower = -10, upper = 10)
}
\arguments{
\item{xmat}{design matrix}

\item{prop}{A number between 0 and 1 that represents case/success/outcome==1 proportion}

\item{targetpwr}{target power, default = 0.80}

\item{alpha}{significance level}

\item{reg}{either "uni" for univaraite calculations or "multi" for bivariate calculations}

\item{b2}{hypothesized adjustment variable effect size, if reg == "multi"}

\item{lower}{lower bound of the hypothesized effect size b1}

\item{upper}{upper bound of the hypothesized effect size b1}
}
\value{
calculated effect size from
design matrix \code{xmat}, proportion \code{prop}, target power \code{targetpwr} and significance level \code{alpha}.
If multivariate, also use the hypothesized adjustment variable effect size \code{b2},
and correlation \code{rho} between variable of interest X1 and adjustment variable X2.
Seaching from \code{lower} to \code{upper} of the hypothesized effect size of X1.
}
\description{
calc_effect_size_emp: find effect size satisfying a target power and sample size in a logistic model given design matrix
}
\examples{
findb0(b1 = 0.33, prop = 0.5, dist = "t", ctr = 3, scl = 6, df = 9)
}
