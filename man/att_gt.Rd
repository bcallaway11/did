% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/did.R
\name{att_gt}
\alias{att_gt}
\title{att_gt}
\usage{
att_gt(outcome, data, tname, aggte = TRUE, w = NULL, idname = NULL,
  first.treat.name, alp = 0.05, method = "logit", bstrap = FALSE,
  biters = 1000, clustervars = NULL, cband = FALSE, seedvec = NULL,
  pl = FALSE, cores = 2, printdetails = TRUE, maxe = NULL,
  nevertreated = T)
}
\arguments{
\item{outcome}{The outcome y (in quotations, always!)}

\item{data}{The name of the data.frame that contains the data}

\item{tname}{The name of the column containing the time periods}

\item{aggte}{boolean for whether or not to compute aggregate treatment effect parameters, default TRUE}

\item{w}{The name of the column containing the weights}

\item{idname}{The individual (cross-sectional unit) id name}

\item{first.treat.name}{The name of the variable in \code{data} that contains the first
period when a particular observation is treated.  This should be a positive
number for all observations in treated groups.  It should be 0 for observations
in the untreated group.}

\item{alp}{the significance level, default is 0.05}

\item{method}{The method for estimating the propensity score when covariates
are included}

\item{bstrap}{Boolean for whether or not to compute standard errors using
the multiplier boostrap.  If standard errors are clustered, then one
must set \code{bstrap=TRUE}.}

\item{biters}{The number of boostrap iterations to use.  The default is 100,
and this is only applicable if \code{bstrap=TRUE}.}

\item{clustervars}{A vector of variables to cluster on.  At most, there
can be two variables (otherwise will throw an error) and one of these
must be the same as idname which allows for clustering at the individual
level.}

\item{cband}{Boolean for whether or not to compute a uniform confidence
band that covers all of the group-time average treatment effects
with fixed probability \code{1-alp}.  The default is \code{FALSE}
and the resulting standard errors will be pointwise.}

\item{seedvec}{Optional value to set random seed; can possibly be used
in conjunction with bootstrapping standard errors#' (not implemented)}

\item{pl}{Boolean for whether or not to use parallel processing}

\item{cores}{The number of cores to use for parallel processing}

\item{printdetails}{Boolean for showing detailed results or not}

\item{maxe}{maximum values of periods ahead to be computed in event study}

\item{nevertreated}{Boolean for using the group which is never treated in the sample as the comparison unit. Default is TRUE.}
}
\value{
\code{MP} object
}
\description{
\code{att_gt} computes the ATT in the case where there are more
 than two periods of data and allowing for treatment to occur at different points in time
 extending the method of Abadie (2005).  This method relies on once individuals are treated
 they remain in the treated state for the duration.
}
\references{
Callaway, Brantly and Sant'Anna, Pedro.  "Difference-in-Differences with Multiple Time Periods and an Application on the Minimum Wage and Employment." Working Paper <https://ssrn.com/abstract=3148250> (2018).
}