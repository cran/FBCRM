% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetFBCRM.R
\name{GetFBCRM}
\alias{GetFBCRM}
\title{Provides the optimal dose level closest to the mtd where the next cohort of patients should be allotted based on the data.}
\usage{
GetFBCRM(X, Y, Cohort, mu, p_rho, sigma, mtd, B, p_u)
}
\arguments{
\item{X}{Vector of patients allotted to each dose level.}

\item{Y}{Vector of toxicity events in each dose.}

\item{Cohort}{Number of patients within each cohort.}

\item{mu}{Prior expected toxicity probability at each dose.}

\item{p_rho}{Prior probability that two dose-toxicity probabilities will not cluster together.}

\item{sigma}{Prior standard deviation for the parameter alpha.}

\item{mtd}{Maximum Tolerated dose toxicity probability (pre defined).}

\item{B}{Number of Iterations to run for MCMC.}

\item{p_u}{Cut-off toxicity probability for first dose.}
}
\value{
A list containing (1) Design parameters and prior hyperparameters used for running the trials and (2) a posterior summary of the resuls, including the next dose to assign patients to.
}
\description{
Provides the optimal dose level closest to the mtd where the next cohort of patients should be allotted based on the data.
}
\examples{
X=c(3, 6, 3, 3, 3, 9, 15, 6)
Y=c(1, 0, 1, 0, 0, 2,  4, 5)
Cohort=3
mu=seq(0.1,0.8,0.1)
p_rho=0.9
sigma = 2
mtd = 0.3
B=2000 ##Number of iterations
p_u=0.9
Z=GetFBCRM(X, Y, Cohort, mu, p_rho, sigma, mtd, B, p_u)
Z
}
