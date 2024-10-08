% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{Simulate_Enroll}
\alias{Simulate_Enroll}
\title{Survival data simulation}
\usage{
Simulate_Enroll(
  n,
  lambda,
  event,
  M,
  group,
  maxt,
  accrual,
  censor,
  followup,
  partition = "Even"
)
}
\arguments{
\item{n}{Integer. Sample size of patients}

\item{lambda}{Numerical range 0 and 1. Hazard rate of expoential distribution}

\item{event}{Numerical range 0 and 1. Event rate}

\item{M}{Integer. Number of trials generated for multiple streams of MCMC}

\item{group}{Integer. Number of subgroup for patient enrollment}

\item{maxt}{Numerical. The maximum time length of entire trial}

\item{accrual}{Numerical. The duration of patient enrolment}

\item{censor}{Numerical range 0 and 1. The censoring rate of patients leaving before trial ends.}

\item{followup}{Integer. The time length of follow up.}

\item{partition}{Logical. If \verb{partition == "Even", the trial recruits equal numbers of patients in each stage; and if }partition == "Uneven", the trial recruits unequal numbers of patients in each stage.}
}
\value{
Simulated survival data contain both survival time, censoring status, and enrollment time.
}
\description{
\code{Simulate_Enroll} generates multiple streams of data sets with survival time, censoring status, and enrollment time.
}
\examples{
data <- Simulate_Enroll(n = c(50,20,20), lambda = 0.03, event = 0.1, M = 3, group = 3, maxt = 5,
                  accrual = 3,  censor = 0.9, followup = 2,partition = "Uneven")
head(data[[1]])
head(data[[2]])
head(data[[3]])
}
