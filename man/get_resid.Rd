% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NetGAM_Func.R
\name{get_resid}
\alias{get_resid}
\title{# Name: get_resid}
\usage{
get_resid(data, cols, cors)
}
\arguments{
\item{data}{Dataframe with species abundance values and time-series predictor columns.}

\item{cols}{Columns of data that will be modeled using the GAM.}

\item{cors}{Correlation structure class.}
}
\value{
: Vector of GAM residuals.
}
\description{
# Description: Gets the residuals of GAM model
}
