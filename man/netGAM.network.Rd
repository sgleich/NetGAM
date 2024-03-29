% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NetGAM_Func.R
\name{netGAM.network}
\alias{netGAM.network}
\title{Name: netGAM.network}
\usage{
netGAM.network(gam_df, method = "glasso", pvalue = NULL)
}
\arguments{
\item{gam_df}{GAM-transformed species abundance dataframe with samples as rows species as columns (i.e. output of netGAM.df function)}

\item{method}{Networking method to use (default is glasso). "glasso" = graphical lasso network constructed with the "batch.pulsar" function in the pulsar package with StARS selection; "scc" = spearman correlation network constructed with the "corr.test" function in the psych package; "pcc" = pearson correlation network constructed with the "corr.test" function in the psych package.}

\item{pvalue}{P-value cutoff for deciding whether or not an edge exists (default is NULL). P-values in corrleation networks are bonferroni-adjusted prior to declaring cutoff. P-value only needed for scc and pcc networks.}
}
\value{
Adjacency matrix of network predicitons (1 = edge, 0 = no edge) for glasso networks and an edgelist with p-values for correlation networks
}
\description{
Description: Takes a GAM-transformed species abundance dataframe as input (output of netGAM.df function), runs a network analysis on the gamm residuals, and returns an adjacency matrix of network-predicted associations.
}
