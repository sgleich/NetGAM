# NetGAM Method Validation
# Make mock ASV/OTU data with time-series features (i.e. monthly trend and long-term trend)
# By: Samantha J. Gleich, Jacob A. Cram, Jake L. Weissman, and David A. Caron
# Last updated: October 16, 2021

### Libraries ###
library(Rlab)
library(MASS)
library(igraph)
library(Matrix)
library(tidyverse)
library(compositions)
library(igraph)
library(SpiecEasi)

### Create ASV df with underlying ASV-ASV relationships ###
# Set interaction strength, number of ASVs, and probability of interaction
interaction_strength <- 100
num_ASVs <- 400
prob <- 0.01

# Simulate ASV-ASV interactions (output is a covariance matrix)

# 1.) Barabasi-Albert method
ASVxASV <- sample_pa(400,directed=FALSE)
ASVxASV <- get.adjacency(ASVxASV)
ASVxASV<- as.matrix(ASVxASV)
diag(ASVxASV) <- prob
ASVxASV[ASVxASV==1] <- interaction_strength

# 2.) Erdos-Renyi method
# ASVxASV <- matrix(rbern(num_ASVs^2,prob),nrow=num_ASVs)*interaction_strength
# diag(ASVxASV) <- 0.01

# 3.) Network topology based on real data
# data(amgut1.filt)
# depths <- rowSums(amgut1.filt)
# amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
# amgut1.filt.cs <- round(amgut1.filt.n * min(depths))
# e <- ncol(amgut1.filt.cs)
# d <- e
# n <- nrow(amgut1.filt.cs)
# graph <- SpiecEasi::make_graph('cluster', d, e)
# ASVxASV  <- graph2prec(graph)
# ASVxASV[ASVxASV !=0] <- interaction_strength
# diag(ASVxASV) <- prob

#Simulate mean ASV abundances from a normal distribution
mu_vec <- rnorm(num_ASVs, mean=10,sd=1)
sigma_mat <- ASVxASV

# Pull abundance data for all ASVs across 200 samples from a multivariate normal distribution
sim_samples <- mvrnorm(n=200, mu=mu_vec, sigma_mat,tol=1)
# Make all values greater than or equal to 0
sim_samples[sim_samples < 0] <- 0

# Get edgelist from the covariance matrix - this is the edgelist of true ASV-ASV associations
graph <- graph.adjacency(sigma_mat,mode="undirected")
edgelist <- data.frame(get.edgelist(graph))
edgelist <- edgelist%>%group_by(val_1 = pmin(X1, X2), val_2 = pmax(X1,X2))%>%distinct(.keep_all=TRUE)
edgelist$X3 <- paste(edgelist$val_1,edgelist$val_2,sep="_")
edgelist <- edgelist%>%distinct(X3, .keep_all = TRUE)
edgelist$val_1 <- NULL
edgelist$val_2<-NULL
colnames(edgelist)<- c("e1","e2","e3")

### Add monthly time-series signals ###

# Add monthly trend to mock ASV data
# Specify what percent of the ASVs in the df will have gradual and abrupt seasonal signals
sim_samples <- as.data.frame(sim_samples)
percent_szn_abrupt <- 0.2
percent_szn_gradual <- 0.2
remove_abrupt <- sample(1:ncol(sim_samples), percent_szn_abrupt*ncol(sim_samples),replace=FALSE)
no_szn <- subset(sim_samples, select=-c(remove_abrupt))
szn_abrupt <- subset(sim_samples, select=c(remove_abrupt))
remove_gradual <- sample(1:ncol(no_szn), percent_szn_gradual*ncol(no_szn),replace=FALSE)
szn_gradual <- subset(no_szn, select=c(remove_gradual))
no_szn <- subset(no_szn, select=-c(remove_gradual))

no_szn_namez <- colnames(no_szn)
gradual_namez <- colnames(szn_gradual)
abrupt_namez <- colnames(szn_abrupt)

# Simulate abrupt monthly trend
szn_names <- colnames(szn_abrupt)
sim_new_abrupt <- NULL
szn_data_abrupt <- NULL
for (ii in 1:ncol(szn_abrupt)){
  v <- data.frame(szn_abrupt[,ii])
  num_start <- round(runif(1,0,11))
  tot <- (nrow(szn_abrupt)-1)+num_start
  num_list <- num_start:tot
  for (jj in 1:nrow(v)){
    szn_data_abrupt[jj]<- ((cos(num_list[jj]*2*pi/(12))/2) + 0.5)^10}
  sim_new_new <- v[,1]*szn_data_abrupt
  sim_new_abrupt <- cbind(sim_new_abrupt,sim_new_new)}

colnames(sim_new_abrupt)<- c(szn_names)
sim_new_abrupt <- as.data.frame(sim_new_abrupt)

# Simulate gradual monthly trend
szn_names <- colnames(szn_gradual)
sim_new_gradual <- NULL
szn_data_gradual <- NULL
for (ii in 1:ncol(szn_gradual)){
  v <- data.frame(szn_gradual[,ii])
  num_start <- round(runif(1,0,11))
  tot <- (nrow(szn_gradual)-1)+num_start
  num_list <- num_start:tot
  for (jj in 1:nrow(v)){
    szn_data_gradual[jj]<- ((cos(num_list[jj]*2*pi/(12))/2) + 0.5)^1}
  sim_new_new <- v[,1]*szn_data_gradual
  sim_new_gradual <- cbind(sim_new_gradual,sim_new_new)}

colnames(sim_new_gradual)<- c(szn_names)
sim_new_gradual <- as.data.frame(sim_new_gradual)

# Combine ASVs with gradual monthly trend, abrupt monthly trend, and no monthly trend
sim_szn_fin <- cbind(sim_new_gradual,sim_new_abrupt,no_szn)

### Add long-term time-series trend ###

# Add long-term trend to mock ASV data
# Specify what percent of the ASVs in the df will have a long-term trend
percent_time <- 0.5
remove <- sample(1:ncol(sim_szn_fin), percent_time*ncol(sim_szn_fin),replace=FALSE)
no_trend <- subset(sim_szn_fin, select=-c(remove))
trend <- subset(sim_szn_fin, select=c(remove))

# Add long term trend to ASVs in the "trend" df
namez <- colnames(trend)
for (timez in 1:ncol(trend)){
  coef <- rnorm(1,0.01,0.0009)
  line <- seq(from=0.01, by=coef,length.out=200)
  line2 <- rev(line)
  if (timez>100){trend[,timez] <- trend[,timez]+line}
  if (timez < 100){trend[,timez] <- trend[,timez]+line2}}
colnames(trend) <- namez
sim_time_fin <- cbind(trend,no_trend)

### Add time-series variables and remove samples ###
# Add time-series variables (i.e. month of year and month of time-series)
numz <- rep(1:12,100)
sim_time_fin$MonthOfYear <- numz[1:200]
sim_time_fin$MonthCount <- 1:200
df.new <- sim_time_fin

# Remove time-series variables from df temporarily
df.new2 <- subset(df.new,select=-c(MonthOfYear,MonthCount))

### Make ASV abundance data resemble sequence data ###

# Exponentiate mock ASV data
df_out <- exp(df.new2)

# Calculate relative abundance of mock ASVs
df_out <- as.data.frame(df_out)
df_out$sum <- rowSums(df_out)
df_out_rel <- df_out/df_out$sum
df_out_rel$sum <- NULL
df_out_rel <- as.matrix(df_out_rel)

# CLR transform mock ASV data
output.asvs<- as.data.frame(df_out_rel)

### Finalize mock ASV df and save ###

# Add back time-series variables to df (i.e. month of year and month of time-series)
# output.asvs$MonthOfYear <- df.new$MonthOfYear
# output.asvs$MonthCount <- df.new$MonthCount

# Save mock ASV df to use the GAMGlasso network analysis
write.csv(output.asvs,"output_asvs.csv")
