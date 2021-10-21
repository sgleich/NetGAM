#' Name: works_fun
#'
#' Description: Checks to see if the GAM can be fit to each species abundance vector in the input dataframe
#'
#' @param data Dataframe with species abundance values and time-series predictor columns.
#' @param cols Columns of data that will be modeled using the GAM.
#' @param cors Correlation structure class.
#' @return T/F list of whether errors popped up during GAM fitting
works_fun <- function(data,cols,cors) {
  options(warn = -1)
  formula.tmp<- paste(cols, '~ s(MonthOfYear, bs="cc", k=4) + s(MonthCount,bs="ts",k=4)')
  formula.tmp <- formula(formula.tmp)
  cors.tmp <- paste(cors)
  cors.tmp <- match.fun(cors.tmp)
  x <- try(mgcv::gamm(formula.tmp, data = data, correlation=cors.tmp(form = ~MonthCount),knots=list(MonthOfYear = c(1,12))))
  works <- inherits(x,"try-error")
  options(warn = 0)
  return(works)
}

#' Name: ac_remove_fun
#'
#' Description: Fits a GAM to each species and then checks to see if significant AC is still prevalent in GAM residuals
#'
#' @param data Dataframe with species abundance values and time-series predictor columns.
#' @param cols Columns of data that will be modeled using the GAM.
#' @param cors Correlation structure class.
#' @return List of whether significant AC was detected in the GAM residuals of each species.
ac_remove_fun <- function(data,cols,cors,cutoff){
  options(warn = -1)
  cutoff_neg <- cutoff*-1
  formula.tmp<- paste(cols, '~ s(MonthOfYear, bs="cc", k=4) + s(MonthCount,bs="ts",k=4)')
  formula.tmp <- formula(formula.tmp)
  cors.tmp <- paste(cors)
  cors.tmp <- match.fun(cors.tmp)
  m <- mgcv::gamm(formula.tmp, data = data, correlation = cors.tmp(form = ~MonthCount),knots=list(MonthOfYear = c(1,12)))
  res <- stats::resid(m$lme, type = "normalized")
  myacf<-stats::pacf(res, plot=FALSE, lag.max = 1)
  myacf2 <- (myacf$acf)
  sum <- sum(myacf2< cutoff_neg | myacf2 > cutoff)
  options(warn = 0)
  return(sum)}

#' Name: get_resid
#'
#' Description: Gets the residuals of GAM model
#'
#' @param data Dataframe with species abundance values and time-series predictor columns.
#' @param cols Columns of data that will be modeled using the GAM.
#' @param cors Correlation structure class.
#' @return: Vector of GAM residuals.
get_resid <- function(data,cols,cors) {
  options(warn = -1)
  formula.tmp<- paste(cols, '~ s(MonthOfYear, bs="cc", k=4) + s(MonthCount,bs="ts",k=4)')
  formula.tmp <- formula(formula.tmp)
  cors.tmp <- paste(cors)
  cors.tmp <- match.fun(cors.tmp)
  m <- mgcv::gamm(formula.tmp, data = data, correlation = cors.tmp(form = ~MonthCount),knots=list(MonthOfYear = c(1,12)))
  tmp <- stats::resid(m$gam)
  options(warn =0)
  return(tmp)}


#' Name: netGAM.df
#'
#' Description: Fits a gamm to each species in a species abundance dataset and returns a dataframe of the residuals of each gamm. The returned dataframe can be thought of as a species abundance dataframe with a reduced influence of time on the species abundance values. This dataframe can be used as input in downstream network analyses.
#'
#' @param df Species abundance dataframe with samples as rows species as columns.
#' @param MOY Vector that specifies the month of year for each sample (row) in the dataframe (i.e. 1-12).
#' @param MCount Vector that specifies the day of the time-series for each sample (row) in the dataframe (e.g. 1-200).
#' @param clrt If TRUE, the clr transformation in the compositions package is used to clr transform the input species abundance dataframe prior to GAM transformation (default is TRUE). If FALSE, the species abundance data are not clr transformed prior to GAM transformation.
#' @return GAM-transformed dataframe (rows = samples, columns = species)
#'
#' @export
  netGAM.df <- function(df, MOY, MCount, clrt=TRUE){
  #' Get dataframe set up properly. Do CLR transformation if clr == TRUE
  MOY <- as.data.frame(MOY)
  MCount <- as.data.frame(MCount)
  colnames(MOY)<- "MonthOfYear"
  colnames(MCount)<- "MonthCount"

  if (clrt==TRUE){
    df.clr <- compositions::clr(df)
  }

  if(clrt==FALSE){
    df.clr <- df
  }

  df.pred <- cbind(df.clr,MOY,MCount)
  df.pred <- as.data.frame(df.pred)

  # Get column names (species names)
  n <- ncol(df.pred)
  col_names <- colnames(df.pred[1:(n-2)])
  # Use the works_fun function to see if the GAM can be resolved for the species in the df
  works_fun_out <- lapply(c(col_names), works_fun, data=df.pred,cors="corCAR1")
  works_fun_out <- data.frame(unlist(works_fun_out, use.names=FALSE))
  works_fun_out$cols <- col_names
  colnames(works_fun_out) <- c("TF","cols")

  # Separate species for which GAM worked (keep) from species for which GAM did not work (try_again)
  keep_list <- subset(works_fun_out, TF=="FALSE")
  keep_list_names <- keep_list$cols
  keep <- subset(df.pred, select=c(keep_list_names,"MonthOfYear","MonthCount"))

  remove_list <- subset(works_fun_out, TF=="TRUE")
  remove_list_names <- remove_list$cols
  try_again<- subset(df.pred, select=c(remove_list_names))
  try_again$MonthOfYear <- df.pred$MonthOfYear
  try_again$MonthCount <- df.pred$MonthCount


  # We have 2 predictor columns in the df (MonthOfYear, MonthCount)
  n <- ncol(keep)-2
  dates <- data.frame(rownames(keep))
  names <- data.frame(colnames(keep))
  col_names <- colnames(keep[1:n])

  # AC cutoff values
  cut <- 1.96/sqrt(nrow(keep)-1)

  # Use ac_remove_fun to see if significant AC is removed with the GAM
  ac_remove_fun_out <- lapply(c(col_names), ac_remove_fun, data=keep,cors="corCAR1",cutoff=cut)
  ac_remove_fun_out <- data.frame(unlist(ac_remove_fun_out, use.names=FALSE))
  ac_remove_fun_out$cols <- col_names
  colnames(ac_remove_fun_out) <- c("AC","cols")

  # If there are species with sinificant AC in the GAM residuals, remove them and add them to the try_again df
  if(sum(ac_remove_fun_out$AC==1)>=1){
    remove <- subset(ac_remove_fun_out, AC==1)
    remove <- subset(keep,select=c(remove$cols))
    try_again <- cbind(remove,try_again)
    remove <- subset(ac_remove_fun_out, AC==1)
    nums <- which(colnames(keep) %in% remove$cols)
    keep <- subset(keep, select=-c(nums))}

  # We have 2 predictor columns in the df (MonthOfYear, MonthCount)
  n <- ncol(try_again)-2
  dates <- data.frame(rownames(try_again))
  names <- data.frame(colnames(try_again))

  # Try corAR1 cor structure
  if(n>0){
    col_names <- colnames(try_again[1:n])
    workz_fun_out <- lapply(c(col_names), works_fun, data=try_again, cors="corAR1")
    works_fun_out <- data.frame(unlist(works_fun_out, use.names=FALSE))
    works_fun_out$cols <- col_names
    colnames(works_fun_out) <- c("TF","cols")
    not_work <- sum(works_fun_out$TF=="TRUE")
    dates <- data.frame(rownames(try_again))
    names <- data.frame(colnames(try_again))

    # See if corAR1 removes significant AC
    cut <- 1.96/sqrt(nrow(try_again)-1)

    if(not_work==0){
      col_names <- colnames(try_again[1:n])
      ac_remove_fun_out <- lapply(c(col_names), ac_remove_fun, data=try_again,cors="corAR1",cutoff=cut)
      ac_remove_fun_out <- data.frame(unlist(ac_remove_fun_out, use.names=FALSE))
      ac_remove_fun_out$cols <- col_names
      colnames(ac_remove_fun_out) <- c("AC","cols")
      AR1_AC <- sum(ac_remove_fun_out$AC==1)}
    if (not_work>0){AR1_AC <- NA}}

  # Try corCompSymm cor structure
  # We have 2 predictor columns in the df (MonthOfYear, MonthCount)
  n <- ncol(try_again)-2
  if(n>0){
    col_names <- colnames(try_again[1:n])
    works_fun_out <- lapply(c(col_names), works_fun, data=try_again,cors="corCompSymm")
    works_fun_out <- data.frame(unlist(works_fun_out, use.names=FALSE))
    works_fun_out$cols <- col_names
    colnames(works_fun_out) <- c("TF","cols")
    not_work <- sum(works_fun_out$TF=="TRUE")
    dates <- data.frame(rownames(try_again))
    names <- data.frame(colnames(try_again))

    # AC cutoff values
    cut <- 1.96/sqrt(nrow(try_again)-1)

    # See if corCompSymm removes significant AC
    if(not_work==0){
      col_names <- colnames(try_again[1:n])
      ac_remove_fun_out <- lapply(c(col_names), ac_remove_fun, data=try_again,cors="corCompSymm",cutoff=cut)
      ac_remove_fun_out <- data.frame(unlist(ac_remove_fun_out, use.names=FALSE))
      ac_remove_fun_out$cols <- col_names
      colnames(ac_remove_fun_out) <- c("AC","cols")
      CompSymm_AC <- sum(ac_remove_fun_out$AC==1)}
    if (not_work>0){CompSymm_AC <- NA}}

  # Try corExp cor structure
  # We have 2 predictor columns in the df (MonthOfYear, MonthCount)
  n <- ncol(try_again)-2

  if(n>0) {
    col_names <- colnames(try_again[1:n])
    workz_fun_out <- lapply(c(col_names), works_fun, data=try_again,cors="corExp")
    works_fun_out <- data.frame(unlist(works_fun_out, use.names=FALSE))
    works_fun_out$cols <- col_names
    colnames(works_fun_out) <- c("TF","colz")
    not_work <- sum(works_fun_out$TF=="TRUE")
    dates <- data.frame(rownames(try_again))
    names <- data.frame(colnames(try_again))

    # AC cutoff values
    cut <- 1.96/sqrt(nrow(try_again)-1)

    # See if corExp removes significant AC
    if(not_work==0){
      col_names <- colnames(try_again[1:n])
      ac_remove_fun_out <- lapply(c(col_names), ac_remove_fun, data=try_again,cors="corExp",cutoff=cut)
      ac_remove_fun_out <- data.frame(unlist(ac_remove_fun_out, use.names=FALSE))
      ac_remove_fun_out$cols <- col_names
      colnames(ac_remove_fun_out) <- c("AC","colz")
      Exp_AC <- sum(ac_remove_fun_out$AC==1)}
    if (not_work>0){Exp_AC <- NA}}

  # Try corGaus cor structure
  # We have 2 predictor columns in the df (MonthOfYear, MonthCount)
  n <- ncol(try_again)-2

  if(n>0){
    col_names <- colnames(try_again[1:n])
    works_fun_out <- lapply(c(col_names), works_fun, data=try_again,cors="corGaus")
    works_fun_out <- data.frame(unlist(works_fun_out, use.names=FALSE))
    works_fun_out$cols <- col_names
    colnames(works_fun_out) <- c("TF","colz")
    not_work <- sum(works_fun_out$TF=="TRUE")
    dates <- data.frame(rownames(try_again))
    names <- data.frame(colnames(try_again))

    # AC cutoff values
    cut <- 1.96/sqrt(nrow(try_again)-1)

    # See if corGaus removes significant AC
    if(not_work==0){
      col_names <- colnames(try_again[1:n])
      ac_remove_fun_out <- lapply(c(col_names), ac_remove_fun, data=try_again,cors="corGaus",cutoff=cut)
      ac_remove_fun_out <- data.frame(unlist(ac_remove_fun_out, use.names=FALSE))
      ac_remove_fun_out$cols <- col_names
      colnames(ac_remove_fun_out) <- c("AC","cols")
      Gaus_AC <- sum(ac_remove_fun_out$AC==1)}

    if (not_work>0){Gaus_AC <- NA}}

  # Determine which GAM equation will be used for the species that could not be GAM resolved and/or that had signficant AC left in corCAR1 model residuals
  if(n>0){
    x <- data.frame(a=c(AR1_AC,CompSymm_AC, Exp_AC,Gaus_AC))
    x$b <- c("corAR1","corCompSymm","corExp","corGaus")
    x <- subset(x, !is.na(a))
    x <- dplyr::arrange(x,a)
    x <- as.data.frame(x)
    method <- x[1,2]}

  # Get residuals from keep df for network using the corCAR1 model
  # We have 2 predictor columns in the df (MonthOfYear, MonthCount)
  n <- ncol(keep)-2
  col_names <- colnames(keep[1:n])
  dates <- data.frame(rownames(keep))
  names <- data.frame(colnames(keep))

  output.asvs.fin <- lapply(c(col_names), get_resid, data=keep,cors="corCAR1")
  output.asvs.fin <- data.frame(matrix(unlist(output.asvs.fin), nrow=length(output.asvs.fin), byrow=TRUE))
  output.asvs.fin <- as.data.frame(t(output.asvs.fin))
  colnames(output.asvs.fin)<-names$colnames.keep.[1:n]
  rownames(output.asvs.fin)<-dates$rownames.keep.

  # Get residuals from try_again df for network using the best GAM cor class determined above
  # We have 2 predictor columns in the df (MonthOfYear, MonthCount)
  n <- ncol(try_again)-2
  col_names <- colnames(try_again[1:n])
  if(n>0){
    dates <- data.frame(rownames(try_again))
    names <- data.frame(colnames(try_again))
    output.asvs.fin2 <- lapply(c(col_names), get_resid, data=try_again,cors=method)
    output.asvs.fin2 <- data.frame(matrix(unlist(output.asvs.fin2), nrow=length(output.asvs.fin2), byrow=TRUE))
    output.asvs.fin2 <- as.data.frame(t(output.asvs.fin2))
    colnames(output.asvs.fin2)<-col_names
    rownames(output.asvs.fin2)<-dates$rownames.try_again.

    output.asvs.gam <- cbind(output.asvs.fin,output.asvs.fin2)}
  if(n==0){output.asvs.gam <- output.asvs.fin}
  return(output.asvs.gam)}


#' Name: netGAM.network
#'
#' Description: Fits a gamm to each species in a species abundance dataset, extracts the gamm residuals for each species, runs a network analysis on the gamm residuals, and returns an adjacency matrix of network-predicted associations.
#'
#' @param df Species abundance dataframe with samples as rows species as columns.
#' @param MOY Vector that specifies the month of year for each sample (row) in the dataframe (i.e. 1-12).
#' @param MCount Vector that specifies the day of the time-series for each sample (row) in the dataframe (e.g. 1-200).
#' @param clrt If TRUE, the clr transformation in the compositions package is used to clr transform the input species abundance dataframe prior to GAM transformation (default is TRUE). If FALSE, the species abundance data are not clr transformed prior to GAM transformation.
#' @param method Networking method to use (default is glasso). "glasso" = graphical lasso network constructed with the "batch.pulsar" function in the pulsar package with StARS selection; "scc" = spearman correlation network constructed with the "corr.test" function in the psych package; "pcc" = pearson correlation network constructed with the "corr.test" function in the psych package.
#' @param pvalue P-value cutoff for deciding whether or not an edge exists (default is NULL). P-values in corrleation networks are bonferroni-adjusted prior to declaring cutoff. P-value only needed for scc and pcc networks.
#' @return Adjacency matrix of network predicitons (1 = edge, 0 = no edge)
#'
#' @export
netGAM.network <- function(df,MOY,MCount,clrt=TRUE,method="glasso",pvalue=NULL){
  #' Get dataframe set up properly. Do CLR transformation if clr == TRUE
  MOY <- as.data.frame(MOY)
  MCount <- as.data.frame(MCount)
  colnames(MOY)<- "MonthOfYear"
  colnames(MCount)<- "MonthCount"

  if (clrt==TRUE){
    df.clr <- compositions::clr(df)
  }

  if(clrt==FALSE){
    df.clr <- df
  }

  df.pred <- cbind(df.clr,MOY,MCount)
  df.pred <- as.data.frame(df.pred)

  # Get column names (species names)
  n <- ncol(df.pred)
  col_names <- colnames(df.pred[1:(n-2)])

  # Use the works_fun function to see if the GAM can be resolved for the species in the df
  works_fun_out <- lapply(c(col_names), works_fun, data=df.pred,cors="corCAR1")
  works_fun_out <- data.frame(unlist(works_fun_out, use.names=FALSE))
  works_fun_out$cols <- col_names
  colnames(works_fun_out) <- c("TF","cols")

  # Separate species for which GAM worked (keep) from species for which GAM did not work (try_again)
  keep_list <- subset(works_fun_out, TF=="FALSE")
  keep_list_names <- keep_list$cols
  keep <- subset(df.pred, select=c(keep_list_names,"MonthOfYear","MonthCount"))

  remove_list <- subset(works_fun_out, TF=="TRUE")
  remove_list_names <- remove_list$cols
  try_again<- subset(df.pred, select=c(remove_list_names))
  try_again$MonthOfYear <- df.pred$MonthOfYear
  try_again$MonthCount <- df.pred$MonthCount

  # We have 2 predictor columns in the df (MonthOfYear, MonthCount)
  n <- ncol(keep)-2
  dates <- data.frame(rownames(keep))
  names <- data.frame(colnames(keep))
  col_names <- colnames(keep[1:n])

  # AC cutoff values
  cut <- 1.96/sqrt(nrow(keep)-1)

  # Use ac_remove_fun to see if significant AC is left over in the GAM residuals
  ac_remove_fun_out <- lapply(c(col_names), ac_remove_fun, data=keep,cors="corCAR1",cutoff=cut)
  ac_remove_fun_out <- data.frame(unlist(ac_remove_fun_out, use.names=FALSE))
  ac_remove_fun_out$cols <- col_names
  colnames(ac_remove_fun_out) <- c("AC","cols")

  # If there are species with sinificant AC, remove them and add them to the try_again df
  if(sum(ac_remove_fun_out$AC==1)>=1){
    remove <- subset(ac_remove_fun_out, AC==1)
    remove <- subset(keep,select=c(remove$cols))
    try_again <- cbind(remove,try_again)
    remove <- subset(ac_remove_fun_out, AC==1)
    nums <- which(colnames(keep) %in% remove$cols)
    keep <- subset(keep, select=-c(nums))}

  # We have 2 predictor columns in the df (MonthOfYear, MonthCount)
  n <- ncol(try_again)-2
  dates <- data.frame(rownames(try_again))
  names <- data.frame(colnames(try_again))

  # Try corAR1 cor structure
  if(n>0){
    col_names <- colnames(try_again[1:n])
    workz_fun_out <- lapply(c(col_names), works_fun, data=try_again, cors="corAR1")
    works_fun_out <- data.frame(unlist(works_fun_out, use.names=FALSE))
    works_fun_out$cols <- col_names
    colnames(works_fun_out) <- c("TF","cols")
    not_work <- sum(works_fun_out$TF=="TRUE")
    dates <- data.frame(rownames(try_again))
    names <- data.frame(colnames(try_again))

    # See if corAR1 removes significant AC
    cut <- 1.96/sqrt(nrow(try_again)-1)

    if(not_work==0){
      col_names <- colnames(try_again[1:n])
      ac_remove_fun_out <- lapply(c(col_names), ac_remove_fun, data=try_again,cors="corAR1",cutoff=cut)
      ac_remove_fun_out <- data.frame(unlist(ac_remove_fun_out, use.names=FALSE))
      ac_remove_fun_out$cols <- col_names
      colnames(ac_remove_fun_out) <- c("AC","cols")
      AR1_AC <- sum(ac_remove_fun_out$AC==1)}
    if (not_work>0){AR1_AC <- NA}}

  # Try corCompSymm cor structure
  # We have 2 predictor columns in the df (MonthOfYear, MonthCount)
  n <- ncol(try_again)-2
  if(n>0){
    col_names <- colnames(try_again[1:n])
    works_fun_out <- lapply(c(col_names), works_fun, data=try_again,cors="corCompSymm")
    works_fun_out <- data.frame(unlist(works_fun_out, use.names=FALSE))
    works_fun_out$cols <- col_names
    colnames(works_fun_out) <- c("TF","cols")
    not_work <- sum(works_fun_out$TF=="TRUE")
    dates <- data.frame(rownames(try_again))
    names <- data.frame(colnames(try_again))

    # AC cutoff values
    cut <- 1.96/sqrt(nrow(try_again)-1)

    if(not_work==0){
      col_names <- colnames(try_again[1:n])
      ac_remove_fun_out <- lapply(c(col_names), ac_remove_fun, data=try_again,cors="corCompSymm",cutoff=cut)
      ac_remove_fun_out <- data.frame(unlist(ac_remove_fun_out, use.names=FALSE))
      ac_remove_fun_out$cols <- col_names
      colnames(ac_remove_fun_out) <- c("AC","cols")
      CompSymm_AC <- sum(ac_remove_fun_out$AC==1)}
    if (not_work>0){CompSymm_AC <- NA}}

  # Try corExp cor structure
  # We have 2 predictor columns in the df (MonthOfYear, MonthCount)
  n <- ncol(try_again)-2

  if(n>0) {
    col_names <- colnames(try_again[1:n])
    workz_fun_out <- lapply(c(col_names), works_fun, data=try_again,cors="corExp")
    works_fun_out <- data.frame(unlist(works_fun_out, use.names=FALSE))
    works_fun_out$cols <- col_names
    colnames(works_fun_out) <- c("TF","colz")
    not_work <- sum(works_fun_out$TF=="TRUE")
    dates <- data.frame(rownames(try_again))
    names <- data.frame(colnames(try_again))

    # AC cutoff values
    cut <- 1.96/sqrt(nrow(try_again)-1)

    if(not_work==0){
      col_names <- colnames(try_again[1:n])
      ac_remove_fun_out <- lapply(c(col_names), ac_remove_fun, data=try_again,cors="corExp",cutoff=cut)
      ac_remove_fun_out <- data.frame(unlist(ac_remove_fun_out, use.names=FALSE))
      ac_remove_fun_out$cols <- col_names
      colnames(ac_remove_fun_out) <- c("AC","colz")
      Exp_AC <- sum(ac_remove_fun_out$AC==1)}
    if (not_work>0){Exp_AC <- NA}}

  # Try corGaus cor structure
  # We have 2 predictor columns in the df (MonthOfYear, MonthCount)
  n <- ncol(try_again)-2

  if(n>0){
    col_names <- colnames(try_again[1:n])
    works_fun_out <- lapply(c(col_names), works_fun, data=try_again,cors="corGaus")
    works_fun_out <- data.frame(unlist(works_fun_out, use.names=FALSE))
    works_fun_out$cols <- col_names
    colnames(works_fun_out) <- c("TF","colz")
    not_work <- sum(works_fun_out$TF=="TRUE")
    dates <- data.frame(rownames(try_again))
    names <- data.frame(colnames(try_again))

    # AC cutoff values
    cut <- 1.96/sqrt(nrow(try_again)-1)

    if(not_work==0){
      col_names <- colnames(try_again[1:n])
      ac_remove_fun_out <- lapply(c(col_names), ac_remove_fun, data=try_again,cors="corGaus",cutoff=cut)
      ac_remove_fun_out <- data.frame(unlist(ac_remove_fun_out, use.names=FALSE))
      ac_remove_fun_out$cols <- col_names
      colnames(ac_remove_fun_out) <- c("AC","cols")
      Gaus_AC <- sum(ac_remove_fun_out$AC==1)}

    if (not_work>0){Gaus_AC <- NA}}

  # Determine which GAM equation will be used for the species that the GAM could not be resolved for and for which signficant AC was left in the corCAR1 GAM residuals
  if(n>0){
    x <- data.frame(a=c(AR1_AC,CompSymm_AC, Exp_AC,Gaus_AC))
    x$b <- c("corAR1","corCompSymm","corExp","corGaus")
    x <- subset(x, !is.na(a))
    x <- dplyr::arrange(x,a)
    x <- as.data.frame(x)
    method <- x[1,2]}

  # Get residuals from keep df for network using the corCAR1 model
  # We have 2 predictor columns in the df (MonthOfYear, MonthCount)
  n <- ncol(keep)-2
  col_names <- colnames(keep[1:n])
  dates <- data.frame(rownames(keep))
  names <- data.frame(colnames(keep))

  output.asvs.fin <- lapply(c(col_names), get_resid, data=keep,cors="corCAR1")
  output.asvs.fin <- data.frame(matrix(unlist(output.asvs.fin), nrow=length(output.asvs.fin), byrow=TRUE))
  output.asvs.fin <- as.data.frame(t(output.asvs.fin))
  colnames(output.asvs.fin)<-names$colnames.keep.[1:n]
  rownames(output.asvs.fin)<-dates$rownames.keep.

  # Get residuals from try_again df for network using the best GAM cor class determined above
  # We have 2 predictor columns in the df (MonthOfYear, MonthCount)
  n <- ncol(try_again)-2
  col_names <- colnames(try_again[1:n])
  if(n>0){
    dates <- data.frame(rownames(try_again))
    names <- data.frame(colnames(try_again))
    output.asvs.fin2 <- lapply(c(col_names), get_resid, data=try_again,cors=method)
    output.asvs.fin2 <- data.frame(matrix(unlist(output.asvs.fin2), nrow=length(output.asvs.fin2), byrow=TRUE))
    output.asvs.fin2 <- as.data.frame(t(output.asvs.fin2))
    colnames(output.asvs.fin2)<-col_names
    rownames(output.asvs.fin2)<-dates$rownames.try_again.

    output.asvs.gam <- cbind(output.asvs.fin,output.asvs.fin2)}
  if(n==0){output.asvs.gam <- output.asvs.fin}

  # Run the networks on the GAM-transformed  df
  # GAM-graphical lasso
  output.asvs.gam <- huge::huge.npn(output.asvs.gam)
  output.asvs.gam <- as.matrix(output.asvs.gam)

  if (method=="glasso"){
    lams  <- pulsar::getLamPath(pulsar::getMaxCov(output.asvs.gam), .01, len=30)
    hugeargs <- list(lambda=lams, verbose=FALSE)
    out.p <- pulsar::batch.pulsar(output.asvs.gam, fargs=hugeargs,rep.num=50, criterion = "stars")
    opt <- out.p$stars
    n <- opt$opt.index
    # Get output adjacency matrix from graphical lasso model
    fit <- refit(out.p)
    fit2 <- fit$refit
    fit.fin <- fit2$stars
    fit.fin <- as.matrix(fit.fin)
    fit.fin <- as.data.frame(fit.fin)
    colnames(fit.fin) <- colnames(output.asvs.gam)
    rownames(fit.fin)<- colnames(output.asvs.gam)
    fit.fin <- as.matrix(fit.fin)
    return(fit.fin)}

  # GAM-SCC
  if (method=="scc"){
    output.asvs.cor <- psych::corr.test(output.asvs.gam, method="spearman",adjust="bonferroni")
    p.val <- as.matrix(output.asvs.cor$p)
    diag(p.val) <- 1
    p.val[p.val==1] <- 0.99
    p.val[p.val<pvalue] <- 1
    p.val[p.val !=1] <- 0
    fit.fin <- as.matrix(p.val)
    return(fit.fin)}

  # GAM-PCC
  if (method=="pcc"){
    output.asvs.cor <- psych::corr.test(output.asvs.gam, method="pearson",adjust="bonferroni")
    p.val <- as.matrix(output.asvs.cor$p)
    diag(p.val) <- 1
    p.val[p.val==1] <- 0.99
    p.val[p.val<pvalue] <- 1
    p.val[p.val !=1] <- 0
    fit.fin <- as.matrix(p.val)}
  return(fit.fin)}








