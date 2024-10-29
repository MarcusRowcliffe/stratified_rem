library(tidyverse)
library(activity)
library(circular)
library(lubridate)
library(camtrapDensity)
library(scales)
library(ggplot2)

# Von Mises mixture probability density

# INPUT
# x: vector of values at which to evaluate
# prm: named vector of distribution parameters;
#      names ki, mi, pi respectively k, m, p for the ith component,
#      p with 1 fewer than k/m, last is inferred since sum(pi) = 1

# OUTPUT
# A vector of probability densities
dvmsmix <- function(x, prm){
  k <- exp(prm[grep("k",names(prm))])
  m <- prm[grep("m",names(prm))]
  p <- prm[grep("p",names(prm))]
  if(length(k)==1) p  <- 1 else p <- c(1-sum(p), p)
  res <- rep(0, length(x))
  for(i in 1:length(k))
    res <- res + p[i]*dvonmises(circular(x),circular(m[i]),k[i])
  res
}

# Generate random radian time observations from an empirical distribution function

# INPUT
# n: number of observations to generate
# fit: dataframe describing the empirical distribution, with columns:
#      x: a regular sequence of times from 0 to 2*pi
#      y: the (relative) probabilities of observing times

# OUTPUT
# A vector of radian times
rand_time <- function (n, fit){
  if (sum(c("x", "y") %in% names(fit)) != 2) 
    stop("fit must be a dataframe with (at least) columns named x and y")
  if (diff(range(diff(fit$x))) > 1e-04) 
    stop("x doesn't seem to be a regular sequence")
  cdf <- c(0, cumsum(head(fit$y, -1)) / sum(head(fit$y, -1)))
  rn <- stats::runif(n)
  res <- stats::approx(cdf, fit$x, rn)$y
  res <- res - mean(fit$x[1:2])
  res <- ifelse(res<0, res+2*pi, res)
  res
}

# Make a time x stratum matrix from a table of parameter values
#
# INPUT
#   table: a dataframe of parameter values with (at least) fields:
#     est: parameter estimates
#     timeID: (otional) integer time period IDs
#     stratumID: (optional) stratum IDs
#   str_data: dataframe of stratum data with (at least) a stratumID field
#   time_thresholds: vector of radian times demarcating time periods
#   nt: number of time points to generate for matrices
# OUTPUT
#   A matrix with dimensions nt rows and nrow(str_data) columns.
#
# DETAILS
# parameters$timeID values must be in 1:length(time_thresholds)
# parameters$stratumID values must match str_data$stratumID
# Each row of parameters must represent a unique combination of timeID 
# and stratumID (if both are present).
# time_thresholds must be strictly increasing, defining period 1 as
# tt[1]< <=tt[2] and so on.
pmat_from_table <- function(table, str_data, time_thresholds=NULL, nt=513){
  if(!"est" %in% names(table)) stop("table must contain est column")
  if("stratumID" %in% names(table) &
     !all(str_data$stratumID %in% unique(table$stratumID)))
    stop("Not all str_data$stratumID can be found in table$stratumID")
  if("timeID" %in% names(table)){
    if(is.null(time_thresholds)) 
      stop("time_thresholds must be provided if table contains timeID")
    if(min(time_thresholds)<0 | 
       max(time_thresholds)>2*pi | 
       any(diff(time_thresholds) <= 0))
      stop("time_thresholds must be an increasing sequence of radian times between 0 and 2*pi")
    tID <- unique(table$timeID)
    ntt <- length(time_thresholds)
    if(length(tID) != ntt | !all(range(tID)==c(1,ntt)))
      stop("table$timeID must contain integers from 1 to length(time_threshold)")
  }
  
  ni <- length(time_thresholds)
  ns <- nrow(str_data)
  sq <- seq(0, 2*pi, len=nt)
  
  if("timeID" %in% names(table)){
    i <- expand.grid(sq, time_thresholds) %>%
      apply(1, function(x) sign(diff(x))) %>%
      matrix(ncol=ni) %>%
      apply(1, sum) %>%
      (function(x) (ni-x)/2) %>%
      floor() %>%
      (function(x){
        x[x==0] <- ni
        x
      })
    if("stratumID" %in% names(table)){
      i <- i %>%
        expand.grid(str_data$stratumID) %>%
        apply(1, paste, collapse="") %>%
        match(paste0(table$timeID, table$stratumID))
    } else
      i <- match(i, table$timeID)
  } else{
    if("stratumID" %in% names(table)){
      i <- match(str_data$stratumID, table$stratumID) %>%
        rep(each=nt)
    } else
      i <- 1
  }
  res <- matrix(table$est[i], nrow=nt, ncol=ns)
  colnames(res) <- str_data$stratumID
  res
}

# Make speed/radius/angle parameter array from a function

# INPUT
# str_data: dataframe with a row per stratum and stratumID column
# f: a function with arguments time and stratum, used to generate time- and/or stratum-specific parameter values

# OUTPUT
# Returns a 513 times by strata matrix of parameter values
pmat_from_func <- function(func, str_data, nt=513){
  if(!any(class(func) %in% "function"))
    stop("func must be a function")
  args <- formalArgs(func)
  if(!all(args %in% c("time", "stratum")))
    stop("Function func cannot have arguments other than time and stratum")
  if(!"stratumID" %in% names(str_data)) 
    stop("str_data must contain a stratumID column")  
  
  narg <- length(args)
  nstr <- nrow(str_data)
  sq <- seq(0, 2*pi, len=nt)
  res <- if(narg==0) func() else
    if(narg == 1){
      if(args=="time") func(sq) else
        if(args=="stratum") rep(func(str_data$stratumID), each=nt)
    } else
      func(time = rep(sq, nstr),
           stratum = rep(str_data$stratumID, each=nt))
  res <- matrix(res, nrow=nt, ncol=nstr)
  colnames(res) <- str_data$stratumID
  res
}

# Generate randomised dataset

# INPUT
# dep_dat: dataframe of deployment data, columns deploymentID, stratumID, effort
# str_dat: dataframe of stratum data, columns stratumID, area
# N: scalar population size
# a: matrix of activity levels, eval times (wrapped) by strata
# p: matrix of proportion in stratum, eval times (wrapped) by strata
# sfunc/rfunc/afunc: NULL or a function with arguments time and stratum, passed to make_parmat
# size: negative binomial size parameter for random observations per deployment generation
# NB strata in above vectors/matrices must occur in a consistent order (there is no name matching)

# OUTPUT
# List holding input dep_dat and str_dat plus obs_dat, a dataframe of
# observations with columns deploymentID and time (radian time of day of 
# each observation)
data_gen <- function(dep_dat, str_dat, 
                     N, a, p, 
                     sfunc=NULL, rfunc=NULL, afunc=NULL,
                     size=10){
  nd <- nrow(dep_dat)
  ns <- nrow(str_dat) # number of strata
  
  spd <- pmat_from_func(sfunc, str_dat, 513)
  rad <- pmat_from_func(rfunc, str_dat, 513)
  ang <- pmat_from_func(afunc, str_dat, 513)
  q <- spd * rad * (2+ang)
  if(length(q) > 1) q <- t(q)  # parameter weight, stratum x time
  D <- t(p) * N / str_dat$area # density, stratum x time
  E <- dep_dat %>%
    group_by(stratumID) %>%
    summarise(E = sum(effort)) %>%
    pull(E) # effort, stratum
  P <- q * t(a) * D * E / (512 * pi) # expected obs, stratum x time
  ncams <- with(dep_dat, tapply(stratumID, stratumID, length)) # cameras, stratum
  Pdep <- (apply(P[,-1], 1, sum) / ncams)[dep_dat$stratumID] # expected observations, deployment
  Ydep <- rnbinom(nd, mu=Pdep, size=size) # realised observations, deployment
  Y <- tapply(Ydep, dep_dat$stratumID, sum) # realised observations, stratum
  
  sq <- (0:512)*2*pi/512 # generate random observation times
  f <- function(i){
    df <- data.frame(x=sq, y=P[i,])
    rand_time(Y[i], df)
  }
  obs_dat <- data.frame(deploymentID = rep(dep_dat$deploymentID, Ydep),
                        time = unlist(sapply(1:ns, f)))
  list(dep=dep_dat, str=str_dat, obs=obs_dat)
}

# Checks and modifies stratified activity data ready for input to sa function
# INPUT
# obs_data: dataframe of observations with fields deploymentID, time
# dep_data: dataframe of deployments with fields deploymentID, stratumID, effort
# str_data: dataframe of strata with fields stratumID, area
# OUTPUT
# list of dataframes:
#   obs: observations with stratumID field added
#   str: strata with area, effort and observations fields added
make_sa_data <- function(obs_data, dep_data, str_data){
  # Check for required fields and IDs)
  orf <- c("deploymentID", "time")
  drf <- c("deploymentID", "stratumID", "effort")
  srf <- c("stratumID", "area")
  msg <- function(df, ff) paste(df, "must contain columns:", paste(ff, collapse=", "))
  if(!all(orf %in% names(obs_data))) stop(msg("obs_data", orf))
  if(!all(drf %in% names(dep_data))) stop(msg("dep_data", drf))
  if(!all(srf %in% names(str_data))) stop(msg("str_data", srf))
  if(!all(obs_data$deploymentID %in% dep_data$deploymentID)) 
    stop("Cant find all observation deploymentIDs in dep_data")
  if(!all(dep_data$stratumID %in% str_data$stratumID)) 
    stop("Cant find all deployment stratumIDs in str_data")
  
  # add stratumID to obs_dat
  if(!"stratumID" %in% names(obs_data)) 
    obs_data <- dep_data %>%
    select(deploymentID, stratumID) %>%
    right_join(obs_data, by="deploymentID")
  
  # add n observations to str_dat, remove strata with no observations
  str_data <- obs_data %>%
    group_by(stratumID) %>%
    summarise(observations=n()) %>%
    right_join(str_data, by="stratumID") %>%
    mutate(observations = ifelse(is.na(observations), 0, observations))
  # add effort to str_dat
  str_data <- dep_data %>%
    group_by(stratumID) %>%
    summarise(effort=sum(effort)) %>%
    right_join(str_data, by="stratumID")
  
  list(obs=obs_data, str=str_data)
}

# function fits a stratified activity pattern  with or without bootstrapping
# INPUT
# dat: list of dataframes
#     obs with required fields: stratumID, time
#     str with required fields: stratumID, area, effort, observations
# wt: times by strata matrix of weights
sa <- function(dat, wt, adj=1, boot=FALSE){
  nt <- nrow(wt)
  odat <- lapply(dat$str$stratumID, function(s) subset(dat$obs, stratumID==s)$time)
  if(boot) 
    odat <- lapply(odat, function(x) sample(x, length(x), replace=TRUE))
  sq <- seq(0,2*pi,len=nt)
  f <- sapply(odat, function(d) 
    if(length(d)==0) rep(0,nt) else dvmkern(sq, d, adj=adj))
  colnames(f) <- dat$str$stratumID
  fw <- f * wt
  sum_fw <- apply(fw, 1, sum)
  pdf <- sum_fw / (2*pi*mean(head(sum_fw,-1))) #activity pdf over time
  act <- 1 / (2 * pi * max(pdf)) # overall activity level
  popdist <- fw / sum_fw # population distribution between strata
  pa <- popdist * pdf / max(pdf) # population distribution x activity level
  psum <- apply(head(popdist, -1), 2, sum)
  act_stratum <- apply(head(pa, -1), 2, sum) / ifelse(psum==0, 1, psum)
  list(act = act,
       act_stratum = act_stratum,
       pdf = pdf,
       popdist = popdist)
}

# gets summary error list from bootstrap results (boots)
get_errlist <- function(boots){
  lim <- c(0.025, 0.975)
  # overall activity level errors
  boot_act <- unlist(boots["act", ])
  se_act <- sd(boot_act)
  ci_act <- quantile(boot_act, lim)
  
  # stratum-specific activity level errors
  strata <- names(boots[[2]])
  ns <- length(strata)
  nt <- length(boots[[3]])
  reps <- ncol(boots)
  boot_act_stratum <- matrix(unlist(boots["act_stratum", ]), nrow=ns)
  se_act_stratum <- apply(boot_act_stratum, 1, sd)
  ci_act_stratum <- apply(boot_act_stratum, 1, quantile, lim)
  colnames(ci_act_stratum) <- names(se_act_stratum) <- strata
  
  # activity pattern pdf errors
  boot_pdf <- matrix(unlist(boots["pdf", ]), nrow=nt)
  se_pdf <- apply(boot_pdf, 1, sd)
  ci_pdf <- t(apply(boot_pdf, 1, quantile, lim))
  
  # population distribution errors
  boot_popdist <- array(unlist(boots["popdist", ]), dim=c(nt, ns, reps))
  se_popdist <- apply(boot_popdist, 1:2, sd)
  ci <- apply(boot_popdist, 1:2, quantile, lim)
  lcl_popdist <- ci[1,,]
  ucl_popdist <- ci[2,,]
  colnames(lcl_popdist) <- colnames(ucl_popdist) <- colnames(se_popdist) <- strata
  
  list(se = list(act=se_act, act_stratum=se_act_stratum, pdf=se_pdf, 
                 popdist=se_popdist),
       ci = list(act=ci_act, act_stratum=ci_act_stratum, pdf=ci_pdf,
                 popdist=list(lcl=lcl_popdist, ucl=ucl_popdist)))
}

# Fit a stratified activity distribution

# INPUT
# obs_dat: observations data with columns deploymentID and time
# dep_dat: deployment data with columns deploymentID and stratumID and effort
# str_dat: stratum data with columns stratumID and area

# OUTPUT
# A list with elements:
# data: input data, components dep_dat, str_dat, obs_dat
# actmods: stratum-specific activity models
# popdist: time by stratum matrix, proportion of population in each stratum time by time of day
# pdf: probability density function for overall activity pattern
# act: overall activity level
# act_stratum: stratum-specific activity levels
#obs_data <- dat$obs
#str_data <- dat$str
fitact_strat <- function(obs_data, dep_data, str_data, 
                         speed=NULL, radius=NULL, angle=NULL,
                         adj=1, reps=100, nt=513){
  
  check_param <- function(prm){
    strata <- as.character(str_data$stratumID)
    if(is.null(prm)) return(1) else
      if(class(prm)[1]=="matrix" & mode(prm)=="numeric" & !any(is.na(prm))){
        if(nrow(prm)==nt & ncol(prm)==nrow(str_data) & 
           identical(sort(colnames(prm)), sort(strata)))
          return(prm[, strata]) else
            stop("Parameter matrix must have dimensions nt by nrow(str_data), and have column names matching str_data$stratumID")
      } else
        stop("Parameter must be a matrix with no missing values")
  }
  
  dat <- make_sa_data(obs_data, dep_data, str_data)  
  
  speed <- check_param(speed)
  radius <- check_param(radius)
  angle <- check_param(angle)
  q <- speed * radius * (2+angle)
  A <- matrix(rep(dat$str$area, each=nt), nrow=nt)
  E <- rep(dat$str$effort, each=nt)
  Y <- rep(dat$str$observations, each=nt)
  wt <- Y * A / (q * E)
  
  est <- sa(dat, wt, adj)
  if(reps==0) err <- NULL else{
    boots <- pbapply::pbreplicate(reps, sa(dat, wt, adj, TRUE))
    err <- get_errlist(boots)
  }
  res <- list(data=list(obs=dat$obs, dep=dep_data, str=dat$str), 
              est=est, se=err$se, ci=err$ci)
}

time_of_day <- function(x, scale=c("radian", "hour", "proportion")){
  scale <- match.arg(scale)
  m <- switch(scale, 
              radian = 2*pi, 
              hour = 24, 
              proportion = 1)
  
  x <- lubridate::ymd_hms(x) %>% 
    lubridate::local_time() %>% 
    as.numeric()
  m * x  / (24 * 60^2) 
}

# Calculate SE of the mean of several estimates and their SEs
se_from_ses <- function(est, se){
  if(length(est) != length(se))
    stop("All input vectors must have the same length")
  sqrt(mean((mean(est) - est)^2 + se^2))
}

# Stratified REM
# INPUT
#   parameters: dataframe with columns:
#     stratumID: stratum ID
#     param: parameter name
#     estimate: parameter estimate
#     se: parameter standard error
#     unit: character unit (see camtrapDensity::convert_units for allowable values)
#   str_data: stratum data with columns:
#     stratumID: character or numeric stratum ID
#     area: numeric stratum area
#     unit: character area unit
rem_strat <- function(parameters, str_data){
  strfields <- c("stratumID", "area")
  parfields <- c("stratumID", "param", "estimate", "se", "unit")
  params <- c("traprate", "speed", "activity", "radius", "angle")
  ptabval <- unique(with(parameters, table(stratumID, param)))
  if(!all(strfields %in% names(str_data)))
    stop(paste("str_data must contain fields:", paste(depfields, collapse=", ")))
  if(!all(parfields %in% names(parameters)))
    stop(paste("parameters must contain fields:", paste(parfields, collapse=", ")))
  if(length(ptabval) > 1 | ptabval[1]!=1)
    stop("parameters must have exactly one unique combination of each stratumID and param")
  if(!all(str_data$stratumID %in% parameters$stratumID))
    stop("Some stratumID's in str_data are missing from stratumID in parameters")
  if(!all(unique(parameters$param) %in% params))
    
    strata <- str_data$stratumID
  
  # add overall_speed to parameters
  parameters <- parameters %>%
    filter(param %in% c("speed", "activity")) %>%
    group_by(stratumID) %>%
    summarise(est = prod(estimate),
              se = est * sqrt(sum((se/estimate)^2)),
              unit = unit[param=="speed"],
              param = "overall_speed") %>%
    rename(estimate = est) %>%
    bind_rows(parameters)
  
  # stratum specific estimates
  stratum_rem <- lapply(strata, function(s)
    parameters %>%
      subset(stratumID == s) %>%
      column_to_rownames("param") %>%
      camtrapDensity::rem() %>%
      mutate(stratumID = s,
             se = ifelse(is.nan(se), 0, se)) %>%
      rownames_to_column("param")) %>%
    bind_rows() %>%
    mutate(stratumID = as.character(stratumID))
  
  # weighted mean function
  wtd_mean <- function(mn, se, wt){
    c(mean = sum(wt*mn) / sum(wt),
      se = sqrt(sum(wt^2 * se^2) / sum(wt)^2))
  }
  
  prm_unit <- stratum_rem %>%
    filter(stratumID ==strata[1]) %>%
    select(param, unit)
  
  # stratified estimates across all strata
  stratified_rem <- sapply(prm_unit$param, function(prm)
    stratum_rem %>%
      filter(param == prm) %>%
      reframe(x = wtd_mean(estimate, se, str_data$area)) %>%
      .$x) %>%
    t() %>%
    as.data.frame() %>%
    rename(estimate = "V1", se = "V2") %>%
    rownames_to_column("param") %>%
    mutate(stratumID = "stratified",
           unit = prm_unit$unit) %>%
    bind_rows(stratum_rem) %>%
    relocate(stratumID)
  
  return(stratified_rem)
}
