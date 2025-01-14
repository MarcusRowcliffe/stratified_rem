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
