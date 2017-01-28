#genv <- .GlobalEnv

# Name: getOS
# Tasks: Returns the guessed operating system from a variety of functions in R.
#        Uses Sys.Info() and .Platform finally returns "osx", "windows" or "linux"
# Example: os <- getOS()
getOS <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

# Name: loadPkgs
# Task: A function that takes a vector of names of packages, installs the ones
#       that are missing, loads all the packages and returns a list of the 
#       missing packages that had to be installed
# Example: Create the vector of packages
#       packages <- c("ggplot2", "dplyr", "Hmisc", "lme4", "arm", "lattice")
#       Call the function
#       loadPkgs(packages)
loadPkgs <- function(packages){
  missingPkgs <- setdiff(packages, rownames(installed.packages()))
  if (length(missingPkgs) > 0) {
    install.packages(missingPkgs)  
  }
  lapply(packages, require, character.only = TRUE)
  return(missingPkgs)
}

# Name: corstar
# Task: Uses corr.test from psych package to calculate correlations and
#       statistical significance (at p< 0.001, p< 0.01 and p< 0.05) and returns
#       a matrix of correlations, p-values and stars corresponding to 
#       significance
corstar <- function(x, y = NULL, use = "pairwise", method = "pearson", 
                    round = 3, row.labels, col.labels, ...) {
  require(psych)
  
  ct <- corr.test(x, y, use, method)    # calculate correlation
  r <- ct$r                             # get correlation coefs
  p <- ct$p                             # get p-values
  
  # generate significance stars
  stars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "*  ", "   ")))
  
  m <- matrix(NA, nrow = nrow(r) * 2, ncol = ncol(r) + 1) # create empty matrix
  
  rlab <- if(missing(row.labels)) rownames(r) else row.labels # add row labels
  clab <- if(missing(col.labels)) {
    if(is.null(colnames(r)))
      deparse(substitute(y))
    else
      colnames(r)
  } else {
    col.labels # add column labels
  }
  
  rows <- 1:nrow(m)                     # row indices
  cols <- 2:ncol(m)                     # column indices
  
  odd <- rows %% 2 == 1                 # odd rows
  even <- rows %% 2 == 0                # even rows
  m[odd, 1] <- rlab                     # add variable names
  m[even, 1] <- rep("", sum(even))      # add blank
  
  # add r coefs
  m[odd, cols] <- paste(format(round(r, round), nsmall = round, ...), stars, sep = "")
  
  # add p-values
  m[even, cols] <- paste("(", format(round(p, round), nsmall = round, ...), ")", sep = "")
  
  colnames(m) <- c(" ", clab)           # add colnames
  m                                     # return matrix
}

# Table of correlations with significance stars
corstars <- function(x){
  require(psych)
  x <- as.matrix(x)
  R <- corr.test(x)$r
  p <- corr.test(x)$p
  
  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "*  ", "   ")))
  
  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew)
  
  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew)
}

# Name: fixNames
# Task: In imported data the names of columns is often the question that was
#       asked. This can cause a lot of issues when printing reports etc. This
#       function takes a dataframe, and if any name is larger than maxlength
#       (default 80 chars), all names are replace with a pseudo name
fixNames <- function(x, maxlength = 80){
  require(stringr)
  numLongNames <- sum(str_length(names(x)) > maxlength)
  if(numLongNames > 0){
    c(paste0("x",seq(1,length(x))))
  } else {
    names(x)
  }  
}

# Returns a named vector of number of NAs in each column of a matrix or data
# frame
countNA <- function(x){
  sum(is.na(x))
}
checkNA <- function(x){
  NADat <- apply(x,2,countNA)
  return(NADat)
}

# Check the autocorrelation of variables to decide which variables might be
# sequences or ids and therefore not data collected from respondents
autoCorOneLag <- function(x){
  AcFV <- acf(x, lag.max = 1, type = "correlation", plot = FALSE, na.action = na.pass)
  AcFV$acf[2]
}
checkACFVar <- function(x, maxacf = 0.80){
  bVecACF <- (apply(x,2,autoCorOneLag) > maxacf)
  a <- names(x[bVecACF])
  return(a)
}

# Checking the SD of responses given by a single respondent. If there is very
# little variance, the respondent probably filled the data carelessly.
# Returns a vector of the respondents that have sd < minsd (default 0.35)
# This might need to be modified depending upon the scale used. For a 5 point
# scale (1 - 5), this value is suitable
RespSD <- function(x,minsd = 0.35){
  sdResp <- apply(x,1,sd,na.rm = TRUE)
  if(sum((sdResp < minsd)) != 0) {
    respsWithLowSD <- order(sdResp)[1:sum(sdResp < minsd)]
  } else {
    respsWithLowSD <- "None"
  }
  return(list("RespsWithLowSD" = respsWithLowSD, "Threshold" = minsd))
}

# Survey software sometimes have columns for start time and end time. These
# could be used to confirm that respondents have filled after reading the
# questions.
# Returns a list of vectors:
# $StartTime   : Vector of start times formatted as POSIXlt
# $EndTime     : Vector of end times formatted as POSIXlt
# $TimeTaken   : Amount of time taken in minutes
# $LessThanMin : Respondents that took less than mintime (default = 3)
# $MinTime     : Defined minimum time
GetTimeTaken <- function(x,y,z, format = "%Y-%m-%d %H:%M:%S", mintime = 3){
  require(pander)
  TimeStarted <- strptime(x[,y], format = format)
  TimeEnded   <- strptime(x[,z], format = format)
  TimeTaken   <- as.numeric(TimeEnded - TimeStarted, units = "mins")
  numLessThanMin <- sum(TimeTaken < mintime)
  if(numLessThanMin != 0) {
    LessThanMin <- order(TimeTaken)[1:numLessThanMin]
  } else {
    LessThanMin <- "None"
  }
  catlimit("Average time taken:", mean(TimeTaken))
  catlimit("Minimum time taken:", min(TimeTaken))
  if(numLessThanMin != 0) {
    catlimit(paste0("Some respondents have taken less than ", mintime, " minutes to complete the survey. This could be an anomaly or random responses. You could consider removing these respondents."))
    catlimit("Time taken (in minutes) by these respondents was: ", pandoc.list.return(TimeTaken[LessThanMin]))
  }

  return(list("StartTime" = TimeStarted, "EndTime" = TimeEnded,
              "TimeTaken" = TimeTaken, "LessThanMin" = LessThanMin,
              "MinTime" = mintime))
}

catlimit <- function(..., limit = TRUE){
  cat(..., fill = limit)
}

# Return dataframe with only the numeric data
# Useful for using scoreItems from psych package
GetNumericData <- function(x){
    return(x[sapply(x,is.numeric)])
}

# Write comments
numAnalysisDF <- function(x){
  require(pander)
  require(psych)

  numData <- GetNumericData(x)
  names(numData) <- fixNames(numData)
  colsWithNA <- checkNA(numData)
  cat("\nMISSING VALUE TEST")
  cat("\n------------------")
  if(sum(colsWithNA) == 0){
    catlimit("\nNumeric Data do not contain any missing values")
  } else {
    catlimit("\nNumeric Data contains missing values. Missing values would be ignored for subsequent analysis. Number of missing values for each variable is given below:")
    pander(colsWithNA, style = "rmarkdown")
  }
  cat("\nAUTOCORRELATION TEST")
  cat("\n--------------------")
  drops <- checkACFVar(numData)
  catlimit("\nThe following variable(s) have an autocorrelation at a lag of 1 of more than 0.8. These might be sequence or ID type of variables, they will be excluded from further analysis.")
  catlimit(paste0("Variables excluded: ", pandoc.list.return(drops)))
  numData <- numData[!names(numData) %in% drops]

  cat("\nVARIANCE TEST")
  cat("\n-------------")
  lowVarResp <- RespSD(numData)
  catlimit(paste0("\nThe following respondents have very little variance (< ", as.character(lowVarResp$Threshold), ")", " in their responses, without accounting for missing values. This could be a case of disinterest, consider removing these respondents."))
  catlimit(paste0("Respondents with low variance: ", pandoc.list.return(lowVarResp$RespsWithLowSD)))

  cat("\nDATA DESCRIPTION")
  cat("\n----------------")
  DescData <- describe(numData)
  catlimit(paste0("\nRange of skewness of variables: ", range(DescData$skew)))
  catlimit(paste0("\nRange of kurtosis of variables: ", range(DescData$kurtosis)))
  pander(DescData, style = "rmarkdown")
  #rm(list = c("drops","numData","DescData"))
  return(lowVarResp$RespsWithLowSD)
}
