#' @name summary.batss
#' @title Summary function for 'BATSS' outputs
#' @description Summary method function for objects of class 'batss'.
#' @param object An object of class 'batss' (i.e., output of the function [batss.glm]).
#' @param extended A logical indicating if a standard (extended = FALSE, default) or extended output (extended = TRUE) should be returned. Default to `NULL` in which case the input of the argument `extended` chosen when generating `object` with [batss.glm()] is used.
#' @param ... Additional arguments affecting the summary produced.
#' @returns Object of class 'summary.batss'.
#' @returns The function [summary.batss] returns an S3 list of class 'summary.batss' with available print functions. The list elements are
#' \itemize{
#'   \item beta - A data frame providing information related to the beta parameter vector, such as parameter names and values, for example.
#'   \item look - A data frame providing information related to looks, like sample size of a given interim (m) and cumulative sample size at a given interim (n), for example.
#'   \item par - A list providing different information, like the used seeds (seed) and the groups (group), for example.
#'   \item H1 - A list providing trial aggregated results under the alternative, like the probability of efficacy, futility, or both, per arm or globally (`object$H1$target`), the probability of stopping early for efficacy (`object$H1$efficacy`) and futility (`object$H1$futility`), the sample size expectation, standard deviation, and quantiles 0.1, 0.5 and 0.9, per group and overall (`object$H1$summary.sample.sizes`), the probabilities associated to each combination of efficacy and futility per group (scenario).
#'   \item H0 - A list providing trial aggregated results under the global null hypothesis (same structure as H1).
#'   \item call - The matched call.
#'   \item type - The type of 'BATSS' analysis (only 'glm' is currently available).
#' }
#' @seealso [batss.glm()], the function generating S3 objects of class 'batss'. 
#' @export
summary.batss = function(object, extended=NULL, ...){

  res          <- list()
  res$call     <- object$call
  res$par      <- object$par
  res$type     <- object$type
  res$beta     <- object$beta
  res$sample   <- object$look
  colnames(res$sample)[1] = ""
  if(is.null(extended)){
    res$extended   <- ifelse(is.null(object$call$extended),FALSE,object$call$extended>0)
  }else{
    res$extended   <- as.numeric(extended) 
  }

  if (object$par$H0) {
    
    # sample sizes
    res$H0$sample.sizes <- object$H0$sample[,1:sum(object$beta$target)]    
    # target parameters
    res$H0$target <- object$H0$target
    temp = rbind(object$H0$target$par,object$H0$target$global)
    if(all(temp$both==0)){temp = temp[,colnames(temp)!="both"]}       
    res$H0$target <- temp
    # efficacy    
    temp = rbind(object$H0$efficacy$par,object$H0$efficacy$global)
    temp[,1][-(1:nrow(object$H0$efficacy$par))] = ""
    res$H0$efficacy <- temp
    # futility
    temp = rbind(object$H0$futility$par,object$H0$futility$global)
    temp[,1][-(1:nrow(object$H0$futility$par))] = ""
    res$H0$futility <- temp
    
    if(res$extended>0){    
      # ss
      temp = object$H0$sample[,1:nrow(object$par$group)]
      temp = cbind(temp,Total=apply(temp,1,sum))
      temp = data.frame(pos=1:ncol(temp),
                        id=colnames(temp),
                        'ESS'=round(apply(temp,2,mean),2),
                        'St.Dev'=round(sqrt(apply(temp,2,var)),2),
                        'q10'=round(apply(temp,2,quantile,probs=0.1),2),
                        'q50'=round(apply(temp,2,quantile,probs=0.5),2),
                        'q90'=round(apply(temp,2,quantile,probs=0.9),2))
      colnames(temp)[c(1,5:7)] = c("",paste0("q(",c(0.1,0.5,0.9),")"))
      temp[,1][nrow(temp)] = ""
      res$H0$summary.sample.sizes <- temp
      # scenarios    
      res$H0$scenario <- object$H0$scenario
    }    
  }else{
    res$H0 <- NULL 
  } 
  
  if(object$par$H1){

    # sample sizes
    res$H1$sample.sizes <- object$H1$sample[,1:sum(object$beta$target)]
    # target parameters
    res$H1$target <- object$H1$target
    temp = rbind(object$H1$target$par,object$H1$target$global)
    if(all(temp$both==0)){temp = temp[,colnames(temp)!="both"]}       
    res$H1$target <- temp
    # efficacy    
    temp = rbind(object$H1$efficacy$par,object$H1$efficacy$global)
    temp[,1][-(1:nrow(object$H1$efficacy$par))] = ""
    res$H1$efficacy <- temp
    # futility
    temp = rbind(object$H1$futility$par,object$H1$futility$global)
    temp[,1][-(1:nrow(object$H1$futility$par))] = ""    
    res$H1$futility <- temp

    if(res$extended>0){    
      # ss
      temp = object$H1$sample[,1:nrow(object$par$group)]
      temp = cbind(temp,Total=apply(temp,1,sum))
      temp = data.frame(pos=1:ncol(temp),
                        id=colnames(temp),
                        'ESS'=round(apply(temp,2,mean),2),
                        'St.Dev'=round(sqrt(apply(temp,2,var)),2),
                        'q10'=round(apply(temp,2,quantile,probs=0.1),2),
                        'q50'=round(apply(temp,2,quantile,probs=0.5),2),
                        'q90'=round(apply(temp,2,quantile,probs=0.9),2))
      colnames(temp)[c(1,5:7)] = c("",paste0("q(",c(0.1,0.5,0.9),")"))
      temp[,1][nrow(temp)] = ""
      
      res$H1$summary.sample.sizes <- temp
      # scenarios    
      res$H1$scenario <- object$H1$scenario
    } 
  } else {
    res$H1 <- NULL 
  } 
  class(res) = "summary.batss"
  return(res)
}
  
#' @name print.summary.batss
#' @title Print function for objects of class 'summary.batss' 
#' @description Print function for objects of class 'summary.batss' 
#' @param x An object of class 'summary.batss' (i.e., output of the function [summary] used on an output of the function [batss.glm]).
#' @param ... Additional arguments affecting the summary produced.
#' @returns Prints a summary for objects of class 'batss'.
#' @seealso [batss.glm()], the function generating S3 objects of class 'batss'. 
#' @export 
print.summary.batss = function(x, ...){
  object = x
  cat("\n")
  if(!is.null(object$par$RAR)){
    cli_h1("Bayesian Adaptive Design with Laplace Approx.")
  }else{
    cli_h1("MAMS with Laplace Approx.")
  }
  cat("  (",length(object$par$seed)," Monte Carlo samples)\n",sep="")
  cat("\n")
  cli_h3("Variables:")
  for(i in 2:length(object$call$var)){
    cat("  *",names(object$call$var)[i],":",as.character(object$call$var)[i],"\n")
  }
  if(!is.null(object$par$RAR)){
    cat("\n")
    cli_h3("Group randomisation:")
    cat("  *",object$call$RAR,"\n")
  }
  if(!(is.null(object$par$eff.arm) && is.null(object$par$fut.arm))){
    cat("\n")
    cli_h3("Decision rules:")
    if(!is.null(object$par$eff.arm)) cat("  * Efficacy: ",format(object$call$eff.arm),"\n")
    if(!is.null(object$par$fut.arm)) cat("  * Futility: ",format(object$call$fut.arm),"\n")
  }
  cat("\n")
  cli_h3("Model: ")
  cat("  *",format(object$call$model),paste("(with",ifelse(is.null(object$call$link),"identity",object$call$link), "link)\n"))
  cat("\n")
  cli_h3("Fixed effect parameters:\n")
  objectw = object$beta
  colnames(objectw)[1] = ""
  print(objectw,row.names=FALSE)
  cat("\n")
  cli_h3("Sample size per interim analyis:\n")
  objectw = object$sample
  colnames(objectw)[1] = ""
  print(objectw,row.names=FALSE)

  # H0
  if(object$par$H0){
    cat("\n\n")
    #cat(paste0(rep("-",floor(options()$width/2)),collapse=""))
    cli_h2("\n H0: Under the null hypothesis\n")
    #cat(paste0(rep("-",floor(options()$width/2)),collapse=""))
    #
    cat("\n")
    cli_h3("Target parameters:\n")
    print(object$H0$target,row.names=FALSE)
    #
    if(object$extended>0){
      cat("\n")
      cli_h3("Efficacy:\n")
      print(object$H0$efficacy,row.names=FALSE)
      #
      cat("\n")
      cli_h3("Futility:\n")
      print(object$H0$futility,row.names=FALSE)
      #
      cat("\n")
      cli_h3("Sample size per group:\n")
      print(object$H0$summary.sample.sizes,row.names=FALSE)
      #cat("\n")
      #
      cat("\n")
      cli_h3("Scenarios:\n")
      print(object$H0$scenario,row.names=FALSE)
      cat(" where 0 = no stop, 1 = efficacy stop, 2 = futility stop\n")
      if(any(object$H0$scenario[,object$H1$target$id[-nrow(object$H0$target)+c(0:1)]]==3)){
        cat(",\n       3 = simultaneous efficacy and futility stops")
      }else{cat("\n")}
    }
  }
  # H1
  if(object$par$H1){
    cat("\n\n")
    #cat(paste0(rep("-",floor(options()$width/2)),collapse=""))
    cli_h2("\n H1: Under the alternative hypothesis\n")#cat("\n H1: Under the alternative hypothesis\n")
    #cat(paste0(rep("-",floor(options()$width/2)),collapse=""))
    #
    cat("\n")
    cli_h3("Target parameters:\n")
    print(object$H1$target,row.names=FALSE)
    #
    if(object$extended>0){
      cat("\n")
      cli_h3("Efficacy:\n")
      print(object$H1$efficacy,row.names=FALSE)
      #
      cat("\n")
      cli_h3("Futility:\n")
      print(object$H1$futility,row.names=FALSE)
      #
      cat("\n")
      cli_h3("Sample size per group:\n")
      print(object$H1$summary.sample.sizes,row.names=FALSE)
      #cat("\n")
      cat("\n")
      cli_h3("Scenarios:\n")
      print(object$H1$scenario,row.names=FALSE)
      cat(" where 0 = no stop, 1 = efficacy stop, 2 = futility stop")
      if(any(object$H1$scenario[,object$H1$target$id[-nrow(object$H1$target)+c(0:1)]]==3)){
        cat(",\n       3 = simultaneous efficacy and futility stops")
      }else{cat("\n")}
    }
  }
  cli_h1("")#cat(paste0(rep("-",options()$width),collapse=""))
  cat("\n")
}

