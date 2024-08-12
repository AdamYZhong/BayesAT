#' @title Bayesian inference for survival analysis
#' @description \code{Bayes_test} conduct hypothesis test through Bayesian survival model
#' @param data Matrix. The data contains both survival time and event status.
#' @param alpha Numerical. Gamma distribution alpha parameter.
#' @param beta Numerical. Gamma distribution beta parameter (rate = 1/scale).
#' @param test Categorical. Three types of hypothesis includes "greater", "less", or "two_sided".
#' @param threshold Numerical. The value tested against hypothesis or evidence.
#' @param type Categorical. The types of Bayesian inference include "Posterior" for estimation of parameters or "Predictive" for predicted survival rate.
#' @param pred Numerical. The time point for predicted survival rate, for example, 2 years, or 5 years survival probability.
#' @param diagnosis Logical. If `diagnosis == TRUE`, the Bayes factor is calculated, and the formulation of Bayesian factors is given in details.
#' @returns Bayesian test provide `mean`, `sd`, `CI`, `z_score`, `prob`, and `bf`.
#' @returns `mean` Posterior mean is estimated by calculating the mean of MCMC outputs.
#' @returns `sd` Posterior standard deviation is estimated as the standard deviation of MCMC outputs.
#' @returns `CI`Summary statistics provides the credible intervals and specific quantile.
#' @returns `z_score` Standardized test of statistics is calculated based on MCMC outputs. For example,
#' @returns  \deqn{ \frac{\hat{\lambda} - \lambda_0}{SD( \hat{\lambda} )} \text{ or } \frac{ \hat{S} - S_0}{SD( \hat{S} )},}
#' @returns where \eqn{\hat{\lambda}} is the estimated posterior mean of hazard rate, and \eqn{\hat{S}} is the predicted survival probability. Both \eqn{\lambda_0} and \eqn{S_0} are threshold used for test against hypothesis or evidence.
#' @returns `prob` Posterior probability: \eqn{P(\hat{\lambda} > \lambda_0)} if `test` is "greater", \eqn{P(\hat{\lambda} \le \lambda_0)} if `test` is "less", and \eqn{2 min( P(\hat{\lambda} > \lambda_0),P(\hat{\lambda} \le \lambda_0))} if `test` is "two-sided".
#' @returns `bf` Bayes Factor is calculated if `diagnosis = TRUE`, and the comparison model is non-informative prior, Jeffreys prior, \eqn{\pi \propto 1/\lambda}.
#' @references Jeffreys, H. (1946). An invariant form for the prior probability in estimation problems. Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences, 186(1007), 453-461.
#' @references Kass, R. E., & Raftery, A. E. (1995). Bayes factors. Journal of the american statistical association, 90(430), 773-795.
#' @examples
#' data <- Simulate_Enroll(n = c(50,20,20), lambda = 0.03, event = 0.1, M = 1, group = 3, maxt = 5,
#'                   accrual = 3, censor = 0.9, followup = 2,partition = "Uneven")
#' test <- Bayes_test(data, alpha = 3, beta = 82, test = "greater", pred = 2, threshold = 0.9,
#'                 type = "Predictive",diagnosis = T)
#' print(test)
#' @export
Bayes_test <- function(data,alpha,beta, test, threshold, type, pred, diagnosis = FALSE ){
  bf = NULL
  posterior1 <- .gamma_exponential(100000,data,alpha,beta)
  mean <- posterior1$mean
  posterior_CI <- posterior1$credible_interval
  MCMC <- posterior1$MCMC
  sd <- posterior1$sd
  if( diagnosis == TRUE){
    norm1 <- .norm(data,alpha,beta)
    norm0 <- .norm(data,0,0)
    bf <- norm1/norm0
  }
  ## lambda
  if (type == "Posterior"){

  CI <- posterior_CI

  z <- (mean(MCMC) - threshold )/sd

  if (test == "greater") output <-  mean(threshold<=MCMC)
  if (test == "less")  output <- mean(threshold>=MCMC)

  upper <- qnorm(0.975)*sd + threshold
  lower <- qnorm(0.025)*sd + threshold

  if (test == "two_sided") output <-  (mean(upper<=MCMC)+mean(lower>=MCMC))


  }
  ## survival
  if (type == "Predictive"){

  CI <-  sort( exp(-pred*posterior_CI)*100 , decreasing = F)
  names( CI ) <- c("2.5%", "5%", "25%", "50%", "75%", "95%", "97.5%")

  MCMC <- exp(-pred*MCMC)
  mean <- mean(MCMC)
  sd <- sd(MCMC)
  z <- (mean(MCMC) - threshold )/sd

  if (test == "greater") output <-  mean(threshold<=MCMC)
  if (test == "less")  output <- mean(threshold>=MCMC)

  upper <- qnorm(0.975)*sd + threshold
  lower <- qnorm(0.025)*sd + threshold

  if (test == "two_sided") output <-  (mean(upper<=MCMC)+mean(lower>=MCMC))

  }
  if(diagnosis == TRUE){
    output = list(mean = mean,
                  sd = sd,
                  CI = CI,
                  z_score = z,
                  prob = output,
                  bf = bf,
                  test = test,
                  type = type,
                  threshold = threshold)
  }else{
    output = list(mean = mean,
                  sd = sd,
                  CI = CI,
                  z_score = z,
                  prob = output,
                  test = test,
                  type = type,
                  threshold = threshold)
  }
  class(output) <- "Bayes_test"

  return(output)


}
#' @title Survival data simulation
#' @description \code{Simulate_Enroll} generates multiple streams of data sets with survival time, censoring status, and enrollment time.
#' @param n Integer. Sample size of patients
#' @param lambda Numerical range 0 and 1. Hazard rate of expoential distribution
#' @param event Numerical range 0 and 1. Event rate
#' @param M Integer. Number of trials generated for multiple streams of MCMC
#' @param group Integer. Number of subgroup for patient enrollment
#' @param maxt Numerical. The maximum time length of entire trial
#' @param accrual Numerical. The duration of patient enrolment
#' @param censor Numerical range 0 and 1. The censoring rate of patients leaving before trial ends.
#' @param followup Integer. The time length of follow up.
#' @param partition Logical. If `partition == "Even", the trial recruits equal numbers of patients in each stage; and if `partition == "Uneven", the trial recruits unequal numbers of patients in each stage.
#' @return Simulated survival data contain both survival time, censoring status, and enrollment time.
#' @examples
#' data <- Simulate_Enroll(n = c(50,20,20), lambda = 0.03, event = 0.1, M = 3, group = 3, maxt = 5,
#'                   accrual = 3,  censor = 0.9, followup = 2,partition = "Uneven")
#' head(data[[1]])
#' head(data[[2]])
#' head(data[[3]])
#' @export
Simulate_Enroll <- function(n,lambda,event,M,group,maxt,accrual,censor,followup,partition = "Even"){

  if(M != 1){
  data <- vector(mode = "list",length = M)
  for(m in 1:M){
    data[[m]] <- .simul(n,lambda,event,group,maxt,censor,accrual,followup,partition )
  }
  return(data)
  }
  if(M ==1){
    data <- .simul(n,lambda,event,group,maxt,censor,accrual,followup,partition )
  }
  return(data)
}
#' @title Bayesian adaptive trial interim analysis
#' @aliases BayesAT summary.BayesAT plot.BayesAT
#' @description \code{BayesAT} conducts Bayesian adaptive trials through multiple-stage interim analysis.
#' @param data Matrix. The data contains both survival time and event status.
#' @param D Numerical. The duration of interim analysis, matching the length of enrollment time.
#' @param stage Integer. Numbers of interim analysis stages.
#' @param threshold Numerical. The value tested against hypothesis or evidence.
#' @param start Numerical. The time point when the interim analysis starts.
#' @param objective Numerical. The time point for predicted survival rate, for example, 2 years, or 5 years survival probability.
#' @param alpha Numerical. Gamma distribution alpha parameter.
#' @param beta Numerical. Gamma distribution beta parameter (rate = 1/scale).
#' @param boundary The stopping criterion for interim analysis, and the default sets at 5% significance level and calculate quantiles by `qnorm()` for each stages.
#' @returns Interim analysis reporting Bayesian adaptive trial results.
#' @returns If there is one data set applied to `BayesAT`, the result will provide a table containing:
#' @returns `Upper bound` can be used as stopping criterion for efficacy;
#' @returns `Lower bound` can be used as stopping criterion for futility;
#' @returns `Z score` Z statistic is calculated based on the predicted survival probability:
#' @returns \deqn{\frac{\hat{S} - S_0}{SD( \hat{S} )}}
#' @returns with predicted mean survival rate \eqn{\hat{S}} and test evidence or threshold \eqn{S_0}.
#' @returns `Efficacy Prob` and `Futility Prob` Predictive probability measures the efficacy or futility, such as \eqn{P(\hat{S} > \text{Efficacy})} and \eqn{P(\hat{S} < \text{Futility})}.
#' @returns `Efficacy` and `Futility` indicate the interim analysis results: `+` means the trial reach the stopping criterion, otherwise it is `-`.
#' @examples
#' data <- Simulate_Enroll(n = c(30,20,20,15,30), lambda = 0.03, event = 0.1, M = 3, group = 5, maxt = 5,
#'                  accrual = 3,  censor = 0.9, followup = 2,partition = "Uneven")
#' ## assign patients in each group analyzed at each stage of time points
#' IA <- BayesAT(data,D = 3,stage = 5,threshold = 0.9, start = 1.5, objective = 2, alpha = 3, beta = 82)
#' summary(IA)
#' plot(IA)
#' @export
BayesAT <- function(data,D,stage,threshold,start,objective,alpha,beta,boundary = NULL) {
  if( is.null(boundary) ){
    boundary = matrix(0,nrow = 2, ncol = stage)
    boundary[1,] = qnorm(seq(1,0.95,length.out= stage) )
    boundary[1,1] = 4.3
    boundary[2,] = qnorm(seq(0,0.95,length.out= stage) )
    boundary[2,1] = -4.3
  }

  boundary <- round(boundary,4)

  rownames(boundary) <- c("Upper bound","Lower bound")

   if(class(data)=="data.frame"){

     M = 1
     output <- vector(mode = "list", length = M)
     IA <- .interim_test(data,D,stage,threshold, start,obj = objective,alpha,beta,boundary)
     z_score <- round(IA$z_score,4)
     #prob <- round(IA$prob,4)
     CP_efficay <- round(IA$eff,4)
     CP_futility <- round(IA$fut,4)
     Efficacy <- z_score - boundary[1,]
     Futility <- z_score - boundary[2,]
     Efficacy <- ifelse(Efficacy>=0,"+","-")
     Futility <- ifelse(Futility>=0,"-","+")
     output_m = rbind.data.frame(boundary,z_score,CP_efficay,CP_futility,
                                 Efficacy,Futility)
     colnames(output_m) <- paste("Stage", 1:stage)
     rownames(output_m) <- c("Upper bound","Lower bound","Z score",
                           "Efficacy Prob","Futility Prob",
                           "Efficacy", "Futility")
     output[[1]] <- output_m
     names(output)[1] <- paste("Interim analysis result for trial", M, sep = " ")
   }

   if(class(data)=="list"){
    M <- length(data)
    output <- vector(mode = "list", length = M)
    z_score <- prob <- NULL
    for(m in 1:M){
      data_m <- data[[m]]
      IA <- .interim_test(data_m,D,stage,threshold,start,obj = objective,alpha,beta,boundary)
      z_score <- round(c(IA$z_score),4)
      #prob  <- round(c(IA$prob),4)

      CP_efficay <- round(IA$eff,4)
      CP_futility <- round(IA$fut,4)
      Efficacy <- z_score - boundary[1,]
      Futility <- z_score - boundary[2,]
      Efficacy <- ifelse(Efficacy>=0,"+","-")
      Futility <- ifelse(Futility>=0,"-","+")
      output_m = rbind.data.frame(boundary,z_score,CP_efficay,CP_futility,
                                  Efficacy,Futility)
      colnames(output_m) <- paste("Stage", 1:stage)
      rownames(output_m) <- c("Upper bound","Lower bound","Z score",
                              "Efficacy Prob","Futility Prob",
                              "Efficacy", "Futility")
      output[[m]] <- output_m
      names(output)[m] <- paste("Interim analysis result for trial", m, sep = " ")
    }
   }


  class(output) <- "BayesAT"

  return(output)
}
#' @method summary BayesAT
#' @export
summary.BayesAT <- function(object){
  if (length(object) == 1){
    name <- names(object)
    assign(name, object[[1]])
    cat("\nInterim analysis results:\n")
    print(get(name[[1]]))
  }

  if (length(object) != 1){
    name <- names(object)
    cat("\nInterim analysis results:\n")
    for( m in 1:length(object)){

      assign(name[m], object[[m]])
      cat(paste( "\nTrial", m, ":\n"))
      print(get(name[[m]]))
    }
  }

}
#' @method plot BayesAT
#' @export
plot.BayesAT <- function(object){
  if (length(object) == 1){


    stage <- 1:ncol(object[[1]])
    x_poly <- c(stage, max(stage), min(stage))
    y_poly1 <-  c(as.numeric(object[[1]][1,]),5,5)
    y_poly2 <-  c(as.numeric(object[[1]][2,]),-5,-5)
    plot( stage , as.numeric(object[[1]][1,]), type = "l", ylim = c(-5,5),
          xlab = "Stage", ylab = "Z-scores"  )
    lines( stage,  as.numeric(object[[1]][2,]), type = "l")
    polygon(x_poly, y_poly1, col = "lightgreen", border = NA )
    polygon(x_poly, y_poly2, col = "pink", border = NA )
    lines(stage,as.numeric(object[[1]][3,]),type = "b",
          col = "black",pch=1,cex = 0.9)
    text(mean(stage),-4,"Futility")
    text(mean(stage),4,"Efficacy")

  }

  if (length(object) != 1){

    stage <- 1:ncol(object[[1]])
    x_poly <- c(stage, max(stage), min(stage))
    y_poly1 <-  c(as.numeric(object[[1]][1,]),5,5)
    y_poly2 <-  c(as.numeric(object[[1]][2,]),-5,-5)
    plot( stage , as.numeric(object[[1]][1,]), type = "l", ylim = c(-5,5),
          xlab = "Stages", ylab = "Normal critical values"  )
    lines( stage,  as.numeric(object[[1]][2,]), type = "l")
    polygon(x_poly, y_poly1, col = "lightgreen", border = NA )
    polygon(x_poly, y_poly2, col = "pink", border = NA )

    for (k in 1:length(object)){
      lines(stage,as.numeric(object[[k]][3,]),type = "b",
            col = k,pch=1,cex = 0.9)
    }

    trobjectl <- names(object)

    legend("bottomright",legend = trobjectl, cex = 0.7,
           lty = rep(1,length(object)), col = 1:length(object) ,pch= rep(1,length(object)) )


  }

}
#' @method summary Bayes_test
#' @export
summary.Bayes_test <- function(object){

  cat(  "\nBayesian inference can conclude that: \n"  )
  if( object$type == "Predictive" ){
    cat("Predictive survival mean is ", round(object$mean,4), "and standard deviation is ", round(object$sd,4),"\n")
    cat("Credible interval and summary statistics is given by \n")

    print(round(object$CI,2))
    cat("\n")
    if (object$test == "greater"){
      cat("Predictive probability P( S^ >=", object$threshold, ") = ", round(object$prob,4), "\n")
      cat("Standardized z score is", round(object$z_score,4), "\n")
    }
    if (object$test == "less"){
      cat("Predictive probability P( S^ < ", object$threshold, ") = ", round(object$prob,4), "\n")
      cat("Standardized z score is", round(object$z_score,4), "\n")
    }
    if (object$test == "two_sided"){
      cat("Predictive probability P( S^ != ", object$threshold, ") = ", round(object$prob,4), "\n")
      cat("Standardized z score is", round(object$z_score,4), "\n")
    }


  }
  if( object$type == "Posterior" ){
    cat("Posterior mean is ", round(object$mean,4), "and standard deviation is ", round(object$prob,4),"\n")
    cat("Credible interval and summary statistics is given by \n")

    print(round(object$CI,4))
    cat("\n")
    if (object$test == "greater"){
      cat("Posterior probability P( hr^ >= ", object$threshold, ") = ", round(object$prob,4), "\n")
      cat("Standardized z score is", round(object$z_score,4), "\n")
    }
    if (object$test == "less"){
      cat("Posterior probability P( hr^ < ", object$threshold, ") = ", round(object$prob,4), "\n")
      cat("Standardized z score is", round(object$z_score,4), "\n")
    }
    if (object$test == "two_sided"){
      cat("Posterior probability P( hr^ != ", object$threshold, ") = ", round(object$prob,4), "\n")
      cat("Standardized z score is", round(object$z_score,4), "\n")
    }


  }
  if ( !is.null(object$bf) ){
    cat("The Bayes factor of chosen prior is ", round(object$bf,4), " compared with non-informative prior. \n")
  }
}

