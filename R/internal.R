
BayesAT <- function(data,D,stage,threshold,start,objective,alpha,beta,boundary = NULL) {
  if( is.null(boundary) ){
    boundary = matrix(0,nrow = 2, ncol = stage)
    boundary[1,] = qnorm(seq(1,0.95,length.out= stage) )
    boundary[1,1] = 3.9
    boundary[2,] = qnorm(seq(0,0.95,length.out= stage) )
    boundary[2,1] = -3.9
  }

  boundary <- round(boundary,4)

  rownames(boundary) <- c("Upper bound","Lower bound")

  if(class(data)=="data.frame"){

    M = 1
    output <- vector(mode = "list", length = M)
    IA <-  .interim_test(data,D,stage,threshold, start,obj = objective,alpha,beta,boundary)
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
      IA <-  .interim_test(data_m,D,stage,threshold,start,obj = objective,alpha,beta,boundary)
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



.gamma_exponential <- function(N,data,alpha,beta){

  ## data should be data.frame
  cumul_time <- sum(data$Time)
  cumul_event <- sum(data$Censor)

  ## posterior parameters
  post_alpha <- cumul_event + alpha
  post_beta <- cumul_time + beta

  ## MCMC
  MCMC_posterior <- rgamma(N,post_alpha,rate = post_beta)
  posterior_mean <- mean(MCMC_posterior)
  posterior_sd <- sd(MCMC_posterior)
  posterior_CI <-  sort(quantile(MCMC_posterior, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975) ),decreasing = F)
  names(posterior_CI ) <- c("2.5%", "5%", "25%", "50%", "75%", "95%", "97.5%")

  return(list(mean = c(posterior_mean),
              sd = c(posterior_sd),
              credible_interval = posterior_CI,
              MCMC = sort(MCMC_posterior, decreasing = F)))

}
.surv_data <- function(Time, Censoring){
  data <- cbind.data.frame(Time, Censoring)
  colnames(data) <- c("Time","Censor")

  if( length(unique(data$Censor) ) > 2 ) warning("Incorrect category of censoring status")
  if( length(unique(data$Censor) ) ==1 ) warning("Only one category of censoring status")

  return(data)
}
.norm <- function(data,alpha,beta,N =10^5 ){

  cumul_time <- sum(data$Time)
  cumul_event <- sum(data$Censor)
  if(alpha*beta !=0){
    z <-  beta^alpha/(cumul_time+beta)^(cumul_event+alpha)*gamma(cumul_event+alpha)/gamma(alpha)
  }
  if(alpha*beta==0){
    lambda <- seq(1/N,2, length.out = N )
    dens <- 1/lambda
    dens <- dens/sum(dens)
    z <- sum( lambda^cumul_event*exp(-lambda*cumul_time)*dens ) #gamma(cumul_event)/(cumul_time)^(cumul_event)
  }
  return(z)
}

.simul <- function(n,lambda,event,group,maxt,censor,accrual,followup,partition ="Even"){

  Time <- NULL
  C <- NULL
  Group <- NULL

  if(partition == "Even"){
    sub_sample_size <- rep(round(n/group), group)
    sub_sample_size[group] <- n - sum(sub_sample_size[1:(group-1)])
  }
  if(partition == "Uneven"){
    if (length(n) != group){
      warning("sub sample size not match number of group")
      sub_sample_size  <- n
    }
    if (length(n) == group){
      sub_sample_size  <- n
    }
  }


  for ( h in 1:(group-1)){

    t <- sort(rexp( round(sub_sample_size[h],0),lambda))

    ### this method applies fixed event rate
    c <- rep(0, length(t))
    c[1:round(event*length(t) )] <- 1

    ### truncate the maximum time
    t[t>=maxt] <- maxt  + c(rep(0, round(sum(t>=maxt)*censor)) ,
                            - runif( sum(t>=maxt) - round(sum(t>=maxt)*censor)))

    ### set enrollment to each group

    Group <- c(Group, runif(length(t),0,1)+h-1)
    Time <- c(Time, t)
    C <- c(C,c)
    maxt <-  followup + accrual /group*(group-h)
  }


  t <- sort(rexp( sub_sample_size[group],lambda))

  ### this method applies fixed event rate
  c <- rep(0,  length(t))
  c[1:round(event*length(t) )] <- 1

  ### truncate the maximum time
  t[t>=maxt] <- maxt  + c(rep(0, round(sum(t>=maxt)*censor)) ,
                          - runif( sum(t>=maxt) - round(sum(t>=maxt)*censor)))
  Group <- c(Group, runif(length(t),0,1)+h)
  Time <- c(Time, t)
  C <- c(C,c)
  Group <- Group*accrual/group

  data <- cbind.data.frame(Time, C,Group)
  colnames(data) <- c("Time","Censor","Enroll")
  return(data)
}

.interim_test <- function(data,D,stage,threshold,start,obj,alpha,beta,boundary) {

  time_point <- seq(0,D/(stage - 1)*stage,length.out = stage+1)
  index <- vector(mode = "list",length = stage)
  sub_n <- rep(0,stage)
  for ( h in 1:(stage)){
    index[[h]] <-  which(data$Enroll >= time_point[h] & data$Enroll < time_point[h+1])
    sub_n[h] <- length(which(data$Enroll >= time_point[h] & data$Enroll < time_point[h+1]))
  }

  stage_time <- stage_event <- matrix(0,nrow = stage,ncol = max(sub_n))

  for( h in 1:stage){
    stage_time[h,1:sub_n[h]] <- data$Time[index[[h]]]
    stage_event[h,1:sub_n[h]] <- data$Censor[index[[h]]]
  }

  alpha_k <- alpha
  beta_k <- beta
  ia <- prob <- eff_p <- fut_p <- numeric(stage)
  for ( k in 1:stage){
    updata_t <- stage_time[k,]
    updata_c <- stage_event[k,]

    updata_t[which(updata_t >=start)] <- start
    updata_c[which(updata_t >=start)] <- 0

    interim_data <- cbind.data.frame(updata_t,updata_c)
    colnames( interim_data) <- c("Time","Censor")

    interim_result <- Bayes_test(data = interim_data,alpha = alpha_k, beta = beta_k, test = "greater",
                                 threshold = threshold, type = "Predictive",pred = obj, diagnosis = FALSE)
    ia[k] <- interim_result$z_score
    prob[k] <-   interim_result$prob

    eff_boundary <- interim_result$sd*boundary[1,k] + threshold

    eff_result <- Bayes_test(data = interim_data,alpha = alpha_k, beta = beta_k, test = "greater",
                             threshold = eff_boundary, type = "Predictive",pred = obj, diagnosis = FALSE)
    eff_p[k] <- eff_result$prob

    fut_boundary <- interim_result$sd*boundary[2,k] + threshold
    fut_result <- Bayes_test(data = interim_data,alpha = alpha_k, beta = beta_k, test = "less",
                             threshold = fut_boundary, type = "Predictive",pred = obj, diagnosis = FALSE)
    fut_p[k] <- fut_result$prob

    bound <- sort( start + D/(stage-1)*(0:(k-1)), decreasing = T)
    history_t <- matrix(stage_time[1:k,],nrow = k )
    history_c <- matrix(stage_event[1:k,],nrow = k )

    for (step in 1:k){
      temp <- history_t[step, ]
      temp[which(temp >= bound[step])] <- bound[step]
      temp1 <- history_c[step,]
      temp1[which(temp >= bound[step])] <- 0
      history_t[step,] <- temp
      history_c[step,] <- temp1
    }

    alpha_k <- alpha + sum(history_c[1:k,])
    beta_k <- beta + sum(history_t[1:k,])


  }
  return(list( z_score = ia, prob = prob, eff = eff_p, fut =  fut_p))

}


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

