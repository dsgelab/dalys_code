#Shrinkage procedure to assign given log hazard ratio and it's standard error into NULL or ALT models
#where NULL is point mass at zero
#ALT is a 50%:50% mixture of Normals centered at +/- mu.1 and with sd of sd.1

###inputs:
##data, rows correspond to 80 genetic exposure-disease combinations, includes columns:
#"rsid": identifier of the genetic exposure
#"cause": the disease
#"loghr_m": the log hazard ratio (meta-analyzed for FinnGen and UK Biobank)
#"loghr_se_m": the standard error for of log hazard ratios
##mu.1, sd.1, alpha, beta: the prior distribution parameters
##niter: number of Markov chain Monte Carlo iterations

###main output:
##"prob" column of results_df is the posterior probability of the ALT hypothesis being true

get_null_alt_data <- function(data, mu.1, sd.1, alpha, beta, niter) {
    

    b <- data$loghr_m
    se <- data$loghr_se_m
    
    K = length(b)

    #prior on prop.1 is Beta(a.1, a.0)
    a.1 = alpha
    a.0 = beta

    #initialize
    burnin = 50 #will be discarded as the initial phase
    prop.1 = rbeta(1, a.1, a.0) #sample initial value
    z = sample(c(0,1), prob = c(1-prop.1, prop.1), size = K, replace = T) #initial assignment

    z.count = rep(0,K) #to collect how many times each trait will be in ALT during MCMC
    prop.1.res = rep(NA, niter - burnin) #to collect posterior for prop.1

    for(iter in 1:niter){
      prop.1 = rbeta(1, a.1 + sum(z==1), a.0 + sum(z==0)) #sample prop.1

      #To sample the group assignments, compute marginal likelihoods for ALT:
      A = log(0.5) + cbind(dnorm(b, mu.1, sqrt(se^2 + sd.1^2), log = T),
                           dnorm(b, -mu.1, sqrt(se^2 + sd.1^2), log = T))
      A.max = apply(A, 1, max)
      logmlk.1 = log(rowSums(exp(A - A.max))) + A.max #log of marg. likelihood for ALT

      #compute also NULL model marg. likelihood and account for prop.1:
      B = cbind(log(prop.1) + logmlk.1, log(1-prop.1) + dnorm(b, 0, se, log = T))
      B.max = apply(B, 1, max)
      C = exp(B - B.max)
      C = C / rowSums(C) #C has the probabilities for ALT (in col1) and for NULL (in col2)

      z = as.numeric(runif(K) < C[,1]) #sample new group assignment

      if(iter > burnin){ #if > burnin, then collect results
        z.count[z == 1] = z.count[z == 1] + 1
        prop.1.res[iter-burnin] = prop.1
      }

    }

    z.res = z.count/(niter-burnin)
    
    results_df <- data[,c("rsid", "cause")]
    results_df$loghr_m <- data$loghr_m
    results_df$loghr_se_m <- data$loghr_se_m
    results_df$loghr_p_m <- data$loghr_p_m
    results_df$mu.1 <- mu.1
    results_df$sd.1 <- sd.1
    results_df$alpha <- alpha
    results_df$beta <- beta
    
    results_df$prob <- z.res
    
    return(results_df)
    
}