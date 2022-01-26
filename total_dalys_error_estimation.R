#implementation of uncertainty estimation for total DALYs for each genetic exposure

###inputs:
##data, rows correspond to genetic exposure-disease combinations for all disease that had unshrunk HR, includes columns:
#"rsid": identifier of the genetic exposure
#"cause": the disease
#"loghr_m": the log hazard ratio (meta-analyzed for FinnGen and UK Biobank)
#"loghr_se_m": the standard error for of log hazard ratios
##n_sample, how many times to resample
##C: the correlation matrix of betas, see supplementary information and get_covariance_matrix_of_betas.R for how this is estimated
##dalys_df, data frame containing population DALYs for Finland 2019 by "cause"
##geno, contains frequencies of the genetic exposures:
#"p0": proportion of unexposed
#"p1": proportion of exposed (e.g., 1 copy of common variant or in top 10% of a PGS)
#"p2": proportion of those with 2 copies of common variant (exposed), otherwise 0
##exposure_name, the name of the genetic exposure
##sex, sex to perform analysis on (both, male or female)

###main output:
##"prob" column of results_df is the posterior probability of the ALT hypothesis being true

resample_sum_daly <- function(data,
                              n_sample=10000,
                              C=corr_mat,
                              dalys_df = gbd_dalys_all,
                              geno=geno_freqs,
                              exposure_name,
                              sex) {
    # get causes that had shrunk coefficient
    causes <- data$cause
    # only keep correlation matrix of betas for those causes
    temp_C <- C[causes, causes]
    # calculate diagonal standard error matrix
    diag_se_mat <- data[,"loghr_se_m"] * diag(nrow(data))
    # calculate covariance matrix of multivariate normal to sample from
    cov_mat <- diag_se_mat  %*% as.matrix(temp_C) %*% diag_se_mat 
    # resample betas from multivariate normal where mean vector is the point estimates of coefficients
    # and covariance matrix is derived from standard errors for the diagonal, and from the correlation matrix
    # of null association coefficients for off-diagonal elements
    resample_df <- as.data.frame(mvrnorm(n=n_sample, mu=data[,"loghr_m"], Sigma=cov_mat))
    names(resample_df) <- data[,"cause"]
    resample_df$B <- 1:nrow(resample_df)
    resample_df <- resample_df %>%
    pivot_longer(names_to="cause", values_to="loghr", c(-(ncol(resample_df))))
    resample_df$p0 <- geno$p0
    resample_df$p1 <- geno$p1
    resample_df$p2 <- geno$p2
    resample_df$maf <- geno$maf
    resample_df$sex <- sex
    resample_df <- left_join(resample_df, dalys_df[,c("cause", "dalys", "sex")], by=c("cause", "sex"))
    
    #population numbers from GBD
    fin_2019_males <- 2730110
    fin_2019_females <- 2803985
    fin_2019_both <- 5534095
    fin_le_2019_males <- 79.16
    fin_le_2019_females <- 84.53
    fin_le_2019_both <- (fin_le_2019_males+fin_le_2019_females)/2

    pop_df <- data.frame(sex = c("both", "male", "female"),
                         population = c(fin_2019_both , fin_2019_males, fin_2019_females),
                        life_expectancy = c(fin_le_2019_both, fin_le_2019_males, fin_le_2019_females),
                        stringsAsFactors=FALSE)
    
    resample_df <- left_join(resample_df, pop_df, by="sex")

    resample_df <- resample_df %>%
        filter(!is.na(loghr)) %>%
        mutate(hr01 = exp(loghr),
                                hr02 = exp(2*loghr),
                                af1 = (p1*hr01) / (p0 + p1*hr01 + p2*hr02),
                                af2 = (p2*hr02) / (p0 + p1*hr01 + p2*hr02),
                                af0 = 1 - af1 - af2,
                                dalys0 = af0*dalys,
                                dalys1 = af1*dalys,
                                dalys2 = af2*dalys,
                                n0 = population*p0,
                                n1 = population*p1,
                                n2 = population*p2,
                                person_dalys0 = dalys0/(population*p0),
                                person_dalys1 = dalys1/(population*p1),
                                person_dalys2 = dalys2/(population*p2),
                                attr_person_dalys01 = (person_dalys1 - person_dalys0)*life_expectancy,
                                attr_person_dalys02 = (person_dalys2 - person_dalys0)*life_expectancy,
                                attr_person_dalys12 = (person_dalys2 - person_dalys1)*life_expectancy,
                                attr_population_dalys01 = attr_person_dalys01*n1/life_expectancy,
                                attr_population_dalys10 = -1*attr_person_dalys01*n0/life_expectancy,
                                attr_population_dalys02 = attr_person_dalys02*n2/life_expectancy,
                                attr_population_dalys20 = -1*attr_person_dalys02*n0/life_expectancy,
                                attr_population_dalys12 = attr_person_dalys12*n2/life_expectancy,
                                attr_population_dalys21 = -1*attr_person_dalys12*n1/life_expectancy,
                                attr_population_dalys0_12 = attr_population_dalys01 + attr_population_dalys02,
                                attr_population_dalys2_10 = attr_population_dalys21 + attr_population_dalys20)



    resample_estimates <- resample_df %>%
        group_by(B) %>%
        summarize(attr_person_dalys01=sum(attr_person_dalys01), .groups = 'drop') %>%
        ungroup() %>%
        pull(attr_person_dalys01)
    
     # calculate original estimate
    data$sex <- sex
    data <- left_join(data, dalys_df[,c("cause", "dalys", "sex")], by=c("cause", "sex"))
    data$p0 <- geno_freqs$p0
    data$p1 <- geno_freqs$p1
    data$p2 <- geno_freqs$p2
    data$maf <- geno_freqs$maf
    
    
    data <- left_join(data, pop_df, by="sex")

    data <- data %>%
        filter(!is.na(loghr_m)) %>%
        mutate(hr01 = exp(loghr_m),
                                hr02 = exp(2*loghr_m),
                                af1 = (p1*hr01) / (p0 + p1*hr01 + p2*hr02),
                                af2 = (p2*hr02) / (p0 + p1*hr01 + p2*hr02),
                                af0 = 1 - af1 - af2,
                                dalys0 = af0*dalys,
                                dalys1 = af1*dalys,
                                dalys2 = af2*dalys,
                                n0 = population*p0,
                                n1 = population*p1,
                                n2 = population*p2,
                                person_dalys0 = dalys0/(population*p0),
                                person_dalys1 = dalys1/(population*p1),
                                person_dalys2 = dalys2/(population*p2),
                                attr_person_dalys01 = (person_dalys1 - person_dalys0)*life_expectancy,
                                attr_person_dalys02 = (person_dalys2 - person_dalys0)*life_expectancy,
                                attr_person_dalys12 = (person_dalys2 - person_dalys1)*life_expectancy,
                                attr_population_dalys01 = attr_person_dalys01*n1/life_expectancy,
                                attr_population_dalys10 = -1*attr_person_dalys01*n0/life_expectancy,
                                attr_population_dalys02 = attr_person_dalys02*n2/life_expectancy,
                                attr_population_dalys20 = -1*attr_person_dalys02*n0/life_expectancy,
                                attr_population_dalys12 = attr_person_dalys12*n2/life_expectancy,
                                attr_population_dalys21 = -1*attr_person_dalys12*n1/life_expectancy,
                                attr_population_dalys0_12 = attr_population_dalys01 + attr_population_dalys02,
                                attr_population_dalys2_10 = attr_population_dalys21 + attr_population_dalys20)
    point_estimate <- data %>%
        summarize(attr_person_dalys01=sum(attr_person_dalys01)) %>%
        pull(attr_person_dalys01)
    
    df <- data.frame(rsid=exposure_name,
                     t0=point_estimate,
                     mean_t=mean(resample_estimates),
                     median_t=median(resample_estimates),
                     sd_t=sd(resample_estimates),
                     p025_t=quantile(resample_estimates, probs=c(0.025)),
                     p975_t=quantile(resample_estimates, probs=c(0.975)))
    
    rownames(df) <- NULL
    
    resample <- list(df=df, t=resample_estimates)
    
    return(resample)
}