### Simulated stock composition data
## Generate stock composition data with various structural properties and then
# recover with models 
## Oct. 4 2019

require(tidyverse)
require(rethinking)


## First simulate five year's composition with three stocks
N <- 500 #total sample number
yrs <- seq(1, 5, by = 1)
N_yr <- max(yrs) #number of years
ppnStocks <- c(0.2, 0.4, 0.6)
K <- length(ppnStocks) #number of stocks

## Generate proportion for each year
dat <- data.frame(
  id = seq(1, N, by = 1),
  year = rep(yrs, each = 100),
  stock = NA
)
trueMat <- matrix(NA, nrow = N_yr, ncol = K + 1)
trueMat[ , 1] <- yrs
for (i in seq_along(yrs)) {
  error <- runif(nStocks, 0.0001, 0.9999)
  # generate error
  trueMat[i, 2:4] <- round(samSim::ppnAgeErr(ppnStocks, 0.1, error), digits = 2)
  simStocks <- sample(1:3, size = 100, prob = ppnI, replace = T)
  dat[dat$year == yrs[i], ]$stock <- simStocks
}

# Model with random effects by year
dat_list_yr <- list(
  K = K,
  N = N,
  N_yr = N_yr,
  y = dat$stock,
  yr = dat$year
)

# Model code for STAN
model_code_yr <- "
data{
    int N;
    int N_yr;
    int y[N];
    int yr[N];
    int K;
}
parameters{
    real a[K-1]; 						// intercepts for each stock, minus reference category
    matrix[K-1,N_yr] z_yr;      		// matrix of standardized random effects
    vector<lower=0>[K-1] sigma_yr;   	// stddev of random effects
    cholesky_factor_corr[K-1] L_Rho_yr; // correlation matrix of random effects, Choleskey decomposition
}
transformed parameters{
    matrix[N_yr,K-1] v_yr;     					// matrix of scaled random effects
    v_yr = (diag_pre_multiply(sigma_yr,L_Rho_yr) * z_yr)';   // note transpose in this transformation
}
model{

    // priors for fixed effects, mean followed by standard deviation
    a ~ normal(0,1);

	// hyper-priors
    to_vector(z_yr) ~ normal(0,1);
    sigma_yr ~ exponential(1);
    L_Rho_yr ~ lkj_corr_cholesky(2);

    // Likelihood function
    // This code sets up a function for each of the K-1 responses.
    // For each function (k), an intercept (a) is paramaterized along with
    // a subject-level varying intercept (v_yr). We use STAN's built-in categorical_logit
    // function for multinomial logistic regression.
    for ( i in 1:N ) {
        vector[K] p;
        for ( k in 1:(K-1) )
            p[k] = a[k] +
            v_yr[yr[i],k];
        p[K] = 0;
        y[i] ~ categorical_logit( p );
    }
}

	// In this block, we generate the variance-covariance matrix of year-level
	// random effects for the K-1 stocks. We then calculate the correlation between
	// these effects, Rho_yr, via a recomposition from the Cholesky matrix.
	// We also define a vector of length N for the log likelihood values, subsequently calling
	// STAN's categorical_logit_lpmf to generate the likelihood of each observation, conditional
	// on the model. Note that this step requires a repetition of the likelihood function, as above.
generated quantities{
    matrix[K-1,K-1] Rho_yr;
    vector[N] log_lik;
    Rho_yr = L_Rho_yr * L_Rho_yr';

    for ( i in 1:N ) {
        vector[K] p;
        for ( k in 1:(K-1) )
            p[k] = a[k] +
            v_yr[yr[i],k];
        p[K] = 0;
        log_lik[i] = categorical_logit_lpmf( y[i] | p );
    }
}
"

start_yr <- list (
  a = rep(0, K-1),
  sigma_yr = rep(1, K-1),
  L_Rho_yr = diag(K-1),
  z_yr = matrix(0, nrow=K-1, ncol=N_yr)
)

n_chains <- 3
init_yr <- list()
for ( i in 1:n_chains ) init_yr[[i]] <- start_yr

## We define a model fit object (mfit_i in this case), as is common with other model fitting functions
## in R.
mfit_yr <- stan(model_code=model_code_yr, data=dat_list_yr, chains=n_chains,
               cores= n_chains, warmup=1000, iter=2000, init=init_yr, 
               control = list(adapt_delta = 0.95))
precis(mfit_yr, depth = 2, prob = .96)

post <- extract.samples(mfit_yr)
stockB <- inv_logit(post$a[,1])
stockC <- inv_logit(post$a[,2])

precis( list(stockB = stockB, stockC = stockC) )
postcheck(mfit_yr, n = 1e4)



## brms example
dd <- data.frame(
  y1 = rbinom(N, 10, 0.3), y2 = rbinom(N, 10, 0.5), 
  y3 = rbinom(N, 10, 0.7), x = rnorm(N)
)
dd$size <- with(dd, y1 + y2 + y3)
dd$y <- with(dd, cbind(y1, y2, y3))
