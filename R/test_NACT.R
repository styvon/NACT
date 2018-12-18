require("gee")
require("gmm")

# ---- Assisting Functions ------
get_mu <- function(k, x, g,
                   x_effect, g_effect){
  invlogit(
    logit(k) + 
      x_effect*(x>60) + 
      g_effect*g + 
      rnorm(length(g),0,0.01) # model with covariate effect (for X >60), genetic effect and residual error
  )
}

get_pheno <- function(mu){
  # function to generate phenotype from a list of prob(pheno=1)
  unlist(lapply(mu, rbinom, n=1, size=1))
}

invlogit <- function(x){
  # function for inverse logit
  exp(x)/(1+exp(x))
}

logit <- function(x){
  # function for logit
  log(x/(1-x))
}

gendata_prep <- function(f, q=1, k, maf, 
                    g_effect, x_effect){
  # Generate simulated data for pedigree with two founders and two off-springs
  # Used in gendata()
  # Input:
  #   f: int. number of independent pedigrees (families)
  #   q: int. number of genotypes (only q=1 implemented for now)
  #   k: num. param in Sara paper
  #   maf: num 0-1. allele frequency
  #   g_effect: num. genetic effect
  #   x_effect: num. covariate effect if covariate >60
  # Output: list
  #   G: matrix. genotype (0,1,2) of 4*f samples
  #   X: matrix. covariate binary (e.g.age>60) from norm(60, 15)
  #   pheno: matrix. binary phenotype of n+m samples
  #   K: matrix. kinship matrix of n+m samples.
  
  # number of subjects (generation2, offsprings)
  n <- f*2
  # number of relatives (generation 1, parents)
  m <- f*2
  ALLZERO <- TRUE
  FIRSTRUN <- TRUE
  while(ALLZERO){
    geno_g1 <- replicate(2*q, rbinom(m, 1, maf))
    geno_g2 <- matrix(rep(-1,n*2*q), ncol = 2*q)
    if(FIRSTRUN){ # if first run of the procedure, generate kinship matrix
      K <- diag(n+m)
      for(i in 1:f){
        for(j in 1:q){
          # sample one allele from paternal (odd row in geno_g1)
          geno_g2[(2*i-1):(2*i),2*q-1] = replicate(2,sample(geno_g1[2*i-1,(2*q-1):(2*q)],1, prob = c(0.5,0.5)))
          # sample one allele from maternal (even row in geno_g1)
          geno_g2[(2*i-1):(2*i),2*q] = replicate(2,sample(geno_g1[2*i,(2*q-1):(2*q)],1, prob = c(0.5,0.5)))
        }
        K[(2*i-1),c(2*i, n+2*i-1, n+2*i)] <- 0.5
        K[(2*i),c(2*i-1, n+2*i-1, n+2*i)] <- 0.5
        K[(n+2*i-1),c(2*i-1, 2*i, n+2*i)] <- 0.5
        K[(n+2*i),c(2*i-1, 2*i, n+2*i-1)] <- 0.5
      }
    }
    G_g1 <- matrix(apply(geno_g1, 1, sum),ncol=q)
    G_g2 <- matrix(apply(geno_g2, 1, sum),ncol=q)
    
    ALLZERO <- (sum(G_g2)==0)
    if(ALLZERO){
      print("No minor allele in sample, resampling...")
    }
  }
  
  G <- rbind(G_g2, G_g1)
  
  X <- matrix(round(rnorm(n+m, 60, 15),0),ncol=1)
  # X_g1 <- matrix(round(rtruncnorm(m, 45, 80, 60, 15),0),ncol=1) # covariate for relatives
  # X_g2 <- matrix(X_g1-30,ncol=1) # covariate for subjects
  # X <- rbind(X_g2, X_g1)
  
  mu <- get_mu(k, X, G[,1], x_effect, g_effect) # only consider one SNP scenario for now
  pheno <- matrix(get_pheno(mu),ncol=1)
  list(G=G,
       X=X>60,
       pheno=pheno,
       K=K)
  
}

gendata <- function(N, nsnp, k, maf, g_effect, x_effect, xi=0.9){
  # Generate simulated data for a total of 4*N subjects
  
  # Input:
  #   N: int. number of independent pedigrees (families)
  #   nsnp: int. number of total genotypes for generating GRM
  #   k: num. param in Sara paper
  #   maf: num 0-1. allele frequency
  #   g_effect: num. genetic effect of snp of interest (last column of data)
  #   x_effect: num. covariate effect if covariate >60
  #   xi: num. parameter for adjustment, another (1-xi) diagonal component will be added to  covariance matrix
  # Output: matrix with columns
  #   1: binary phenotype of 4*N samples
  #   2: covariate binary (e.g.age>60) from norm(60, 15)
  #   3-last: genotypes, with last column the genotype of interest
  
  n <- N*4 # total number of subjects
  
  # genotypes for generating GRM
  H1<-matrix(rbinom(N*(nsnp-1),1, maf), ncol=(nsnp-1))
  H2<-matrix(rbinom(N*(nsnp-1),1, maf),  ncol=(nsnp-1))
  H3<-matrix(rbinom(N*(nsnp-1),1, maf),  ncol=(nsnp-1))
  H4<-matrix(rbinom(N*(nsnp-1),1, maf),  ncol=(nsnp-1))
  H5<-matrix(rbinom(N*(nsnp-1),1, maf),  ncol=(nsnp-1))
  H6<-matrix(rbinom(N*(nsnp-1),1, maf),  ncol=(nsnp-1))
  
  G1<-H1+H2
  G2<-H3+H4
  G3<-H1+H3
  G4<-H1+H4
  
  Ga1<-rbind(G1,G2,G3,G4) # 1,2 parents; 3,4 child
  rm(list=c("H1","H2","H3","H4","H5","H6","G1","G2","G3","G4"))
  
  simdat <- gendata_prep(N,k=k, maf=maf, g_effect = g_effect,x_effect = x_effect)
  # Sigma1 <- cov(t(Ga1))* xi + (1-xi)*diag(n) # true matrix
  dat1 <- data.matrix(cbind(simdat$pheno, simdat$X,Ga1,simdat$G))
  rm(list=c("maf","n","k","N","nsnp","Ga1"))
  dat1
}

get_mu_gmm <- function(x, g,
                       x_effect, g_effect){
  # function for calculating mu with binary outcome
  # used in nact()
  # x: vector of covariates
  # x_effect: vector of covariate effects
  temp <- rowSums(cbind(x_effect[1], x_effect[2]*x, g_effect*g))
  invlogit(temp) 
}

moments <- function(beta, data) {
  # function to calculate moments
  # used in nact()
  # beta: coefficient of interest to be estimated
  # params: list (global)
  #   1: null estimate of nuisance coefficients
  #   2: list of inverse approximation matrices 
  # !!! paralell
  
  dim <- dim(data)
  y <- as.numeric(data[, 1])
  covariate <- as.numeric(data[, 2])
  g <- data.matrix(data[, dim[2]]) # last covariate as genotype
  mu <- get_mu_gmm(x=covariate, g=g, x_effect=params[[1]], g_effect=beta)
  mu0 <- get_mu_gmm(x=covariate, g=g, x_effect=params[[1]], g_effect=0) # null model estimates
  gamma <- sqrt(mu0 * (1-mu0))
  m <- c()
  for(i in 1:length(params[[2]])){ # global variable params
    # combine all formula with approximation matrix
    m <- c(m, t(g * gamma) %*% params[[2]][[i]] %*% as.matrix((y - mu)/gamma))
  }
  return(as.matrix(cbind(m))) # change to column vector
}

# ---- Major Functions ------
nact <- function(dat, n_components, n_approx, mfunc=moments){
  # Main function of NACT.
  # Inputs
  #   dat: matrix with rows as subjects, 1st column as phenotype, 2nd column as covariate, 3-last as genotype, last as genotype of interest
  #   n_components: int, number of sampled columns for Nystrom approximation of GRM
  #   n_approx: int, number of approximation matrices
  #   mfunc: function, e.g. see moments()
  # returns a gmm object
  
  # nuisance param estimate in null model
  params <<- list()
  
  id <- 1:dim(dat)[1]
  model0_gee <- gee(y1 ~ x, data = as.data.frame(cbind(id,y1=dat[,1],x=dat[,2],g=dat[,dim(dat)[2]])), 
                    family = binomial)
  params[[1]] <<- model0_gee$coefficients
  params[[2]] <<- list()
  rm(list=c("model0_gee","id"))
  
  # Nystrom approximation
  mymat <- cov(t(dat[,3:dim(dat)[2]]))* xi  # without adding diagonal
  inv_approxs <- list()
  n_samples <- dim(mymat)[1]
  n_snps <- dim(mymat)[2] 
  if(n_components > n_snps){
    cat("Number of components should not exceed number of SNPs, value is set to number of SNPs.")
    n_components <- n_snps
  }
  n_components <- min(n_components,n_samples)
  
  for(k in 1:n_approx){
    # get approximation matrices 
    # !!! parallel instead?
    ALLZERO <- TRUE
    while(ALLZERO){
      basis_inds <-  sample(1:n_snps,n_components, replace = FALSE)
      basis <-  mymat[,basis_inds]
      ALLZERO <- (sum(abs(mymat[basis_inds,basis_inds]))==0)
      if(ALLZERO){
        print(head(basis))
        print("All zeros in sample matrix. Resampling...")
      }
    }
   
    s <- svd(mymat[basis_inds,basis_inds])
    s$d[s$d<1e-12] <- 1e-12 # calibrate extreme small numbers
    w_k <-  (s$u) %*% (t(s$v) / (s$d))
    # mymat_approx <- basis  %*% w_k %*% t(basis) # reconstruct matrix
    
    # inverse using Woodbury
    params[[2]][[k]] <<- 1/(1-xi)*(diag(n_samples) - 
                                    basis %*% solve( (1-xi) * diag(n_components) +
                                                       w_k %*% t(basis) %*% basis) %*%
                                    w_k %*% t(basis) )
  }
  rm(list = c("basis_inds","basis","s","w_k"))
  
  # gmm
  init <- glm(dat[,1] ~ dat[,dim(dat)[2]],family = binomial)$coefficients[2] # glm estimate as initial value
  my_gmm <- gmm(g=mfunc, x = dat, 
                t0 = init, type = "iterative", crit = 1e-25, wmatrix = "ident",
                vcov="iid", prewhite=FALSE,
                method = "Brent",  lower=-100,upper=100,
                onlyCoefficients=FALSE,
                control = list(reltol = 1e-25, maxit = 20000))
  return(my_gmm)
}