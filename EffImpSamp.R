#Install and load relevant packages
install.packages("partitions")
library(partitions)


#### CODE FOR COUNTING HISTORIES#############

GroupSeq <- function(D)
{
  X <- list()
  cnt <- 0
  for ( i in 1:nrow(D[[1]]))
  {
    k <- cnt + 1
    cnt <- cnt + D[[2]][i]
    X[[i]] <- k:cnt
    }
 return(X)
}

GroupAl <- function(D)
{
  X <- list()
  for(i in 1:ncol(D[[1]]))
  {
    X[[i]] <- which(D[[1]][,i] == 1)
  }
  
  if(anyDuplicated(X))
  {
    X <- unique(X)
  }
  
  Y <- rep(1, length(X))
  for (i in 1:length(X))
  {
    for(j in setdiff(1:length(X), i))
    {
      if( all(X[[i]] %in% X[[j]]) )
      {
        Y[i] <- 0
        break
      }
    }
  }
  X[Y == 0] <- NULL
  return(X)
}

#Returns sequence IDs
seqId <- function(D)
{
  A <- 1:nrow(D[[1]])
  #A list to store subtrees
  ST <- list()
  GS <- GroupSeq(D)
  GA <- GroupAl(D)
  if(any(rowSums(D[[1]]) == 0))
  {
    GA[[length(GA) + 1]] <- setdiff(A, unlist(GA))
  }
  return(GA)
}

#Returns subtrees of a tree
SubTree <- function(D)
{
  ST <- list()
  if(class(D[[1]]) == "numeric")
  {
    ST[[1]] <- D
    return(ST)
  }
  A <- 1:nrow(D[[1]])
  #A list to store subtrees
  
  GS <- GroupSeq(D)
  GA <- seqId(D)
  
  for(i in 1:length(GA))
  {
    v <- GA[[i]]
    ST[[i]] <- list(D[[1]][v,], D[[2]][v])
  }
  
  for(j in 1:length(ST))
  {
    if(class(ST[[j]][[1]]) == "numeric")
    {
      del <- which(ST[[j]][[1]] == 0)
      if (length(del) == 0)
      {
        next
      }
      if(length(del) == length(ST[[j]][[1]]))
      {
        ST[[j]][[1]] <- 0
      }else{
        ST[[j]][[1]] <- ST[[j]][[1]][-del]
      }
      next
    }
    del <- which(colSums(ST[[j]][[1]]) == 0)
    if(length(del) == 0)
    {
      next
    }
    ST[[j]][[1]] <- ST[[j]][[1]][,-del]
  }
    return(ST)
}

#Returns number of histories for a sequence of n ancestral alleles
anc <- function(n)
{
  if(n==1)
  {
    return(1)
  }
  return(((n*(n-1))/2)*anc(n-1))
}

#Number of binary hierarchies
NBH <- function(r)
{
  if(r == 2)
  {
    return (1)
  }
  return((2*r-3)*NBH(r-1))
}

#Returns binary hierarchies
BH <- function(x)
{
  if(length(x) == 3)
  {
    xx <- restrictedparts(length(x), 2)
    parts <- listParts(xx)
    out <- rapply(parts, function(ii) x[ii], how="replace")
    out[[1]] <- NULL
    
  }else{
    xx <- restrictedparts(length(x), 2, include.zero = FALSE)
    parts <- listParts(xx)
    out <- rapply(parts, function(ii) x[ii], how="replace")
  }
  return(out)
}

#Returns number of nodes of a hierarchy
nodes <- function(H)
{
  n <- 0
  for(i in 1:length(H))
  {
    if(class(H[[i]][[1]]) == "numeric")
    {
      if(length(H[[i]][[1]]) == 1)
      {
        if(H[[i]][[1]] == 1)
        {
          n <- n + 1 + sum(H[[i]][[2]])
        }else{
          n <- n + sum(H[[i]][[2]])
        }
      }else{
        n <- n + sum(H[[i]][[1]]) + sum(H[[i]][[2]])
      }
    }else{
      n <- n + ncol(H[[i]][[1]]) + sum(H[[i]][[2]])
    }
  }
  return(n)
}

#Counts histories of trees with degree 1
count_deg1 <- function(D)
{
  #Subtrees
  ST <- SubTree(D)
  #Degree of root
  r <- length(ST)
  if(!is.matrix(ST[[1]][[1]]))
  {
    n <- ST[[1]][[2]]
    return(anc(n))
  }
  if(all(ST[[1]][[1]] == 1))
  {
    n <- ST[[1]][[2]]
    return(anc(n))
  }
  del <- which(colSums(ST[[1]][[1]]) == nrow(ST[[1]][[1]]))
  ST[[1]][[1]] <- ST[[1]][[1]][,-del]
  if(class(ST[[1]][[1]]) != "matrix")
  {
    ST[[1]][[1]] <- as.matrix(ST[[1]][[1]])
  }
  return(count(ST[[1]]))
}

#Counts hierarchies
count_H <- function(H)
{
  if(length(H) == 1)
  {
    return(count_deg1(H[[1]]))
  }
  R <- 1:length(H)
  B <- BH(R)
  v <- 0
  for(i in 1:length(B))
  {
    h1 <- B[[i]][[1]]
    h2 <- B[[i]][[2]]
    H1 <- H[h1]
    H2 <- H[h2]
    c_H1 <- count_H(H1) 
    c_H2 <- count_H(H2) 
    n1 <- nodes(H1)
    n2 <- nodes(H2)
    nodefactor <- choose(n1+n2-2, n1-1)
    v <- v + (c_H1*c_H2*nodefactor)
  }
  return(v)
}

#Counts histories for the S = (1,0), N = (1, f) case
spcl <- function(f)
{
  if(f == 1)
  {
    return (1)
  }
  return(anc(f+1) + spcl(f-1))
}

#Counts histories of a data set D consistent with the ISM
count <- function(D)
{
  #only ancestral alleles
  if(class(D[[1]]) == "numeric")
  {
    if(length(D[[1]]) == 1)
    {
      if(D[[1]] == 0)
      {
        a <- D[[2]]
        return(anc(a))
      }
    }
  }
  #base case(s)
  a <- as.matrix(c(1,0))
  b <- as.matrix(c(0,1))
  c <- c(1,1)
  #spcl
  if( (nrow(D[[1]]) == 2) & (ncol(D[[1]]) == 1))
  {
    if((all(D[[1]] == a) & (D[[2]][1] == 1)) |
    (all(D[[1]] == b) & (D[[2]][2] == 1)))
    {
      if( (all(D[[1]] == a)) )
      {
        f <- D[[2]][2]
      }else{
        f <- D[[2]][1]
      }
       return (spcl(f))
    }      
  }
  ST <- SubTree(D)
  seqST <- seqId(D)
  #Degree of root
  r <- length(ST)
  if(r==1)
  {
    return(count_deg1(D))
  }
  R <- 1:r
  B <- BH(R)
  v <- 0
  for(i in 1:length(B))
  {
    h1 <- B[[i]][[1]]
    h2 <- B[[i]][[2]]
    H1 <- ST[h1]
    H2 <- ST[h2]
    n1 <- nodes(H1)
    n2 <- nodes(H2)
    nodefactor <- choose(n1+n2-2, n1-1)
    c_H1 <- count_H(H1)  
    c_H2 <- count_H(H2) 
    v <- v + (c_H1*c_H2*nodefactor)
  }
  return(v)
}
#############################

#p(d_m) for Hobolth Scheme
p <- function(D, theta)
{
  n <- sum(D[[2]])
  s <- ncol(D[[1]])
  a <- nrow(D[[1]])
  ret <- rep(0, s)
  u <- matrix(0, a, s)
  for(i in 1:s)
  {
    f <- which(D[[1]][,i] == 1)
    d <- sum(D[[2]][f])
    k <- n-d+1
    
    if(d>1)
    {
      p_num <- function(theta, n, d, k)
      {
        if(k == 2)
        {
          return((d-1)/((1+theta)*(n-2)*(n-1)))
        }
        
        r <-  ((d-1)/(n-k))*
        (1/(k-1+theta))*(factorial(n-d-1)/
        (factorial(k-2)*factorial(n-d-k+1)))
        *(factorial(k-1)
        *factorial(n-k)/ (factorial(n-1)))
        return(r + p_num(theta, n, d, k-1))  
      }
      
      p_den <- function(theta, n, d, k)
      {
        if(k == 2)
        {
          return( 1/((n-1)*(1+theta)) )
        }
        
        r <- (1/(k-1+theta))*
        (factorial(n-d-1)/(factorial(k-2)*
        factorial(n-d-k+1)))*(factorial(k-1)*
        factorial(n-k)/(factorial(n-1)))
        return(r + p_den(theta, n, d, k-1))
      }
      
      ret[i] <- p_num(theta, n, d, k)
     /p_den(theta, n, d, k)
    }else
    {
      num <- 1/(n-1+theta)
      p1_den <- function(theta, n, k)
      {
        if(k==2)
        {
          return( (1/((1+theta)*(n-1))) )
        }
        
        r <- (k-1)/((k-1+theta)*(n-1))
        return(r + p1_den(theta, n, k-1))
      }
      
      den <- p1_den(theta, n, k = n)
      ret[i] <- num/den
    }
    
    for (j in 1:a)
    {
      if(D[[1]][j, i] == 1)
      {
        u[j, i] <- (ret[i]*D[[2]][j])/(d)
      }else {
        u[j, i] <- ((1 - ret[i])*D[[2]][j])/(n-d)
      }
    }
  }
  return(rowSums(u))
}

#Watterson's estimator for theta
tht <- function(D)
{
  s <- ncol(D[[1]]) #no. of segragating sites
  n <- sum(D[[2]]) #sample size
  
  den <- function(n)
  {
    if (n==1)
    {
      return(1)
    }
    return(1/n + den(n-1))
  }
  
  return(s/den(n-1))
}

#Identifies candidates eligible for
#involvement in previous event
candidates <- function(D)
{
  S <- D[[1]]
  N <- D[[2]]
  n <- sum(N)
  # no. of types of alleles
  k <- nrow(S)
  A <- 1:k
  # no. of seg sites
  s <- ncol(S)
  #Identifying allelic canidates for last event
  C <- which(N>1) l
  dm <- which(colSums(S) == 1) 
  dm2 <- (which(S[,dm] == 1))%%k 
  dm2[which(dm2==0)] = k
  M <- setdiff(dm2, C) 
  M1 <- rep(0, length(M))
  #in the sample
  nj <- rep(0, length(M))
  for(i in 1:length(M))
  {
    Dtemp <- D[[1]][,-(dm[which(D[[1]][M[i], dm] == 1)])]
    Dtemp <- as.matrix(Dtemp)
    if(length(dm[which(D[[1]][M[i], dm] == 1)]) > 1)
    {
      next
    }
    if(any(duplicated(Dtemp, MARGIN = 1)))
    {
      M1[i] <- M[i] 
      
      srch <- (1:k)[-M[i]]
      
      for(j in srch) #.. correct indexing
      {
        if(sum(Dtemp[j,] == Dtemp[M[i],]) == s-1)
        {
          nj[i] <- j
          break
        }
      }
    }
  }
  if(length(nj) > 1)
  {
    if(any(nj == 0))
    {
      nj <- nj[-which(nj==0)]
    }
  }
  M2 <- setdiff(M, M1)
  M1 <- setdiff(M,M2)
  
  Cand <- list()
  #Candidates
  Cand[[1]] <- C
  Cand[[2]] <- M1
  Cand[[3]] <- M2
  #Alleles corresponding to M1 
  Cand[[4]] <- nj
  #for use in parents() function
  Cand[[5]] <- dm
  
  return(Cand)
}

#Identifies parents of a given configuration
parents <- function(D, C, M1, M2, nj, dm)
{
  S <- D[[1]]
  N <- D[[2]]
  n <- sum(N)
  
  # no. of types of alleles
  k <- nrow(S)
  A <- 1:k
  # no. of seg sites
  s <- ncol(S)
  
  pnts <- list()
  D_C <- list()
  if(length(C) > 0)
  {
    N_C <- list()
    for(i in 1:length(C))
    {
      N_C[[i]] <- N
      N_C[[i]][C[i]] = N_C[[i]][C[i]] - 1
      D_C[[i]] <- list(S, N_C[[i]])
    }
  }
  
  D_M1 <- list()
  if(length(M1) > 0)
  {
    N_M1 <- list()
    S_M1 <- list()
    for(i in 1:length(M1))
    {
      N_M1[[i]] <- N
      S_M1[[i]] <- S
      mr <- dm[which(S[M1[i], dm] == 1)]
      
      S_M1[[i]] <- S_M1[[i]][,-mr] 
      #Not letting matrix to collapse into a vector
      if(class(S_M1[[i]])=="numeric")
      {
        S_M1[[i]] <- as.matrix(S_M1[[i]])
      }
      w <- nj[i]
      N_M1[[i]][w] <- N_M1[[i]][w] + 1
      N_M1[[i]] <- N_M1[[i]][-M1[i]]
      S_M1[[i]] <- S_M1[[i]][-M1[i],]
      if(class(S_M1[[i]])=="numeric")
      {
        S_M1[[i]] <- as.matrix(S_M1[[i]])
      }
      
      D_M1[[i]] <- list(S_M1[[i]], N_M1[[i]])
    }
  }
  
  
  D_M2 <- list()
  if(length(M2) > 0)
  {
    S_M2 <- list()
    for(i in 1:length(M2))
    {
      
      S_M2[[i]] <- S
      mr <- dm[which(S[M2[i], dm] == 1)]
      mr <- min(mr)
      
      S_M2[[i]] <- S[,-mr]
      if(class(S_M2[[i]])=="numeric")
      {
        S_M2[[i]] <- as.matrix(S_M2[[i]])
      }
      D_M2[[i]] <- list(S_M2[[i]], N) 
    }
  }
  
  
  pnts[[1]] <- D_C
  pnts[[2]] <- D_M1
  pnts[[3]] <- D_M2
  
  return(pnts)
}
#Probability of a sequence of only ancestral alleles
anc2 <- function(f, theta)
{
  if(f==1)
  {
    return(1)
  }
  return( ((f-1)/(f-1+theta))*anc2(f-1, theta) )
}

#Importance Sampling schemes
#Sampler = "GT", "SD", "Hob", or "CIS"
SIS <- function(D, theta, sampler)
{
  if ( !(sampler == "GT" || sampler == "SD" ||
        sampler == "HOB" || sampler == "CIS") )
    stop("Sampler must be either GT, SD, HOB, or CIS")
  if (theta == 0) {
    return(0)
  }
  
  type = sampler
  
  #Base case(s)
  
  #Only ancestral alleles remaining
  if(class(D[[1]]) == "numeric")
  {
    if(length(D[[1]]) == 1)
    {
      if(D[[1]] == 0)
      {
        return(anc2(D[[2]], theta))
      }
    }
  }
  
  #special case
  if( (nrow(D[[1]]) == 2) & (ncol(D[[1]]) == 1))
  {
    a <- as.matrix(c(1,0))
    b <- as.matrix(c(0,1))
    if((all(D[[1]] == a) & (D[[2]][1] == 1)) |
       (all(D[[1]] == b) & (D[[2]][2] == 1)))
    {
      if( (all(D[[1]] == a)) )
      {
        f <- D[[2]][2]
      }else{
        f <- D[[2]][1]
      }
      
      #Base case
      if(f == 1)
      {
        return(theta/((1+theta)^2))
      }
      
      PAR_C <- list()
      PAR_M1 <- list()
      PAR_C[[1]] <- a
      PAR_C[[2]] <- c(1, f-1)
      PAR_M1[[1]] <- 0
      PAR_M1[[2]] <- f+1
      Q <- rep(0, 2)
      
      if(sampler == "CIS")
      {
        Q[1] <- count(PAR_C)
        Q[2] <- count(PAR_M1)
      }
      
      if(sampler == "SD")
      {
        Q[1] <- f 
        Q[2] <- 1
      }
      
      if(sampler == "GT")
      {
        theta0 <- tht(D)
        Q[1] <- f-1
        Q[2] <- theta0
      }
      
      if(sampler == "HOB")
      {
        theta0 <- tht(D)
        U <- p(D, theta0)
        Q[1] <-  U[(which(D[[2]] > 1))]
        Q[2] <-  setdiff(U, Q[1])
      }
      
      Q <- Q/(sum(Q))
      smp <- sample(c(1,2), 1, prob = Q)
      Qw <- Q[smp]
      
      if(smp ==1)
      {
        P <- (f-1)/(f+theta)
        PAR <- PAR_C
      }else{
        P <- (theta)/(f+theta)
        PAR <- PAR_M1
      }
      
      return((P/Qw)*SIS(PAR, theta, type))
    }      
  }
  
  #Sample parameters
  S <- D[[1]]
  N <- D[[2]]
  n <- sum(N)
  ks <- nrow(S)
  A <- 1:ks
  
  #Eligibile candidates for
  #involvment in previous event
  cand <- candidates(D)
  C <- cand[[1]]
  M1 <- cand[[2]]
  M2 <- cand[[3]]
  nj <- cand[[4]]
  dm <- cand[[5]]
  #Parental configurations of D
  ps <- parents(D, C, M1, M2, nj, dm)
  D_C <- ps[[1]]
  D_M1 <- ps[[2]]
  D_M2 <- ps[[3]]
  #Proposal distribution
  Q <- rep(0, ks)
  PAR <- list()
  if(length(D_C)>0)
  {
    for(i in 1:length(D_C))
    {
      PAR[[C[i]]] <- D_C[[i]]
    }
  }
  
  if(length(D_M1)>0)
  {
    for(j in 1:length(D_M1))
    {
      PAR[[M1[j]]] <- D_M1[[j]]
    }
  }
  
  if(length(D_M2)>0)
  {
    for(k in 1:length(D_M2))
    {
      PAR[[M2[k]]] <- D_M2[[k]]
    }
  }
  
  if(sampler == "CIS")
  {
    if(length(D_C)>0)
      for(i in 1:length(D_C))
      {
        Q[C[i]] <- count(D_C[[i]])*choose(N[C[i]], 2)
      }
    
    if(length(D_M1)>0)
      for(j in 1:length(D_M1))
      {
        Q[M1[j]] <- count(D_M1[[j]])
      }
    
    if(length(D_M2)>0)
      for(k in 1:length(D_M2))
      {
        Q[M2[k]] <- count(D_M2[[k]])
      }
  }
  
  if(sampler == "SD")
  {
    if(length(C) > 0)
      Q[C] <- N[C]
    if(length(M1) > 0)
      Q[M1] <- 1
    if(length(M2) > 0)
      Q[M2] <- 1
  }
  
  if(sampler == "GT")
  {
    theta0 <- tht(D)
    if(length(C) > 0)
      Q[C] <- N[C] - 1
    if(length(M1) > 0)
      Q[M1] <- (theta0/n)*(1+nj)
    if(length(M2) > 0)
      Q[M2] <- (theta0/n)
  }
  
  if(sampler == "HOB")
  {
    theta0 <- tht(D) 
    U <- p(D, theta0)
    if(length(C) > 0)
      Q[C] <- U[C]
    if(length(M1) > 0)
      Q[M1] <- U[M1]
    if(length(M2) > 0)
      Q[M2] <- U[M2]
  }
  
  Q <- Q/(sum(Q))
  if (sum(is.na(Q)) > 0) {
    return(0)
  }
  smp <- sample(A, 1, prob = Q)
  Qw <- Q[smp]
  
  #Nominal distribution
  if(smp %in% C)
  {
    P <- ((N[smp]-1)/(n-1+theta))
  }
  
  if(smp %in% M1)
  {
    sc <- which(M1 == smp)
    P  <- (theta/(n*(n-1+theta)))
    *(N[nj[sc]] + 1)
  }
  
  if(smp %in% M2)
  {
    P <-  (theta/(n*(n-1+theta)))
  }
  return((P/Qw)*SIS(PAR[[smp]], theta, type))
}

#Parallelizing the importance sampling schemes
D <- list()
D[[1]] <- matrix(c(1,1,0,0,1,1,
0,1,0,0,1,0,0,0,0,0), 4,4, byrow = T)
D[[2]] <- c(1,1,2,1)
thetas <- c(0.5,2.12,5)
L_GT <- rep(0, length(thetas))
SAM <- c("GT", "SD", "HOB", "CIS")

 ########################################
# Experiments 1, 2
########################################

library(foreach)
library(doParallel)
library(Rmisc)

# Setup processor cores
cores <- detectCores()
n.cores <-  detectCores() - 1 # Number of cores to use
n.cores <- ifelse(n.cores > cores,
cores[1] - 1, n.cores)
cluster <- makeCluster(n.cores,
outfile="C:/Users/Drona/Desktop/Dissertation/tst1.txt") 
registerDoParallel(cluster)

tic <- proc.time()
par.m.2 <- foreach(a.index=1:length(SAM),
  .combine=rbind) %do% {
  library(partitions)  # Must be loaded on every core
  library(Rmisc)
  library(foreach)
  a <- SAM[a.index]
  foreach(i=1:n.cores) %dopar% {
    library(partitions)
    library(Rmisc)
  }
  results <- foreach(i=1:length(thetas),
  .combine=rbind) %dopar% {
    L_GT <- rep(0, 50000)
    L_GT <- sapply(1:50000, function(x)
    SIS(D, theta=thetas[i], sampler=a))
    
    print(paste0(
      "SAMPLER: ", a, "; ",
      "THETA: ", thetas[i], "; ",
      "APPROX LIK: ", mean(L_GT), ";",
      "ESS: ",
    ((50000)/(1+((sd(L_GT)/mean(L_GT))^2))), ";"
      
    ))
    
    data.frame(SAMPLER=a, THETA=thetas[i],
  APPROX_LIK=mean(L_GT), ESS = ((50000)/
  (1+((sd(L_GT)/mean(L_GT))^2))) )
  }
  results
}
toc2 <- proc.time() - tic

stopCluster(cluster)

#########################################################################
# SAMPLING PLOTS
#########################################################################

library(ggplot2)

par.m.2 <- par.m.2[par.m.2$APPROX_LIK > 0, ]
sampling.plot <- ggplot(par.m.2, 
aes(x=THETA, y=APPROX_LIK, col=SAMPLER)) +
geom_line() + geom_line(aes(y=LOWER, col = SAMPLER),
linetype="dotted") + geom_line(aes(y=UPPER, col = SAMPLER),
linetype="dotted") 

sampling.plot


#############################################
######## EGT Recursion ##############

#.............................................................
#The S = (1,0) special case 
spc <- function(f, theta)
{
  if(f==1)
  {
    return (theta/((1+theta)^2))
  }
  n <- f+1
  return( ((theta)/(n-1+theta))*anc2(f+1, theta)
+ ((n-1)/(n-1+theta))*((f-1)/f)*spc(f-1, theta)) #check
}

#Creates a unique id for data
dat.id <- function(D, theta) {
  id <- ""
  for (i in 1:length(D)) {
    Di <- paste(D[[i]], collapse=",")
    id <- paste0(id, "D[[", i, "]]=", Di, ";")
  }
  id <- paste0(id, "theta=", theta)
  return(id)
}

# Set global variable
lookup_lik <- list()


# Likelihood function
lik <- function(D, theta)
{
  
  #Base case(s)
  #Case when either M1, M2 or C is empty
  if(length(D) == 0)
  {
    return (0)
  }
  
  #Only anc
  if(class(D[[1]]) == "numeric")
  {
    if(length(D[[1]]) == 1)
    {
      if(D[[1]] == 0)
      {
        return(anc2(D[[2]], theta))
      }
    }
  }
  #spc
  if( (nrow(D[[1]]) == 2) & (ncol(D[[1]]) == 1))
  {
    a <- as.matrix(c(1,0))
    b <- as.matrix(c(0,1))
    if((all(D[[1]] == a) & (D[[2]][1] == 1))
    | (all(D[[1]] == b) & (D[[2]][2] == 1)))
    {
      if( (all(D[[1]] == a)) )
      {
        f <- D[[2]][2]
      }else{
        f <- D[[2]][1]
      }
      return (spc(f, theta))
      
    }      
  }
  
  #Sample parameters
  S <- D[[1]]
  N <- D[[2]]
  n <- sum(N)
  #Candidate
  cand <- candidates(D)
  C <- cand[[1]]
  M1 <- cand[[2]]
  M2 <- cand[[3]]
  nj <- cand[[4]]
  dm <- cand[[5]]
  #Parents
  ps <- parents(D, C, M1, M2, nj, dm)
  D_C <- ps[[1]]
  D_M1 <- ps[[2]]
  D_M2 <- ps[[3]]
  #Likelihood terms
  return( ((1)/(n-1+theta))*sum((N[C]-1)*
          vapply(D_C, function(x)
        {likM(x, theta)},
          FUN.VALUE = numeric(1))) + 
        (theta/(n*(n-1+theta)))*sum((N[nj]+1)*
          vapply(D_M1, function(y) {likM(y, theta)},
        FUN.VALUE = numeric(1))) + 
        (theta/(n*(n-1+theta)))*sum(vapply(D_M2,
        function(z){likM(z, theta)},
        FUN.VALUE = numeric(1))) )
}

#Memoized likelihood function
likM <- function(D, theta, lookup_lik.=lookup_lik) {
  
  # Get D, theta id
  D.theta.id <- dat.id(D, theta)
  print(D.theta.id)
  
  #check if likelihood of D already computed
  lookup_likelihood <- lookup_lik.[[D.theta.id]]
  if(length(lookup_likelihood) > 0) {
    print('likelihood memo')
    return(lookup_likelihood)
  }
  print('not likelihood memo')
  
  # Compute likelihood
  likelihood <- lik(D, theta)
  
  # Store value
  lookup_lik[[D.theta.id]] <<- likelihood
  return(likelihood)
}

#.............................................................................
#########################################################
#Experiment 3
#########################################################
cnt <- 0
lookup_pth <- list()
#Returns true probabilities of paths and identifies paths
cnt <- 0
pathProb <- function(D, theta, path='')
{
  cnt <<- cnt + 1
  
  pathTrue <- list()
  path <- paste0(path, 
'{{ ', dat.id(D, theta), ' }} -- ')
  #base cases
  #Only ancestral alleles remaining
  if(class(D[[1]]) == "numeric")
  {
    if(length(D[[1]]) == 1)
    {
      if(D[[1]] == 0)
      { 
        results <- list()
        results[[path]] <- anc2(D[[2]], theta)
        return(results)
      }
    }
  }
  
  #The S = (1,0), N = (1, f) case
  if( (nrow(D[[1]]) == 2) & (ncol(D[[1]]) == 1))
  {
    a <- as.matrix(c(1,0))
    b <- as.matrix(c(0,1))
    if((all(D[[1]] == a)
    & (D[[2]][1] == 1)) |
    (all(D[[1]] == b)
    & (D[[2]][2] == 1)))
    {
      if( (all(D[[1]] == a)) )
      {
        f <- D[[2]][2]
      }else{
        f <- D[[2]][1]
      }
      
      if(f==1)
      {
        results <- list()
        results[[path]] <- 
        theta/((1+theta)^2)
        return(results)
      }
      
      PAR_C <- list()
      PAR_M1 <- list()
      PAR_C[[1]] <- a
      PAR_C[[2]] <- c(1, f-1)
      PAR_M1[[1]] <- 0
      PAR_M1[[2]] <- f+1
      
      p.1 <- lapply(pathProbM(PAR_C, theta),
      function(x) {x * ((f-1)/(f+theta))})
      names(p.1) <- paste0(names(p.1), path)
      p.2 <- lapply(pathProbM(PAR_M1, theta),
      function(x) {x * (theta)/(f+theta)})
      names(p.2) <- paste0(names(p.2), path)
      results <- c(p.1, p.2)
      return(results)
    }
  }
  
  #Sample parameters
  S <- D[[1]]
  N <- D[[2]]
  n <- sum(N)
  ks <- nrow(S)
  A <- 1:ks
  #Eligibile candidates for involvment in previous event
  cand <- candidates(D)
  C <- cand[[1]]
  M1 <- cand[[2]]
  M2 <- cand[[3]]
  nj <- cand[[4]]
  dm <- cand[[5]]
  #Parental configurations of D
  ps <- parents(D, C, M1, M2, nj, dm)
  D_C <- ps[[1]]
  D_M1 <- ps[[2]]
  D_M2 <- ps[[3]]
  PAR <- list()
  i <- 0
  if(length(D_C)>0)
    for(i in 1:length(D_C))
    {
      PAR[[C[i]]] <- D_C[[i]]
    }
  j <- 0
  if(length(D_M1)>0)
    for(j in 1:length(D_M1))
    {
      PAR[[M1[j]]] <- D_M1[[j]]
    }
  k <- 0
  if(length(D_M2)>0)
    for(k in 1:length(D_M2))
    {
      PAR[[M2[k]]] <- D_M2[[k]]
    }
  
  i1 <- 0
  results_D_C <- list()
  if(length(D_C)>0)
    for(i1 in 1:length(D_C))
    {
      p <- lapply(pathProbM(PAR[[C[i1]]],
      theta), function(x)
      {x * ((N[C[i1]]-1)/(n-1+theta))})
      names(p) <- paste0(names(p), path)
      results_D_C[[i1]] <- p
    }
  results_D_C <- do.call(c, results_D_C)
  
  results_D_M1 <- list()
  j1 <- 0
  if(length(D_M1)>0)
    for(j1 in 1:length(D_M1))
    {
      p <- lapply(pathProbM(PAR[[M1[j1]]],
      theta), function(x)
      {x * ((N[nj[j1]]+1)/(n-1+theta))*
      (theta/n)})
      names(p) <- paste0(names(p), path)
      results_D_M1[[j1]] <- p
    }
  results_D_M1 <- do.call(c, results_D_M1)
  
  results_D_M2 <- list()
  k1 <- 0
  if(length(D_M2)>0)
    for(k1 in 1:length(D_M2))
    {
      p <- lapply(pathProbM(PAR[[M2[k1]]], theta), 
      function(x)
      {x * (theta/(n*(n-1+theta)))})
      names(p) <- paste0(names(p), path)
      results_D_M2[[k1]] <- p
    }
  results_D_M2 <- do.call(c, results_D_M2)
  
  
  results <- c(results_D_C, results_D_M1, results_D_M2)
  return(results)
  
}


#Memoized path prob function
pathProbM <- function(D, theta,
lookup_pth.=lookup_pth) {
  
  # Get D, theta id
  D.theta.id <- dat.id(D, theta)
  print(D.theta.id)
  
  #check if path prob of D already computed
  lookup_pth <- lookup_pth.[[D.theta.id]]
  if(length(lookup_pth) > 0) {
    print('path memo')
    return(lookup_pth)
  }
  print('not path memo')
  
  # Compute pathProb
  pP <- pathProb(D, theta)
  
  # Store value
  lookup_pth[[D.theta.id]] <<- pP
  return(pP)
}



lookup_Q <- list()
#Returns proposal based path probabilities
pathQ <- function(D, theta, sampler)
{
  if ( !(sampler == "GT" || sampler
== "SD" || sampler == "HOB" || 
         sampler == "CIS") )
    stop("Sampler must be either
  GT, SD, HOB, or CIS")
  pathTrue <- list()
  
  #base cases
  #Only ancestral alleles remaining
  if(class(D[[1]]) == "numeric")
  {
    if(length(D[[1]]) == 1)
    {
      if(D[[1]] == 0)
      { 
        return(1)
      }
    }
  }
  
  #The S = (1,0), N = (1, f) case
  if( (nrow(D[[1]]) == 2) &
  (ncol(D[[1]]) == 1))
  {
    a <- as.matrix(c(1,0))
    b <- as.matrix(c(0,1))
    if((all(D[[1]] == a) & 
    (D[[2]][1] == 1)) |
    (all(D[[1]] == b) 
    & (D[[2]][2] == 1)))
    {
      if( (all(D[[1]] == a)) )
      {
        f <- D[[2]][2]
      }else{
        f <- D[[2]][1]
      }
      
      if(f==1)
      {
        return(1)
      }
      
      PAR_C <- list()
      PAR_M1 <- list()
      PAR_C[[1]] <- a
      PAR_C[[2]] <- c(1, f-1)
      PAR_M1[[1]] <- 0
      PAR_M1[[2]] <- f+1
      Q <- rep(0, 2)
      
      if(sampler == "CIS")
      {
        Q[1] <- count(PAR_C)
        Q[2] <- count(PAR_M1)
      }
      
      if(sampler == "SD")
      {
        Q[1] <- f 
        Q[2] <- 1
      }
      
      if(sampler == "GT")
      {
        theta0 <- tht(D)
        Q[1] <- f-1
        Q[2] <- theta0
      }
      
      if(sampler == "HOB")
      {
        theta0 <- tht(D)
        U <- p(D, theta0)
        Q[1] <-  U[(which(D[[2]] > 1))]
        Q[2] <-  setdiff(U, Q[1])
      }
      
      Q <- Q/(sum(Q))
      pathTrue[[1]] <- Q[1]*
      pathQ(PAR_C, theta, sampler)
      pathTrue[[2]] <- Q[2]*
      pathQ(PAR_M1, theta, sampler)
      return(unlist(pathTrue))
    }
  }
  
  #Sample parameters
  S <- D[[1]]
  N <- D[[2]]
  n <- sum(N)
  ks <- nrow(S)
  A <- 1:ks
  #Eligibile candidates for involvment in previous event
  cand <- candidates(D)
  C <- cand[[1]]
  M1 <- cand[[2]]
  M2 <- cand[[3]]
  nj <- cand[[4]]
  dm <- cand[[5]]
  #Parental configurations of D
  ps <- parents(D, C, M1, M2, nj, dm)
  D_C <- ps[[1]]
  D_M1 <- ps[[2]]
  D_M2 <- ps[[3]]
  Q <- rep(0, ks)
  
  PAR <- list()
  i <- 0
  if(length(D_C)>0)
    for(i in 1:length(D_C))
    {
      PAR[[C[i]]] <- D_C[[i]]
    }
  j <- 0
  if(length(D_M1)>0)
    for(j in 1:length(D_M1))
    {
      PAR[[M1[j]]] <- D_M1[[j]]
    }
  k <- 0
  if(length(D_M2)>0)
    for(k in 1:length(D_M2))
    {
      PAR[[M2[k]]] <- D_M2[[k]]
    }
  
  if(sampler == "CIS")
  {
    if(length(D_C)>0)
      for(i in 1:length(D_C))
      {
        Q[C[i]] <- count(D_C[[i]])*choose(N[C[i]], 2)
      }
    
    if(length(D_M1)>0)
      for(j in 1:length(D_M1))
      {
        Q[M1[j]] <- count(D_M1[[j]])
      }
    
    if(length(D_M2)>0)
      for(k in 1:length(D_M2))
      {
        Q[M2[k]] <- count(D_M2[[k]])
      }
  }
  
  if(sampler == "SD")
  {
    if(length(C) > 0)
      Q[C] <- N[C]
    if(length(M1) > 0)
      Q[M1] <- 1
    if(length(M2) > 0)
      Q[M2] <- 1
  }
  
  if(sampler == "GT")
  {
    theta0 <- tht(D) #Updates theta 
    if(length(C) > 0)
      Q[C] <- N[C] - 1
    if(length(M1) > 0)
      Q[M1] <- (theta0/n)*(1+nj)
    if(length(M2) > 0)
      Q[M2] <- (theta0/n)
  }
  
  if(sampler == "HOB")
  {
    theta0 <- tht(D) #Updates theta
    U <- p(D, theta0)
    if(length(C) > 0)
      Q[C] <- U[C]
    if(length(M1) > 0)
      Q[M1] <- U[M1]
    if(length(M2) > 0)
      Q[M2] <- U[M2]
  }
  
  Q <- Q/(sum(Q))
  
  i1 <- 0
  if(length(D_C)>0)
    for(i1 in 1:length(D_C))
    {
      pathTrue[[i1]] <- Q[C[i1]]*
      pathQ(PAR[[C[i1]]], theta,sampler)
    }
  j1 <- 0
  if(length(D_M1)>0)
    for(j1 in 1:length(D_M1))
    {
      pathTrue[[i1 + j1]] <- Q[M1[j1]]*
    pathQ(PAR[[M1[j1]]], theta, sampler)
    }
  k1 <- 0
  if(length(D_M2)>0)
    for(k1 in 1:length(D_M2))
    {
      pathTrue[[i1 + j1 + k1]] <- 
    Q[M2[k1]]*pathQ(PAR[[M2[k1]]], 
    theta,sampler)
    }
  
  
  for(ol in 1:length(pathTrue))
  {
    if(is.null(pathTrue[[ol]]))
    {
      pathTrue[[ol]] <- NULL 
    }
  }
  
  return(unlist(pathTrue))
  
}

#Memoized path sample prob function
pathQM <- function(D, theta, sampler,
lookup_Q.=lookup_Q) {
  
  # Get D, theta id
  D.theta.id <- dat.id(D, theta)
  print(D.theta.id)
  
  #check if path prob of D already computed
  lookup_Q <- lookup_Q.[[D.theta.id]]
  if(length(lookup_Q) > 0) {
    print('path memo')
    return(lookup_Q)
  }
  print('not path memo')
  
  # Compute pathProb
  pQ <- pathQ(D, theta, sampler)
  
  # Store value
  lookup_Q[[D.theta.id]] <<- pQ
  return(pQ)
}
#plots
lookup_Q <- list()
par(mfrow=c(1,1))
pT <- pathProbM(D, theta)
pT <- pT/sum(pT)
plot( pT, pch = 16, xlab = 
"Path index", ylab = "GT and True Probabilities")
pGT <- pathQM(D, theta, "GT")
points(pGT, pch = 16, col = "purple")
plot(pT, pch = 16, xlab = "Path index",
ylab = "SD and True Probabilities")
lookup_Q <- list()
pSD <- pathQM(D, theta, "SD")
points(pSD, pch = 16, col = "red")
plot(pT, pch = 16, xlab = "Path index",
ylab = "Hob and True Probabilities")
lookup_Q <- list()
pHob <- pathQM(D, theta, "HOB")
points(pHob, pch = 16, col = "blue")
plot(pT, pch = 16, xlab = "Path index",
ylab = "CIS and True Probabilities")
lookup_Q <- list()
pCIS <- pathQM(D,theta, "CIS")
points(pCIS, pch = 16, col = "orange")
#KL divergence
lookup_pth <- list()
pT <- pathProbM(D, theta)
pT <- pT/sum(pT)
#GT
lookup_Q <- list()
pGT <- pathQM(D, theta, "GT")
KLGT<-sum((pT)*(log(pT) - log(pGT)))
#SD
lookup_Q <- list()
pSD <- pathQM(D, theta, "SD")
KLSD<-sum((pT)*(log(pT) - log(pSD)))
#Hob
lookup_Q <- list()
pHOB <- pathQM(D, theta, "HOB")
KLHOB<-sum((pT)*(log(pT) - log(pHOB)))
#CIS
lookup_Q <- list()
pCIS <- pathQM(D, theta, "CIS")
KLCIS<-sum((pT)*(log(pT) - log(pCIS)))
#L1 distance
#GT
L1GT <- sum(abs(pT-pGT))
#SD
L1SD <- sum(abs(pT-pSD))
#Hob
L1HOB <- sum(abs(pT-pHOB))
#CIS
L1CIS <- sum(abs(pT-pCIS))



