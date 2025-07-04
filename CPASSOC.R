library(data.table)
library(dplyr)
library(tibble)
args <- commandArgs(trailingOnly=TRUE)
exp_path <- args[1]
data3<-fread(exp_path,header=T,sep="\t",check.names = FALSE)
data3<-data3%>%distinct(snpid,.keep_all=T)
data3<-data.frame(column_to_rownames(data3, var = "snpid"))
CorrMatrix=cor(data3)
# 
P1 <- as.numeric(args[2])
P2 <- as.numeric(args[3])
# 
#cat(sprintf("SampleSize = c(%d, %d)\n", P1, P2))
SampleSize <- c(P1, P2)
cat(sprintf("SampleSize = c(%d, %d)\n", SampleSize[1], SampleSize[2]))
library("MASS")
library("Matrix")
require(compiler)
enableJIT(3)
Non_Trucated_TestScore <- function(X, SampleSize, CorrMatrix)
{
  Wi = matrix(SampleSize, nrow = 1);
  sumW = sqrt(sum(Wi^2));
  W = Wi / sumW;
  
  Sigma = ginv(CorrMatrix);
  XX = apply(X, 1, function(x) {
    x1 <- matrix(x, ncol = length(x), nrow = 1);
    T = W %*% Sigma %*% t(x1);
    T = (T*T) / (W %*% Sigma %*% t(W));
    return(T[1,1]);
  }
  );
  return(XX);
}
SHom <- cmpfun(Non_Trucated_TestScore);

Trucated_TestScore <- function(X, SampleSize, CorrMatrix, correct = 1, startCutoff = 0, endCutoff = 1, CutoffStep = 0.05, isAllpoosible = T)
{
  N = dim(X)[2];
  
  Wi = matrix(SampleSize, nrow = 1);
  sumW = sqrt(sum(Wi^2));
  W = Wi / sumW;
  
  XX = apply(X, 1, function(x) {
    TTT = -1;
    
    if (isAllpoosible == T ) {
      cutoff = sort(unique(abs(x)));    ## it will filter out any of them.
    } else {
      cutoff = seq(startCutoff, endCutoff, CutoffStep);
    }
    
    
    for (threshold in cutoff) {
      x1 = x;
      index = which(abs(x1) < threshold);
      
      if (length(index) == N) break;
      
      A = CorrMatrix;
      
      W1 = W;
      if (length(index) !=0 ) {
        x1 = x1[-index];
        A  = A[-index, -index];   ## update the matrix
        W1 = W[-index];
      }
      
      if (correct == 1)
      {
        index = which(x1 < 0);
        if (length(index) != 0) {
          W1[index] = -W1[index];    ## update the sign
        }
      }
      
      A = ginv(A);
      x1 = matrix(x1, nrow = 1);
      W1 = matrix(W1, nrow = 1);
      T = W1 %*% A %*% t(x1);
      T = (T*T) / (W1 %*% A %*% t(W1));
      
      if (TTT < T[1,1]) TTT = T[1,1];
    }
    return(TTT);
  }
  );
  return(XX);
}
SHet <- cmpfun(Trucated_TestScore);

EstimateGamma <- function (N = 1E4, SampleSize, CorrMatrix, correct = 1, startCutoff = 0, endCutoff = 1, CutoffStep = 0.05, isAllpossible = T) {
  
  Wi = matrix(SampleSize, nrow = 1);
  sumW = sqrt(sum(Wi^2));
  W = Wi / sumW;
  
  Permutation = mvrnorm(n = N, mu = c(rep(0, length(SampleSize))), Sigma = CorrMatrix, tol = 1e-6, empirical = F);
  
  Stat =  TrucatedTestScore(X = Permutation, SampleSize = SampleSize, CorrMatrix = CorrMatrix,
                            correct = correct, startCutoff = startCutoff, endCutoff = endCutoff,
                            CutoffStep = CutoffStep, isAllpossible = isAllpossible);
  a = min(Stat)*3/4
  ex3 = mean(Stat*Stat*Stat)
  V =     var(Stat);
  
  for (i in 1:100){
    E = mean(Stat)-a;
    k = E^2/V
    theta = V/E
    a = (-3*k*(k+1)*theta**2+sqrt(9*k**2*(k+1)**2*theta**4-12*k*theta*(k*(k+1)*(k+2)*theta**3-ex3)))/6/k/theta
  }
  
  para = c(k,theta,a);
  return(para);
}
#shom=SHom(X=data3,SampleSize=SampleSize,CorrMatrix=CorrMatrix)
#shom=as.matrix(shom)
x=SHet(X=data3,SampleSize=SampleSize,CorrMatrix=CorrMatrix)#correct=1,isAllpossible=T
SHet=as.matrix(x)
write.table(SHet,'456.txt', quote = FALSE,row.names = T,sep="\t")
