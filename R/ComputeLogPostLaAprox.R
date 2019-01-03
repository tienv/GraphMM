# require package trust
ComputeLF = function(Label, E, ListPara)
{
    
    dat1 = ListPara$dat1
    dat2 = ListPara$dat2
    matA = ListPara$matA
    matB = ListPara$matB
    mu0 = ListPara$mu0
    tau = ListPara$tau
    d0 = ListPara$d0
    delta = ListPara$delta
    deg = ListPara$deg
    PriorPara = ListPara$PriorPara
    PriorType = ListPara$PriorType
    
    LabelDiff = which(E==0)
    BlockSize = rowsum(rep(1,length(Label)), group = Label)
    m1 = dim(dat1)[1]
    m2 = dim(dat2)[1]
    NVer = length(Label)
    no_of_clusters = length(BlockSize)
    no_of_diff = length(LabelDiff)
    df = deg + NVer
    
    datmean1 = apply(dat1, MARGIN = 2, mean)
    datmean2 = apply(dat2, MARGIN = 2, mean)
    blocksize = rowsum(rep(1,length(Label)), group = Label)
    blockmean1 = rowsum(datmean1, group = Label)/BlockSize
    blockmean2 = rowsum(datmean2, group = Label)/BlockSize
    phi = blockmean1
    d = blockmean2 - blockmean1
    phi[which(E==1)] = (blockmean1[which(E==1)] + blockmean2[which(E==1)])/2
    d[which(E==1)] = 0
    InitValue = c(phi,d[LabelDiff])
    
    OptObj2 = trust::trust(OptimFunc_trust, InitValue, rinit=5, rmax = 5, minimize = F,
          Label = Label, E=E, ListPara=ListPara)
    x_max2 = OptObj2$argument
    f_max2 = OptimFunc(x_max2, Label, E, ListPara)

    eigval = eigen(-Hess_OptimFunc(x_max2, Label, E, ListPara))$values
    if (sum(eigval<=0)>0) 
    {
      print(eigval[which(eigval<=0)])
      print("Error: Hessian matrix has negative eigen values")
    }
    detMinusHessAprox = prod(abs(eigval))
    
    Loglkh = df/2*(log(det(matA))+log(det(matB))) - 
            no_of_clusters*log(tau*sqrt(2*pi)) - no_of_diff*log(delta*sqrt(2*pi)) + 
            (df+m1)*f_max2 + 0.5*(no_of_clusters+length(LabelDiff))*log(2*pi/(df+m1)) - 
            0.5*log(detMinusHessAprox)
    
    f = Loglkh;
    if (PriorType == 1)
    {
        for (k in 1:no_of_clusters)
        {
            f = f + lgamma(BlockSize[k])
        }
        f = f + no_of_clusters*log(PriorPara[1])
        
        w = sum(BlockSize[which(E==1)])
        f = f + w*log(PriorPara[2]) + (NVer-w)*log(1-PriorPara[2])
     }
    if (PriorType == 2)
    {
        for (k in 1:no_of_clusters)
        {
            f = f + lgamma(blocksize[k])
        }
        f = f + no_of_clusters*log(PriorPara[1])
        f = f + sum(E==1)*log(PriorPara[2]) + sum(E==0)*log(1-PriorPara[2])
    }
    if (PriorType == 3)
    {
        w = sum(BlockSize[which(E==1)])
        f = f + w*log(PriorPara[1]) + (NVer-w)*log(1-PriorPara[1])
    }
    if (PriorType == 4)
    {
        f = f + sum(E==1)*log(PriorPara[1]) + sum(E==0)*log(1-PriorPara[1])
    }
    return(f)
}
ComputeLogPost = function(partition, ListPara)
{
  return(ComputeLF(partition$GraphLabel, partition$E, ListPara))
  
}

