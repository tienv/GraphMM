OptimFunc = function(x, Label, E, ListPara)
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
    
    matS1 = cov(dat1)
    matT1 = cov(dat2)
    matS4 = 1/(m1-1)*matA + matS1
    matT4 = 1/(m2-1)*matB + matT1
    datmean1 = apply(dat1, MARGIN = 2, mean)
    datmean2 = apply(dat2, MARGIN = 2, mean)
    
    phi = x[1:no_of_clusters]
    d = rep(0, no_of_clusters)
    d[LabelDiff] = x[no_of_clusters+(1:no_of_diff)]
    
    priorSS2 =  sum((phi-mu0)^2)/2/tau^2/(df+m1) + sum((d[LabelDiff]-d0)^2)/2/delta^2/(df+m1) 
    muX = phi[Label]
    muY = phi[Label] + d[Label]
    A_tild = (matA + (m1-1)*matS1 + m1*(datmean1-muX)%*%t(datmean1-muX))
    B_tild = (matB + (m2-1)*matT1 + m2*(datmean2-muY)%*%t(datmean2-muY))
    f = -0.5*log(det(A_tild)) -0.5*(df+m2)/(df+m1)*log(det(B_tild)) - priorSS2
    return(f)
}

Grad_OptimFunc = function(x, Label, E, ListPara)
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
    
    matS1 = cov(dat1)
    matT1 = cov(dat2)
    matS4 = 1/(m1-1)*matA + matS1
    matT4 = 1/(m2-1)*matB + matT1
    C = matrix(0, NVer, no_of_clusters)
    for (i in 1:NVer)
        C[i,Label[i]] = 1
    datmean1 = apply(dat1, MARGIN = 2, mean)
    datmean2 = apply(dat2, MARGIN = 2, mean)
    
    phi = x[1:no_of_clusters]
    d = rep(0, no_of_clusters)
    if (no_of_diff>0)
        d[LabelDiff] = x[no_of_clusters+(1:no_of_diff)]
    
    vX = solve(matS4, datmean1)
    muX = phi[Label]
    wX = solve(matS4, muX)
    qX = t(C)%*%(wX-vX)
    vY = solve(matT4, datmean2)
    muY = phi[Label] + d[Label]
    wY = solve(matT4, muY)
    qY = t(C)%*%(wY-vY)
    A_tild = (matA + (m1-1)*matS1 + m1*(datmean1-muX)%*%t(datmean1-muX))
    B_tild = (matB + (m2-1)*matT1 + m2*(datmean2-muY)%*%t(datmean2-muY))
    
    rdetX = m1/(m1-1)*det(matA + (m1-1)*matS1)/det(A_tild)
    rdetY = m2/(m2-1)*det(matB + (m2-1)*matT1)/det(B_tild)
    
    grad = rep(0, no_of_clusters+no_of_diff)
    grad[1:no_of_clusters] = -rdetX*qX - (df+m2)/(df+m1)*rdetY*qY - 1/tau^2/(df+m1)*(phi-mu0)
    if (no_of_diff>0)
    {
        grad[no_of_clusters+(1:length(LabelDiff))] = -(df+m2)/(df+m1)*rdetY*qY[LabelDiff] - 
                                                    1/delta^2/(df+m1)*(d[LabelDiff]-d0)
    }
    
    return(grad)
}

Hess_OptimFunc = function(x, Label, E, ListPara)
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
    
    matS1 = cov(dat1)
    matT1 = cov(dat2)
    matS4 = 1/(m1-1)*matA + matS1
    matT4 = 1/(m2-1)*matB + matT1
    C = matrix(0, NVer, no_of_clusters)
    for (i in 1:NVer)
        C[i,Label[i]] = 1
    datmean1 = apply(dat1, MARGIN = 2, mean)
    datmean2 = apply(dat2, MARGIN = 2, mean)
    
    phi = x[1:no_of_clusters]
    d = rep(0, no_of_clusters)
    if (no_of_diff>0)
        d[LabelDiff] = x[no_of_clusters+(1:no_of_diff)]
    
    vX = solve(matS4, datmean1)
    muX = phi[Label]
    wX = solve(matS4, muX)
    qX = t(C)%*%(wX-vX)
    vY = solve(matT4, datmean2)
    muY = phi[Label] + d[Label]
    wY = solve(matT4, muY)
    qY = t(C)%*%(wY-vY)
    A_tild = (matA + (m1-1)*matS1 + m1*(datmean1-muX)%*%t(datmean1-muX))
    B_tild = (matB + (m2-1)*matT1 + m2*(datmean2-muY)%*%t(datmean2-muY))
    
    rdetX = m1/(m1-1)*det(matA + (m1-1)*matS1)/det(A_tild)
    rdetY = m2/(m2-1)*det(matB + (m2-1)*matT1)/det(B_tild)
    HX = t(C)%*%solve(matS4)%*%C
    HY = t(C)%*%solve(matT4)%*%C
    
    MinusHessMat = matrix(0, no_of_clusters+length(LabelDiff), no_of_clusters+length(LabelDiff))
    MinusHessMat[1:no_of_clusters, 1:no_of_clusters] = rdetX*(HX-2*rdetX*qX%*%t(qX)) +
      (df+m2)/(df+m1)*rdetY*(HY-2*rdetY*qY%*%t(qY)) + diag(1/tau^2/(df+m1), no_of_clusters, no_of_clusters)
    
    HYY = (df+m2)/(df+m1)*rdetY*(HY - 2*rdetY*qY%*%t(qY))
    if (no_of_diff>0)
    {
        MinusHessMat[no_of_clusters+(1:length(LabelDiff)), no_of_clusters+(1:length(LabelDiff))] = HYY[LabelDiff, LabelDiff] + 
                                                                                                    diag(1/delta^2/(df+m1), length(LabelDiff), length(LabelDiff))
        MinusHessMat[1:no_of_clusters, no_of_clusters+(1:length(LabelDiff))] = HYY[1:no_of_clusters,LabelDiff]
        MinusHessMat[no_of_clusters+(1:length(LabelDiff)), 1:no_of_clusters] = HYY[LabelDiff, 1:no_of_clusters]
    }
    return(-MinusHessMat)
    
}

OptimFunc_trust = function(x, Label, E, ListPara)
{
    f = OptimFunc(x, Label, E, ListPara)
    g = Grad_OptimFunc(x, Label, E, ListPara)
    HessMat = Hess_OptimFunc(x, Label, E, ListPara)
    
    return(list(value=f, gradient=g, hessian = HessMat))
}
