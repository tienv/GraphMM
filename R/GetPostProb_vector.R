GetPostProb_vector = function(v, dat1.vec, dat2.vec, nbh_size, 
                              folder, est_hyper, ListPara)
{ 
  #######################
  # Setup neighborhoods #
  #######################
  nbh_rad = floor(nbh_size/2)
  # nbs <- vector( mode="list", length=nrow(dat1.vec) )
  # center_gr_id = rep(NA, nrow(dat1.vec))
  nbs = c()
  if (v <= nbh_rad) {
    nbs = c(1:(v+nbh_rad))
    center_gr_id = v
  } else if ((v>nbh_rad)&&(v <= (ncol(dat1.vec)-nbh_rad))) {
    nbs = ((v-nbh_rad):(v+nbh_rad))
    center_gr_id = nbh_rad + 1
  } else if (v > (ncol(dat1.vec)-nbh_rad)) {
    nbs = ((v-nbh_rad):ncol(dat1.vec))
    center_gr_id = nbh_rad + 1
  }
  
  ###########################
  #  Setup graph structure  #
  ###########################
  Gr = igraph::make_graph(c(rbind(1:(length(nbs)-1),
                                  2:(length(nbs)))), directed = F)
  G.nver = igraph::vcount(Gr)
  G.ned = igraph::ecount(Gr)
  G = igraph::get.edges(Gr, 1:G.ned)
  G.edgelist = cbind(G-1,(0:(G.ned-1)))
  n.ver = G.nver
  n.ed = G.ned
  
  seed1 = 65713501
  set.seed(seed1)
  seed = .Random.seed[2:626]
  
  datafolder = paste(folder, "/Neigh",v, sep = "")
  dir.create(datafolder, recursive = T)
  WriteVec(seed, paste(datafolder, "Seed.txt", sep = "/"))
  WriteMat(G.edgelist, paste(datafolder, "Graph.txt", sep = "/"))
  
  #####################################################
  #  Sampling graph-respecting partition for 1 group  #
  #####################################################
  iteration = 100*(2^n.ver)         
  prior.type = 2   # uniform prior
  a = 2 # this is ignored in case of uniform prior
  L = SamplingPrior(iter = iteration, NVer = n.ver, 
                    NEd = n.ed, PriorType = prior.type, a=a,
                    datafolder = datafolder)
  UniqueStates = L$UniqueStates
  UniqueStates = UniqueStates+1
  
  # UniqueStates = as.matrix(data.table::fread(paste(datafolder, 
  #                 "UniqueStates.txt", sep = "/"), sep = " ", header = F,
  #                 skip = 1)) +1
  
  
  ################################
  #  Get partitions for 2 groups  #
  ################################
  count = 0
  for (i in 1:dim(UniqueStates)[1])
  {
    count = count + 2^max(UniqueStates[i,])
  }
  AllPartitions = vector(mode = "list", count)
  EEmat = matrix(0,count, n.ver)
  count1 = 0
  for (i in 1:dim(UniqueStates)[1])
  {
    count2 = 0
    
    count2 = count2+1
    AllPartitions[[count1+count2]]$GraphLabel = UniqueStates[i,]
    AllPartitions[[count1+count2]]$E = rep(1, max(UniqueStates[i,]))
    EEmat[count1+count2,] =  AllPartitions[[count1+count2]]$E[AllPartitions[[count1+count2]]$GraphLabel]
    for (j in max(UniqueStates[i,]):1)
    {
      for (k in 1:count2)
      {
        AllPartitions[[count1+count2+k]]$GraphLabel = UniqueStates[i,]
        AllPartitions[[count1+count2+k]]$E = AllPartitions[[count1+k]]$E
        AllPartitions[[count1+count2+k]]$E[j] = 0
        EEmat[count1+count2+k,] =  AllPartitions[[count1+count2+k]]$E[AllPartitions[[count1+count2+k]]$GraphLabel]
      }
      count2 = count2*2
    }
    count1 = count1+count2
  }
  
  ###########################
  #  Setup hyperparameters  #
  ###########################
  ListPara$dat1 = dat1.vec[, nbs]
  ListPara$dat2 = dat2.vec[, nbs]
  EstMatA = matrix(ListPara$mean.cor1, n.ver, n.ver)
  diag(EstMatA) = rep(1, n.ver)
  EstMatA = ListPara$mean.sd1^2*EstMatA
  
  EstMatB = matrix(ListPara$mean.cor2, n.ver, n.ver)
  diag(EstMatB) = rep(1, n.ver)
  EstMatB = ListPara$mean.sd2^2*EstMatB
  
  if (est_hyper == "global") {
    ListPara$matA = EstMatA
    ListPara$matB = EstMatB
  } else {
    dat.mean1 = apply(ListPara$dat1,2,mean)
    dat.mean2 = apply(ListPara$dat2,2,mean)
    dat.sum = c(apply(ListPara$dat1,2,sum),apply(ListPara$dat2,2,sum))
    dat.sd1 = apply(ListPara$dat1,2,sd)
    dat.sd2 = apply(ListPara$dat2,2,sd)
    
    ListPara$mu0 = mean(dat.mean1)
    ListPara$d0 = mean(dat.mean2 - dat.mean1)
    ListPara$tau = sd(dat.mean1)
    ListPara$delta = sd(dat.mean2 - dat.mean1)
    
    if (est_hyper == "local") {
      ListPara$matA = cov(ListPara$dat1)
      ListPara$matB = cov(ListPara$dat2)
    } else if (est_hyper == "mixed") {
      ListPara$matA = EstMatA
      ListPara$matB = EstMatB
    }
  } 
  
  logpostprob = unlist(lapply(AllPartitions, ComputeLogPost, ListPara))
  postprob = exp(logpostprob-max(logpostprob))
  postprob = postprob/sum(postprob)
  postEE = colSums(EEmat*postprob)
  return(postEE[center_gr_id])
}
