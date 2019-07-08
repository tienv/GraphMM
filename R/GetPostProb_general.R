GetPostProb_general = function(v, dat1.vec, dat2.vec, graph, 
                              folder, est_hyper, ListPara)
{ 
  #######################
  # Setup neighborhoods #
  #######################
  Nunits = graph[[1]]
  FullG = graph[[2]]
  
  # Finding neighbors of node 'v'
  ind1 = which(FullG[,1]==v)
  ind2 = which(FullG[,2]==v)
  neighbors = c(FullG[ind1,2],
                FullG[ind2,1])
  nbs = unique(neighbors)
  
  ###########################
  #  Setup graph structure  #
  ###########################
  v_list = sort(c(nbs, v))
  A = matrix((FullG %in% v_list), nrow = dim(FullG)[1], 
             ncol = dim(FullG)[2])
  ind3 = which(rowSums(A)==2)
  if (length(ind3)>1) SubG = FullG[ind3,] else 
  SubG = matrix(FullG[ind3,], nrow = length(ind3), ncol = dim(FullG)[2])

  SubG_reindex = matrix(match(SubG, v_list), 
                        nrow = dim(SubG)[1], 
                        ncol = dim(SubG)[2])
  center_gr_id = match(v, v_list)
  n.ver = length(v_list)
  n.ed = dim(SubG)[1]
  G.edgelist = cbind(SubG_reindex-1, (0:(n.ed-1)))
  
  if (n.ver > 12) {
      stop(paste("Error at node ", v, 
          "The number of vertices in the",
          "neighborhood must be smaller than 12. \n\n"))
  }
  
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
  Bell_number = c(1, 2, 5, 15, 52, 203, 877, 4140, 21147,
                  115975, 678570, 4213597)
  iteration = Bell_number[n.ver]*10       
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
  ListPara$dat1 = dat1.vec[, v_list]
  ListPara$dat2 = dat2.vec[, v_list]
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
