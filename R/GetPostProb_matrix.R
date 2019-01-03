GetPostProb_matrix = function(v, dat1.vec, dat2.vec, size.im, folder, 
                              est_hyper, ListPara, IndicatorVector, IndexList, mccores)
{ 
  #######################
  # Setup neighborhoods #
  #######################
  dim = size.im;
  NbhSize = 3;
  CenterId = as.matrix(IndexList[v,])
  d = c(NbhSize, NbhSize)
  G1 = igraph::make_lattice(dimvector = d)
  NVer1 = igraph::vcount(G1)
  igraph::V(G1)$OrgId = (1:NVer1)
  c = ceiling(NbhSize/2)
  CenterArrId = c(c, c)
  CenterNodeId = (c-1)*d[1]+c
  rad = floor(NbhSize/2)
  NodeImageCoor = matrix(0, NVer1, length(dim))
  for (i in -rad:rad)
    for (j in -rad:rad)
    {
      arr_id = CenterArrId + c(i,j)
      node_id = arr_id[1] + (arr_id[2]-1)*d[1]
      NodeImageCoor[node_id,] = CenterId + c(i, j)
    }
  NodeIdRetain = c()
  for (i in 1:NVer1)
  {
    id = NodeImageCoor[i,]
    if((id[1]>=1)&&(id[2]>=1)&&(id[1]<=dim[1])&&(id[2]<=dim[2]))
    {
      vec_id = id[1] + (id[2]-1)*dim[1]
      if (IndicatorVector[vec_id]>0)
        NodeIdRetain = c(NodeIdRetain, i)
      
    }
  }
  G2 = igraph::induced_subgraph(G1, NodeIdRetain)
  cl = igraph::components(G2)
  if (cl$no > 1)
  {
    VCen = which(igraph::V(G2)$OrgId==CenterNodeId)
    Vlist = which(cl$membership==cl$membership[VCen])
    G = igraph::induced_subgraph(G2, Vlist)
  } else {
    G = G2 }
  NodeImageCoorRetain = NodeImageCoor[igraph::V(G)$OrgId,]
  # WriteMat(NodeImageCoorRetain-1, paste(foldername, "NodeImageCoor.txt", sep="/"))
  # WriteMat(get.edgelist(G, 1:ecount(G))-1, 
  #          paste(foldername, "SingleGraph.txt", sep="/"))
  node.coor = NodeImageCoorRetain
  G.nver = igraph::vcount(G)
  G.ned = igraph::ecount(G)
  Gr = igraph::get.edges(G, 1:G.ned)
  G.edgelist = cbind(Gr-1,(0:(G.ned-1)))
  n.ver = G.nver
  n.ed = G.ned
  
  center_gr_id = which((node.coor[,1] == CenterId[1])&(node.coor[,2] == CenterId[2]))
  
  seed1 = 65713501
  set.seed(seed1)
  seed = .Random.seed[2:626]
  
  datafolder = paste(folder, "/Neigh_X", CenterId[1], "_Y", CenterId[2], sep = "")
  dir.create(datafolder, recursive = T)
  WriteVec(seed, paste(datafolder, "Seed.txt", sep = "/"))
  WriteMat(G.edgelist, paste(datafolder, "Graph.txt", sep = "/"))
  
  #####################################################
  #  Sampling graph-respecting partition for 1 group  #
  #####################################################
  iteration = 20000         
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
  indlist = rep(-1, dim(node.coor)[1])
  for (i in 1:dim(node.coor)[1])
  {
    indlist[i] = which((IndexList[,1] == node.coor[i,1])&
                         (IndexList[,2] == node.coor[i,2]))
  }
  
  ListPara$dat1 = dat1.vec[, indlist]
  ListPara$dat2 = dat2.vec[, indlist]
  
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
  
  logpostprob = unlist(parallel::mclapply(AllPartitions, ComputeLogPost, ListPara, mc.cores = mccores))
  postprob = exp(logpostprob-max(logpostprob))
  postprob = postprob/sum(postprob)
  postEE = colSums(EEmat*postprob)
  return(postEE[center_gr_id])
}
