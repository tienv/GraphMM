#include "Functions.h"

int Sample(int size, arma::vec& prob_raw, double sum_prob, arma::uvec& index, Rand* rand)
{
  int n = prob_raw.n_elem;
  int nOrig_1 = n-1;
  double rT, mass, totalmass = 1.0;
  arma::vec prob(prob_raw.n_elem);
  for (int i=0; i<prob_raw.n_elem; i++) prob[i] = prob_raw[i]/sum_prob;
  arma::uvec perm = arma::sort_index(prob, 1); /*descending sort of index*/
  prob = arma::sort(prob, 1);  /*descending sort of prob*/ 
  int ii, jj, kk;
  for (ii = 0; ii < size; ii++, nOrig_1--)
  {
    double random = rand->runif();
    rT = totalmass * random;
    mass = 0;
    for (jj = 0; jj < nOrig_1; jj++) 
    {
      mass += prob[jj];
      if (rT <= mass)
        break;
    }
    index[ii] = perm[jj];
    totalmass -= prob[jj];
    for (kk = jj; kk < nOrig_1; kk++) 
    {
      prob[kk] = prob[kk+1];
      perm[kk] = perm[kk+1];
    }
  }
  return(size); 
}

int FindComponents(const arma::umat& G, int NVer, const std::vector<int>& EdValue, 
                   std::vector<int>& membership, std::vector<int>& csize)
{
  arma::uvec EdVal = arma::conv_to<arma::uvec>::from(EdValue) ;
  arma::umat G1 = G.rows(find(EdVal>0));
  
  
  std::vector<int> next_nei(NVer,0);
  std::vector<int> q;
  int Nq = 0;
  int no_of_clusters = 1;
  int act_cluster_size;
  std::vector<int> out;
  int Nout=0;
  out.reserve(NVer);
  
  /* Finding neighbors of each node */
  std::vector<std::vector<int> > neigh(NVer);
  int E = G1.n_rows;
  for (int i=0; i<NVer; i++)
  {
    neigh[i].reserve(std::min(E,50));
  }
  
  for (int j=0; j<E; j++)
  {
    neigh[G1(j,0)].push_back(G1(j,1));
    neigh[G1(j,1)].push_back(G1(j,0));
  }
  
  std::vector<int> tmp;
  for (int i=0; i<NVer; i++)
  {
    /* Finding neighbors of node i */
      //tmp = sort(neigh[i], "descend");
      tmp = neigh[i];  
      //std::sort(tmp.begin(), tmp.end());
      
      if (next_nei[i] > tmp.size()) continue;
      q.push_back(i);
      Nq++;
      while (Nq > 0)
      {
        int act_node = q[Nq-1];
        
        /* Finding neighbors of node act_node */
          // tmp = sort(neigh[act_node], "descend");
          tmp = neigh[act_node];
          //std::sort(tmp.begin(), tmp.end());
          
          if (next_nei[act_node] == 0)
          {
            next_nei[act_node]++;
          } else if (next_nei[act_node]<=tmp.size()) {
            int neighbor = tmp[next_nei[act_node]-1];
            if (next_nei[neighbor]==0)
            {
              q.push_back(neighbor);
              Nq++;
            }
            next_nei[act_node]++;
          } else {
            out.push_back(act_node);
            Nout++;
            q.pop_back();
            Nq--;
          }
      }
  }
  
  
  std::fill(next_nei.begin(),next_nei.end(),0); 
  while (Nout > 0)
  {
    int grandfather = out[Nout-1];
    out.pop_back();
    Nout--;
    if (next_nei[grandfather]!=0) continue;
    next_nei[grandfather] = 1;
    act_cluster_size = 1;
    membership[grandfather] = no_of_clusters-1;
    q.push_back(grandfather);
    Nq++;
    while (Nq>0)
    {
      int act_node = q[Nq-1];
      q.pop_back();
      Nq--;
      
      /* Finding neighbors of node act_node */
        //tmp = sort(neigh[act_node], "descend");
        tmp = neigh[act_node];
        //std::sort(tmp.begin(), tmp.end());
        for (int j=0; j<tmp.size(); j++)
        {
          int neighbor = tmp[j];
          if (next_nei[neighbor]!=0) continue;
          q.push_back(neighbor);
          Nq++;
          next_nei[neighbor] = 1;
          act_cluster_size++;
          membership[neighbor] = no_of_clusters-1;
        }
    }
    no_of_clusters++;
    csize.push_back(act_cluster_size);
    
  }
  no_of_clusters--;
  
  LabelNormalize(membership, csize);
  
  return(no_of_clusters);
}

void LabelNormalize (vector<int> & label, vector<int> & csize)
{
    vector<int> csize_new(csize.size(),-1);
    map<int, int> IDtoSeq;
    int index = 0;
    for (int i=0; i<label.size(); i++)
    {
        if (IDtoSeq.find(label[i]) == IDtoSeq.end())
        {
            csize_new[index] = csize[label[i]];
            IDtoSeq[label[i]] = index++;
        }
        label[i] = IDtoSeq[label[i]];
    }
    csize = csize_new;
}

void Exchange (std::vector<int>& EdValue1, std::vector<int> & EdValue2, std::vector<int> & EdValue_tmp,
               std::vector<int>& s1, std::vector<int>& s2, std::vector<int>& s_tmp,
               std::vector<int>& tree1, std::vector<int>& tree2, std::vector<int>& tree_tmp,
               double& f1, double& f2, double& f_tmp,
               double& Loglkh1, double& Loglkh2, double& Loglkh_tmp,
               double& log_n_tree1, double& log_n_tree2, double& log_n_tree_tmp,
               int& no_of_clusters1, int& no_of_clusters2, int& no_of_clusters_tmp, int& flag)
{
  if (flag == 1)
  {
    for (int j=0; j<EdValue1.size(); j++) EdValue1[j] = EdValue_tmp[j];
    for (int j=0; j<s1.size(); j++) s1[j] = s_tmp[j];
    for (int j=0; j<tree1.size(); j++) tree1[j] = tree_tmp[j];
    f1 = f_tmp;
    Loglkh1 = Loglkh_tmp;
    log_n_tree1 = log_n_tree_tmp;
    no_of_clusters1 = no_of_clusters_tmp;  
  }
  if (flag == 0)
  {
    /* Save the value to _tmp */
      for (int j=0; j<EdValue1.size(); j++) EdValue_tmp[j] = EdValue1[j];
      for (int j=0; j<s1.size(); j++) s_tmp[j] = s1[j];
      for (int j=0; j<tree1.size(); j++) tree_tmp[j] = tree1[j];
      f_tmp = f1;
      Loglkh_tmp = Loglkh1;
      log_n_tree_tmp = log_n_tree1;
      no_of_clusters_tmp = no_of_clusters1; 
      
      /* Assign new values */
        for (int j=0; j<EdValue1.size(); j++) EdValue1[j] = EdValue2[j];
        for (int j=0; j<s1.size(); j++) s1[j] = s2[j];
        for (int j=0; j<tree1.size(); j++) tree1[j] = tree2[j];
        f1 = f2;
        Loglkh1 = Loglkh2;
        log_n_tree1 = log_n_tree2;
        no_of_clusters1 = no_of_clusters2;  
        
        flag = 1;
  }
}

void SampleTree(const arma::umat& Gr, int n_ver, std::vector<int> & tr, Rand* rand)
{
   if (n_ver <= 1) return;
  
  /* Finding neighbors and incident edges of each node */
  std::vector<std::vector<int> > inci(n_ver);
  std::vector<std::vector<int> > neigh(n_ver);
  int E = Gr.n_rows;
  for (int i=0; i<n_ver; i++)
  {
    neigh[i].reserve(std::min(E,50));
    inci[i].reserve(std::min(E,50));
  }
  
  for (int j=0; j<Gr.n_rows; j++)
  {
    neigh[Gr(j,0)].push_back(Gr(j,1));
    neigh[Gr(j,1)].push_back(Gr(j,0));
    inci[Gr(j,0)].push_back(j);
    inci[Gr(j,1)].push_back(j);
  }
  
  
  /* Finding neighbors and incident edges of each node */
  // std::vector<std::vector<int> > inci(n_ver);
  // std::vector<std::vector<int> > neigh(n_ver);
  // 
  // for (int j=0; j<n_ver; j++)
    // {
      //   arma::uvec v0 = find(Gr.col(0)==j);
      //   arma::uvec v1 = find(Gr.col(1)==j);
      //   std::vector<int> tmp_neigh(v0.n_elem+v1.n_elem);
      //   std::vector<int> tmp_inci(v0.n_elem+v1.n_elem);
      //   for (int l=0; l<v0.n_elem; l++)
        //   {
          //     tmp_neigh[l] = Gr.col(1)[v0[l]];
          //     tmp_inci[l] = v0[l];
          //   }
      //   for (int l=0; l<v1.n_elem; l++)
        //   {
          //     tmp_neigh[l+v0.n_elem] = Gr.col(0)[v1[l]];
          //     tmp_inci[l+v0.n_elem] = v1[l];
          //   }
      // 
        //   neigh[j] = tmp_neigh;
        //   inci[j] = tmp_inci;
        // 
          // }
  
  int root = rand->rint(n_ver);

  /* Wilson algorithm */
    
    std::vector<int> NextVer(n_ver,-1);
  std::vector<int> NextEd(n_ver,-1);
  std::vector<bool> intree(n_ver,false);
  intree[root] = true;
  int uu;
  
  for (int i = 0; i<n_ver; i++)
  {
    /*Find a loop-erased random walk from vertex i to the current tree*/
      uu = i;
      while (intree[uu]==false)
      {
        int n = inci[uu].size();
        if (n>1)
        {
          int ind = rand->rint(n);
          NextEd[uu] = inci[uu][ind];
          NextVer[uu] = neigh[uu][ind];
        }
        else {
          NextEd[uu] = inci[uu][0];
          NextVer[uu] = neigh[uu][0];
        }
        uu = NextVer[uu];
      }
      
      /* Add the loop-erased random walk to the current tree */
        uu = i;
        while (intree[uu]==false)
        {
          intree[uu] = true;
          uu = NextVer[uu];
        }
  }
  
  /* The spanning tree is contained in NextEd */
    for (int j=0; j<NextEd.size(); j++)
    {
      if (NextEd[j]==-1) continue;
      tr.push_back(NextEd[j]);
    } 
  
}

double ComputeLogNTree(const arma::umat & Gr, int n_ver)
{
  if (n_ver<=1) return(0.0);
  double log_n_tree = 0;
  MatrixXd L = MatrixXd::Zero(n_ver,n_ver);
  int ind1 = -1;
  int ind2 = -1;
  
  for (int j=0; j<Gr.n_rows; j++)
  {
    ind1 = Gr.col(0)[j];
    ind2 = Gr.col(1)[j];
    L(ind1,ind2) -= 1;
    L(ind2,ind1) -= 1;
    L(ind1,ind1) += 1;
    L(ind2,ind2) += 1;
  }
  
  VectorXd Dvec(L.topLeftCorner(n_ver-1,n_ver-1).ldlt().vectorD());
  log_n_tree = Dvec.array().log().sum();
  return(log_n_tree);
}

void CreateSubgraphs(const arma::umat & Gr, const std::vector<int> & label,
                     int no_of_clusters, vector<arma::umat>& GB, 
                     std::vector<int>& NVer_B)
{
  vector<vector<int> > v_list1(no_of_clusters);
  vector<vector<int> > v_list2(no_of_clusters);
  vector<vector<int> > e_list(no_of_clusters);
  vector<int> counter(no_of_clusters,0);
  for (int i=0; i<Gr.n_rows ; i++) {
    unsigned int node1 = Gr(i,0);
    unsigned int node2 = Gr(i,1);
    if ( label[node1] == label[node2] )
    {
      e_list[label[node1]].push_back(Gr(i,2));
      v_list1[label[node1]].push_back(node1);
      v_list2[label[node1]].push_back(node2);
      counter[label[node1]]++;
    }
  }
  
  GB.resize(no_of_clusters);
  NVer_B.resize(no_of_clusters,0);
  for (int i=0; i<no_of_clusters ; i++) {
    arma::umat G_b;
    G_b.set_size(counter[i],3);
    map<int, int> v2v;
    int index1 = 0;
    for (int j=0; j<counter[i]; j++)
    {
      if (v2v.find(v_list1[i][j]) == v2v.end()) {
        v2v[v_list1[i][j]] = index1++;
      }
      if (v2v.find(v_list2[i][j]) == v2v.end()) {
        v2v[v_list2[i][j]] = index1++;
      }
      G_b(j,0) = v2v[v_list1[i][j]];
      G_b(j,1) = v2v[v_list2[i][j]];
      G_b(j,2) = e_list[i][j];
    }
    GB[i] = G_b;
    NVer_B[i] = index1;
  }
  
}


// int CreateSubgraph(const arma::umat & Gr, const std::vector<int> & membership,
                      //                    int b, arma::umat & G_b)
// {
  //   std::vector<int> v_list;
  //   for (int j=0; j<membership.size(); j++)
    //       if (membership[j] == b) v_list.push_back(j);
  //   int n_b = v_list.size();
  //   if (n_b == 1)
    //     return (n_b);
  //   
    //   
    //   std::vector<int> e_list;
  //   for (int j=0; j<Gr.n_rows; j++)
    //     if ((membership[Gr.col(0)[j]]==b) && (membership[Gr.col(1)[j]]==b)) e_list.push_back(j);
  //     
    //   if (e_list.size()>0)
      //   {
        //     G_b.set_size(e_list.size(),3);
        //     for (int j=0; j<e_list.size(); j++)
          //     {
            //       int ind1 = -1;
            //       int ind2 = -1;
            //       for (int l=0; l<v_list.size(); l++)
              //       {
                //         if (v_list[l] == Gr.col(0)[e_list[j]]) ind1 = l;
                //         if (v_list[l] == Gr.col(1)[e_list[j]]) ind2 = l;
                //         if ((ind1 != -1)&&(ind2 != -1)) break;
                //       }
            //       G_b.col(0)[j] = ind1;
            //       G_b.col(1)[j] = ind2;
            //       G_b.col(2)[j] = Gr.col(2)[e_list[j]];
            //     }
        //   }
  //   return(n_b);    
  // }

int CreateHyperGraph(const arma::umat& Gr, const std::vector<int> & membership, 
                     int no_of_clusters, arma::umat& G_H)
{
  int NVer_h = no_of_clusters;
  std::vector<int> endpoint_class1(Gr.n_rows);
  std::vector<int> endpoint_class2(Gr.n_rows);
  for (int i=0; i<Gr.n_rows; i++)
  {
    endpoint_class1[i] = membership[Gr.col(0)[i]];
    endpoint_class2[i] = membership[Gr.col(1)[i]];
  }
  std::vector<int> edge_list;
  for (int j=0; j<Gr.n_rows; j++)
    if (endpoint_class1[j]!=endpoint_class2[j]) edge_list.push_back(j);
  
  if (edge_list.size()==0) return(NVer_h);
  if (edge_list.size()>0)
  {
    G_H.set_size(edge_list.size(),3);
    for (int j=0; j<edge_list.size(); j++)
    {
      G_H.col(0)[j] = endpoint_class1[edge_list[j]];
      G_H.col(1)[j] = endpoint_class2[edge_list[j]];
      G_H.col(2)[j] = Gr.col(2)[edge_list[j]];
    }
  }
  return(NVer_h);
  
}


double Sim(const arma::vec & dat1, const arma::vec & dat2, int simtype)
{
  if (simtype == 1)
  {
    double sum1 = 0;
    double sum2 = 0;
    for (int i=0; i<dat1.n_elem; i++)
    {
      sum1 += dat1[i];
      sum2 += dat2[i];
    }
    double sim = 1/(fabs(sum1-sum2)+1);
    return (sim);
  }
  if (simtype == 2)
  {
    double sum1 = 0;
    double sum2 = 0;
    double sum1sq = 0;
    double sum2sq = 0;
    double sumxy = 0;
    for (int i=0; i<dat1.n_elem; i++)
    {
      sum1 += dat1[i];
      sum2 += dat2[i];
      sum1sq += pow(dat1[i],2);
      sum2sq += pow(dat2[i],2);
      sumxy += dat1[i]*dat2[i];
    }
    double cor;
    cor = (sumxy-sum1*sum2/dat1.n_elem)/sqrt((sum1sq-pow(sum1,2)/dat1.n_elem)*(sum2sq-pow(sum2,2)/dat1.n_elem));
    return (cor);
  }
  if (simtype == 3)
  {
    double sumsq = 0;
    for (int i=0; i<dat1.n_elem; i++)
      sumsq += pow((dat1[i]-dat2[i]),2);
      double sim = 1/(sqrt(sumsq)+1);
      return (sim);
  }
}

double ComputeSimilarity(const arma::mat & Dat1, const arma::mat & Dat2, const arma::umat & G, 
                         std::vector<double> & Similarity, int NEd, int NVer, int SimType)
{
  for (int i=0; i<NEd; i++)
  {
    arma::vec dat1;
    arma::vec dat2;
    unsigned int dat_point1;
    unsigned int dat_point2;
    if (G(i,0)>=(NVer/2))
    {
      dat_point1 = (G(i,0)%(NVer/2));
      dat1 = Dat2.col(dat_point1);
    } else {
      dat_point1 = G(i,0);
      dat1 = Dat1.col(dat_point1);
    }
    
    if (G(i,1)>=(NVer/2))
    {
      dat_point2 = (G(i,1)%(NVer/2));
      dat2 = Dat2.col(dat_point2);
    } else {
      dat_point2 = G(i,1);
      dat2 = Dat1.col(dat_point2);
    }
    
    Similarity[i] = Sim(dat1, dat2, SimType);
  }
  
}

int is_connected_component(const arma::umat& Gr, int NVer, vector<int> v_list)
{
    if (v_list.size()==1) return(1);
    
    // Create subgraph containing vertices in v_list
    vector<int> v_list1;
    vector<int> v_list2;
    vector<int> e_list;
    int counter = 0;
    for (int i=0; i<Gr.n_rows ; i++) {
        int node1 = Gr(i,0);
        int node2 = Gr(i,1);
        vector<int>::iterator it1, it2;
        it1 = std::find(v_list.begin(), v_list.end(), node1);
        it2 = std::find(v_list.begin(), v_list.end(), node2);
        if ((it1 != v_list.end())&&(it2 != v_list.end()))
        {
            e_list.push_back(Gr(i,2));
            v_list1.push_back(node1);
            v_list2.push_back(node2);
            counter++;
        }
    }
    
    if (counter == 0) return(0);
    
    arma::umat G_b;
    G_b.set_size(counter,3);
    map<int, int> v2v;
    int index1 = 0;
    for (int j=0; j<counter; j++)
    {
        if (v2v.find(v_list1[j]) == v2v.end()) {
            v2v[v_list1[j]] = index1++;
        }
        if (v2v.find(v_list2[j]) == v2v.end()) {
            v2v[v_list2[j]] = index1++;
        }
        G_b(j,0) = v2v[v_list1[j]];
        G_b(j,1) = v2v[v_list2[j]];
        G_b(j,2) = e_list[j];
    }
    int NVer_b = v_list.size();
    if (index1 < NVer_b) return(0);

    vector<int> EdVal(counter, 1);
    vector<int> membership(NVer_b, -1);
    vector<int> csize;
    int no_of_clusters = FindComponents(G_b, NVer_b, EdVal, membership, csize);
    if (no_of_clusters == 1) return(1);
    else return(0);
}