#include "MCMCOneStep.h"

void CreateInitialState(std::vector<int> & TreeOld, std::vector<int> & EdValueOld,
                        std::vector<int> & Label_old, std::vector<int> & BlockSize_old,
                        double & f_old, double & Loglkh_old, double & log_n_tree_old,
                        int & no_of_clusters_old, Rand* randgen,
                        const arma::umat & G, int NVer, int NEd,  double a, int PriorType)
{
    /* Generate spanning tree for graph G using Wilson algorithm */
    
    TreeOld.reserve(NVer-1);
    SampleTree(G, NVer, TreeOld, randgen);
    
     /* Generate edgevalue */
    EdValueOld.resize(NEd, 0);

    /* Derive Label_old */
    Label_old.resize(NVer, -1);
    no_of_clusters_old = FindComponents(G, NVer, EdValueOld, Label_old, BlockSize_old);

    /* Derive f_old */
    f_old = ComputeLogPrior(Label_old, BlockSize_old, a, Loglkh_old, PriorType);

    /* Derive log_n_tree_old */
    std::vector<arma::umat> GB;
    vector<int> NVer_B;
    CreateSubgraphs(G, Label_old, no_of_clusters_old, GB, NVer_B);
    
    log_n_tree_old = 0.0;
    for (int i=0; i<NVer_B.size(); i++)
    {
        if (NVer_B[i] > 1) log_n_tree_old += ComputeLogNTree(GB[i],NVer_B[i]);
    }
    if (no_of_clusters_old >1)
    {
      arma::umat GH;
      int NVer_GH;
      NVer_GH = CreateHyperGraph(G, Label_old, no_of_clusters_old, GH);
      log_n_tree_old += ComputeLogNTree(GH,NVer_GH);
     }
}

int MCMCOneStep(const arma::umat & G, int NVer, int NEd, double a, 
                  Rand* rand_item, int PriorType,
                  vector<int> & TreeOld, vector<int> & EdValueOld, vector<int> & Label_old, 
                  vector<int> & BlockSize_old, unsigned int & NoChangeCountInChain, 
                  double & f_old, double & Loglkh_old, 
                  int &  no_of_clusters_old, double & log_n_tree_old)
{
    int NEdChange = 1;
    double n_tree_change_ratio = 0.8;
    int ans = 1;
        MCMCOneStep1(NEdChange, n_tree_change_ratio, G, NVer, NEd, a, rand_item, PriorType,  
                    TreeOld,  EdValueOld, Label_old, BlockSize_old, 
                    NoChangeCountInChain, f_old,  Loglkh_old, 
                    no_of_clusters_old, log_n_tree_old);
     return(ans); 
}


    

