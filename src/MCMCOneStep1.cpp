#include "MCMCOneStep1.h"

void  MCMCOneStep1(int NEdChange, double n_tree_change_ratio, const arma::umat & G, int NVer, int NEd, 
                   double a, Rand* rand_item, int PriorType,
                   vector<int> & TreeOld, vector<int> & EdValueOld, 
                  vector<int> & Label_old, vector<int> & BlockSize_old, 
                  unsigned int & NoChangeCountInChain, 
                  double & f_old, double & Loglkh_old, 
                  int &  no_of_clusters_old, double & log_n_tree_old)
{
    vector<int> EdValueNew;
    EdValueNew = EdValueOld;
    double proposal_new2old;
    double proposal_old2new;
    double LMHratio;
    vector<int> Label_new(Label_old), BlockSize_new;
    double f_new, Loglkh_new;
    int no_of_clusters_new;
    
    arma::vec WithinProb_old2new(TreeOld.size());
    WithinProb_old2new.fill(1.0);
    double WithinSum_old2new = double(TreeOld.size());
    
    int EdChange;
    arma::uvec WithinEdChange_Id(NEdChange);
    Sample(NEdChange, WithinProb_old2new, WithinSum_old2new, WithinEdChange_Id, rand_item);     
    for (int i = 0; i<WithinEdChange_Id.n_elem; i++)
    {
        EdChange = TreeOld[WithinEdChange_Id[i]];
        EdValueNew[EdChange] = 1-EdValueOld[EdChange];
     }    
    
    proposal_new2old = log(1.0/double(TreeOld.size()));
    proposal_old2new = log(1.0/double(TreeOld.size()));;  
    

    LMHratio = proposal_new2old - proposal_old2new;

    no_of_clusters_new = FindComponents(G, NVer, EdValueNew, Label_new, BlockSize_new);
    f_new = ComputeLogPrior(Label_new, BlockSize_new, a, Loglkh_new, PriorType);
    LMHratio += f_new - f_old;

    /* Calculating number of spanning tree given partition using parallel */
    std::vector<arma::umat> GB;
    vector<int> NVer_B;
    CreateSubgraphs(G, Label_new, no_of_clusters_new, GB, NVer_B);
    
    vector<double> log_n_tree_new1(GB.size(),0.0);
    // #pragma omp parallel for num_threads(4)
        for (int i=0; i<GB.size(); i++)
        {
            log_n_tree_new1[i] = ComputeLogNTree(GB[i], NVer_B[i]);
        }
    double log_n_tree_new = std::accumulate(log_n_tree_new1.begin(),log_n_tree_new1.end(),0.0);
    
    if (no_of_clusters_new >1)
    {
        arma::umat GH;
        int NVer_GH;
        NVer_GH = CreateHyperGraph(G, Label_new, no_of_clusters_new, GH);
        log_n_tree_new += ComputeLogNTree(GH,NVer_GH);
    }
    
    /* Calculating MH acceptance probability */
    LMHratio += log_n_tree_old - log_n_tree_new;
    double accept = std::min(0.0, LMHratio);
    double logU = log(rand_item->runif());
    bool IsAccept;
    if (logU <= accept) {
        IsAccept = 1; 
    } else IsAccept = 0; 
    if (logU<=accept) 
    {
        no_of_clusters_old = no_of_clusters_new;
        EdValueOld = EdValueNew;
        Label_old = Label_new;
        BlockSize_old = BlockSize_new;
        f_old = f_new;
        Loglkh_old = Loglkh_new;
        log_n_tree_old = log_n_tree_new;
    } else NoChangeCountInChain++;
    
        /* Update tree, fixed partition */
        /* Update tree*/
        std::vector<int> tree;
        tree.reserve(NVer-1);
        
        /* Sample the trees to be changed */
        int n_b_change = std::ceil(no_of_clusters_old*n_tree_change_ratio);
         arma::uvec b_change(n_b_change);
        arma::vec prob(no_of_clusters_old);
        for (int i=0; i<no_of_clusters_old; i++) prob[i] = 1;
        Sample(n_b_change, prob, no_of_clusters_old, b_change, rand_item);

        /* Adding trees that don't change */
        for (int i=0; i<TreeOld.size(); i++) 
        {
            int z1 = Label_old[G.col(0)[TreeOld[i]]];
            int z2 = Label_old[G.col(1)[TreeOld[i]]];
            if (z1 != z2) continue;
            int flag1 = 0;
            for (int j=0; j<b_change.size(); j++)
                if (b_change[j] == z1)
                {
                    flag1 = 1;
                    break;
                }
                if (flag1 == 0)
                {
                    tree.push_back(TreeOld[i]);
                }
        }
        
        /* Tree to be changed */
        
        {
            std::vector<std::vector<int> > tr(b_change.size());
            std::vector<arma::umat> GB1;
            vector<int> NVer_B1;
            CreateSubgraphs(G, Label_old, no_of_clusters_old, GB1, NVer_B1);
            for (int k=0; k<b_change.size(); k++)
            {
                /* Create subgraph */
                arma::umat Gb;
                int NVer_Gb;
                Gb = GB1[b_change[k]];
                NVer_Gb = NVer_B1[b_change[k]];
                
                
                
                /* Find spanning tree for subgraph */
                if (NVer_Gb <= 1) continue;

                vector<int> tr_sub;
                tr_sub.reserve(NVer_Gb-1);
                
                SampleTree(Gb, NVer_Gb, tr_sub, rand_item);
                
                vector<int> tr_sub_orgind(NVer_Gb-1,-1);
                for (int j=0; j<tr_sub.size(); j++)
                {
                    tr_sub_orgind[j] = Gb.col(2)(tr_sub[j]);
                }
                tr[k] = tr_sub_orgind;
            }
            
            for (int k=0; k<b_change.size(); k++)
            {
                if (tr[k].size()==0) continue; 
                for (int j=0; j<tr[k].size(); j++)
                {
                    if (tr[k][j]==-1) continue;
                    tree.push_back(tr[k][j]);
                }
            }
        }
        
        /*Spanning tree for Hyper-graph */
        if (no_of_clusters_old > 1)
        {
            arma::umat GH;
            int NVer_GH;
            NVer_GH = CreateHyperGraph(G, Label_old, no_of_clusters_old, GH);
            
            /* Finding spanning tree for Hyper-graph */
            std::vector<int> tr_H;
            tr_H.reserve(NVer_GH-1);
            SampleTree(GH, NVer_GH, tr_H, rand_item);
            
            /* Add this to the overall tree */
            if (tr_H.size() > 0)
            {
                for (int j=0; j<tr_H.size(); j++)
                {
                    tree.push_back(GH.col(2)(tr_H[j]));
                }
            }
        }
        TreeOld = tree;
        
        /* Update binary edges so that the partition is unchanged */
        std::fill(EdValueNew.begin(),EdValueNew.end(),0);
        for (int i=0; i<TreeOld.size(); i++)
            if (Label_old[G.col(0)[TreeOld[i]]] == Label_old[G.col(1)[TreeOld[i]]]) EdValueNew[TreeOld[i]] = 1;    
            EdValueOld = EdValueNew;
            
}
    