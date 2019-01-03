#include <ComputeLogPrior.h>

double ComputeLogPrior(const std::vector<int>& Label, const vector<int> & BlockSize,
                      double a, double& Loglkh, int PriorType)
{
    int no_of_clusters = *max_element(Label.begin(), Label.end())+1;
    int NVer = Label.size();
    Loglkh = 0.0;
    double f = Loglkh;
    if (PriorType == 1)
    {
      for (int k=0; k<no_of_clusters; k++)
      {
        f += lgamma(BlockSize[k]);
      }
      f = f + no_of_clusters*log(a);
    }
    return(f);
}

