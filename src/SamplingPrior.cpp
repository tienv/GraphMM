#include <RcppArmadillo.h>

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <numeric>

#include "rand.h"
#include "ReadWrite.h"
#include "Functions.h"
#include "ComputeLogPrior.h"
#include "MCMCOneStep.h"

using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

std::vector<int> GetCsize (std::vector<int> membership)
{
  int no_of_clusters = *std::max_element(membership.begin(),membership.end());
  no_of_clusters ++;
  std::vector<int> csize(no_of_clusters, 0);
  for (int i=0; i<membership.size(); i++)
  {
    csize[membership[i]]++;
  }
  return(csize);
}
double GetTheoPrior (std::vector<int> membership, double a, int PriorType)
{
  std::vector<int> csize;
  csize = GetCsize(membership);
  double prior;
  double Loglkh = 0;
  prior = ComputeLogPrior(membership, csize, a, Loglkh, PriorType);
  return(prior);
}

Rcpp::List CountUniqueStates(int burn, int NSave, int NVer, int iter, 
                      double a, int PriorType, arma::Mat<int> states)
{
  
  // Read all states
  ifstream fin;
  ofstream fout;
  std::vector<std::vector<int> > AS(iter);
  for (int i=0; i<iter; i++)
  {
      AS[i].reserve(NVer);
      for (int j=0; j<NVer; j++)
      {
        AS[i].push_back(states(i,j));
      }
  }

  std::sort(AS.begin(), AS.end());
  std::vector<int> freq;
  freq.reserve(iter);
  std::vector<double> theo_prior;
  theo_prior.reserve(iter);
  std::vector<int> z = AS[0];
  freq.push_back(1);
  theo_prior.push_back(GetTheoPrior(z, a, PriorType));
  
  std::vector<int> z1;
  for (int i=1; i<AS.size(); i++)
  {
    z1 = AS[i];
    if (z1 != z)
    {
      freq.push_back(0);
      theo_prior.push_back(GetTheoPrior(z1, a, PriorType));
      z = z1;
    }
    ++freq.back();
  }
  AS.erase(std::unique(AS.begin(), AS.end()), AS.end());
  arma::Mat<int> UniqueStates(AS.size(), NVer);
  for (int i=0; i<AS.size(); i++)
    for (int j=0; j<NVer; j++)
      UniqueStates(i,j) = AS[i][j];
  
  return Rcpp::List::create(Rcpp::Named("UniqueStates") = Rcpp::as<Rcpp::IntegerMatrix>(Rcpp::wrap(UniqueStates)),
                            Rcpp::Named("Frequency") = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(freq))); 
  
} 


Rcpp::List MCMCOneChain (int iter, int NSave, arma::umat G, int NVer, 
                  int NEd, int burn, double a, int PriorType,
                  vector<int> & TreeOld, vector<int> & EdValueOld,
                  vector<int> & Label_old, vector<int> & BlockSize_old,
                  double & f_old, double & Loglkh_old, double & log_n_tree_old, 
                  int & no_of_clusters_old, Rand * rand_item)
{
 
  unsigned int  NoChangeCountInChain = 0;
  int CurrentIter = 0;
  
  bool IsBurning;
  int t = CurrentIter;
  vector<double> LogPostProb;
  LogPostProb.reserve(iter);
  /* vector<double> LogLikelihood; */
  /* LogLikelihood.reserve(n_iter_checkpoint); */
  arma::Mat<int> states(iter, NVer);
  states.fill(-1);
  while (t<iter)
  {
    if (t>=burn) IsBurning = false;
    else IsBurning = true;
    int success = 1;
    success = MCMCOneStep(G, NVer, NEd, a, rand_item, PriorType,  
                TreeOld,  EdValueOld, Label_old, BlockSize_old, 
                NoChangeCountInChain, f_old,  Loglkh_old, 
                no_of_clusters_old, log_n_tree_old);
    
    if (success)
    {
        LogPostProb.push_back(f_old);
        for (int j=0; j<Label_old.size(); j++)
        {
          states(t,j) = Label_old[j];
        }
        t++;  
    }
  }
  delete rand_item;

  Rcpp::List RL = CountUniqueStates(burn, NSave, NVer, iter, a, PriorType, states);

  return Rcpp::List::create(Rcpp::Named("AllStates") = Rcpp::as<Rcpp::IntegerMatrix>(Rcpp::wrap(states)),
                            Rcpp::Named("LogPostProb") =  Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(LogPostProb)),
                            Rcpp::Named("UniqueStates") = RL["UniqueStates"],
                            Rcpp::Named("Frequency") = RL["Frequency"]); 
  
}

// [[Rcpp::export]]    
Rcpp::List SamplingPrior(int iter, int NVer, int NEd, int PriorType, double a,
                         std::string datafolder)
{

  int NSave = iter;
  int burn = 0;
  /* Read data */
  vector<unsigned int> seeds;
  arma::umat G;
  ReadData(seeds, G, datafolder);

  // initialize random number generators
  Rand * randgen = new Rand(seeds);

  /* Initial state */
  std::vector<int> TreeOld;
  std::vector<int> EdValueOld;
  double log_n_tree_old = 0.0;
  int  no_of_clusters_old = 0;
  std::vector<int> Label_old;
  std::vector<int> BlockSize_old;
  double f_old = 0.0;
  double Loglkh_old = 0.0;
  

  CreateInitialState(TreeOld, EdValueOld, Label_old,
                     BlockSize_old, f_old, Loglkh_old, log_n_tree_old, 
                     no_of_clusters_old, randgen, G, NVer, NEd, a, PriorType);

  /* Run MCMC chain */
  Rcpp::List L = MCMCOneChain(iter, NSave, G, NVer, NEd, burn, a, PriorType,
               TreeOld, EdValueOld, Label_old, BlockSize_old, f_old, Loglkh_old, 
               log_n_tree_old, no_of_clusters_old, randgen);
  return(L);
}  
