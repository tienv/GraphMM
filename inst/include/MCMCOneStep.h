#ifndef MCMCONESTEP_H_
#define MCMCONESTEP_H_

#include <iostream>
#include <fstream>
#include <math.h>
#include <numeric>

using namespace std;

#include "Functions.h"
#include "rand.h"
#include "ReadWrite.h"
#include "ComputeLogPrior.h"
#include "MCMCOneStep1.h"

void CreateInitialState(std::vector<int> & TreeOld, std::vector<int> & EdValueOld,
                        std::vector<int> & Label_old, std::vector<int> & BlockSize_old,
                        double & f_old, double & Loglkh_old, double & log_n_tree_old,
                        int & no_of_clusters_old, Rand* randgen,
                        const arma::umat & G, int NVer, int NEd,  double a, int PriorType);
int MCMCOneStep(const arma::umat & G, int NVer, int NEd, double a, 
                Rand* rand_item, int PriorType,
                vector<int> & TreeOld, vector<int> & EdValueOld, vector<int> & Label_old, 
                vector<int> & BlockSize_old, unsigned int & NoChangeCountInChain, 
                double & f_old, double & Loglkh_old, 
                int &  no_of_clusters_old, double & log_n_tree_old);

#endif

