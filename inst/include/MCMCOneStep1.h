#ifndef MCMCONESTEP1_H_
#define MCMCONESTEP1_H_

#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <numeric>
#include <RcppArmadillo.h>
#include <RcppEigen.h>


using namespace std;
using Eigen::MatrixXd;                  // variable size matrix, double  precision
using Eigen::VectorXd;  // variable size vector, double precision

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]


#include "Functions.h"
#include "rand.h"
#include "ReadWrite.h"
#include "ComputeLogPrior.h"

void  MCMCOneStep1(int NEdChange, double n_tree_change_ratio, const arma::umat & G, int NVer, int NEd, 
                   double a, Rand* rand_item, int PriorType,
                   vector<int> & TreeOld, vector<int> & EdValueOld, 
                   vector<int> & Label_old, vector<int> & BlockSize_old, 
                   unsigned int & NoChangeCountInChain, 
                   double & f_old, double & Loglkh_old, 
                   int &  no_of_clusters_old, double & log_n_tree_old);

#endif