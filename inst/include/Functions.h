#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <RcppArmadillo.h>
#include <math.h>
#include <time.h>
#include <map>
#include "rand.h"
#include <RcppEigen.h>

using namespace std;
using Eigen::MatrixXd;                  // variable size matrix, double  precision
using Eigen::VectorXd;  // variable size vector, double precision

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

int Sample(int size, arma::vec& prob_raw, double sum_prob, arma::uvec& index, Rand*);
void LabelNormalize (vector<int> & label, vector<int> & csize);
int FindComponents(const arma::umat& G, int NVer, const std::vector<int>& EdValue, 
                   std::vector<int>& membership, std::vector<int>& csize);
void Exchange (std::vector<int>& EdValue1, std::vector<int> & EdValue2, std::vector<int> & EdValue_tmp,
               std::vector<int>& s1, std::vector<int>& s2, std::vector<int>& s_tmp,
               std::vector<int>& tree1, std::vector<int>& tree2, std::vector<int>& tree_tmp,
               double& f1, double& f2, double& f_tmp,
               double& Loglkh1, double& Loglkh2, double& Loglkh_tmp,
               double& log_n_tree1, double& log_n_tree2, double& log_n_tree_tmp,
               int& no_of_clusters1, int& no_of_clusters2, int& no_of_clusters_tmp, int& flag);
void SampleTree(const arma::umat& Gr, int n_ver, std::vector<int> & tr, Rand*);
double ComputeLogNTree(const arma::umat & Gr, int n_ver);
void CreateSubgraphs(const arma::umat & Gr, const std::vector<int> & label,
               int no_of_clusters, vector<arma::umat>& GB, 
               std::vector<int>& NVer_B);
//int CreateSubgraph(const arma::umat & Gr, const std::vector<int> & membership,
//                   int b, arma::umat & G_b);
int CreateHyperGraph(const arma::umat& Gr, const std::vector<int> & membership, 
                     int no_of_clusters, arma::umat& G_H);
double Sim(const arma::vec & dat1, const arma::vec & dat2, int simtype);
double ComputeSimilarity(const arma::mat & Dat1, const arma::mat & Dat2, const arma::umat & G, 
                         std::vector<double> & Similarity, int NEd, int NVer, int SimType);
int is_connected_component(const arma::umat& Gr, int NVer, vector<int> v_list);

#endif
