#ifndef COMPUTELOGPRIOR_H_
#define COMPUTELOGPRIOR_H_


#include <math.h>
#include <float.h>
#include <numeric>  //to use std::accummulate
#include <Functions.h>

#include "ReadWrite.h"


using Eigen::MatrixXd;                  // variable size matrix, double  precision
using Eigen::VectorXd;  // variable size vector, double precision

using namespace std;

using Eigen::MatrixXd;                  // variable size matrix, double  precision
using Eigen::VectorXd;  // variable size vector, double precision
using namespace std;

double ComputeLogPrior(const std::vector<int>& Label, const vector<int> & BlockSize,
                    double a, double& Loglkh, int PriorType);
#endif
