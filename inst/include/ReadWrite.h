#ifndef READWRITE_H_
#define READWRITE_H_

#include <iostream>
#include<fstream>
#include <vector>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;

template <typename T> void ReadMat(arma::Mat<T> & M, ifstream &fin)
{
    int nrow, ncol;
    fin >> nrow >> ncol;
    M.set_size(nrow,ncol);
    for (int i=0; i<nrow; i++)
      for (int j=0; j<ncol; j++)
	fin >> M(i,j);
};

template <typename T> void WriteMat(const arma::Mat<T> & M, ofstream &fout, bool write_no_of_rows_cols = true)
{
  int nrow, ncol;
  nrow = M.n_rows;
  ncol = M.n_cols;
  fout.precision(16); 
  if (write_no_of_rows_cols) fout << nrow << " " << ncol << "\n";
  for (int i=0; i<nrow; i++)
  {
     for (int j=0; j<ncol; j++) fout << M(i,j) << " ";
     fout << "\n";
  }
};


template <typename T> void ReadVec(vector<T> & V, ifstream & fin)
{
  size_t n;
  fin >> n;
  V.reserve(n);
  for (int i=0; i<n; i++)
  { 
    T x;
    fin >> x;
    V.push_back(x);
  }
}

template <typename T> void WriteVec(const vector<T> & V, ofstream & fout, bool write_no_of_rows = true)
{
  int n;
  n = V.size();
  fout.precision(16);
  if (write_no_of_rows)
    {
      fout << n << "\n";
    }
  
  for (int i=0; i<n; i++)
  {
    fout << V[i] << "\n";
  }
}

template <typename T> void PrintVec(const vector<T> & V, const char *s)
{
  int n;
  n = V.size();
  cout << s << "\n";
  for (int i=0; i<n; i++)
  {
    cout << V[i] << " ";
  }
  cout << "\n";
}

template <typename T> void PrintVec1(const arma::Col<T> & V, const char *s)
{
  int n;
  n = V.size();
  cout << s << "\n";
  for (int i=0; i<n; i++)
  {
    cout << V[i] << " ";
  }
  cout << "\n";
}
 
 
void ReadData(vector<unsigned int> & seeds, arma::umat & G, std::string datafolder);

#endif
