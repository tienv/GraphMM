#include <RcppArmadillo.h>
#include"ReadWrite.h"

void ReadData(vector<unsigned int> & seeds, arma::umat & G, std::string datafolder)
{
    ifstream fin;
    ofstream fout;
    std::string filename1, filename2;
    filename1 = datafolder + "/Seed.txt";
    filename2 = datafolder + "/Graph.txt";
    
    fin.open(filename1.c_str());
    if (fin.is_open())
    {
        ReadVec(seeds, fin);
        fin.close();
    } else Rcpp::Rcout << "Error: Cannot open/read file Seed.txt" << endl;
    
    
    fin.open(filename2.c_str());
    if (fin.is_open())
    {
        ReadMat(G,fin);
        fin.close();
    } else Rcpp::Rcout << "Error: Cannot open/read file Graph.txt" << endl;
    
}

