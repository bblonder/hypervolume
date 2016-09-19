#include <vector>
#include <iostream>
#include <stdexcept>
#include <Rcpp.h>
#include "KDTree.h"

using namespace std;

vector<vector< double > > convertMatrixToVector(double array[], int nrow, int ncol)
{
    vector<vector< double > > v(nrow, vector< double >(ncol));
    
    for(int i = 0; i < nrow*ncol; ++i) {
        v[i / ncol][i % ncol] = array[i];
    }
    
    return v;
}

using namespace Rcpp;

// [[Rcpp::export]]
SEXP kdtree_build_intl(SEXP d, SEXP nr, SEXP nc) // d is the numeric data, nr is the rows, nc is the cols
{
  int nrow = as<int>(nr);
  int ncol = as<int>(nc);
  NumericVector data(d);
  
  if(data.size() != nrow*ncol)
  {
    throw(length_error("Data not same size as product of nrow and ncol"));
  }
  
  vector<vector< double > > dataMatrix
        = convertMatrixToVector(data.begin(), nrow, ncol);

  KDTree *t = new KDTree(dataMatrix);
  
  XPtr<KDTree> p(t, TRUE);
  
  return(p);
}

// [[Rcpp::export]]
SEXP kdtree_ball_query_multiple(SEXP tr, SEXP ptlist, SEXP nr, SEXP nc, SEXP r, SEXP verb)
{
  XPtr<KDTree> tree = as<XPtr<KDTree> >(tr);
  int nrow = as<int>(nr);
  int ncol = as<int>(nc);
  NumericVector data(ptlist);
  double radius = as<double>(r);
  bool verbose = as<int>(verb);
  
  vector<vector< double > > dataMatrix
    = convertMatrixToVector(data.begin(), nrow, ncol);
    
  vector<int> finalCounts;
    
  if (ncol != tree->ndims())
  {
    throw(length_error("Points not same dimensionality as data in kdtree"));  
  }
  
  for (int i=0; i<nrow; i++)
  {
    vector<int> thisIndices;
    vector<double> thisDistances;
    
    vector<double> thisPoint = dataMatrix[i];
    tree->ball_query(thisPoint, radius, thisIndices, thisDistances);
    
    // store the number of points within the ball for each point
    finalCounts.push_back(thisIndices.size());
    
    if (verbose==1 && (i%10000)==1)
    {
      Rcpp::Rcout << " " << 1.0*i/nrow << " ";
      if ((i%50000)==1)
      {
        Rcpp::Rcout << "\n";
      }
    }
  }
  
  return(wrap(finalCounts));
}

// [[Rcpp::export]]
SEXP kdtree_ball_query_id_multiple(SEXP tr, SEXP ptlist, SEXP nr, SEXP nc, SEXP r, SEXP verb)
{
  XPtr<KDTree> tree = as<XPtr<KDTree> >(tr);
  int nrow = as<int>(nr);
  int ncol = as<int>(nc);
  NumericVector data(ptlist);
  double radius = as<double>(r);
  bool verbose = as<int>(verb);
  
  vector<vector< double > > dataMatrix
    = convertMatrixToVector(data.begin(), nrow, ncol);
  
  vector< vector< int > > finalIDs;
  
  if (ncol != tree->ndims())
  {
    throw(length_error("Points not same dimensionality as data in kdtree"));  
  }
  
  for (int i=0; i<nrow; i++)
  {
    vector<int> thisIndices;
    
    vector<double> thisDistances;
    
    vector<double> thisPoint = dataMatrix[i];
    tree->ball_query(thisPoint, radius, thisIndices, thisDistances);
    
    // store the number of points within the ball for each point
    if (thisIndices.size() > 0)
    {
      finalIDs.push_back(thisIndices);
    }
    else
    {
      vector<int> empty;
      empty.push_back(-1);
      finalIDs.push_back(empty);
    }
    
    if (verbose==1 && (i%10000)==1)
    {
      Rcpp::Rcout << " " << 1.0*i/nrow << " ";
      if ((i%50000)==1)
      {
        Rcpp::Rcout << "\n";
      }
    }
  }
  
  return(wrap(finalIDs));
}

// [[Rcpp::export]]
SEXP kdtree_range_query_multiple(SEXP tr, SEXP pminlist, SEXP pmaxlist, SEXP nr, SEXP nc, SEXP verb)
{
  XPtr<KDTree> tree = as<XPtr<KDTree> >(tr);
  int nrow = as<int>(nr);
  int ncol = as<int>(nc);
  NumericVector datamin(pminlist);
  NumericVector datamax(pmaxlist);
  bool verbose = as<int>(verb);

  if (ncol != tree->ndims())
  {
    throw(length_error("pmin or pmax not same dimensionality as data in kdtree"));  
  }
  
  vector<vector< double > > dataminMatrix
    = convertMatrixToVector(datamin.begin(), nrow, ncol);
  vector<vector< double > > datamaxMatrix
    = convertMatrixToVector(datamax.begin(), nrow, ncol);
    
  vector<int> finalCounts;
    
  if (ncol != tree->ndims())
  {
    throw(length_error("Points not same dimensionality as data in kdtree"));  
  }
  
  for (int i=0; i<nrow; i++)
  {
    vector<int> thisIndices;
    
    vector<double> thisPointMin = dataminMatrix[i];
    vector<double> thisPointMax = datamaxMatrix[i];

    tree->range_query(thisPointMin, thisPointMax, thisIndices);
    // store the number of points within the ball for each point
    finalCounts.push_back(thisIndices.size());
    
    if (verbose==1 && (i%10000)==1)
    {
      Rcpp::Rcout << " " << 1.0*i/nrow << " ";
      if ((i%50000)==1)
      {
        Rcpp::Rcout << "\n";
      }
    }
  }
  
  return(wrap(finalCounts));
}