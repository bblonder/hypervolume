#include <vector>
#include <iostream>
#include <stdexcept>
#include <Rcpp.h>
#include "KDTree.h"
//#include <RProgress.h>

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
SEXP kdtree_build_intl(SEXP d, SEXP nr, SEXP nc) // d is the numeric col-major R matrix, nr rows, nc cols
{
  int nrow = as<int>(nr);
  int ncol = as<int>(nc);
  NumericVector data(d);

  if (data.size() != (R_xlen_t) nrow * ncol) {
    throw(length_error("Data not same size as product of nrow and ncol"));
  }

  // Transpose R's col-major buffer into a row-major flat vector<double>
  // (KDTree stores points as points[i*ncol + j]). This is one pass through
  // the buffer; previously this codepath allocated a vector<vector<double>>
  // (one heap allocation per row) via convertMatrixToVector and then the
  // KDTree constructor copied that into its own field -- three copies total
  // for one build. Now: one linear transpose + one std::move.
  const double* src = data.begin();
  vector<double> flat((size_t) nrow * ncol);
  for (int j = 0; j < ncol; j++) {
    for (int i = 0; i < nrow; i++) {
      flat[(size_t) i * ncol + j] = src[(size_t) i + (size_t) j * nrow];
    }
  }

  KDTree* t = new KDTree(std::move(flat), nrow, ncol);
  XPtr<KDTree> p(t, TRUE);
  return p;
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
  
 //RProgress::RProgress pb("[:bar]", nrow);
  if (verbose==1)
  {
    Rcpp::Rcout << "Ball query... \n";
    //pb.tick(0);
  }
  
  for (int i=0; i<nrow; i++)
  {
    
    vector<int> thisIndices;
    vector<double> thisDistances;
    
    vector<double> thisPoint = dataMatrix[i];
    tree->ball_query(thisPoint, radius, thisIndices, thisDistances);
    
    // store the number of points within the ball for each point
    finalCounts.push_back(thisIndices.size());
    
    if (i%10==0 && verbose==1)
    {
      //pb.update(1.0*(i+1)/nrow);
    }
    
  }
  
  if (verbose==1)
  {
    Rcpp::Rcout << "\ndone.\n";
  }
  
  //pb.update(1);
  
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
  
  //RProgress::RProgress pb("[:bar]", nrow);
  if (verbose==1)
  {
    Rcpp::Rcout << "Ball query... \n";
    //pb.tick(0);
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
    
    if (i%10==0 && verbose==1)
    {
      //pb.update(1.0*(i+1)/nrow);
    }
  }
  
  //pb.update(1);
  
  if (verbose==1)
  {
    Rcpp::Rcout << "\ndone.\n";
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
  
  //RProgress::RProgress pb("[:bar]", nrow);
  if (verbose==1)
  {
    Rcpp::Rcout << "Ball query... \n";
    //pb.tick(0);
  }
  
  for (int i=0; i<nrow; i++)
  {
    vector<int> thisIndices;
    
    vector<double> thisPointMin = dataminMatrix[i];
    vector<double> thisPointMax = datamaxMatrix[i];
    
    tree->range_query(thisPointMin, thisPointMax, thisIndices);
    // store the number of points within the ball for each point
    finalCounts.push_back(thisIndices.size());
    
    if (i%10==0 && verbose==1)
    {
      //pb.update(1.0*(i+1)/nrow);
    }
    
    
  }
  
  //pb.update(1);
  
  if (verbose==1)
  {
    Rcpp::Rcout << "\ndone.\n";
  }
  
  return(wrap(finalCounts));
}