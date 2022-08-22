#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix dist_wtbinary(Rcpp::NumericMatrix ww) {
  
  //ww is the weight matrix, after getWeights
  int ncols = ww.ncol();
  int nrows = ww.nrow();
  
  //Rcout << "No of rows " << nrows; 
  //Rcout << "No of cols " << ncols;
  
  //create variables
  double dist_value;
  double w11; double w10; double w01; double fract; double fracb;  
  
  Rcpp::NumericMatrix dist_matrix(nrows, nrows);
  // http://stackoverflow.com/questions/23748572/initializing-a-matrix-to-na-in-rcpp
  std::fill( dist_matrix.begin(), dist_matrix.end(), Rcpp::NumericVector::get_na() );
  
  for (int i=0; i<nrows; i++){
    for (int j=0; j<nrows; j++){
      
      //initializa variables
      dist_value = w11 = w10 = w01 = 0.0; 
      
      //we fill matrix for pairs - i,j and j,i - NA means that pair is empty, otherwise skip
      if(Rcpp::NumericVector::is_na(dist_matrix(i,j))){
        
        NumericVector v1= ww.row(i);
        NumericVector v2= ww.row(j);
        //a third loop to go ove rall pairs of 01, 10 and 11. Can improve this? 
          for (int k=0; k<ncols; k++){
            //Rcout << "k= " << k;
            //Rcout << "\n"; 
            if(  (v1[k] != 0. && v2[k] != 0.) ){ w11 += (v1[k] + v2[k]);}
            //Rcout << "w11= " << w11; 
            if(  (v1[k] != 0. && v2[k] == 0.) ){ w10 += (v1[k] + v2[k]);}
            //Rcout << "w10=" << w10; }
        if(  (v1[k] == 0. && v2[k] != 0.) ){ w01 += (v1[k] + v2[k]);}
        //Rcout << "w01= " << w01; }
    }
    
    fract = w01 + w10; 
    fracb = w01 + w10 + w11; 
    
    if(fract ==0 && fracb == 0){dist_value = 0.0;}
    //computing weighted distance 
    else{dist_value = (w01 + w10)/(w01 + w10 + w11);}
    dist_matrix(i,j) = dist_value;
    dist_matrix(j,i) = dist_value;
  }
}
}
return dist_matrix;    
}




