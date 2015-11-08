// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 
#include <iostream> 
#include <Rmath.h> 

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

// fct 1 : inline double medianIV_hpp(IntegerVector data);
// fct 2 : inline double medianNV_hpp(NumericVector data);
// fct 3 : inline double medianVecInteger_hpp(std::vector < int > data);
// fct 4 : inline double medianVecDouble_hpp(std::vector < double > data);
// fct 5 : inline IntegerVector rank_hpp(NumericVector x);
inline IntegerVector order_hpp(const NumericVector& x); // fct 6 :

#ifndef __UTILITIES__
#define __UTILITIES__

const double CST_EPSILON = std::numeric_limits < double>::epsilon();
const double CST_PROGRESS = 51;

// fct 1 ////////////////////////////////////////////////////////////
inline double medianIV_hpp(IntegerVector data){
    
    // taille du vecteur
    int taille = data.size() ;
    if (taille == 0) throw "Empty";
    
    // rangement du vecteur
    std::sort(data.begin(), data.end()) ;
    
    // calcul de la mediane
    double mediane ;
    if (taille  % 2 == 0)
    {   mediane = (data[ taille / 2 - 1 ] + data[taille / 2]) / 2.0; }
    else 
    {   mediane = data[taille / 2] ;  }
    
    return mediane;
  }


// fct 2 ////////////////////////////////////////////////////////////
inline double medianNV_hpp(NumericVector data){
    
    // taille du vecteur
    int taille = data.size() ;
    if (taille == 0) throw "Empty";
    
    // rangement du vecteur
    std::sort(data.begin(), data.end()) ;
    
    // calcul de la mediane
    double mediane ;
    if (taille  % 2 == 0)
    {   mediane = (data[ taille / 2 - 1 ] + data[taille / 2]) / 2.0; }
    else 
    {   mediane = data[taille / 2] ;  }
    
    return mediane;
  }


// fct 3 ////////////////////////////////////////////////////////////
inline double medianVecInteger_hpp(std::vector < int > data){
    
    // taille du vecteur
    size_t taille = data.size() ;
    if (taille == 0) throw "Empty";
    
    // rangement du vecteur
    std::sort(data.begin(), data.end()) ;
    
    // calcul de la mediane
    double mediane ;
    if (taille  % 2 == 0)
    {   mediane = (data[ taille / 2 - 1 ] + data[taille / 2]) / 2.0; }
    else 
    {   mediane = data[taille / 2] ;  }
    
    return mediane;
  }   

// fct 4 ////////////////////////////////////////////////////////////
inline double medianVecDouble_hpp(std::vector < double > data){
    
    // taille du vecteur
    size_t taille = data.size() ;
    if (taille == 0) throw "Empty";
    
    // rangement du vecteur
    std::sort(data.begin(), data.end()) ;
    
    // calcul de la mediane
    double mediane ;
    if (taille  % 2 == 0)
    {   mediane = (data[ taille / 2 - 1 ] + data[taille / 2]) / 2.0; }
    else 
    {   mediane = data[taille / 2] ;  }
    
    return mediane;
  }   

// fct 5 ////////////////////////////////////////////////////////////
inline IntegerVector rank_hpp(NumericVector x) {
  
  NumericVector rank;
  rank = order_hpp(x);
  return(order_hpp(rank));
  
}

// fct 6 ////////////////////////////////////////////////////////////
inline IntegerVector order_hpp(const NumericVector& x) {
  NumericVector sorted = clone(x).sort();
  return match(sorted, x);
}

#endif //__UTILITIES__
