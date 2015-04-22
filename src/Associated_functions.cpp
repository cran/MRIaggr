#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

// List calcHemi_cpp(const arma::mat& px_hemiL, const arma::mat& px_hemiR, double sd_data, int p, bool symetrie)
// List calcContro_cpp(const arma::mat& contrast, const arma::mat& coords_px, const IntegerVector& index_k, const IntegerVector& index_k_contro,
//                    double d_lim, double lambda, int param_ref, double var_ref ,
//                    bool type_moy, bool type_med, bool type_NN, bool trace)
// arma::mat filtrage2D_cpp(const arma::mat& M_data, const arma::mat& M_operateur, const arma::mat& index_data, bool w_contrast, bool na_rm){
// arma::mat filtrage2Dmed_cpp(const arma::mat& M_data, const arma::mat& M_operateur, const arma::mat& index_data, bool na_rm)
// List filtrage3D_cpp(const NumericVector& Vec_data, const IntegerVector& p_data,
//                    const NumericVector& Vec_operateur, const IntegerVector& p_operateur, 
//                    const arma::mat& index_data, bool w_contrast, bool na_rm)
// arma::mat filtrage3Dmed_cpp(const arma::mat& M_data, const arma::mat& M_operateur, const arma::mat& index_data, bool na_rm)
double medianIV_cpp(IntegerVector data);
double medianNV_cpp(NumericVector data);
double medianVecInteger_cpp(std::vector<int> data);
double medianVecDouble_cpp(std::vector<double> data);
// List calcGroupsCoords_cpp(const arma::mat& coords_NNA, const IntegerVector& index_NNA, const arma::mat& Neighborhood, IntegerVector coords_max,
//                           int max_groups, bool trace)
// List calcGroupsW_cpp(const S4& W, const vector<int>& subset, int max_groups)
 List subsetW_cpp(const IntegerVector& W_i, const IntegerVector& W_p, IntegerVector subset);
// List calcRadius_cpp(const arma::mat& coords, const NumericVector& sample, double threshold, const LogicalVector& subset_bary, bool trace)

//  1 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List calcHemi_cpp(const IntegerVector& coordsI, const IntegerVector& coordsJ, List ls_indexK, int n_num, const NumericVector& value, int n,
                    double i_pos, double j_pos, double angle_pos,
                    double sd_data, double p, bool symetrie){ 
  
  // critere de symetrie
  double valeur_vois, dist_vois, asymetrie_valeur=0, poids_vois;
  int asymetrie_nb=0;
  double distI, distJ ;
  double sdp_data = pow(sd_data,p);
  
  // parametrage du nouveau de repere
  int n_indexK, n_hemiL, n_hemiR, iter_vxAll ;
  double cos_alpha = cos(angle_pos);
  double sin_alpha = sin(angle_pos);
  
  vector<int> indexK;
  vector<double> i_nouveauL, i_nouveauR, j_nouveauL, j_nouveauR, valueL, valueR;
  double i_tempo;
  
  for(int iter_k=0 ; iter_k<n_num ; iter_k++){
    
    indexK = ls_indexK[iter_k];
    n_indexK = indexK.size();
    n_hemiL = 0;
    n_hemiR = 0;
    
    i_nouveauL.resize(0);
    i_nouveauR.resize(0);
    j_nouveauL.resize(0);
    j_nouveauR.resize(0);
    valueL.resize(0);
    valueR.resize(0);
    
    for(int iter_vx=0 ; iter_vx<n_indexK ; iter_vx++){
      iter_vxAll = indexK[iter_vx];  
      
      // nouvelles coordonnees
      i_tempo = cos_alpha*(coordsI[iter_vxAll]-i_pos) - sin_alpha*(coordsJ[iter_vxAll]-j_pos);
      
      if(i_tempo>0){ // left hemisphere (radiological convention)
      
      i_nouveauL.push_back(i_tempo);
      j_nouveauL.push_back(sin_alpha*(coordsI[iter_vxAll]-i_pos) + cos_alpha*(coordsJ[iter_vxAll]-j_pos));
      valueL.push_back(value[iter_vxAll]);
      n_hemiL++;
      
      }else{ // right hemisphere
      
      i_nouveauR.push_back(i_tempo);
      j_nouveauR.push_back(sin_alpha*(coordsI[iter_vxAll]-i_pos) + cos_alpha*(coordsJ[iter_vxAll]-j_pos));
      valueR.push_back(value[iter_vxAll]);
      n_hemiR++;
      }
      
    }
    
    // evalutation du critere de symetrie vis a vis de l hemisphere gauche
    for(int iter_pxL=0 ; iter_pxL<n_hemiL ; iter_pxL++){
      poids_vois = 0;
      valeur_vois = 0;        
      
      // iteration sur les voisins : moyenne ponderee par les distance des valeurs voisines
      for(int iter_pxR=0 ; iter_pxR<n_hemiR ; iter_pxR++){
        
        distI = abs(i_nouveauL[iter_pxL]+i_nouveauR[iter_pxR]) ;
        distJ = abs(j_nouveauL[iter_pxL]-j_nouveauR[iter_pxR]) ;
        
        if(distI < 1 && distJ < 1){
          dist_vois = sqrt(pow(distI,2) + pow(distJ,2));
          valeur_vois = valeur_vois + valueR[iter_pxR] * (2-dist_vois)         ;           
          poids_vois = poids_vois + (2-dist_vois) ;
        }
        
      }
      
      if(poids_vois>0){
        
        valeur_vois = valeur_vois/poids_vois;
        if(symetrie){
          asymetrie_valeur =  asymetrie_valeur + 1/(1+pow(abs(valueL[iter_pxL]-valeur_vois),p)/sdp_data);
        }else{        
          asymetrie_valeur = asymetrie_valeur + pow(abs(valueL[iter_pxL]-valeur_vois),p)/sdp_data   ;
        }
         
        asymetrie_nb = asymetrie_nb + 1;
      }
    }
  }
  
  // export
  return(List::create(Named("asymetrie_valeur")  = asymetrie_valeur,
                      Named("asymetrie_moy")  = asymetrie_valeur/asymetrie_nb,
                      Named("asymetrie_nb")  = asymetrie_nb));
  
}


//  2 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List calcContro_cpp(const arma::mat& contrast, const arma::mat& coords_px, const IntegerVector& index_k, const IntegerVector& index_k_contro,
                    double d_lim, double lambda, int param_ref, double var_ref ,
                    bool type_moy, bool type_med, bool type_NN, bool trace){
                
                // preparation
                int n_k = index_k.size();
                int n_k_contro = index_k_contro.size();
                int D = contrast.n_cols;
                double dist1, dist2 ;
                bool voisin_proche ;
                int index_px, index_contro, nb_vois;       
                arma::mat valeur_NormContro(n_k,D);
                std::fill(valeur_NormContro.begin(),valeur_NormContro.end(),NA_REAL);
                vector<bool> index_plot_k(n_k,false),index_plot_k_contro(n_k_contro,false); 
                vector<int> index_vois;
                vector<int> iter_vois;
                vector<double> valeurs_miroir;
                
                double poids_vois, norm , valeur_moyenne ; 
                int index_ref_contro; 
                double diff_ref, diff_tempo, dist_vois;     
                  
                // iteration sur les pixels
                for(int iter_px=0 ; iter_px < n_k ; iter_px++){
                  index_px = index_k[iter_px];
                  index_vois.resize(0);
                  iter_vois.resize(0);
                  voisin_proche = false;
                                
                  for(int iter_contro=0 ; iter_contro < n_k_contro ; iter_contro++ ){
                  index_contro = index_k_contro[iter_contro];
                  // evaluation des distances
                  dist1 = abs(coords_px(index_contro,0)+coords_px(index_px,0)) ;
                  dist2 = abs(coords_px(index_contro,1)-coords_px(index_px,1)) ;
                  
                 
                  // si le point contro est assez proche
                  if( dist1 <= d_lim && dist2 <= d_lim )
                  {   
                      if(sqrt(pow(dist1,2) + pow(dist2,2)) <= 1)
                     {voisin_proche = true ;}
                  
                    // maj des vois utilises pour la normalisation
                    index_vois.push_back(index_contro)  ;             
                    if(trace){iter_vois.push_back(iter_contro);} 
                  }
                  
                 }
                // si l on trouve des voisins dont au moins un proche
                nb_vois = index_vois.size(); 
                if(voisin_proche && nb_vois>0){
                
                // maj des px utilises
                if(trace){
                index_plot_k[iter_px] = true;
                for(int iter=0 ; iter < nb_vois ; iter++){
                index_plot_k_contro[iter_vois[iter]] = true;
                }
                }
                
                if(type_moy){
                                     
                  for(int iter_d=0 ; iter_d<D ; iter_d++){
                    norm=0;
                    valeur_moyenne=0;

                    for(int iter_vois=0 ;  iter_vois<nb_vois ; iter_vois++){
                          
                          index_contro = index_vois[iter_vois];
                          poids_vois= 2*d_lim-sqrt(pow(coords_px(index_contro,0)+coords_px(index_px,0),2) + pow(coords_px(index_contro,1)-coords_px(index_px,1),2));

                          valeur_moyenne = valeur_moyenne + poids_vois*contrast(index_contro,iter_d);
                          norm = norm + poids_vois;
                    }

                     valeur_NormContro(iter_px,iter_d) = contrast(index_px,iter_d)-valeur_moyenne/norm;
                  }
                }
                
                if(type_med){
                  
                  for(int iter_d=0 ; iter_d<D ; iter_d++){
                    
                    // vecteur des donnees
                    valeurs_miroir.resize(nb_vois);
                     for(int iter_vois=0 ;  iter_vois<nb_vois ; iter_vois++){
                          valeurs_miroir[iter_vois] = contrast(index_vois[iter_vois],iter_d);
                     }
                    
                     valeur_NormContro(iter_px,iter_d) = contrast(index_px,iter_d)-medianVecDouble_cpp(valeurs_miroir); //medianNV_cpp(valeurs_miroir);
                  }
                }
                
                if(type_NN){
                  index_contro = index_vois[0];
                  dist_vois = sqrt(pow(coords_px(index_contro,0)+coords_px(index_px,0),2) + pow(coords_px(index_contro,1)-coords_px(index_px,1),2));
                  diff_ref = abs(contrast(index_px,param_ref) - contrast(index_contro,param_ref))/var_ref + lambda*dist_vois;
                  index_ref_contro = index_contro;

                  // identification du pixel le plus proche en intensite
                  for(int iter_vois=1 ;  iter_vois<nb_vois ; iter_vois++){
                    index_contro = index_vois[iter_vois];
                    dist_vois = sqrt(pow(coords_px(index_contro,0)+coords_px(index_px,0),2) + pow(coords_px(index_contro,1)-coords_px(index_px,1),2));
                         
                    diff_tempo = abs(contrast(index_px,param_ref) - contrast(index_contro,param_ref))/var_ref + lambda*dist_vois;
                    
                    if(diff_tempo - diff_ref < - 0.0000000000001){
                    diff_ref = diff_tempo;
                    index_ref_contro = index_contro;
                    }                  
                  }
                    for(int iter_d=0 ; iter_d<D ; iter_d++){                  
                    valeur_NormContro(iter_px,iter_d) = contrast(index_px,iter_d)-contrast(index_ref_contro,iter_d);                 
                    }
                   
                 }
              
                }
                }
                
                return(List::create(Named("valeur_NormContro")  = valeur_NormContro,
                          Named("index_plot_k")  = index_plot_k,
                          Named("index_plot_k_contro")  = index_plot_k_contro
                ));
}

//  3 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List filtrage2D_cpp(const arma::mat& M_data, const arma::mat& M_operateur, const arma::mat& index_data,
                    bool w_contrast, bool na_rm){

  int const n_data = index_data.n_rows;
  vector<int>  p_data(2);
  p_data[0]=M_data.n_rows ; p_data[1] = M_data.n_cols ;
  vector<int>  p_operateur(2);
  p_operateur[0] = M_operateur.n_rows ; p_operateur[1] = M_operateur.n_cols;

  vector<int> v0(p_operateur[0]);
  for(int iter=0 ; iter<p_operateur[0] ; iter++){v0[iter] = iter;}
  vector<int> v1(p_operateur[1]);
  for(int iter=0 ; iter<p_operateur[1] ; iter++){v1[iter] = iter;}
  vector<double> p_ref(2); // vector<int> p_ref(2);
  p_ref[0] = medianVecInteger_cpp(v0) ; p_ref[1] = medianVecInteger_cpp(v1); 

  NumericVector M_vector(n_data);
  for(int iter_px=0 ; iter_px<n_data ; iter_px++)
  {M_vector(iter_px) = M_data(index_data(iter_px,0),index_data(iter_px,1)) ;}
  
  double const sigma = sqrt(var(M_vector));
  double const dnorm_sigma = R::dnorm(0.0,0.0,sigma,false); 
 
 
  arma::mat Mres(p_data[0],p_data[1]);  
  std::fill(Mres.begin(),Mres.end(),NA_REAL);
  arma::mat Wres(p_data[0],p_data[1]);  
  std::fill(Wres.begin(),Wres.end(),NA_REAL);
  bool test1a, test1b, test2a, test2b, break_loop ; 
  double tempo, sumNNA,Mintensite, Mvoisin ;
  
  for(int iter_px=0 ; iter_px<n_data ; iter_px++){
    // initialisation
    Mres(index_data(iter_px,0),index_data(iter_px,1)) = 0;
    sumNNA = 0;
    break_loop = false;
      
    // mise en place des matrices
    for(int iter_l=0 ; iter_l<p_operateur[0] ; iter_l++){

      for(int iter_c=0 ; iter_c<p_operateur[1] ; iter_c++){   

         Mintensite=1;
         Mvoisin=NA_REAL;
        
         test1a = index_data(iter_px,0)+iter_l-p_ref[0] >= 0;
         test1b = index_data(iter_px,0)+iter_l-p_ref[0] < p_data[0];
         test2a = index_data(iter_px,1)+iter_c-p_ref[1] >= 0;
         test2b = index_data(iter_px,1)+iter_c-p_ref[1] < p_data[1];
    
    if(test1a*test1b*test2a*test2b == 1 && R_IsNA(M_data(index_data(iter_px,0)+iter_l-p_ref[0],index_data(iter_px,1)+iter_c-p_ref[1]))==0){      
      // elements voisins
      if(R_IsNA(M_operateur(iter_l,iter_c))==0){
        Mvoisin = M_data(index_data(iter_px,0)+iter_l-p_ref[0],
        index_data(iter_px,1)+iter_c-p_ref[1]);
        
        // intensite des pixels
        if(w_contrast==1 && M_operateur(iter_l,iter_c)!=0){
          tempo = Mvoisin-M_data(index_data(iter_px,0),index_data(iter_px,1));
          Mintensite = R::dnorm(tempo,0.0,sigma,false)/dnorm_sigma;
        }
        
      }
    }
     
    // ajustement des NA et maj
    if(R_IsNA(Mvoisin)==0){
      sumNNA = sumNNA + abs(M_operateur(iter_l,iter_c)*Mintensite);
 
      Mres(index_data(iter_px,0),index_data(iter_px,1)) = Mres(index_data(iter_px,0),index_data(iter_px,1)) + M_operateur(iter_l,iter_c)*Mintensite*Mvoisin;
    }else{
      
      if(na_rm){
        Mres(index_data(iter_px,0),index_data(iter_px,1)) = NA_REAL;
        break_loop = true;
        break;
      }
            
    }
    }
    if(break_loop){break;}
    }    

    // correction des NA
    Wres(index_data(iter_px,0),index_data(iter_px,1)) = sumNNA;
//    if(break_loop==false && sumNNA>0)
//    {Mres(index_data(iter_px,0),index_data(iter_px,1)) = Mres(index_data(iter_px,0),index_data(iter_px,1))/sumNNA;}
 }  

    return(List::create(Named("Mres")  = Mres,
                        Named("Wres")  = Wres                          
                ));
}

//  4 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List filtrage3D_cpp(const NumericVector& Vec_data, const IntegerVector& p_data,
                    const NumericVector& Vec_operateur, const IntegerVector& p_operateur, 
                    const arma::mat& index_data, bool w_contrast, bool na_rm){


  arma::cube A_data(Vec_data.begin(), p_data[0], p_data[1], p_data[2]);
  arma::cube A_operateur(Vec_operateur.begin(), p_operateur[0], p_operateur[1], p_operateur[2]);
  int const n_data = index_data.n_rows;
  
  vector<int> v0(p_operateur[0]);
  for(int iter=0 ; iter<p_operateur[0] ; iter++){v0[iter] = iter;}
  vector<int> v1(p_operateur[1]);
  for(int iter=0 ; iter<p_operateur[1] ; iter++){v1[iter] = iter;}
  vector<int> v2(p_operateur[2]);
  for(int iter=0 ; iter<p_operateur[2] ; iter++){v2[iter] = iter;}
  vector<double> p_ref(3); // vector<int> p_ref(3);
  p_ref[0] = medianVecInteger_cpp(v0) ; p_ref[1] = medianVecInteger_cpp(v1);  p_ref[2] = medianVecInteger_cpp(v2);

  // variance des donnees pour le lissage
  NumericVector M_vector(n_data);
  for(int iter_px=0 ; iter_px<n_data ; iter_px++)
  {M_vector(iter_px) = A_data(index_data(iter_px,0),index_data(iter_px,1),index_data(iter_px,2)) ;}
  
  double const sigma = sqrt(var(M_vector));
  double const dnorm_sigma = R::dnorm(0.0,0.0,sigma,false); 
 
  // stockage resultats
  arma::cube Mres(p_data(0),p_data(1),p_data(2));
  std::fill(Mres.begin(),Mres.end(),NA_REAL);
  arma::cube Wres(p_data(0),p_data(1),p_data(2));
  std::fill(Wres.begin(),Wres.end(),NA_REAL);
  
  bool test1a, test1b, test2a, test2b, test3a, test3b, break_loop ; 
  double tempo, sumNNA, Mvoisin, Mintensite;
    
  for(int iter_px=0 ; iter_px<n_data ; iter_px++){
  
    // initialisation
    Mres(index_data(iter_px,0),index_data(iter_px,1),index_data(iter_px,2)) = 0;
    sumNNA = 0;
    break_loop = false;
       
    // mise en place des matrices
    for(int iter_l=0 ; iter_l<p_operateur[0] ; iter_l++){
      
      for(int iter_c=0 ; iter_c<p_operateur[1] ; iter_c++){      
        
         for(int iter_d=0 ; iter_d<p_operateur[2] ; iter_d++){ 
           Mintensite=1;
           Mvoisin=NA_REAL;
    
         test1a = index_data(iter_px,0)+iter_l-p_ref[0] >= 0;
         test1b = index_data(iter_px,0)+iter_l-p_ref[0] < p_data[0];
         test2a = index_data(iter_px,1)+iter_c-p_ref[1] >= 0;
         test2b = index_data(iter_px,1)+iter_c-p_ref[1] < p_data[1];
         test3a = index_data(iter_px,2)+iter_d-p_ref[2] >= 0;
         test3b = index_data(iter_px,2)+iter_d-p_ref[2] < p_data[2];
      
    
    if(test1a*test1b*test2a*test2b*test3a*test3b == 1 && R_IsNA(A_data(index_data(iter_px,0)+iter_l-p_ref[0],index_data(iter_px,1)+iter_c-p_ref[1],index_data(iter_px,2)+iter_d-p_ref[2]))==0){
      // elements voisins
       if(R_IsNA(A_operateur(iter_l,iter_c,iter_d))==0){
              Mvoisin = A_data(index_data(iter_px,0)+iter_l-p_ref[0],
                               index_data(iter_px,1)+iter_c-p_ref[1],
                               index_data(iter_px,2)+iter_d-p_ref[2]);
            
      // intensite des pixels
      if(w_contrast==1){
      tempo = Mvoisin-A_data(index_data(iter_px,0),index_data(iter_px,1),index_data(iter_px,2));
      Mintensite = R::dnorm(tempo,0.0,sigma,false)/dnorm_sigma;
      }
       
      }
      
    }

    // ajustement des NA et maj
    if(R_IsNA(Mvoisin)==0){
        sumNNA = sumNNA + abs(A_operateur(iter_l,iter_c,iter_d)*Mintensite);

        Mres(index_data(iter_px,0),index_data(iter_px,1),index_data(iter_px,2)) = Mres(index_data(iter_px,0),index_data(iter_px,1),index_data(iter_px,2)) + A_operateur(iter_l,iter_c,iter_d)*Mintensite*Mvoisin;
                     
    }else{
      
      if(na_rm)
      { Mres(index_data(iter_px,0),index_data(iter_px,1),index_data(iter_px,2)) = NA_REAL;
        break_loop = true;
        break;
      }
            
    }
        }
         if(break_loop){break;}
    }
    if(break_loop){break;}
    }    

    // correction des NA
    Wres(index_data(iter_px,0),index_data(iter_px,1),index_data(iter_px,2)) = sumNNA;
//    if(break_loop==false && sumNNA>0)
//    {Mres(index_data(iter_px,0),index_data(iter_px,1),index_data(iter_px,2)) = Mres(index_data(iter_px,0),index_data(iter_px,1),index_data(iter_px,2))/sumNNA;}

 }  

    return(List::create(Named("Mres")  = Mres,
                        Named("Wres")  = Wres                          
                ));
}

//  5 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat filtrage2Dmed_cpp(const arma::mat& M_data, const arma::mat& M_operateur, const arma::mat& index_data,
                            bool na_rm){
  
  int const n_data = index_data.n_rows;
  vector<int>  p_data(2);
  p_data[0]=M_data.n_rows ; p_data[1] = M_data.n_cols ;
  vector<int>  p_operateur(2);
  p_operateur[0] = M_operateur.n_rows ; p_operateur[1] = M_operateur.n_cols;

  vector<int> v0(p_operateur[0]);
  for(int iter=0 ; iter<p_operateur[0] ; iter++){v0[iter] = iter;}
  vector<int> v1(p_operateur[1]);
  for(int iter=0 ; iter<p_operateur[1] ; iter++){v1[iter] = iter;}
  vector<double> p_ref(2); //vector<int> p_ref(2);
  p_ref[0] = medianVecInteger_cpp(v0) ; p_ref[1] = medianVecInteger_cpp(v1);
  vector<double> Mvoisin ;

  arma::mat Mres(p_data[0],p_data[1]);
  std::fill(Mres.begin(),Mres.end(),NA_REAL);
  bool test1a, test1b, test2a, test2b, break_loop ; 
 
  for(int iter_px=0 ; iter_px<n_data ; iter_px++){
    
    // initialisation
    Mres(index_data(iter_px,0),index_data(iter_px,1)) = 0;
    Mvoisin.resize(0);
    break_loop = false;
      
      
    // mise en place des matrices
    for(int iter_l=0 ; iter_l<p_operateur[0] ; iter_l++){
      
      for(int iter_c=0 ; iter_c<p_operateur[1] ; iter_c++){      
         test1a = index_data(iter_px,0)+iter_l-p_ref[0] >= 0;
         test1b = index_data(iter_px,0)+iter_l-p_ref[0] < p_data[0];
         test2a = index_data(iter_px,1)+iter_c-p_ref[1] >= 0;
         test2b = index_data(iter_px,1)+iter_c-p_ref[1] < p_data[1];
                          
    // elements voisins
    if( test1a*test1b*test2a*test2b == 1 && R_IsNA(M_data(index_data(iter_px,0)+iter_l-p_ref[0],index_data(iter_px,1)+iter_c-p_ref[1]))==0){
          
      if(R_IsNA(M_operateur(iter_l,iter_c))==0  && M_operateur(iter_l,iter_c)!=0){
      Mvoisin.push_back(M_data(index_data(iter_px,0)+iter_l-p_ref[0],
                          index_data(iter_px,1)+iter_c-p_ref[1]));
      }
      
    }else{
      if(na_rm){
        Mres(index_data(iter_px,0),index_data(iter_px,1))= NA_REAL;
        break_loop = true;
        break;                       
      }
    }
        
    }
    if(break_loop){break;}
    }    
  
    // calcul de la mediane
    if(R_IsNA(Mres(index_data(iter_px,0),index_data(iter_px,1)))==0){
    Mres(index_data(iter_px,0),index_data(iter_px,1)) = medianVecDouble_cpp(Mvoisin);
    }
    
  
 }  

 return Mres;
}

//  6 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
NumericVector filtrage3Dmed_cpp(const NumericVector& Vec_data, const IntegerVector& p_data,
                       const NumericVector& Vec_operateur, const IntegerVector& p_operateur, 
                       const arma::mat& index_data, bool na_rm){
  
    // conversion en array
   arma::cube A_data(Vec_data.begin(), p_data[0], p_data[1], p_data[2]);
   arma::cube A_operateur(Vec_operateur.begin(), p_operateur[0], p_operateur[1], p_operateur[2]);
   int const n_data = index_data.n_rows;

  vector<int> v0(p_operateur[0]);
  for(int iter=0 ; iter<p_operateur[0] ; iter++){v0[iter] = iter;}
  vector<int> v1(p_operateur[1]);
  for(int iter=0 ; iter<p_operateur[1] ; iter++){v1[iter] = iter;}
  vector<int> v2(p_operateur[2]);
  for(int iter=0 ; iter<p_operateur[2] ; iter++){v2[iter] = iter;}

  vector<double> p_ref(3);  //vector<int> p_ref(3);
   p_ref[0] = medianVecInteger_cpp(v0) ; p_ref[1] = medianVecInteger_cpp(v1);  p_ref[2] = medianVecInteger_cpp(v2);

  vector<double> Mvoisin ;
 
  arma::cube Mres(p_data[0],p_data[1],p_data[2]);
  std::fill(Mres.begin(),Mres.end(),NA_REAL);
    
    bool test1a, test1b, test2a, test2b, test3a, test3b, break_loop ; 
   
      for(int iter_px=0 ; iter_px<n_data ; iter_px++){
          
            // initialisation
          Mres(index_data(iter_px,0),index_data(iter_px,1),index_data(iter_px,2)) = 0;
          Mvoisin.resize(0);
          break_loop = false;
          
            // mise en place des matrices
            for(int iter_l=0 ; iter_l<p_operateur[0] ; iter_l++){
               
                for(int iter_c=0 ; iter_c<p_operateur[1] ; iter_c++){      
                  
                    for(int iter_d=0 ; iter_d<p_operateur[2] ; iter_d++){ 
                       
                     test1a = index_data(iter_px,0)+iter_l-p_ref[0] >= 0;
                     test1b = index_data(iter_px,0)+iter_l-p_ref[0] < p_data[0];
                     test2a = index_data(iter_px,1)+iter_c-p_ref[1] >= 0;
                     test2b = index_data(iter_px,1)+iter_c-p_ref[1] < p_data[1];
                     test3a = index_data(iter_px,2)+iter_d-p_ref[2] >= 0;
                     test3b = index_data(iter_px,2)+iter_d-p_ref[2] < p_data[2];
                                      
                  // elements voisins
                if( test1a*test1b*test2a*test2b*test3a*test3b == 1 && R_IsNA(A_data(index_data(iter_px,0)+iter_l-p_ref[0],index_data(iter_px,1)+iter_c-p_ref[1],index_data(iter_px,2)+iter_d-p_ref[2]))==0){
                        
                      if(R_IsNA(A_operateur(iter_l,iter_c,iter_d))==0  && A_operateur(iter_l,iter_c,iter_d)!=0){
                        Mvoisin.push_back(A_data(index_data(iter_px,0)+iter_l-p_ref[0],
                                            index_data(iter_px,1)+iter_c-p_ref[1],
                                            index_data(iter_px,2)+iter_d-p_ref[2]));
                        }
                    
                    }else{
                        if(na_rm){
                            Mres(index_data(iter_px,0),index_data(iter_px,1),index_data(iter_px,2))= NA_REAL;
                            break_loop = true;
                            break;                       
                          }
                      }
                     }                    
                    if(break_loop){break;}
                   }                  
            if(break_loop){break;}
            }   
            
        
            // calcul de la mediane
          if(R_IsNA(Mres(index_data(iter_px,0),index_data(iter_px,1),index_data(iter_px,2)))==0)
            {Mres(index_data(iter_px,0),index_data(iter_px,1),index_data(iter_px,2)) = medianVecDouble_cpp(Mvoisin);}
        
         }  
   
   return(wrap(Mres));
}

//  7 ////////////////////////////////////////////////////////////
double medianIV_cpp(IntegerVector data){
  
  // taille du vecteur
  int taille = data.size() ;
  if (taille==0) throw "Empty";
  
  // rangement du vecteur
  std::sort(data.begin(), data.end()) ;
  
  // calcul de la mediane
  double mediane ;
  if (taille  % 2 == 0)
  {   mediane = (data[ taille / 2 - 1 ] + data[taille / 2]) / 2.0; }
  else 
  {   mediane = data[taille/2] ;  }
  
  return mediane;
}


//  8 ////////////////////////////////////////////////////////////
double medianNV_cpp(NumericVector data){
   
    // taille du vecteur
    int taille = data.size() ;
    if (taille==0) throw "Empty";
    
    // rangement du vecteur
    std::sort(data.begin(), data.end()) ;
    
    // calcul de la mediane
    double mediane ;
    if (taille  % 2 == 0)
    {   mediane = (data[ taille / 2 - 1 ] + data[taille / 2]) / 2.0; }
    else 
    {   mediane = data[taille/2] ;  }
    
    return mediane;
  }
  
  
//  9 ////////////////////////////////////////////////////////////
double medianVecInteger_cpp(std::vector<int> data){
    
   // taille du vecteur
    unsigned int taille = data.size() ;
    if (taille==0) throw "Empty";
    
    // rangement du vecteur
    std::sort(data.begin(), data.end()) ;
    
    // calcul de la mediane
    double mediane ;
    if (taille  % 2 == 0)
    {   mediane = (data[ taille / 2 - 1 ] + data[taille / 2]) / 2.0; }
    else 
    {   mediane = data[taille/2] ;  }
    
    return mediane;
}   

//  10 ////////////////////////////////////////////////////////////
double medianVecDouble_cpp(std::vector<double> data){
    
   // taille du vecteur
    unsigned int taille = data.size() ;
    if (taille==0) throw "Empty";
    
    // rangement du vecteur
    std::sort(data.begin(), data.end()) ;
    
    // calcul de la mediane
    double mediane ;
    if (taille  % 2 == 0)
    {   mediane = (data[ taille / 2 - 1 ] + data[taille / 2]) / 2.0; }
    else 
    {   mediane = data[taille/2] ;  }
    
    return mediane;
}   

//  11 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List calcGroupsCoords_cpp(const arma::mat& coords_NNA, const IntegerVector& index_NNA, const arma::mat& Neighborhood, IntegerVector coords_max,
                          int max_groups, bool trace){
  
   int n=coords_NNA.n_rows; // nb de vx
   int p_Mneighbors=Neighborhood.n_cols; // dim du voisinage
   int n_Mneighbors=Neighborhood.n_rows; // taille du voisinage
  vector< int> cum_coords_max(coords_max.size());
  cum_coords_max[0] = 1;
  for (int iter=1; iter<coords_max.size(); iter++){cum_coords_max[iter]=cum_coords_max[iter-1]*coords_max[iter-1];}  
 
  int iter_group=0; // iteration sur les groupes
  int iter_res; // iteration pour le stockage
  int n_neighbors,n_neighbors_new; // nb de nouveaux voisins
  int Sum_group_size=0 ;
  int index_tempo ;
  double sum_tempo ; // int sum_tempo ;
  unsigned int iter_vois_min ;
  
  vector<int> group_size(0); // taille de chaque groupe
  vector<int> residual(n); // index des voxels a classer
  vector<int> coords_newVoisin(0);
  for(int iter=0 ; iter<n ; iter++){residual[iter] = iter;}
  int n_affected=0; // iteration sur les voxels
  arma::mat group(n,2);  // definition des groupes  
  bool test_newVoisin;
  
  // trace
  double seq_trace = n*0.1;
       
  while(n_affected<n && iter_group<=max_groups){ // tant qu il reste des vx a classer et que le nb max de groupe n est pas atteint
    iter_group = iter_group + 1 ;
    
    // ajout d un nouveau groupe
    group(n_affected,0) = residual[0];
    group(n_affected,1) = iter_group ;
    n_affected = n_affected + 1 ;  
    residual.erase(residual.begin());    
      
    group_size.push_back(group_size.size());
    group_size[group_size.size()-1] = 1 ;
    n_neighbors = 1;
               
    while(n_neighbors>0 && n_affected<n){ // tant qu il reste des nouveaux voisins
      n_neighbors_new=0;
      coords_newVoisin.resize(0);
      
       for(int iter_vois=Sum_group_size+group_size[iter_group-1]-n_neighbors ; iter_vois < Sum_group_size+group_size[iter_group-1] ; iter_vois++){ // pour chaque voisin
    
            for(int iter_l=0 ; iter_l < n_Mneighbors ; iter_l++){ // identifier les nouveaux voisins potentiels
            test_newVoisin = true;    
            index_tempo=0;
            
                for(int iter_p=0 ; iter_p < p_Mneighbors ; iter_p++){
                
                  sum_tempo = coords_NNA(group(iter_vois,0),iter_p) + Neighborhood(iter_l,iter_p);
                 
                  if(sum_tempo <0 || sum_tempo >= coords_max[iter_p] ){ 
                    test_newVoisin = false;
                    break;                     
                  }else{
                  index_tempo = index_tempo + sum_tempo*cum_coords_max[iter_p] ;
                  }                 
                  
                }
                
                if(test_newVoisin){  
                  coords_newVoisin.push_back(coords_newVoisin.size());  
                  coords_newVoisin[coords_newVoisin.size()-1] = index_tempo;
                 } 
            }
        }  

        std::sort(coords_newVoisin.begin(), coords_newVoisin.end());
        iter_vois_min=0;       
        
      iter_res = 0 ;
      while(iter_res<n-n_affected && iter_vois_min<coords_newVoisin.size()){
              
              while(index_NNA[residual[iter_res]]>coords_newVoisin[iter_vois_min] && iter_vois_min<coords_newVoisin.size()){
                coords_newVoisin.erase(coords_newVoisin.begin() + iter_vois_min); // iter_vois_min++;                
              }
              
              if(index_NNA[residual[iter_res]]==coords_newVoisin[iter_vois_min]){ // matching               
               // attribution du groupe au voxel
                 group(n_affected,0) = residual[iter_res];
                 group(n_affected,1) = iter_group ;
                 n_affected++ ;  
               
                // maj de residual et neighbors_new
                residual.erase(residual.begin() + iter_res);
                n_neighbors_new++ ;
                iter_vois_min++;                                
              }else{
                iter_res++;
              }              
      }
      
      n_neighbors = n_neighbors_new;
      group_size[iter_group-1] = group_size[iter_group-1] + n_neighbors;
      
      if(trace && n_affected<n && (n_affected > seq_trace) ){
        Rcout << "*" ;
        seq_trace = seq_trace + n*0.1;
      }
      
    }
    
    Sum_group_size = Sum_group_size + group_size[iter_group-1];
    
  }  
  if(trace){Rcout << endl ;}
  
  
  
 //  export
  return(List::create(Named("group")  = group,
                     Named("group_size")  = group_size)                     
      );
  
}

//  12 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List calcGroupsW_cpp(const S4& W, const IntegerVector& subset, int max_groups){
  
  int n = subset.size(); 

  // cas pathologique
  if(n==1){    
    return(List::create(Named("group_subset")  = 1,
                        Named("group_size")  = 1)                     
    );
  }
  
  // subset
  List res_subset=subsetW_cpp(W.slot("i"),W.slot("p"),subset);
  
  vector<int> W_i=res_subset[0];  
  vector<int> W_p=res_subset[1]; 
  
  // initialization
  vector<int> group(n,-1); 
  vector<int> groupeSize(0);  
  vector<int> indexObs(n); 
  
  // loop
  for(int iter=0 ; iter<n ; iter++){indexObs[iter]=iter;}
  
  vector<int> voisin_0; 
  vector<int> voisin_1; 
  
  unsigned int iter_newvois, iter_index;
    
  // iteration sur les groupes
  int iter_group=0; 
  
  while(iter_group < max_groups && indexObs.size()>0){
    
    group[indexObs[0]] = iter_group + 1 ; 
    groupeSize.push_back(1);
    
    voisin_0.resize(1);
    voisin_0[0]=indexObs[0]; 
    indexObs.erase( indexObs.begin() );
    
    while(voisin_0.size()>0 && indexObs.size()>0){ // tant qu il y a des nouveaux voisins
    voisin_1.resize(0);
    
    for(unsigned int iter_vois=0 ; iter_vois<voisin_0.size() ; iter_vois++){ // chercher les nouveaux voisins
    
    if(abs(W_p[voisin_0[iter_vois]] - W_p[voisin_0[iter_vois]+1])>0.1){ // si le point a des voisins
    
    for(int iter_newvois=W_p[voisin_0[iter_vois]]; iter_newvois<W_p[voisin_0[iter_vois]+1]; iter_newvois++){
      voisin_1.push_back(W_i[iter_newvois]);  
    }
    }
    
    }    
    
    std::sort(voisin_1.begin(), voisin_1.end()); // rangement par ordre croissant
    
    iter_newvois=0;
    iter_index=0;
    
    while(iter_newvois<voisin_1.size() && iter_index < indexObs.size()){
      
      while(voisin_1[iter_newvois] < indexObs[iter_index] && iter_newvois<voisin_1.size()){
        voisin_1.erase(voisin_1.begin() + iter_newvois); 
      }     
      
      if(voisin_1[iter_newvois] == indexObs[iter_index]){ 
        
        group[voisin_1[iter_newvois]] = iter_group + 1;
        indexObs.erase(indexObs.begin() + iter_index);   
        iter_newvois++; 
        
      }else{
        iter_index++;
      }
      
    }
    
    if(iter_newvois<voisin_1.size()){    
      voisin_1.erase(voisin_1.begin() + iter_newvois,voisin_1.begin() + voisin_1.size());          
    }
    
    voisin_0 = voisin_1;
    groupeSize[iter_group] = groupeSize[iter_group] + voisin_0.size();
    }
    
    iter_group = iter_group + 1 ;
    
  }
  
  
  // export
  int group_length=groupeSize.size();
  int nbObsRes=indexObs.size();
  
  return(List::create(
  Named("group_subset")  = group,                    
  Named("group_size")  = groupeSize,
  Named("group_length") = group_length,
  Named("nbObsRes")  = nbObsRes)                     
  );
}

//  13 ////////////////////////////////////////////////////////////
List subsetW_cpp(const IntegerVector& W_i, const IntegerVector& W_p, IntegerVector subset){
  
  std::sort(subset.begin(), subset.end());
  
  int n_subset = subset.size();
  int index_min,index_max;
  vector<int> correspondance(W_p.size()-1,0);
  for(int iter_subset=0; iter_subset<n_subset ; iter_subset++){
    correspondance[subset[iter_subset]] = iter_subset;
  }
  vector<int> index_i,index_i_sort;
  vector<int> Wnew_p(n_subset+1,0);
  vector<int> Wnew_i;  
  
  for(int iter_subset=0; iter_subset<n_subset ; iter_subset++){
    
    index_min = W_p[subset[iter_subset]];
    index_max = W_p[subset[iter_subset]+1];
    
    index_i.resize(index_max-index_min);
    for(int iter_i=index_min ; iter_i<index_max ; iter_i++){
      index_i[iter_i-index_min]= W_i[iter_i];
    }
    
    index_i_sort.resize(0);
    std::sort(index_i.begin(), index_i.end());
    
    set_intersection(index_i.begin(),index_i.end(),subset.begin(),subset.end(),back_inserter(index_i_sort));
    
    if(index_i_sort.size()>0){
      for(unsigned int iter_i_sort=0 ; iter_i_sort < index_i_sort.size() ; iter_i_sort++){
        Wnew_i.push_back(correspondance[index_i_sort[iter_i_sort]]);
        index_i_sort[iter_i_sort] = correspondance[index_i_sort[iter_i_sort]];
      }
    }
    
    Wnew_p[iter_subset+1] = Wnew_p[iter_subset]+index_i_sort.size();
    
  }
  
  return(List::create(Named("W_i")  = Wnew_i,
  Named("W_p")  = Wnew_p));       
  
}

//  13 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List calcRadius_cpp(const arma::mat& coords, const NumericVector& sample, double threshold, const LogicalVector& subset_bary, bool trace){
  
  // preparation
  int n = sample.size();
  vector<int> Groupe_potentiel; 
  for(int iter_n=0; iter_n<n; iter_n++)
  { if(sample(iter_n)>threshold && subset_bary(iter_n))
  {Groupe_potentiel.push_back(iter_n);}       
  }
  int n_Gp = Groupe_potentiel.size();
  
  // cas particulier
  if(n_Gp<2){
  if(trace==true){Rcout << "Rayon : " << 0 << endl;}  
  return(0);
  }
  
  int index_Gp;     
  double Rayon_medAVC=0, norm=0;
  
  // calcul du rayon - moyenne ponderee des distances au barycentre
  int D = coords.n_cols;
 
  vector<double> coord_bary(D);
  std::fill(coord_bary.begin(),coord_bary.end(),0);  
  // calcul du barycentre
  for(int iter_Gp=0 ; iter_Gp<n_Gp ; iter_Gp++){ 
    index_Gp = Groupe_potentiel[iter_Gp];  

    for( int iter_d=0; iter_d < D ; iter_d++){
      coord_bary[iter_d] = coord_bary[iter_d] + coords(index_Gp,iter_d)*sample[index_Gp];
    }
    norm = norm + sample[index_Gp];
  }
  for( int iter_d=0; iter_d < D ; iter_d++){
    coord_bary[iter_d] = coord_bary[iter_d]/norm;
  }

  // moyenne des distances
  norm=0;
  vector<double> distance(n,0);
  
  for(int iter_Gp=0 ; iter_Gp<n_Gp ; iter_Gp++){ 
    index_Gp = Groupe_potentiel[iter_Gp];
    
    for( int iter_d=0; iter_d < D ; iter_d++){
      distance[index_Gp] = distance[index_Gp]+pow(coord_bary[iter_d]-coords(index_Gp,iter_d),2);      
    }
    distance[index_Gp] = sqrt(distance[index_Gp]);
    Rayon_medAVC = Rayon_medAVC + distance[index_Gp]*sample[index_Gp];
    norm = norm + sample(index_Gp);
  }
  
  
  Rayon_medAVC = Rayon_medAVC/norm ;
  if(trace==true){Rcout << "Rayon : " << Rayon_medAVC << endl;}
  
   return(List::create(Named("radius") = Rayon_medAVC,
                       Named("distance") = distance,
                       Named("coord_bary") = coord_bary));   
  
}
