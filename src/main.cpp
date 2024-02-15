// this code contains functions to estimate the parameters of six types of LCPLFMs
// model M_1 with constant classification of object parameters and class-specific object parameters (disjunctive mapping rule)
// model M_2 with constant classification of object parameters and class-specific attribute parameters (disjunctive mapping rule)
// model M_3 with constant classification of object parameters and class-specific object- and attribute parameters (disjunctive mapping rule)

// model M_1 with constant classification of object parameters and class-specific object parameters  (additive mapping rule)
// model M_2 with constant classification of object parameters and class-specific attribute parameters (additive mapping rule)
// model M_3 with constant classification of object parameters and class-specific object- and attribute parameters (additive mapping rule)


// LCPLFMs are applied to a binary data array with  nO objects, nA attributes and nR persons
// the object, attribute and class size parameters of LCPLFMs are indicated as ro, ta and ga, respectively


#include <iostream>
#include <cstring>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> // for NULL
#include <fstream>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <math.h>
#include <cmath>
#include <R_ext/Rdynload.h> 

using namespace std;

// data structures
typedef bool binary;
typedef bool ** binary2;
typedef bool *** binary3;
typedef long double **** extended4;
typedef long double *** extended3;
typedef long double ** extended2;
typedef long double * extended1;

// data types
int seed1, seed2, seed3;
int _nO; // number of objects
int _nA; // number of attributes
int _nR; // number of persons
int _nruns; 
long double emcrit;
long double delta; // The precision used to compute standard errors
unsigned short int _nF,_nT,_nS; // number of latent features _nF, number of latent classes _nT, number of feature patterns _nS=2^_nF
long double logold,lognew,conv,logpost_n,loglik_n;
binary3 dataB; // _nO x _nA x _nR binary array of observed data
binary2 patS; // _nS x _nF matrix of binary latent data patterns 
extended4 omega;
// bootstrap vars
int _Nbootstrap;


//C++ functions
long double pow_1(long double, binary);
long double pow_2(long double, binary);
long double pow_3(unsigned short int, unsigned short int, extended2, binary);
long double pow_3(unsigned short int, unsigned short int,unsigned short int, extended3, binary);
long double logposterior (extended3, extended2, extended1, extended2, extended3,extended3, extended2);
long double logposterior (extended2, extended3, extended1, extended3, extended2, extended2, extended3);
long double logposterior (extended3, extended3, extended1, extended3, extended3, extended3, extended3);

void calculate_dims(void);
void calculate_condprobX_DC(extended2 , extended2 );
void calculate_condprobX_DC(extended3 , extended3 );
void calculate_condprobX_ADD(extended2 , extended2 );
void calculate_condprobX_ADD(extended3 , extended3 );
void calculate_margprobX(extended3 , extended3 );
void calculate_margprobX(extended2 , extended2 );
void replace(extended3 ,extended3 ,extended2 ,extended2 ,extended1 ,extended1 , extended3, extended2, extended2 ,extended2 ,extended3 ,extended3 );
void replace(extended2 ,extended2 ,extended3 ,extended3 ,extended1 ,extended1 , extended2, extended3, extended3 ,extended3 ,extended2 ,extended2 );
void replace(extended3 ,extended3 ,extended3 ,extended3 ,extended1 ,extended1 , extended3, extended3, extended3 ,extended3 ,extended3 ,extended3 );
void calculate_probmat(extended1 , extended2 , extended3,extended2  );
void calculate_probmat(extended1 , extended3 , extended2,extended2  );
void calculate_probmat(extended1 , extended3 , extended3,extended2  );
void emgamma( extended1,extended2  );
void calculate_condprob_pattern(extended2 , extended3,extended3 );
void calculate_condprob_pattern(extended3 , extended2,extended3 );
void calculate_condprob_pattern(extended3 , extended3,extended3 );
void update_emro( extended3,extended2, extended3 );
void update_emro( extended2,extended2, extended2 );
void update_emta_DC(extended2, extended2, extended2, extended2, extended2);
void update_emta_DC(extended3, extended3, extended3, extended2, extended3);
void update_emta_ADD(extended2, extended2, extended2, extended2, extended2);
void update_emta_ADD(extended3, extended3, extended3, extended2, extended3);
void calculate_conv(extended3 ,extended3 ,extended2 ,extended2 ,extended1 ,extended1, extended3, extended2);
void calculate_conv(extended2 ,extended2 ,extended3 ,extended3 ,extended1 ,extended1, extended2, extended3 );
void calculate_conv(extended3 ,extended3 ,extended3 ,extended3 ,extended1 ,extended1, extended3, extended3 );
extern "C" void PlFm_XZ_Y_intern_DC(binary3  R_data , int R_nR , int R_nO , int R_nA , int R_nF, int R_nT, binary2 R_patS,							
		  long double R_emcrit, extended3 R_ro_o, extended2 R_ta_o, extended1 R_ga_o, 
          extended3 R_ro_n, extended2 R_ta_n, extended1 R_ga_n, extended3 R_ro_update, extended2 R_ta_update, bool flag, double R_delta,
		  extended3 margprobx_o,extended3 margprobx_n,extended3 somega,extended3 gradient_ro,extended3 se_ro,extended3 p_o_r_t,
		  extended2 probmat,extended2 condprobx_o,extended2 condprobx_n,extended2 gradient_ta,extended2 se_ta,extended2  p_r_t,
		 extended1 gradient_ga,extended1 p_r,extended1 se_ga, 
		 int R_Nbootstrap, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att, extended2 report_common_obj);
extern "C" void PlFm_X_YZ_intern_DC(binary3 R_data, int R_nR, int R_nO, int R_nA, int R_nF, int R_nT, binary2 R_patS, 
		  long double R_emcrit,extended2 R_ro_o, extended3 R_ta_o, extended1 R_ga_o,
		  extended2 R_ro_n, extended3 R_ta_n, extended1 R_ga_n, extended2 R_ro_update, extended3 R_ta_update, bool flag, double R_delta,
		 extended2 margprobx_o,extended2 margprobx_n,extended3 somega,extended2 gradient_ro,extended2  se_ro,extended3  p_o_r_t,
		 extended2 probmat,extended3 condprobx_o,extended3 condprobx_n,extended3 gradient_ta,extended3 se_ta,extended2 p_r_t,
		 extended1 gradient_ga,extended1 p_r,extended1 se_ga,
		 int R_Nbootstrap, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att, extended2 report_common_obj);
extern "C" void PlFm_XZ_YZ_intern_DC(binary3 R_data, int R_nR, int R_nO, int R_nA, int R_nF, int R_nT, binary2 R_patS, 
		  long double R_emcrit,extended3 R_ro_o, extended3 R_ta_o, extended1 R_ga_o,
		  extended3 R_ro_n, extended3 R_ta_n, extended1 R_ga_n, extended3 R_ro_update, extended3 R_ta_update, bool flag, double R_delta,
		 extended3 margprobx_o,extended3 margprobx_n,extended3 somega,extended3 gradient_ro,extended3  se_ro,extended3  p_o_r_t,
		 extended2 probmat,extended3 condprobx_o,extended3 condprobx_n,extended3 gradient_ta,extended3 se_ta,extended2 p_r_t,
		 extended1 gradient_ga,extended1 p_r,extended1 se_ga,
		 int R_Nbootstrap, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att, extended2 report_common_obj);
extern "C" void PlFm_XZ_Y_intern_ADD(binary3  R_data , int R_nR , int R_nO , int R_nA , int R_nF, int R_nT, binary2 R_patS,							
		  long double R_emcrit, extended3 R_ro_o, extended2 R_ta_o, extended1 R_ga_o, 
          extended3 R_ro_n, extended2 R_ta_n, extended1 R_ga_n, extended3 R_ro_update, extended2 R_ta_update, bool flag, double R_delta,
		  extended3 margprobx_o,extended3 margprobx_n,extended3 somega,extended3 gradient_ro,extended3 se_ro,extended3 p_o_r_t,
		  extended2 probmat,extended2 condprobx_o,extended2 condprobx_n,extended2 gradient_ta,extended2 se_ta,extended2  p_r_t,
		 extended1 gradient_ga,extended1 p_r,extended1 se_ga, 
		 int R_Nbootstrap, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att, extended2 report_common_obj);
extern "C" void PlFm_X_YZ_intern_ADD(binary3 R_data, int R_nR, int R_nO, int R_nA, int R_nF, int R_nT, binary2 R_patS, 
		  long double R_emcrit,extended2 R_ro_o, extended3 R_ta_o, extended1 R_ga_o,
		  extended2 R_ro_n, extended3 R_ta_n, extended1 R_ga_n, extended2 R_ro_update, extended3 R_ta_update, bool flag, double R_delta,
		 extended2 margprobx_o,extended2 margprobx_n,extended3 somega,extended2 gradient_ro,extended2  se_ro,extended3  p_o_r_t,
		 extended2 probmat,extended3 condprobx_o,extended3 condprobx_n,extended3 gradient_ta,extended3 se_ta,extended2 p_r_t,
		 extended1 gradient_ga,extended1 p_r,extended1 se_ga,
		 int R_Nbootstrap, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att, extended2 report_common_obj);
extern "C" void PlFm_XZ_YZ_intern_ADD(binary3 R_data, int R_nR, int R_nO, int R_nA, int R_nF, int R_nT, binary2 R_patS, 
		  long double R_emcrit,extended3 R_ro_o, extended3 R_ta_o, extended1 R_ga_o,
		  extended3 R_ro_n, extended3 R_ta_n, extended1 R_ga_n, extended3 R_ro_update, extended3 R_ta_update, bool flag, double R_delta,
		 extended3 margprobx_o,extended3 margprobx_n,extended3 somega,extended3 gradient_ro,extended3  se_ro,extended3  p_o_r_t,
		 extended2 probmat,extended3 condprobx_o,extended3 condprobx_n,extended3 gradient_ta,extended3 se_ta,extended2 p_r_t,
		 extended1 gradient_ga,extended1 p_r,extended1 se_ga,
		 int R_Nbootstrap, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att, extended2 report_common_obj);
void CreateVariables(extended3 ,extended2 ,extended2 ,extended3 ,extended3 ,extended2 , extended3  ,extended2 ,extended1 );
void CreateVariables( extended3 , extended3 , extended3 , extended2 , extended2 , extended2 , extended3 , extended2 , extended1 );
void CreateVariables( extended3 , extended3 , extended3 , extended3 , extended3 , extended2 , extended3 , extended2 , extended1 );
	
void C_destructor(extended3 ,extended2 ,extended2 ,extended3 ,extended3 ,extended2 , extended3  ,extended2 ,extended1);
void C_destructor(extended3 ,extended3 ,extended3 ,extended2 ,extended2 ,extended2 , extended3  ,extended2 ,extended1);
void C_destructor(extended3 ,extended3 ,extended3 ,extended3 ,extended3 ,extended2 , extended3  ,extended2 ,extended1);

extern "C" void R2CMapper(int * , int * , int * , int * , int * , int * , int * , 
						  double * , double * , double * , double * ,
						  double * , double * , double * );
void calculate_probmat_gradient(extended1 , extended2 , extended3,extended3  ,extended2 ,extended1 ,extended2  );
void calculate_probmat_gradient(extended1 , extended3 , extended2,extended3  ,extended2 ,extended1 ,extended2  );
void calculate_probmat_gradient(extended1 , extended3 , extended3,extended3  ,extended2 ,extended1 ,extended2  );
void calculate_gradient_ro(extended3 , extended1 , extended3 , extended2,extended3,extended3  ,extended2 ,extended1, extended3 );
void calculate_gradient_ro(extended2 , extended1 , extended2 , extended3,extended2,extended3  ,extended2 ,extended1, extended2 );
void calculate_gradient_ro(extended3 , extended1 , extended3 , extended3,extended3,extended3  ,extended2 ,extended1, extended3 );
long double element_gradient_ro(unsigned short int , unsigned short int , unsigned short int ,  extended3 , extended1 ,extended3 , extended2 ,extended3  ,extended2 ,extended1  );
long double element_gradient_ro(unsigned short int , unsigned short int ,  extended2, extended1 ,extended2 , extended3 ,extended3  ,extended2 , extended1 );
long double element_gradient_ro(unsigned short int , unsigned short int , unsigned short int , extended3, extended1 ,extended3 , extended3 ,extended3  ,extended2 ,extended1  );
void calculate_gradient_ta_DC(extended2 , extended1 , extended2 , extended3,extended2 ,extended3  ,extended2 ,extended1, extended2  );
void calculate_gradient_ta_DC(extended3 , extended1 , extended3 , extended2,extended3 ,extended3  ,extended2 ,extended1, extended3  );
void calculate_gradient_ta_DC(extended3 , extended1 , extended3 , extended3,extended3 ,extended3  ,extended2 ,extended1, extended3  );
void calculate_gradient_ta_ADD(extended2 , extended1 , extended2 , extended3,extended2 ,extended3  ,extended2 ,extended1, extended2  );
void calculate_gradient_ta_ADD(extended3 , extended1 , extended3 , extended2,extended3 ,extended3  ,extended2 ,extended1, extended3  );
void calculate_gradient_ta_ADD(extended3 , extended1 , extended3 , extended3,extended3 ,extended3  ,extended2 ,extended1, extended3  );
long double element_gradient_ta_DC(unsigned short int , unsigned short int , extended2 , extended1 , extended2 , extended3,extended3  ,extended2 ,extended1   );
long double element_gradient_ta_DC(unsigned short int , unsigned short int ,unsigned short int, extended3 , extended1 , extended3 , extended2,extended3  ,extended2 ,extended1   );
long double element_gradient_ta_DC(unsigned short int , unsigned short int ,unsigned short int, extended3 , extended1 , extended3 , extended3,extended3  ,extended2 ,extended1   );
long double element_gradient_ta_ADD(unsigned short int , unsigned short int , extended2 , extended1 , extended2 , extended3,extended3  ,extended2 ,extended1   );
long double element_gradient_ta_ADD(unsigned short int , unsigned short int ,unsigned short int, extended3 , extended1 , extended3 , extended2,extended3  ,extended2 ,extended1   );
long double element_gradient_ta_ADD(unsigned short int , unsigned short int ,unsigned short int, extended3 , extended1 , extended3 , extended3,extended3  ,extended2 ,extended1   );
void calculate_gradient_ga(extended1,extended1,extended3  ,extended2 ,extended1    );
long double element_gradient_ga(unsigned short int , extended1,extended3  ,extended2 ,extended1   );
void calculate_se_ro_DC(extended3 ,extended2 ,extended1 ,extended2 ,extended3 ,extended3 ,extended3  ,extended2 ,extended1 ,extended2, extended3 );
void calculate_se_ro_DC(extended2 ,extended3 ,extended1 ,extended3 ,extended2 ,extended2 ,extended3  ,extended2 ,extended1 ,extended2, extended2 );
void calculate_se_ro_DC(extended3 ,extended3 ,extended1 ,extended3 ,extended3 ,extended3 ,extended3  ,extended2 ,extended1 ,extended2, extended3 );
void calculate_se_ta_DC(extended3 ,extended2 ,extended1 ,extended2 ,extended3 ,extended2 ,extended3  ,extended2 ,extended1 ,extended2, extended2 );
void calculate_se_ta_DC(extended2 ,extended3 ,extended1 ,extended3 ,extended2 ,extended3 ,extended3  ,extended2 ,extended1 ,extended2, extended3 );
void calculate_se_ta_DC(extended3 ,extended3 ,extended1 ,extended3 ,extended3 ,extended3 ,extended3  ,extended2 ,extended1 ,extended2, extended3 );
void calculate_se_ga_DC(extended3 ,extended2 ,extended1 ,extended2 ,extended3 ,extended1 ,extended3  ,extended2 ,extended1 ,extended2 );
void calculate_se_ga_DC(extended2 ,extended3 ,extended1 ,extended3 ,extended2 ,extended1 ,extended3  ,extended2 ,extended1 ,extended2 );
void calculate_se_ga_DC(extended3 ,extended3 ,extended1 ,extended3 ,extended3 ,extended1 ,extended3  ,extended2 ,extended1 ,extended2 );
void calculate_se_ro_ADD(extended3 ,extended2 ,extended1 ,extended2 ,extended3 ,extended3 ,extended3  ,extended2 ,extended1 ,extended2, extended3 );
void calculate_se_ro_ADD(extended2 ,extended3 ,extended1 ,extended3 ,extended2 ,extended2 ,extended3  ,extended2 ,extended1 ,extended2, extended2 );
void calculate_se_ro_ADD(extended3 ,extended3 ,extended1 ,extended3 ,extended3 ,extended3 ,extended3  ,extended2 ,extended1 ,extended2, extended3 );
void calculate_se_ta_ADD(extended3 ,extended2 ,extended1 ,extended2 ,extended3 ,extended2 ,extended3  ,extended2 ,extended1 ,extended2, extended2 );
void calculate_se_ta_ADD(extended2 ,extended3 ,extended1 ,extended3 ,extended2 ,extended3 ,extended3  ,extended2 ,extended1 ,extended2, extended3 );
void calculate_se_ta_ADD(extended3 ,extended3 ,extended1 ,extended3 ,extended3 ,extended3 ,extended3  ,extended2 ,extended1 ,extended2, extended3 );
void calculate_se_ga_ADD(extended3 ,extended2 ,extended1 ,extended2 ,extended3 ,extended1 ,extended3  ,extended2 ,extended1 ,extended2 );
void calculate_se_ga_ADD(extended2 ,extended3 ,extended1 ,extended3 ,extended2 ,extended1 ,extended3  ,extended2 ,extended1 ,extended2 );
void calculate_se_ga_ADD(extended3 ,extended3 ,extended1 ,extended3 ,extended3 ,extended1 ,extended3  ,extended2 ,extended1 ,extended2 );


// bootstrapping functions
void bootstrap_dependency_XZ_Y_DC(binary3, extended3,extended2,extended1, extended3, extended3, extended2, extended2);
void bootstrap_dependency_X_YZ_DC(binary3, extended2,extended3,extended1, extended3, extended3, extended2, extended2);
void bootstrap_dependency_XZ_YZ_DC(binary3, extended3,extended3,extended1, extended3, extended3, extended2, extended2);
void bootstrap_dependency_XZ_Y_ADD(binary3, extended3,extended2,extended1, extended3, extended3, extended2, extended2);
void bootstrap_dependency_X_YZ_ADD(binary3, extended2,extended3,extended1, extended3, extended3, extended2, extended2);
void bootstrap_dependency_XZ_YZ_ADD(binary3, extended3,extended3,extended1, extended3, extended3, extended2, extended2);
void compute_OR_common_attribute_obs(binary3, extended2);
void compute_OR_common_object_obs(binary3, extended2);
long double or_att_pair(unsigned short int,unsigned short int,unsigned short int, binary3);
long double or_obj_pair(unsigned short int,unsigned short int,unsigned short int, binary3);
void generate_data_fixedobject_XZ_Y_DC(extended1, extended3, extended2, binary3);
void generate_data_fixedobject_X_YZ_DC(extended1, extended2, extended3, binary3);
void generate_data_fixedobject_XZ_YZ_DC(extended1, extended3, extended3, binary3);
void generate_data_fixedobject_XZ_Y_ADD(extended1, extended3, extended2, binary3);
void generate_data_fixedobject_X_YZ_ADD(extended1, extended2, extended3, binary3);
void generate_data_fixedobject_XZ_YZ_ADD(extended1, extended3, extended3, binary3);
void compute_OR_common_attribute_rep(binary3, unsigned short int, extended3);
void compute_OR_common_object_rep(binary3, unsigned short int, extended3);
unsigned short int trek_latente_klasse( unsigned short int, extended1, extended1);
void compute_mean_OR( extended3, unsigned short int, unsigned short int, unsigned short int, extended2 );
void create_report_common_att(extended2, extended2, extended2);
void create_report_common_obj(extended2, extended2, extended2);



// C++  functions to generate random numbers
double runif(int , int , int );
double rstnorm(void);
double rgamma_best(double );
double rbeta(double , double );

long double loglikelihood (extended1, extended2,extended3);
long double loglikelihood (extended1, extended3,extended2);
long double loglikelihood (extended1, extended3,extended3);

// R functions
void Create_ro_ta(extended3, extended2);
void Create_ro_ta(extended2, extended3);
void Create_ro_ta(extended3, extended3);

void generate_true_parameters(extended3 , extended2 , extended1 );
void generate_true_parameters(extended2 , extended3 , extended1 );
void generate_true_parameters(extended3 , extended3 , extended1 );
void R_destructor(binary3,binary2,extended3,extended2,extended1,extended3,extended2,extended1);





int main(void)
{
	return 0;
}


// function to estimate disjunctive LCPLFM with constant classification of objects and heterogeneous object parameters (model M_1)
extern "C" void PlFm_XZ_Y_intern_DC(binary3 R_data, int R_nR, int R_nO, int R_nA, int R_nF, int R_nT, binary2 R_patS, 
		  long double R_emcrit,extended3 R_ro_o, extended2 R_ta_o, extended1 R_ga_o,
		  extended3 R_ro_n, extended2 R_ta_n, extended1 R_ga_n, extended3 R_ro_update, extended2 R_ta_update, bool flag, double R_delta,
		 extended3 margprobx_o,extended3 margprobx_n,extended3 somega,extended3 gradient_ro,extended3  se_ro,extended3  p_o_r_t,
		 extended2 probmat,extended2 condprobx_o,extended2 condprobx_n,extended2 gradient_ta,extended2 se_ta,extended2 p_r_t,
		 extended1 gradient_ga,extended1 p_r,extended1 se_ga,
		 int R_Nbootstrap, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att, extended2 report_common_obj)
{


	extended3 ro_o; extended3 ro_n; extended3 ro_update;
	extended2 ta_o; extended2 ta_n; extended2 ta_update;
	extended1 ga_o; extended1 ga_n;
	

	// pass all data from R to C
	dataB = R_data;
	_nR = R_nR;
	_nO = R_nO;
	_nA = R_nA;
	_nF = R_nF;
	_nT = R_nT;
	_Nbootstrap = R_Nbootstrap;
	patS = R_patS;
	emcrit = R_emcrit;
	ro_o = R_ro_o;
	ta_o = R_ta_o;
	ga_o = R_ga_o;
	ro_n = R_ro_n;
	ta_n = R_ta_n;
	ga_n = R_ga_n;
	ro_update=R_ro_update;
	ta_update=R_ta_update;
	
	delta = R_delta;

	///////////////////// main program /////////////////////
	
    //compute conditional probability of observing an object-attribute association given object classification using old parameter values
	calculate_condprobX_DC(ta_o,condprobx_o); 
	// compute marginal probability of object classification given person classification using old parameter values
	calculate_margprobX(ro_o,margprobx_o); 
	// compute posterior density using old parameter values
	logold= logposterior(ro_o,ta_o,ga_o,condprobx_o,margprobx_o,ro_update,ta_update);  

    //compute conditional probability of observing an object-attribute association given object classification using new parameter values
	calculate_condprobX_DC(ta_n,condprobx_n);
    // compute marginal probability of object classification given person classification using new parameter values
	calculate_margprobX(ro_n,margprobx_n);
	// compute posterior density using new parameter values
	lognew= logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n,ro_update,ta_update);

	// compute value of convergence criterion
    calculate_conv( ro_o, ro_n, ta_o, ta_n, ga_o, ga_n, ro_update, ta_update);

	while( conv>emcrit)
	{ 
		// replace old parameter values by new parameter values
        replace( ro_o, ro_n, ta_o, ta_n, ga_o, ga_n,ro_update,ta_update, condprobx_o, condprobx_n, margprobx_o, margprobx_n);
		// compute posterior probabilities of classification person
        calculate_probmat(ga_o,condprobx_o,margprobx_o,probmat);
		// update the class size parameters
        emgamma(ga_n, probmat);
		// compute probability of object classification given data and classification of person 
        calculate_condprob_pattern(condprobx_o,margprobx_o,somega);
        // update object parameters
		update_emro(ro_n,probmat,ro_update);
		// update attribute parameters
		update_emta_DC(condprobx_o,ta_o, ta_n,probmat,ta_update);
		// compute conditional probability of data given object classification and given person classification using new parameter values
		calculate_condprobX_DC(ta_n,condprobx_n);
	    // compute marginal probability of object classification given person classification using new parameter values
		calculate_margprobX(ro_n,margprobx_n);
		// compute posterior density using new parameter values
        lognew= logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n,ro_update,ta_update);
       	// compute value of convergence criterion
		calculate_conv(ro_o, ro_n, ta_o, ta_n, ga_o, ga_n,ro_update,ta_update);
	} 

    // compute posterior density using new parameter values
	logpost_n = logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n,ro_update,ta_update);
	// compute log likelihood using new parameter values
	loglik_n  = loglikelihood (ga_n, condprobx_n, margprobx_n);
    
    // if flag==0 computation of gradient and standard errors can be skipped
    // if flag==1 computation of gradient and standard errors is needed
	if (flag == true)
	{
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t , p_r_t, p_r, probmat);
		calculate_gradient_ro(ro_n, ga_n, margprobx_n, condprobx_n, gradient_ro,p_o_r_t , p_r_t, p_r,ro_update);
		calculate_gradient_ta_DC(ta_n, ga_n,condprobx_n, margprobx_n, gradient_ta,p_o_r_t , p_r_t, p_r,ta_update);
		calculate_gradient_ga(ga_n, gradient_ga,p_o_r_t , p_r_t, p_r);

		calculate_se_ro_DC( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ro, p_o_r_t , p_r_t, p_r, probmat, ro_update);
		calculate_se_ta_DC( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ta, p_o_r_t , p_r_t, p_r, probmat, ta_update);
		calculate_se_ga_DC( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ga, p_o_r_t , p_r_t, p_r, probmat);

		calculate_probmat(ga_n,condprobx_n,margprobx_n,probmat);
	
		// bootstrap procedure
		bootstrap_dependency_XZ_Y_DC(dataB, ro_n, ta_n, ga_n, OR_common_att_rep, OR_common_obj_rep, report_common_att, report_common_obj);

	}

	C_destructor( somega, condprobx_o, condprobx_n, margprobx_o, margprobx_n, probmat,   p_o_r_t, p_r_t, p_r);
}

// function to estimate disjunctive LCPLFM with fixed object classification and heterogeneous attribute parameters (Model M_2)
extern "C" void PlFm_X_YZ_intern_DC(binary3 R_data, int R_nR, int R_nO, int R_nA, int R_nF, int R_nT, binary2 R_patS, 
		  long double R_emcrit,extended2 R_ro_o, extended3 R_ta_o, extended1 R_ga_o,
		  extended2 R_ro_n, extended3 R_ta_n, extended1 R_ga_n,extended2 R_ro_update, extended3 R_ta_update, bool flag, double R_delta,
		 extended2 margprobx_o,extended2 margprobx_n,extended3 somega,extended2 gradient_ro,extended2  se_ro,extended3  p_o_r_t,
		 extended2 probmat,extended3 condprobx_o,extended3 condprobx_n,extended3 gradient_ta,extended3 se_ta,extended2 p_r_t,
		 extended1 gradient_ga,extended1 p_r,extended1 se_ga,
		  int R_Nbootstrap, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att, extended2 report_common_obj)
{

	extended2 ro_o; extended2 ro_n; extended2 ro_update;
	extended3 ta_o; extended3 ta_n; extended3 ta_update;
	extended1 ga_o; extended1 ga_n;

	// pass all data from R to C
	dataB = R_data;
	_nR = R_nR;
	_nO = R_nO;
	_nA = R_nA;
	_nF = R_nF;
	_nT = R_nT;
	_Nbootstrap = R_Nbootstrap;
	patS = R_patS;
	emcrit = R_emcrit;
	ro_o = R_ro_o;
	ta_o = R_ta_o;
	ga_o = R_ga_o;
	ro_n = R_ro_n;
	ta_n = R_ta_n;
	ga_n = R_ga_n;
	ro_update=R_ro_update;
	ta_update=R_ta_update;
	delta = R_delta;

	///////////////////// main program /////////////////////
	
    //compute conditional probability of observing an object-attribute association given object and person classification using old parameter values
	calculate_condprobX_DC(ta_o,condprobx_o);
	// compute marginal probability of object classification using old parameter values
	calculate_margprobX(ro_o,margprobx_o);
	// compute posterior density using old parameter values
    logold= logposterior(ro_o,ta_o,ga_o,condprobx_o,margprobx_o,ro_update,ta_update);

    //compute conditional probability of observing an object-attribute association given object and person classification using new parameter values
	calculate_condprobX_DC(ta_n,condprobx_n);
	// compute marginal probability of object classification using new parameter values
	calculate_margprobX(ro_n,margprobx_n);
	// compute posterior density using new parameter values
    lognew= logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n,ro_update,ta_update);
    
    // compute value of convergence criterion
	calculate_conv( ro_o, ro_n, ta_o, ta_n, ga_o, ga_n, ro_update, ta_update);


	while( conv>emcrit)
	{ 
        // replace old parameter values by new parameter values
		replace( ro_o, ro_n, ta_o, ta_n, ga_o, ga_n, ro_update,ta_update, condprobx_o, condprobx_n, margprobx_o, margprobx_n);
		// compute posterior probabilities of classification person
		calculate_probmat(ga_o,condprobx_o,margprobx_o,probmat);
		// update the class size parameters
		emgamma(ga_n, probmat);
		// compute probability of object classification given data and classification of person 
		calculate_condprob_pattern(condprobx_o,margprobx_o,somega);
		// update object parameters
		update_emro(ro_n,probmat, ro_update);
		// update attribute parameters
		update_emta_DC(condprobx_o,ta_o, ta_n,probmat, ta_update);
		// compute conditional probability of data given object classification and given person classification using new parameter values
		calculate_condprobX_DC(ta_n,condprobx_n);
		// compute marginal probability of object classification given person classification using new parameter values
		calculate_margprobX(ro_n,margprobx_n);
		// compute posterior density using new parameter values
		lognew= logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n, ro_update, ta_update);
		// compute value of convergence criterion
		calculate_conv(ro_o, ro_n, ta_o, ta_n, ga_o, ga_n, ro_update, ta_update);

	} 
     // compute posterior density using new parameter values
	logpost_n = logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n, ro_update, ta_update);
	// compute loglikelihood using new parameter values
	loglik_n  = loglikelihood (ga_n, condprobx_n, margprobx_n);

    // if flag==0 computation of gradient and standard errors can be skipped
    // if flag==1 computation of gradient and standard errors is needed
	if (flag == true)
	{
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t , p_r_t, p_r, probmat);
		calculate_gradient_ro(ro_n, ga_n, margprobx_n, condprobx_n, gradient_ro,p_o_r_t , p_r_t, p_r, ro_update);
		calculate_gradient_ta_DC(ta_n, ga_n,condprobx_n, margprobx_n, gradient_ta,p_o_r_t , p_r_t, p_r, ta_update);
		calculate_gradient_ga(ga_n, gradient_ga,p_o_r_t , p_r_t, p_r);

		calculate_se_ro_DC( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ro, p_o_r_t , p_r_t, p_r, probmat, ro_update);
		calculate_se_ta_DC( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ta, p_o_r_t , p_r_t, p_r, probmat, ta_update);
		calculate_se_ga_DC( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ga, p_o_r_t , p_r_t, p_r, probmat);

		calculate_probmat(ga_n,condprobx_n,margprobx_n,probmat);
		
		// bootstrap procedure
		bootstrap_dependency_X_YZ_DC(dataB, ro_n, ta_n, ga_n, OR_common_att_rep, OR_common_obj_rep, report_common_att, report_common_obj);

	}

	C_destructor( somega, condprobx_o, condprobx_n, margprobx_o, margprobx_n, probmat,   p_o_r_t, p_r_t, p_r);
}

// function to estimate disjunctive LCPLFM with fixed object classification and heterogeneous object- and attribute parameters (model M_3)
extern "C" void PlFm_XZ_YZ_intern_DC(binary3 R_data, int R_nR, int R_nO, int R_nA, int R_nF, int R_nT, binary2 R_patS, 
		  long double R_emcrit,extended3 R_ro_o, extended3 R_ta_o, extended1 R_ga_o,
		  extended3 R_ro_n, extended3 R_ta_n, extended1 R_ga_n, extended3 R_ro_update, extended3 R_ta_update, bool flag, double R_delta,
		 extended3 margprobx_o,extended3 margprobx_n,extended3 somega,extended3 gradient_ro,extended3  se_ro,extended3  p_o_r_t,
		 extended2 probmat,extended3 condprobx_o,extended3 condprobx_n,extended3 gradient_ta,extended3 se_ta,extended2 p_r_t,
		 extended1 gradient_ga,extended1 p_r,extended1 se_ga,
		 int R_Nbootstrap, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att, extended2 report_common_obj)
{

	extended3 ro_o; extended3 ro_n; extended3 ro_update;
	extended3 ta_o; extended3 ta_n; extended3 ta_update;
	extended1 ga_o; extended1 ga_n;

	// pass all data from R to C
	dataB = R_data;
	_nR = R_nR;
	_nO = R_nO;
	_nA = R_nA;
	_nF = R_nF;
	_nT = R_nT;
	patS = R_patS;
	emcrit = R_emcrit;
	ro_o = R_ro_o;
	ta_o = R_ta_o;
	ga_o = R_ga_o;
	ro_n = R_ro_n;
	ta_n = R_ta_n;
	ga_n = R_ga_n;
	ro_update=R_ro_update;
	ta_update=R_ta_update;
	delta = R_delta;

	///////////////////// main program /////////////////////
	
    //compute conditional probability of observing an object-attribute association given object and person classification using old parameter values
	calculate_condprobX_DC(ta_o,condprobx_o);
	// compute marginal probability of object classification given person classification using old parameter values
	calculate_margprobX(ro_o,margprobx_o);
	// compute posterior density using old parameter values
    logold= logposterior(ro_o,ta_o,ga_o,condprobx_o,margprobx_o,ro_update,ta_update);

    //compute conditional probability of observing an object-attribute association given object and person classification using new parameter values
	calculate_condprobX_DC(ta_n,condprobx_n);
	// compute marginal probability of object classification given person classification using new parameter values
	calculate_margprobX(ro_n,margprobx_n);
	// compute posterior density using new parameter values
    lognew= logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n,ro_update,ta_update);

    // compute value of convergence criterion
	calculate_conv( ro_o, ro_n, ta_o, ta_n, ga_o, ga_n,ro_update,ta_update);


	while( conv>emcrit)
	{ 
		// replace old parameter values by new parameter values
        replace( ro_o, ro_n, ta_o, ta_n, ga_o, ga_n, ro_update,ta_update,  condprobx_o, condprobx_n, margprobx_o, margprobx_n);
        // compute posterior probabilities of classification person
		calculate_probmat(ga_o,condprobx_o,margprobx_o,probmat);
		// update the class size parameters
		emgamma(ga_n, probmat);
		// compute probability of object classification given data and classification of person 
		calculate_condprob_pattern(condprobx_o,margprobx_o,somega);
		// update object parameters
		update_emro(ro_n,probmat,ro_update);
		// update attribute parameters
		update_emta_DC(condprobx_o,ta_o, ta_n,probmat,ta_update);
		// compute conditional probability of data given object classification and given person classification using new parameter values
		calculate_condprobX_DC(ta_n,condprobx_n);
		// compute marginal probability of object classification given person classification using new parameter values
		calculate_margprobX(ro_n,margprobx_n);
		// compute posterior density using new parameter values
		lognew= logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n,ro_update,ta_update);
		// compute value of convergence criterion
		calculate_conv(ro_o, ro_n, ta_o, ta_n, ga_o, ga_n,ro_update,ta_update);

	} 
    // compute posterior density using new parameter values
	logpost_n = logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n,ro_update,ta_update);
	// compute loglikelihood using new parameter values
	loglik_n  = loglikelihood (ga_n, condprobx_n, margprobx_n);

    // if flag==0 computation of gradient and standard errors can be skipped
    // if flag==1 computation of gradient and standard errors is needed
    
	if (flag == true)
	{
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t , p_r_t, p_r, probmat);
		calculate_gradient_ro(ro_n, ga_n, margprobx_n, condprobx_n, gradient_ro,p_o_r_t , p_r_t, p_r,ro_update);
		calculate_gradient_ta_DC(ta_n, ga_n,condprobx_n, margprobx_n, gradient_ta,p_o_r_t , p_r_t, p_r,ta_update);
		calculate_gradient_ga(ga_n, gradient_ga,p_o_r_t , p_r_t, p_r);

		calculate_se_ro_DC( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ro, p_o_r_t , p_r_t, p_r, probmat,ro_update);
		calculate_se_ta_DC( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ta, p_o_r_t , p_r_t, p_r, probmat,ta_update);
		calculate_se_ga_DC( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ga, p_o_r_t , p_r_t, p_r, probmat);

		calculate_probmat(ga_n,condprobx_n,margprobx_n,probmat);
		
		// bootstrap procedure
		bootstrap_dependency_XZ_YZ_DC(dataB, ro_n, ta_n, ga_n, OR_common_att_rep, OR_common_obj_rep, report_common_att, report_common_obj);

	}

	C_destructor( somega, condprobx_o, condprobx_n, margprobx_o, margprobx_n, probmat,   p_o_r_t, p_r_t, p_r);
}




// estimate disjunctive model M_1
extern "C" void PlFm_XZ_Y_DC(int * _R_data, int * _R_nR, int * _R_nO, int * _R_nA, int * _R_nF, int * _R_nT, int * _R_patS, 
						  double * _R_emcrit,
						  double * _R_ro_n, double * _R_ta_n, double * _R_ga_n,
						  double * _R_ro_update, double * _R_ta_update,
						  double * _R_gradient_ro, double * _R_gradient_ta, double * _R_gradient_ga,
						  double * _R_se_ro, double * _R_se_ta, double * _R_se_ga, double * _R_delta, double * _R_logpost_n, double * _R_loglik_n, double * _R_probmat, int * _R_flag,
						  int * _R_Nbootstrap, double * _R_OR_common_att_rep, double * _R_OR_common_obj_rep, double * _R_report_common_att, double * _R_report_common_obj)
{
	// variables
	
	extended3 margprobx_o,margprobx_n,somega,gradient_ro, se_ro, p_o_r_t;
	extended2 probmat,condprobx_o,condprobx_n,gradient_ta, se_ta, p_r_t;
	extended1 gradient_ga,p_r,se_ga;


	binary3 R_data; int R_nR; int R_nO; int R_nA; int R_nF; int R_nT; int R_Nbootstrap; binary2 R_patS; 
	long double R_emcrit;extended3 R_ro_o; extended2 R_ta_o; extended1 R_ga_o;
	extended3 R_ro_n; extended2 R_ta_n; extended1 R_ga_n; extended3 R_ro_update; extended2 R_ta_update;
	

	int R_flag;
	double R_delta;

	R_nR = *_R_nR; 
	R_nO = *_R_nO; 
	R_nA = *_R_nA; 
	R_nF = *_R_nF; 
	R_nT = *_R_nT;
	R_Nbootstrap = *_R_Nbootstrap;

	R_flag = *_R_flag;
	R_delta = *_R_delta;
	
	bool flag;
	if (R_flag == 0)
		flag = false;
	else
		flag = true;

	_nO = R_nO;
        _nA = R_nA;
        _nR = R_nR;
	_nT = R_nT;
	_nF = R_nF;

	R_data = new bool **[R_nO+1]; 
	for (int i = 0; i < R_nO+1; i++) 
	{
		R_data[i] = new bool *[R_nA+1];
		for (int ii = 0; ii < R_nA+1; ii++)
		{
			R_data[i][ii] = new bool [R_nR+1];
		}
	}
	for (int i = 0; i < R_nO; i++)
	{
		for (int j = 0; j < R_nA; j++)
		{
			for (int k = 0; k < R_nR; k++)
			{
				R_data[i+1][j+1][k+1] = (bool)_R_data[i * (R_nA) * (R_nR) + j * (R_nR) + k];
			}
		}
	}

	calculate_dims();

	somega = new long double ** [_nR+1]; 
	condprobx_o = new long double * [_nS+1]; 
	margprobx_o = new long double ** [_nS+1]; 
	condprobx_n = new long double * [_nS+1]; 
	margprobx_n = new long double ** [_nS+1]; 
	probmat = new long double * [_nT+1]; 
	p_o_r_t = new long double ** [_nO+1]; 
	p_r_t = new long double * [_nR+1];
	p_r = new long double [_nR+1];
	CreateVariables( somega, condprobx_o, condprobx_n, margprobx_o, margprobx_n, probmat, p_o_r_t, p_r_t, p_r);

	R_patS = new bool * [_nS+1]; 
	for (int i = 0; i < _nS+1; i++)
	{
		R_patS[i] = new bool [R_nF+1];
	}
	for (int i = 0; i < _nS; i++)
	{
		for (int j = 0; j < R_nF; j++)
		{
			R_patS[i+1][j+1] = (bool)_R_patS[i * (R_nF) + j ];			
		}
	}
	R_emcrit = *_R_emcrit; 
	R_ro_o = new long double ** [_nO+1] ; R_ro_n = new long double ** [_nO+1]; R_ro_update = new long double ** [_nO+1];
	R_ta_o = new long double * [_nA+1] ; R_ta_n = new long double * [_nA+1] ; R_ta_update = new long double * [_nA+1] ;
	R_ga_o = new long double [_nT+1]; R_ga_n = new long double [_nT+1];
	Create_ro_ta(R_ro_o,R_ta_o);
	Create_ro_ta(R_ro_n,R_ta_n);
	Create_ro_ta(R_ro_update,R_ta_update);
	seed1= 0; 
	seed2= 1; 
	seed3= 2; 
	generate_true_parameters(R_ro_o,R_ta_o,R_ga_o);

	for (int i = 0; i < _nO; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			for (int k = 0; k < _nT; k++)
			{

				R_ro_n[i+1][j+1][k+1] =  _R_ro_n[  i * (_nF) * (_nT) + j * (_nT) + k];
			}
		}
	}

    for (int i = 0; i < _nO; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			for (int k = 0; k < _nT; k++)
			{

				R_ro_update[i+1][j+1][k+1] =  _R_ro_update[  i * (_nF) * (_nT) + j * (_nT) + k];
			}
		}
	}
	
    for (int i = 0; i < _nA; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
		
			R_ta_n[i+1][j+1] =  _R_ta_n[  i * (_nF) + j];	
		}
	}

    for (int i = 0; i < _nA; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
		
			R_ta_update[i+1][j+1] =  _R_ta_update[  i * (_nF) + j];	
		}
	}

	for (int i = 0; i < _nT; i++)
	{
	
		R_ga_n[i+1] =  _R_ga_n[i];	
	}

	// bootstrapping variables
	extended3 OR_common_att_rep, OR_common_obj_rep;
	extended2 report_common_att, report_common_obj;
	int N_attpair, N_objpair;
	N_attpair = (_nA*(_nA-1)) / 2;
	N_objpair = (_nO*(_nO-1)) / 2;


	if (flag == true)
	{

		gradient_ta = new long double * [_nA+1]; 
		for (int i = 0; i < _nA+1; i++)
		{
			gradient_ta[i] = new long double [_nF+1];
		}


		gradient_ro = new long double ** [_nO+1]; 
		for (int i = 0; i < _nO+1; i++)
		{
			gradient_ro[i] = new long double * [_nF+1];
			for (int ii = 0; ii < _nF+1; ii++)
			{
				gradient_ro[i][ii] = new long double [_nT+1];
			}
		}


		gradient_ga = new long double [_nT+1];



		se_ta = new long double * [_nA+1]; 
		for (int i = 0; i < _nA+1; i++)
		{
			se_ta[i] = new long double [_nF+1];
		}


		se_ro = new long double ** [_nO+1]; 
		for (int i = 0; i < _nO+1; i++)
		{
			se_ro[i] = new long double * [_nF+1];
			for (int ii = 0; ii < _nF+1; ii++)
			{
				se_ro[i][ii] = new long double [_nT+1];
			}
		}


		se_ga = new long double [_nT+1];

		// bootstrapping variables

		

		// construct OR_common_att_rep,OR_common_obj_rep
		OR_common_att_rep = new long double ** [_Nbootstrap+1] ;
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			OR_common_att_rep[i] = new long double * [_nO+1];
			for (int ii = 0; ii < _nO+1; ii++)
			{
				OR_common_att_rep[i][ii] = new long double [N_attpair+1];
			}
		}
		OR_common_obj_rep = new long double ** [_Nbootstrap+1] ;
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			OR_common_obj_rep[i] = new long double * [_nA+1];
			for (int ii = 0; ii < _nA+1; ii++)
			{
				OR_common_obj_rep[i][ii] = new long double [N_objpair+1];
			}
		}
		//construct report_common_att & report_common_obj
		report_common_att = new long double * [_nO*N_attpair+1] ;
		for (int i = 0; i < _nO*N_attpair+1; i++)
		{
			report_common_att[i] = new long double [5+1];
		}
		report_common_obj = new long double * [_nA*N_objpair+1] ;
		for (int i = 0; i < _nA*N_objpair+1; i++)
		{
			report_common_obj[i] = new long double [5+1];
		}

	}
	else
	{
		se_ta = NULL;
		se_ga = NULL;
		se_ro = NULL;
		gradient_ta = NULL; 
		gradient_ga = NULL; 
		gradient_ro = NULL; 
		OR_common_att_rep = NULL;
		OR_common_obj_rep = NULL;
		report_common_att = NULL;
		report_common_obj = NULL;

	}

	//////////////////// call C function model M_1 ///////////////////////
	PlFm_XZ_Y_intern_DC(R_data,R_nR,R_nO,R_nA,R_nF,R_nT,R_patS,R_emcrit,R_ro_o,R_ta_o,R_ga_o,R_ro_n,R_ta_n,R_ga_n,R_ro_update,R_ta_update,flag,R_delta,
		margprobx_o,margprobx_n,somega,gradient_ro, se_ro, p_o_r_t,
		probmat,condprobx_o,condprobx_n,gradient_ta, se_ta, p_r_t,
		gradient_ga,p_r,se_ga,R_Nbootstrap, OR_common_att_rep, OR_common_obj_rep, report_common_att, report_common_obj);
	////////////////////////////////////////////////////////////



	for (int i = 1; i < _nO+1; i++)
	{
		for (int j = 1; j < _nF+1; j++)
		{
			for (int k = 1; k < _nT+1; k++)
			{

				_R_ro_n[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = R_ro_n[i][j][k];
			}
		}
	}
	for (int i = 1; i < _nA+1; i++)
	{
		for (int j = 1; j < _nF+1; j++)
		{
	
			_R_ta_n[  (i-1) * (_nF) + (j-1)] = R_ta_n[i][j];	
		}
	}
	for (int i = 1; i < _nT+1; i++)
	{

		_R_ga_n[i-1] = R_ga_n[i];		
	}

	*_R_logpost_n = logpost_n;
	*_R_loglik_n = loglik_n;


	if (flag == true)
	{
		
		for (int i = 1; i < _nA+1; i++)
		{
			for (int j = 1; j < _nF+1; j++)
			{
				_R_gradient_ta[  (i-1) * (_nF) + (j-1)] = gradient_ta[i][j];
				_R_se_ta[  (i-1) * (_nF) + (j-1)] = se_ta[i][j];
			}
		}

		
		for (int i = 1; i < _nO+1; i++)
		{
			for (int j = 1; j < _nF+1; j++)
			{
				for (int k = 1; k < _nT+1; k++)
				{
					_R_gradient_ro[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = gradient_ro[i][j][k];
					_R_se_ro[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = se_ro[i][j][k];
				}
			}
		}

		for (int i = 1; i < _nT+1; i++)
		{
			_R_gradient_ga[i-1] = gradient_ga[i];
			_R_se_ga[i-1] = se_ga[i];
		}

		for (int i = 1; i < _nT+1; i++)
		{
			for (int j = 1; j < _nR+1; j++)
			{
				_R_probmat[ (i-1) * (_nR) + (j-1)] = probmat[i][j];
			}
		}

		// bootstrapping variables
		for (int i = 1; i < _Nbootstrap+1; i++)
		{
			for (int j = 1; j < _nO+1; j++)
			{
				for (int k = 1; k < N_attpair+1; k++)
				{
					_R_OR_common_att_rep[ (i-1) * (_nO) * (N_attpair) + (j-1) * (N_attpair) + (k-1)] = OR_common_att_rep[i][j][k];
					
				}
			}
		}
		for (int i = 1; i < _Nbootstrap+1; i++)
		{
			for (int j = 1; j < _nA+1; j++)
			{
				for (int k = 1; k < N_objpair+1; k++)
				{
					_R_OR_common_obj_rep[ (i-1) * (_nA) * (N_objpair) + (j-1) * (N_objpair) + (k-1)] = OR_common_obj_rep[i][j][k];
					
				}
			}
		}
		for (int i = 1; i < _nO*N_attpair+1; i++)
		{
			for (int j = 1; j < 5+1; j++)
			{
				_R_report_common_att[  (i-1) * (5) + (j-1)] = report_common_att[i][j];
			}
		}
		for (int i = 1; i < _nA*N_objpair+1; i++)
		{
			for (int j = 1; j < 5+1; j++)
			{
				_R_report_common_obj[  (i-1) * (5) + (j-1)] = report_common_obj[i][j];
			}
		}

	}


	// destructor

	for (int i = 0; i < _nO+1; i++) 
	{
		for (int ii = 0; ii < _nA+1; ii++)
		{
			delete R_data[i][ii];
		}
		delete R_data[i];
	}
	delete[] R_data ;
	for (int i = 0; i < _nS+1; i++)
	{
		delete R_patS[i];
	}
	delete[] R_patS ;

	for (int i = 0; i < _nO+1; i++)
	{
		for (int ii = 0; ii < _nF+1; ii++)
		{
			delete R_ro_o[i][ii];
			delete R_ro_n[i][ii];
			delete R_ro_update[i][ii];
		}
		delete R_ro_o[i];
		delete R_ro_n[i];
		delete R_ro_update[i];
	}
	delete[] R_ro_o ;
	delete[] R_ro_n ;
	delete[] R_ro_update ;
	
	for (int i = 0; i < _nA+1; i++)
	{
		delete R_ta_o[i];
		delete R_ta_n[i];
		delete R_ta_update[i];
	}
	delete[] R_ta_o;
	delete[] R_ta_n;
	delete[] R_ta_update;	
	delete[] R_ga_o;
	delete[] R_ga_n;

	for (int i = 0; i < _nT+1; i++)
	{
		delete probmat[i];
	}
	delete[] probmat;

	if (flag == true)
	{
		for (int i = 0; i < _nA+1; i++)
		{
			delete gradient_ta[i];
		}
		delete[] gradient_ta; 

		for (int i = 0; i < _nO+1; i++)
		{
			for (int ii = 0; ii < _nF+1; ii++)
			{
				delete gradient_ro[i][ii];
			}
			delete gradient_ro[i];
		}
		delete[] gradient_ro;

		delete[] gradient_ga;

		for (int i = 0; i < _nA+1; i++)
		{
			delete se_ta[i];
		}
		delete[] se_ta;

		for (int i = 0; i < _nO+1; i++)
		{
			for (int ii = 0; ii < _nF+1; ii++)
			{
				delete se_ro[i][ii];
			}
			delete se_ro[i];
		}
		delete[] se_ro;

		delete[] se_ga;
		
		// OR_common_att_rep
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			for (int ii = 0; ii < _nO+1; ii++)
			{
				delete OR_common_att_rep[i][ii];
			}
			delete OR_common_att_rep[i];
		}
		delete[]  OR_common_att_rep;
		
		// OR_common_obj_rep
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			for (int ii = 0; ii < _nA+1; ii++)
			{
				delete OR_common_obj_rep[i][ii];
			}
			delete OR_common_obj_rep[i];
		}
		delete[] OR_common_obj_rep;
				
		// report_common_att
		
		for (int i = 0; i < _nO*N_attpair+1; i++)
		{
			delete report_common_att[i];
		}
		delete[] report_common_att;
		
	   // report_common_obj	
	   for (int i = 0; i < _nA*N_objpair+1; i++)
		{
			delete report_common_obj[i];
		}
		delete[] report_common_obj;		
		
	}
}

// estimate disjunctive model M_2
extern "C" void PlFm_X_YZ_DC(int * _R_data, int * _R_nR, int * _R_nO, int * _R_nA, int * _R_nF, int * _R_nT, int * _R_patS, 
						  double * _R_emcrit,
						  double * _R_ro_n, double * _R_ta_n, double * _R_ga_n,
						  double * _R_ro_update, double * _R_ta_update,
						  double * _R_gradient_ro, double * _R_gradient_ta, double * _R_gradient_ga,
						  double * _R_se_ro, double * _R_se_ta, double * _R_se_ga, double * _R_delta, double * _R_logpost_n, double * _R_loglik_n, double * _R_probmat, int * _R_flag,
						  int * _R_Nbootstrap, double * _R_OR_common_att_rep, double * _R_OR_common_obj_rep, double * _R_report_common_att, double * _R_report_common_obj)
{
	// variables
	extended3 condprobx_o,condprobx_n,somega,p_o_r_t,gradient_ta,se_ta;
	extended2 probmat, margprobx_o,margprobx_n,p_r_t,gradient_ro,se_ro;
    extended1 p_r,gradient_ga,se_ga;

	binary3 R_data; int R_nR; int R_nO; int R_nA; int R_nF; int R_nT; int R_Nbootstrap; binary2 R_patS;
	long double R_emcrit; extended2 R_ro_o; extended3 R_ta_o; extended1 R_ga_o;
	extended2 R_ro_n; extended3 R_ta_n; extended1 R_ga_n; extended2 R_ro_update; extended3 R_ta_update;

	int R_flag;
	double R_delta;

	R_nR = *_R_nR; 
	R_nO = *_R_nO; 
	R_nA = *_R_nA; 
	R_nF = *_R_nF; 
	R_nT = *_R_nT;
	R_Nbootstrap = *_R_Nbootstrap;

	R_flag = *_R_flag;
	R_delta = *_R_delta;
	
	bool flag;
	if (R_flag == 0)
		flag = false;
	else
		flag = true;

	_nO = R_nO;
        _nA = R_nA;
        _nR = R_nR;
	_nT = R_nT;
	_nF = R_nF;

	R_data = new bool **[R_nO+1]; 
	for (int i = 0; i < R_nO+1; i++) 
	{
		R_data[i] = new bool *[R_nA+1];
		for (int ii = 0; ii < R_nA+1; ii++)
		{
			R_data[i][ii] = new bool [R_nR+1];
		}
	}
	for (int i = 0; i < R_nO; i++)
	{
		for (int j = 0; j < R_nA; j++)
		{
			for (int k = 0; k < R_nR; k++)
			{
				R_data[i+1][j+1][k+1] = (bool)_R_data[i * (R_nA) * (R_nR) + j * (R_nR) + k];
			}
		}
	}

	calculate_dims();

	somega = new long double ** [_nR+1]; 
	condprobx_o = new long double ** [_nS+1]; 
	margprobx_o = new long double * [_nS+1]; 
	condprobx_n = new long double ** [_nS+1]; 
	margprobx_n = new long double * [_nS+1]; 
	probmat = new long double * [_nT+1]; 
	p_o_r_t = new long double ** [_nO+1]; 
	p_r_t = new long double * [_nR+1];
	p_r = new long double [_nR+1];
	CreateVariables( somega, condprobx_o, condprobx_n, margprobx_o, margprobx_n, probmat, p_o_r_t, p_r_t, p_r);
	

	R_patS = new bool * [_nS+1]; 
	for (int i = 0; i < _nS+1; i++)
	{
		R_patS[i] = new bool [R_nF+1];
	}
	for (int i = 0; i < _nS; i++)
	{
		for (int j = 0; j < R_nF; j++)
		{
			R_patS[i+1][j+1] = (bool)_R_patS[i * (R_nF) + j ];			
		}
	}
	R_emcrit = *_R_emcrit; 
	
	R_ro_o = new long double * [_nO+1] ; R_ro_n = new long double * [_nO+1]; R_ro_update = new long double * [_nO+1];
	R_ta_o = new long double ** [_nA+1] ; R_ta_n = new long double ** [_nA+1] ; R_ta_update = new long double ** [_nA+1] ;
	R_ga_o = new long double [_nT+1]; R_ga_n = new long double [_nT+1];
	Create_ro_ta(R_ro_o,R_ta_o);
	Create_ro_ta(R_ro_n,R_ta_n);
	Create_ro_ta(R_ro_update,R_ta_update);
	seed1= 0; 
	seed2= 1; 
	seed3= 2; 
	generate_true_parameters(R_ro_o,R_ta_o,R_ga_o);

	for (int i = 0; i < _nO; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			R_ro_n[i+1][j+1] =  _R_ro_n[  i * (_nF) + j ];
		}
	}

	for (int i = 0; i < _nO; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			R_ro_update[i+1][j+1] =  _R_ro_update[  i * (_nF) + j ];
		}
	}


	for (int i = 0; i < _nA; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			for (int k = 0; k < _nT; k++)
			{
				R_ta_n[i+1][j+1][k+1] =  _R_ta_n[  i * (_nF) * (_nT) + j * (_nT) + k];
			}

		}
	}

	for (int i = 0; i < _nA; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			for (int k = 0; k < _nT; k++)
			{
				R_ta_update[i+1][j+1][k+1] =  _R_ta_update[  i * (_nF) * (_nT) + j * (_nT) + k];
			}

		}
	}


	for (int i = 0; i < _nT; i++)
	{	
		R_ga_n[i+1] =  _R_ga_n[i];	
	}
	
		// bootstrapping variables
	extended3 OR_common_att_rep, OR_common_obj_rep;
	extended2 report_common_att, report_common_obj;
	int N_attpair, N_objpair;
	N_attpair = (_nA*(_nA-1)) / 2;
	N_objpair = (_nO*(_nO-1)) / 2;
	
	

	if (flag == true)
	{

		gradient_ta = new long double ** [_nA+1]; 
		for (int i = 0; i < _nA+1; i++)
		{
			gradient_ta[i] = new long double * [_nF+1];
			for (int ii = 0; ii < _nF+1; ii++)
			{
				gradient_ta[i][ii] = new long double [_nT+1];
			}
		}


		gradient_ro = new long double * [_nO+1]; 
		for (int i = 0; i < _nO+1; i++)
		{
			gradient_ro[i] = new long double  [_nF+1];
		}


		gradient_ga = new long double [_nT+1];



		se_ta = new long double ** [_nA+1]; 
		for (int i = 0; i < _nA+1; i++)
		{
			se_ta[i] = new long double *[_nF+1];
			for (int ii = 0; ii < _nF+1; ii++)
			{
				se_ta[i][ii] = new long double [_nT+1];
			}
		}


		se_ro = new long double * [_nO+1]; 
		for (int i = 0; i < _nO+1; i++)
		{
			se_ro[i] = new long double [_nF+1];
		}


		se_ga = new long double [_nT+1];
		
		
		
		// bootstrapping variables

		// construct OR_common_att_rep,OR_common_obj_rep
		OR_common_att_rep = new long double ** [_Nbootstrap+1] ;
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			OR_common_att_rep[i] = new long double * [_nO+1];
			for (int ii = 0; ii < _nO+1; ii++)
			{
				OR_common_att_rep[i][ii] = new long double [N_attpair+1];
			}
		}
		OR_common_obj_rep = new long double ** [_Nbootstrap+1] ;
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			OR_common_obj_rep[i] = new long double * [_nA+1];
			for (int ii = 0; ii < _nA+1; ii++)
			{
				OR_common_obj_rep[i][ii] = new long double [N_objpair+1];
			}
		}
		//construct report_common_att & report_common_obj
		report_common_att = new long double * [_nO*N_attpair+1] ;
		for (int i = 0; i < _nO*N_attpair+1; i++)
		{
			report_common_att[i] = new long double [5+1];
		}
		report_common_obj = new long double * [_nA*N_objpair+1] ;
		for (int i = 0; i < _nA*N_objpair+1; i++)
		{
			report_common_obj[i] = new long double [5+1];
		}


		
	}
		else
	{
		se_ta = NULL;
		se_ga = NULL;
		se_ro = NULL;
		gradient_ta = NULL; 
		gradient_ga = NULL; 
		gradient_ro = NULL;
		OR_common_att_rep = NULL;
		OR_common_obj_rep = NULL;
		report_common_att = NULL;
		report_common_obj = NULL; 
	}

	//////////////////// call C function model M_2///////////////////////
	PlFm_X_YZ_intern_DC(R_data,R_nR,R_nO,R_nA,R_nF,R_nT,R_patS,R_emcrit,R_ro_o,R_ta_o,R_ga_o,R_ro_n,R_ta_n,R_ga_n,R_ro_update,R_ta_update,flag,R_delta,
		margprobx_o,margprobx_n,somega,gradient_ro, se_ro, p_o_r_t,
		probmat,condprobx_o,condprobx_n,gradient_ta, se_ta, p_r_t,
		gradient_ga,p_r,se_ga,R_Nbootstrap, OR_common_att_rep, OR_common_obj_rep, report_common_att, report_common_obj);
	////////////////////////////////////////////////////////////



	for (int i = 1; i < _nO+1; i++)
	{
		for (int j = 1; j < _nF+1; j++)
		{
			_R_ro_n[  (i-1) * (_nF)  + (j-1) ] = R_ro_n[i][j];
		}
	}
	for (int i = 1; i < _nA+1; i++)
	{
		for (int j = 1; j < _nF+1; j++)
		{
			for (int k = 1; k < _nT+1; k++)
			{
				_R_ta_n[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = R_ta_n[i][j][k];	
			}	
		}
	}
	for (int i = 1; i < _nT+1; i++)
	{
		_R_ga_n[i-1] = R_ga_n[i];		
	}

	*_R_logpost_n = logpost_n;
	*_R_loglik_n = loglik_n;


	if (flag == true)
	{
		for (int i = 1; i < _nA+1; i++)
		{
			for (int j = 1; j < _nF+1; j++)
			{
				for (int k = 1; k < _nT+1; k++)
				{
					_R_gradient_ta[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = gradient_ta[i][j][k];
					_R_se_ta[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = se_ta[i][j][k];
				}	
			}
		}

		
		for (int i = 1; i < _nO+1; i++)
		{
			for (int j = 1; j < _nF+1; j++)
			{
				_R_gradient_ro[  (i-1) * (_nF) + (j-1) ] = gradient_ro[i][j];
				_R_se_ro[  (i-1) * (_nF)  + (j-1) ] = se_ro[i][j];
			}
		}

		for (int i = 1; i < _nT+1; i++)
		{
			_R_gradient_ga[i-1] = gradient_ga[i];
			_R_se_ga[i-1] = se_ga[i];
		}

		for (int i = 1; i < _nT+1; i++)
		{
			for (int j = 1; j < _nR+1; j++)
			{
				_R_probmat[ (i-1) * (_nR) + (j-1)] = probmat[i][j];
			}
		}
		
		
			// bootstrapping variables
		for (int i = 1; i < _Nbootstrap+1; i++)
		{
			for (int j = 1; j < _nO+1; j++)
			{
				for (int k = 1; k < N_attpair+1; k++)
				{
					_R_OR_common_att_rep[ (i-1) * (_nO) * (N_attpair) + (j-1) * (N_attpair) + (k-1)] = OR_common_att_rep[i][j][k];
					
				}
			}
		}
		for (int i = 1; i < _Nbootstrap+1; i++)
		{
			for (int j = 1; j < _nA+1; j++)
			{
				for (int k = 1; k < N_objpair+1; k++)
				{
					_R_OR_common_obj_rep[ (i-1) * (_nA) * (N_objpair) + (j-1) * (N_objpair) + (k-1)] = OR_common_obj_rep[i][j][k];
					
				}
			}
		}
		for (int i = 1; i < _nO*N_attpair+1; i++)
		{
			for (int j = 1; j < 5+1; j++)
			{
				_R_report_common_att[  (i-1) * (5) + (j-1)] = report_common_att[i][j];
			}
		}
		for (int i = 1; i < _nA*N_objpair+1; i++)
		{
			for (int j = 1; j < 5+1; j++)
			{
				_R_report_common_obj[  (i-1) * (5) + (j-1)] = report_common_obj[i][j];
			}
		}

		

	}


	// destructor

	for (int i = 0; i < _nO+1; i++) 
	{
		for (int ii = 0; ii < _nA+1; ii++)
		{
			delete R_data[i][ii];
		}
		delete R_data[i];
	}
	delete[] R_data;
	for (int i = 0; i < _nS+1; i++)
	{
		delete R_patS[i];
	}
	delete[] R_patS;
	for (int i = 0; i < _nO+1; i++)
	{
		delete R_ro_o[i];
		delete R_ro_n[i];
		delete R_ro_update[i];
	}
	delete[] R_ro_o;
	delete[] R_ro_n;
	delete[] R_ro_update;
	
	for (int i = 0; i < _nA+1; i++)
	{
		for (int ii = 0; ii < _nF+1; ii++)
		{
			delete R_ta_o[i][ii];
			delete R_ta_n[i][ii];
			delete R_ta_update[i][ii];
		}
		delete R_ta_o[i];
		delete R_ta_n[i];
		delete R_ta_update[i];
	}
	delete[] R_ta_o;
	delete[] R_ta_n;
	delete[] R_ta_update;
	delete[] R_ga_o;
	delete[] R_ga_n;

	for (int i = 0; i < _nT+1; i++)
	{
		delete probmat[i];
	}
	delete[] probmat;

	if (flag == true)
	{
		for (int i = 0; i < _nA+1; i++)
		{
			for (int ii = 0; ii < _nF+1; ii++)
			{
				delete gradient_ta[i][ii];
			}
			delete gradient_ta[i];
		}
		delete[] gradient_ta; 

		for (int i = 0; i < _nO+1; i++)
		{
			delete gradient_ro[i];
		}
		delete[] gradient_ro;

		delete[] gradient_ga;

		for (int i = 0; i < _nA+1; i++)
		{
			for (int ii = 0; ii < _nF+1; ii++)
			{
				delete se_ta[i][ii];
			}
			delete se_ta[i];
		}
		delete[] se_ta;

		for (int i = 0; i < _nO+1; i++)
		{
			delete se_ro[i];
		}
		delete[] se_ro;

		delete[] se_ga;
		
		
		
		
		// OR_common_att_rep
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			for (int ii = 0; ii < _nO+1; ii++)
			{
				delete OR_common_att_rep[i][ii];
			}
			delete OR_common_att_rep[i];
		}
		delete[] OR_common_att_rep;
		
		// OR_common_obj_rep
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			for (int ii = 0; ii < _nA+1; ii++)
			{
				delete OR_common_obj_rep[i][ii];
			}
			delete OR_common_obj_rep[i];
		}
		delete[] OR_common_obj_rep;
				
		// report_common_att
		
		for (int i = 0; i < _nO*N_attpair+1; i++)
		{
			delete report_common_att[i];
		}
		delete[] report_common_att;
		
	   // report_common_obj	
	   for (int i = 0; i < _nA*N_objpair+1; i++)
		{
			delete report_common_obj[i];
		}
		delete[] report_common_obj;		
		
	
	}
}

// estimate disjunctive model M_3
extern "C" void PlFm_XZ_YZ_DC(int * _R_data, int * _R_nR, int * _R_nO, int * _R_nA, int * _R_nF, int * _R_nT, int * _R_patS, 
						  double * _R_emcrit,
						  double * _R_ro_n, double * _R_ta_n, double * _R_ga_n,
						  double * _R_ro_update, double * _R_ta_update,
						  double * _R_gradient_ro, double * _R_gradient_ta, double * _R_gradient_ga,
						  double * _R_se_ro, double * _R_se_ta, double * _R_se_ga, double * _R_delta, double * _R_logpost_n, double * _R_loglik_n, double * _R_probmat, int * _R_flag,
						  int * _R_Nbootstrap, double * _R_OR_common_att_rep, double * _R_OR_common_obj_rep, double * _R_report_common_att, double * _R_report_common_obj)
{
	// variables
	extended3 condprobx_o,condprobx_n,somega,margprobx_o,margprobx_n,p_o_r_t,gradient_ta,se_ta,gradient_ro,se_ro;
	extended2 probmat, p_r_t;
    extended1 p_r,gradient_ga,se_ga;

	binary3 R_data; int R_nR; int R_nO; int R_nA; int R_nF; int R_nT; int R_Nbootstrap; binary2 R_patS;
	long double R_emcrit; extended3 R_ro_o; extended3 R_ta_o; extended1 R_ga_o;
	extended3 R_ro_n; extended3 R_ta_n; extended1 R_ga_n; extended3 R_ro_update; extended3 R_ta_update;

	int R_flag;
	double R_delta;

	R_nR = *_R_nR; 
	R_nO = *_R_nO; 
	R_nA = *_R_nA; 
	R_nF = *_R_nF; 
	R_nT = *_R_nT;
	R_Nbootstrap = *_R_Nbootstrap;

	R_flag = *_R_flag;
	R_delta = *_R_delta;
	
	bool flag;
	if (R_flag == 0)
		flag = false;
	else
		flag = true;

	_nO = R_nO;
        _nA = R_nA;
        _nR = R_nR;
	_nT = R_nT;
	_nF = R_nF;

	R_data = new bool **[R_nO+1]; 
	for (int i = 0; i < R_nO+1; i++) 
	{
		R_data[i] = new bool *[R_nA+1];
		for (int ii = 0; ii < R_nA+1; ii++)
		{
			R_data[i][ii] = new bool [R_nR+1];
		}
	}
	for (int i = 0; i < R_nO; i++)
	{
		for (int j = 0; j < R_nA; j++)
		{
			for (int k = 0; k < R_nR; k++)
			{
				R_data[i+1][j+1][k+1] = (bool)_R_data[i * (R_nA) * (R_nR) + j * (R_nR) + k];
			}
		}
	}


	calculate_dims();

	somega = new long double ** [_nR+1]; 
	condprobx_o = new long double ** [_nS+1]; 
	margprobx_o = new long double ** [_nS+1]; 
	condprobx_n = new long double ** [_nS+1]; 
	margprobx_n = new long double ** [_nS+1]; 
	probmat = new long double * [_nT+1]; 
	p_o_r_t = new long double ** [_nO+1]; 
	p_r_t = new long double * [_nR+1];
	p_r = new long double [_nR+1];
	CreateVariables( somega, condprobx_o, condprobx_n, margprobx_o, margprobx_n, probmat, p_o_r_t, p_r_t, p_r);
	

	R_patS = new bool * [_nS+1]; 
	for (int i = 0; i < _nS+1; i++)
	{
		R_patS[i] = new bool [R_nF+1];
	}
	for (int i = 0; i < _nS; i++)
	{
		for (int j = 0; j < R_nF; j++)
		{
			R_patS[i+1][j+1] = (bool)_R_patS[i * (R_nF) + j ];			
		}
	}
	R_emcrit = *_R_emcrit; 
	
	R_ro_o = new long double ** [_nO+1] ; R_ro_n = new long double ** [_nO+1]; R_ro_update = new long double ** [_nO+1];
	R_ta_o = new long double ** [_nA+1] ; R_ta_n = new long double ** [_nA+1] ; R_ta_update = new long double ** [_nA+1] ;
	R_ga_o = new long double [_nT+1]; R_ga_n = new long double [_nT+1];
	Create_ro_ta(R_ro_o,R_ta_o);
	Create_ro_ta(R_ro_n,R_ta_n);
	Create_ro_ta(R_ro_update,R_ta_update);
	seed1= 0; 
	seed2= 1; 
	seed3= 2; 
	generate_true_parameters(R_ro_o,R_ta_o,R_ga_o);

	for (int i = 0; i < _nO; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			for (int k = 0; k < _nT; k++)
			{

				R_ro_n[i+1][j+1][k+1] =  _R_ro_n[  i * (_nF) * (_nT) + j * (_nT) + k];
			}
		}
	}

	for (int i = 0; i < _nO; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			for (int k = 0; k < _nT; k++)
			{

				R_ro_update[i+1][j+1][k+1] =  _R_ro_update[  i * (_nF) * (_nT) + j * (_nT) + k];
			}
		}
	}


	for (int i = 0; i < _nA; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			for (int k = 0; k < _nT; k++)
			{
				R_ta_n[i+1][j+1][k+1] =  _R_ta_n[  i * (_nF) * (_nT) + j * (_nT) + k];
			}

		}
	}

	for (int i = 0; i < _nA; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			for (int k = 0; k < _nT; k++)
			{
				R_ta_update[i+1][j+1][k+1] =  _R_ta_update[  i * (_nF) * (_nT) + j * (_nT) + k];
			}

		}
	}

	for (int i = 0; i < _nT; i++)
	{	
		R_ga_n[i+1] =  _R_ga_n[i];	
	}
	
	// bootstrapping variables
	extended3 OR_common_att_rep, OR_common_obj_rep;
	extended2 report_common_att, report_common_obj;
	int N_attpair, N_objpair;
	N_attpair = (_nA*(_nA-1)) / 2;
	N_objpair = (_nO*(_nO-1)) / 2;

	if (flag == true)
	{

		gradient_ta = new long double ** [_nA+1]; 
		for (int i = 0; i < _nA+1; i++)
		{
			gradient_ta[i] = new long double * [_nF+1];
			for (int ii = 0; ii < _nF+1; ii++)
			{
				gradient_ta[i][ii] = new long double [_nT+1];
			}
		}


		gradient_ro = new long double ** [_nO+1]; 
		for (int i = 0; i < _nO+1; i++)
		{
			gradient_ro[i] = new long double * [_nF+1];
			for (int ii = 0; ii < _nF+1; ii++)
			{
				gradient_ro[i][ii] = new long double [_nT+1];
			}
		}


		gradient_ga = new long double [_nT+1];



		se_ta = new long double ** [_nA+1]; 
		for (int i = 0; i < _nA+1; i++)
		{
			se_ta[i] = new long double *[_nF+1];
			for (int ii = 0; ii < _nF+1; ii++)
			{
				se_ta[i][ii] = new long double [_nT+1];
			}
		}


		se_ro = new long double ** [_nO+1]; 
		for (int i = 0; i < _nO+1; i++)
		{
			se_ro[i] = new long double * [_nF+1];
			for (int ii = 0; ii < _nF+1; ii++)
			{
				se_ro[i][ii] = new long double [_nT+1];
			}
		}


		se_ga = new long double [_nT+1];
		
		
		// bootstrapping variables
		

		// construct OR_common_att_rep,OR_common_obj_rep
		OR_common_att_rep = new long double ** [_Nbootstrap+1] ;
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			OR_common_att_rep[i] = new long double * [_nO+1];
			for (int ii = 0; ii < _nO+1; ii++)
			{
				OR_common_att_rep[i][ii] = new long double [N_attpair+1];
			}
		}
		OR_common_obj_rep = new long double ** [_Nbootstrap+1] ;
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			OR_common_obj_rep[i] = new long double * [_nA+1];
			for (int ii = 0; ii < _nA+1; ii++)
			{
				OR_common_obj_rep[i][ii] = new long double [N_objpair+1];
			}
		}
		//construct report_common_att & report_common_obj
		report_common_att = new long double * [_nO*N_attpair+1] ;
		for (int i = 0; i < _nO*N_attpair+1; i++)
		{
			report_common_att[i] = new long double [5+1];
		}
		report_common_obj = new long double * [_nA*N_objpair+1] ;
		for (int i = 0; i < _nA*N_objpair+1; i++)
		{
			report_common_obj[i] = new long double [5+1];
		}

		
	}
		else
	{
		se_ta = NULL;
		se_ga = NULL;
		se_ro = NULL;
		gradient_ta = NULL; 
		gradient_ga = NULL; 
		gradient_ro = NULL;
		OR_common_att_rep = NULL;
		OR_common_obj_rep = NULL;
		report_common_att = NULL;
		report_common_obj = NULL; 
	}

	//////////////////// call C function model M_3 ///////////////////////
	PlFm_XZ_YZ_intern_DC(R_data,R_nR,R_nO,R_nA,R_nF,R_nT,R_patS,R_emcrit,R_ro_o,R_ta_o,R_ga_o,R_ro_n,R_ta_n,R_ga_n,R_ro_update,R_ta_update,flag,R_delta,
		margprobx_o,margprobx_n,somega,gradient_ro, se_ro, p_o_r_t,
		probmat,condprobx_o,condprobx_n,gradient_ta, se_ta, p_r_t,
		gradient_ga,p_r,se_ga,R_Nbootstrap, OR_common_att_rep, OR_common_obj_rep, report_common_att, report_common_obj);
	////////////////////////////////////////////////////////////



	for (int i = 1; i < _nO+1; i++)
	{
		for (int j = 1; j < _nF+1; j++)
		{
			for (int k = 1; k < _nT+1; k++)
			{

				_R_ro_n[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = R_ro_n[i][j][k];
			}
		}
	}
	for (int i = 1; i < _nA+1; i++)
	{
		for (int j = 1; j < _nF+1; j++)
		{
			for (int k = 1; k < _nT+1; k++)
			{
				_R_ta_n[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = R_ta_n[i][j][k];	
			}	
		}
	}
	for (int i = 1; i < _nT+1; i++)
	{
		_R_ga_n[i-1] = R_ga_n[i];		
	}

	*_R_logpost_n = logpost_n;
	*_R_loglik_n = loglik_n;


	if (flag == true)
	{
		for (int i = 1; i < _nA+1; i++)
		{
			for (int j = 1; j < _nF+1; j++)
			{
				for (int k = 1; k < _nT+1; k++)
				{
					_R_gradient_ta[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = gradient_ta[i][j][k];
					_R_se_ta[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = se_ta[i][j][k];
				}	
			}
		}

		
		for (int i = 1; i < _nO+1; i++)
		{
			for (int j = 1; j < _nF+1; j++)
			{
				for (int k = 1; k < _nT+1; k++)
				{
					_R_gradient_ro[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = gradient_ro[i][j][k];
					_R_se_ro[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = se_ro[i][j][k];
				}
			}
		}

		for (int i = 1; i < _nT+1; i++)
		{
			_R_gradient_ga[i-1] = gradient_ga[i];
			_R_se_ga[i-1] = se_ga[i];
		}

		for (int i = 1; i < _nT+1; i++)
		{
			for (int j = 1; j < _nR+1; j++)
			{
				_R_probmat[ (i-1) * (_nR) + (j-1)] = probmat[i][j];
			}
		}
		
		
		// bootstrapping variables
		for (int i = 1; i < _Nbootstrap+1; i++)
		{
			for (int j = 1; j < _nO+1; j++)
			{
				for (int k = 1; k < N_attpair+1; k++)
				{
					_R_OR_common_att_rep[ (i-1) * (_nO) * (N_attpair) + (j-1) * (N_attpair) + (k-1)] = OR_common_att_rep[i][j][k];
					
				}
			}
		}
		for (int i = 1; i < _Nbootstrap+1; i++)
		{
			for (int j = 1; j < _nA+1; j++)
			{
				for (int k = 1; k < N_objpair+1; k++)
				{
					_R_OR_common_obj_rep[ (i-1) * (_nA) * (N_objpair) + (j-1) * (N_objpair) + (k-1)] = OR_common_obj_rep[i][j][k];
					
				}
			}
		}
		for (int i = 1; i < _nO*N_attpair+1; i++)
		{
			for (int j = 1; j < 5+1; j++)
			{
				_R_report_common_att[  (i-1) * (5) + (j-1)] = report_common_att[i][j];
			}
		}
		for (int i = 1; i < _nA*N_objpair+1; i++)
		{
			for (int j = 1; j < 5+1; j++)
			{
				_R_report_common_obj[  (i-1) * (5) + (j-1)] = report_common_obj[i][j];
			}
		}
		

	}


	// destructor

	for (int i = 0; i < _nO+1; i++) 
	{
		for (int ii = 0; ii < _nA+1; ii++)
		{
			delete R_data[i][ii];
		}
		delete R_data[i];
	}
	delete[] R_data;
	for (int i = 0; i < _nS+1; i++)
	{
		delete R_patS[i];
	}
	delete[] R_patS;
	for (int i = 0; i < _nO+1; i++)
	{
		for (int ii = 0; ii < _nF+1; ii++)
		{
			delete R_ro_o[i][ii];
			delete R_ro_n[i][ii];
			delete R_ro_update[i][ii];
		}
		delete R_ro_o[i];
		delete R_ro_n[i];
		delete R_ro_update[i];
	}
	delete[] R_ro_o;
	delete[] R_ro_n;
	delete[] R_ro_update;
	for (int i = 0; i < _nA+1; i++)
	{
		for (int ii = 0; ii < _nF+1; ii++)
		{
			delete R_ta_o[i][ii];
			delete R_ta_n[i][ii];
			delete R_ta_update[i][ii];
		}
		delete R_ta_o[i];
		delete R_ta_n[i];
		delete R_ta_update[i];
	}
	delete[] R_ta_o;
	delete[] R_ta_n;
	delete[] R_ta_update;
	delete[] R_ga_o;
	delete[] R_ga_n;

	for (int i = 0; i < _nT+1; i++)
	{
		delete probmat[i];
	}
	delete[] probmat;

	if (flag == true)
	{
		for (int i = 0; i < _nA+1; i++)
		{
			for (int ii = 0; ii < _nF+1; ii++)
			{
				delete gradient_ta[i][ii];
			}
			delete gradient_ta[i];
		}
		delete[] gradient_ta; 

		for (int i = 0; i < _nO+1; i++)
		{
			for (int ii = 0; ii < _nF+1; ii++)
			{
				delete gradient_ro[i][ii];
			}
			delete gradient_ro[i];
		}
		delete[]  gradient_ro;

		delete[] gradient_ga;

		for (int i = 0; i < _nA+1; i++)
		{
			for (int ii = 0; ii < _nF+1; ii++)
			{
				delete se_ta[i][ii];
			}
			delete se_ta[i];
		}
		delete[] se_ta;

		for (int i = 0; i < _nO+1; i++)
		{
			for (int ii = 0; ii < _nF+1; ii++)
			{
				delete se_ro[i][ii];
			}
			delete se_ro[i];
		}
		delete[] se_ro;

		delete[] se_ga;
		
		
	// OR_common_att_rep
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			for (int ii = 0; ii < _nO+1; ii++)
			{
				delete OR_common_att_rep[i][ii];
			}
			delete OR_common_att_rep[i];
		}
		delete[] OR_common_att_rep;
		
		// OR_common_obj_rep
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			for (int ii = 0; ii < _nA+1; ii++)
			{
				delete OR_common_obj_rep[i][ii];
			}
			delete OR_common_obj_rep[i];
		}
		delete[] OR_common_obj_rep;
				
		// report_common_att
		
		for (int i = 0; i < _nO*N_attpair+1; i++)
		{
			delete report_common_att[i];
		}
		delete[] report_common_att;
		
	   // report_common_obj	
	   for (int i = 0; i < _nA*N_objpair+1; i++)
		{
			delete report_common_obj[i];
		}
		delete[] report_common_obj;			
	
	}
}


// function to estimate additive LCPLFM with constant classification of objects and heterogeneous object parameters (model M_1)
extern "C" void PlFm_XZ_Y_intern_ADD(binary3 R_data, int R_nR, int R_nO, int R_nA, int R_nF, int R_nT, binary2 R_patS, 
		  long double R_emcrit,extended3 R_ro_o, extended2 R_ta_o, extended1 R_ga_o,
		  extended3 R_ro_n, extended2 R_ta_n, extended1 R_ga_n, extended3 R_ro_update, extended2 R_ta_update, bool flag, double R_delta,
		 extended3 margprobx_o,extended3 margprobx_n,extended3 somega,extended3 gradient_ro,extended3  se_ro,extended3  p_o_r_t,
		 extended2 probmat,extended2 condprobx_o,extended2 condprobx_n,extended2 gradient_ta,extended2 se_ta,extended2 p_r_t,
		 extended1 gradient_ga,extended1 p_r,extended1 se_ga,
		 int R_Nbootstrap, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att, extended2 report_common_obj)
{


	extended3 ro_o; extended3 ro_n; extended3 ro_update;
	extended2 ta_o; extended2 ta_n; extended2 ta_update;
	extended1 ga_o; extended1 ga_n;
	

	// pass all data from R to C
	dataB = R_data;
	_nR = R_nR;
	_nO = R_nO;
	_nA = R_nA;
	_nF = R_nF;
	_nT = R_nT;
	_Nbootstrap = R_Nbootstrap;
	patS = R_patS;
	emcrit = R_emcrit;
	ro_o = R_ro_o;
	ta_o = R_ta_o;
	ga_o = R_ga_o;
	ro_n = R_ro_n;
	ta_n = R_ta_n;
	ga_n = R_ga_n;
	ro_update=R_ro_update;
	ta_update=R_ta_update;
	
	delta = R_delta;

	///////////////////// main program /////////////////////
	
    //compute conditional probability of observing an object-attribute association given object classification using old parameter values
	calculate_condprobX_ADD(ta_o,condprobx_o); 
	// compute marginal probability of object classification given person classification using old parameter values
	calculate_margprobX(ro_o,margprobx_o); 
	// compute posterior density using old parameter values
	logold= logposterior(ro_o,ta_o,ga_o,condprobx_o,margprobx_o,ro_update,ta_update);  

    //compute conditional probability of observing an object-attribute association given object classification using new parameter values
	calculate_condprobX_ADD(ta_n,condprobx_n);
    // compute marginal probability of object classification given person classification using new parameter values
	calculate_margprobX(ro_n,margprobx_n);
	// compute posterior density using new parameter values
	lognew= logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n,ro_update,ta_update);

	// compute value of convergence criterion
    calculate_conv( ro_o, ro_n, ta_o, ta_n, ga_o, ga_n, ro_update, ta_update);

	while( conv>emcrit)
	{ 
		// replace old parameter values by new parameter values
        replace( ro_o, ro_n, ta_o, ta_n, ga_o, ga_n,ro_update,ta_update, condprobx_o, condprobx_n, margprobx_o, margprobx_n);
		// compute posterior probabilities of classification person
        calculate_probmat(ga_o,condprobx_o,margprobx_o,probmat);
		// update the class size parameters
        emgamma(ga_n, probmat);
		// compute probability of object classification given data and classification of person 
        calculate_condprob_pattern(condprobx_o,margprobx_o,somega);
        // update object parameters
		update_emro(ro_n,probmat,ro_update);
		// update attribute parameters
		update_emta_ADD(condprobx_o,ta_o, ta_n,probmat,ta_update);
		// compute conditional probability of data given object classification and given person classification using new parameter values
		calculate_condprobX_ADD(ta_n,condprobx_n);
	    // compute marginal probability of object classification given person classification using new parameter values
		calculate_margprobX(ro_n,margprobx_n);
		// compute posterior density using new parameter values
        lognew= logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n,ro_update,ta_update);
       	// compute value of convergence criterion
		calculate_conv(ro_o, ro_n, ta_o, ta_n, ga_o, ga_n,ro_update,ta_update);
	} 

    // compute posterior density using new parameter values
	logpost_n = logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n,ro_update,ta_update);
	// compute log likelihood using new parameter values
	loglik_n  = loglikelihood (ga_n, condprobx_n, margprobx_n);
    
    // if flag==0 computation of gradient and standard errors can be skipped
    // if flag==1 computation of gradient and standard errors is needed
	if (flag == true)
	{
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t , p_r_t, p_r, probmat);
		calculate_gradient_ro(ro_n, ga_n, margprobx_n, condprobx_n, gradient_ro,p_o_r_t , p_r_t, p_r,ro_update);
		calculate_gradient_ta_ADD(ta_n, ga_n,condprobx_n, margprobx_n, gradient_ta,p_o_r_t , p_r_t, p_r,ta_update);
		calculate_gradient_ga(ga_n, gradient_ga,p_o_r_t , p_r_t, p_r);

		calculate_se_ro_ADD( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ro, p_o_r_t , p_r_t, p_r, probmat, ro_update);
		calculate_se_ta_ADD( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ta, p_o_r_t , p_r_t, p_r, probmat, ta_update);
		calculate_se_ga_ADD( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ga, p_o_r_t , p_r_t, p_r, probmat);

		calculate_probmat(ga_n,condprobx_n,margprobx_n,probmat);
	
		// bootstrap procedure
		bootstrap_dependency_XZ_Y_ADD(dataB, ro_n, ta_n, ga_n, OR_common_att_rep, OR_common_obj_rep, report_common_att, report_common_obj);

	}

	C_destructor( somega, condprobx_o, condprobx_n, margprobx_o, margprobx_n, probmat,   p_o_r_t, p_r_t, p_r);
}

// function to estimate additive LCPLFM with fixed object classification and heterogeneous attribute parameters (Model M_2)
extern "C" void PlFm_X_YZ_intern_ADD(binary3 R_data, int R_nR, int R_nO, int R_nA, int R_nF, int R_nT, binary2 R_patS, 
		  long double R_emcrit,extended2 R_ro_o, extended3 R_ta_o, extended1 R_ga_o,
		  extended2 R_ro_n, extended3 R_ta_n, extended1 R_ga_n,extended2 R_ro_update, extended3 R_ta_update, bool flag, double R_delta,
		 extended2 margprobx_o,extended2 margprobx_n,extended3 somega,extended2 gradient_ro,extended2  se_ro,extended3  p_o_r_t,
		 extended2 probmat,extended3 condprobx_o,extended3 condprobx_n,extended3 gradient_ta,extended3 se_ta,extended2 p_r_t,
		 extended1 gradient_ga,extended1 p_r,extended1 se_ga,
		  int R_Nbootstrap, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att, extended2 report_common_obj)
{

	extended2 ro_o; extended2 ro_n; extended2 ro_update;
	extended3 ta_o; extended3 ta_n; extended3 ta_update;
	extended1 ga_o; extended1 ga_n;

	// pass all data from R to C
	dataB = R_data;
	_nR = R_nR;
	_nO = R_nO;
	_nA = R_nA;
	_nF = R_nF;
	_nT = R_nT;
	_Nbootstrap = R_Nbootstrap;
	patS = R_patS;
	emcrit = R_emcrit;
	ro_o = R_ro_o;
	ta_o = R_ta_o;
	ga_o = R_ga_o;
	ro_n = R_ro_n;
	ta_n = R_ta_n;
	ga_n = R_ga_n;
	ro_update=R_ro_update;
	ta_update=R_ta_update;
	delta = R_delta;

	///////////////////// main program /////////////////////
	
    //compute conditional probability of observing an object-attribute association given object and person classification using old parameter values
	calculate_condprobX_ADD(ta_o,condprobx_o);
	// compute marginal probability of object classification using old parameter values
	calculate_margprobX(ro_o,margprobx_o);
	// compute posterior density using old parameter values
    logold= logposterior(ro_o,ta_o,ga_o,condprobx_o,margprobx_o,ro_update,ta_update);

    //compute conditional probability of observing an object-attribute association given object and person classification using new parameter values
	calculate_condprobX_ADD(ta_n,condprobx_n);
	// compute marginal probability of object classification using new parameter values
	calculate_margprobX(ro_n,margprobx_n);
	// compute posterior density using new parameter values
    lognew= logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n,ro_update,ta_update);
    
    // compute value of convergence criterion
	calculate_conv( ro_o, ro_n, ta_o, ta_n, ga_o, ga_n, ro_update, ta_update);


	while( conv>emcrit)
	{ 
        // replace old parameter values by new parameter values
		replace( ro_o, ro_n, ta_o, ta_n, ga_o, ga_n, ro_update,ta_update, condprobx_o, condprobx_n, margprobx_o, margprobx_n);
		// compute posterior probabilities of classification person
		calculate_probmat(ga_o,condprobx_o,margprobx_o,probmat);
		// update the class size parameters
		emgamma(ga_n, probmat);
		// compute probability of object classification given data and classification of person 
		calculate_condprob_pattern(condprobx_o,margprobx_o,somega);
		// update object parameters
		update_emro(ro_n,probmat, ro_update);
		// update attribute parameters
		update_emta_ADD(condprobx_o,ta_o, ta_n,probmat, ta_update);
		// compute conditional probability of data given object classification and given person classification using new parameter values
		calculate_condprobX_ADD(ta_n,condprobx_n);
		// compute marginal probability of object classification given person classification using new parameter values
		calculate_margprobX(ro_n,margprobx_n);
		// compute posterior density using new parameter values
		lognew= logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n, ro_update, ta_update);
		// compute value of convergence criterion
		calculate_conv(ro_o, ro_n, ta_o, ta_n, ga_o, ga_n, ro_update, ta_update);

	} 
     // compute posterior density using new parameter values
	logpost_n = logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n, ro_update, ta_update);
	// compute loglikelihood using new parameter values
	loglik_n  = loglikelihood (ga_n, condprobx_n, margprobx_n);

    // if flag==0 computation of gradient and standard errors can be skipped
    // if flag==1 computation of gradient and standard errors is needed
	if (flag == true)
	{
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t , p_r_t, p_r, probmat);
		calculate_gradient_ro(ro_n, ga_n, margprobx_n, condprobx_n, gradient_ro,p_o_r_t , p_r_t, p_r, ro_update);
		calculate_gradient_ta_ADD(ta_n, ga_n,condprobx_n, margprobx_n, gradient_ta,p_o_r_t , p_r_t, p_r, ta_update);
		calculate_gradient_ga(ga_n, gradient_ga,p_o_r_t , p_r_t, p_r);

		calculate_se_ro_ADD( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ro, p_o_r_t , p_r_t, p_r, probmat, ro_update);
		calculate_se_ta_ADD( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ta, p_o_r_t , p_r_t, p_r, probmat, ta_update);
		calculate_se_ga_ADD( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ga, p_o_r_t , p_r_t, p_r, probmat);

		calculate_probmat(ga_n,condprobx_n,margprobx_n,probmat);
		
		// bootstrap procedure
		bootstrap_dependency_X_YZ_ADD(dataB, ro_n, ta_n, ga_n, OR_common_att_rep, OR_common_obj_rep, report_common_att, report_common_obj);

	}

	C_destructor( somega, condprobx_o, condprobx_n, margprobx_o, margprobx_n, probmat,   p_o_r_t, p_r_t, p_r);
}

// function to estimate additive LCPLFM with fixed object classification and heterogeneous object- and attribute parameters (model M_3)
extern "C" void PlFm_XZ_YZ_intern_ADD(binary3 R_data, int R_nR, int R_nO, int R_nA, int R_nF, int R_nT, binary2 R_patS, 
		  long double R_emcrit,extended3 R_ro_o, extended3 R_ta_o, extended1 R_ga_o,
		  extended3 R_ro_n, extended3 R_ta_n, extended1 R_ga_n, extended3 R_ro_update, extended3 R_ta_update, bool flag, double R_delta,
		 extended3 margprobx_o,extended3 margprobx_n,extended3 somega,extended3 gradient_ro,extended3  se_ro,extended3  p_o_r_t,
		 extended2 probmat,extended3 condprobx_o,extended3 condprobx_n,extended3 gradient_ta,extended3 se_ta,extended2 p_r_t,
		 extended1 gradient_ga,extended1 p_r,extended1 se_ga,
		 int R_Nbootstrap, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att, extended2 report_common_obj)
{

	extended3 ro_o; extended3 ro_n; extended3 ro_update;
	extended3 ta_o; extended3 ta_n; extended3 ta_update;
	extended1 ga_o; extended1 ga_n;

	// pass all data from R to C
	dataB = R_data;
	_nR = R_nR;
	_nO = R_nO;
	_nA = R_nA;
	_nF = R_nF;
	_nT = R_nT;
	patS = R_patS;
	emcrit = R_emcrit;
	ro_o = R_ro_o;
	ta_o = R_ta_o;
	ga_o = R_ga_o;
	ro_n = R_ro_n;
	ta_n = R_ta_n;
	ga_n = R_ga_n;
	ro_update=R_ro_update;
	ta_update=R_ta_update;
	delta = R_delta;

	///////////////////// main program /////////////////////
	
    //compute conditional probability of observing an object-attribute association given object and person classification using old parameter values
	calculate_condprobX_ADD(ta_o,condprobx_o);
	// compute marginal probability of object classification given person classification using old parameter values
	calculate_margprobX(ro_o,margprobx_o);
	// compute posterior density using old parameter values
    logold= logposterior(ro_o,ta_o,ga_o,condprobx_o,margprobx_o,ro_update,ta_update);

    //compute conditional probability of observing an object-attribute association given object and person classification using new parameter values
	calculate_condprobX_ADD(ta_n,condprobx_n);
	// compute marginal probability of object classification given person classification using new parameter values
	calculate_margprobX(ro_n,margprobx_n);
	// compute posterior density using new parameter values
    lognew= logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n,ro_update,ta_update);

    // compute value of convergence criterion
	calculate_conv( ro_o, ro_n, ta_o, ta_n, ga_o, ga_n,ro_update,ta_update);


	while( conv>emcrit)
	{ 
		// replace old parameter values by new parameter values
        replace( ro_o, ro_n, ta_o, ta_n, ga_o, ga_n, ro_update,ta_update,  condprobx_o, condprobx_n, margprobx_o, margprobx_n);
        // compute posterior probabilities of classification person
		calculate_probmat(ga_o,condprobx_o,margprobx_o,probmat);
		// update the class size parameters
		emgamma(ga_n, probmat);
		// compute probability of object classification given data and classification of person 
		calculate_condprob_pattern(condprobx_o,margprobx_o,somega);
		// update object parameters
		update_emro(ro_n,probmat,ro_update);
		// update attribute parameters
		update_emta_ADD(condprobx_o, ta_o, ta_n,probmat,ta_update);
		// compute conditional probability of data given object classification and given person classification using new parameter values
		calculate_condprobX_ADD(ta_n,condprobx_n);
		// compute marginal probability of object classification given person classification using new parameter values
		calculate_margprobX(ro_n,margprobx_n);
		// compute posterior density using new parameter values
		lognew= logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n,ro_update,ta_update);
		// compute value of convergence criterion
		calculate_conv(ro_o, ro_n, ta_o, ta_n, ga_o, ga_n,ro_update,ta_update);

	} 
    // compute posterior density using new parameter values
	logpost_n = logposterior(ro_n,ta_n,ga_n,condprobx_n,margprobx_n,ro_update,ta_update);
	// compute loglikelihood using new parameter values
	loglik_n  = loglikelihood (ga_n, condprobx_n, margprobx_n);

    // if flag==0 computation of gradient and standard errors can be skipped
    // if flag==1 computation of gradient and standard errors is needed
    
	if (flag == true)
	{
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t , p_r_t, p_r, probmat);
		calculate_gradient_ro(ro_n, ga_n, margprobx_n, condprobx_n, gradient_ro,p_o_r_t , p_r_t, p_r,ro_update);
		calculate_gradient_ta_ADD(ta_n, ga_n,condprobx_n, margprobx_n, gradient_ta,p_o_r_t , p_r_t, p_r,ta_update);
		calculate_gradient_ga(ga_n, gradient_ga,p_o_r_t , p_r_t, p_r);

		calculate_se_ro_ADD( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ro, p_o_r_t , p_r_t, p_r, probmat,ro_update);
		calculate_se_ta_ADD( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ta, p_o_r_t , p_r_t, p_r, probmat,ta_update);
		calculate_se_ga_ADD( ro_n, ta_n, ga_n, condprobx_n, margprobx_n, se_ga, p_o_r_t , p_r_t, p_r, probmat);

		calculate_probmat(ga_n,condprobx_n,margprobx_n,probmat);
		
		// bootstrap procedure
		bootstrap_dependency_XZ_YZ_ADD(dataB, ro_n, ta_n, ga_n, OR_common_att_rep, OR_common_obj_rep, report_common_att, report_common_obj);

	}

	C_destructor( somega, condprobx_o, condprobx_n, margprobx_o, margprobx_n, probmat,   p_o_r_t, p_r_t, p_r);
}




// estimate additive model M_1
extern "C" void PlFm_XZ_Y_ADD(int * _R_data, int * _R_nR, int * _R_nO, int * _R_nA, int * _R_nF, int * _R_nT, int * _R_patS, 
						  double * _R_emcrit,
						  double * _R_ro_n, double * _R_ta_n, double * _R_ga_n,
						  double * _R_ro_update, double * _R_ta_update,
						  double * _R_gradient_ro, double * _R_gradient_ta, double * _R_gradient_ga,
						  double * _R_se_ro, double * _R_se_ta, double * _R_se_ga, double * _R_delta, double * _R_logpost_n, double * _R_loglik_n, double * _R_probmat, int * _R_flag,
						  int * _R_Nbootstrap, double * _R_OR_common_att_rep, double * _R_OR_common_obj_rep, double * _R_report_common_att, double * _R_report_common_obj)
{
	// variables
	
	extended3 margprobx_o,margprobx_n,somega,gradient_ro, se_ro, p_o_r_t;
	extended2 probmat,condprobx_o,condprobx_n,gradient_ta, se_ta, p_r_t;
	extended1 gradient_ga,p_r,se_ga;


	binary3 R_data; int R_nR; int R_nO; int R_nA; int R_nF; int R_nT; int R_Nbootstrap; binary2 R_patS; 
	long double R_emcrit;extended3 R_ro_o; extended2 R_ta_o; extended1 R_ga_o;
	extended3 R_ro_n; extended2 R_ta_n; extended1 R_ga_n; extended3 R_ro_update; extended2 R_ta_update;
	

	int R_flag;
	double R_delta;

	R_nR = *_R_nR; 
	R_nO = *_R_nO; 
	R_nA = *_R_nA; 
	R_nF = *_R_nF; 
	R_nT = *_R_nT;
	R_Nbootstrap = *_R_Nbootstrap;

	R_flag = *_R_flag;
	R_delta = *_R_delta;
	
	bool flag;
	if (R_flag == 0)
		flag = false;
	else
		flag = true;

	_nO = R_nO;
        _nA = R_nA;
        _nR = R_nR;
	_nT = R_nT;
	_nF = R_nF;

	R_data = new bool **[R_nO+1]; 
	for (int i = 0; i < R_nO+1; i++) 
	{
		R_data[i] = new bool *[R_nA+1];
		for (int ii = 0; ii < R_nA+1; ii++)
		{
			R_data[i][ii] = new bool [R_nR+1];
		}
	}
	for (int i = 0; i < R_nO; i++)
	{
		for (int j = 0; j < R_nA; j++)
		{
			for (int k = 0; k < R_nR; k++)
			{
				R_data[i+1][j+1][k+1] = (bool)_R_data[i * (R_nA) * (R_nR) + j * (R_nR) + k];
			}
		}
	}

	calculate_dims();

	somega = new long double ** [_nR+1]; 
	condprobx_o = new long double * [_nS+1]; 
	margprobx_o = new long double ** [_nS+1]; 
	condprobx_n = new long double * [_nS+1]; 
	margprobx_n = new long double ** [_nS+1]; 
	probmat = new long double * [_nT+1]; 
	p_o_r_t = new long double ** [_nO+1]; 
	p_r_t = new long double * [_nR+1];
	p_r = new long double [_nR+1];
	CreateVariables( somega, condprobx_o, condprobx_n, margprobx_o, margprobx_n, probmat, p_o_r_t, p_r_t, p_r);

	R_patS = new bool * [_nS+1]; 
	for (int i = 0; i < _nS+1; i++)
	{
		R_patS[i] = new bool [R_nF+1];
	}
	for (int i = 0; i < _nS; i++)
	{
		for (int j = 0; j < R_nF; j++)
		{
			R_patS[i+1][j+1] = (bool)_R_patS[i * (R_nF) + j ];			
		}
	}
	R_emcrit = *_R_emcrit; 
	R_ro_o = new long double ** [_nO+1] ; R_ro_n = new long double ** [_nO+1]; R_ro_update = new long double ** [_nO+1];
	R_ta_o = new long double * [_nA+1] ; R_ta_n = new long double * [_nA+1] ; R_ta_update = new long double * [_nA+1] ;
	R_ga_o = new long double [_nT+1]; R_ga_n = new long double [_nT+1];
	Create_ro_ta(R_ro_o,R_ta_o);
	Create_ro_ta(R_ro_n,R_ta_n);
	Create_ro_ta(R_ro_update,R_ta_update);
	seed1= 0; 
	seed2= 1; 
	seed3= 2; 
	generate_true_parameters(R_ro_o,R_ta_o,R_ga_o);

	for (int i = 0; i < _nO; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			for (int k = 0; k < _nT; k++)
			{

				R_ro_n[i+1][j+1][k+1] =  _R_ro_n[  i * (_nF) * (_nT) + j * (_nT) + k];
			}
		}
	}

    for (int i = 0; i < _nO; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			for (int k = 0; k < _nT; k++)
			{

				R_ro_update[i+1][j+1][k+1] =  _R_ro_update[  i * (_nF) * (_nT) + j * (_nT) + k];
			}
		}
	}
	
    for (int i = 0; i < _nA; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
		
			R_ta_n[i+1][j+1] =  _R_ta_n[  i * (_nF) + j];	
		}
	}

    for (int i = 0; i < _nA; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
		
			R_ta_update[i+1][j+1] =  _R_ta_update[  i * (_nF) + j];	
		}
	}

	for (int i = 0; i < _nT; i++)
	{
	
		R_ga_n[i+1] =  _R_ga_n[i];	
	}

	// bootstrapping variables
	extended3 OR_common_att_rep, OR_common_obj_rep;
	extended2 report_common_att, report_common_obj;
	int N_attpair, N_objpair;
	N_attpair = (_nA*(_nA-1)) / 2;
	N_objpair = (_nO*(_nO-1)) / 2;


	if (flag == true)
	{

		gradient_ta = new long double * [_nA+1]; 
		for (int i = 0; i < _nA+1; i++)
		{
			gradient_ta[i] = new long double [_nF+1];
		}


		gradient_ro = new long double ** [_nO+1]; 
		for (int i = 0; i < _nO+1; i++)
		{
			gradient_ro[i] = new long double * [_nF+1];
			for (int ii = 0; ii < _nF+1; ii++)
			{
				gradient_ro[i][ii] = new long double [_nT+1];
			}
		}


		gradient_ga = new long double [_nT+1];



		se_ta = new long double * [_nA+1]; 
		for (int i = 0; i < _nA+1; i++)
		{
			se_ta[i] = new long double [_nF+1];
		}


		se_ro = new long double ** [_nO+1]; 
		for (int i = 0; i < _nO+1; i++)
		{
			se_ro[i] = new long double * [_nF+1];
			for (int ii = 0; ii < _nF+1; ii++)
			{
				se_ro[i][ii] = new long double [_nT+1];
			}
		}


		se_ga = new long double [_nT+1];

		// bootstrapping variables

		

		// construct OR_common_att_rep,OR_common_obj_rep
		OR_common_att_rep = new long double ** [_Nbootstrap+1] ;
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			OR_common_att_rep[i] = new long double * [_nO+1];
			for (int ii = 0; ii < _nO+1; ii++)
			{
				OR_common_att_rep[i][ii] = new long double [N_attpair+1];
			}
		}
		OR_common_obj_rep = new long double ** [_Nbootstrap+1] ;
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			OR_common_obj_rep[i] = new long double * [_nA+1];
			for (int ii = 0; ii < _nA+1; ii++)
			{
				OR_common_obj_rep[i][ii] = new long double [N_objpair+1];
			}
		}
		//construct report_common_att & report_common_obj
		report_common_att = new long double * [_nO*N_attpair+1] ;
		for (int i = 0; i < _nO*N_attpair+1; i++)
		{
			report_common_att[i] = new long double [5+1];
		}
		report_common_obj = new long double * [_nA*N_objpair+1] ;
		for (int i = 0; i < _nA*N_objpair+1; i++)
		{
			report_common_obj[i] = new long double [5+1];
		}

	}
	else
	{
		se_ta = NULL;
		se_ga = NULL;
		se_ro = NULL;
		gradient_ta = NULL; 
		gradient_ga = NULL; 
		gradient_ro = NULL; 
		OR_common_att_rep = NULL;
		OR_common_obj_rep = NULL;
		report_common_att = NULL;
		report_common_obj = NULL;

	}

	//////////////////// call C function model M_1 ///////////////////////
	PlFm_XZ_Y_intern_ADD(R_data,R_nR,R_nO,R_nA,R_nF,R_nT,R_patS,R_emcrit,R_ro_o,R_ta_o,R_ga_o,R_ro_n,R_ta_n,R_ga_n,R_ro_update,R_ta_update,flag,R_delta,
		margprobx_o,margprobx_n,somega,gradient_ro, se_ro, p_o_r_t,
		probmat,condprobx_o,condprobx_n,gradient_ta, se_ta, p_r_t,
		gradient_ga,p_r,se_ga,R_Nbootstrap, OR_common_att_rep, OR_common_obj_rep, report_common_att, report_common_obj);
	////////////////////////////////////////////////////////////



	for (int i = 1; i < _nO+1; i++)
	{
		for (int j = 1; j < _nF+1; j++)
		{
			for (int k = 1; k < _nT+1; k++)
			{

				_R_ro_n[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = R_ro_n[i][j][k];
			}
		}
	}
	for (int i = 1; i < _nA+1; i++)
	{
		for (int j = 1; j < _nF+1; j++)
		{
	
			_R_ta_n[  (i-1) * (_nF) + (j-1)] = R_ta_n[i][j];	
		}
	}
	for (int i = 1; i < _nT+1; i++)
	{

		_R_ga_n[i-1] = R_ga_n[i];		
	}

	*_R_logpost_n = logpost_n;
	*_R_loglik_n = loglik_n;


	if (flag == true)
	{
		
		for (int i = 1; i < _nA+1; i++)
		{
			for (int j = 1; j < _nF+1; j++)
			{
				_R_gradient_ta[  (i-1) * (_nF) + (j-1)] = gradient_ta[i][j];
				_R_se_ta[  (i-1) * (_nF) + (j-1)] = se_ta[i][j];
			}
		}

		
		for (int i = 1; i < _nO+1; i++)
		{
			for (int j = 1; j < _nF+1; j++)
			{
				for (int k = 1; k < _nT+1; k++)
				{
					_R_gradient_ro[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = gradient_ro[i][j][k];
					_R_se_ro[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = se_ro[i][j][k];
				}
			}
		}

		for (int i = 1; i < _nT+1; i++)
		{
			_R_gradient_ga[i-1] = gradient_ga[i];
			_R_se_ga[i-1] = se_ga[i];
		}

		for (int i = 1; i < _nT+1; i++)
		{
			for (int j = 1; j < _nR+1; j++)
			{
				_R_probmat[ (i-1) * (_nR) + (j-1)] = probmat[i][j];
			}
		}

		// bootstrapping variables
		for (int i = 1; i < _Nbootstrap+1; i++)
		{
			for (int j = 1; j < _nO+1; j++)
			{
				for (int k = 1; k < N_attpair+1; k++)
				{
					_R_OR_common_att_rep[ (i-1) * (_nO) * (N_attpair) + (j-1) * (N_attpair) + (k-1)] = OR_common_att_rep[i][j][k];
					
				}
			}
		}
		for (int i = 1; i < _Nbootstrap+1; i++)
		{
			for (int j = 1; j < _nA+1; j++)
			{
				for (int k = 1; k < N_objpair+1; k++)
				{
					_R_OR_common_obj_rep[ (i-1) * (_nA) * (N_objpair) + (j-1) * (N_objpair) + (k-1)] = OR_common_obj_rep[i][j][k];
					
				}
			}
		}
		for (int i = 1; i < _nO*N_attpair+1; i++)
		{
			for (int j = 1; j < 5+1; j++)
			{
				_R_report_common_att[  (i-1) * (5) + (j-1)] = report_common_att[i][j];
			}
		}
		for (int i = 1; i < _nA*N_objpair+1; i++)
		{
			for (int j = 1; j < 5+1; j++)
			{
				_R_report_common_obj[  (i-1) * (5) + (j-1)] = report_common_obj[i][j];
			}
		}

	}


	// destructor

	for (int i = 0; i < _nO+1; i++) 
	{
		for (int ii = 0; ii < _nA+1; ii++)
		{
			delete R_data[i][ii];
		}
		delete R_data[i];
	}
	delete[] R_data;
	for (int i = 0; i < _nS+1; i++)
	{
		delete R_patS[i];
	}
	delete[] R_patS;

	for (int i = 0; i < _nO+1; i++)
	{
		for (int ii = 0; ii < _nF+1; ii++)
		{
			delete R_ro_o[i][ii];
			delete R_ro_n[i][ii];
			delete R_ro_update[i][ii];
		}
		delete R_ro_o[i];
		delete R_ro_n[i];
		delete R_ro_update[i];
	}
	delete[] R_ro_o;
	delete[] R_ro_n;
	delete[] R_ro_update;
	
	for (int i = 0; i < _nA+1; i++)
	{
		delete R_ta_o[i];
		delete R_ta_n[i];
		delete R_ta_update[i];
	}
	delete[] R_ta_o;
	delete[] R_ta_n;
	delete[] R_ta_update;	
	delete[] R_ga_o;
	delete[] R_ga_n;

	for (int i = 0; i < _nT+1; i++)
	{
		delete probmat[i];
	}
	delete[] probmat;

	if (flag == true)
	{
		for (int i = 0; i < _nA+1; i++)
		{
			delete gradient_ta[i];
		}
		delete[] gradient_ta; 

		for (int i = 0; i < _nO+1; i++)
		{
			for (int ii = 0; ii < _nF+1; ii++)
			{
				delete gradient_ro[i][ii];
			}
			delete gradient_ro[i];
		}
		delete[]  gradient_ro;

		delete[] gradient_ga;

		for (int i = 0; i < _nA+1; i++)
		{
			delete se_ta[i];
		}
		delete[] se_ta;

		for (int i = 0; i < _nO+1; i++)
		{
			for (int ii = 0; ii < _nF+1; ii++)
			{
				delete se_ro[i][ii];
			}
			delete se_ro[i];
		}
		delete[] se_ro;

		delete[] se_ga;
		
		// OR_common_att_rep
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			for (int ii = 0; ii < _nO+1; ii++)
			{
				delete OR_common_att_rep[i][ii];
			}
			delete OR_common_att_rep[i];
		}
		delete[] OR_common_att_rep;
		
		// OR_common_obj_rep
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			for (int ii = 0; ii < _nA+1; ii++)
			{
				delete OR_common_obj_rep[i][ii];
			}
			delete OR_common_obj_rep[i];
		}
		delete[] OR_common_obj_rep;
				
		// report_common_att
		
		for (int i = 0; i < _nO*N_attpair+1; i++)
		{
			delete report_common_att[i];
		}
		delete[] report_common_att;
		
	       // report_common_obj	
	       for (int i = 0; i < _nA*N_objpair+1; i++)
		{
			delete report_common_obj[i];
		}
		delete[] report_common_obj;		
		
	}
}

// estimate additive model M_2
extern "C" void PlFm_X_YZ_ADD(int * _R_data, int * _R_nR, int * _R_nO, int * _R_nA, int * _R_nF, int * _R_nT, int * _R_patS, 
						  double * _R_emcrit,
						  double * _R_ro_n, double * _R_ta_n, double * _R_ga_n,
						  double * _R_ro_update, double * _R_ta_update,
						  double * _R_gradient_ro, double * _R_gradient_ta, double * _R_gradient_ga,
						  double * _R_se_ro, double * _R_se_ta, double * _R_se_ga, double * _R_delta, double * _R_logpost_n, double * _R_loglik_n, double * _R_probmat, int * _R_flag,
						  int * _R_Nbootstrap, double * _R_OR_common_att_rep, double * _R_OR_common_obj_rep, double * _R_report_common_att, double * _R_report_common_obj)
{
	// variables
	extended3 condprobx_o,condprobx_n,somega,p_o_r_t,gradient_ta,se_ta;
	extended2 probmat, margprobx_o,margprobx_n,p_r_t,gradient_ro,se_ro;
    extended1 p_r,gradient_ga,se_ga;

	binary3 R_data; int R_nR; int R_nO; int R_nA; int R_nF; int R_nT; int R_Nbootstrap; binary2 R_patS;
	long double R_emcrit; extended2 R_ro_o; extended3 R_ta_o; extended1 R_ga_o;
	extended2 R_ro_n; extended3 R_ta_n; extended1 R_ga_n; extended2 R_ro_update; extended3 R_ta_update;

	int R_flag;
	double R_delta;

	R_nR = *_R_nR; 
	R_nO = *_R_nO; 
	R_nA = *_R_nA; 
	R_nF = *_R_nF; 
	R_nT = *_R_nT;
	R_Nbootstrap = *_R_Nbootstrap;

	R_flag = *_R_flag;
	R_delta = *_R_delta;
	
	bool flag;
	if (R_flag == 0)
		flag = false;
	else
		flag = true;

	_nO = R_nO;
        _nA = R_nA;
        _nR = R_nR;
   	_nT = R_nT;
	_nF = R_nF;

	R_data = new bool **[R_nO+1]; 
	for (int i = 0; i < R_nO+1; i++) 
	{
		R_data[i] = new bool *[R_nA+1];
		for (int ii = 0; ii < R_nA+1; ii++)
		{
			R_data[i][ii] = new bool [R_nR+1];
		}
	}
	for (int i = 0; i < R_nO; i++)
	{
		for (int j = 0; j < R_nA; j++)
		{
			for (int k = 0; k < R_nR; k++)
			{
				R_data[i+1][j+1][k+1] = (bool)_R_data[i * (R_nA) * (R_nR) + j * (R_nR) + k];
			}
		}
	}

	calculate_dims();

	somega = new long double ** [_nR+1]; 
	condprobx_o = new long double ** [_nS+1]; 
	margprobx_o = new long double * [_nS+1]; 
	condprobx_n = new long double ** [_nS+1]; 
	margprobx_n = new long double * [_nS+1]; 
	probmat = new long double * [_nT+1]; 
	p_o_r_t = new long double ** [_nO+1]; 
	p_r_t = new long double * [_nR+1];
	p_r = new long double [_nR+1];
	CreateVariables( somega, condprobx_o, condprobx_n, margprobx_o, margprobx_n, probmat, p_o_r_t, p_r_t, p_r);
	

	R_patS = new bool * [_nS+1]; 
	for (int i = 0; i < _nS+1; i++)
	{
		R_patS[i] = new bool [R_nF+1];
	}
	for (int i = 0; i < _nS; i++)
	{
		for (int j = 0; j < R_nF; j++)
		{
			R_patS[i+1][j+1] = (bool)_R_patS[i * (R_nF) + j ];			
		}
	}
	R_emcrit = *_R_emcrit; 
	
	R_ro_o = new long double * [_nO+1] ; R_ro_n = new long double * [_nO+1]; R_ro_update = new long double * [_nO+1];
	R_ta_o = new long double ** [_nA+1] ; R_ta_n = new long double ** [_nA+1] ; R_ta_update = new long double ** [_nA+1] ;
	R_ga_o = new long double [_nT+1]; R_ga_n = new long double [_nT+1];
	Create_ro_ta(R_ro_o,R_ta_o);
	Create_ro_ta(R_ro_n,R_ta_n);
	Create_ro_ta(R_ro_update,R_ta_update);
	seed1= 0; 
	seed2= 1; 
	seed3= 2; 
	generate_true_parameters(R_ro_o,R_ta_o,R_ga_o);

	for (int i = 0; i < _nO; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			R_ro_n[i+1][j+1] =  _R_ro_n[  i * (_nF) + j ];
		}
	}

	for (int i = 0; i < _nO; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			R_ro_update[i+1][j+1] =  _R_ro_update[  i * (_nF) + j ];
		}
	}


	for (int i = 0; i < _nA; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			for (int k = 0; k < _nT; k++)
			{
				R_ta_n[i+1][j+1][k+1] =  _R_ta_n[  i * (_nF) * (_nT) + j * (_nT) + k];
			}

		}
	}

	for (int i = 0; i < _nA; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			for (int k = 0; k < _nT; k++)
			{
				R_ta_update[i+1][j+1][k+1] =  _R_ta_update[  i * (_nF) * (_nT) + j * (_nT) + k];
			}

		}
	}


	for (int i = 0; i < _nT; i++)
	{	
		R_ga_n[i+1] =  _R_ga_n[i];	
	}
	
		// bootstrapping variables
	extended3 OR_common_att_rep, OR_common_obj_rep;
	extended2 report_common_att, report_common_obj;
	int N_attpair, N_objpair;
	N_attpair = (_nA*(_nA-1)) / 2;
	N_objpair = (_nO*(_nO-1)) / 2;
	
	

	if (flag == true)
	{

		gradient_ta = new long double ** [_nA+1]; 
		for (int i = 0; i < _nA+1; i++)
		{
			gradient_ta[i] = new long double * [_nF+1];
			for (int ii = 0; ii < _nF+1; ii++)
			{
				gradient_ta[i][ii] = new long double [_nT+1];
			}
		}


		gradient_ro = new long double * [_nO+1]; 
		for (int i = 0; i < _nO+1; i++)
		{
			gradient_ro[i] = new long double  [_nF+1];
		}


		gradient_ga = new long double [_nT+1];



		se_ta = new long double ** [_nA+1]; 
		for (int i = 0; i < _nA+1; i++)
		{
			se_ta[i] = new long double *[_nF+1];
			for (int ii = 0; ii < _nF+1; ii++)
			{
				se_ta[i][ii] = new long double [_nT+1];
			}
		}


		se_ro = new long double * [_nO+1]; 
		for (int i = 0; i < _nO+1; i++)
		{
			se_ro[i] = new long double [_nF+1];
		}


		se_ga = new long double [_nT+1];
		
		
		
		// bootstrapping variables

		// construct OR_common_att_rep,OR_common_obj_rep
		OR_common_att_rep = new long double ** [_Nbootstrap+1] ;
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			OR_common_att_rep[i] = new long double * [_nO+1];
			for (int ii = 0; ii < _nO+1; ii++)
			{
				OR_common_att_rep[i][ii] = new long double [N_attpair+1];
			}
		}
		OR_common_obj_rep = new long double ** [_Nbootstrap+1] ;
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			OR_common_obj_rep[i] = new long double * [_nA+1];
			for (int ii = 0; ii < _nA+1; ii++)
			{
				OR_common_obj_rep[i][ii] = new long double [N_objpair+1];
			}
		}
		//construct report_common_att & report_common_obj
		report_common_att = new long double * [_nO*N_attpair+1] ;
		for (int i = 0; i < _nO*N_attpair+1; i++)
		{
			report_common_att[i] = new long double [5+1];
		}
		report_common_obj = new long double * [_nA*N_objpair+1] ;
		for (int i = 0; i < _nA*N_objpair+1; i++)
		{
			report_common_obj[i] = new long double [5+1];
		}

		
		
		
	}
		else
	{
		se_ta = NULL;
		se_ga = NULL;
		se_ro = NULL;
		gradient_ta = NULL; 
		gradient_ga = NULL; 
		gradient_ro = NULL;
		OR_common_att_rep = NULL;
		OR_common_obj_rep = NULL;
		report_common_att = NULL;
		report_common_obj = NULL; 
	}

	//////////////////// call C function model M_2///////////////////////
	PlFm_X_YZ_intern_ADD(R_data,R_nR,R_nO,R_nA,R_nF,R_nT,R_patS,R_emcrit,R_ro_o,R_ta_o,R_ga_o,R_ro_n,R_ta_n,R_ga_n,R_ro_update,R_ta_update,flag,R_delta,
		margprobx_o,margprobx_n,somega,gradient_ro, se_ro, p_o_r_t,
		probmat,condprobx_o,condprobx_n,gradient_ta, se_ta, p_r_t,
		gradient_ga,p_r,se_ga,R_Nbootstrap, OR_common_att_rep, OR_common_obj_rep, report_common_att, report_common_obj);
	////////////////////////////////////////////////////////////



	for (int i = 1; i < _nO+1; i++)
	{
		for (int j = 1; j < _nF+1; j++)
		{
			_R_ro_n[  (i-1) * (_nF)  + (j-1) ] = R_ro_n[i][j];
		}
	}
	for (int i = 1; i < _nA+1; i++)
	{
		for (int j = 1; j < _nF+1; j++)
		{
			for (int k = 1; k < _nT+1; k++)
			{
				_R_ta_n[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = R_ta_n[i][j][k];	
			}	
		}
	}
	for (int i = 1; i < _nT+1; i++)
	{
		_R_ga_n[i-1] = R_ga_n[i];		
	}

	*_R_logpost_n = logpost_n;
	*_R_loglik_n = loglik_n;


	if (flag == true)
	{
		for (int i = 1; i < _nA+1; i++)
		{
			for (int j = 1; j < _nF+1; j++)
			{
				for (int k = 1; k < _nT+1; k++)
				{
					_R_gradient_ta[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = gradient_ta[i][j][k];
					_R_se_ta[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = se_ta[i][j][k];
				}	
			}
		}

		
		for (int i = 1; i < _nO+1; i++)
		{
			for (int j = 1; j < _nF+1; j++)
			{
				_R_gradient_ro[  (i-1) * (_nF) + (j-1) ] = gradient_ro[i][j];
				_R_se_ro[  (i-1) * (_nF)  + (j-1) ] = se_ro[i][j];
			}
		}

		for (int i = 1; i < _nT+1; i++)
		{
			_R_gradient_ga[i-1] = gradient_ga[i];
			_R_se_ga[i-1] = se_ga[i];
		}

		for (int i = 1; i < _nT+1; i++)
		{
			for (int j = 1; j < _nR+1; j++)
			{
				_R_probmat[ (i-1) * (_nR) + (j-1)] = probmat[i][j];
			}
		}
		
		
			// bootstrapping variables
		for (int i = 1; i < _Nbootstrap+1; i++)
		{
			for (int j = 1; j < _nO+1; j++)
			{
				for (int k = 1; k < N_attpair+1; k++)
				{
					_R_OR_common_att_rep[ (i-1) * (_nO) * (N_attpair) + (j-1) * (N_attpair) + (k-1)] = OR_common_att_rep[i][j][k];
					
				}
			}
		}
		for (int i = 1; i < _Nbootstrap+1; i++)
		{
			for (int j = 1; j < _nA+1; j++)
			{
				for (int k = 1; k < N_objpair+1; k++)
				{
					_R_OR_common_obj_rep[ (i-1) * (_nA) * (N_objpair) + (j-1) * (N_objpair) + (k-1)] = OR_common_obj_rep[i][j][k];
					
				}
			}
		}
		for (int i = 1; i < _nO*N_attpair+1; i++)
		{
			for (int j = 1; j < 5+1; j++)
			{
				_R_report_common_att[  (i-1) * (5) + (j-1)] = report_common_att[i][j];
			}
		}
		for (int i = 1; i < _nA*N_objpair+1; i++)
		{
			for (int j = 1; j < 5+1; j++)
			{
				_R_report_common_obj[  (i-1) * (5) + (j-1)] = report_common_obj[i][j];
			}
		}

		

	}


	// destructor

	for (int i = 0; i < _nO+1; i++) 
	{
		for (int ii = 0; ii < _nA+1; ii++)
		{
			delete R_data[i][ii];
		}
		delete R_data[i];
	}
	delete[] R_data;
	for (int i = 0; i < _nS+1; i++)
	{
		delete R_patS[i];
	}
	delete[] R_patS;
	for (int i = 0; i < _nO+1; i++)
	{
		delete R_ro_o[i];
		delete R_ro_n[i];
		delete R_ro_update[i];
	}
	delete[] R_ro_o;
	delete[] R_ro_n;
	delete[] R_ro_update;
	
	for (int i = 0; i < _nA+1; i++)
	{
		for (int ii = 0; ii < _nF+1; ii++)
		{
			delete R_ta_o[i][ii];
			delete R_ta_n[i][ii];
			delete R_ta_update[i][ii];
		}
		delete R_ta_o[i];
		delete R_ta_n[i];
		delete R_ta_update[i];
	}
	delete[]  R_ta_o;
	delete[]  R_ta_n;
	delete[]  R_ta_update;
	delete[]  R_ga_o;
	delete[]  R_ga_n;

	for (int i = 0; i < _nT+1; i++)
	{
		delete probmat[i];
	}
	delete[] probmat;

	if (flag == true)
	{
		for (int i = 0; i < _nA+1; i++)
		{
			for (int ii = 0; ii < _nF+1; ii++)
			{
				delete gradient_ta[i][ii];
			}
			delete gradient_ta[i];
		}
		delete[] gradient_ta; 

		for (int i = 0; i < _nO+1; i++)
		{
			delete gradient_ro[i];
		}
		delete[] gradient_ro;

		delete[] gradient_ga;

		for (int i = 0; i < _nA+1; i++)
		{
			for (int ii = 0; ii < _nF+1; ii++)
			{
				delete se_ta[i][ii];
			}
			delete se_ta[i];
		}
		delete[] se_ta;

		for (int i = 0; i < _nO+1; i++)
		{
			delete se_ro[i];
		}
		delete[] se_ro;

		delete[] se_ga;
		
		
		
		
		// OR_common_att_rep
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			for (int ii = 0; ii < _nO+1; ii++)
			{
				delete OR_common_att_rep[i][ii];
			}
			delete OR_common_att_rep[i];
		}
		delete[] OR_common_att_rep;
		
		// OR_common_obj_rep
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			for (int ii = 0; ii < _nA+1; ii++)
			{
				delete OR_common_obj_rep[i][ii];
			}
			delete OR_common_obj_rep[i];
		}
		delete[] OR_common_obj_rep;
				
		// report_common_att
		
		for (int i = 0; i < _nO*N_attpair+1; i++)
		{
			delete report_common_att[i];
		}
		delete[] report_common_att;
		
	       // report_common_obj	
	        for (int i = 0; i < _nA*N_objpair+1; i++)
		{
			delete report_common_obj[i];
		}
		delete[] report_common_obj;		
		
	
	}
}

// estimate additive model M_3
extern "C" void PlFm_XZ_YZ_ADD(int * _R_data, int * _R_nR, int * _R_nO, int * _R_nA, int * _R_nF, int * _R_nT, int * _R_patS, 
						  double * _R_emcrit,
						  double * _R_ro_n, double * _R_ta_n, double * _R_ga_n,
						  double * _R_ro_update, double * _R_ta_update,
						  double * _R_gradient_ro, double * _R_gradient_ta, double * _R_gradient_ga,
						  double * _R_se_ro, double * _R_se_ta, double * _R_se_ga, double * _R_delta, double * _R_logpost_n, double * _R_loglik_n, double * _R_probmat, int * _R_flag,
						  int * _R_Nbootstrap, double * _R_OR_common_att_rep, double * _R_OR_common_obj_rep, double * _R_report_common_att, double * _R_report_common_obj)
{
	// variables
	extended3 condprobx_o,condprobx_n,somega,margprobx_o,margprobx_n,p_o_r_t,gradient_ta,se_ta,gradient_ro,se_ro;
	extended2 probmat, p_r_t;
    extended1 p_r,gradient_ga,se_ga;

	binary3 R_data; int R_nR; int R_nO; int R_nA; int R_nF; int R_nT; int R_Nbootstrap; binary2 R_patS;
	long double R_emcrit; extended3 R_ro_o; extended3 R_ta_o; extended1 R_ga_o;
	extended3 R_ro_n; extended3 R_ta_n; extended1 R_ga_n; extended3 R_ro_update; extended3 R_ta_update;

	int R_flag;
	double R_delta;

	R_nR = *_R_nR; 
	R_nO = *_R_nO; 
	R_nA = *_R_nA; 
	R_nF = *_R_nF; 
	R_nT = *_R_nT;
	R_Nbootstrap = *_R_Nbootstrap;

	R_flag = *_R_flag;
	R_delta = *_R_delta;
	
	bool flag;
	if (R_flag == 0)
		flag = false;
	else
		flag = true;

	_nO = R_nO;
    _nA = R_nA;
    _nR = R_nR;
	_nT = R_nT;
	_nF = R_nF;

	R_data = new bool **[R_nO+1]; 
	for (int i = 0; i < R_nO+1; i++) 
	{
		R_data[i] = new bool *[R_nA+1];
		for (int ii = 0; ii < R_nA+1; ii++)
		{
			R_data[i][ii] = new bool [R_nR+1];
		}
	}
	for (int i = 0; i < R_nO; i++)
	{
		for (int j = 0; j < R_nA; j++)
		{
			for (int k = 0; k < R_nR; k++)
			{
				R_data[i+1][j+1][k+1] = (bool)_R_data[i * (R_nA) * (R_nR) + j * (R_nR) + k];
			}
		}
	}


	calculate_dims();

	somega = new long double ** [_nR+1]; 
	condprobx_o = new long double ** [_nS+1]; 
	margprobx_o = new long double ** [_nS+1]; 
	condprobx_n = new long double ** [_nS+1]; 
	margprobx_n = new long double ** [_nS+1]; 
	probmat = new long double * [_nT+1]; 
	p_o_r_t = new long double ** [_nO+1]; 
	p_r_t = new long double * [_nR+1];
	p_r = new long double [_nR+1];
	CreateVariables( somega, condprobx_o, condprobx_n, margprobx_o, margprobx_n, probmat, p_o_r_t, p_r_t, p_r);
	

	R_patS = new bool * [_nS+1]; 
	for (int i = 0; i < _nS+1; i++)
	{
		R_patS[i] = new bool [R_nF+1];
	}
	for (int i = 0; i < _nS; i++)
	{
		for (int j = 0; j < R_nF; j++)
		{
			R_patS[i+1][j+1] = (bool)_R_patS[i * (R_nF) + j ];			
		}
	}
	R_emcrit = *_R_emcrit; 
	
	R_ro_o = new long double ** [_nO+1] ; R_ro_n = new long double ** [_nO+1]; R_ro_update = new long double ** [_nO+1];
	R_ta_o = new long double ** [_nA+1] ; R_ta_n = new long double ** [_nA+1] ; R_ta_update = new long double ** [_nA+1] ;
	R_ga_o = new long double [_nT+1]; R_ga_n = new long double [_nT+1];
	Create_ro_ta(R_ro_o,R_ta_o);
	Create_ro_ta(R_ro_n,R_ta_n);
	Create_ro_ta(R_ro_update,R_ta_update);
	seed1= 0; 
	seed2= 1; 
	seed3= 2; 
	generate_true_parameters(R_ro_o,R_ta_o,R_ga_o);

	for (int i = 0; i < _nO; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			for (int k = 0; k < _nT; k++)
			{

				R_ro_n[i+1][j+1][k+1] =  _R_ro_n[  i * (_nF) * (_nT) + j * (_nT) + k];
			}
		}
	}

	for (int i = 0; i < _nO; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			for (int k = 0; k < _nT; k++)
			{

				R_ro_update[i+1][j+1][k+1] =  _R_ro_update[  i * (_nF) * (_nT) + j * (_nT) + k];
			}
		}
	}


	for (int i = 0; i < _nA; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			for (int k = 0; k < _nT; k++)
			{
				R_ta_n[i+1][j+1][k+1] =  _R_ta_n[  i * (_nF) * (_nT) + j * (_nT) + k];
			}

		}
	}

	for (int i = 0; i < _nA; i++)
	{
		for (int j = 0; j < _nF; j++)
		{
			for (int k = 0; k < _nT; k++)
			{
				R_ta_update[i+1][j+1][k+1] =  _R_ta_update[  i * (_nF) * (_nT) + j * (_nT) + k];
			}

		}
	}

	for (int i = 0; i < _nT; i++)
	{	
		R_ga_n[i+1] =  _R_ga_n[i];	
	}
	
	// bootstrapping variables
	extended3 OR_common_att_rep, OR_common_obj_rep;
	extended2 report_common_att, report_common_obj;
	int N_attpair, N_objpair;
	N_attpair = (_nA*(_nA-1)) / 2;
	N_objpair = (_nO*(_nO-1)) / 2;

	if (flag == true)
	{

		gradient_ta = new long double ** [_nA+1]; 
		for (int i = 0; i < _nA+1; i++)
		{
			gradient_ta[i] = new long double * [_nF+1];
			for (int ii = 0; ii < _nF+1; ii++)
			{
				gradient_ta[i][ii] = new long double [_nT+1];
			}
		}


		gradient_ro = new long double ** [_nO+1]; 
		for (int i = 0; i < _nO+1; i++)
		{
			gradient_ro[i] = new long double * [_nF+1];
			for (int ii = 0; ii < _nF+1; ii++)
			{
				gradient_ro[i][ii] = new long double [_nT+1];
			}
		}


		gradient_ga = new long double [_nT+1];



		se_ta = new long double ** [_nA+1]; 
		for (int i = 0; i < _nA+1; i++)
		{
			se_ta[i] = new long double *[_nF+1];
			for (int ii = 0; ii < _nF+1; ii++)
			{
				se_ta[i][ii] = new long double [_nT+1];
			}
		}


		se_ro = new long double ** [_nO+1]; 
		for (int i = 0; i < _nO+1; i++)
		{
			se_ro[i] = new long double * [_nF+1];
			for (int ii = 0; ii < _nF+1; ii++)
			{
				se_ro[i][ii] = new long double [_nT+1];
			}
		}


		se_ga = new long double [_nT+1];
		
		
		// bootstrapping variables
		

		// construct OR_common_att_rep,OR_common_obj_rep
		OR_common_att_rep = new long double ** [_Nbootstrap+1] ;
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			OR_common_att_rep[i] = new long double * [_nO+1];
			for (int ii = 0; ii < _nO+1; ii++)
			{
				OR_common_att_rep[i][ii] = new long double [N_attpair+1];
			}
		}
		OR_common_obj_rep = new long double ** [_Nbootstrap+1] ;
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			OR_common_obj_rep[i] = new long double * [_nA+1];
			for (int ii = 0; ii < _nA+1; ii++)
			{
				OR_common_obj_rep[i][ii] = new long double [N_objpair+1];
			}
		}
		//construct report_common_att & report_common_obj
		report_common_att = new long double * [_nO*N_attpair+1] ;
		for (int i = 0; i < _nO*N_attpair+1; i++)
		{
			report_common_att[i] = new long double [5+1];
		}
		report_common_obj = new long double * [_nA*N_objpair+1] ;
		for (int i = 0; i < _nA*N_objpair+1; i++)
		{
			report_common_obj[i] = new long double [5+1];
		}

		
	}
		else
	{
		se_ta = NULL;
		se_ga = NULL;
		se_ro = NULL;
		gradient_ta = NULL; 
		gradient_ga = NULL; 
		gradient_ro = NULL;
		OR_common_att_rep = NULL;
		OR_common_obj_rep = NULL;
		report_common_att = NULL;
		report_common_obj = NULL; 
	}

	//////////////////// call C function model M_3 ///////////////////////
	PlFm_XZ_YZ_intern_ADD(R_data,R_nR,R_nO,R_nA,R_nF,R_nT,R_patS,R_emcrit,R_ro_o,R_ta_o,R_ga_o,R_ro_n,R_ta_n,R_ga_n,R_ro_update,R_ta_update,flag,R_delta,
		margprobx_o,margprobx_n,somega,gradient_ro, se_ro, p_o_r_t,
		probmat,condprobx_o,condprobx_n,gradient_ta, se_ta, p_r_t,
		gradient_ga,p_r,se_ga,R_Nbootstrap, OR_common_att_rep, OR_common_obj_rep, report_common_att, report_common_obj);
	////////////////////////////////////////////////////////////



	for (int i = 1; i < _nO+1; i++)
	{
		for (int j = 1; j < _nF+1; j++)
		{
			for (int k = 1; k < _nT+1; k++)
			{

				_R_ro_n[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = R_ro_n[i][j][k];
			}
		}
	}
	for (int i = 1; i < _nA+1; i++)
	{
		for (int j = 1; j < _nF+1; j++)
		{
			for (int k = 1; k < _nT+1; k++)
			{
				_R_ta_n[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = R_ta_n[i][j][k];	
			}	
		}
	}
	for (int i = 1; i < _nT+1; i++)
	{
		_R_ga_n[i-1] = R_ga_n[i];		
	}

	*_R_logpost_n = logpost_n;
	*_R_loglik_n = loglik_n;


	if (flag == true)
	{
		for (int i = 1; i < _nA+1; i++)
		{
			for (int j = 1; j < _nF+1; j++)
			{
				for (int k = 1; k < _nT+1; k++)
				{
					_R_gradient_ta[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = gradient_ta[i][j][k];
					_R_se_ta[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = se_ta[i][j][k];
				}	
			}
		}

		
		for (int i = 1; i < _nO+1; i++)
		{
			for (int j = 1; j < _nF+1; j++)
			{
				for (int k = 1; k < _nT+1; k++)
				{
					_R_gradient_ro[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = gradient_ro[i][j][k];
					_R_se_ro[  (i-1) * (_nF) * (_nT) + (j-1) * (_nT) + (k-1)] = se_ro[i][j][k];
				}
			}
		}

		for (int i = 1; i < _nT+1; i++)
		{
			_R_gradient_ga[i-1] = gradient_ga[i];
			_R_se_ga[i-1] = se_ga[i];
		}

		for (int i = 1; i < _nT+1; i++)
		{
			for (int j = 1; j < _nR+1; j++)
			{
				_R_probmat[ (i-1) * (_nR) + (j-1)] = probmat[i][j];
			}
		}
		
		
		// bootstrapping variables
		for (int i = 1; i < _Nbootstrap+1; i++)
		{
			for (int j = 1; j < _nO+1; j++)
			{
				for (int k = 1; k < N_attpair+1; k++)
				{
					_R_OR_common_att_rep[ (i-1) * (_nO) * (N_attpair) + (j-1) * (N_attpair) + (k-1)] = OR_common_att_rep[i][j][k];
					
				}
			}
		}
		for (int i = 1; i < _Nbootstrap+1; i++)
		{
			for (int j = 1; j < _nA+1; j++)
			{
				for (int k = 1; k < N_objpair+1; k++)
				{
					_R_OR_common_obj_rep[ (i-1) * (_nA) * (N_objpair) + (j-1) * (N_objpair) + (k-1)] = OR_common_obj_rep[i][j][k];
					
				}
			}
		}
		for (int i = 1; i < _nO*N_attpair+1; i++)
		{
			for (int j = 1; j < 5+1; j++)
			{
				_R_report_common_att[  (i-1) * (5) + (j-1)] = report_common_att[i][j];
			}
		}
		for (int i = 1; i < _nA*N_objpair+1; i++)
		{
			for (int j = 1; j < 5+1; j++)
			{
				_R_report_common_obj[  (i-1) * (5) + (j-1)] = report_common_obj[i][j];
			}
		}
		

	}


	// destructor

	for (int i = 0; i < _nO+1; i++) 
	{
		for (int ii = 0; ii < _nA+1; ii++)
		{
			delete R_data[i][ii];
		}
		delete R_data[i];
	}
	delete[] R_data;
	for (int i = 0; i < _nS+1; i++)
	{
		delete R_patS[i];
	}
	delete[] R_patS;
	for (int i = 0; i < _nO+1; i++)
	{
		for (int ii = 0; ii < _nF+1; ii++)
		{
			delete R_ro_o[i][ii];
			delete R_ro_n[i][ii];
			delete R_ro_update[i][ii];
		}
		delete R_ro_o[i];
		delete R_ro_n[i];
		delete R_ro_update[i];
	}
	delete[] R_ro_o;
	delete[] R_ro_n;
	delete[] R_ro_update;
	for (int i = 0; i < _nA+1; i++)
	{
		for (int ii = 0; ii < _nF+1; ii++)
		{
			delete R_ta_o[i][ii];
			delete R_ta_n[i][ii];
			delete R_ta_update[i][ii];
		}
		delete R_ta_o[i];
		delete R_ta_n[i];
		delete R_ta_update[i];
	}
	delete[]  R_ta_o;
	delete[]  R_ta_n;
	delete[]  R_ta_update;
	delete[]  R_ga_o;
	delete[]  R_ga_n;

	for (int i = 0; i < _nT+1; i++)
	{
		delete probmat[i];
	}
	delete[] probmat;

	if (flag == true)
	{
		for (int i = 0; i < _nA+1; i++)
		{
			for (int ii = 0; ii < _nF+1; ii++)
			{
				delete gradient_ta[i][ii];
			}
			delete gradient_ta[i];
		}
		delete[] gradient_ta; 

		for (int i = 0; i < _nO+1; i++)
		{
			for (int ii = 0; ii < _nF+1; ii++)
			{
				delete gradient_ro[i][ii];
			}
			delete gradient_ro[i];
		}
		delete[] gradient_ro;

		delete[] gradient_ga;

		for (int i = 0; i < _nA+1; i++)
		{
			for (int ii = 0; ii < _nF+1; ii++)
			{
				delete se_ta[i][ii];
			}
			delete se_ta[i];
		}
		delete[] se_ta;

		for (int i = 0; i < _nO+1; i++)
		{
			for (int ii = 0; ii < _nF+1; ii++)
			{
				delete se_ro[i][ii];
			}
			delete se_ro[i];
		}
		delete[] se_ro;

		delete[] se_ga;
		
		
	// OR_common_att_rep
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			for (int ii = 0; ii < _nO+1; ii++)
			{
				delete OR_common_att_rep[i][ii];
			}
			delete OR_common_att_rep[i];
		}
		delete[] OR_common_att_rep;
		
		// OR_common_obj_rep
		for (int i = 0; i < _Nbootstrap+1; i++)
		{
			for (int ii = 0; ii < _nA+1; ii++)
			{
				delete OR_common_obj_rep[i][ii];
			}
			delete OR_common_obj_rep[i];
		}
		delete[] OR_common_obj_rep;
				
		// report_common_att
		
		for (int i = 0; i < _nO*N_attpair+1; i++)
		{
			delete report_common_att[i];
		}
		delete[] report_common_att;
		
	   // report_common_obj	
	   for (int i = 0; i < _nA*N_objpair+1; i++)
		{
			delete report_common_obj[i];
		}
		delete[] report_common_obj;			
	
	}
}



long double pow_1(long double grondtal, binary exponent)
{
	if (exponent==0)    
		return (1-grondtal); 
	else 
		return grondtal;
}

long double pow_2(long double grondtal, binary exponent)
{
	if (exponent==0)   
		return 1;
	else 
		return grondtal;
}

long double pow_3(unsigned short int x, unsigned short int a, extended2 condprobx, binary exponent)
{
	// XZ_Y
	if ((exponent==1) &&  (x==1))    return 1.0;
	else if ((exponent==1) &&  (x>1))  return condprobx[x][a];
	else if (exponent==0)  return 1.0 -condprobx[x][a];
	return 0.0;
}

long double pow_3(unsigned short int x, unsigned short int a,unsigned short int t, extended3 condprobx, binary exponent)
{
	// X_YZ
	if ((exponent==1) &&  (x==1))    return 1.0;
	else if ((exponent==1) &&  (x>1))   return condprobx[x][a][t];
	else if (exponent==0)    return 1.0 -condprobx[x][a][t];
	return 0.0;
}

// function to compute loglikelihood model M_1
long double loglikelihood (extended1 ga, extended2 condprobx,extended3 margprobx)
{

	unsigned short int o,a,r,t,s; 
	long double totsum,sumpattern,sumclas,prod,prodobject;
  
	totsum= 0;
	for( r= 1; r<=_nR; r++)  
	{ 
		sumclas= 0;
		for( t= 1; t<= _nT; t++)  
		{ 
			prodobject= 1;
			for( o= 1; o<=_nO; o++)  
			{ 
				sumpattern= 0;
				for( s= 1; s<= _nS; s++)  
				{ 
					prod= 1;
					for( a= 1; a<=_nA; a++)  
					{ 
						prod= prod*pow_1(condprobx[s][a],dataB[o][a][r]);
					} 
					prod= prod*margprobx[s][o][t];
					sumpattern= sumpattern+prod;
				} 
				prodobject= prodobject*sumpattern;
			} 
			sumclas= sumclas+prodobject*ga[t];
		} 
		totsum= totsum+log(sumclas);
	} 
	return totsum;
}

// function to compute loglikelihood model M_2
long double loglikelihood (extended1 ga, extended3 condprobx,extended2 margprobx)
{

	unsigned short int o,a,r,t,s; 
	long double totsum,sumpattern,sumclas,prod,prodobject;
  
	totsum= 0;
	for( r= 1; r<=_nR; r++)  { 
	sumclas= 0;
	for( t= 1; t<=_nT; t++)  { 
	prodobject= 1;
	for( o= 1; o<=_nO; o++)  { 
		sumpattern= 0;
		for( s= 1; s<=_nS; s++)  { 
		prod= 1;
		for( a= 1; a<=_nA; a++)  { 
		prod= prod*pow_1(condprobx[s][a][t],dataB[o][a][r]);
		} 
		prod= prod*margprobx[s][o];
		sumpattern= sumpattern+prod;
		} 
	prodobject= prodobject*sumpattern;
		} 
	sumclas= sumclas+prodobject*ga[t];
	} 
	totsum= totsum+log(sumclas);
	} 
	return totsum;

}

// function to compute loglikelihood model M_3
long double loglikelihood (extended1 ga, extended3 condprobx,extended3 margprobx)
{

	unsigned short int o,a,r,t,s; 
	long double totsum,sumpattern,sumclas,prod,prodobject;
  
	totsum= 0;
	for( r= 1; r<=_nR; r++)  { 
	sumclas= 0;
	for( t= 1; t<=_nT; t++)  { 
	prodobject= 1;
	for( o= 1; o<=_nO; o++)  { 
		sumpattern= 0;
		for( s= 1; s<=_nS; s++)  { 
		prod= 1;
		for( a= 1; a<=_nA; a++)  { 
		prod= prod*pow_1(condprobx[s][a][t],dataB[o][a][r]);
		} 
		prod= prod*margprobx[s][o][t];
		sumpattern= sumpattern+prod;
		} 
	prodobject= prodobject*sumpattern;
		} 
	sumclas= sumclas+prodobject*ga[t];
	} 
	totsum= totsum+log(sumclas);
	} 
	return totsum;


}

// function to compute log posterior density model M_1
long double logposterior (extended3 ro, extended2 ta, extended1 ga, extended2 condprobx, extended3 margprobx, extended3 ro_update, extended2 ta_update)
{

	unsigned short int o,a,r,t,f,s;	
	long double totsum,sumpattern,sumclas,prod,prodobject,logprior;

	totsum= 0;
	for( r= 1; r<=_nR; r++)  
	{ 
		sumclas= 0;
		for( t= 1; t<=_nT; t++)  
		{ 
			prodobject= 1;
			for( o= 1; o<=_nO; o++)  
			{ 
				sumpattern= 0;
				for( s= 1; s<=_nS; s++)  
				{ 
					prod= 1;
					for( a= 1; a<=_nA; a++)  
					{ 
						prod= prod*pow_1(condprobx[s][a],dataB[o][a][r]);
					} 
					prod= prod*margprobx[s][o][t];
					sumpattern= sumpattern+prod;
				} 
				prodobject= prodobject*sumpattern;
			} 
			sumclas= sumclas+prodobject*ga[t];
		} 
		totsum= totsum+log(sumclas);
	} 

	logprior= 0;
	for( t= 1; t<=_nT; t++) logprior= logprior+(2.0/((long double)_nT))*log(ga[t]);

	for( o= 1; o<=_nO; o++)  
	{ 
		for( f= 1; f<=_nF; f++)  
		{ 
			for( t= 1; t<=_nT; t++)  
			{ 
                if (ro_update[o][f][t]==1)
                { 
				   logprior= logprior+(1.0/(((long double)_nO)*((long double)_nT)))*log(ro[o][f][t])+(1.0/(((long double)_nO)*((long double)_nT)))*log(1-ro[o][f][t]);
                 }
			} 
		} 
	} 

	for( a= 1; a<=_nA; a++)  
	{ 
		for( f= 1; f<=_nF; f++)  
		{
            if (ta_update[a][f]==1)
            {  
               logprior= logprior+(1.0/((long double)_nA))*log(ta[a][f])+(1.0/((long double)_nA))*log(1-ta[a][f]);
            }
		} 
	} 

	return totsum+logprior; 
	
}

// function to compute log posterior density model M_2
long double logposterior (extended2 ro, extended3 ta, extended1 ga, extended3 condprobx, extended2 margprobx, extended2 ro_update, extended3 ta_update)
{

	unsigned short int o,a,r,t,f,s;	
	long double totsum,sumpattern,sumclas,prod,prodobject,logprior;

	totsum= 0;
	for( r= 1; r<=_nR; r++)  { 
	sumclas= 0;
	for( t= 1; t<=_nT; t++)  { 
	prodobject= 1;
	for( o= 1; o<=_nO; o++)  { 
		sumpattern= 0;
		for( s= 1; s<=_nS; s++)  { 
		prod= 1;
		for( a= 1; a<=_nA; a++)  { 
		prod= prod*pow_1(condprobx[s][a][t],dataB[o][a][r]);
		} 
		prod= prod*margprobx[s][o];
		sumpattern= sumpattern+prod;
		} 
	prodobject= prodobject*sumpattern;
		} 
	sumclas= sumclas+prodobject*ga[t];
	} 
	totsum= totsum+log(sumclas);
	} 

	logprior= 0;
	for( t= 1; t<=_nT; t++) logprior= logprior+(2.0/((long double)_nT))*log(ga[t]); // same as Latent Gold

	for( o= 1; o<=_nO; o++)  
      { 
	     for( f= 1; f<=_nF; f++)  
          { 
            if (ro_update[o][f]==1)
            {  
               logprior= logprior+(1.0/((long double)_nO))*log(ro[o][f])+(1.0/((long double)_nO))*log(1-ro[o][f]);  // same as Latent Gold
            }
          } 
      } 

	for( a= 1; a<=_nA; a++)  { 
	for( f= 1; f<=_nF; f++)  { 
	for( t= 1; t<=_nT; t++)  { 
    if (ta_update[a][f][t]==1)
    {     
	logprior= logprior+(1.0/(((long double)_nT)*((long double)_nA)))*log(ta[a][f][t])+(1.0/(((long double)_nT)*((long double)_nA)))*log(1-ta[a][f][t]);  // same as in Latent Gold
    }  
    } 
	} 
	} 

	return ( totsum+logprior);

	
}

// function to compute log posterior density model M_3
long double logposterior (extended3 ro, extended3 ta, extended1 ga, extended3 condprobx, extended3 margprobx, extended3 ro_update, extended3 ta_update)
{

	unsigned short int o,a,r,t,f,s;	
	long double totsum,sumpattern,sumclas,prod,prodobject,logprior;

	totsum= 0;
	for( r= 1; r<=_nR; r++)  { 
	sumclas= 0;
	for( t= 1; t<=_nT; t++)  { 
	prodobject= 1;
	for( o= 1; o<=_nO; o++)  { 
		sumpattern= 0;
		for( s= 1; s<=_nS; s++)  { 
		prod= 1;
		for( a= 1; a<=_nA; a++)  { 
		prod= prod*pow_1(condprobx[s][a][t],dataB[o][a][r]);
		} 
		prod= prod*margprobx[s][o][t];
		sumpattern= sumpattern+prod;
		} 
	prodobject= prodobject*sumpattern;
		} 
	sumclas= sumclas+prodobject*ga[t];
	} 
	totsum= totsum+log(sumclas);
	} 

	logprior= 0;
	for( t= 1; t<=_nT; t++) logprior= logprior+(2.0/((long double)_nT))*log(ga[t]); // same prior as Latent Gold

	for( o= 1; o<=_nO; o++)  { 
	for( f= 1; f<=_nF; f++)  { 
	for( t= 1; t<=_nT; t++)  {
    if (ro_update[o][f][t]==1)
    {      
	logprior= logprior+(1.0/((long double)(_nO*_nT)))*log(ro[o][f][t])+(1.0/((long double)(_nO*_nT)))*log(1-ro[o][f][t]);   // same prior as latent gold
    }
    } 
	} 
	} 

	for( a= 1; a<=_nA; a++)  { 
	for( f= 1; f<=_nF; f++)  { 
	for( t= 1; t<=_nT; t++)  { 
    if (ta_update[a][f][t]==1)
    { 
	logprior= logprior+(1.0/((long double)(_nA*_nT)))*log(ta[a][f][t])+(1.0/((long double)(_nA*_nT)))*log(1.0-ta[a][f][t]); //same prior as latent gold
    }
	} 
	} 
	} 

	return (totsum+logprior);


	
}

// function to compute number of data latent data patterns for an object/attribute 
// as a function of the number of latent features
void calculate_dims(void)
{
	if (_nF==1   )
	{ 
	_nS= 2;
	}  
	else if (_nF==2   )
	{ 
	_nS= 4;
	}  
	else if (_nF==3   )
	{ 
	_nS= 8;
	}  
	else if (_nF==4   )
	{ 
	_nS= 16;
	}  
	else if (_nF==5   )
	{ 
	_nS= 32;
	}  
	else if (_nF==6   )
	{ 
	_nS= 64;
	}  
	else if (_nF==7   )
	{ 
	_nS= 128;
	}
    else if (_nF==8   )
	{ 
	_nS= 256;
	}
    else if (_nF==9   )
	{ 
	_nS= 512;
	} 
    else if (_nF==10   )
	{ 
	_nS= 1024;
	} 
}

// generate parameters model M_1
void generate_true_parameters(extended3 ro, extended2 ta, extended1 ga)
{

	unsigned short int o,a,t,f;
	long double     sum;
  
	for( o= 1; o<=_nO; o++)  
	{ 
		for( f= 1; f<=_nF; f++)  
		{ 
			for( t= 1; t<=_nT; t++)  
			{ 
				ro[o][f][t]= runif(seed1,seed2,seed3);
			} 
		} 
	} 
	for( a= 1; a<=_nA; a++)  
	{ 
		for( f= 1; f<=_nF; f++)  
		{ 
			ta[a][f]= runif(seed1,seed2,seed3);
		} 
	} 
	for( t= 1; t<=_nT; t++) 
		ga[t]= 4+runif(seed1,seed2,seed3);
	sum= 0;
	for( t= 1; t<=_nT; t++) 
		sum= sum+ga[t];
	for( t= 1; t<=_nT; t++) 
		ga[t]= ga[t]/sum;
}

// generate parameters model M_2
void generate_true_parameters(extended2 ro, extended3 ta, extended1 ga)
{

	unsigned short int o,a,t,f;
	long double     sum;
  
	for( o= 1; o<=_nO; o++)  
	{ 
		for( f= 1; f<=_nF; f++)  
		{ 
			ro[o][f]= runif(seed1,seed2,seed3);
		} 
	} 
	for( a= 1; a<=_nA; a++)  
	{ 
		for( f= 1; f<=_nF; f++)  
		{ 
			for( t= 1; t<=_nT; t++)  
			{ 
				ta[a][f][t]= runif(seed1,seed2,seed3);
			} 
		} 
	} 
	for( t= 1; t<=_nT; t++) 
		ga[t]= 4+runif(seed1,seed2,seed3);
	sum= 0;
	for( t= 1; t<=_nT; t++) 
		sum= sum+ga[t];
	for( t= 1; t<=_nT; t++) 
		ga[t]= ga[t]/sum;
}

// generate parameters model M_3
void generate_true_parameters(extended3 ro, extended3 ta, extended1 ga)
{

	unsigned short int o,a,t,f;
	long double     sum;
  
	for( o= 1; o<=_nO; o++)  
	{ 
		for( f= 1; f<=_nF; f++)  
		{ 
			for( t= 1; t<=_nT; t++)  
			{ 
				ro[o][f][t]= runif(seed1,seed2,seed3);
			} 
		} 
	} 
	for( a= 1; a<=_nA; a++)  
	{ 
		for( f= 1; f<=_nF; f++)  
		{ 
			for( t= 1; t<=_nT; t++)  
			{ 
				ta[a][f][t]= runif(seed1,seed2,seed3);
			} 
		} 
	} 
	for( t= 1; t<=_nT; t++) 
		ga[t]= 4+runif(seed1,seed2,seed3);
	sum= 0;
	for( t= 1; t<=_nT; t++) 
		sum= sum+ga[t];
	for( t= 1; t<=_nT; t++) 
		ga[t]= ga[t]/sum;
}

// function to compute conditional probability of observing an object-attribute association given object classification (model M_1)
// disjunctive mapping rule
void calculate_condprobX_DC(extended2 ta, extended2 condprobx)
{

	unsigned short int s,a,f;
	long double     prod;
  
	for( s= 1; s<=_nS; s++)  { 
	for( a= 1; a<=_nA; a++)  { 
	prod= 1;
	for( f= 1; f<=_nF; f++)  { 
	prod= prod*(1-patS[s][f]*ta[a][f]);
	} 
	condprobx[s][a]= 1-prod;
	} 
	} 
}


// function to compute conditional probability of observing an object-attribute association given object classification and given person classification (model M_2, M_3)
// disjunctive mapping rule
void calculate_condprobX_DC(extended3 ta, extended3 condprobx)
{

	unsigned short int s,a,f,t;
	long double     prod;
  
	for( s= 1; s<=_nS; s++)  { 
	for( a= 1; a<=_nA; a++)  { 
	for( t= 1; t<=_nT; t++)  { 
	prod= 1;
	for( f= 1; f<=_nF; f++)  { 
	prod= prod*(1-patS[s][f]*ta[a][f][t]);
	} 
	condprobx[s][a][t]= 1-prod;
	} 
	} 
	} 
}

// function to compute conditional probability of observing an object-attribute association given object classification (model M_1)
// additive mapping rule
void calculate_condprobX_ADD(extended2 ta, extended2 condprobx)
{

	unsigned short int s,a,f;
	long double  sum;
  
	for( s= 1; s<=_nS; s++)  { 
	for( a= 1; a<=_nA; a++)  { 
	sum=0;
	for( f= 1; f<=_nF; f++)  { 
	sum= sum + (patS[s][f]*ta[a][f]);
	} 
	condprobx[s][a]=sum/_nF;
	} 
	} 
}


// function to compute conditional probability of observing an object-attribute association given object classification and given person classification (model M_2, M_3)
// additive mapping rule
void calculate_condprobX_ADD(extended3 ta, extended3 condprobx)
{

	unsigned short int s,a,f,t;
	long double  sum;
  
	for( s= 1; s<=_nS; s++)  { 
	for( a= 1; a<=_nA; a++)  { 
	for( t= 1; t<=_nT; t++)  { 
	sum= 0;
	for( f= 1; f<=_nF; f++)  { 
	sum= sum + (patS[s][f]*ta[a][f][t]);
	} 
	condprobx[s][a][t]=sum/_nF;
	} 
	} 
	} 
}


//function to compute the marginal probability of an object classification  given person classification (model M_1, M3)
void calculate_margprobX(extended3 ro, extended3 margprobx)
{
	// XZ_Y
	unsigned short int  s,o,t,f;
	long double      prod;
  
	for( s= 1; s<=_nS; s++)  { 
	for( t= 1; t<=_nT; t++)  { 
	for( o= 1; o<=_nO; o++)  { 
	prod= 1;
	for( f= 1; f<=_nF; f++)  { 
		prod= prod*pow_1(ro[o][f][t],patS[s][f]);
		} 
	margprobx[s][o][t]= prod;
	} 
	} 
	} 
}

//function to compute the marginal probability of an object classification  (model M_2)
void calculate_margprobX(extended2 ro, extended2 margprobx)
{
	// X_YZ
	unsigned short int  s,o,f;
	long double      prod;
  
	for( s= 1; s<=_nS; s++)  { 
	for( o= 1; o<=_nO; o++)  { 
	prod= 1;
	for( f= 1; f<=_nF; f++)  { 
		prod= prod*pow_1(ro[o][f],patS[s][f]);
		} 
	margprobx[s][o]= prod;
	} 
	} 

}

// function to replace old parameter values by new parameter values (model M_1)
void replace(extended3 ro_o,extended3 ro_n,extended2 ta_o,extended2 ta_n,extended1 ga_o,extended1 ga_n, extended3 ro_update, extended2 ta_update, extended2 condprobx_o,extended2 condprobx_n,extended3 margprobx_o,extended3 margprobx_n)
{

	unsigned short int o,a,t,f,s;

	for( o= 1; o<=_nO; o++)  { 
	for( f= 1; f<=_nF; f++)  { 
	for( t= 1; t<=_nT; t++)  { 
	 if (ro_update[o][f][t]==1){ro_o[o][f][t]= ro_n[o][f][t];}
	} 
	} 
	} 



	for( a= 1; a<=_nA; a++)  { 
	for( f= 1; f<=_nF; f++)  { 
	if (ta_update[a][f]==1) {ta_o[a][f]= ta_n[a][f];}
	} 
	} 



	for( t= 1; t<=_nT; t++) ga_o[t]= ga_n[t];



	for( s= 1; s<=_nS; s++)  { 
	for( a= 1; a<=_nA; a++)  { 
	condprobx_o[s][a]= condprobx_n[s][a];
	} 
	} 



	for( s= 1; s<=_nS; s++)  { 
	for( o= 1; o<=_nO; o++)  { 
	for( t= 1; t<=_nT; t++)  { 
	margprobx_o[s][o][t]= margprobx_n[s][o][t];
	} 
	} 
	} 

	logold= lognew;
}

// function to replace old parameter values by new parameter values (model M_2)
void replace(extended2 ro_o,extended2 ro_n,extended3 ta_o,extended3 ta_n,extended1 ga_o,extended1 ga_n, extended2 ro_update, extended3 ta_update, extended3 condprobx_o,extended3 condprobx_n,extended2 margprobx_o,extended2 margprobx_n)
{

	unsigned short int o,a,t,f,s;


	for( o= 1; o<=_nO; o++)  { 
	for( f= 1; f<=_nF; f++)  { 
	if (ro_update[o][f]==1) {ro_o[o][f]= ro_n[o][f];}
	} 
	} 

	for( a= 1; a<=_nA; a++)  { 
	for( f= 1; f<=_nF; f++)  { 
	for( t= 1; t<=_nT; t++)  { 
	if (ta_update[a][f][t]==1) {ta_o[a][f][t]= ta_n[a][f][t];}
	} 
	} 
	} 

	for( t= 1; t<=_nT; t++) ga_o[t]= ga_n[t];

	for( s= 1; s<=_nS; s++)  { 
	for( a= 1; a<=_nA; a++)  { 
	for( t= 1; t<=_nT; t++)  { 
		condprobx_o[s][a][t]= condprobx_n[s][a][t];
	} 
	} 
	} 

	for( s= 1; s<=_nS; s++)  { 
	for( o= 1; o<=_nO; o++)  { 
	margprobx_o[s][o]= margprobx_n[s][o];
	} 
	} 

	logold= lognew;


}

// function to replace old parameter values by new parameter values (model M_3)
void replace(extended3 ro_o,extended3 ro_n,extended3 ta_o,extended3 ta_n,extended1 ga_o,extended1 ga_n, extended3 ro_update, extended3 ta_update, 
              extended3 condprobx_o,extended3 condprobx_n,extended3 margprobx_o,extended3 margprobx_n)
{

	unsigned short int o,a,t,f,s;


	for( o= 1; o<=_nO; o++)  { 
	for( f= 1; f<=_nF; f++)  { 
	for( t= 1; t<=_nT; t++)  { 
	 if (ro_update[o][f][t]==1) {ro_o[o][f][t]= ro_n[o][f][t];}
	} 
	} 
	} 

	for( a= 1; a<=_nA; a++)  { 
	for( f= 1; f<=_nF; f++)  { 
	for( t= 1; t<=_nT; t++)  { 
	if (ta_update[a][f][t]==1) {ta_o[a][f][t]= ta_n[a][f][t];}
	} 
	} 
	} 

	for( t= 1; t<=_nT; t++) ga_o[t]= ga_n[t];

	for( s= 1; s<=_nS; s++)  { 
	for( a= 1; a<=_nA; a++)  { 
	for( t= 1; t<=_nT; t++)  { 
		condprobx_o[s][a][t]= condprobx_n[s][a][t];
	} 
	} 
	} 

	for( s= 1; s<=_nS; s++)  { 
	for( o= 1; o<=_nO; o++)  { 
	for( t= 1; t<=_nT; t++)  { 
	margprobx_o[s][o][t]= margprobx_n[s][o][t];
	} 
	} 
	} 

	logold= lognew;


}

// function to compute posterior probabilities of data given person classification (model M_1)
void calculate_probmat(extended1 ga, extended2 condprobx, extended3 margprobx,extended2 probmat)
{

	unsigned short int r,t,o,s,a;
	long double   prodobject,sumpattern,prod,sum;
  
	for( t= 1; t<=_nT; t++)  { 
	for( r= 1; r<=_nR; r++)  { 
	prodobject= 1;
	for( o= 1; o<=_nO; o++)  { 
		sumpattern= 0;
		for( s= 1; s<=_nS; s++)  { 
		prod= 1;
		for( a= 1; a<=_nA; a++)  { 
		prod= prod*pow_1(condprobx[s][a],dataB[o][a][r]);
		} 
		prod= prod*margprobx[s][o][t];
		sumpattern= sumpattern+prod;
		} 
	prodobject= prodobject*sumpattern;
		} 
	probmat[t][r]= prodobject*ga[t];
	} 
	} 

	for( r= 1; r<=_nR; r++)  { 
	sum= 0;
	for( t= 1; t<=_nT; t++) sum= sum+probmat[t][r];
	for( t= 1; t<=_nT; t++) probmat[t][r]= probmat[t][r]/sum;
	} 
}

// function to compute posterior probabilities of data given person classification (model M_2)
void calculate_probmat(extended1 ga, extended3 condprobx, extended2 margprobx,extended2 probmat)
{

	unsigned short int r,t,o,s,a;
	long double     prodobject,sumpattern,prod,sum;
  
	for( t= 1; t<=_nT; t++)  { 
	for( r= 1; r<=_nR; r++)  { 
	prodobject= 1;
	for( o= 1; o<=_nO; o++)  { 
		sumpattern= 0;
		for( s= 1; s<=_nS; s++)  { 
		prod= 1;
		for( a= 1; a<=_nA; a++)  { 
		prod= prod*pow_1(condprobx[s][a][t],dataB[o][a][r]);
		} 
		prod= prod*margprobx[s][o];
		sumpattern= sumpattern+prod;
		} 
	prodobject= prodobject*sumpattern;
		} 
	probmat[t][r]= prodobject*ga[t];
	} 
	} 

	for( r= 1; r<=_nR; r++)  { 
	sum= 0;
	for( t= 1; t<=_nT; t++) sum= sum+probmat[t][r];
	for( t= 1; t<=_nT; t++) probmat[t][r]= probmat[t][r]/sum;
	} 

}

// function to compute posterior probabilities of data given person classification (model M_3)
void calculate_probmat(extended1 ga, extended3 condprobx, extended3 margprobx,extended2 probmat)
{

	unsigned short int r,t,o,s,a;
	long double     prodobject,sumpattern,prod,sum;
  
	for( t= 1; t<=_nT; t++)  { 
	for( r= 1; r<=_nR; r++)  { 
	prodobject= 1;
	for( o= 1; o<=_nO; o++)  { 
		sumpattern= 0;
		for( s= 1; s<=_nS; s++)  { 
		prod= 1;
		for( a= 1; a<=_nA; a++)  { 
		prod= prod*pow_1(condprobx[s][a][t],dataB[o][a][r]);
		} 
		prod= prod*margprobx[s][o][t];
		sumpattern= sumpattern+prod;
		} 
	prodobject= prodobject*sumpattern;
		} 
	probmat[t][r]= prodobject*ga[t];
	} 
	} 

	for( r= 1; r<=_nR; r++)  { 
	sum= 0;
	for( t= 1; t<=_nT; t++) sum= sum+probmat[t][r];
	for( t= 1; t<=_nT; t++) probmat[t][r]= probmat[t][r]/sum;
	} 


}

// function to update class-size parameters
void emgamma(extended1 ga,extended2 probmat)
{
	unsigned short int t,r;
	long double     sum;
  
	for( t= 1; t<=_nT; t++)  { 
	sum= 0;
	for( r= 1; r<=_nR; r++) sum= sum+probmat[t][r];
	ga[t]= (sum+2.0/(double(_nT)))/(   (double(_nR)) + 2.0 *(double(_nT))/(double(_nT))  );
	} 
}

// function to compute probability of object classification given data and classification of person (model M_1) 
void calculate_condprob_pattern(extended2 condprobx, extended3 margprobx,extended3 somega)
{
	unsigned short int r,o,s,t,a;
	long double      prod,sum;

	for( r= 1; r<=_nR; r++)  { 
	for( o= 1; o<=_nO; o++)  { 
	for( t= 1; t<=_nT; t++)  { 
	for( s= 1; s<=_nS; s++)  { 
		prod= 1;
		for( a= 1; a<=_nA; a++)  { 
		prod= prod*pow_1(condprobx[s][a],dataB[o][a][r]);
		} 
		prod= prod*margprobx[s][o][t];
		omega[r][o][s][t]= prod;
		} 
	} 
	} 
	} 

	for( r= 1; r<=_nR; r++)  { 
	for( o= 1; o<=_nO; o++)  { 
	for( t= 1; t<=_nT; t++)  { 
	sum= 0;
	for( s= 1; s<=_nS; s++)  { 
		sum= sum+omega[r][o][s][t];
		} 
	somega[r][o][t]= sum;
	} 
	} 
	} 

	for( r= 1; r<=_nR; r++)  { 
	for( o= 1; o<=_nO; o++)  { 
	for( t= 1; t<=_nT; t++)  { 
	for( s= 1; s<=_nS; s++)  { 
		omega[r][o][s][t]= omega[r][o][s][t]/somega[r][o][t];
		} 
	} 
	} 
	} 
}

// function to compute probability of object classification given data and classification of person (model M_2) 
void calculate_condprob_pattern(extended3 condprobx, extended2 margprobx,extended3 somega)
{
	unsigned short int r,o,s,t,a;
	long double      prod,sum;

	for( r= 1; r<=_nR; r++)  { 
	for( o= 1; o<=_nO; o++)  { 
	for( t= 1; t<=_nT; t++)  { 
	for( s= 1; s<=_nS; s++)  { 
		prod= 1;
		for( a= 1; a<=_nA; a++)  { 
		prod= prod*pow_1(condprobx[s][a][t],dataB[o][a][r]);
		} 
		prod= prod*margprobx[s][o];
		omega[r][o][s][t]= prod;
		} 
	} 
	} 
	} 

	for( r= 1; r<=_nR; r++)  { 
	for( o= 1; o<=_nO; o++)  { 
	for( t= 1; t<=_nT; t++)  { 
	sum= 0;
	for( s= 1; s<=_nS; s++)  { 
		sum= sum+omega[r][o][s][t];
		} 
	somega[r][o][t]= sum;
	} 
	} 
	} 

	for( r= 1; r<=_nR; r++)  { 
	for( o= 1; o<=_nO; o++)  { 
	for( t= 1; t<=_nT; t++)  { 
	for( s= 1; s<=_nS; s++)  { 
		omega[r][o][s][t]= omega[r][o][s][t]/somega[r][o][t];
		} 
	} 
	} 
	} 

}

// function to compute probability of object classification given data and classification of person (model M_3) 
void calculate_condprob_pattern(extended3 condprobx, extended3 margprobx,extended3 somega)
{
	unsigned short int r,o,s,t,a;
	long double      prod,sum;

	for( r= 1; r<=_nR; r++)  { 
	for( o= 1; o<=_nO; o++)  { 
	for( t= 1; t<=_nT; t++)  { 
	for( s= 1; s<=_nS; s++)  { 
		prod= 1;
		for( a= 1; a<=_nA; a++)  { 
		prod= prod*pow_1(condprobx[s][a][t],dataB[o][a][r]);
		} 
		prod= prod*margprobx[s][o][t];
		omega[r][o][s][t]= prod;
		} 
	} 
	} 
	} 

	for( r= 1; r<=_nR; r++)  { 
	for( o= 1; o<=_nO; o++)  { 
	for( t= 1; t<=_nT; t++)  { 
	sum= 0;
	for( s= 1; s<=_nS; s++)  { 
		sum= sum+omega[r][o][s][t];
		} 
	somega[r][o][t]= sum;
	} 
	} 
	} 

	for( r= 1; r<=_nR; r++)  { 
	for( o= 1; o<=_nO; o++)  { 
	for( t= 1; t<=_nT; t++)  { 
	for( s= 1; s<=_nS; s++)  { 
		omega[r][o][s][t]= omega[r][o][s][t]/somega[r][o][t];
		} 
	} 
	} 
	} 


}

// function to update object parameters (model M_1, M_3)
void update_emro(extended3 ro,extended2 probmat, extended3 ro_update)
{

	unsigned short int r,s,f,o,t;
	long double     teller,sum,noemer;
  
	for( o= 1; o<=_nO; o++)  { 
	for( f= 1; f<=_nF; f++)  { 
	for( t= 1; t<=_nT; t++)  { 
          
    if (ro_update[o][f][t]==1){
                               
	// calculate numerator;
	teller= 0;
	for( r= 1; r<=_nR; r++)  { 
	sum= 0;
	for( s= 1; s<=_nS; s++)  { 
	sum= sum+patS[s][f]*omega[r][o][s][t];
	} 
	teller= teller+sum*probmat[t][r];
	} 

	// calculate denominator
	noemer= 0;
	for( r= 1; r<=_nR; r++)  { 
	noemer= noemer+probmat[t][r];
	} 

	ro[o][f][t]= (teller+1.0/((long double)(_nO*_nT)))/(noemer+2.0/((long double)(_nO*_nT))); //prior latent gold
    }
    
	} 
	} 
	} 

}

// function to update object parameters (model M_2)
void update_emro(extended2 ro,extended2 probmat, extended2 ro_update)
{

	unsigned short int r,s,f,o,t;
	long double     teller,sum,noemer;
  
	for( o= 1; o<=_nO; o++)  { 
	for( f= 1; f<=_nF; f++)  { 
    
    if (ro_update[o][f]==1)
    {

	// bereken numerator;
	teller= 0;
	for( r= 1; r<=_nR; r++)  { 
	for( t= 1; t<=_nT; t++)  { 
	sum= 0;
	for( s= 1; s<=_nS; s++)  { 
	sum= sum+patS[s][f]*omega[r][o][s][t];
	} 
	teller= teller+sum*probmat[t][r];
	} 
	} 

	// bereken denominator
	noemer= _nR;


	ro[o][f] = (teller+1.0/((long double)_nO))/(noemer+2.0/((long double)_nO)); //prior latent gold
    }
    
	} 
	} 

}

// function to update attribute parameters (model M_1)
// disjunctive mapping rule
void update_emta_DC(extended2 condprobx, extended2 ta_o, extended2 ta_n, extended2 probmat, extended2 ta_update)
{

	unsigned short int  o,a,t,f,r,s;
	long double     teller,noemer,sum1,sum2,py;
  

	for( a= 1; a<=_nA; a++)  { 
	for( f= 1; f<=_nF; f++)  { 
    
    if (ta_update[a][f]==1)      
    {
	// calculate denominator;

	noemer= 0;
	for( t= 1; t<=_nT; t++)  { 
	for( r= 1; r<=_nR; r++)  { 
	sum1= 0;
	for( o= 1; o<=_nO; o++)  { 
	for( s= 1; s<=_nS; s++)  { 
		sum1= sum1+patS[s][f]*omega[r][o][s][t];
		} 
	} 
	noemer= noemer+sum1*probmat[t][r];
	} 
	} 

	// calculate numerator;

	teller= 0;
	for( t= 1; t<=_nT; t++)  { 
	for( r= 1; r<=_nR; r++)  { 
	sum2= 0;
	for( o= 1; o<=_nO; o++)  { 
	for( s= 1; s<=_nS; s++)  { 
		if (patS[s][f]==0)
			py= 0; 
		else  { 
		if (dataB[o][a][r]==0)    
			py= 0;
		else py= ta_o[a][f]/(condprobx[s][a]);
		} 
		sum2= sum2+patS[s][f]*omega[r][o][s][t]*py;
		} 
	} 
	teller= teller+sum2*probmat[t][r];
	} 
	} 

	ta_n[a][f]= (teller+1.0/((long double)(_nA)))/(noemer+2.0/((long double)(_nA))); //prior latent gold

	
    }
    } 
	} 
}

// function to update attribute parameteters (model M_2, M_3)
// disjunctive mapping rule
void update_emta_DC(extended3 condprobx, extended3 ta_o, extended3 ta_n, extended2 probmat, extended3 ta_update)
{

	unsigned short int  o,a,t,f,r,s;
	long double     teller,noemer,sum1,sum2,py;
  

	for( a= 1; a<=_nA; a++)  { 
	for( f= 1; f<=_nF; f++)  { 
	for( t= 1; t<=_nT; t++)  { 
    
    if (ta_update[a][f][t]==1)
    { 
	// bereken denominator;

	noemer= 0;

	for( r= 1; r<=_nR; r++)  { 
	sum1= 0;
	for( o= 1; o<=_nO; o++)  { 
	for( s= 1; s<=_nS; s++)  { 
		sum1= sum1+patS[s][f]*omega[r][o][s][t];
		} 
	} 
	noemer= noemer+sum1*probmat[t][r];
	} 


	// bereken numerator;

	teller= 0;

	for( r= 1; r<=_nR; r++)  { 
	sum2= 0;
	for( o= 1; o<=_nO; o++)  { 
	for( s= 1; s<=_nS; s++)  { 
		if (patS[s][f]==0)
			py= 0;
		else  { 
		if (dataB[o][a][r]==0) py= 0; 
		else py= ta_o[a][f][t]/(condprobx[s][a][t]);
		} 
		sum2= sum2+patS[s][f]*omega[r][o][s][t]*py;
		} 
	} 
	teller= teller+sum2*probmat[t][r];
	} 


	ta_n[a][f][t]= (teller+1.0/((long double)(_nA*_nT)))/(noemer+2.0/((long double)(_nA*_nT))); //prior latent gold

    }
	} 
	} 
	} 

}


// function to update attribute parameters (model M_1)
// additive mapping rule
void update_emta_ADD(extended2 condprobx, extended2 ta_o, extended2 ta_n, extended2 probmat, extended2 ta_update)
{

	unsigned short int  o,a,t,f,r,s;
	long double     teller,noemer,sum1,sum2,py;
  

	for( a= 1; a<=_nA; a++)  { 
	for( f= 1; f<=_nF; f++)  { 
    
    if (ta_update[a][f]==1)      
    {
	// calculate denominator;

	noemer= 0;
	for( t= 1; t<=_nT; t++)  { 
	for( r= 1; r<=_nR; r++)  { 
	sum1= 0;
	for( o= 1; o<=_nO; o++)  { 
	for( s= 1; s<=_nS; s++)  { 
		sum1= sum1+patS[s][f]*omega[r][o][s][t];
		} 
	} 
	noemer= noemer+sum1*probmat[t][r];
	} 
	} 

	// calculate numerator;

	teller= 0;
	for( t= 1; t<=_nT; t++)  { 
	for( r= 1; r<=_nR; r++)  { 
	sum2= 0;
	for( o= 1; o<=_nO; o++)  { 
	for( s= 1; s<=_nS; s++)  { 
		if (patS[s][f]==0)
			py= 0; 
		else  { 
		if (dataB[o][a][r]==0) 
			py= (ta_o[a][f])*(_nF-1- _nF*condprobx[s][a]+patS[s][f]*ta_o[a][f])/(_nF-ta_o[a][f]- _nF*condprobx[s][a]+patS[s][f]*ta_o[a][f]);
		else py= (ta_o[a][f]*(1+ _nF*condprobx[s][a]-patS[s][f]*ta_o[a][f]))/(ta_o[a][f]+ _nF*condprobx[s][a]-patS[s][f]*ta_o[a][f]);
		} 
		sum2= sum2+patS[s][f]*omega[r][o][s][t]*py;
		} 
	} 
	teller= teller+sum2*probmat[t][r];
	} 
	} 

	ta_n[a][f]= (teller+1.0/((long double)(_nA)))/(noemer+2.0/((long double)(_nA))); //prior latent gold

	
    }
    } 
	} 
}


// function to update attribute parameteters (model M_2, M_3)
// additive mapping rule
void update_emta_ADD(extended3 condprobx, extended3 ta_o, extended3 ta_n, extended2 probmat, extended3 ta_update)
{

	unsigned short int  o,a,t,f,r,s;
	long double     teller,noemer,sum1,sum2,py;
  

	for( a= 1; a<=_nA; a++)  { 
	for( f= 1; f<=_nF; f++)  { 
	for( t= 1; t<=_nT; t++)  { 
    
    if (ta_update[a][f][t]==1)
    { 
	// bereken denominator;

	noemer= 0;

	for( r= 1; r<=_nR; r++)  { 
	sum1= 0;
	for( o= 1; o<=_nO; o++)  { 
	for( s= 1; s<=_nS; s++)  { 
		sum1= sum1+patS[s][f]*omega[r][o][s][t];
		} 
	} 
	noemer= noemer+sum1*probmat[t][r];
	} 


	// bereken numerator;

	teller= 0;

	for( r= 1; r<=_nR; r++)  { 
	sum2= 0;
	for (o= 1; o<=_nO; o++)  { 
	for (s= 1; s<=_nS; s++)  { 
		if (patS[s][f]==0)
			py= 0;
		else  { 
		if (dataB[o][a][r]==0) 
		     py=(ta_o[a][f][t])*(_nF-1- _nF*condprobx[s][a][t]+patS[s][f]*ta_o[a][f][t])/(_nF-ta_o[a][f][t]-_nF*condprobx[s][a][t]+patS[s][f]*ta_o[a][f][t]); 
		else py=(ta_o[a][f][t]*(1+_nF*condprobx[s][a][t]-patS[s][f]*ta_o[a][f][t]))/(ta_o[a][f][t]+_nF*condprobx[s][a][t]-patS[s][f]*ta_o[a][f][t]);
		} 
		sum2= sum2+patS[s][f]*omega[r][o][s][t]*py;
		} 
	} 
	teller= teller+sum2*probmat[t][r];
	} 


	ta_n[a][f][t]= (teller+1.0/((long double)(_nA*_nT)))/(noemer+2.0/((long double)(_nA*_nT))); //prior latent gold

    }
	} 
	} 
	} 

}




// function to compute value of convergence criterion (model M_1)
void calculate_conv(extended3 ro_o, extended3 ro_n,extended2 ta_o,extended2 ta_n,extended1 ga_o,extended1 ga_n, extended3 ro_update, extended2 ta_update)
{

	unsigned short int o,a,t,f;
  
	conv= 0;
	for( o= 1; o<=_nO; o++)  { 
	for( f= 1; f<=_nF; f++)  { 
	for( t= 1; t<=_nT; t++)  { 
	if (ro_update[o][f][t]==1) {conv= conv+ abs(ro_n[o][f][t]-ro_o[o][f][t]);}

	} 
	} 
	} 
	for( a= 1; a<=_nA; a++)  { 
	for( f= 1; f<=_nF; f++)  { 
	if (ta_update[a][f]==1) {conv= conv+abs(ta_n[a][f]-ta_o[a][f]);}
	} 
	} 
	for( t= 1; t<=_nT; t++) conv= conv+abs(ga_o[t]-ga_n[t]);

	conv = conv/(((long double)_nO)*((long double)_nF)*((long double)_nT)+((long double)_nA)*((long double)_nF)+((long double)_nT));
	

}

// function to compute value of convergence criterion (model M_2)
void calculate_conv(extended2 ro_o,extended2 ro_n,extended3 ta_o,extended3 ta_n,extended1 ga_o,extended1 ga_n, extended2 ro_update, extended3 ta_update)
{
	unsigned short int o,a,t,f;
  
	conv= 0;
	for( o= 1; o<=_nO; o++)  { 
	for( f= 1; f<=_nF; f++)  { 
	if (ro_update[o][f]==1) {conv= conv+abs(ro_n[o][f]-ro_o[o][f]);}
	} 
	} 
	for( a= 1; a<=_nA; a++)  { 
	for( f= 1; f<=_nF; f++)  { 
	for( t= 1; t<=_nT; t++)  { 
	if	(ta_update[a][f][t]==1) {conv= conv+abs(ta_n[a][f][t]-ta_o[a][f][t]);}
	} 
	} 
	} 
	for( t= 1; t<=_nT; t++) conv= conv+abs(ga_o[t]-ga_n[t]);
	conv = conv/(((long double)_nO)*((long double)_nF)+((long double)_nA)*((long double)_nF)*((long double)_nT)+((long double)_nT));
	


}

// function to compute value of convergence criterion (model M_3)
void calculate_conv(extended3 ro_o,extended3 ro_n,extended3 ta_o,extended3 ta_n,extended1 ga_o,extended1 ga_n, extended3 ro_update, extended3 ta_update)
{

	unsigned short int o,a,t,f;
  
	conv= 0;
	for( o= 1; o<=_nO; o++)  { 
	for( f= 1; f<=_nF; f++)  { 
	for( t= 1; t<=_nT; t++)  { 
	if (ro_update[o][f][t]==1) {conv= conv+ abs(ro_n[o][f][t]-ro_o[o][f][t]);}

	} 
	} 
	} 
	for( a= 1; a<=_nA; a++)  { 
	for( f= 1; f<=_nF; f++)  { 
	for( t= 1; t<=_nT; t++)  { 
	 if (ta_update[a][f][t]==1)	{conv= conv+abs(ta_n[a][f][t]-ta_o[a][f][t]);}
	} 
	} 
	} 
	for( t= 1; t<=_nT; t++) conv= conv+abs(ga_o[t]-ga_n[t]);
	conv = conv/(((long double)_nO)*((long double)_nF)*((long double)_nT)+((long double)_nA)*((long double)_nF)*((long double)_nT)+((long double)_nT));
	

}

// create arrays to store object/attribute parameters (Model M_1)
void Create_ro_ta(extended3 ro, extended2 ta)
{

	for (int i = 0; i < _nO+1; i++)
	{
		ro[i] = new long double * [_nF+1];
		for (int ii = 0; ii < _nF+1; ii++)
		{
			ro[i][ii] = new long double [_nT+1];
		}
	}

 
	for (int i = 0; i < _nA+1; i++)
	{
		ta[i] = new long double [_nF+1];
	}


}

// create arrays to store object/attribute parameters (Model M_2)
void Create_ro_ta(extended2 ro, extended3 ta)
{


	for (int i = 0; i < _nO+1; i++)
	{
		ro[i] = new long double [_nF+1];
	}


	for (int i = 0; i < _nA+1; i++)
	{
		ta[i] = new long double * [_nF+1];
		for (int ii = 0; ii < _nF+1; ii++)
		{
			ta[i][ii] = new long double [_nT+1];
		}
	}


}

// create arrays to store object/attribute parameters (Model M_3)
void Create_ro_ta(extended3 ro, extended3 ta)
{

	for (int i = 0; i < _nO+1; i++)
	{
		ro[i] = new long double * [_nF+1];
		for (int ii = 0; ii < _nF+1; ii++)
		{
			ro[i][ii] = new long double [_nT+1];
		}
	}


	for (int i = 0; i < _nA+1; i++)
	{
		ta[i] = new long double * [_nF+1];
		for (int ii = 0; ii < _nF+1; ii++)
		{
			ta[i][ii] = new long double [_nT+1];
		}
	}

}


// create arrays to store marginal en conditional probabilities needed for computations (Model M_1)
void CreateVariables(extended3 somega,extended2 condprobx_o,extended2 condprobx_n,extended3 margprobx_o,extended3 margprobx_n,extended2 probmat, extended3  p_o_r_t,extended2 p_r_t,extended1 p_r)
{

	omega = new long double *** [_nR+1]; 
	for (int i = 0; i < _nR+1; i++)
	{
		omega[i] = new long double ** [_nO+1];
		for (int ii = 0; ii < _nO+1; ii++)
		{
			omega[i][ii] = new long double * [_nS+1];
			for (int iii = 0; iii < _nS+1; iii++)
			{
				omega[i][ii][iii] = new long double [_nT+1];
			}
		}
	}


	
	for (int i = 0; i < _nR+1; i++)
	{
		somega[i] = new long double * [_nO+1];
		for (int ii = 0; ii < _nO+1; ii++)
		{
			somega[i][ii] = new long double [_nT+1];
		}
	}


	
	for (int i = 0; i < _nS+1; i++)
	{
		condprobx_o[i] = new long double [_nA+1];
	}


	
	for (int i = 0; i < _nS+1; i++)
	{
		margprobx_o[i] = new long double * [_nO+1];
		for (int ii = 0; ii < _nO+1; ii++)
		{
			margprobx_o[i][ii] = new long double [_nT+1];
		}
	}
	

	
	for (int i = 0; i < _nS+1; i++)
	{
		condprobx_n[i] = new long double [_nA+1];
	}


	
	for (int i = 0; i < _nS+1; i++)
	{
		margprobx_n[i] = new long double * [_nO+1];
		for (int ii = 0; ii < _nO+1; ii++)
		{
			margprobx_n[i][ii] = new long double [_nT+1];
		}
	}


	
	for (int i = 0; i < _nT+1; i++)
	{
		probmat[i] = new long double [_nR+1];
	}


	
	for (int i = 0; i < _nO+1; i++)
	{
		p_o_r_t[i] = new long double * [_nR+1];
		for (int ii = 0; ii < _nR+1; ii++)
		{
			p_o_r_t[i][ii] = new long double [_nT+1];
		}
	}


	 
	for (int i = 0; i < _nR+1; i++)
	{
		p_r_t[i] = new long double [_nT+1];
	}



}

// create arrays to store marginal en conditional probabilities needed for computations (Model M_2)
void CreateVariables(extended3 somega, extended3 condprobx_o, extended3 condprobx_n, extended2 margprobx_o, extended2 margprobx_n, extended2 probmat, extended3 p_o_r_t, extended2 p_r_t, extended1 p_r)
{

	omega = new long double *** [_nR+1]; 
	for (int i = 0; i < _nR+1; i++)
	{
		omega[i] = new long double ** [_nO+1];
		for (int ii = 0; ii < _nO+1; ii++)
		{
			omega[i][ii] = new long double * [_nS+1];
			for (int iii = 0; iii < _nS+1; iii++)
			{
				omega[i][ii][iii] = new long double [_nT+1];
			}
		}
	}


	
	for (int i = 0; i < _nR+1; i++)
	{
		somega[i] = new long double * [_nO+1];
		for (int ii = 0; ii < _nO+1; ii++)
		{
			somega[i][ii] = new long double [_nT+1];
		}
	}


	
	for (int i = 0; i < _nS+1; i++)
	{
		condprobx_o[i] = new long double * [_nA+1];
		for (int ii = 0; ii < _nA+1; ii++)
		{
			condprobx_o[i][ii] = new long double [_nT+1];
		}
	}

	for (int i = 0; i < _nS+1; i++)
	{
		condprobx_n[i] = new long double * [_nA+1];
		for (int ii = 0; ii < _nA+1; ii++)
		{
			condprobx_n[i][ii] = new long double [_nT+1];
		}
	}


	for (int i = 0; i < _nS+1; i++)
	{
		margprobx_o[i] = new long double  [_nO+1];
	}
	for (int i = 0; i < _nS+1; i++)
	{
		margprobx_n[i] = new long double  [_nO+1];
	}
	


	
	for (int i = 0; i < _nT+1; i++)
	{
		probmat[i] = new long double [_nR+1];
	}

	
	for (int i = 0; i < _nO+1; i++)
	{
		p_o_r_t[i] = new long double * [_nR+1];
		for (int ii = 0; ii < _nR+1; ii++)
		{
			p_o_r_t[i][ii] = new long double [_nT+1];
		}
	}


	 
	for (int i = 0; i < _nR+1; i++)
	{
		p_r_t[i] = new long double [_nT+1];
	}


	

}

// create arrays to store marginal en conditional probabilities needed for computations (Model M_3)
void CreateVariables(extended3 somega, extended3 condprobx_o, extended3 condprobx_n, extended3 margprobx_o, extended3 margprobx_n, extended2 probmat, extended3 p_o_r_t, extended2 p_r_t, extended1 p_r)
{

	omega = new long double *** [_nR+1]; 
	for (int i = 0; i < _nR+1; i++)
	{
		omega[i] = new long double ** [_nO+1];
		for (int ii = 0; ii < _nO+1; ii++)
		{
			omega[i][ii] = new long double * [_nS+1];
			for (int iii = 0; iii < _nS+1; iii++)
			{
				omega[i][ii][iii] = new long double [_nT+1];
			}
		}
	}


	
	for (int i = 0; i < _nR+1; i++)
	{
		somega[i] = new long double * [_nO+1];
		for (int ii = 0; ii < _nO+1; ii++)
		{
			somega[i][ii] = new long double [_nT+1];
		}
	}


	
	for (int i = 0; i < _nS+1; i++)
	{
		condprobx_o[i] = new long double * [_nA+1];
		for (int ii = 0; ii < _nA+1; ii++)
		{
			condprobx_o[i][ii] = new long double [_nT+1];
		}
	}

	for (int i = 0; i < _nS+1; i++)
	{
		condprobx_n[i] = new long double * [_nA+1];
		for (int ii = 0; ii < _nA+1; ii++)
		{
			condprobx_n[i][ii] = new long double [_nT+1];
		}
	}

	for (int i = 0; i < _nS+1; i++)
	{
		margprobx_o[i] = new long double * [_nO+1];
		for (int ii = 0; ii < _nO+1; ii++)
		{
			margprobx_o[i][ii] = new long double [_nT+1];
		}
	}

	for (int i = 0; i < _nS+1; i++)
	{
		margprobx_n[i] = new long double * [_nO+1];
		for (int ii = 0; ii < _nO+1; ii++)
		{
			margprobx_n[i][ii] = new long double [_nT+1];
		}
	}


	
	for (int i = 0; i < _nT+1; i++)
	{
		probmat[i] = new long double [_nR+1];
	}


	
	for (int i = 0; i < _nO+1; i++)
	{
		p_o_r_t[i] = new long double * [_nR+1];
		for (int ii = 0; ii < _nR+1; ii++)
		{
			p_o_r_t[i][ii] = new long double [_nT+1];
		}
	}


	 
	for (int i = 0; i < _nR+1; i++)
	{
		p_r_t[i] = new long double [_nT+1];
	}


	

}

// function to delete variables 
void R_destructor(binary3 R_data,binary2 R_patS,extended3 ro,extended2 ta,extended1 ga,extended3 ro2,extended2 ta2,extended1 ga2)
{
	for (int i = 0; i < _nO+1; i++) 
	{
		for (int ii = 0; ii < _nA+1; ii++)
		{
			delete R_data[i][ii];
		}
		delete R_data[i];
	}
	delete[] R_data;
	for (int i = 0; i < _nS+1; i++)
	{
		delete R_patS[i];
	}
	delete[] R_patS;
	for (int i = 0; i < _nO+1; i++)
	{
		for (int ii = 0; ii < _nF+1; ii++)
		{
			delete ro[i][ii];
		}
		delete ro[i];
	}
	delete[] ro;
	for (int i = 0; i < _nA+1; i++)
	{
		delete ta[i];
	}
	delete[] ta;
	delete[] ga;
	for (int i = 0; i < _nO+1; i++)
	{
		for (int ii = 0; ii < _nF+1; ii++)
		{
			delete ro2[i][ii];
		}
		delete ro2[i];
	}
	delete[] ro2;
	for (int i = 0; i < _nA+1; i++)
	{
		delete ta2[i];
	}
	delete[] ta2;
	delete[] ga2;

}

// function to delete arrays with probabilities needed for computations (model M_1)
void C_destructor(extended3 somega,extended2 condprobx_o,extended2 condprobx_n,extended3 margprobx_o,extended3 margprobx_n,extended2 probmat, extended3  p_o_r_t,extended2 p_r_t,extended1 p_r)
{

	for (int i = 0; i < _nR+1; i++)
	{
		for (int ii = 0; ii < _nO+1; ii++)
		{
			for (int iii = 0; iii < _nS+1; iii++)
			{
				delete omega[i][ii][iii];
			}
			delete omega[i][ii];
		}
		delete omega[i];
	}
	delete[] omega; 

	for (int i = 0; i < _nR+1; i++)
	{
		for (int ii = 0; ii < _nO+1; ii++)
		{
			delete somega[i][ii];
		}
		delete somega[i];
	}
	delete[] somega; 

	for (int i = 0; i < _nS+1; i++)
	{
		delete condprobx_o[i];
	}
	delete[] condprobx_o; 

	for (int i = 0; i < _nS+1; i++)
	{
		for (int ii = 0; ii < _nO+1; ii++)
		{
			delete margprobx_o[i][ii];
		}
		delete margprobx_o[i];
	}
	delete[] margprobx_o;
	
	for (int i = 0; i < _nS+1; i++)
	{
		delete condprobx_n[i];
	}
	delete[] condprobx_n; 

	for (int i = 0; i < _nS+1; i++)
	{
		for (int ii = 0; ii < _nO+1; ii++)
		{
			delete margprobx_n[i][ii];
		}
		delete margprobx_n[i];
	}
	delete[] margprobx_n;

	for (int i = 0; i < _nO+1; i++)
	{
		for (int ii = 0; ii < _nR+1; ii++)
		{
			delete p_o_r_t[i][ii];
		}
		delete p_o_r_t[i];
	}
	delete[] p_o_r_t;

	for (int i = 0; i < _nR+1; i++)
	{
		delete p_r_t[i];
	}
	delete[] p_r_t;

	delete[] p_r;

	

}

// function to delete arrays with probabilities needed for computations (model M_2)
void C_destructor(extended3 somega,extended3 condprobx_o,extended3 condprobx_n,extended2 margprobx_o,extended2 margprobx_n,extended2 probmat, extended3  p_o_r_t,extended2 p_r_t,extended1 p_r)
{

	for (int i = 0; i < _nR+1; i++)
	{
		for (int ii = 0; ii < _nO+1; ii++)
		{
			for (int iii = 0; iii < _nS+1; iii++)
			{
				delete omega[i][ii][iii];
			}
			delete omega[i][ii];
		}
		delete omega[i];
	}
	delete[] omega; 

	for (int i = 0; i < _nR+1; i++)
	{
		for (int ii = 0; ii < _nO+1; ii++)
		{
			delete somega[i][ii];
		}
		delete somega[i];
	}
	delete[] somega; 

	for (int i = 0; i < _nS+1; i++)
	{
		delete margprobx_o[i];
		delete margprobx_n[i];
	}
	delete[] margprobx_o;
	delete[] margprobx_n;

	for (int i = 0; i < _nS+1; i++)
	{
		for (int ii = 0; ii < _nA+1; ii++)
		{
			delete condprobx_o[i][ii];
			delete condprobx_n[i][ii];
		}
		delete condprobx_o[i];
		delete condprobx_n[i];
	}
	delete[] condprobx_o;
	delete[] condprobx_n;
	
	for (int i = 0; i < _nO+1; i++)
	{
		for (int ii = 0; ii < _nR+1; ii++)
		{
			delete p_o_r_t[i][ii];
		}
		delete p_o_r_t[i];
	}
	delete[] p_o_r_t;

	for (int i = 0; i < _nR+1; i++)
	{
		delete p_r_t[i];
	}
	delete[] p_r_t;

	delete[] p_r;

	

}

// function to delete arrays with probabilities needed for computations (model M_3)
void C_destructor(extended3 somega,extended3 condprobx_o,extended3 condprobx_n,extended3 margprobx_o,extended3 margprobx_n,extended2 probmat, extended3  p_o_r_t,extended2 p_r_t,extended1 p_r)
{

	for (int i = 0; i < _nR+1; i++)
	{
		for (int ii = 0; ii < _nO+1; ii++)
		{
			for (int iii = 0; iii < _nS+1; iii++)
			{
				delete omega[i][ii][iii];
			}
			delete omega[i][ii];
		}
		delete omega[i];
	}
	delete[] omega; 

	for (int i = 0; i < _nR+1; i++)
	{
		for (int ii = 0; ii < _nO+1; ii++)
		{
			delete somega[i][ii];
		}
		delete somega[i];
	}
	delete[] somega; 


	for (int i = 0; i < _nS+1; i++)
	{
		for (int ii = 0; ii < _nO+1; ii++)
		{
			delete margprobx_o[i][ii];
			delete margprobx_n[i][ii];
		}
		delete margprobx_o[i];
		delete margprobx_n[i];

	}
	delete[] margprobx_o;
	delete[] margprobx_n;

	for (int i = 0; i < _nS+1; i++)
	{
		for (int ii = 0; ii < _nA+1; ii++)
		{
			delete condprobx_o[i][ii];
			delete condprobx_n[i][ii];
		}
		delete condprobx_o[i];
		delete condprobx_n[i];
	}
	delete[] condprobx_o;
	delete[] condprobx_n;
	
	for (int i = 0; i < _nO+1; i++)
	{
		for (int ii = 0; ii < _nR+1; ii++)
		{
			delete p_o_r_t[i][ii];
		}
		delete p_o_r_t[i];
	}
	delete[] p_o_r_t;

	for (int i = 0; i < _nR+1; i++)
	{
		delete p_r_t[i];
	}
	delete[] p_r_t;

	delete[] p_r;

	

}

// function to draw from a uniform distribution
double runif(int xp1, int xp2, int xp3)
{   
	double intpart = 0.0;
	seed1= 171*(xp1  %  177) - 2*(xp1  /  177);
	seed2= 172*(xp2  %  176) - 35*(xp2  /  176);
	seed3= 170*(xp3  %  178) - 63*(xp3  /  178);
	if(seed1 < 0   ) seed1= seed1 + 30269;
	if(seed2 < 0   ) seed2= seed2 + 30307;
	if(seed3 < 0   ) seed3= seed3 + 30323;
	return modf(((double)seed1/30269.0 + (double)seed2/30307.0 + (double)seed3/30323.0), &intpart);
} 

// function to draw from a standard normal distribution
double rstnorm(void)
{
	double tp1,tp2,tp3,tp4,tp5;
	tp1= runif(seed1,seed2,seed3);
	tp2= 2.0*M_PI*tp1;
	tp3= runif(seed1,seed2,seed3);
	tp4= sqrt(2.0*(-log(tp3)));    
	tp5= tp4*cos(tp2);
	return tp5;
}

// function to draw from a gamma distribution
double rgamma_best(double a)
{
	double b,c,U,V,W,Y,X,Z;
	binary accept;

	accept= 0;
	b= a-1;
	c= 3.0*a-0.75;
	do
	{
		U= runif(seed1,seed2,seed3);
		V= runif(seed1,seed2,seed3);
		W= U*(1-U);
		Y= sqrt(c/W)*(U-0.5);
		X= b+Y;
		if((X>=0)   )
		{ 
			Z= 64*W*W*W*V*V;
			if((Z<=1-(2.0*sqrt(Y)/X))   ) 
				accept= 1 ;
			else
			{ 
				if((log(Z)<=2.0*(b*log(X/b)-Y))   ) 
					accept= 1;
			}  
		} 
	}
	while (accept==1);
	return X;

} 

// function to draw from a beta distribution
double rbeta(double x1, double x2)
{
	double u1,u2;
	long double v1,v2,sum;
	double gamma_x1,gamma_x2;

  
	if(((x1<=1) || (x2<=1))   )
	{ 
		do
		{
			u1= runif(seed1,seed2,seed3);
			u2= runif(seed1,seed2,seed3);

			v1= (1.0/x1)*log(u1);
			if((v1>1E4)   ) 
				v1= 1E4;
			if((v1<-1E4)   ) 
				v1= -1E4;
			v1= exp(v1);

			v2= (1.0/x2)*log(u2);
			if((v2>1E4)   ) 
				v2= 1E4;
			if((v2<-1E4)   ) 
				v2= -1E4;
			v2= exp(v2);
			sum= v1+v2;
		}
		while( (sum<=1));
		return v1/sum;
	}  
	else
	{ 
		gamma_x1= rgamma_best(x1);
		gamma_x2= rgamma_best(x2);
		return gamma_x1/(gamma_x1+gamma_x2);
	} 
} 

// compute the gradient of all class size parameters using function element_gradient_ga
void calculate_gradient_ga(extended1 ga,extended1 gradient_ga,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r)
{
	unsigned short int u;

	for( u= 1; u<=_nT; u++)  { 
    gradient_ga[u]= element_gradient_ga(u,ga, p_o_r_t , p_r_t, p_r);
  } 
}

// compute the gradient of a specific class size parameter
long double element_gradient_ga(unsigned short int u, extended1 ga,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r)
{
	unsigned short int r;
	long double sum;

	sum= 0;
	for( r= 1; r<=_nR; r++)  { 
	sum = sum +(1.0/p_r[r])*p_r_t[r][u];
	} 
	return sum+(2.0/(double(_nT)))*(1.0/ga[u])-(_nR+2.0);

}

// compute the gradient of all attribute parameters for model M_1 using function element_gradient_ta
// disjunctive mapping rule
void calculate_gradient_ta_DC(extended2 ta, extended1 ga, extended2 condprobx, extended3 margprobx,extended2 gradient_ta,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r, extended2 ta_update)
{

	unsigned short int u,v;

	for( u= 1; u<=_nA; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
	if (ta_update[u][v]==1) {gradient_ta[u][v] = element_gradient_ta_DC(u,v,ta,ga,condprobx,margprobx, p_o_r_t , p_r_t, p_r);}
	} 
	} 
 } 

// compute the gradient of all attribute parameters for model M_2 using function element_gradient_ta
// disjunctive mapping rule
void calculate_gradient_ta_DC(extended3 ta, extended1 ga, extended3 condprobx, extended2 margprobx,extended3 gradient_ta,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r, extended3 ta_update)
{

	unsigned short int u,v,w;

	for( u= 1; u<=_nA; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
	for( w= 1; w<=_nT; w++)  { 

	if (ta_update[u][v][w]==1) {gradient_ta[u][v][w]= element_gradient_ta_DC(u,v,w,ta,ga,condprobx,margprobx,p_o_r_t,p_r_t, p_r );}

	} 
	} 
	} 

 } 

// compute the gradient of all attribute parameters for model M_3 using function element_gradient_ta
// disjunctive mapping rule
void calculate_gradient_ta_DC(extended3 ta, extended1 ga, extended3 condprobx, extended3 margprobx,extended3 gradient_ta,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r, extended3 ta_update)
{
	
	unsigned short int u,v,w;

	for( u= 1; u<=_nA; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
	for( w= 1; w<=_nT; w++)  { 

	if (ta_update[u][v][w]==1) {gradient_ta[u][v][w]= element_gradient_ta_DC(u,v,w,ta,ga,condprobx,margprobx,p_o_r_t,p_r_t, p_r );}

	} 
	} 
	} 

 } 


// compute the gradient of all attribute parameters for model M_1 using function element_gradient_ta
// additive mapping rule
void calculate_gradient_ta_ADD(extended2 ta, extended1 ga, extended2 condprobx, extended3 margprobx,extended2 gradient_ta,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r, extended2 ta_update)
{

	unsigned short int u,v;

	for( u= 1; u<=_nA; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
	if (ta_update[u][v]==1) {gradient_ta[u][v] = element_gradient_ta_ADD(u,v,ta,ga,condprobx,margprobx, p_o_r_t , p_r_t, p_r);}
	} 
	} 
 } 

// compute the gradient of all attribute parameters for model M_2 using function element_gradient_ta
// additive mapping rule
void calculate_gradient_ta_ADD(extended3 ta, extended1 ga, extended3 condprobx, extended2 margprobx,extended3 gradient_ta,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r, extended3 ta_update)
{

	unsigned short int u,v,w;

	for( u= 1; u<=_nA; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
	for( w= 1; w<=_nT; w++)  { 

	if (ta_update[u][v][w]==1) {gradient_ta[u][v][w]= element_gradient_ta_ADD(u,v,w,ta,ga,condprobx,margprobx,p_o_r_t,p_r_t, p_r );}

	} 
	} 
	} 

 } 


// compute the gradient of all attribute parameters for model M_3 using function element_gradient_ta
// additive mapping rule
void calculate_gradient_ta_ADD(extended3 ta, extended1 ga, extended3 condprobx, extended3 margprobx,extended3 gradient_ta,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r, extended3 ta_update)
{
	
	unsigned short int u,v,w;

	for( u= 1; u<=_nA; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
	for( w= 1; w<=_nT; w++)  { 

	if (ta_update[u][v][w]==1) {gradient_ta[u][v][w]= element_gradient_ta_ADD(u,v,w,ta,ga,condprobx,margprobx,p_o_r_t,p_r_t, p_r );}

	} 
	} 
	} 

 } 


// function to compute the gradient of one attribute parameter for model M_1 (see appendix paper)
// disjunctive mapping rule
long double element_gradient_ta_DC(unsigned short int u, unsigned short int v, extended2 ta, extended1 ga, extended2 condprobx, extended3 margprobx,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r)
{

	unsigned short int r,t,o,a,s;
	long double     sum_s,sum_o,sum_t,sum_r,prod;
  

	sum_r= 0;
	for( r= 1; r<=_nR; r++)  { 
	sum_t= 0;
	for( t= 1; t<=_nT; t++)  { 
	sum_o= 0;
	for( o= 1; o<=_nO; o++)  { 

		sum_s= 0;
		for( s= 1; s<=_nS; s++)  { 
		prod= 1;
			for( a= 1; a<=_nA; a++)  { 
			prod= prod*pow_1(condprobx[s][a],dataB[o][a][r]);
			} 
			prod= prod*margprobx[s][o][t];
			prod= prod/pow_3(s,u,condprobx,dataB[o][u][r]);
			prod= prod*(2.0*dataB[o][u][r]-1.0)*patS[s][v]*(1.0-condprobx[s][u])/(1.0-patS[s][v]*ta[u][v]);
		sum_s= sum_s+prod;
		} 
		sum_o= sum_o+p_r_t[r][t]*(1.0/p_o_r_t[o][r][t])*sum_s;
		} 
	sum_t= sum_t+ga[t]*sum_o;
	} 
	sum_r= sum_r+(1/p_r[r])*sum_t;
	} 

   return sum_r+(1.0/(double(_nA)))*(1.0/ta[u][v])-(1.0/(double(_nA)))*(1.0/(1.0-ta[u][v]));

 } 

// function to compute the gradient of one attribute parameter for model M_2 (see appendix paper)
// disjunctive mapping rule
long double element_gradient_ta_DC(unsigned short int u, unsigned short int v, unsigned short int w, extended3 ta, extended1 ga, extended3 condprobx, extended2 margprobx,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r)
{

	unsigned short int r,o,a,s;
	long double     sum_s,sum_o,sum_r,prod;
  

	sum_r= 0;
	for( r= 1; r<=_nR; r++)  { 

	sum_o= 0;
	for( o= 1; o<=_nO; o++)  { 

		sum_s= 0;
		for( s= 1; s<=_nS; s++)  { 
		prod= 1;
			for( a= 1; a<=_nA; a++)  { 
			prod= prod*pow_1(condprobx[s][a][w],dataB[o][a][r]);
			} 
			prod= prod*margprobx[s][o];
			prod= prod/pow_3(s,u,w,condprobx,dataB[o][u][r]);
			prod= prod*(2.0*dataB[o][u][r]-1.0)*patS[s][v]*(1.0-condprobx[s][u][w])/(1.0-patS[s][v]*ta[u][v][w]);
		sum_s= sum_s+prod;
		} 
		sum_o= sum_o+p_r_t[r][w]*(1.0/p_o_r_t[o][r][w])*sum_s;
		} 

	sum_r= sum_r+ga[w]*(1.0/p_r[r])*sum_o;
	} 

	return (sum_r+(1.0/(double(_nA*_nT)))*(1.0/ta[u][v][w])-(1.0/(double(_nA*_nT)))*(1.0/(1.0-ta[u][v][w])));


 } 

// function to compute the gradient of one attribute parameter for model M_3 (see appendix paper)
// disjunctive mapping rule
long double element_gradient_ta_DC(unsigned short int u, unsigned short int v, unsigned short int w, extended3 ta, extended1 ga, extended3 condprobx, extended3 margprobx,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r)
{

	unsigned short int r,o,a,s;
	long double     sum_s,sum_o,sum_r,prod;
  
	sum_r= 0;
		for( r= 1; r<=_nR; r++)  { 

		sum_o= 0;
		for( o= 1; o<=_nO; o++)  { 

			sum_s= 0;
			for( s= 1; s<=_nS; s++)  { 
			prod= 1;
				for( a= 1; a<=_nA; a++)  { 
				prod= prod*pow_1(condprobx[s][a][w],dataB[o][a][r]);
				} 
				prod= prod*margprobx[s][o][w];
				prod= prod/pow_3(s,u,w,condprobx,dataB[o][u][r]);
				prod= prod*(2.0*dataB[o][u][r]-1.0)*patS[s][v]*(1.0-condprobx[s][u][w])/(1.0-patS[s][v]*ta[u][v][w]);
			sum_s= sum_s+prod;
			} 
			sum_o= sum_o+p_r_t[r][w]*(1.0/p_o_r_t[o][r][w])*sum_s;
			} 

		sum_r= sum_r+ga[w]*(1/p_r[r])*sum_o;
		} 

	return (sum_r+(1.0/(double(_nA*_nT)))*(1.0/ta[u][v][w])-(1.0/(double(_nA*_nT)))*(1.0/(1.0-ta[u][v][w])));


 } 


// function to compute the gradient of one attribute parameter for model M_1 (see appendix paper)
// additive mapping rule
long double element_gradient_ta_ADD(unsigned short int u, unsigned short int v, extended2 ta, extended1 ga, extended2 condprobx, extended3 margprobx,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r)
{

	unsigned short int r,t,o,a,s;
	long double     sum_s,sum_o,sum_t,sum_r,prod;
  

	sum_r= 0;
	for( r= 1; r<=_nR; r++)  { 
	sum_t= 0;
	for( t= 1; t<=_nT; t++)  { 
	sum_o= 0;
	for( o= 1; o<=_nO; o++)  { 

		sum_s= 0;
		for( s= 1; s<=_nS; s++)  { 
		prod= 1;
			for( a= 1; a<=_nA; a++)  { 
			prod= prod*pow_1(condprobx[s][a],dataB[o][a][r]);
			} 
			prod= prod*margprobx[s][o][t];
			prod= prod/pow_3(s,u,condprobx,dataB[o][u][r]);
			prod= prod*(2.0*dataB[o][u][r]-1.0)*patS[s][v]/(_nF);
		sum_s= sum_s+prod;
		} 
		sum_o= sum_o+p_r_t[r][t]*(1.0/p_o_r_t[o][r][t])*sum_s;
		} 
	sum_t= sum_t+ga[t]*sum_o;
	} 
	sum_r= sum_r+(1/p_r[r])*sum_t;
	} 

   return sum_r+(1.0/(double(_nA)))*(1.0/ta[u][v])-(1.0/(double(_nA)))*(1.0/(1.0-ta[u][v]));

 } 



// function to compute the gradient of one attribute parameter for model M_2 (see appendix paper)
// additive mapping rule
long double element_gradient_ta_ADD(unsigned short int u, unsigned short int v, unsigned short int w, extended3 ta, extended1 ga, extended3 condprobx, extended2 margprobx,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r)
{

	unsigned short int r,o,a,s;
	long double     sum_s,sum_o,sum_r,prod;
  
	sum_r= 0;
	for( r= 1; r<=_nR; r++)  { 

	sum_o= 0;
	for( o= 1; o<=_nO; o++)  { 

		sum_s= 0;
		for( s= 1; s<=_nS; s++)  { 
		prod= 1;
			for( a= 1; a<=_nA; a++)  { 
			prod= prod*pow_1(condprobx[s][a][w],dataB[o][a][r]);
			} 
			prod= prod*margprobx[s][o];
			prod= prod/pow_3(s,u,w,condprobx,dataB[o][u][r]);
			prod= prod*(2.0*dataB[o][u][r]-1.0)*patS[s][v]/(_nF);
		sum_s= sum_s+prod;
		} 
		sum_o= sum_o+p_r_t[r][w]*(1.0/p_o_r_t[o][r][w])*sum_s;
		} 

	sum_r= sum_r+ga[w]*(1.0/p_r[r])*sum_o;
	} 

	return (sum_r+(1.0/(double(_nA*_nT)))*(1.0/ta[u][v][w])-(1.0/(double(_nA*_nT)))*(1.0/(1.0-ta[u][v][w])));


 } 

// function to compute the gradient of one attribute parameter for model M_3 (see appendix paper)
// additive mapping rule
long double element_gradient_ta_ADD(unsigned short int u, unsigned short int v, unsigned short int w, extended3 ta, extended1 ga, extended3 condprobx, extended3 margprobx,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r)
{

	unsigned short int r,o,a,s;
	long double     sum_s,sum_o,sum_r,prod;
  
	sum_r= 0;
		for( r= 1; r<=_nR; r++)  { 

		sum_o= 0;
		for( o= 1; o<=_nO; o++)  { 

			sum_s= 0;
			for( s= 1; s<=_nS; s++)  { 
			prod= 1;
				for( a= 1; a<=_nA; a++)  { 
				prod= prod*pow_1(condprobx[s][a][w],dataB[o][a][r]);
				} 
				prod= prod*margprobx[s][o][w];
				prod= prod/pow_3(s,u,w,condprobx,dataB[o][u][r]);
				prod= prod*(2.0*dataB[o][u][r]-1.0)*patS[s][v]/(_nF);
			sum_s= sum_s+prod;
			} 
			sum_o= sum_o+p_r_t[r][w]*(1.0/p_o_r_t[o][r][w])*sum_s;
			} 

		sum_r= sum_r+ga[w]*(1/p_r[r])*sum_o;
		} 

	return (sum_r+(1.0/(double(_nA*_nT)))*(1.0/ta[u][v][w])-(1.0/(double(_nA*_nT)))*(1.0/(1.0-ta[u][v][w])));


 } 


// compute the gradient of all object parameters for model M_1 using function element_gradient_ro
void calculate_gradient_ro(extended3 ro, extended1 ga, extended3 margprobx, extended2 condprobx,extended3 gradient_ro,
                            extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r, extended3 ro_update)
{

	unsigned short int u,v,w;
  
	for( u= 1; u<=_nO; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
	for( w= 1; w<=_nT; w++)  { 

	if (ro_update[u][v][w]==1) {gradient_ro[u][v][w]= element_gradient_ro(u,v,w,ro,ga,margprobx,condprobx, p_o_r_t , p_r_t, p_r);}

	} 
	} 
	} 
	
}

// compute the gradient of all object parameters for model M_2 using function element_gradient_ro
void calculate_gradient_ro(extended2 ro, extended1 ga, extended2 margprobx, extended3 condprobx,extended2 gradient_ro,
                           extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r, extended2 ro_update)
{

	unsigned short int u,v;
  
	for( u= 1; u<=_nO; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
	if (ro_update[u][v]==1) {gradient_ro[u][v] = element_gradient_ro(u,v,ro,ga,margprobx,condprobx,p_o_r_t,p_r_t,p_r);}
	} 
	} 

	
}

// compute the gradient of all object parameters for model M_3 using function element_gradient_ro
void calculate_gradient_ro(extended3 ro, extended1 ga, extended3 margprobx, extended3 condprobx,extended3 gradient_ro,
                            extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r, extended3 ro_update)
{

	unsigned short int u,v,w;
  
	for( u= 1; u<=_nO; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
	for( w= 1; w<=_nT; w++)  { 

	if (ro_update[u][v][w]==1) {gradient_ro[u][v][w]= element_gradient_ro(u,v,w,ro,ga,margprobx,condprobx, p_o_r_t , p_r_t, p_r);}

	} 
	} 
	} 

	
}

// function to compute the gradient of one object parameter for model M_1 (see appendix paper)
long double element_gradient_ro(unsigned short int u, unsigned short int v, unsigned short int w,  extended3 ro, extended1 ga,
                              extended3 margprobx, extended2 condprobx,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r)
{

	unsigned short int s,a,r;
	long double sum_r,prod,sum_s;
    sum_r= 0;
	for( r= 1; r<=_nR; r++)  { 
		sum_s= 0;
		for( s= 1; s<=_nS; s++)  { 
		prod= 1;
		for( a= 1; a<=_nA; a++)  { 
			prod= prod*pow_1(condprobx[s][a],dataB[u][a][r]);
			} 
		prod= prod*margprobx[s][u][w];
		prod= prod/(ro[u][v][w]*patS[s][v]+(1.0-ro[u][v][w])*(1.0 -patS[s][v]));
		prod= prod*(2.0*patS[s][v] - 1.0);
		sum_s= sum_s+prod;
		} 
	sum_r= sum_r+sum_s*ga[w]*(1.0/p_r[r])*p_r_t[r][w]*(1.0/p_o_r_t[u][r][w]);
		} 

    return sum_r+(1.0/(double(_nO*_nT)))*(1.0/ro[u][v][w])-(1.0/(double(_nO*_nT)))*(1.0/(1.0-ro[u][v][w]));
} 

// function to compute the gradient of one object parameter for model M_2 (see appendix paper)
long double element_gradient_ro(unsigned short int u, unsigned short int v,  extended2 ro, extended1 ga,
                              extended2 margprobx, extended3 condprobx,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r)
{

	unsigned short int s,a,r,t;
	long double sum_r,prod,sum_s,sum_t;

    sum_r= 0;
    for( r= 1; r<=_nR; r++)  { 
     sum_t= 0;
     for( t= 1; t<=_nT; t++)  { 
        sum_s= 0;
         for( s= 1; s<=_nS; s++)  { 
          prod= 1;
          for( a= 1; a<=_nA; a++)  { 
           prod= prod*pow_1(condprobx[s][a][t],dataB[u][a][r]);
           } 
          prod= prod*margprobx[s][u];
          prod= prod/(ro[u][v]*patS[s][v]+(1.0-ro[u][v])*(1.0-patS[s][v]));
          prod= prod*(2.0*patS[s][v]-1.0);
          sum_s= sum_s+prod;
          } 
       sum_t= sum_t+sum_s*ga[t]*p_r_t[r][t]*(1.0/p_o_r_t[u][r][t]);
       } 
    sum_r= sum_r+sum_t*(1.0/p_r[r]);
     } 

	return (sum_r+(1.0/(double(_nO)))*(1.0/ro[u][v])-(1.0/(double(_nO)))*(1.0/(1.0-ro[u][v])));


} 

// function to compute the gradient of one object parameter for model M_3 (see appendix paper)
long double element_gradient_ro(unsigned short int u, unsigned short int v, unsigned short int w, extended3 ro, extended1 ga,
                              extended3 margprobx, extended3 condprobx,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r)
{

	unsigned short int s,a,r;
	long double sum_r,prod,sum_s;

	sum_r= 0;
	for( r= 1; r<=_nR; r++)  { 
		sum_s= 0;
		for( s= 1; s<=_nS; s++)  { 
		prod= 1;
		for( a= 1; a<=_nA; a++)  { 
			prod= prod*pow_1(condprobx[s][a][w],dataB[u][a][r]);
			} 
		prod= prod*margprobx[s][u][w];
		prod= prod/(ro[u][v][w]*patS[s][v]+(1.0-ro[u][v][w])*(1.0-patS[s][v]));
		prod= prod*(2.0*patS[s][v]-1.0);
		sum_s= sum_s+prod;
		} 
	sum_r= sum_r+sum_s*ga[w]*(1.0/p_r[r])*p_r_t[r][w]*(1.0/p_o_r_t[u][r][w]);
		} 
	return  (sum_r+(1.0/(double(_nO*_nT)))*(1.0/ro[u][v][w])-(1.0/(double(_nO*_nT)))*(1.0/(1.0-ro[u][v][w])));




} 

// function to compute probabilities needed for computing gradient of parameters (model M_1)
// p_r represents probability of data for person r 
// p_r_t represents probability of data for person r given membership of class t
// p_o_r_t represents probability of data for person r on object o given membership of class t
void calculate_probmat_gradient(extended1 ga, extended2 condprobx, extended3 margprobx,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat)
{

	unsigned short int r,t,o,s,a;
	long double prodobject,sumpattern,prod;
  
	for( t= 1; t<=_nT; t++)  { 
	for( r= 1; r<=_nR; r++)  { 
	prodobject= 1;
	for( o= 1; o<=_nO; o++)  { 
		sumpattern= 0;
		for( s= 1; s<=_nS; s++)  { 
		prod= 1;
		for( a= 1; a<=_nA; a++)  { 
		prod= prod*pow_1(condprobx[s][a],dataB[o][a][r]);
		} 
		prod= prod*margprobx[s][o][t];
		sumpattern= sumpattern+prod;
		} 
		p_o_r_t[o][r][t]= sumpattern;
	prodobject= prodobject*sumpattern;
		} 
	p_r_t[r][t]= prodobject;
	probmat[t][r]= prodobject*ga[t];
	} 
	} 

	for( r= 1; r<=_nR; r++)  { 
	p_r[r]= 0;
	for( t= 1; t<=_nT; t++) p_r[r]= p_r[r]+probmat[t][r];
	} 

 
}

// function to compute probabilities needed for computing gradient of object parameters (model M_2)
// p_r represents probability of data for person r 
// p_r_t represents probability of data for person r given membership of class t
// p_o_r_t represents probability of data for person r on object o given membership of class t
void calculate_probmat_gradient(extended1 ga, extended3 condprobx, extended2 margprobx,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat)
{

	unsigned short int r,t,o,s,a;
	long double prodobject,sumpattern,prod;
  
	for( t= 1; t<=_nT; t++)  { 
	for( r= 1; r<=_nR; r++)  { 
	prodobject= 1;
	for( o= 1; o<=_nO; o++)  { 
		sumpattern= 0;
		for( s= 1; s<=_nS; s++)  { 
		prod= 1;
		for( a= 1; a<=_nA; a++)  { 
		prod= prod*pow_1(condprobx[s][a][t],dataB[o][a][r]);
		} 
		prod= prod*margprobx[s][o];
		sumpattern= sumpattern+prod;
		} 
		p_o_r_t[o][r][t]= sumpattern;
	prodobject= prodobject*sumpattern;
		} 
	p_r_t[r][t]= prodobject;
	probmat[t][r]= prodobject*ga[t];
	} 
	} 

	for( r= 1; r<=_nR; r++)  { 
	p_r[r]= 0;
	for( t= 1; t<=_nT; t++) p_r[r]= p_r[r]+probmat[t][r];
	} 
}

// function to compute probabilities needed for computing gradient of object parameters (model M_3)
// p_r represents probability of data for person r 
// p_r_t represents probability of data for person r given membership of class t
// p_o_r_t represents probability of data for person r on object o given membership of class t
void calculate_probmat_gradient(extended1 ga, extended3 condprobx, extended3 margprobx,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat)
{

	unsigned short int r,t,o,s,a;
	long double prodobject,sumpattern,prod;
  
	for( t= 1; t<=_nT; t++)  { 
	for( r= 1; r<=_nR; r++)  { 
	prodobject= 1;
	for( o= 1; o<=_nO; o++)  { 
		sumpattern= 0;
		for( s= 1; s<=_nS; s++)  { 
		prod= 1;
		for( a= 1; a<=_nA; a++)  { 
		prod= prod*pow_1(condprobx[s][a][t],dataB[o][a][r]);
		} 
		prod= prod*margprobx[s][o][t];
		sumpattern= sumpattern+prod;
		} 
		p_o_r_t[o][r][t]= sumpattern;
	prodobject= prodobject*sumpattern;
		} 
	p_r_t[r][t]= prodobject;
	probmat[t][r]= prodobject*ga[t];
	} 
	} 

	for( r= 1; r<=_nR; r++)  { 
	p_r[r]= 0;
	for( t= 1; t<=_nT; t++) p_r[r]= p_r[r]+probmat[t][r];
	} 

}

// function to compute standard errors of class size parameters for model M_1 using method of finite differences
// disjunctive model
void calculate_se_ga_DC(extended3 ro_n,extended2 ta_n,extended1 ga_n,extended2 condprobx_n,extended3 margprobx_n,extended1 se_ga,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat)
{

	unsigned short int r,u;
	long double     sum;
  
	calculate_margprobX(ro_n,margprobx_n);
	calculate_condprobX_DC(ta_n,condprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n,p_o_r_t , p_r_t, p_r, probmat);

	for( u= 1; u<=_nT; u++)  { 
	sum= 0;
	for( r= 1; r<=_nR; r++)  { 
	sum= sum+(-1.0)*(1.0/((p_r[r])*(p_r[r])))*(p_r_t[r][u])*(p_r_t[r][u]);
	} 
	se_ga[u]= 1.0/sqrt(-sum);
	} 
}

// function to compute standard errors of class size parameters for model M_2 using method of finite differences
// disjunctive model
void calculate_se_ga_DC(extended2 ro_n,extended3 ta_n,extended1 ga_n,extended3 condprobx_n,extended2 margprobx_n,extended1 se_ga,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat)
{

	unsigned short int r,u;
	long double     sum;
  
	calculate_margprobX(ro_n,margprobx_n);
	calculate_condprobX_DC(ta_n,condprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n,p_o_r_t , p_r_t, p_r, probmat);

	for( u= 1; u<=_nT; u++)  { 
	sum= 0;
	for( r= 1; r<=_nR; r++)  { 
	sum= sum+(-1.0)*(1.0/((p_r[r])*(p_r[r])))*(p_r_t[r][u])*(p_r_t[r][u]);
	} 
	se_ga[u]= 1.0/sqrt(-sum);
	} 
}

// function to compute standard errors of class size parameters for model M_3 using method of finite differences
// disjunctive model
void calculate_se_ga_DC(extended3 ro_n,extended3 ta_n,extended1 ga_n,extended3 condprobx_n,extended3 margprobx_n,extended1 se_ga,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat)
{

	unsigned short int r,u;
	long double     sum;
  
	calculate_margprobX(ro_n,margprobx_n);
	calculate_condprobX_DC(ta_n,condprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n,p_o_r_t , p_r_t, p_r, probmat);

	for( u= 1; u<=_nT; u++)  { 
	sum= 0;
	for( r= 1; r<=_nR; r++)  { 
	sum= sum+(-1.0)*(1.0/((p_r[r])*(p_r[r])))*(p_r_t[r][u])*(p_r_t[r][u]);
	} 
	se_ga[u]= 1.0/sqrt(-sum);
	} 
}

// function to compute standard errors of attribute parameters for model M_1 using method of finite differences
// disjunctive model
void calculate_se_ta_DC(extended3 ro_n,extended2 ta_n,extended1 ga_n,extended2 condprobx_n,extended3 margprobx_n,extended2 se_ta,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat, extended2 ta_update)
{

	unsigned short int u,v;
	long double     up,low;

  
	calculate_margprobX(ro_n,margprobx_n);

	for( u= 1; u<=_nA; u++) { 
	for( v= 1; v<=_nF; v++) {
    
    if (ta_update[u][v]==1)
    {
           
	// upper
		ta_n[u][v]= ta_n[u][v]+delta;
		calculate_condprobX_DC(ta_n,condprobx_n);
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n,p_o_r_t , p_r_t, p_r, probmat);
		up= element_gradient_ta_DC(u,v,ta_n,ga_n,condprobx_n,margprobx_n,p_o_r_t , p_r_t, p_r);

	// lower
		ta_n[u][v]= ta_n[u][v]-2.0*delta;
		calculate_condprobX_DC(ta_n,condprobx_n);
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n,p_o_r_t , p_r_t, p_r, probmat);
		low= element_gradient_ta_DC(u,v,ta_n,ga_n,condprobx_n,margprobx_n,p_o_r_t , p_r_t, p_r);

		se_ta[u][v]= 1.0/sqrt(-(up-low)/(2.0*delta));
		ta_n[u][v]= ta_n[u][v]+delta;
    } 
	}
    } 
}

// function to compute standard errors of attribute parameters for model M_2 using method of finite differences
// disjunctive model
void calculate_se_ta_DC(extended2 ro_n,extended3 ta_n,extended1 ga_n,extended3 condprobx_n,extended2 margprobx_n,extended3 se_ta,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat, extended3 ta_update)
{

	unsigned short int u,v,w;
	long double     up,low;

	calculate_margprobX(ro_n,margprobx_n);

	for( u= 1; u<=_nA; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
	for( w= 1; w<=_nT; w++)  { 
    
    if (ta_update[u][v][w]==1)
    {
    
	// upper
		ta_n[u][v][w]= ta_n[u][v][w]+delta;
		calculate_condprobX_DC(ta_n,condprobx_n);
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t, p_r_t, p_r,probmat);
		up= element_gradient_ta_DC(u,v,w,ta_n,ga_n,condprobx_n,margprobx_n, p_o_r_t, p_r_t, p_r);

	// lower
		ta_n[u][v][w]= ta_n[u][v][w]-2.0*delta;
		calculate_condprobX_DC(ta_n,condprobx_n);
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t, p_r_t, p_r,probmat);
		low= element_gradient_ta_DC(u,v,w,ta_n,ga_n,condprobx_n,margprobx_n, p_o_r_t, p_r_t, p_r);

		se_ta[u][v][w]= 1.0/sqrt(-(up-low)/(2.0*delta));
		ta_n[u][v][w]= ta_n[u][v][w]+delta;
    } 
	} 
	}
    } 
		
}

// function to compute standard errors of attribute parameters for model M_3 using method of finite differences
// disjunctive model
void calculate_se_ta_DC(extended3 ro_n,extended3 ta_n,extended1 ga_n,extended3 condprobx_n,extended3 margprobx_n,extended3 se_ta,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat, extended3 ta_update)
{

	unsigned short int u,v,w;
	long double     up,low;

	calculate_margprobX(ro_n,margprobx_n);

	for( u= 1; u<=_nA; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
	for( w= 1; w<=_nT; w++)  { 

    if (ta_update[u][v][w]==1)
    {
    
	// upper
		ta_n[u][v][w]= ta_n[u][v][w]+delta;
		calculate_condprobX_DC(ta_n,condprobx_n);
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t, p_r_t, p_r,probmat);
		up= element_gradient_ta_DC(u,v,w,ta_n,ga_n,condprobx_n,margprobx_n, p_o_r_t, p_r_t, p_r);

	// lower
		ta_n[u][v][w]= ta_n[u][v][w]-2.0*delta;
		calculate_condprobX_DC(ta_n,condprobx_n);
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t, p_r_t, p_r,probmat);
		low= element_gradient_ta_DC(u,v,w,ta_n,ga_n,condprobx_n,margprobx_n, p_o_r_t, p_r_t, p_r);

		se_ta[u][v][w]= 1.0/sqrt(-(up-low)/(2.0*delta));
		ta_n[u][v][w]= ta_n[u][v][w]+delta;
    } 
	} 
	}
    } 
		
}

// function to compute standard errors of object parameters for model M_1 using method of finite differences
// disjunctive model
void calculate_se_ro_DC(extended3 ro_n, extended2 ta_n,extended1 ga_n,extended2 condprobx_n,extended3 margprobx_n,extended3 se_ro,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat, extended3 ro_update)
{

	unsigned short int u,v,w;
	long double     up,low;
	
	calculate_condprobX_DC(ta_n,condprobx_n);

	for( u= 1; u<=_nO; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
	for( w= 1; w<=_nT; w++)  { 
    
    if (ro_update[u][v][w]==1)
    { 
	// upper
	ro_n[u][v][w]= ro_n[u][v][w]+delta;
	calculate_margprobX(ro_n,margprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t , p_r_t, p_r, probmat);
	up= element_gradient_ro(u,v,w,ro_n,ga_n,margprobx_n,condprobx_n,p_o_r_t , p_r_t, p_r);

	// lower
	ro_n[u][v][w]= ro_n[u][v][w]-2.0*delta;
	calculate_margprobX(ro_n,margprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t , p_r_t, p_r, probmat);
	low= element_gradient_ro(u,v,w,ro_n,ga_n,margprobx_n,condprobx_n,p_o_r_t , p_r_t, p_r);

	se_ro[u][v][w]= 1.0/sqrt(-(up-low)/(2.0*delta));
	ro_n[u][v][w]= ro_n[u][v][w]+delta;
    }
	} 
	} 
	} 
}

// function to compute standard errors of object parameters for model M_2 using method of finite differences
// disjunctive model
void calculate_se_ro_DC(extended2 ro_n,extended3 ta_n,extended1 ga_n,extended3 condprobx_n,extended2 margprobx_n,extended2 se_ro,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat, extended2 ro_update)
{

	unsigned short int u,v;
	long double     up,low;
	
	calculate_condprobX_DC(ta_n,condprobx_n);

	for( u= 1; u<=_nO; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
       
     if (ro_update[u][v]==1)
     {
      
	// upper
	ro_n[u][v]= ro_n[u][v]+delta;
	calculate_margprobX(ro_n,margprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t, p_r_t, p_r,probmat);
	up= element_gradient_ro(u,v,ro_n,ga_n,margprobx_n,condprobx_n,p_o_r_t,p_r_t,p_r);

	// lower
	ro_n[u][v]= ro_n[u][v]-2.0*delta;
	calculate_margprobX(ro_n,margprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t, p_r_t, p_r,probmat);
	low= element_gradient_ro(u,v,ro_n,ga_n,margprobx_n,condprobx_n,p_o_r_t,p_r_t,p_r);

	se_ro[u][v]= 1.0/sqrt(-(up-low)/(2.0*delta));
	ro_n[u][v]= ro_n[u][v]+delta;
	} 
	} 
    }

}

// function to compute standard errors of object parameters for model M_3 using method of finite differences
// disjunctive model
void calculate_se_ro_DC(extended3 ro_n, extended3 ta_n,extended1 ga_n,extended3 condprobx_n,extended3 margprobx_n,extended3 se_ro,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat, extended3 ro_update)
{

	unsigned short int u,v,w;
	long double     up,low;
	
	calculate_condprobX_DC(ta_n,condprobx_n);

	for( u= 1; u<=_nO; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
	for( w= 1; w<=_nT; w++)  { 
    
    if (ro_update[u][v][w]==1) 
    {
	// upper
	ro_n[u][v][w]= ro_n[u][v][w]+delta;
	calculate_margprobX(ro_n,margprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t, p_r_t, p_r,probmat);
	up= element_gradient_ro(u,v,w,ro_n,ga_n,margprobx_n,condprobx_n, p_o_r_t, p_r_t, p_r);

	// lower
	ro_n[u][v][w]= ro_n[u][v][w]-2.0*delta;
	calculate_margprobX(ro_n,margprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t, p_r_t, p_r,probmat);
	low= element_gradient_ro(u,v,w,ro_n,ga_n,margprobx_n,condprobx_n, p_o_r_t, p_r_t, p_r);

	se_ro[u][v][w]= 1.0/sqrt(-(up-low)/(2.0*delta));
	ro_n[u][v][w]= ro_n[u][v][w]+delta;
	} 
	} 
	}
    } 


}


// function to compute standard errors of class size parameters for model M_1 using method of finite differences
// additive model
void calculate_se_ga_ADD(extended3 ro_n,extended2 ta_n,extended1 ga_n,extended2 condprobx_n,extended3 margprobx_n,extended1 se_ga,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat)
{

	unsigned short int r,u;
	long double     sum;
  
	calculate_margprobX(ro_n,margprobx_n);
	calculate_condprobX_ADD(ta_n,condprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n,p_o_r_t , p_r_t, p_r, probmat);

	for( u= 1; u<=_nT; u++)  { 
	sum= 0;
	for( r= 1; r<=_nR; r++)  { 
	sum= sum+(-1.0)*(1.0/((p_r[r])*(p_r[r])))*(p_r_t[r][u])*(p_r_t[r][u]);
	} 
	se_ga[u]= 1.0/sqrt(-sum);
	} 
}

// function to compute standard errors of class size parameters for model M_2 using method of finite differences
// additive model
void calculate_se_ga_ADD(extended2 ro_n,extended3 ta_n,extended1 ga_n,extended3 condprobx_n,extended2 margprobx_n,extended1 se_ga,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat)
{

	unsigned short int r,u;
	long double     sum;
  
	calculate_margprobX(ro_n,margprobx_n);
	calculate_condprobX_ADD(ta_n,condprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n,p_o_r_t , p_r_t, p_r, probmat);

	for( u= 1; u<=_nT; u++)  { 
	sum= 0;
	for( r= 1; r<=_nR; r++)  { 
	sum= sum+(-1.0)*(1.0/((p_r[r])*(p_r[r])))*(p_r_t[r][u])*(p_r_t[r][u]);
	} 
	se_ga[u]= 1.0/sqrt(-sum);
	} 
}

// function to compute standard errors of class size parameters for model M_3 using method of finite differences
// additive model
void calculate_se_ga_ADD(extended3 ro_n,extended3 ta_n,extended1 ga_n,extended3 condprobx_n,extended3 margprobx_n,extended1 se_ga,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat)
{

	unsigned short int r,u;
	long double     sum;
  
	calculate_margprobX(ro_n,margprobx_n);
	calculate_condprobX_ADD(ta_n,condprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n,p_o_r_t , p_r_t, p_r, probmat);

	for( u= 1; u<=_nT; u++)  { 
	sum= 0;
	for( r= 1; r<=_nR; r++)  { 
	sum= sum+(-1.0)*(1.0/((p_r[r])*(p_r[r])))*(p_r_t[r][u])*(p_r_t[r][u]);
	} 
	se_ga[u]= 1.0/sqrt(-sum);
	} 
}

// function to compute standard errors of attribute parameters for model M_1 using method of finite differences
// additive model
void calculate_se_ta_ADD(extended3 ro_n,extended2 ta_n,extended1 ga_n,extended2 condprobx_n,extended3 margprobx_n,extended2 se_ta,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat, extended2 ta_update)
{

	unsigned short int u,v;
	long double     up,low;

  
	calculate_margprobX(ro_n,margprobx_n);

	for( u= 1; u<=_nA; u++) { 
	for( v= 1; v<=_nF; v++) {
    
    if (ta_update[u][v]==1)
    {
           
	// upper
		ta_n[u][v]= ta_n[u][v]+delta;
		calculate_condprobX_ADD(ta_n,condprobx_n);
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n,p_o_r_t , p_r_t, p_r, probmat);
		up= element_gradient_ta_ADD(u,v,ta_n,ga_n,condprobx_n,margprobx_n,p_o_r_t , p_r_t, p_r);

	// lower
		ta_n[u][v]= ta_n[u][v]-2.0*delta;
		calculate_condprobX_ADD(ta_n,condprobx_n);
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n,p_o_r_t , p_r_t, p_r, probmat);
		low= element_gradient_ta_ADD(u,v,ta_n,ga_n,condprobx_n,margprobx_n,p_o_r_t , p_r_t, p_r);

		se_ta[u][v]= 1.0/sqrt(-(up-low)/(2.0*delta));
		ta_n[u][v]= ta_n[u][v]+delta;
    } 
	}
    } 
}

// function to compute standard errors of attribute parameters for model M_2 using method of finite differences
// additive model
void calculate_se_ta_ADD(extended2 ro_n,extended3 ta_n,extended1 ga_n,extended3 condprobx_n,extended2 margprobx_n,extended3 se_ta,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat, extended3 ta_update)
{

	unsigned short int u,v,w;
	long double     up,low;

	calculate_margprobX(ro_n,margprobx_n);

	for( u= 1; u<=_nA; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
	for( w= 1; w<=_nT; w++)  { 
    
    if (ta_update[u][v][w]==1)
    {
    
	// upper
		ta_n[u][v][w]= ta_n[u][v][w]+delta;
		calculate_condprobX_ADD(ta_n,condprobx_n);
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t, p_r_t, p_r,probmat);
		up= element_gradient_ta_ADD(u,v,w,ta_n,ga_n,condprobx_n,margprobx_n, p_o_r_t, p_r_t, p_r);

	// lower
		ta_n[u][v][w]= ta_n[u][v][w]-2.0*delta;
		calculate_condprobX_ADD(ta_n,condprobx_n);
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t, p_r_t, p_r,probmat);
		low= element_gradient_ta_ADD(u,v,w,ta_n,ga_n,condprobx_n,margprobx_n, p_o_r_t, p_r_t, p_r);

		se_ta[u][v][w]= 1.0/sqrt(-(up-low)/(2.0*delta));
		ta_n[u][v][w]= ta_n[u][v][w]+delta;
    } 
	} 
	}
    } 
		
}

// function to compute standard errors of attribute parameters for model M_3 using method of finite differences
// additive model
void calculate_se_ta_ADD(extended3 ro_n,extended3 ta_n,extended1 ga_n,extended3 condprobx_n,extended3 margprobx_n,extended3 se_ta,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat, extended3 ta_update)
{

	unsigned short int u,v,w;
	long double     up,low;

	calculate_margprobX(ro_n,margprobx_n);

	for( u= 1; u<=_nA; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
	for( w= 1; w<=_nT; w++)  { 

    if (ta_update[u][v][w]==1)
    {
    
	// upper
		ta_n[u][v][w]= ta_n[u][v][w]+delta;
		calculate_condprobX_ADD(ta_n,condprobx_n);
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t, p_r_t, p_r,probmat);
		up= element_gradient_ta_ADD(u,v,w,ta_n,ga_n,condprobx_n,margprobx_n, p_o_r_t, p_r_t, p_r);

	// lower
		ta_n[u][v][w]= ta_n[u][v][w]-2.0*delta;
		calculate_condprobX_ADD(ta_n,condprobx_n);
		calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t, p_r_t, p_r,probmat);
		low= element_gradient_ta_ADD(u,v,w,ta_n,ga_n,condprobx_n,margprobx_n, p_o_r_t, p_r_t, p_r);

		se_ta[u][v][w]= 1.0/sqrt(-(up-low)/(2.0*delta));
		ta_n[u][v][w]= ta_n[u][v][w]+delta;
    } 
	} 
	}
    } 
		
}

// function to compute standard errors of object parameters for model M_1 using method of finite differences
// additive model
void calculate_se_ro_ADD(extended3 ro_n, extended2 ta_n,extended1 ga_n,extended2 condprobx_n,extended3 margprobx_n,extended3 se_ro,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat, extended3 ro_update)
{

	unsigned short int u,v,w;
	long double     up,low;
	
	calculate_condprobX_ADD(ta_n,condprobx_n);

	for( u= 1; u<=_nO; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
	for( w= 1; w<=_nT; w++)  { 
    
    if (ro_update[u][v][w]==1)
    { 
	// upper
	ro_n[u][v][w]= ro_n[u][v][w]+delta;
	calculate_margprobX(ro_n,margprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t , p_r_t, p_r, probmat);
	up= element_gradient_ro(u,v,w,ro_n,ga_n,margprobx_n,condprobx_n,p_o_r_t , p_r_t, p_r);

	// lower
	ro_n[u][v][w]= ro_n[u][v][w]-2.0*delta;
	calculate_margprobX(ro_n,margprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t , p_r_t, p_r, probmat);
	low= element_gradient_ro(u,v,w,ro_n,ga_n,margprobx_n,condprobx_n,p_o_r_t , p_r_t, p_r);

	se_ro[u][v][w]= 1.0/sqrt(-(up-low)/(2.0*delta));
	ro_n[u][v][w]= ro_n[u][v][w]+delta;
    }
	} 
	} 
	} 
}

// function to compute standard errors of object parameters for model M_2 using method of finite differences
// additive model
void calculate_se_ro_ADD(extended2 ro_n,extended3 ta_n,extended1 ga_n,extended3 condprobx_n,extended2 margprobx_n,extended2 se_ro,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat, extended2 ro_update)
{

	unsigned short int u,v;
	long double     up,low;
	
	calculate_condprobX_ADD(ta_n,condprobx_n);

	for( u= 1; u<=_nO; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
       
     if (ro_update[u][v]==1)
     {
      
	// upper
	ro_n[u][v]= ro_n[u][v]+delta;
	calculate_margprobX(ro_n,margprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t, p_r_t, p_r,probmat);
	up= element_gradient_ro(u,v,ro_n,ga_n,margprobx_n,condprobx_n,p_o_r_t,p_r_t,p_r);

	// lower
	ro_n[u][v]= ro_n[u][v]-2.0*delta;
	calculate_margprobX(ro_n,margprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t, p_r_t, p_r,probmat);
	low= element_gradient_ro(u,v,ro_n,ga_n,margprobx_n,condprobx_n,p_o_r_t,p_r_t,p_r);

	se_ro[u][v]= 1.0/sqrt(-(up-low)/(2.0*delta));
	ro_n[u][v]= ro_n[u][v]+delta;
	} 
	} 
    }

}

// function to compute standard errors of object parameters for model M_3 using method of finite differences
// additive model
void calculate_se_ro_ADD(extended3 ro_n, extended3 ta_n,extended1 ga_n,extended3 condprobx_n,extended3 margprobx_n,extended3 se_ro,extended3 p_o_r_t ,extended2 p_r_t,extended1 p_r,extended2 probmat, extended3 ro_update)
{

	unsigned short int u,v,w;
	long double     up,low;
	
	calculate_condprobX_ADD(ta_n,condprobx_n);

	for( u= 1; u<=_nO; u++)  { 
	for( v= 1; v<=_nF; v++)  { 
	for( w= 1; w<=_nT; w++)  { 
    
    if (ro_update[u][v][w]==1) 
    {
	// upper
	ro_n[u][v][w]= ro_n[u][v][w]+delta;
	calculate_margprobX(ro_n,margprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t, p_r_t, p_r,probmat);
	up= element_gradient_ro(u,v,w,ro_n,ga_n,margprobx_n,condprobx_n, p_o_r_t, p_r_t, p_r);

	// lower
	ro_n[u][v][w]= ro_n[u][v][w]-2.0*delta;
	calculate_margprobX(ro_n,margprobx_n);
	calculate_probmat_gradient(ga_n, condprobx_n, margprobx_n, p_o_r_t, p_r_t, p_r,probmat);
	low= element_gradient_ro(u,v,w,ro_n,ga_n,margprobx_n,condprobx_n, p_o_r_t, p_r_t, p_r);

	se_ro[u][v][w]= 1.0/sqrt(-(up-low)/(2.0*delta));
	ro_n[u][v][w]= ro_n[u][v][w]+delta;
	} 
	} 
	}
    } 


}



// bootstrap model fixed_object_XZ_Y
// disjunctive model
void bootstrap_dependency_XZ_Y_DC(binary3 dataobs, extended3 ro, extended2 ta, extended1 ga, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att,extended2 report_common_obj)
{
	// bootstrap vars
	extended2 OR_common_att_obs,OR_common_obj_obs;
	extended2 mean_OR_common_att;
	extended2 mean_OR_common_obj;
	binary3 datarep;
	int N_attpair, N_objpair;

	N_attpair = (_nA*(_nA-1)) / 2; 
	N_objpair = (_nO*(_nO-1)) / 2;

	// construct OR_common_att_obs
	OR_common_att_obs = new long double * [_nO+1]; 
	for (int i = 0; i < _nO+1; i++)
	{
		OR_common_att_obs[i] = new long double [N_attpair+1];
	}
	// construct OR_common_obj_obs
	OR_common_obj_obs = new long double * [_nA+1]; 
	for (int i = 0; i < _nA+1; i++)
	{
		OR_common_obj_obs[i] = new long double [N_objpair+1];
	}



	// construct datarep
	datarep = new bool ** [_nO+1];
	for (int i = 0; i < _nO+1; i++)
	{
		datarep[i] = new bool * [_nA+1];
		for (int ii = 0; ii < _nA+1; ii++)
		{
			datarep[i][ii] = new bool [_nR+1];
		}
	}

	//construct mean_OR_common_att & mean_OR_common_obj
	mean_OR_common_att = new long double * [_nO+1] ;
	for (int i = 0; i < _nO+1; i++)
	{
		mean_OR_common_att[i] = new long double [N_attpair+1];
	}
	mean_OR_common_obj = new long double * [_nA+1] ;
	for (int i = 0; i < _nA+1; i++)
	{
		mean_OR_common_obj[i] = new long double [N_objpair+1];
	}



	//  --------------------------------------  STAP 1  ----------------------------------------------

	compute_OR_common_attribute_obs(dataobs, OR_common_att_obs);
	compute_OR_common_object_obs(dataobs, OR_common_obj_obs);


	//  --------------------------------------  STAP 2  ----------------------------------------------

	
	for (int iter = 1; iter <= _Nbootstrap; iter++)
	{ 
		//if (iter mod 500)=0 then writeln(iter:5);
		generate_data_fixedobject_XZ_Y_DC( ga, ro, ta, datarep); 
		compute_OR_common_attribute_rep( datarep, iter, OR_common_att_rep );
		compute_OR_common_object_rep( datarep, iter, OR_common_obj_rep );
	}
	

	//  --------------------------------------  STAP 3  ----------------------------------------------

	compute_mean_OR(OR_common_att_rep, _Nbootstrap,_nO, N_attpair,mean_OR_common_att);
	compute_mean_OR(OR_common_obj_rep, _Nbootstrap,_nA, N_objpair,mean_OR_common_obj);

	
	create_report_common_att(OR_common_att_obs, mean_OR_common_att,  report_common_att);
	create_report_common_obj(OR_common_obj_obs, mean_OR_common_obj,  report_common_obj);


 // destruct
 
     // OR_common_att_obs
 		for (int i = 0; i < _nO+1; i++)
		{
			delete OR_common_att_obs[i];
		}
		delete[] OR_common_att_obs;
	
	// OR_common_obj_obs	
		for (int i = 0; i < _nA+1; i++)
		{
			delete OR_common_obj_obs[i];
		}
		delete[] OR_common_obj_obs;
		
		//delete cumul;
		//delete latklas;	
		
	// datarep	
	  for (int i = 0; i < _nO+1; i++)
		{
			for (int ii = 0; ii < _nA+1; ii++)
			{
				delete datarep[i][ii];
			}
			delete datarep[i];
		}
		delete[] datarep;

      // mean_OR_common_att
	  for (int i = 0; i < _nO+1; i++)
		{
			delete mean_OR_common_att[i];
		}
		delete[] mean_OR_common_att;
		
	  // mean_OR_common_obj	
	  for (int i = 0; i < _nA+1; i++)
		{
			delete mean_OR_common_obj[i];
		}
		delete[] mean_OR_common_obj;
 

}


// bootstrap model fixed_object_X_YZ
// disjunctive model
void bootstrap_dependency_X_YZ_DC(binary3 dataobs, extended2 ro, extended3 ta, extended1 ga, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att,extended2 report_common_obj)
{
	// bootstrap vars
	extended2 OR_common_att_obs,OR_common_obj_obs;
	extended2 mean_OR_common_att;
	extended2 mean_OR_common_obj;
	binary3 datarep;
	int N_attpair, N_objpair;

	N_attpair = (_nA*(_nA-1)) / 2; 
	N_objpair = (_nO*(_nO-1)) / 2;

	// construct OR_common_att_obs
	OR_common_att_obs = new long double * [_nO+1]; 
	for (int i = 0; i < _nO+1; i++)
	{
		OR_common_att_obs[i] = new long double [N_attpair+1];
	}
	// construct OR_common_obj_obs
	OR_common_obj_obs = new long double * [_nA+1]; 
	for (int i = 0; i < _nA+1; i++)
	{
		OR_common_obj_obs[i] = new long double [N_objpair+1];
	}


	// construct datarep
	datarep = new bool ** [_nO+1];
	for (int i = 0; i < _nO+1; i++)
	{
		datarep[i] = new bool * [_nA+1];
		for (int ii = 0; ii < _nA+1; ii++)
		{
			datarep[i][ii] = new bool [_nR+1];
		}
	}

	//construct mean_OR_common_att & mean_OR_common_obj
	mean_OR_common_att = new long double * [_nO+1] ;
	for (int i = 0; i < _nO+1; i++)
	{
		mean_OR_common_att[i] = new long double [N_attpair+1];
	}
	mean_OR_common_obj = new long double * [_nA+1] ;
	for (int i = 0; i < _nA+1; i++)
	{
		mean_OR_common_obj[i] = new long double [N_objpair+1];
	}



	//  --------------------------------------  STAP 1  ----------------------------------------------

	compute_OR_common_attribute_obs(dataobs, OR_common_att_obs);
	compute_OR_common_object_obs(dataobs, OR_common_obj_obs);


	//  --------------------------------------  STAP 2  ----------------------------------------------

	
	for (int iter = 1; iter <= _Nbootstrap; iter++)
	{ 
		//if (iter mod 500)=0 then writeln(iter:5);
		generate_data_fixedobject_X_YZ_DC( ga, ro, ta, datarep); 
		compute_OR_common_attribute_rep( datarep, iter, OR_common_att_rep );
		compute_OR_common_object_rep( datarep, iter, OR_common_obj_rep );
	}
	

	//  --------------------------------------  STAP 3  ----------------------------------------------

	compute_mean_OR(OR_common_att_rep, _Nbootstrap,_nO, N_attpair,mean_OR_common_att);
	compute_mean_OR(OR_common_obj_rep, _Nbootstrap,_nA, N_objpair,mean_OR_common_obj);

	
	create_report_common_att(OR_common_att_obs, mean_OR_common_att,  report_common_att);
	create_report_common_obj(OR_common_obj_obs, mean_OR_common_obj,  report_common_obj);


 // destruct
 
     // OR_common_att_obs
 		for (int i = 0; i < _nO+1; i++)
		{
			delete OR_common_att_obs[i];
		}
		delete[] OR_common_att_obs;
	
	// OR_common_obj_obs	
		for (int i = 0; i < _nA+1; i++)
		{
			delete OR_common_obj_obs[i];
		}
		delete[] OR_common_obj_obs;
		
		//delete cumul;
		//delete latklas;	
		
	// datarep	
	  for (int i = 0; i < _nO+1; i++)
		{
			for (int ii = 0; ii < _nA+1; ii++)
			{
				delete datarep[i][ii];
			}
			delete datarep[i];
		}
		delete[] datarep;

      // mean_OR_common_att
	  for (int i = 0; i < _nO+1; i++)
		{
			delete mean_OR_common_att[i];
		}
		delete[] mean_OR_common_att;
		
	  // mean_OR_common_obj	
	  for (int i = 0; i < _nA+1; i++)
		{
			delete mean_OR_common_obj[i];
		}
		delete[] mean_OR_common_obj;
 

}




// bootstrap model fixed_object_XZ_YZ
// disjunctive model
void bootstrap_dependency_XZ_YZ_DC(binary3 dataobs, extended3 ro, extended3 ta, extended1 ga, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att,extended2 report_common_obj)
{
	// bootstrap vars
	extended2 OR_common_att_obs,OR_common_obj_obs;
	extended2 mean_OR_common_att;
	extended2 mean_OR_common_obj;
	binary3 datarep;
	int N_attpair, N_objpair;

	N_attpair = (_nA*(_nA-1)) / 2; 
	N_objpair = (_nO*(_nO-1)) / 2;

	// construct OR_common_att_obs
	OR_common_att_obs = new long double * [_nO+1]; 
	for (int i = 0; i < _nO+1; i++)
	{
		OR_common_att_obs[i] = new long double [N_attpair+1];
	}
	// construct OR_common_obj_obs
	OR_common_obj_obs = new long double * [_nA+1]; 
	for (int i = 0; i < _nA+1; i++)
	{
		OR_common_obj_obs[i] = new long double [N_objpair+1];
	}


	// construct datarep
	datarep = new bool ** [_nO+1];
	for (int i = 0; i < _nO+1; i++)
	{
		datarep[i] = new bool * [_nA+1];
		for (int ii = 0; ii < _nA+1; ii++)
		{
			datarep[i][ii] = new bool [_nR+1];
		}
	}

	//construct mean_OR_common_att & mean_OR_common_obj
	mean_OR_common_att = new long double * [_nO+1] ;
	for (int i = 0; i < _nO+1; i++)
	{
		mean_OR_common_att[i] = new long double [N_attpair+1];
	}
	mean_OR_common_obj = new long double * [_nA+1] ;
	for (int i = 0; i < _nA+1; i++)
	{
		mean_OR_common_obj[i] = new long double [N_objpair+1];
	}



	//  --------------------------------------  STAP 1  ----------------------------------------------

	compute_OR_common_attribute_obs(dataobs, OR_common_att_obs);
	compute_OR_common_object_obs(dataobs, OR_common_obj_obs);


	//  --------------------------------------  STAP 2  ----------------------------------------------

	
	for (int iter = 1; iter <= _Nbootstrap; iter++)
	{ 
		//if (iter mod 500)=0 then writeln(iter:5);
		generate_data_fixedobject_XZ_YZ_DC( ga, ro, ta, datarep); 
		compute_OR_common_attribute_rep( datarep, iter, OR_common_att_rep );
		compute_OR_common_object_rep( datarep, iter, OR_common_obj_rep );
	}
	

	//  --------------------------------------  STAP 3  ----------------------------------------------

	compute_mean_OR(OR_common_att_rep, _Nbootstrap,_nO, N_attpair,mean_OR_common_att);
	compute_mean_OR(OR_common_obj_rep, _Nbootstrap,_nA, N_objpair,mean_OR_common_obj);

	
	create_report_common_att(OR_common_att_obs, mean_OR_common_att,  report_common_att);
	create_report_common_obj(OR_common_obj_obs, mean_OR_common_obj,  report_common_obj);


 // destruct
 
     // OR_common_att_obs
 		for (int i = 0; i < _nO+1; i++)
		{
			delete OR_common_att_obs[i];
		}
		delete[] OR_common_att_obs;
	
	// OR_common_obj_obs	
		for (int i = 0; i < _nA+1; i++)
		{
			delete OR_common_obj_obs[i];
		}
		delete[] OR_common_obj_obs;
		
	
		
	// datarep	
	  for (int i = 0; i < _nO+1; i++)
		{
			for (int ii = 0; ii < _nA+1; ii++)
			{
				delete datarep[i][ii];
			}
			delete datarep[i];
		}
		delete[] datarep;

      // mean_OR_common_att
	  for (int i = 0; i < _nO+1; i++)
		{
			delete mean_OR_common_att[i];
		}
		delete[] mean_OR_common_att;
		
	  // mean_OR_common_obj	
	  for (int i = 0; i < _nA+1; i++)
		{
			delete mean_OR_common_obj[i];
		}
		delete[] mean_OR_common_obj;
 

}



// bootstrap model fixed_object_XZ_Y
// additive model
void bootstrap_dependency_XZ_Y_ADD(binary3 dataobs, extended3 ro, extended2 ta, extended1 ga, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att,extended2 report_common_obj)
{
	// bootstrap vars
	extended2 OR_common_att_obs,OR_common_obj_obs;
	extended2 mean_OR_common_att;
	extended2 mean_OR_common_obj;
	binary3 datarep;
	int N_attpair, N_objpair;

	N_attpair = (_nA*(_nA-1)) / 2; 
	N_objpair = (_nO*(_nO-1)) / 2;

	// construct OR_common_att_obs
	OR_common_att_obs = new long double * [_nO+1]; 
	for (int i = 0; i < _nO+1; i++)
	{
		OR_common_att_obs[i] = new long double [N_attpair+1];
	}
	// construct OR_common_obj_obs
	OR_common_obj_obs = new long double * [_nA+1]; 
	for (int i = 0; i < _nA+1; i++)
	{
		OR_common_obj_obs[i] = new long double [N_objpair+1];
	}



	// construct datarep
	datarep = new bool ** [_nO+1];
	for (int i = 0; i < _nO+1; i++)
	{
		datarep[i] = new bool * [_nA+1];
		for (int ii = 0; ii < _nA+1; ii++)
		{
			datarep[i][ii] = new bool [_nR+1];
		}
	}

	//construct mean_OR_common_att & mean_OR_common_obj
	mean_OR_common_att = new long double * [_nO+1] ;
	for (int i = 0; i < _nO+1; i++)
	{
		mean_OR_common_att[i] = new long double [N_attpair+1];
	}
	mean_OR_common_obj = new long double * [_nA+1] ;
	for (int i = 0; i < _nA+1; i++)
	{
		mean_OR_common_obj[i] = new long double [N_objpair+1];
	}



	//  --------------------------------------  STAP 1  ----------------------------------------------

	compute_OR_common_attribute_obs(dataobs, OR_common_att_obs);
	compute_OR_common_object_obs(dataobs, OR_common_obj_obs);


	//  --------------------------------------  STAP 2  ----------------------------------------------

	
	for (int iter = 1; iter <= _Nbootstrap; iter++)
	{ 
		//if (iter mod 500)=0 then writeln(iter:5);
		generate_data_fixedobject_XZ_Y_ADD( ga, ro, ta, datarep); 
		compute_OR_common_attribute_rep( datarep, iter, OR_common_att_rep );
		compute_OR_common_object_rep( datarep, iter, OR_common_obj_rep );
	}
	

	//  --------------------------------------  STAP 3  ----------------------------------------------

	compute_mean_OR(OR_common_att_rep, _Nbootstrap,_nO, N_attpair,mean_OR_common_att);
	compute_mean_OR(OR_common_obj_rep, _Nbootstrap,_nA, N_objpair,mean_OR_common_obj);

	
	create_report_common_att(OR_common_att_obs, mean_OR_common_att,  report_common_att);
	create_report_common_obj(OR_common_obj_obs, mean_OR_common_obj,  report_common_obj);


 // destruct
 
     // OR_common_att_obs
 		for (int i = 0; i < _nO+1; i++)
		{
			delete OR_common_att_obs[i];
		}
		delete[] OR_common_att_obs;
	
	// OR_common_obj_obs	
		for (int i = 0; i < _nA+1; i++)
		{
			delete OR_common_obj_obs[i];
		}
		delete[] OR_common_obj_obs;
		
		//delete cumul;
		//delete latklas;	
		
	// datarep	
	  for (int i = 0; i < _nO+1; i++)
		{
			for (int ii = 0; ii < _nA+1; ii++)
			{
				delete datarep[i][ii];
			}
			delete datarep[i];
		}
		delete[] datarep;

      // mean_OR_common_att
	  for (int i = 0; i < _nO+1; i++)
		{
			delete mean_OR_common_att[i];
		}
		delete[] mean_OR_common_att;
		
	  // mean_OR_common_obj	
	  for (int i = 0; i < _nA+1; i++)
		{
			delete mean_OR_common_obj[i];
		}
		delete[] mean_OR_common_obj;
 

}


// bootstrap model fixed_object_X_YZ
// additive model
void bootstrap_dependency_X_YZ_ADD(binary3 dataobs, extended2 ro, extended3 ta, extended1 ga, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att,extended2 report_common_obj)
{
	// bootstrap vars
	extended2 OR_common_att_obs,OR_common_obj_obs;
	extended2 mean_OR_common_att;
	extended2 mean_OR_common_obj;
	binary3 datarep;
	int N_attpair, N_objpair;

	N_attpair = (_nA*(_nA-1)) / 2; 
	N_objpair = (_nO*(_nO-1)) / 2;

	// construct OR_common_att_obs
	OR_common_att_obs = new long double * [_nO+1]; 
	for (int i = 0; i < _nO+1; i++)
	{
		OR_common_att_obs[i] = new long double [N_attpair+1];
	}
	// construct OR_common_obj_obs
	OR_common_obj_obs = new long double * [_nA+1]; 
	for (int i = 0; i < _nA+1; i++)
	{
		OR_common_obj_obs[i] = new long double [N_objpair+1];
	}


	// construct datarep
	datarep = new bool ** [_nO+1];
	for (int i = 0; i < _nO+1; i++)
	{
		datarep[i] = new bool * [_nA+1];
		for (int ii = 0; ii < _nA+1; ii++)
		{
			datarep[i][ii] = new bool [_nR+1];
		}
	}

	//construct mean_OR_common_att & mean_OR_common_obj
	mean_OR_common_att = new long double * [_nO+1] ;
	for (int i = 0; i < _nO+1; i++)
	{
		mean_OR_common_att[i] = new long double [N_attpair+1];
	}
	mean_OR_common_obj = new long double * [_nA+1] ;
	for (int i = 0; i < _nA+1; i++)
	{
		mean_OR_common_obj[i] = new long double [N_objpair+1];
	}



	//  --------------------------------------  STAP 1  ----------------------------------------------

	compute_OR_common_attribute_obs(dataobs, OR_common_att_obs);
	compute_OR_common_object_obs(dataobs, OR_common_obj_obs);


	//  --------------------------------------  STAP 2  ----------------------------------------------

	
	for (int iter = 1; iter <= _Nbootstrap; iter++)
	{ 
		//if (iter mod 500)=0 then writeln(iter:5);
		generate_data_fixedobject_X_YZ_ADD( ga, ro, ta, datarep); 
		compute_OR_common_attribute_rep( datarep, iter, OR_common_att_rep );
		compute_OR_common_object_rep( datarep, iter, OR_common_obj_rep );
	}
	

	//  --------------------------------------  STAP 3  ----------------------------------------------

	compute_mean_OR(OR_common_att_rep, _Nbootstrap,_nO, N_attpair, mean_OR_common_att);
	compute_mean_OR(OR_common_obj_rep, _Nbootstrap,_nA, N_objpair, mean_OR_common_obj);

	
	create_report_common_att(OR_common_att_obs, mean_OR_common_att,  report_common_att);
	create_report_common_obj(OR_common_obj_obs, mean_OR_common_obj,  report_common_obj);


 // destruct
 
     // OR_common_att_obs
 		for (int i = 0; i < _nO+1; i++)
		{
			delete OR_common_att_obs[i];
		}
		delete[] OR_common_att_obs;
	
	// OR_common_obj_obs	
		for (int i = 0; i < _nA+1; i++)
		{
			delete OR_common_obj_obs[i];
		}
		delete[] OR_common_obj_obs;
		
		//delete cumul;
		//delete latklas;	
		
	// datarep	
	  for (int i = 0; i < _nO+1; i++)
		{
			for (int ii = 0; ii < _nA+1; ii++)
			{
				delete datarep[i][ii];
			}
			delete datarep[i];
		}
		delete[] datarep;

      // mean_OR_common_att
	  for (int i = 0; i < _nO+1; i++)
		{
			delete mean_OR_common_att[i];
		}
		delete[] mean_OR_common_att;
		
	  // mean_OR_common_obj	
	  for (int i = 0; i < _nA+1; i++)
		{
			delete mean_OR_common_obj[i];
		}
		delete[] mean_OR_common_obj;
 

}




// bootstrap model fixed_object_XZ_YZ
// additive model
void bootstrap_dependency_XZ_YZ_ADD(binary3 dataobs, extended3 ro, extended3 ta, extended1 ga, extended3 OR_common_att_rep, extended3 OR_common_obj_rep, extended2 report_common_att,extended2 report_common_obj)
{
	// bootstrap vars
	extended2 OR_common_att_obs,OR_common_obj_obs;
	extended2 mean_OR_common_att;
	extended2 mean_OR_common_obj;
	binary3 datarep;
	int N_attpair, N_objpair;

	N_attpair = (_nA*(_nA-1)) / 2; 
	N_objpair = (_nO*(_nO-1)) / 2;

	// construct OR_common_att_obs
	OR_common_att_obs = new long double * [_nO+1]; 
	for (int i = 0; i < _nO+1; i++)
	{
		OR_common_att_obs[i] = new long double [N_attpair+1];
	}
	// construct OR_common_obj_obs
	OR_common_obj_obs = new long double * [_nA+1]; 
	for (int i = 0; i < _nA+1; i++)
	{
		OR_common_obj_obs[i] = new long double [N_objpair+1];
	}


	// construct datarep
	datarep = new bool ** [_nO+1];
	for (int i = 0; i < _nO+1; i++)
	{
		datarep[i] = new bool * [_nA+1];
		for (int ii = 0; ii < _nA+1; ii++)
		{
			datarep[i][ii] = new bool [_nR+1];
		}
	}

	//construct mean_OR_common_att & mean_OR_common_obj
	mean_OR_common_att = new long double * [_nO+1] ;
	for (int i = 0; i < _nO+1; i++)
	{
		mean_OR_common_att[i] = new long double [N_attpair+1];
	}
	mean_OR_common_obj = new long double * [_nA+1] ;
	for (int i = 0; i < _nA+1; i++)
	{
		mean_OR_common_obj[i] = new long double [N_objpair+1];
	}



	//  --------------------------------------  STAP 1  ----------------------------------------------

	compute_OR_common_attribute_obs(dataobs, OR_common_att_obs);
	compute_OR_common_object_obs(dataobs, OR_common_obj_obs);


	//  --------------------------------------  STAP 2  ----------------------------------------------

	
	for (int iter = 1; iter <= _Nbootstrap; iter++)
	{ 
		//if (iter mod 500)=0 then writeln(iter:5);
		generate_data_fixedobject_XZ_YZ_ADD( ga, ro, ta, datarep); 
		compute_OR_common_attribute_rep( datarep, iter, OR_common_att_rep );
		compute_OR_common_object_rep( datarep, iter, OR_common_obj_rep );
	}
	

	//  --------------------------------------  STAP 3  ----------------------------------------------

	compute_mean_OR(OR_common_att_rep, _Nbootstrap,_nO, N_attpair,mean_OR_common_att);
	compute_mean_OR(OR_common_obj_rep, _Nbootstrap,_nA, N_objpair,mean_OR_common_obj);

	
	create_report_common_att(OR_common_att_obs, mean_OR_common_att,  report_common_att);
	create_report_common_obj(OR_common_obj_obs, mean_OR_common_obj,  report_common_obj);


 // destruct
 
     // OR_common_att_obs
 		for (int i = 0; i < _nO+1; i++)
		{
			delete OR_common_att_obs[i];
		}
		delete[] OR_common_att_obs;
	
	// OR_common_obj_obs	
		for (int i = 0; i < _nA+1; i++)
		{
			delete OR_common_obj_obs[i];
		}
		delete[] OR_common_obj_obs;
		
	
		
	// datarep	
	  for (int i = 0; i < _nO+1; i++)
		{
			for (int ii = 0; ii < _nA+1; ii++)
			{
				delete datarep[i][ii];
			}
			delete datarep[i];
		}
		delete[] datarep;

      // mean_OR_common_att
	  for (int i = 0; i < _nO+1; i++)
		{
			delete mean_OR_common_att[i];
		}
		delete[] mean_OR_common_att;
		
	  // mean_OR_common_obj	
	  for (int i = 0; i < _nA+1; i++)
		{
			delete mean_OR_common_obj[i];
		}
		delete[] mean_OR_common_obj;
 

}



void compute_OR_common_attribute_obs(binary3 dataB, extended2 OR_common_att_obs) 
{
	unsigned short int countvar = 0;
	for (int j = 1; j <= _nO; j++)                
	{
		countvar = 0;
		for (int u = 1; u <= ( _nA - 1 ); u++)
		{
			for (int v = u + 1; v <= _nA; v++)
			{
				countvar = countvar + 1;
				OR_common_att_obs[j][ countvar] = or_att_pair( j, u, v, dataB ); 
			}
		}
	}
}

void compute_OR_common_object_obs(binary3 dataB, extended2 OR_common_obj_obs) 
{
	unsigned short int countvar = 0;
	for (int k = 1; k <= _nA; k++)
	{
		countvar = 0;
		for (int u = 1; u <= ( _nO - 1 ); u++)
		{
			for (int v = u + 1; v <= _nO; v++)
			{
				countvar = countvar + 1;
				OR_common_obj_obs[k][ countvar] = or_obj_pair( k, u, v, dataB );
			}
		}
	}
}

long double or_att_pair(unsigned short int j, unsigned short int u, unsigned short int v, binary3 data2)
{
	long double result = 0.0;
	long double som00 = 0.0, som01 = 0.0, som10 = 0.0, som11 = 0.0;
	som00 = 0;
	som01 = 0;
	som10 = 0;
	som11 = 0;
	for (int i = 1; i <= _nR; i++)
	{
		if ( ( data2[j][u][i]== 0 ) && ( data2[j][v][i] == 0 ) )
			som00 = som00 + 1; 
		if ( ( data2[j][u][i] == 0 ) && ( data2[j][v][i] == 1 ) )
			som01 = som01 + 1;
		if ( ( data2[j][u][i] == 1 ) && ( data2[j][v][i] == 0 ) )
			som10 = som10 + 1;
		if ( ( data2[j][u][i] == 1 ) && ( data2[j][v][i]== 1 ) )
			som11 = som11 + 1;
	}

	result = logl(som11+0.5) + logl( som00 + 0.5 ) - logl( som10 + 0.5 ) - logl( som01 + 0.5 ); // adapted log long double
	return result;

}


long double or_obj_pair(unsigned short int k,unsigned short int u, unsigned short int v, binary3 data2)
{
	long double result = 0.0;
	long double som00 = 0.0, som01 = 0.0, som10 = 0.0, som11 = 0.0;
	som00 = 0;
	som01 = 0;
	som10 = 0;
	som11 = 0;
	for (int i = 1; i <= _nR; i++)
	{
		if ( ( data2[u][k][i] == 0 ) && ( data2[v][k][i]== 0 ) )
			som00 = som00 + 1;
		if ( ( data2[u][k][i] == 0 ) && ( data2[v][k][i] == 1 ) )
			som01 = som01 + 1;
		if ( ( data2[u][k][i] == 1 ) && ( data2[v][k][i] == 0 ) )
			som10 = som10 + 1;
		if ( ( data2[u][k][i] == 1 ) && ( data2[v][k][i] == 1 ) )
			som11 = som11 + 1;
	}
	result = logl( som11 + 0.5 ) + logl( som00 + 0.5 ) - logl( som10 + 0.5 ) - logl( som01 + 0.5 );
	return result;
}

// generate data fixedobject_XZ_Y
// disjunctive model
void generate_data_fixedobject_XZ_Y_DC(extended1 ga, extended3 ro , extended2 ta, binary3 datarep) 
{
	bool prod;
	bool * latx;
	bool * laty;
	
	unsigned short int * latklas;
	extended1 cumul;
	
	// construct latx & laty
	latx = new bool [_nF+1];
	laty = new bool [_nF+1];
	
	latklas= new unsigned short int [_nR+1];
	cumul=new long double [_nT+1];

	for (int i = 1; i <= _nR; i++)
	{
		latklas[i] = trek_latente_klasse( _nT, ga, cumul);
	}

	for (int j = 1; j <= _nO; j++)
	{
		for (int i = 1; i <= _nR; i++)
		{
			for (int f = 1; f <= _nF; f++)
			{
				if ( runif( seed1, seed2, seed3 ) <= ro[j][ f][ latklas[i]] )
					latx[f] = 1;
				else
					latx[f] = 0;    
			}
			for (int k = 1; k <= _nA; k++)
			{
				for (int f = 1; f <= _nF; f++)
				{
					if ( runif( seed1, seed2, seed3 ) <= ta[k][ f] )
						laty[f] = 1;
					else
						laty[f] = 0; 
				}
				prod = 1;
				for (int f = 1; f <= _nF; f++)
					prod = prod && ( 1 - (latx[f] && laty[f]));
				datarep[j][k][i] = 1 - prod;
			}
		}
	}

	// destruct latx & laty
	delete[] latx;
	delete[] laty;
	
	delete[] cumul;
	delete[] latklas;

}



// generate data fixedobject_X_YZ
// disjunctive model
void generate_data_fixedobject_X_YZ_DC(extended1 ga, extended2 ro , extended3 ta, binary3 datarep) 
{
	bool prod;
	bool * latx;
	bool * laty;
	
	unsigned short int * latklas;
	extended1 cumul;
	
	// construct latx & laty
	latx = new bool [_nF+1];
	laty = new bool [_nF+1];
	
	latklas= new unsigned short int [_nR+1];
	cumul=new long double [_nT+1];

	for (int i = 1; i <= _nR; i++)
	{
		latklas[i] = trek_latente_klasse( _nT, ga, cumul);
	}

	for (int j = 1; j <= _nO; j++)
	{
		for (int i = 1; i <= _nR; i++)
		{
			for (int f = 1; f <= _nF; f++)
			{
				if ( runif( seed1, seed2, seed3 ) <= ro[j][f])
					latx[f] = 1;
				else
					latx[f] = 0;    
			}
			for (int k = 1; k <= _nA; k++)
			{
				for (int f = 1; f <= _nF; f++)
				{
					if ( runif( seed1, seed2, seed3 ) <= ta[k][f][latklas[i]])
						laty[f] = 1;
					else
						laty[f] = 0; 
				}
				prod = 1;
				for (int f = 1; f <= _nF; f++)
					prod = prod && ( 1 - (latx[f] && laty[f]));
				datarep[j][k][i] = 1 - prod;
			}
		}
	}


	// destruct latx & laty
	delete[] latx;
	delete[] laty;
	
	delete[] cumul;
	delete[] latklas;

}




// generate data fixedobject_XZ_YZ
// disjunctive model
void generate_data_fixedobject_XZ_YZ_DC(extended1 ga, extended3 ro , extended3 ta, binary3 datarep) 
{
	bool prod;
	bool * latx;
	bool * laty;
	
	unsigned short int * latklas;
	extended1 cumul;
	
	// construct latx & laty
	latx = new bool [_nF+1];
	laty = new bool [_nF+1];
	
	latklas= new unsigned short int [_nR+1];
	cumul=new long double [_nT+1];

	for (int i = 1; i <= _nR; i++)
	{
		latklas[i] = trek_latente_klasse( _nT, ga, cumul);
	}

	for (int j = 1; j <= _nO; j++)
	{
		for (int i = 1; i <= _nR; i++)
		{
			for (int f = 1; f <= _nF; f++)
			{
				if ( runif( seed1, seed2, seed3 ) <= ro[j][f][latklas[i]])
					latx[f] = 1;
				else
					latx[f] = 0;    
			}
			for (int k = 1; k <= _nA; k++)
			{
				for (int f = 1; f <= _nF; f++)
				{
					if ( runif( seed1, seed2, seed3 ) <= ta[k][f][latklas[i]])
						laty[f] = 1;
					else
						laty[f] = 0; 
				}
				prod = 1;
				for (int f = 1; f <= _nF; f++)
					prod = prod && ( 1 - (latx[f] && laty[f]));
				datarep[j][k][i] = 1 - prod;
			}
		}
	}

	// destruct latx & laty
	delete[] latx;
	delete[] laty;
	
	delete[] cumul;
	delete[] latklas;

}


// generate data fixedobject_XZ_Y
// additive model
void generate_data_fixedobject_XZ_Y_ADD(extended1 ga, extended3 ro , extended2 ta, binary3 datarep) 
{
	long double sum=0.0;
	bool * latx;
	bool * laty;
	
	unsigned short int * latklas;
	extended1 cumul;
	
	// construct latx & laty
	latx = new bool [_nF+1];
	laty = new bool [_nF+1];
	
	latklas= new unsigned short int [_nR+1];
	cumul=new long double [_nT+1];

	for (int i = 1; i <= _nR; i++)
	{
		latklas[i] = trek_latente_klasse( _nT, ga, cumul);
	}

	for (int j = 1; j <= _nO; j++)
	{
		for (int i = 1; i <= _nR; i++)
		{
			for (int f = 1; f <= _nF; f++)
			{
				if (runif( seed1, seed2, seed3) <= ro[j][f][latklas[i]])
					latx[f] = 1;
				else
					latx[f] = 0;    
			}
			for (int k = 1; k <= _nA; k++)
			{
				for (int f = 1; f <= _nF; f++)
				{
					if ( runif( seed1, seed2, seed3 ) <= ta[k][f] )
						laty[f] = 1;
					else
						laty[f] = 0; 
				}
				sum = 0;
				for (int f = 1; f <= _nF; f++) {
					sum = sum + (latx[f] * laty[f]);
				}
				sum=sum/_nF;
				if (runif(seed1, seed2, seed3) <= sum) datarep[j][k][i] = 1;
				else datarep[j][k][i] = 0;
			}
		}
	}

	// destruct latx & laty
	delete[] latx;
	delete[] laty;
	
	delete[] cumul;
	delete[] latklas;

}


// generate data fixedobject_X_YZ
// additive model
void generate_data_fixedobject_X_YZ_ADD(extended1 ga, extended2 ro , extended3 ta, binary3 datarep) 
{
	long double sum=0.0;
	bool * latx;
	bool * laty;
	
	unsigned short int * latklas;
	extended1 cumul;
	
	// construct latx & laty
	latx = new bool [_nF+1];
	laty = new bool [_nF+1];
	
	latklas= new unsigned short int [_nR+1];
	cumul=new long double [_nT+1];

	for (int i = 1; i <= _nR; i++)
	{
		latklas[i] = trek_latente_klasse( _nT, ga, cumul);
	}

	for (int j = 1; j <= _nO; j++)
	{
		for (int i = 1; i <= _nR; i++)
		{
			for (int f = 1; f <= _nF; f++)
			{
				if ( runif( seed1, seed2, seed3 ) <= ro[j][f])
					latx[f] = 1;
				else
					latx[f] = 0;    
			}
			for (int k = 1; k <= _nA; k++)
			{
				for (int f = 1; f <= _nF; f++)
				{
					if ( runif( seed1, seed2, seed3 ) <= ta[k][f][latklas[i]])
						laty[f] = 1;
					else
						laty[f] = 0; 
				}
                sum = 0;
				for (int f = 1; f <= _nF; f++) {
					sum = sum + (latx[f] * laty[f]);
				}
				sum=sum/_nF;
				if (runif(seed1, seed2, seed3) <= sum) datarep[j][k][i] = 1;
				else datarep[j][k][i] = 0;
			}
		}
	}


	// destruct latx & laty
	delete[] latx;
	delete[] laty;
	
	delete[] cumul;
	delete[] latklas;

}




// generate data fixedobject_XZ_YZ
//additive model
void generate_data_fixedobject_XZ_YZ_ADD(extended1 ga, extended3 ro , extended3 ta, binary3 datarep) 
{
	long double sum=0.0;
	bool * latx;
	bool * laty;
	
	unsigned short int * latklas;
	extended1 cumul;
	
	// construct latx & laty
	latx = new bool [_nF+1];
	laty = new bool [_nF+1];
	
	latklas= new unsigned short int [_nR+1];
	cumul=new long double [_nT+1];

	for (int i = 1; i <= _nR; i++)
	{
		latklas[i] = trek_latente_klasse( _nT, ga, cumul);
	}

	for (int j = 1; j <= _nO; j++)
	{
		for (int i = 1; i <= _nR; i++)
		{
			for (int f = 1; f <= _nF; f++)
			{
				if ( runif( seed1, seed2, seed3 ) <= ro[j][f][latklas[i]])
					latx[f] = 1;
				else
					latx[f] = 0;    
			}
			for (int k = 1; k <= _nA; k++)
			{
				for (int f = 1; f <= _nF; f++)
				{
					if ( runif( seed1, seed2, seed3 ) <= ta[k][f][latklas[i]])
						laty[f] = 1;
					else
						laty[f] = 0; 
				}
                sum = 0;
				for (int f = 1; f <= _nF; f++) {
					sum = sum + (latx[f] * laty[f]);
				}
				sum=sum/_nF;
				if (runif(seed1, seed2, seed3) <= sum) datarep[j][k][i] = 1;
				else datarep[j][k][i] = 0;
			}
		}
	}

	// destruct latx & laty
	delete[] latx;
	delete[] laty;
	
	delete[] cumul;
	delete[] latklas;

}





unsigned short int trek_latente_klasse( unsigned short int lengte, extended1 vec, extended1 cumul)
{

	unsigned short int teller = 0;
	double getal = 0.0;
	cumul[0] = 0;
	for (int i = 1; i <= lengte; i++)
		cumul[i] = vec[i];
	for (int i = 1; i <= lengte; i++)
		cumul[i] = cumul[i - 1] + cumul[i];
	cumul[lengte] = 1;
	teller = 0;
	getal = runif( seed1, seed2, seed3 );
	do
	{
		teller = teller + 1;
	}
	while ( ! ( ( cumul[teller - 1] < getal ) && ( cumul[teller] >= getal ) ) );
	return teller;
}

void compute_OR_common_attribute_rep(binary3 datarep, unsigned short int row, extended3 OR_common_att_rep)
{
	unsigned short int countvar = 0;
	for (int j = 1; j <= _nO; j++)
	{
		countvar = 0;
		for (int u = 1; u <= ( _nA - 1 ); u++)
		{
			for (int v = u + 1; v <= _nA; v++)
			{
				countvar = countvar + 1;
				OR_common_att_rep[row][ j][ countvar] = or_att_pair( j, u, v, datarep );
			}
		}
	}
}

void compute_OR_common_object_rep(binary3 datarep, unsigned short int row, extended3 OR_common_obj_rep)
{
	unsigned short int countvar = 0;
	for (int k = 1; k <= _nA; k++)
	{
		countvar = 0;
		for (int u = 1; u <= ( _nO - 1 ); u++)
		{
			for (int v = u + 1; v <= _nO; v++)
			{
				countvar = countvar + 1;
				OR_common_obj_rep[row][ k][ countvar] = or_obj_pair( k, u, v, datarep );
			}
		}
	}
}


void compute_mean_OR( const extended3 datarep, unsigned short int Nbootstrap, unsigned short int Nelement, unsigned short int Npair, extended2 mean_OR)
{
	long double som = 0.0;
	for (int j = 1; j <= Nelement; j++)
	{
		for (int k = 1; k <= Npair; k++)
		{
			som = 0;
			for (int i = 1; i <= Nbootstrap; i++) 
			som = som + datarep[i][j][k];
			som = double( som ) / Nbootstrap;
			mean_OR[j][ k] = som;
		}
	}
}

void create_report_common_att(extended2 OR_common_att_obs, extended2 mean_OR_common_att, extended2 report_common_att)
{
	unsigned short int countvar1 = 0;
	unsigned short int countvar2 = 0;
	for (int j = 1; j <= _nO; j++)
	{
		countvar2 = 0;
		for (int k = 1; k <= _nA - 1; k++)
		{
			for (int l = k + 1; l <= _nA; l++)
			{
				countvar1 = countvar1 + 1;
				countvar2 = countvar2 + 1;
				report_common_att[countvar1][ 1] = j;
				report_common_att[countvar1][ 2] = k;
				report_common_att[countvar1][ 3] = l;
				report_common_att[countvar1][ 4] = OR_common_att_obs[j][ countvar2];
				report_common_att[countvar1][ 5] = mean_OR_common_att[j][ countvar2];
			}
		}
	}
}

void create_report_common_obj(extended2 OR_common_obj_obs, extended2 mean_OR_common_obj, extended2 report_common_obj)
{
	unsigned short int countvar1 = 0;
	unsigned short int countvar2 = 0;
	for (int k = 1; k <= _nA; k++)
	{
		countvar2 = 0;
		for (int j = 1; j <= _nO - 1; j++)
		{
			for (int l = j + 1; l <= _nO; l++)
			{
				countvar1 = countvar1 + 1;
				countvar2 = countvar2 + 1;
				report_common_obj[countvar1][ 1] = k;
				report_common_obj[countvar1][ 2] = j;
				report_common_obj[countvar1][ 3] = l;
				report_common_obj[countvar1][ 4] = OR_common_obj_obs[k][ countvar2];
				report_common_obj[countvar1][ 5] = mean_OR_common_obj[k][ countvar2];
			}
		}
	}
}
























