/*##############################################################################################*/
// source:       rInterface.c
//------------------------------------------------------------------------------------------------
// type:         Function
//------------------------------------------------------------------------------------------------
// title:        C interface for R functions
//------------------------------------------------------------------------------------------------
// description:  Interface to make connections between R and the called C functions
//                |                                      |                                       |
//                |                R                     |                   C                   |
//                |                                      |                                       |
//
//                c_blas_dgemm                          -> _r_cblas_dgemm    -> gsl_blas_dgemm()
//
//                c_int_det         -> c_call_integrate -> _r_int_det        -> int_det()
//              * c_int_det_hemodyn -> c_call_integrate -> _r_int_hemodyn    -> int_det_hemodyn()
//                c_int_dfdp        -> c_call_intdiff   -> _r_int_dfdp       -> int_dfdp()
//
//                c_bireduce        -> c_call_reduce    -> _r_bireduce       -> bireduce()
//                                                                           -> reduce_out()
//                c_soreduce        -> c_call_reduce    -> _r_soreduce       -> soreduce()
//                                                                           -> reduce_out()
//
//                c_fx_fmri         -> c_call_fx        -> _r_fx_fmri        -> fx_fmri()
//                c_dfdx            -> c_call_fx        -> _r_dfdx           -> dfdx()
//                c_dfdu            -> c_call_fx        -> _r_dfdu           -> dfdu()
//                c_dfdxdu          -> c_call_fx        -> _r_dfdxdu         -> dfdxdu()
//                c_dfdxdx          -> c_call_fx        -> _r_dfdxdx         -> dfdxdx()
//
//                c_gx_fmri         -> c_call_gx        -> _r_gx_fmri        -> gx_fmri()
//                c_dgdx            -> c_call_gx        -> _r_dgdx           -> dgdx()
//                c_dgdxdx          -> c_call_gx        -> _r_dgdxdx         -> dgdxdx()
//
//    * obsolete function
//   ** not tested
//  *** missing interface
//------------------------------------------------------------------------------------------------
// based on:     ---
//------------------------------------------------------------------------------------------------
// reference:    ---
//------------------------------------------------------------------------------------------------
// remarks:      jacobians outputs are read column-wise in R by default
//------------------------------------------------------------------------------------------------
// TODO:         ---
/*##############################################################################################*/

//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_multifit_nlin.h>
//#include <gsl/gsl_blas.h>

//#include "arrays.h"
#include "utils.h"
#include "integrate.h"
#include "reduce.h"
#include "fx_fmri.h"
#include "gx_fmri.h"

#include "redcmc_export.h"




// Interface for utility functions
// -----------------------------------------------------------------------------------------------------------------
REDCMC_EXPORT void _r_cblas_dgemm( double* A, double* B, int* d, double* Y )
{

    gsl_matrix_view Ax = gsl_matrix_view_array( A, d[0], d[1] );
    gsl_matrix_view Bx = gsl_matrix_view_array( B, d[2], d[3] );
    gsl_matrix_view Yx = gsl_matrix_view_array( Y, d[0], d[3] );

    gsl_blas_dgemm ( CblasNoTrans, CblasNoTrans, 1.0, &Ax.matrix, &Bx.matrix, 0.0, &Yx.matrix );

}




// Interface for model integration functions
// -----------------------------------------------------------------------------------------------------------------

REDCMC_EXPORT void _r_int_det( double* x, double* timing, double* P, int* M, double* U, double* H, double* Y )
{

    red_par* Model = (red_par*) malloc( sizeof(red_par) );
    unsigned int i;
    int nu = M[0];
    int ny = M[1];
    int nyny = ny*ny;
    int postD;
    int A_offs;

    double* u = (double*) calloc( nu, sizeof(double) );

    int_par* ModelInt = (int_par*) malloc( sizeof(int_par) );


    Model->nu = M[0];							// number of exogenous inputs
    Model->ny = M[1];							// number of response variables
    Model->nonlin = (char) M[2];					// nonlinear DCM switch
    Model->twostate = (char) M[3];					// two-state DCM switch
    Model->symmetry = (char) M[4];					// symmetric connectivity matrix switch
    Model->backward = (char) M[5];					// backward connections
    Model->A_eigen = (char) M[6];					// average connections encodes eigenvalue
	A_offs = (!Model->backward) ? 0 : nyny;
	postD = (!Model->nonlin) ? ((!Model->backward) ? nyny+nu*nyny+nu*ny : 2*nyny+nu*nyny+nu*ny) :
			    ((!Model->backward) ? nyny+nu*nyny+nu*ny+ny*nyny : 2*nyny+nu*nyny+nu*ny+ny*nyny );

    Model->x = x;							// state variables (neuronal+hemodynamic)

    Model->u = u;							// time-point of input function

    Model->A = &P[0];							// linear parameters
    Model->B = &P[nyny+A_offs];						// bilinear parameters
    Model->C = &P[nyny+A_offs+nu*nyny];					// exogenous parameters
    Model->D = (!Model->nonlin) ? NULL : &P[nyny+A_offs+nu*nyny+nu*ny];	// nonlinear parameters
    Model->transit = &P[postD];						// transit time coefficient
    Model->decay = &P[postD+ny];					// signal decay coefficient
    Model->epsilon = P[postD+ny+ny];
    Model->H = H;							// hemodynamic constants


    printf( "decay: " );
    for ( i = 0; i < Model->ny; ++i )
    {
	printf( "%f ", Model->decay[i] );
    }
    printf( "\ntransit: " );
    for ( i = 0; i < Model->ny; ++i )
    {
	printf( "%f ", Model->transit[i] );
    }
    printf( "\nepsilon: %f\n", Model->epsilon );


    ModelInt->Model = Model;						// Model structure
    ModelInt->TR = timing[0];						// TR (session)
    ModelInt->U_ext = U;						// external input
    ModelInt->delays = &timing[4];					// TR delay (regional)
    ModelInt->dt = timing[1];						// microtime intervals
    ModelInt->u = (int)timing[2];					// microtime bins
    ModelInt->ns = (int)timing[3];					// time samples

    //printf( "delays: %f %f %f %f %f\n", ModelInt->delays[0], ModelInt->delays[1], ModelInt->delays[2], ModelInt->delays[3], ModelInt->delays[4]);

    if ( !Model->nonlin )
    {
	int_det( ModelInt, Y );
    }
    else
    {
	int_det_D( ModelInt, Y );
    }

    free( u );
    free( Model );
    free( ModelInt );

}


REDCMC_EXPORT void _r_int_det_hemodyn( double* x, double* timing, double* P, int* M, double* U, double* H, double* Y, double* X )
{

    red_par* Model = (red_par*) malloc( sizeof(red_par) );
    int nu = M[0];
    int ny = M[1];
    int nyny = ny*ny;
    int postD;
    int A_offs;

    double* u = (double*) calloc( nu, sizeof(double) );

    int_par* ModelInt = (int_par*) malloc( sizeof(int_par) );


    Model->nu = M[0];							// number of exogenous inputs
    Model->ny = M[1];							// number of response variables
    Model->nonlin = (char) M[2];					// nonlinear DCM switch
    Model->twostate = (char) M[3];					// two-state DCM switch
    Model->symmetry = (char) M[4];					// symmetric connectivity matrix switch
    Model->backward = (char) M[5];					// backward connections
    Model->A_eigen = (char) M[6];					// average connections encodes eigenvalue
	A_offs = (!Model->backward) ? 0 : nyny;
	postD = (!Model->nonlin) ? ((!Model->backward) ? nyny+nu*nyny+nu*ny : 2*nyny+nu*nyny+nu*ny) :
			    ((!Model->backward) ? nyny+nu*nyny+nu*ny+ny*nyny : 2*nyny+nu*nyny+nu*ny+ny*nyny );

    Model->x = x;							// state variables (neuronal+hemodynamic)

    Model->u = u;							// time-point of input function

    Model->A = &P[0];							// linear parameters
    Model->B = &P[nyny+A_offs];						// bilinear parameters
    Model->C = &P[nyny+A_offs+nu*nyny];					// exogenous parameters
    Model->D = (!Model->nonlin) ? NULL : &P[nyny+A_offs+nu*nyny+nu*ny];	// nonlinear parameters
    Model->transit = &P[postD];						// transit time coefficient
    Model->decay = &P[postD+ny];					// signal decay coefficient
    Model->epsilon = P[postD+ny+ny];
    Model->H = H;							// hemodynamic constants


    ModelInt->Model = Model;						// Model structure
    ModelInt->TR = timing[0];						// TR (session)
    ModelInt->U_ext = U;						// external input
    ModelInt->delays = &timing[4];					// TR delay (regional)
    ModelInt->dt = timing[1];						// microtime intervals
    ModelInt->u = (int)timing[2];					// microtime bins
    ModelInt->ns = (int)timing[3];					// time samples


    int_det_hemodyn( ModelInt, Y, X );

    free( u );
    free( Model );
    free( ModelInt );

}


REDCMC_EXPORT void _r_int_dfdp( double* x, double* timing, double* Pr, int* M, double* U, double* H, double* V, double* Y )
{

    red_par* Model = (red_par*) malloc( sizeof(red_par) );
    int_par* ModelInt = (int_par*) malloc( sizeof(int_par) );
    gsl_matrix_view J;

    int nu = M[0];
    int ny = M[1];
    int nyny = ny*ny;
    int postD;
    int A_offs;

    double* u = (double*) calloc( nu, sizeof(double) );

    Model->nu = M[0];							// number of exogenous inputs
    Model->ny = M[1];							// number of response variables
    Model->nonlin = (char) M[2];					// nonlinear DCM switch
    Model->twostate = (char) M[3];					// two-state DCM switch
    Model->symmetry = (char) M[4];					// symmetric connectivity matrix switch
    Model->backward = (char) M[5];					// backward connections
    Model->A_eigen = (char) M[6];					// average connections encodes eigenvalue
	A_offs = (!Model->backward) ? 0 : nyny;
	postD = (!Model->nonlin) ? ((!Model->backward) ? nyny+nu*nyny+nu*ny : 2*nyny+nu*nyny+nu*ny) :
			    ((!Model->backward) ? nyny+nu*nyny+nu*ny+ny*nyny : 2*nyny+nu*nyny+nu*ny+ny*nyny );

    Model->x = x;							// state variables (neuronal+hemodynamic)

    Model->u = u;							// time-point of input function

    Model->H = H;							// hemodynamic constants

    /*
    Model->A = &P[0];							// linear parameters
    Model->B = &P[nyny+A_offs];						// bilinear parameters
    Model->C = &P[nyny+A_offs+nu*nyny];					// exogenous parameters
    Model->D = (!Model->nonlin) ? NULL : &P[nyny+A_offs+nu*nyny+nu*ny];	// nonlinear parameters
    Model->transit = &P[postD];						// transit time coefficient
    Model->decay = &P[postD+ny];					// signal decay coefficient
    Model->epsilon = P[postD+ny+ny];
    */


    ModelInt->Model = Model;						// Model structure
    ModelInt->TR = timing[0];						// TR (session)
    ModelInt->U_ext = U;						// external input
    ModelInt->delays = &timing[4];					// TR delay (regional)
    ModelInt->dt = timing[1];						// microtime intervals
    ModelInt->u = (int)timing[2];					// microtime bins
    ModelInt->ns = (int)timing[3];					// time samples
    ModelInt->Pr = Pr;							// parameters in reduced space
    ModelInt->V = V;							// transformation matrix
    ModelInt->Vr = M[7];						// V_rows
    ModelInt->Vc = M[8];						// V_columns



    J = gsl_matrix_view_array( Y, ModelInt->ns*ny, ModelInt->Vc );

    int_dfdp( ModelInt, &J.matrix );

    free( u );
    free( Model );
    free( ModelInt );

}



// Interface for nonlinear MIMO reduction functions
// -----------------------------------------------------------------------------------------------------------------

REDCMC_EXPORT void _r_bireduce( double* x, double* H, int* M, double* P, double* M0, double* M1, double* L1, double* L2, int* out )
{

    unsigned int i;

    red_par* Model = (red_par*) malloc( sizeof(red_par) );
    int nu = M[0];
    int ny = M[1];
    int nyny = ny*ny;
    int postD;
    int A_offs;
    int mxsize;
    gsl_matrix_view mxv1, mxv2;

    double* u = (double*) calloc( nu, sizeof(double) );

    Model->nu = M[0];							// number of exogenous inputs
    Model->ny = M[1];							// number of response variables
    Model->nonlin = (char) M[2];					// nonlinear DCM switch
    Model->twostate = (char) M[3];					// two-state DCM switch
    Model->symmetry = (char) M[4];					// symmetric connectivity matrix switch
    Model->backward = (char) M[5];					// backward connections
    Model->A_eigen = (char) M[6];					// average connections encodes eigenvalue
	A_offs = (!Model->backward) ? 0 : nyny;
	postD = (!Model->nonlin) ? ((!Model->backward) ? nyny+nu*nyny+nu*ny : 2*nyny+nu*nyny+nu*ny) :
			    ((!Model->backward) ? nyny+nu*nyny+nu*ny+ny*nyny : 2*nyny+nu*nyny+nu*ny+ny*nyny );
	mxsize = (!Model->twostate) ? 5*ny+1 : 6*ny+1;

    Model->x = x;							// state variables (neuronal+hemodynamic)

    Model->u = u;							// time-point of input function

    Model->A = &P[0];							// linear parameters
    Model->B = &P[nyny+A_offs];						// bilinear parameters
    Model->C = &P[nyny+A_offs+nu*nyny];					// exogenous parameters
    Model->D = (!Model->nonlin) ? NULL : &P[nyny+A_offs+nu*nyny+nu*ny];	// nonlinear parameters
    Model->transit = &P[postD];						// transit time coefficient
    Model->decay = &P[postD+ny];					// signal decay coefficient
    Model->epsilon = P[postD+ny+ny];
    Model->H = H;							// hemodynamic constants


    //printf( "decay: %f %f %f %f\n", Model->decay[0], Model->decay[1], Model->decay[2], Model->decay[3] );
    //printf( "transit: %f %f %f %f\n", Model->transit[0], Model->transit[1], Model->transit[2], Model->transit[3] );
    //printf( "epsilon: %f\n", Model->epsilon );


    //printf( "bireduce: model derivatives\n" );
    bireduce( Model, M0, M1, L1, L2 );
    if ( *out == 1 )
    {
	//printf( "bireduce: output derivatives\n" );
	reduce_out( Model, M0, M1, L1, L2 );
    }


    mxv1 = gsl_matrix_view_array( M0, mxsize, mxsize );
    gsl_matrix_transpose( &mxv1.matrix );
    mxv2 = gsl_matrix_view_array( L1, mxsize, mxsize );
    gsl_matrix_transpose( &mxv2.matrix );
    for ( i = 0; i < nu; ++i )
    {
	mxv1 = gsl_matrix_view_array( &M1[i*mxsize*mxsize], mxsize, mxsize );
	gsl_matrix_transpose( &mxv1.matrix );
    }
    for ( i = 0; i < ny; ++i )
    {
	mxv2 = gsl_matrix_view_array( &L2[i*mxsize*mxsize], mxsize, mxsize );
	gsl_matrix_transpose( &mxv2.matrix );
    }

    free( u );
    free( Model );

}


REDCMC_EXPORT void _r_soreduce( double* x, double* H, int* M, double* P, double* M0, double* M1, double* M2, double* L1, double* L2, int* out )
{

    unsigned int i;

    red_par* Model = (red_par*) malloc( sizeof(red_par) );
    int nu = M[0];
    int ny = M[1];
    int nyny = ny*ny;
    int postD;
    int A_offs;
    int mxsize;
    gsl_matrix_view mxv1, mxv2;

    double* u = (double*) calloc( nu, sizeof(double) );

    Model->nu = M[0];							// number of exogenous inputs
    Model->ny = M[1];							// number of response variables
    Model->nonlin = (char) M[2];					// nonlinear DCM switch
    Model->twostate = (char) M[3];					// two-state DCM switch
    Model->symmetry = (char) M[4];					// symmetric connectivity matrix switch
    Model->backward = (char) M[5];					// backward connections
    Model->A_eigen = (char) M[6];					// average connections encodes eigenvalue
	A_offs = (!Model->backward) ? 0 : nyny;
	postD = (!Model->nonlin) ? ((!Model->backward) ? nyny+nu*nyny+nu*ny : 2*nyny+nu*nyny+nu*ny) :
			    ((!Model->backward) ? nyny+nu*nyny+nu*ny+ny*nyny : 2*nyny+nu*nyny+nu*ny+ny*nyny );
	mxsize = (!Model->twostate) ? 5*ny+1 : 6*ny+1;

    Model->x = x;							// state variables (neuronal+hemodynamic)

    Model->u = u;							// time-point of input function

    Model->A = &P[0];							// linear parameters
    Model->B = &P[nyny+A_offs];						// bilinear parameters
    Model->C = &P[nyny+A_offs+nu*nyny];					// exogenous parameters
    Model->D = (!Model->nonlin) ? NULL : &P[nyny+A_offs+nu*nyny+nu*ny];	// nonlinear parameters
    Model->transit = &P[postD];						// transit time coefficient
    Model->decay = &P[postD+ny];					// signal decay coefficient
    Model->epsilon = P[postD+ny+ny];
    Model->H = H;							// hemodynamic constants


    //printf( "decay: %f %f %f %f\n", Model->decay[0], Model->decay[1], Model->decay[2], Model->decay[3] );
    //printf( "transit: %f %f %f %f\n", Model->transit[0], Model->transit[1], Model->transit[2], Model->transit[3] );
    //printf( "epsilon: %f\n", Model->epsilon );


    //printf( "bireduce: model derivatives\n" );
    soreduce( Model, M0, M1, M2, L1, L2 );
    if ( *out == 1 )
    {
	//printf( "bireduce: output derivatives\n" );
	reduce_out( Model, M0, M1, L1, L2 );
    }


    mxv1 = gsl_matrix_view_array( M0, mxsize, mxsize );
    gsl_matrix_transpose( &mxv1.matrix );
    mxv2 = gsl_matrix_view_array( L1, mxsize, mxsize );
    gsl_matrix_transpose( &mxv2.matrix );
    for ( i = 0; i < nu; ++i )
    {
	mxv1 = gsl_matrix_view_array( &M1[i*mxsize*mxsize], mxsize, mxsize );
	gsl_matrix_transpose( &mxv1.matrix );
    }
    for ( i = 0; i < ny; ++i )
    {
	mxv2 = gsl_matrix_view_array( &L2[i*mxsize*mxsize], mxsize, mxsize );
	gsl_matrix_transpose( &mxv2.matrix );
    }
    for ( i = 0; i < mxsize - 1; ++i )
    {
	mxv1 = gsl_matrix_view_array( &M2[i*mxsize*mxsize], mxsize, mxsize );
	gsl_matrix_transpose( &mxv2.matrix );
    }

    free( u );
    free( Model );

}




// Interface for fx_fmri, gx_fmri functions
// -----------------------------------------------------------------------------------------------------------------

REDCMC_EXPORT void _r_fx_fmri( double* x, double* u, double* P, int* M, double* H, double* f_out )
{

    fx_par* Model = (fx_par*) malloc( sizeof(fx_par) );
    int nu = M[0];
    int ny = M[1];
    int nyny = ny*ny;
    int postD;
    int A_offs;
    int i;

    Model->nu = M[0];							// number of exogenous inputs
    Model->ny = M[1];							// number of response variables
    Model->nonlin = (char) M[2];					// nonlinear DCM switch
    Model->twostate = (char) M[3];					// two-state DCM switch
    Model->symmetry = (char) M[4];					// symmetric connectivity matrix switch
    Model->backward = (char) M[5];					// backward connections
    Model->A_eigen = (char) M[6];					// average connections encodes eigenvalue
	A_offs = (!Model->backward) ? 0 : nyny;
	postD = (!Model->nonlin) ? ((!Model->backward) ? nyny+nu*nyny+nu*ny : 2*nyny+nu*nyny+nu*ny) :
			    ((!Model->backward) ? nyny+nu*nyny+nu*ny+ny*nyny : 2*nyny+nu*nyny+nu*ny+ny*nyny );

    Model->x = x;							// state variables (neuronal+hemodynamic)

    Model->u = u;							// time-point of input function

    Model->A = &P[0];							// linear parameters
    Model->B = &P[nyny+A_offs];						// bilinear parameters
    Model->C = &P[nyny+A_offs+nu*nyny];					// exogenous parameters
    Model->D = (!Model->nonlin) ? NULL : &P[nyny+A_offs+nu*nyny+nu*ny];	// nonlinear parameters
    Model->transit = &P[postD];						// transit time coefficient
    Model->decay = &P[postD+ny];					// signal decay coefficient
    Model->epsilon = P[postD+ny+ny];
    Model->H = H;							// hemodynamic constants


    //printf( "decay: %f %f %f %f\n", Model->decay[0], Model->decay[1], Model->decay[2], Model->decay[3] );
    //printf( "transit: %f %f %f %f\n", Model->transit[0], Model->transit[1], Model->transit[2], Model->transit[3] );
    //printf( "epsilon: %f\n", Model->epsilon );


    fx_fmri( Model, f_out );

    free( Model );

}


REDCMC_EXPORT void _r_dfdx( double* x, double* u, double* P, int* M, double* H, double* f_out )
{

    fx_par* Model = (fx_par*) malloc( sizeof(fx_par) );
    gsl_matrix_view J;

    int nu = M[0];
    int ny = M[1];
    int nyny = ny*ny;
    int postD;
    int A_offs;
    int states;

    Model->nu = M[0];							// number of exogenous inputs
    Model->ny = M[1];							// number of response variables
    Model->nonlin = (char) M[2];					// nonlinear DCM switch
    Model->twostate = (char) M[3];					// two-state DCM switch
    Model->symmetry = (char) M[4];					// symmetric connectivity matrix switch
    Model->backward = (char) M[5];					// backward connections
    Model->A_eigen = (char) M[6];					// average connections encodes eigenvalue
	A_offs = (!Model->backward) ? 0 : nyny;
	postD = (!Model->nonlin) ? ((!Model->backward) ? nyny+nu*nyny+nu*ny : 2*nyny+nu*nyny+nu*ny) :
			    ((!Model->backward) ? nyny+nu*nyny+nu*ny+ny*nyny : 2*nyny+nu*nyny+nu*ny+ny*nyny );
	states = (!Model->twostate) ? 5 : 6;

    Model->x = x;							// state variables (neuronal+hemodynamic)

    Model->u = u;							// time-point of input function

    Model->A = &P[0];							// linear parameters
    Model->B = &P[nyny+A_offs];						// bilinear parameters
    Model->C = &P[nyny+A_offs+nu*nyny];					// exogenous parameters
    Model->D = (!Model->nonlin) ? NULL : &P[nyny+A_offs+nu*nyny+nu*ny];	// nonlinear parameters
    Model->transit = &P[postD];						// transit time coefficient
    Model->decay = &P[postD+ny];					// signal decay coefficient
    Model->epsilon = P[postD+ny+ny];
    Model->H = H;							// hemodynamic constants


    //printf( "decay: %f %f %f %f\n", Model->decay[0], Model->decay[1], Model->decay[2], Model->decay[3] );
    //printf( "transit: %f %f %f %f\n", Model->transit[0], Model->transit[1], Model->transit[2], Model->transit[3] );
    //printf( "epsilon: %f\n", Model->epsilon );


    J = gsl_matrix_view_array( f_out, Model->ny*states, Model->ny*states );

    dfdx( Model, &J.matrix );

    free( Model );

}


REDCMC_EXPORT void _r_dfdu( double* x, double* u, double* P, int* M, double* H, double* f_out )
{

    fx_par* Model = (fx_par*) malloc( sizeof(fx_par) );
    gsl_matrix_view J;

    int nu = M[0];
    int ny = M[1];
    int nyny = ny*ny;
    int postD;
    int A_offs;
    int states;

    Model->nu = M[0];							// number of exogenous inputs
    Model->ny = M[1];							// number of response variables
    Model->nonlin = (char) M[2];					// nonlinear DCM switch
    Model->twostate = (char) M[3];					// two-state DCM switch
    Model->symmetry = (char) M[4];					// symmetric connectivity matrix switch
    Model->backward = (char) M[5];					// backward connections
    Model->A_eigen = (char) M[6];					// average connections encodes eigenvalue
	A_offs = (!Model->backward) ? 0 : nyny;
	postD = (!Model->nonlin) ? ((!Model->backward) ? nyny+nu*nyny+nu*ny : 2*nyny+nu*nyny+nu*ny) :
			    ((!Model->backward) ? nyny+nu*nyny+nu*ny+ny*nyny : 2*nyny+nu*nyny+nu*ny+ny*nyny );
	states = (!Model->twostate) ? 5 : 6;

    Model->x = x;							// state variables (neuronal+hemodynamic)

    Model->u = u;							// time-point of input function

    Model->A = &P[0];							// linear parameters
    Model->B = &P[nyny+A_offs];						// bilinear parameters
    Model->C = &P[nyny+A_offs+nu*nyny];					// exogenous parameters
    Model->D = (!Model->nonlin) ? NULL : &P[nyny+A_offs+nu*nyny+nu*ny];	// nonlinear parameters
    Model->transit = &P[postD];						// transit time coefficient
    Model->decay = &P[postD+ny];					// signal decay coefficient
    Model->epsilon = P[postD+ny+ny];
    Model->H = H;							// hemodynamic constants


    //printf( "decay: %f %f %f %f\n", Model->decay[0], Model->decay[1], Model->decay[2], Model->decay[3] );
    //printf( "transit: %f %f %f %f\n", Model->transit[0], Model->transit[1], Model->transit[2], Model->transit[3] );
    //printf( "epsilon: %f\n", Model->epsilon );


    J = gsl_matrix_view_array( f_out, Model->ny*states, Model->nu );

    dfdu( Model, &J.matrix );

    free( Model );

}


REDCMC_EXPORT void _r_dfdxdu( double* x, double* u, double* P, int* M, double* H, double* f_out )
{

    fx_par* Model = (fx_par*) malloc( sizeof(fx_par) );
    gsl_matrix_view J;

    int nu = M[0];
    int ny = M[1];
    int nyny = ny*ny;
    int postD;
    int A_offs;
    int states;

    Model->nu = M[0];							// number of exogenous inputs
    Model->ny = M[1];							// number of response variables
    Model->nonlin = (char) M[2];					// nonlinear DCM switch
    Model->twostate = (char) M[3];					// two-state DCM switch
    Model->symmetry = (char) M[4];					// symmetric connectivity matrix switch
    Model->backward = (char) M[5];					// backward connections
    Model->A_eigen = (char) M[6];					// average connections encodes eigenvalue
	A_offs = (!Model->backward) ? 0 : nyny;
	postD = (!Model->nonlin) ? ((!Model->backward) ? nyny+nu*nyny+nu*ny : 2*nyny+nu*nyny+nu*ny) :
			    ((!Model->backward) ? nyny+nu*nyny+nu*ny+ny*nyny : 2*nyny+nu*nyny+nu*ny+ny*nyny );
	states = (!Model->twostate)? 5 : 6;

    Model->x = x;							// state variables (neuronal+hemodynamic)

    Model->u = u;							// time-point of input function

    Model->A = &P[0];							// linear parameters
    Model->B = &P[nyny+A_offs];						// bilinear parameters
    Model->C = &P[nyny+A_offs+nu*nyny];					// exogenous parameters
    Model->D = (!Model->nonlin) ? NULL : &P[nyny+A_offs+nu*nyny+nu*ny];	// nonlinear parameters
    Model->transit = &P[postD];						// transit time coefficient
    Model->decay = &P[postD+ny];					// signal decay coefficient
    Model->epsilon = P[postD+ny+ny];
    Model->H = H;							// hemodynamic constants


    //printf( "decay: %f %f %f %f\n", Model->decay[0], Model->decay[1], Model->decay[2], Model->decay[3] );
    //printf( "transit: %f %f %f %f\n", Model->transit[0], Model->transit[1], Model->transit[2], Model->transit[3] );
    //printf( "epsilon: %f\n", Model->epsilon );


    J = gsl_matrix_view_array( f_out, Model->ny*states * Model->ny*states, Model->nu );

    dfdxdu( Model, &J.matrix );

    free( Model );

}


REDCMC_EXPORT void _r_dfdxdx( double* x, double* u, double* P, int* M, double* H, double* f_out )
{

    fx_par* Model = (fx_par*) malloc( sizeof(fx_par) );
    gsl_matrix_view J;

    int nu = M[0];
    int ny = M[1];
    int nyny = ny*ny;
    int postD;
    int A_offs;
    int states;

    Model->nu = M[0];							// number of exogenous inputs
    Model->ny = M[1];							// number of response variables
    Model->nonlin = (char) M[2];					// nonlinear DCM switch
    Model->twostate = (char) M[3];					// two-state DCM switch
    Model->symmetry = (char) M[4];					// symmetric connectivity matrix switch
    Model->backward = (char) M[5];					// backward connections
    Model->A_eigen = (char) M[6];					// average connections encodes eigenvalue
	A_offs = (!Model->backward) ? 0 : nyny;
	postD = (!Model->nonlin) ? ((!Model->backward) ? nyny+nu*nyny+nu*ny : 2*nyny+nu*nyny+nu*ny) :
			    ((!Model->backward) ? nyny+nu*nyny+nu*ny+ny*nyny : 2*nyny+nu*nyny+nu*ny+ny*nyny );
	states = (!Model->twostate)? 5 : 6;

    Model->x = x;							// state variables (neuronal+hemodynamic)

    Model->u = u;							// time-point of input function

    Model->A = &P[0];							// linear parameters
    Model->B = &P[nyny+A_offs];						// bilinear parameters
    Model->C = &P[nyny+A_offs+nu*nyny];					// exogenous parameters
    Model->D = (!Model->nonlin) ? NULL : &P[nyny+A_offs+nu*nyny+nu*ny];	// nonlinear parameters
    Model->transit = &P[postD];						// transit time coefficient
    Model->decay = &P[postD+ny];					// signal decay coefficient
    Model->epsilon = P[postD+ny+ny];
    Model->H = H;							// hemodynamic constants


    //printf( "decay: %f %f %f %f\n", Model->decay[0], Model->decay[1], Model->decay[2], Model->decay[3] );
    //printf( "transit: %f %f %f %f\n", Model->transit[0], Model->transit[1], Model->transit[2], Model->transit[3] );
    //printf( "epsilon: %f\n", Model->epsilon );


    J = gsl_matrix_view_array( f_out, Model->ny*states * Model->ny*states, Model->ny*states );

    dfdxdx( Model, &J.matrix );

    free( Model );

}



REDCMC_EXPORT void _r_gx_fmri( double* x, double* u, double* P, int* M, double* f_out )
{

    gx_par* Response = (gx_par*) malloc( sizeof(gx_par) );
    int ny = M[0];
    int nyny = ny*ny;
    int xlen = M[2];

    Response->ny = M[0];						// number of response variables
    Response->twostate = (char) M[1];					// two-state DCM switch
    Response->x = x;							// state variables (neuronal+hemodynamic)
    Response->epsilon = P[0];


    gx_fmri( Response, f_out );

    free( Response );

}


REDCMC_EXPORT void _r_dgdx( double* x, double* u, double* P, int* M, double* f_out )
{

    gx_par* Response = (gx_par*) malloc( sizeof(gx_par) );
    gsl_matrix_view J;

    int ny = M[0];
    int nyny = ny*ny;
    int xlen = M[2];

    Response->ny = M[0];						// number of response variables
    Response->twostate = (char) M[1];					// two-state DCM switch
    Response->x = x;							// state variables (neuronal+hemodynamic)
    Response->epsilon = P[0];


    J = gsl_matrix_view_array( f_out, Response->ny, xlen );

    dgdx( Response, &J.matrix );

    free( Response );

}


REDCMC_EXPORT void _r_dgdxdx( double* x, double* u, double* P, int* M, double* f_out )
{

    gx_par* Response = (gx_par*) malloc( sizeof(gx_par) );
    gsl_matrix_view J;

    int ny = M[0];
    int nyny = ny*ny;
    int xlen = M[2];

    Response->ny = M[0];						// number of response variables
    Response->twostate = (char) M[1];					// two-state DCM switch
    Response->x = x;							// state variables (neuronal+hemodynamic)
    Response->epsilon = P[0];


    J = gsl_matrix_view_array( f_out, Response->ny*xlen, xlen );

    dgdxdx( Response, &J.matrix );

    free( Response );

}
