/*##############################################################################################*/
// source:       fx_fmri.c
//------------------------------------------------------------------------------------------------
// type:         Function
//------------------------------------------------------------------------------------------------
// title:        State equation for a dynamic model of fMRI responses
//------------------------------------------------------------------------------------------------
// description:  Solves state equations at time point u for a dynamic [bilinear/nonlinear/Balloon]
//               model of fMRI responses
//
//               Input:    x     - state variables
//                                   x[,1]   - excitatory neuronal activity      ue
//                                   x[,2]   - vascular signal                   s
//                                   x[,3]   - rCBF                              ln(f)
//                                   x[,4]   - venous volume                     ln(v)
//                                   x[,5]   - deoxyHb                           ln(q)
//                                   [x[,6]  - inhibitory neuronal activity      ui]
//                         u     - external input time point
//                         P     - model parameters object
//                         M     - model specifiactions object
//
//               Output:   f     - dx/dt
//                         [dfdx - df/dx]
//                         [dfdu - df/du]
//                         [D    - delays]
//------------------------------------------------------------------------------------------------
// based on:     spm_fx_fmri.m
//               DCM12 (http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
//------------------------------------------------------------------------------------------------
// reference:    hemodynamic and neuronal state equation
//                   1.  Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
//                       changes during brain activation: The Balloon model. MRM 39:855-864,
//                       1998.
//                   2.  Friston KJ, Mechelli A, Turner R, Price CJ. Nonlinear responses in
//                       fMRI: the Balloon model, Volterra kernels, and other hemodynamics.
//                       Neuroimage 12:466-477, 2000.
//                   3.  Stephan KE, Kasper L, Harrison LM, Daunizeau J, den Ouden HE,
//                       Breakspear M, Friston KJ. Nonlinear dynamic causal models for fMRI.
//                       Neuroimage 42:649-662, 2008.
//                   4.  Marreiros AC, Kiebel SJ, Friston KJ. Dynamic causal modelling for
//                       fMRI: a two-state model.
//                       Neuroimage. 2008 Jan 1;39(1):269-78.
//------------------------------------------------------------------------------------------------
// TODO:         - implement 2-state branch
//               - implement nonlinear branch
//               - return derivatives
/*##############################################################################################*/

#ifndef FX_FMRI_H
#define FX_FMRI_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

#include "gsl_multifit_nlin_ext.h"
#include "arrays.h"



typedef struct{
    double* x;			// state variables (neuronal+hemodynamic)
    double* u;			// time-point of input function
    double* A;			// linear parameters
    double* B;			// bilinear parameters
    double* C;			// exogenous parameters
    double* D;			// nonlinear parameters
    double* decay;		// signal decay coefficient
    double* transit;		// transit time coefficient
    double epsilon;
    double* H;			// hemodynamic constants
    int nu;			// number of exogenous inputs
    int ny;			// number of response variables
    char nonlin;		// nonlinear DCM switch
    char twostate;		// two-state DCM switch
    char symmetry;		// symmetric connectivity matrix switch
    char backward;		// backward connections
    char A_eigen;		// average connections encodes eigenvalue
} fx_par;



int fx_fmri ( fx_par* data, double* f0 );

int dfdx( fx_par* Model, gsl_matrix* J );
int dfdxdu( fx_par* Model, gsl_matrix* J );
int dfdxdx( fx_par* Model, gsl_matrix* J );
int dfdu( fx_par* Model, gsl_matrix* J );


#endif		// FX_FMRI_H
