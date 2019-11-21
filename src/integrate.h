/*##############################################################################################*/
// source:       integrate.c
//------------------------------------------------------------------------------------------------
// type:         Function
//------------------------------------------------------------------------------------------------
// title:        Integrates a MIMO bilinear system dx/dt = f(x,u) = A*x + B*x*u + Cu + D
//------------------------------------------------------------------------------------------------
// description:  Integrates the bilinear approximation to the MIMO system described by
//
//                   dx/dt = f(x,u,P) = A*x + u*B*x + C*u + D
//                   y     = g(x,u,P) = L*x
//
//               Fast integrator that uses a bilinear approximation to the Jacobian  evaluated
//               using ReDCM_bireduce. This routine will also allow for sparse sampling of the
//               solution and delays in observing outputs, specified by M@delays.
//               It is used primarily for integrating fMRI models (see also spm_int_D)
//
//               Input:    P   - model parameters object
//                         M   - model specification object
//                         U   - external input object
//
//               Output:   y   - response
//------------------------------------------------------------------------------------------------
// based on:     spm_int.m
//               DCM12 (http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
//------------------------------------------------------------------------------------------------
// reference:    ---
//------------------------------------------------------------------------------------------------
// TODO:         - implement other integrators
/*##############################################################################################*/

#ifndef INTEGRATE_MODEL_H
#define INTEGRATE_MODEL_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "fx_fmri.h"
#include "gx_fmri.h"
#include "gsl_multifit_nlin_ext.h"
#include "arrays.h"
#include "reduce.h"


typedef struct{
    red_par* Model;			// Model structure
    double TR;				// TR
    double* delays;			// regional delays
    double* U_ext;			// external input
    double dt;				//
    int u;				// microtime bins
    int ns;				// time samples

    double* Pr;				// model parameters in reduced space
    double* V;				// transformation matrix for parameter space dimension reduction
    int Vr;				// rows of V
    int Vc;				// columns of V
} int_par;


int int_det( int_par* data, double* Y );
int int_det_D( int_par* data, double* Y );

int int_dfdp( int_par* ModelInt, gsl_matrix* J );

int int_det_hemodyn( int_par* data, double* Y, double* X );


#endif		// INTEGRATE_MODEL_H
