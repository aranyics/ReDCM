/*##############################################################################################*/
// source:       gx_fmri.c
//------------------------------------------------------------------------------------------------
// type:         Function
//------------------------------------------------------------------------------------------------
// title:        Simulated BOLD response to input
//------------------------------------------------------------------------------------------------
// description:  Returns a predicted response from state variables passed through the
//               output nonlinearity:
//
//                   g   = V0 * (k1(1 - q) + k2(1 - q/v) + k3(1 - v))
//
//               Input:    x     - state variables
//                         P     - model parameters object
//
//               Output:   g     - BOLD response
//                         [dgdx - dg/dx]
//------------------------------------------------------------------------------------------------
// based on:     spm_gx_fmri.m
//               DCM12 (http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
//------------------------------------------------------------------------------------------------
// reference:    BOLD signal model
//                   1.  Stephan KE, Weiskopf N, Drysdale PM, Robinson PA, Friston KJ (2007)
//                       Comparing hemodynamic models with DCM. NeuroImage 38: 387-401.
//------------------------------------------------------------------------------------------------
// TODO:         - return derivatives
/*##############################################################################################*/

#ifndef GX_FMRI_H
#define GX_FMRI_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

#include "gsl_multifit_nlin_ext.h"
#include "arrays.h"


typedef struct{
    double* x;				// state variables
    int ny;				// number of response variables
    double epsilon;			// ???intra- to extra-vascular signal ratio parameter
    char twostate;
} gx_par;



int gx_fmri( gx_par* data, double* g0 );

int dgdx( gx_par* Bold, gsl_matrix* J );
int dgdxdx( gx_par* Bold, gsl_matrix* J );




#endif		// GX_FMRI_H