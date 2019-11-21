/*##############################################################################################*/
// source:       bireduce.c
//------------------------------------------------------------------------------------------------
// type:         Function
//------------------------------------------------------------------------------------------------
// title:        Reduction of a fully nonlinear MIMO system to Bilinear form
//------------------------------------------------------------------------------------------------
// description:  Returns Matrix operators for the Bilinear approximation to the MIMO
//               system described by
//
//                   dx/dt = f(x,u,P)
//                   y(t)  = g(x,u,P)
//
//               A Bilinear approximation is returned, where the states are
//
//                   q(t)  = [1; x(t) - x(0)]
//
//               Input:    M   - model specification object
//                         P   - model parameters object
//
//               Output:   M0  - Bilinear operator - M0
//                         M1  - Bilinear operator - M1 = dM0/du
//                         [L1 - Output operators  - L1]
//                         [L2 - Output operators  - L2 = dL1/dx]
//------------------------------------------------------------------------------------------------
// based on:     spm_bireduce.m
//               DCM12 (http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
//------------------------------------------------------------------------------------------------
// reference:    ---
//------------------------------------------------------------------------------------------------
// TODO:         ---
/*##############################################################################################*/

#ifndef REDUCE_MIMO_H
#define REDUCE_MIMO_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

#include "fx_fmri.h"
#include "gx_fmri.h"
#include "gsl_multifit_nlin_ext.h"
#include "arrays.h"


typedef fx_par red_par;


int bireduce( red_par* data, double* M0, double* M1, double* L1, double* L2 );
int soreduce( red_par* data, double* M0, double* M1, double* M2, double* L1, double* L2 );
int reduce_out( red_par* data, double* M0, double* M1, double* L1, double* L2 );


#endif		// REDUCE_MIMO_H
