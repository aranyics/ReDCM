/*##############################################################################################*/
// source:       utils.c
//------------------------------------------------------------------------------------------------
// type:         Function
//------------------------------------------------------------------------------------------------
// title:        Utilities
//------------------------------------------------------------------------------------------------
// description:  Utility matrix operations
//------------------------------------------------------------------------------------------------
// based on:     ---
//------------------------------------------------------------------------------------------------
// reference:    ---
//------------------------------------------------------------------------------------------------
// TODO:         ---
/*##############################################################################################*/

#ifndef REDCM_UTILS_H
#define REDCM_UTILS_H


#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>


int kronecker_dm( const gsl_matrix* x, const gsl_matrix* y, gsl_matrix* A );
int invert_dm( const gsl_matrix* X, gsl_matrix* Inv );
double trace_dm( gsl_matrix* X );

#endif		//REDCM_UTILS_H