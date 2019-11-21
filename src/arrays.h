/*##############################################################################################*/
// source:       arrays.c
//------------------------------------------------------------------------------------------------
// type:         Function
//------------------------------------------------------------------------------------------------
// title:        Functions for C array operations
//------------------------------------------------------------------------------------------------
// description:
//------------------------------------------------------------------------------------------------
// based on:
//------------------------------------------------------------------------------------------------
// reference:    ---
//------------------------------------------------------------------------------------------------
// TODO:         ---
/*##############################################################################################*/

#ifndef ARR_FUNC_H
#define ARR_FUNC_H


double* arr_cpyArray( const unsigned int len, double* a, const double* b );

double* arr_addScalar( const unsigned int len, double* a, const double s, const double* b );
double* arr_addScalar2( const unsigned int len, double* a, const double m, const double s, const double* b );
double* arr_mulScalar( const unsigned int len, double* a, const double s, const double* b );
double* arr_divScalar( const unsigned int len, double* a, const double s, const double* b );
double* arr_powScalar( const unsigned int len, double* a, const double s, const double* b );

double* arr_scalarDiv( const unsigned int len, const double s, double* a, const double* b );
double* arr_scalarPow( const unsigned int len, const double s, double* a, const double* b );

double* arr_expArray( const unsigned int len, double* a, const double s, const double* b );

double* arr_addArray( const unsigned int len, double* a, const double* b, const double s );
double* arr_mulArray( const unsigned int len, double* a, const double* b, const double s );
double* arr_divArray( const unsigned int len, double* a, const double* b, const double s );
double* arr_powArray( const unsigned int len, double* a, const double* b );

int arr_unique( const int len, double* a );
double arr_norm1( const double* x, const int cols, const int rows );
double* arr_evalUnary( double (*fn)(double), const unsigned int len, double* a, const double* b, const double s );

int* arr_i_addScalar( const unsigned int len, int* a, const int s, const int* b );
int* arr_i_cpyArray( const unsigned int len, int* a, const int* b );
int* arr_i_mulScalar( const unsigned int len, int* a, const int s, const int* b );
int arr_i_unique( const int len, int* a );

int* arr_dtoi( const int len, int* a, const double* b );



#endif		//ARR_FUNC_H
