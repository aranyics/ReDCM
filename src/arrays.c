#include <stdlib.h>
#include <math.h>

#include "arrays.h"




double* arr_cpyArray( const unsigned int len, double* a, const double* b )
{
    unsigned int i;

    for( i = 0; i < len; ++i )
    {
    	a[i] = b[i];
    }

    return a;
}


double* arr_addScalar( const unsigned int len, double* a, const double s, const double* b )
{
    unsigned int i;

    for ( i = 0; i < len; ++i )
    {
	a[i] = b[i] + s;
    }

    return a;
}


double* arr_addScalar2( const unsigned int len, double* a, const double m, const double s, const double* b )
{
    unsigned int i;

    for ( i = 0; i < len; ++i )
    {
	a[i] = b[i] * m + s;
    }

    return a;
}


double* arr_mulScalar( const unsigned int len, double* a, const double s, const double* b )
{
    unsigned int i;

    for ( i = 0; i < len; ++i )
    {
	a[i] = b[i] * s;
    }

    return a;
}


double* arr_divScalar( const unsigned int len, double* a, const double s, const double* b )
{
    unsigned int i;

    for ( i = 0; i < len; ++i )
    {
	a[i] = b[i] / s;
    }

    return a;
}


double* arr_powScalar( const unsigned int len, double* a, const double s, const double* b )
{
    unsigned int i;

    for ( i = 0; i < len; ++i )
    {
	a[i] = pow( b[i], s );
    }

    return a;
}
 

double* arr_scalarDiv( const unsigned int len, const double s, double* a, const double* b )
{
    unsigned int i;

    for ( i = 0; i < len; ++i )
    {
	a[i] = s / b[i];
    }

    return a;
}


double* arr_scalarPow( const unsigned int len, const double s, double* a, const double* b )
{
    unsigned int i;

    for ( i = 0; i < len; ++i )
    {
	a[i] = pow( s, b[i] );
    }

    return a;
}


double* arr_expArray( const unsigned int len, double* a, const double s, const double* b )
{
    unsigned int i;

    for ( i = 0; i < len; ++i )
    {
	a[i] = exp( b[i] ) * s;
    }

    return a;
}


double* arr_powArray( const unsigned int len, double* a, const double* b )
{
    unsigned int i;

    for ( i = 0; i < len; ++i )
    {
	a[i] = pow( a[i], b[i] );
    }

    return a;
}


double* arr_addArray( const unsigned int len, double* a, const double* b, const double s)
{
    unsigned int i;

    for ( i = 0; i < len; ++i )
    {
	a[i] = a[i] + b[i] * s;
    }

    return a;
}


double* arr_mulArray( const unsigned int len, double* a, const double* b, const double s)
{
    unsigned int i;

    for ( i = 0; i < len; ++i )
    {
	a[i] = a[i] * b[i] * s;
    }

    return a;
}


double* arr_divArray( const unsigned int len, double* a, const double* b, const double s)
{
    unsigned int i;

    for ( i = 0; i < len; ++i )
    {
	a[i] = a[i] / (b[i] * s);
    }

    return a;
}


int arr_unique( const int len, double* a )
{

    int i, d;

    d = 0;
    for ( i = 1; i < len; ++i )
    {
	if ( a[i] - a[i-1] == 0 )
	{
	    d++;
	}

	if ( d )
	{
	    a[i-d] = a[i];
	}
    }

    return d;

}


double arr_norm1( const double* x, const int cols, const int rows )
{

    unsigned int i, j;
    double maxCol = 0;

    for ( i = 0; i < cols; ++i )
    {
	double sum = 0;
	for ( j = 0; j < rows; ++j )
	{
//printf("%f ", x[i*rows + j]);
	    sum += x[i*rows + j];
	}
	if ( sum > maxCol )
	{
	    maxCol = sum;
	}
    }

    return maxCol;

}


double* arr_evalUnary( double (*fn)(double), const unsigned int len, double* a, const double* b, const double s )
{

    unsigned int i;

    for ( i = 0; i < len; ++i )
    {
	//printf("%f ", b[i]);
	a[i] = (*fn)( b[i] * s );
	//printf("%f \n", a[i]);
    }

    return a;

}


int* arr_i_addScalar( const unsigned int len, int* a, const int s, const int* b )
{
    unsigned int i;

    for ( i = 0; i < len; ++i )
    {
	a[i] = b[i] + s;
    }

    return a;
}


int* arr_i_mulScalar( const unsigned int len, int* a, const int s, const int* b )
{
    unsigned int i;

    for ( i = 0; i < len; ++i )
    {
	a[i] = b[i] * s;
    }

    return a;
}


int* arr_i_cpyArray( const unsigned int len, int* a, const int* b )
{
    unsigned int i;

    for( i = 0; i < len; ++i )
    {
    	a[i] = b[i];
    }

    return a;
}


int arr_i_unique( const int len, int* a )
{

    int i, d;

    d = 0;
    for ( i = 1; i < len; ++i )
    {
	if ( a[i] - a[i-1] == 0 )
	{
	    d++;
	}

	if ( d )
	{
	    a[i-d] = a[i];
	}
    }

    return d;

}


int* arr_dtoi( const int len, int* a, const double* b )
{

    int i;

    for ( i = 0; i < len; ++i  )
    {
	a[i] = (int)b[i];
    }

    return a;

}

