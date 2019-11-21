#include "gx_fmri.h"


int _geval( const gsl_vector* x, void* params, gsl_vector* g );
int _dgdxeval( const gsl_vector* x, void* params, gsl_vector* f );


int gx_fmri( gx_par* data, double* g0 )
{

    // declarations
    /* ------------------------------------------------------------------------ */
    unsigned int i;
    int status = 0;

    double TE, V0, ep, r0, nu0, E0;
    double k1, k2, k3;
    double v, q;

    int states = data->twostate ? 6 : 5;
    int c4 = data->ny * 0;	// 0th and...
    int c5 = data->ny * 1;	// 1st column of state submatrix

    double* x = (double*) malloc( sizeof(double) * data->ny * 2 );
    //double* x = (double*) malloc( sizeof(double) * data->ny * states );

    // initialization
    /* ------------------------------------------------------------------------ */
    arr_cpyArray( data->ny*2, x, &data->x[data->ny*3] );
    //arr_cpyArray( data->ny*states, x, data->x );

//for( i=0; i<data->ny*states; ++i )
//{
//    printf("%f ",x[i] );
//}


    /* ======================================================================== */
    // (1) - Biophysical constants for 1.5T
    /* ======================================================================== */

    // time to echo (TE) (default 0.04 sec)
    /* ------------------------------------------------------------------------ */
    TE = 0.04;

    // resting venous volume (%)
    /* ------------------------------------------------------------------------ */
    V0 = 4.0;

    // estimated region-specific ratios of intra- to extra-vascular signal
    /* ------------------------------------------------------------------------ */
    ep = exp( data->epsilon );

    // slope r0 of intravascular relaxation rate R_iv as a function of oxygen
    // saturation S: R_iv = r0*[(1 - S)-(1 - S0)] (Hz)
    /* ------------------------------------------------------------------------ */
    r0 = 25.0;

    // frequency offset at the outer surface of magnetized vessels (Hz)
    /* ------------------------------------------------------------------------ */
    nu0 = 40.3;

    // resting oxygene extraction fraction
    /* ------------------------------------------------------------------------ */
    E0 = 0.4;



    /* ======================================================================== */
    // (2) - Coefficients in BOLD signal model
    /* ======================================================================== */
    k1 = 4.3 * nu0 * E0 * TE;
    k2 = ep * r0 * E0 * TE;
    k3 = 1 - ep;



    /* ======================================================================== */
    // (3) - Output equation of BOLD signal model
    /* ======================================================================== */
    for ( i = 0; i < data->ny; ++i )
    {
    //printf("%.2f %.2f ", x[i*states+3], x[i*states+4]);
	v = exp( x[c4+i] );
	q = exp( x[c5+i] );
	//v = exp( x[i*states+3] );
	//q = exp( x[i*states+4] );

	g0[i] = V0 * ( k1 - k1*q + k2 - k2*q/v + k3 - k3*v );
    }



    free( x );


    return status = 0;

}


int dgdx( gx_par* Bold, gsl_matrix* J )
{

    int status = 0;

    int states = Bold->twostate ? 6 : 5;
    int x_size = Bold->ny * states;

    gsl_vector* g0 = gsl_vector_alloc( Bold->ny );
    gsl_vector_view x = gsl_vector_view_array( Bold->x, x_size );
    gsl_multifit_function_fdf f;

    f.f = &_geval;
    f.df = NULL;
    f.fdf = NULL;
    f.n = Bold->ny;
    f.p = x_size;
    f.params = Bold;


    //status = gsl_multifit_fdfsolver_dif_fdf( &x.vector, &f, g0, J );
    status = gsl_multifit_fdfsolver_dif_fdf_simple( &x.vector, &f, g0, J );


    gsl_vector_free( g0 );

    return status;

}


int dgdxdx( gx_par* Bold, gsl_matrix* J )
{
    int i, j, k;
    int status = 0;
    double tmp;

    int states = Bold->twostate ? 6 : 5;
    int g0_size = Bold->ny * states;

    gsl_vector* gx = gsl_vector_alloc( Bold->ny * g0_size );
    gsl_vector_view x = gsl_vector_view_array( Bold->x, g0_size );
    gsl_multifit_function_fdf f;

    f.f = &_dgdxeval;
    f.df = NULL;
    f.fdf = NULL;
    f.n = Bold->ny * g0_size;
    f.p = g0_size;
    f.params = Bold;


    //status = gsl_multifit_fdfsolver_dif_fdf( &x.vector, &f, gx, J );
    status = gsl_multifit_fdfsolver_dif_fdf_simple( &x.vector, &f, gx, J );


    gsl_vector_free( gx );

    return status;
}




int _geval( const gsl_vector* x, void* params, gsl_vector* g )
{

    unsigned int i;
    int status = 0;

    gx_par* Bold = (gx_par*) params;
    int g0_size = Bold->ny;

    double* g0 = (double*) malloc( sizeof(double) * g0_size );


    status = gx_fmri( Bold, g0 );


    for ( i = 0; i < g0_size; ++i )
    {
	gsl_vector_set( g, i, g0[i] );
    }

    free( g0 );

    return status;

}


int _dgdxeval( const gsl_vector* x, void* params, gsl_vector* f )
{

    int i, j;
    int status = 0;

    gx_par* Bold = (gx_par*) params;

    int states = Bold->twostate ? 6 : 5;
    int g0_size = Bold->ny * states;

    gsl_matrix* J = gsl_matrix_calloc( Bold->ny, g0_size );


    status = dgdx( Bold, J );


    for ( i = 0; i < g0_size; ++i )
    {
	for ( j = 0; j < Bold->ny; ++j )
	{
	    gsl_vector_set( f, i*Bold->ny + j, gsl_matrix_get( J, j, i ) ); // TODO test if coordinate is (j,i) or (i,j)
	}
    }


    gsl_matrix_free( J );

    return status;

}




