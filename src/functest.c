/*********************************************************************/
/* Func test --- /application to test implemented functions for DCM/ */
/* ================================================================= */
/*                                                                   */
/*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "utils.h"
#include "arrays.h"
#include "fx_fmri.h"
#include "gx_fmri.h"
#include "reduce.h"
#include "integrate.h"

#define		SIZE	12



typedef struct{
    double* d;
    int l;
} data_t;


int findMx( const gsl_matrix* fMx, int* r, int* c, double* v )
{

    unsigned int i, j, k;
    double tmp;

    k = 0;
    for ( i = 0; i < fMx->size1; ++i )
    {
	for ( j = 0; j < fMx->size2; ++j )
	{
	    if ( tmp = gsl_matrix_get( fMx, i, j ) )
	    {
		r[k] = i;
		c[k] = j;
		v[k] = tmp;
		++k;
	    }
	}
    }

    return k;

}


double* test_fn( data_t* data, double* o )
{

    arr_cpyArray( SIZE, o, data->d );
    printf("cpy\n");

    arr_powArray( 4, &o[4], &o[0] );
    printf("add 5.0 to x[4:7]\n");

    return o;

}



int main( int argc, char* argv[] )
{

    int i, j, k;

    double* out = (double*) malloc( sizeof(double) * SIZE );
    double* tmp;
    double adat[SIZE] = {1,5,3,4,7,6,5,8,9,7,5,0};


    data_t* data = (data_t*) malloc( sizeof(data_t) );

    int *rows = (int*) malloc( sizeof(int) * SIZE);
    int *cols = (int*) malloc( sizeof(int) * SIZE);
    double *vals = (double*) malloc( sizeof(double) * SIZE);
    gsl_matrix_view fMx = gsl_matrix_view_array( adat, 3, 4 );

    double exp1[16] = {1, 3, 2, 5, 4, 5, 6, 2, 5, 3, 1, 5, 4, 3, 2, 7};
    double exp2[16] = {0};
    gsl_matrix_view mm = gsl_matrix_view_array( exp1, 4, 4 );
    gsl_matrix_view emm = gsl_matrix_view_array( exp2, 4, 4 );

    double svdtest[20] = {1, 3, 2, 5, 4, 5, 6, 2, 5, 3, 1, 5, 4, 3, 2, 7, 3, 6, 2, 1};
    gsl_matrix_view Asvd = gsl_matrix_view_array( svdtest, 5, 4 );
    gsl_matrix* Vsvd = gsl_matrix_calloc( 4, 4 );
    gsl_vector* Ssvd = gsl_vector_calloc( 4 );
    gsl_vector* worksvd = gsl_vector_calloc( 4 );

    double* M0 = (double*) malloc( sizeof(double) * 21 * 21 );
    double* M1 = (double*) malloc( sizeof(double) * 21 * 21 * 3 );
    double* M2 = (double*) malloc( sizeof(double) * 21 * 21 * 20 );
    double* L1 = (double*) malloc( sizeof(double) * 4 * 21 );
    double* L2 = (double*) malloc( sizeof(double) * 21 * 21 * 4 );







    gsl_vector_view x;
    gsl_multifit_function_fdf f;
    gsl_vector* fx = gsl_vector_alloc( 20 );
    gsl_matrix* Jx = gsl_matrix_calloc( 20, 20 );
    gsl_vector* fu = gsl_vector_alloc( 20 );
    gsl_matrix* Ju = gsl_matrix_calloc( 20, 3 );
    gsl_vector* fxu = gsl_vector_alloc( 20*20 );
    gsl_matrix* Jxu = gsl_matrix_calloc( 20*20, 3 );
    gsl_vector* fxx = gsl_vector_alloc( 20*20 );
    gsl_matrix* Jxx = gsl_matrix_calloc( 20*20, 20 );
    gsl_vector* gx = gsl_vector_alloc( 4 );
    gsl_matrix* Jgx = gsl_matrix_calloc( 4, 20 );
    gsl_vector* gxx = gsl_vector_alloc( 4*20 );
    gsl_matrix* Jgxx = gsl_matrix_calloc( 4*20, 20 );
    gsl_matrix* Jip = gsl_matrix_calloc( 12, 23 );

    fx_par* Model = (fx_par*) malloc( sizeof(fx_par) );
    double* f_out = (double*) malloc( sizeof(double) * 20 );
    red_par* MIMO = (red_par*) Model;

    gx_par* Bold = (gx_par*) malloc( sizeof(gx_par) );
    double* g_out = (double*) malloc(sizeof(double) * 4);

    int_par* ModelInt = (int_par*) malloc( sizeof(int_par) );
    //double* b_out = (double*) malloc( sizeof(double) * 4 * 3 );
    //double* x_out = (double*) malloc( sizeof(double) * 20 * 3 );
    double* b_out = (double*) malloc( sizeof(double) * 4 * 48 );
    double* x_out = (double*) malloc( sizeof(double) * 20 * 48 );


    double s[20] = {1, 3, 2, 5, 4, 5, 6, 2, 5, 3, 1, 5, 4, 3, 2, 7, 3, 6, 2, 1};
    double s0[20] = {0};
    double st[20] = {1, 4, 5, 4, 3, 3, 5, 3, 3, 6, 2, 6, 1, 2, 2, 5, 2, 5, 7, 1};
    double s024[24] = {0};
    double s24[24] = {1, 3, 2, 5, 4, 5, 6, 2, 5, 3, 1, 5, 4, 3, 2, 7, 3, 6, 2, 1, 2, 2, 2, 2};
    double u[3] = {0, 0, 0};
    double u2[3] = {0, 1, 0};
    //double A[16] = {0, 0.0078125, 0, 0, 0.0078125, 0, 0, 0.0078125, 0, 0, 0, 0.0078125, 0, 0.0078125, 0.0078125, 0};
    double A[16] = { 0.0813, -0.1826, 0, 0, -0.011, 0.1915, 0, -0.3344, 0, 0, -0.3722, 1.0657, 0, 0.1525, -0.216, 0.0106 };
    double B[48] = {0}; B[30] = -0.00001468; B[36] = 0.7849; B[45] = 1.1467;
    double C[12] = {0}; C[2] = 0.4208;
    double D[64] = {0}; D[0] = 0.87; D[63] = 0.45;
    //double decay[4] = {0};
    double decay[4] = {0.0657, -0.0071, -0.0711, 0.0562};
    //double transit[4] = {0};
    double transit[4] = {0.0781, -0.0175, -0.1160, 0.0742};
    double Hconst[5] = {0.64, 0.32, 2, 0.32, 0.4};

    double U[144] = {0}; U[78] = 1;
    //double delays[4] = {3.22, 3.22, 3.22, 3.22};
    double delays[4] = {2.0, 1.0, 2.0, 3.22};
    //double delays[4] = {0.0, 0.0, 0.0, 0.0};

    double Pr[23] = {-1.468345e-05, 7.848597e-01, 1.146719e+00, 4.207531e-01,-1.825675e-01, 1.915414e-01, 1.525485e-01,-3.722462e-01,
			-2.159942e-01, -3.344438e-01, 1.065666e+00, 1.055359e-02, 8.134213e-02, -1.103227e-02, 7.808828e-02, -1.745144e-02,
			-1.160219e-01,  7.422662e-02, 6.565322e-02, -7.083873e-03, -7.108322e-02, 5.616165e-02, 2.050480e-02};
    double V[1955] = {0}; V[12]=V[36]=V[96]=V[120]=V[167]=V[237]=V[261]=V[308]=V[332]=V[356]=V[989]=V[1128]=V[1267]=V[1659]=V[1762]=V[1786]=V[1810]=V[1834]=V[1858]=V[1882]=V[1906]=V[1930]=V[1954]=1;


    Model->x = s;			// state variables (neuronal+hemodynamic)
    Model->u = u;			// time-point of input function
    Model->A = A;			// linear parameters
    Model->B = B;			// bilinear parameters
    Model->C = C;			// exogenous parameters
    //Model->D = NULL;			// nonlinear parameters
    Model->D = D;			// nonlinear parameters
    Model->decay = decay;		// signal decay coefficient
    Model->transit = transit;		// transit time coefficient
    //Model->epsilon = 0;
    Model->epsilon = 0.0205;
    Model->H = Hconst;			// hemodynamic constants
    Model->nu = 3;			// number of exogenous inputs
    Model->ny = 4;			// number of response variables
    Model->nonlin = 0;			// nonlinear DCM switch
    Model->twostate = 0;		// two-state DCN switch
    Model->symmetry = 0;		// symmetric connectivity matrix switch
    Model->backward = 0;		// backward connections
    Model->A_eigen = 0;			// average connections encodes eigenvalue


    Bold->x = s;
    Bold->ny = 4;
    Bold->epsilon = 0;
    Bold->twostate = 0;


    ModelInt->Model = (red_par*)Model;	// Model structure
    ModelInt->TR = 3.22;		// TR
    //ModelInt->TR = 0.2013;		// TR
    ModelInt->delays = delays;		// regional delays
    ModelInt->U_ext = U;		// external input
    ModelInt->dt = 0.2013;		// microtime intervals
    ModelInt->u = 48;			// #microtime bins
    //ModelInt->ns = 3;			// #time samples
    ModelInt->ns = 48;			// #time samples

    ModelInt->Pr = Pr;
    ModelInt->V = V;
    ModelInt->Vr = 85;
    ModelInt->Vc = 23;



    x = gsl_vector_view_array( u, 3 );
    /*f.f = &_dfdxeval;
    f.df = NULL;
    f.fdf = NULL;
    f.n = 20*20;
    f.p = 3;
    f.params = Model;*/


    double krA[6] = {1, 2, 3, 4, 5, 6};
    double krB[4] = {1, 0, 2, 0};
    double krC[24] = {0};

    double solMx[9] = {2, 3, 6, 1, 2, 3, 5, 4, 2};
    gsl_matrix_view solMxv = gsl_matrix_view_array( solMx, 3, 3 );
    gsl_matrix* invMx = gsl_matrix_alloc( 3, 3 );
    gsl_matrix_view krAm, krBm, krCm;


    printf("start fx_fmri -----\n");
    fx_fmri( Model, f_out );
    printf("finish fx_fmri ----\n");

    printf("\n-----fx----------------------\nk: %d\n", k);
    for ( i = 0; i < 20; ++i )
    {
	printf("%f\t", f_out[i]);
    }
    printf("\n");

//return 0;
    //k = gsl_multifit_fdfsolver_dif_fdf( &x.vector, &f, fx, J );
    k = dfdx( Model, Jx );

    printf("\n-----Jx---------------------\nk: %d\n", k);

    for ( i = 0; i < 20; ++i )
    {
	for ( j = 0; j < 20; ++j )
	{
	    printf("%.3f\t", gsl_matrix_get(Jx, i, j));
	}
	printf("\n");
    }

//return 0;

    //k = gsl_multifit_fdfsolver_dif_fdf( &x.vector, &f, fu, Ju );
    k = dfdu( Model, Ju );


    printf("\n-----Ju---------------------\nk: %d\n", k);

    for ( i = 0; i < 20; ++i )
    {
	for ( j = 0; j < 3; ++j )
	{
	    printf("%.2f\t", gsl_matrix_get(Ju, i, j));
	}
	printf("\n");
    }


    //k = gsl_multifit_fdfsolver_dif_fdf( &x.vector, &f, fxu, Jxu );
    k = dfdxdu( Model, Jxu );


    printf("\n-----Jxu[1:20,]-------------\nk: %d\n", k);

    for ( i = 0; i < 20; ++i )
    {
	for ( j = 0; j < 3; ++j )
	{
	    printf("%.2f\t", gsl_matrix_get(Jxu, i, j));
	}
	printf("\n");
    }


    //k = gsl_multifit_fdfsolver_dif_fdf( &x.vector, &f, fxu, Jxu );
    k = dfdxdx( Model, Jxx );


    printf("\n-----Jxx[1:20,]-------------\nk: %d\n", k);

    for ( i = 0; i < 20; ++i )
    {
	for ( j = 0; j < 20; ++j )
	{
	    printf("%.2f\t", gsl_matrix_get(Jxx, i, j));
	}
	printf("\n");
    }


    // -------gx fmri --------------------

    k = gx_fmri( Bold, g_out );

    printf("\n-----gx---------------------\nk: %d\n", k);

    for ( i = 0; i < 4; ++i )
    {
        printf("%.8f\t", g_out[i]);
    }
    printf("\n");


    //k = gsl_multifit_fdfsolver_dif_fdf( &x.vector, &f, gx, Jgx );
    k = dgdx( Bold, Jgx );


    printf("\n-----Jgx--------------------\nk: %d\n", k);

    for ( i = 0; i < 4; ++i )
    {
	for ( j = 0; j < 20; ++j )
	{
	    printf("%.2f\t", gsl_matrix_get(Jgx, i, j));
	    //printf("%.2f\t", Jgx->data[j]);
	}
	printf("\n");
    }

    //k = gsl_multifit_fdfsolver_dif_fdf( &x.vector, &f, gxx, Jgxx );
    k = dgdxdx( Bold, Jgxx );


    printf("\n-----Jgxx-------------------\nk: %d\n", k);

    for ( i = 0; i < 80; ++i )
    {
	for ( j = 0; j < 20; ++j )
	{
	    printf("%.3f\t", gsl_matrix_get(Jgxx, i, j));
	}
	printf("\n");
    }


    printf("\n-----------------------------------------\n-----unary test-------------------\n");

    arr_evalUnary( ceil, SIZE, out, adat, 1.0 );
    for ( i = 0; i < SIZE; ++i )
    {
	printf("%f\t", out[i]);
    }
    printf("\n");



    k = findMx( &fMx.matrix, rows, cols, vals );
    printf("\n-----findMx test-------------------\nk: %d\n", k);
    for ( i = 0; i < k; ++i )
    {
	printf("%d\t%d\t%f\n", rows[i], cols[i], vals[i]);
    }
    //printf("\n");


    k = arr_unique( 20, st );
    printf("\n-----unique test-------------------\nk: %d\n", k);
    for ( i = 0; i < 20; ++i )
    {
	printf("%f\t", st[i]);
    }
    printf("\n");
    for ( i = 0; i < 20-k; ++i )
    {
	printf("%f\t", st[i]);
    }
    printf("\n");


    k = gsl_linalg_SV_decomp( &Asvd.matrix, Vsvd, Ssvd, worksvd );
    printf("\n-----svd test----------------------\nk: %d\nU\n", k);
    for ( i = 0; i < 5; ++i )
    {
	for ( j = 0; j < 4; ++j )
	{
	    printf("%.3f\t", gsl_matrix_get(&Asvd.matrix, i, j));
	}
	printf("\n");
    }
    printf("V\n");
    for ( i = 0; i < 4; ++i )
    {
	for ( j = 0; j < 4; ++j )
	{
	    printf("%.3f\t", gsl_matrix_get(Vsvd, i, j));
	}
	printf("\n");
    }
    printf("S\n");
    for ( i = 0; i < 4; ++i )
    {
	printf("%f\t", Ssvd->data[i]);
    }
    printf("\n");


    k = 0;
    printf("\n-----norm test-------------------\nk: %d\n", k);
    printf("norm1: %f\n", arr_norm1( s, 5, 4 ));
    //printf("\n");


    k = 0;
    printf("\n-----expm test-------------------\nk: %d\n", k);
    gsl_matrix_transpose( &mm.matrix );
    gsl_linalg_exponential_ss(&mm.matrix, &emm.matrix, GSL_PREC_DOUBLE);

    for ( i = 0; i < 4; ++i )
    {
	for ( j = 0; j < 4; ++j )
	{
	    printf("%.3f\t", emm.matrix.data[i*4+j]);
	}
	printf("\n");
    }


    k = bireduce( MIMO, M0, M1, L1, L2 );
    k = soreduce( MIMO, M0, M1, M2, L1, L2 );
    k = reduce_out( MIMO, M0, M1, L1, L2 );
    printf("\n-----------------------------------------");
    printf("\n-----soreduce test-------------------\nk: %d\n", k);
    for ( i = 0; i < 21; ++i )
    {
	for ( j = 0; j < 21; ++j )
	{
	    printf("%.3f\t", M0[i*21+j]);
	}
	printf("\n");
    }
    for ( i = 0; i < 21; ++i )
    {
	for ( j = 0; j < 21; ++j )
	{
	    printf("%.3f\t", M1[(i*21+j)+(21*21*2)]);
	}
	printf("\n");
    }
    for ( i = 0; i < 21; ++i )
    {
	for ( j = 0; j < 21; ++j )
	{
	    printf("%.3f\t", M2[(i*21+j)+(21*21*9)]);
	}
	printf("\n");
    }
    printf("\n-----reduce_out test-----------------\nk: %d\n", k);
    for ( i = 0; i < 4; ++i )
    {
	for ( j = 0; j < 21; ++j )
	{
	    printf("%.3f\t", L1[i*21+j]);
	}
	printf("\n");
    }
    for ( i = 0; i < 21; ++i )
    {
	for ( j = 0; j < 21; ++j )
	{
	    printf("%.3f\t", L2[(i*21+j)+(21*21*2)]);
	}
	printf("\n");
    }

//printf("size model* %d \n", sizeof(*Model));
    k = int_det( ModelInt, b_out );
    printf("\n-----integrate test------------------\nk: %d\n", k);
    for ( i = 0; i < 3; ++i )
    {
	for ( j = 0; j < 4; ++j )
	{
	    printf("%.3f\t", b_out[i*4+j]);
	}
	printf("\n");
    }
    printf("\n");

//printf("size model* %d \n", sizeof(*Model));
    k = int_det_D( ModelInt, b_out );
    printf("\n-----integrate test---(nonlin)-------\nk: %d\n", k);
    for ( i = 0; i < 3; ++i )
    {
	for ( j = 0; j < 4; ++j )
	{
	    printf("%.3f\t", b_out[i*4+j]);
	}
	printf("\n");
    }
    printf("\n");


/*
    k = int_det_hemodyn( ModelInt, b_out, x_out );
    printf("\n-----integrate hemodynamics test----\nk: %d\n", k);
    for ( i = 0; i < 48; ++i )
    {
	for ( j = 0; j < 20; ++j )
	{
	    printf("%.3f\t", x_out[i*20+j]);
	}
	printf("\n");
    }
    printf("\n");*/

/*
    k = int_dfdp( ModelInt, Jip );
    printf("\n-----Jip--------------------\nk: %d\n", k);

    for ( i = 0; i < 12; ++i )
    {
	for ( j = 0; j < 23; ++j )
	{
	    printf("%.2f\t", gsl_matrix_get(Jip, i, j));
	}
	printf("\n");
    }
*/


    krAm = gsl_matrix_view_array( krA, 2, 3 );
    krBm = gsl_matrix_view_array( krB, 2, 2 );
    krCm = gsl_matrix_view_array( krC, 4, 6 );
    k = kronecker_dm( &krAm.matrix, &krBm.matrix, &krCm.matrix );


    printf("\n-----kron test-------------\nk: %d\n", k);

    for ( i = 0; i < 4; ++i )
    {
	for ( j = 0; j < 6; ++j )
	{
	    printf("%.3f\t", gsl_matrix_get(&krCm.matrix, i, j));
	}
	printf("\n");
    }


    /*k = gsl_linalg_LU_decomp( &solMxv.matrix, perm, &signum );
    // k += gsl_linalg_LU_invert( &solMxv.matrix, perm, invMx );
    gsl_matrix_set_identity( identMx );
    colId = gsl_matrix_column( identMx, 0 );
    colInv = gsl_matrix_column( invMx, 0 );
    k += gsl_linalg_LU_solve( &solMxv.matrix, perm, &colId.vector, &colInv.vector );
    colId = gsl_matrix_column( identMx, 1 );
    colInv = gsl_matrix_column( invMx, 1 );
    k += gsl_linalg_LU_solve( &solMxv.matrix, perm, &colId.vector, &colInv.vector );
    colId = gsl_matrix_column( identMx, 2 );
    colInv = gsl_matrix_column( invMx, 2 );
    k += gsl_linalg_LU_solve( &solMxv.matrix, perm, &colId.vector, &colInv.vector );
    */
    k = invert_dm( &solMxv.matrix, invMx );
    printf("\n-----inverse test-----------\nk: %d\n", k);

    for ( i = 0; i < 3; ++i )
    {
	for ( j = 0; j < 3; ++j )
	{
	    printf("%.3f\t", gsl_matrix_get(invMx, i, j));
	}
	printf("\n");
    }


    double trace;
    trace = trace_dm( invMx );
    printf("\n-----trace test-------------\nk: %f\n", trace);




    free( out );
    free( data );
    free( Model );
    free( f_out );
    free( Bold );
    free( g_out );
    free( ModelInt );
    free( b_out );
    free( x_out );
    gsl_vector_free( fx );
    gsl_matrix_free( Jx );
    gsl_vector_free( fu );
    gsl_matrix_free( Ju );
    gsl_vector_free( fxu );
    gsl_matrix_free( Jxu );
    gsl_vector_free( fxx );
    gsl_matrix_free( Jxx );
    gsl_vector_free( gx );
    gsl_matrix_free( Jgx );
    gsl_vector_free( gxx );
    gsl_matrix_free( Jgxx );
    gsl_matrix_free( Jip );
    gsl_matrix_free( invMx );

    free( rows );
    free( cols );
    free( vals );

    free( M0 );
    free( M1 );
    free( M2 );
    free( L1 );
    free( L2 );

    gsl_matrix_free( Vsvd );
    gsl_vector_free( Ssvd );
    gsl_vector_free( worksvd );

/*
    data->d = adat;
    data->l = SIZE;

    for ( i = 0; i < SIZE; ++i )
    {
	printf("%.3f\t", data->d[i]);
    }
    printf("\n");


    test_fn( data, out );

    tmp = out;

    for ( i = 0; i < SIZE; ++i )
    {
	printf("%.3f\t", tmp[i]);
    }
    printf("\n");
*/

    return 0;

}
