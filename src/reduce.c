#include "reduce.h"



int bireduce( red_par* data, double* M0, double* M1, double* L1, double* L2 )
{

    // Declarations and memory allocation
    // -------------------------------------------------------------------------

    int status;
    unsigned int i, j, k, w, z;

    int m = data->nu;
    int n = (!data->twostate) ? 5 * data->ny : 6 * data->ny;
    int l = data->ny;

    fx_par* fxData = (fx_par*) data;

    double* f0 = (double*) malloc( sizeof(double) * n );
    gsl_matrix* dfdx0 = gsl_matrix_calloc( n, n );
    gsl_matrix* dfdu0 = gsl_matrix_calloc( n, m );
    gsl_matrix* dfdxdu0 = gsl_matrix_calloc( n*n, m );

    gsl_matrix_view M;
    gsl_vector_view VM;
    gsl_vector* V = gsl_vector_alloc( n*n );
    double* U = (double*) malloc( sizeof(double) * n );
    gsl_vector_view xv = gsl_vector_view_array( data->x, n );
    gsl_vector* prodx = gsl_vector_alloc( n );

    // Partial derivatives for 1st order Bilinear operators
    // -------------------------------------------------------------------------

    fx_fmri( fxData, f0 );
    dfdx( fxData, dfdx0 );
    dfdu( fxData, dfdu0 );
    dfdxdu( fxData, dfdxdu0 );


/*printf("dfdx\n");
    for ( i = 0; i < n; ++i )
    {
	for ( j = 0; j < n; ++j )
	{
	    printf( "%.4f\t", gsl_matrix_get( dfdx0, i, j ) );
	}
	printf( "\n" );
    }
printf("dfdu\n");
    for ( i = 0; i < n; ++i )
    {
	for ( j = 0; j < m; ++j )
	{
	    printf( "%.4f\t", gsl_matrix_get(dfdu0, i, j) );
	}
	printf( "\n" );
    }
printf("dfdxdu1\n");
VM = gsl_matrix_column( dfdxdu0, 0 );
printf("  test \n");
    for ( i = 0; i < n; ++i )
    {
	for ( j = 0; j < n; ++j )
	{
	    printf( "%.4f\t", gsl_vector_get( &VM.vector, i*n+j ) );
	}
	printf( "\n" );
    }
printf("dfdxdu2\n");
VM = gsl_matrix_column( dfdxdu0, 1 );
    for ( i = 0; i < n; ++i )
    {
	for ( j = 0; j < n; ++j )
	{
	    printf( "%.4f\t", gsl_vector_get( &VM.vector, i*n+j ) );
	}
	printf( "\n" );
    }
printf("dfdxdu3\n");
VM = gsl_matrix_column( dfdxdu0, 2 );
    for ( i = 0; i < n; ++i )
    {
	for ( j = 0; j < n; ++j )
	{
	    printf( "%.4f\t", gsl_vector_get( &VM.vector, i*n+j ) );
	}
	printf( "\n" );
    }*/


    // Bilinear operators
    // -------------------------------------------------------------------------
    // --- M0 ---

    status = gsl_blas_dgemv (CblasNoTrans, 1.0, dfdx0, &xv.vector, 0.0, prodx);
    arr_addArray( n, f0, prodx->data, -1.0 );    // f0 = f0 - dfdx %*% xv

    k = 0;
    w = 0;
    for ( i = 0; i < n+1; ++i )
    {
	for ( j = 0; j < n+1; ++j )
	{
	    if ( !i )
	    {
		M0[k++] = 0;
	    }
	    else if ( !j )
	    {
		M0[k++] = f0[w++];
	    }
	    else
	    {
		M0[k++] = gsl_matrix_get( dfdx0, i-1, j-1 );
	    }
	}
    }

    // --- M1 = dM0/du ---

    k = 0;
    for ( z = 0; z < m; ++z )
    {

    VM = gsl_matrix_column( dfdxdu0, z );
    for ( i = 0; i < dfdxdu0->size1; ++i )
    {
	V->data[i] = gsl_vector_get( &VM.vector, i );
    }
    M = gsl_matrix_view_vector( V, n, n );     // dfdxdu[,,z]
    //gsl_matrix_transpose( &M.matrix);

    VM = gsl_matrix_column( dfdu0, z );     // dfdu[,z]
    for ( i = 0; i < dfdu0->size1; ++i )
    {
	U[i] = gsl_vector_get( &VM.vector, i );
    }

    status = gsl_blas_dgemv (CblasNoTrans, 1.0, &M.matrix, &xv.vector, 0.0, prodx);
    arr_addArray( n, U, prodx->data, -1.0 );    // dfdu[,i] - dfdxdu[,,i] %*% xv


    w = 0;
    for ( i = 0; i < n+1; ++i )
    {
	for ( j = 0; j < n+1; ++j )
	{
	    if ( !i )
	    {
		M1[k++] = 0;
	    }
	    else if ( !j )
	    {
		M1[k++] = U[w++];
	    }
	    else
	    {
		M1[k++] = gsl_matrix_get( &M.matrix, i-1, j-1 );
	    }
	}
    }


    }


    // Free memory
    // -------------------------------------------------------------------------

    free( f0 );
    gsl_matrix_free( dfdx0 );
    gsl_matrix_free( dfdu0 );
    gsl_matrix_free( dfdxdu0 );
    gsl_vector_free( prodx );
    gsl_vector_free( V );
    free( U );


    return status;

}




int soreduce( red_par* data, double* M0, double* M1, double* M2, double* L1, double* L2 )
{

    // Declarations and memory allocation
    // -------------------------------------------------------------------------

    int status;
    unsigned int i, j, k, w, z;

    int m = data->nu;
    int n = (!data->twostate) ? 5 * data->ny : 6 * data->ny;
    int l = data->ny;

    fx_par* fxData = (fx_par*) data;

    double* f0 = (double*) malloc( sizeof(double) * n );
    gsl_matrix* dfdx0 = gsl_matrix_calloc( n, n );
    gsl_matrix* dfdu0 = gsl_matrix_calloc( n, m );
    gsl_matrix* dfdxdu0 = gsl_matrix_calloc( n*n, m );
    gsl_matrix* dfdxdx0 = gsl_matrix_calloc( n*n, n );

    gsl_matrix_view M;
    gsl_vector_view VM;
    gsl_vector* V = gsl_vector_alloc( n*n );
    double* U = (double*) malloc( sizeof(double) * n );
    gsl_vector_view xv = gsl_vector_view_array( data->x, n );
    gsl_vector* prodx = gsl_vector_alloc( n );

    // Partial derivatives for 1st order Bilinear operators
    // -------------------------------------------------------------------------

    fx_fmri( fxData, f0 );
    dfdx( fxData, dfdx0 );
    dfdu( fxData, dfdu0 );
    dfdxdu( fxData, dfdxdu0 );
    dfdxdx( fxData, dfdxdx0 );


/*printf("dfdx\n");
    for ( i = 0; i < n; ++i )
    {
	for ( j = 0; j < n; ++j )
	{
	    printf( "%.4f\t", gsl_matrix_get( dfdx0, i, j ) );
	}
	printf( "\n" );
    }
printf("dfdu\n");
    for ( i = 0; i < n; ++i )
    {
	for ( j = 0; j < m; ++j )
	{
	    printf( "%.4f\t", gsl_matrix_get(dfdu0, i, j) );
	}
	printf( "\n" );
    }
printf("dfdxdu1\n");
VM = gsl_matrix_column( dfdxdu0, 0 );
printf("  test \n");
    for ( i = 0; i < n; ++i )
    {
	for ( j = 0; j < n; ++j )
	{
	    printf( "%.4f\t", gsl_vector_get( &VM.vector, i*n+j ) );
	}
	printf( "\n" );
    }
printf("dfdxdu2\n");
VM = gsl_matrix_column( dfdxdu0, 1 );
    for ( i = 0; i < n; ++i )
    {
	for ( j = 0; j < n; ++j )
	{
	    printf( "%.4f\t", gsl_vector_get( &VM.vector, i*n+j ) );
	}
	printf( "\n" );
    }
printf("dfdxdu3\n");
VM = gsl_matrix_column( dfdxdu0, 2 );
    for ( i = 0; i < n; ++i )
    {
	for ( j = 0; j < n; ++j )
	{
	    printf( "%.4f\t", gsl_vector_get( &VM.vector, i*n+j ) );
	}
	printf( "\n" );
    }
printf("dfdxdx0\n");
VM = gsl_matrix_column( dfdxdx0, 1 );
    for ( i = 0; i < n; ++i )
    {
	for ( j = 0; j < n; ++j )
	{
	    printf( "%.4f\t", gsl_vector_get( &VM.vector, i*n+j ) );
	}
	printf( "\n" );
    }*/


    // Bilinear operators
    // -------------------------------------------------------------------------
    // --- M0 ---

    status = gsl_blas_dgemv (CblasNoTrans, 1.0, dfdx0, &xv.vector, 0.0, prodx);
    arr_addArray( n, f0, prodx->data, -1.0 );    // f0 = f0 - dfdx %*% xv

    k = 0;
    w = 0;
    for ( i = 0; i < n+1; ++i )
    {
	for ( j = 0; j < n+1; ++j )
	{
	    if ( !i )
	    {
		M0[k++] = 0;
	    }
	    else if ( !j )
	    {
		M0[k++] = f0[w++];
	    }
	    else
	    {
		M0[k++] = gsl_matrix_get( dfdx0, i-1, j-1 );
	    }
	}
    }

    // --- M1 = dM0/du ---

    k = 0;
    for ( z = 0; z < m; ++z )
    {

    VM = gsl_matrix_column( dfdxdu0, z );
    for ( i = 0; i < dfdxdu0->size1; ++i )
    {
	V->data[i] = gsl_vector_get( &VM.vector, i );
    }
    M = gsl_matrix_view_vector( V, n, n );     // dfdxdu[,,z]
    //gsl_matrix_transpose( &M.matrix);

    VM = gsl_matrix_column( dfdu0, z );     // dfdu[,z]
    for ( i = 0; i < dfdu0->size1; ++i )
    {
	U[i] = gsl_vector_get( &VM.vector, i );
    }

    status = gsl_blas_dgemv (CblasNoTrans, 1.0, &M.matrix, &xv.vector, 0.0, prodx);
    arr_addArray( n, U, prodx->data, -1.0 );    // dfdu[,i] - dfdxdu[,,i] %*% xv


    w = 0;
    for ( i = 0; i < n+1; ++i )
    {
	for ( j = 0; j < n+1; ++j )
	{
	    if ( !i )
	    {
		M1[k++] = 0;
	    }
	    else if ( !j )
	    {
		M1[k++] = U[w++];
	    }
	    else
	    {
		M1[k++] = gsl_matrix_get( &M.matrix, i-1, j-1 );
	    }
	}
    }

    }

    // --- M2 = dM0/dx ---

    k = 0;
    for ( z = 0; z < n; ++z )
    {

    VM = gsl_matrix_column( dfdxdx0, z );
    for ( i = 0; i < dfdxdx0->size1; ++i )
    {
	V->data[i] = gsl_vector_get( &VM.vector, i );
    }
    M = gsl_matrix_view_vector( V, n, n );     // dfdxdx[,,z]
    //gsl_matrix_transpose( &M.matrix);

    VM = gsl_matrix_column( dfdx0, z );     // dfdx[,z]
    for ( i = 0; i < dfdx0->size1; ++i )
    {
	U[i] = gsl_vector_get( &VM.vector, i );
    }

    status = gsl_blas_dgemv (CblasNoTrans, 1.0, &M.matrix, &xv.vector, 0.0, prodx);
    arr_addArray( n, U, prodx->data, -1.0 );    // dfdx[,i] - dfdxdx[,,i] %*% xv


    w = 0;
    for ( i = 0; i < n+1; ++i )
    {
	for ( j = 0; j < n+1; ++j )
	{
	    if ( !i )
	    {
		M2[k++] = 0;
	    }
	    else if ( !j )
	    {
		M2[k++] = U[w++];
	    }
	    else
	    {
		M2[k++] = gsl_matrix_get( &M.matrix, i-1, j-1 );
	    }
	}
    }

    }



    // Free memory
    // -------------------------------------------------------------------------

    free( f0 );
    gsl_matrix_free( dfdx0 );
    gsl_matrix_free( dfdu0 );
    gsl_matrix_free( dfdxdu0 );
    gsl_matrix_free( dfdxdx0 );
    gsl_vector_free( prodx );
    gsl_vector_free( V );
    free( U );


    return status;

}


int reduce_out( red_par* data, double* M0, double* M1, double* L1, double* L2 )
{

    // Declarations and memory allocation
    // -------------------------------------------------------------------------

    int status;
    unsigned int h, i, j, k, w, z;

    int m = data->nu;
    int n = (!data->twostate) ? 5 * data->ny : 6 * data->ny;
    int l = data->ny;

    gx_par* gxData = (gx_par*) malloc( sizeof(gx_par) );

    double* g0 = (double*) malloc( sizeof(double) * l );
    gsl_matrix* dgdx0 = gsl_matrix_calloc( l, n );
    gsl_matrix* dgdxdx0 = gsl_matrix_calloc( l*n, n );


    gsl_matrix* D = gsl_matrix_calloc( n, n );
    gsl_vector_view VM;
    gsl_vector_view VMrow;
    gsl_vector_view xv = gsl_vector_view_array( data->x, n );
    gsl_vector* prodx = gsl_vector_alloc( l );


    // Partial derivatives for 1st order Bilinear operators
    // -------------------------------------------------------------------------

    gxData->x = data->x;
    gxData->ny = data->ny;
    gxData->epsilon = data->epsilon;
    gxData->twostate = data->twostate;

    gx_fmri( gxData, g0 );
    dgdx( gxData, dgdx0 );
    dgdxdx( gxData, dgdxdx0 );


/*    printf("test\n");
    for ( i = 0; i < 80; ++i )
    {
	for ( j = 0; j < 20; ++j )
	{
	    printf("%.2f\t", gsl_matrix_get(dgdxdx0, i, j));
	}
	printf("\n");
    }*/



    // Output operators
    // -------------------------------------------------------------------------
    // --- L1 ---

    status = gsl_blas_dgemv (CblasNoTrans, 1.0, dgdx0, &xv.vector, 0.0, prodx);
    arr_addArray( l, g0, prodx->data, -1.0 );    // f0 = f0 - dfdx %*% xv

    k = 0;
    w = 0;
    for ( i = 0; i < l; ++i )
    {
	for ( j = 0; j < n+1; ++j )
	{
	    if ( !j )
	    {
		L1[k++] = g0[w++];
	    }
	    else
	    {
		L1[k++] = gsl_matrix_get( dgdx0, i, j-1 );
	    }
	}
    }

    // --- L2 = dL1/dx ---

    k = 0;
    for ( z = 0; z < l; ++z )
    {

	for ( i = 0; i < n; ++i )
	{
	    gsl_vector* tmpV = gsl_vector_alloc( l * n );

	    VM = gsl_matrix_column( dgdxdx0, i );

	    for ( j = 0; j < l; ++j )
	    {
		for ( h = 0; h < n; ++h )
		{
		    gsl_vector_set( tmpV, j*n+h, gsl_vector_get( &VM.vector, h*l+j ) );
		}
	    }
	    VMrow = gsl_vector_subvector( tmpV, z*n, n );

	    gsl_matrix_set_row( D, i, &VMrow.vector );

	    gsl_vector_free( tmpV );
	}



/*    printf("print D %d\n", z);
    for ( i = 0; i < 20; ++i )
    {
	for ( j = 0; j < 20; ++j )
	{
	    printf("%.2f\t", gsl_matrix_get(D, i, j));
	}
	printf("\n");
    }*/

        w = 0;
        for ( i = 0; i < n+1; ++i )
	{
	    for ( j = 0; j < n+1; ++j )
	    {
		if ( !i )
		{
		    L2[k++] = 0;
		}
		else if ( !j )
		{
		    L2[k++] = 0;
		}
		else
		{
		    L2[k++] = gsl_matrix_get( D, i-1, j-1 );
		}
	    }
	}

    }



    // Free memory
    // -------------------------------------------------------------------------

    free( gxData );
    free( g0 );
    gsl_matrix_free( dgdx0 );
    gsl_matrix_free( dgdxdx0 );
    gsl_vector_free( prodx );
    gsl_matrix_free( D );


    return status;

}
