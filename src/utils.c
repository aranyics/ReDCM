#include<stdio.h>
#include "utils.h"



int kronecker_dm( const gsl_matrix* x, const gsl_matrix* y, gsl_matrix* A )
{

    int status = 0;
    unsigned int i, j;
    gsl_matrix_view subA;

    if ( x->size1 * y->size1 != A->size1 || x->size2 * y->size2 != A->size2 )
    {
	printf( "kronecker: wrong size of matrices\n" );
	status++;
	return ( status );
    }

    for ( i = 0; i < x->size1; ++i )
    {
	for ( j = 0; j < x->size2; ++j )
	{
	    subA = gsl_matrix_submatrix( A, i*y->size1, j*y->size2, y->size1, y->size2 );
	    gsl_matrix_memcpy( &subA.matrix, y );
	    gsl_matrix_scale( &subA.matrix, gsl_matrix_get( x, i, j ) );
	}
    }

    return ( status );

}


int invert_dm( const gsl_matrix* X, gsl_matrix* Inv )
{

    int status = 0;
    int size = X->size1;
    unsigned int i;

    gsl_permutation* perm = gsl_permutation_alloc( size );
    int signum;
    gsl_matrix* LU = gsl_matrix_alloc( size, size );

    gsl_matrix* identMx = gsl_matrix_alloc( size, size );
    gsl_vector_view colId, colInv;

    if ( (X->size1 - X->size2) || (X->size1 - Inv->size1) || ( Inv->size1 - Inv->size2 ) )
    {
	printf( "invert: matrix is not square\n" );
	status++;
	return ( status );
    }

    gsl_matrix_set_identity( identMx );

    gsl_matrix_memcpy( LU, X );
    status += gsl_linalg_LU_decomp( LU, perm, &signum );

    for ( i = 0; i < size; ++i )
    {
	colId = gsl_matrix_column( identMx, i );
	colInv = gsl_matrix_column( Inv, i );
	status += gsl_linalg_LU_solve( LU, perm, &colId.vector, &colInv.vector );
    }

    gsl_matrix_free( LU );
    gsl_matrix_free( identMx );
    gsl_permutation_free( perm );

    return ( status );

}


double trace_dm( gsl_matrix* X )
{

    double trace = 0.0;
    int size = X->size1;
    unsigned int i;

    gsl_vector_view diag;

    if ( X->size1 - X->size2 )
    {
	printf( "trace: matrix is not square\n" );
	return 1.0;
    }

    diag = gsl_matrix_diagonal( X );
    for ( i = 0; i < size; ++i )
    {
	trace += gsl_vector_get( &diag.vector, i );
    }

    return ( trace );

}
