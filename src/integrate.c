#include <stdlib.h>
#include "integrate.h"


int _int_det_eval( const gsl_vector* x, void* params, gsl_vector* f );
int _compare_i( const void* a, const void* b );



int int_det( int_par* data, double* Y )
{

    int status = 0;
    unsigned int i, j, k;
    int* nullVect;

    red_par* Model = data->Model;
    int ns = data->ns;
    double* delays = data->delays;

    //double* U = data->U_ext;		// TODO don't!!! need to be transposed in R
    int m = Model->nu;
    int u = data->u;
    double dt = data->dt;
    double TR = data->TR;
    gsl_matrix_view U = gsl_matrix_view_array( data->U_ext, m, u );

    int states = (!Model->twostate) ? 5 : 6;
    int n = states * Model->ny;
    double* x = (double*) malloc( sizeof(double) * (n+1) );
    gsl_vector_view xview = gsl_vector_view_array( &x[0], n+1 );
    gsl_vector* xtmp = gsl_vector_alloc( n+1 );

    double* M0 = (double*) malloc( sizeof(double) * (n+1) * (n+1) );
    double* M1 = (double*) malloc( sizeof(double) * (n+1) * (n+1) * Model->nu );
    //double* L1 = (double*) malloc( sizeof(double) );
    //double* L2 = (double*) malloc( sizeof(double) );
    double *L1, *L2;
    double* J = (double*) malloc( sizeof(double) * (n+1) * (n+1) );
    double* Jtmp = (double*) malloc( sizeof(double) * (n+1) * (n+1) );
    gsl_matrix_view Jview = gsl_matrix_view_array( Jtmp, n+1, n+1 );
    gsl_matrix* Jexp;
    //gsl_matrix* Jexp = gsl_matrix_calloc( n+1, n+1 );

    double* D = (double*) malloc( sizeof(double) * Model->ny );
    int* Di = (int*) malloc( sizeof(int) * Model->ny );
    int period;
    int pCnt = 0;

    int* su = (int*) malloc( sizeof(int) );
    int* sy = (int*) malloc( sizeof(int) * Model->ny );
    int* t = (int*) malloc( sizeof(int) );
    double* tDiff;
    double tD;
    //int sampMt = ceil( u / ns );
    char sampled = 0;
    int nsu = 2;
    int nsy = 0;
    int nnsy = 0;
    int nt = 2;
    int ru = 0;
    int ry = 0;

    gsl_vector_view Ucol1;
    gsl_vector_view Ucol2;

    gx_par* Bold = (gx_par*) malloc( sizeof(gx_par) );
    double* q = (double*) malloc( sizeof(double) * Model->ny );
    int* toSample = (int*) malloc( sizeof(int) * Model->ny );
    int nSample;
    int nodeCnt = 0;


    // initialize
    // -----------------------------------------------------------------------
    nullVect = (int*) malloc( sizeof(int) * Model->ny );
    arr_i_mulScalar( Model->ny, nullVect, 0, nullVect );

    // get expansion point
    // -----------------------------------------------------------------------
    arr_cpyArray( n, &x[1], Model->x ); x[0] = 1;

    // Bilinear approximation (1st order)
    // -----------------------------------------------------------------------
    status = bireduce( Model, M0, M1, L1, L2 );
    if (status)
    {
	printf( "error bireduce\n" );
	return status;
    }

    // delays TODO: only if delays > dt
    // -----------------------------------------------------------------------
    arr_divScalar( Model->ny, D, dt, delays );
    arr_evalUnary( round, Model->ny, D, D, 1.0 );
    arr_dtoi( Model->ny, Di, D );
    arr_i_addScalar( Model->ny, Di, -1.0, Di );
    period = (int)( round( TR / dt ) );


    // Evaluation times (t) and indicator array for inputs (su) and output (sy)
    // =======================================================================
    for( i = 0; i < u; ++i )
    {
	if ( i == 0 )
	{
	    su[0] = i;
	    t[0] = i;
	    continue;
	}

	// get times that the input changes
	// -------------------------------------------------------------------
	Ucol1 = gsl_matrix_column( &U.matrix, i-1 );
	Ucol2 = gsl_matrix_column( &U.matrix, i );
	for ( j = 0; j < m; ++j )
	{
	    double dU = gsl_vector_get( &Ucol2.vector, j ) - gsl_vector_get( &Ucol1.vector, j );
	    if ( dU )
	    {
		su = (int*) realloc( su, sizeof(int) * nsu++ );
		t = (int*) realloc( t, sizeof(int) * nt++ );

		su[nsu-2] = i;
		t[nt-2] = i;

		//printf("%d ", i);

		break;
	    }
	}

	// get times that the response is sampled (TODO only works with homogene delays)
	// -------------------------------------------------------------------
	if ( (i > 0) && (((i+1) % period) == 0) )
	{
	    //printf("i: %d", i);
	    int* tmpD = (int*) malloc( sizeof(int) * Model->ny );

	    nsy += Model->ny;
	    nnsy++;
	    sy = (int*) realloc( sy, sizeof(int) * nsy );
	    arr_i_addScalar( Model->ny, tmpD, pCnt*period, Di );
	    arr_i_cpyArray( Model->ny, &sy[nsy-Model->ny], tmpD );
	    pCnt++;

	    nt += Model->ny;
	    t = (int*) realloc( t, sizeof(int) * nt );
	    arr_i_cpyArray( Model->ny, &t[nt-1-Model->ny], tmpD );

	    //printf("%d ", i);

	    free( tmpD );
	}

    }
    nt--;
    nsu--;
    //printf("bins: %d\n", i);

    if ( nsy / nnsy != Model->ny )
    {
	printf("error at nsy\n");
	return (nsy%nnsy);
    }
    if ( nnsy-ns != 0 )
    {
	if ( period == 1 && nnsy-ns != -1 )
	{
	    printf("error at ns %d %d %d\n", nnsy, ns, nnsy-ns);
	    return (nnsy);
	}
    }
    if ( nsu + nsy != nt )
    {
	printf("error at nt\n");
	return (nsu+nsy);
    }


    // time in seconds
    // -----------------------------------------------------------------------
    //nt = nt - arr_i_unique( nt, t );
    //t = (int*) realloc( t, sizeof(int) * nt );
    //printf("\n");
    //for ( i = 0; i < nt-1; ++i )
    //{
	//printf(" %d", t[i]);
    //}
    //printf("\nqsort");
/*    qsort( t, nt, sizeof(int), _compare_i );
    nt -= arr_i_unique( nt, t );
    tDiff = (double*) malloc( sizeof(double) * nt );
    for ( i = 0; i < nt-1; ++i )
    {
	tDiff[i] = dt * (t[i+1] - t[i]);
	//printf(" %d", t[i]);
    }*/
    /*printf("\n");
    for ( i = 0; i < nt; ++i )
    {
	printf(" %d", t[i]);
    }
    printf("\nqsort");*/
    qsort( t, nt, sizeof(int), _compare_i );
    /*for ( i = 0; i < nt; ++i )
    {
	printf(" %d", t[i]);
    }
    printf("\nunique");*/
    nt -= arr_i_unique( nt, t );
    /*for ( i = 0; i < nt; ++i )
    {
	printf(" %d", t[i]);
    }*/
    tDiff = (double*) malloc( sizeof(double) * nt );
    for ( i = 0; i < nt-1; ++i )
    {
	tDiff[i] = dt * (t[i+1] - t[i]);
	//printf(" %d", t[i]);
    }

//return 0;


/*printf("\n---sy-%d-%d-\n", nsy, nnsy);
for ( i = 0; i < nsy; ++i )
{
    printf("%d ", sy[i]);
}
printf("\n---su-%d--\n", nsu);
for ( i = 0; i < nsu; ++i )
{
    printf("%d ", su[i]);
}
printf("\n---t-%d--\n", nt);
for ( i = 0; i < nt; ++i )
{
    printf("%d ", t[i]);
}
printf("\n---dt-%d--\n", nt);
for ( i = 0; i < nt; ++i )
{
    printf("%f ", tDiff[i]);
}
printf("\n");*/
// return 0;


    // Integrates
    // =======================================================================
    arr_cpyArray( (n+1)*(n+1), J, M0 );

    for ( i = 0; i < nt; ++i )
    {
	int ryny = 0;
	//printf("t: %d ru: %d ry: %d\n", t[i], ru, ry);

	Jexp = gsl_matrix_calloc( n+1, n+1 );

	// input dependent changes in Jacobian
	// -------------------------------------------------------------------
	if ( t[i] == su[ru] )
	{
	    Ucol1 = gsl_matrix_column( &U.matrix, t[i] );
	    arr_cpyArray( (n+1)*(n+1), J, M0 );
	    for ( j = 0; j < m; ++j )
	    {
		double* tmpJ = (double*) malloc( sizeof(double) * (n+1)*(n+1) );
		arr_mulScalar( (n+1)*(n+1), tmpJ, gsl_vector_get( &Ucol1.vector, j ), &M1[(n+1)*(n+1)*j] );
		arr_addArray( (n+1)*(n+1), J, tmpJ, 1.0 );
		free( tmpJ );
	    }



	    ru++;
	}

	// output sampled
	// ---------------------------------------------------------------
	arr_i_cpyArray( Model->ny, toSample, nullVect );
	ryny = ry*Model->ny;
	nSample = 0;

	for ( j = 0; j < Model->ny; ++j )
	{
	    if ( t[i] == sy[ryny+j] )
	    {
		toSample[nSample++] = j;
	    }
	}

	if ( nSample )
	{
	    Bold->x = &x[1];
	    Bold->ny = Model->ny;
	    Bold->epsilon = Model->epsilon;
	    Bold->twostate = Model->twostate;

	    status += gx_fmri( Bold, q );

	    for ( j = 0; j < nSample; ++j )
	    {
		int offs = toSample[j];
		Y[ryny+offs] = q[offs];
	    }

	    nodeCnt += nSample;
	    if ( nodeCnt == Model->ny )
	    {
		ry++;
		nodeCnt = 0;
	    }
	}

	// compute updated states x = expm(J*dt)*x
	// ---------------------------------------------------------------
	arr_mulScalar( (n+1)*(n+1), Jtmp, tDiff[i], J );
	gsl_linalg_exponential_ss( &Jview.matrix, Jexp, GSL_PREC_DOUBLE );
	gsl_blas_dgemv( CblasNoTrans, 1.0, Jexp, &xview.vector, 0.0, xtmp );
	arr_cpyArray( n+1, x, xtmp->data );


	// check for convergence
	// ---------------------------------------------------------------
	if ( arr_norm1( &x[1], states, Model->ny ) > 1e6 )
	{
	    break;
	}

	gsl_matrix_free( Jexp );

    }






    free( Bold );
    free( M0 );
    free( M1 );
    //free( L1 );
    //free( L2 );
    free( J );
    free( Jtmp );
    free( x );
    free( D );
    free( su );
    free( sy );
    free( t );
    free( tDiff );
    free( nullVect );
    free( q );
    free( toSample );
    //gsl_matrix_free( Jexp );
    gsl_vector_free( xtmp );

    return status;

}



int int_det_D( int_par* data, double* Y )
{

    int status = 0;
    unsigned int i, j, k;
    int* nullVect;

    red_par* Model = data->Model;
    int ns = data->ns;
    double* delays = data->delays;

    //double* U = data->U_ext;		// TODO don't!!! need to be transposed in R
    int m = Model->nu;
    int u = data->u;
    double dt = data->dt;
    double TR = data->TR;
    gsl_matrix_view U = gsl_matrix_view_array( data->U_ext, m, u );

    int states = (!Model->twostate) ? 5 : 6;
    int n = states * Model->ny;
    double* x = (double*) malloc( sizeof(double) * (n+1) );
    gsl_vector_view xview = gsl_vector_view_array( &x[0], n+1 );
    gsl_vector* xtmp = gsl_vector_alloc( n+1 );

    double* M0 = (double*) malloc( sizeof(double) * (n+1) * (n+1) );
    double* M1 = (double*) malloc( sizeof(double) * (n+1) * (n+1) * Model->nu );
    double* M2 = (double*) malloc( sizeof(double) * (n+1) * (n+1) * n );
    //double* L1 = (double*) malloc( sizeof(double) );
    //double* L2 = (double*) malloc( sizeof(double) );
    double *L1, *L2;
    double* J = (double*) malloc( sizeof(double) * (n+1) * (n+1) );
    double* Jtmp = (double*) malloc( sizeof(double) * (n+1) * (n+1) );
    gsl_matrix_view Jview = gsl_matrix_view_array( Jtmp, n+1, n+1 );
    gsl_matrix* Jexp;
    //gsl_matrix* Jexp = gsl_matrix_calloc( n+1, n+1 );

    double* D = (double*) malloc( sizeof(double) * Model->ny );
    int* Di = (int*) malloc( sizeof(int) * Model->ny );
    int period;
    int pCnt = 0;

    int* su = (int*) malloc( sizeof(int) );
    int* sx = (int*) malloc( sizeof(int) );
    int* sy = (int*) malloc( sizeof(int) * Model->ny );
    int* t = (int*) malloc( sizeof(int) );
    double* tDiff;
    double tD;
    //int sampMt = ceil( u / ns );
    char sampled = 0;
    int nsu = 2;
    int nsx = 2;
    int nsy = 0;
    int nnsy = 0;
    int nt = 2;
    int ru = 0;
    int ry = 0;

    gsl_vector_view Ucol1;
    gsl_vector_view Ucol2;

    gx_par* Bold = (gx_par*) malloc( sizeof(gx_par) );
    double* q = (double*) malloc( sizeof(double) * Model->ny );
    int* toSample = (int*) malloc( sizeof(int) * Model->ny );
    int nSample;
    int nodeCnt = 0;


    // initialize
    // -----------------------------------------------------------------------
    nullVect = (int*) malloc( sizeof(int) * Model->ny );
    arr_i_mulScalar( Model->ny, nullVect, 0, nullVect );

    // get expansion point
    // -----------------------------------------------------------------------
    arr_cpyArray( n, &x[1], Model->x ); x[0] = 1;

    // Bilinear approximation (1st order)
    // -----------------------------------------------------------------------
    status = soreduce( Model, M0, M1, M2, L1, L2 );
    if (status)
    {
	printf( "error bireduce\n" );
	return status;
    }

    // delays TODO: only if delays > dt
    // -----------------------------------------------------------------------
    arr_divScalar( Model->ny, D, dt, delays );
    arr_evalUnary( round, Model->ny, D, D, 1.0 );
    arr_dtoi( Model->ny, Di, D );
    arr_i_addScalar( Model->ny, Di, -1.0, Di );
    period = (int)( round( TR / dt ) );


    // Evaluation times (t) and indicator array for inputs (su) and output (sy)
    // =======================================================================
    for( i = 0; i < u; ++i )
    {
	if ( i == 0 )
	{
	    su[0] = i;
	    sx[0] = i;
	    t[0] = i;
	    continue;
	}

	// get times that the input changes
	// -------------------------------------------------------------------
	Ucol1 = gsl_matrix_column( &U.matrix, i-1 );
	Ucol2 = gsl_matrix_column( &U.matrix, i );
	for ( j = 0; j < m; ++j )
	{
	    double dU = gsl_vector_get( &Ucol2.vector, j ) - gsl_vector_get( &Ucol1.vector, j );
	    if ( dU )
	    {
		su = (int*) realloc( su, sizeof(int) * nsu++ );
		t = (int*) realloc( t, sizeof(int) * nt++ );

		su[nsu-2] = i;
		t[nt-2] = i;

		//printf("%d ", i);

		break;
	    }
	}

	// get (N) intervening times
	if ( (i > 0) && ((i+1) % period) == 0 )
	{
		sx = (int*) realloc( sx, sizeof(int) * nsx++ );
		t = (int*) realloc( t, sizeof(int) * nt++ );

		sx[nsx-2] = i;
		t[nt-2] = i;

		//printf("%d ", i);
	}

	// get times that the response is sampled
	// -------------------------------------------------------------------
	if ( (i > 0) && (((i+1) % period) == 0) )
	{
	    //printf("i: %d", i);
	    int* tmpD = (int*) malloc( sizeof(int) * Model->ny );

	    nsy += Model->ny;
	    nnsy++;
	    sy = (int*) realloc( sy, sizeof(int) * nsy );
	    arr_i_addScalar( Model->ny, tmpD, pCnt*period, Di );
	    arr_i_cpyArray( Model->ny, &sy[nsy-Model->ny], tmpD );
	    pCnt++;

	    nt += Model->ny;
	    t = (int*) realloc( t, sizeof(int) * nt );
	    arr_i_cpyArray( Model->ny, &t[nt-1-Model->ny], tmpD );

	    //printf("%d ", i);

	    free( tmpD );
	}


    }
    nt--;
    nsu--;
    nsx--;
    //printf("nt: %d\n", nt);
    //printf("bins: %d\n", i);

    if ( nsy / nnsy != Model->ny )
    {
	printf("error at nsy\n");
	return (nsy%nnsy);
    }
    if ( nnsy-ns != 0 )
    {
	if ( period == 1 && nnsy-ns != -1 )
	{
	    printf("error at ns %d %d %d\n", nnsy, ns, nnsy-ns);
	    return (nnsy);
	}
    }
/*    if ( nsu + nsy != nt )
    {
	printf("error at nt\n");
	return (nsu+nsy);
    }*/


    // time in seconds
    // -----------------------------------------------------------------------
    //nt = nt - arr_i_unique( nt, t );
    //t = (int*) realloc( t, sizeof(int) * nt );
    /*printf("\n");
    for ( i = 0; i < nt; ++i )
    {
	printf(" %d", t[i]);
    }
    printf("\nqsort");*/
    qsort( t, nt, sizeof(int), _compare_i );
    /*for ( i = 0; i < nt; ++i )
    {
	printf(" %d", t[i]);
    }
    printf("\nunique");*/
    nt -= arr_i_unique( nt, t );
    /*for ( i = 0; i < nt; ++i )
    {
	printf(" %d", t[i]);
    }*/
    tDiff = (double*) malloc( sizeof(double) * nt );
    for ( i = 0; i < nt-1; ++i )
    {
	tDiff[i] = dt * (t[i+1] - t[i]);
	//printf(" %d", t[i]);
    }

//return 0;


/*printf("\n---sy-%d-%d-\n", nsy, nnsy);
for ( i = 0; i < nsy; ++i )
{
    printf("%d ", sy[i]);
}
printf("\n---su-%d--\n", nsu);
for ( i = 0; i < nsu; ++i )
{
    printf("%d ", su[i]);
}
printf("\n---t-%d--\n", nt);
for ( i = 0; i < nt; ++i )
{
    printf("%d ", t[i]);
}
printf("\n---dt-%d--\n", nt);
for ( i = 0; i < nt; ++i )
{
    printf("%f ", tDiff[i]);
}
printf("\n");*/
// return 0;


    // Integrates
    // =======================================================================
    arr_cpyArray( (n+1)*(n+1), J, M0 );

    for ( i = 1; i < nt; ++i )
    {
	int ryny = 0;
	//printf("t: %d ru: %d ry: %d\n", t[i], ru, ry);

	Jexp = gsl_matrix_calloc( n+1, n+1 );

	// input dependent changes in Jacobian
	// -------------------------------------------------------------------
	Ucol1 = gsl_matrix_column( &U.matrix, t[i] );
	arr_cpyArray( (n+1)*(n+1), J, M0 );
	for ( j = 0; j < m; ++j )
	{
	    double* tmpJ = (double*) malloc( sizeof(double) * (n+1)*(n+1) );
	    arr_mulScalar( (n+1)*(n+1), tmpJ, gsl_vector_get( &Ucol1.vector, j ), &M1[(n+1)*(n+1)*j] );
	    arr_addArray( (n+1)*(n+1), J, tmpJ, 1.0 );
	    free( tmpJ );
	}

	// state dependent changes in Jacobian
	// -------------------------------------------------------------------
	for ( j = 0; j < n; ++j )
	{
		double* tmpJ = (double*) malloc( sizeof(double) * (n+1)*(n+1) );
		arr_mulScalar( (n+1)*(n+1), tmpJ, x[j+1] - Model->x[j], &M2[(n+1)*(n+1)*j] );
		arr_addArray( (n+1)*(n+1), J, tmpJ, 1.0 );
		free( tmpJ );
	}

	// output sampled
	// ---------------------------------------------------------------
	arr_i_cpyArray( Model->ny, toSample, nullVect );
	ryny = ry*Model->ny;
	nSample = 0;

	for ( j = 0; j < Model->ny; ++j )
	{
	    if ( t[i] == sy[ryny+j] )
	    {
		toSample[nSample++] = j;
	    }
	}

	if ( nSample )
	{
	    Bold->x = &x[1];
	    Bold->ny = Model->ny;
	    Bold->epsilon = Model->epsilon;
	    Bold->twostate = Model->twostate;

	    status += gx_fmri( Bold, q );

	    for ( j = 0; j < nSample; ++j )
	    {
		int offs = toSample[j];
		Y[ryny+offs] = q[offs];
	    }

	    nodeCnt += nSample;
	    if ( nodeCnt == Model->ny )
	    {
		ry++;
		nodeCnt = 0;
	    }
	}

	// compute updated states x = expm(J*dt)*x
	// ---------------------------------------------------------------
	arr_mulScalar( (n+1)*(n+1), Jtmp, tDiff[i], J );
	gsl_linalg_exponential_ss( &Jview.matrix, Jexp, GSL_PREC_DOUBLE );
	gsl_blas_dgemv( CblasNoTrans, 1.0, Jexp, &xview.vector, 0.0, xtmp );
	arr_cpyArray( n+1, x, xtmp->data );


	// check for convergence
	// ---------------------------------------------------------------
	if ( arr_norm1( &x[1], states, Model->ny ) > 1e6 )
	{
	    break;
	}

	gsl_matrix_free( Jexp );

    }






    free( Bold );
    free( M0 );
    free( M1 );
    free( M2 );
    //free( L1 );
    //free( L2 );
    free( J );
    free( Jtmp );
    free( x );
    free( D );
    free( su );
    free( sx );
    free( sy );
    free( t );
    free( tDiff );
    free( nullVect );
    free( q );
    free( toSample );
    //gsl_matrix_free( Jexp );
    gsl_vector_free( xtmp );

    return status;

}




int int_det_hemodyn( int_par* data, double* Y, double* X )
{

    int status = 0;
    unsigned int i, j, k;
    int* nullVect;

    red_par* Model = data->Model;
    int ns = data->ns;
    double* delays = data->delays;

    //double* U = data->U_ext;		// TODO don't!!! need to be transposed in R
    int m = Model->nu;
    int u = data->u;
    double dt = data->dt;
    double TR = data->TR;
    gsl_matrix_view U = gsl_matrix_view_array( data->U_ext, m, u );

    int states = (!Model->twostate) ? 5 : 6;
    int n = states * Model->ny;
    double* x = (double*) malloc( sizeof(double) * (n+1) );
    gsl_vector_view xview = gsl_vector_view_array( &x[0], n+1 );
    gsl_vector* xtmp = gsl_vector_alloc( n+1 );

    double* M0 = (double*) malloc( sizeof(double) * (n+1) * (n+1) );
    double* M1 = (double*) malloc( sizeof(double) * (n+1) * (n+1) * 3 );
    //double* L1 = (double*) malloc( sizeof(double) );
    //double* L2 = (double*) malloc( sizeof(double) );
    double *L1, *L2;
    double* J = (double*) malloc( sizeof(double) * (n+1) * (n+1) );
    double* Jtmp = (double*) malloc( sizeof(double) * (n+1) * (n+1) );
    gsl_matrix_view Jview = gsl_matrix_view_array( Jtmp, n+1, n+1 );
    gsl_matrix_view JviewOrig = gsl_matrix_view_array( J, n+1, n+1 );
    gsl_vector_view JcolOrig = gsl_matrix_column( &JviewOrig.matrix, 0 );
    gsl_matrix* Jexp;
    //gsl_matrix* Jexp = gsl_matrix_calloc( n+1, n+1 );

    double* D = (double*) malloc( sizeof(double) * Model->ny );
    int* Di = (int*) malloc( sizeof(int) * Model->ny );
    int period;
    int pCnt = 0;

    int* su = (int*) malloc( sizeof(int) );
    int* sy = (int*) malloc( sizeof(int) * Model->ny );
    int* t = (int*) malloc( sizeof(int) );
    double* tDiff;
    double tD;
    //int sampMt = ceil( u / ns );
    char sampled = 0;
    int nsu = 2;
    int nsy = 0;
    int nnsy = 0;
    int nt = 2;
    int ru = 0;
    int ry = 0;

    gsl_vector_view Ucol1;
    gsl_vector_view Ucol2;

    gx_par* Bold = (gx_par*) malloc( sizeof(gx_par) );
    double* q = (double*) malloc( sizeof(double) * Model->ny );
    int* toSample = (int*) malloc( sizeof(int) * Model->ny );
    int nSample;
    int nodeCnt = 0;


    // initialize
    // -----------------------------------------------------------------------
    nullVect = (int*) malloc( sizeof(int) * Model->ny );
    arr_i_mulScalar( Model->ny, nullVect, 0, nullVect );

    // get expansion point
    // -----------------------------------------------------------------------
    arr_cpyArray( n, &x[1], Model->x ); x[0] = 1;

    // Bilinear approximation (1st order)
    // -----------------------------------------------------------------------
    status = bireduce( Model, M0, M1, L1, L2 );
    if (status)
    {
	printf( "error bireduce\n" );
	return status;
    }

    // delays TODO: only if delays > dt
    // -----------------------------------------------------------------------
    arr_divScalar( Model->ny, D, dt, delays );
    arr_evalUnary( round, Model->ny, D, D, 1.0 );
    arr_dtoi( Model->ny, Di, D );
    arr_i_addScalar( Model->ny, Di, -1.0, Di );
    period = (int)( round( TR / dt ) );


    // Evaluation times (t) and indicator array for inputs (su) and output (sy)
    // =======================================================================
    for( i = 0; i < u; ++i )
    {
	if ( i == 0 )
	{
	    su[0] = i;
	    t[0] = i;
	    continue;
	}

	// get times that the input changes
	// -------------------------------------------------------------------
	Ucol1 = gsl_matrix_column( &U.matrix, i-1 );
	Ucol2 = gsl_matrix_column( &U.matrix, i );
	for ( j = 0; j < m; ++j )
	{
	    double dU = gsl_vector_get( &Ucol2.vector, j ) - gsl_vector_get( &Ucol1.vector, j );
	    if ( dU )
	    {
		su = (int*) realloc( su, sizeof(int) * nsu++ );
		t = (int*) realloc( t, sizeof(int) * nt++ );

		su[nsu-2] = i;
		t[nt-2] = i;

		break;
	    }
	}

	// get times that the response is sampled (TODO only works with homogene delays)
	// -------------------------------------------------------------------
	if ( (i > 0) && (((i+1) % period) == 0) )
	{
	    //printf("i: %d", i);
	    int* tmpD = (int*) malloc( sizeof(int) * Model->ny );

	    nsy += Model->ny;
	    nnsy++;
	    sy = (int*) realloc( sy, sizeof(int) * nsy );
	    arr_i_addScalar( Model->ny, tmpD, pCnt*period, Di );
	    arr_i_cpyArray( Model->ny, &sy[nsy-Model->ny], tmpD );
	    pCnt++;

	    nt += Model->ny;
	    t = (int*) realloc( t, sizeof(int) * nt );
	    arr_i_cpyArray( Model->ny, &t[nt-1-Model->ny], tmpD );

	    free( tmpD );
	}

    }
    nt--;
    nsu--;

    if ( nsy / nnsy != Model->ny )
    {
	printf("error at nsy\n");
	return (nsy%nnsy);
    }
    if ( nnsy-ns != 0 )
    {
	if ( period == 1 && nnsy-ns != -1 )
	{
	    printf("error at ns %d %d %d\n", nnsy, ns, nnsy-ns);
	    return (nnsy);
	}
    }
    if ( nsu + nsy != nt )
    {
	printf("error at nt\n");
	return (nsu+nsy);
    }


    // time in seconds
    // -----------------------------------------------------------------------
    //nt = nt - arr_i_unique( nt, t );
    //t = (int*) realloc( t, sizeof(int) * nt );
    qsort( t, nt, sizeof(int), _compare_i );
    nt -= arr_i_unique( nt, t );
    tDiff = (double*) malloc( sizeof(double) * nt );
    for ( i = 0; i < nt-1; ++i )
    {
	tDiff[i] = dt * (t[i+1] - t[i]);
    }



/*printf("\n---sy-%d-%d-\n", nsy, nnsy);
for ( i = 0; i < nsy; ++i )
{
    printf("%d ", sy[i]);
}
printf("\n---su-%d--\n", nsu);
for ( i = 0; i < nsu; ++i )
{
    printf("%d ", su[i]);
}
printf("\n---t-%d--\n", nt);
for ( i = 0; i < nt; ++i )
{
    printf("%d ", t[i]);
}
printf("\n---dt-%d--\n", nt);
for ( i = 0; i < nt; ++i )
{
    printf("%f ", tDiff[i]);
}
printf("\n");*/
// return 0;


    // Integrates
    // =======================================================================
    arr_cpyArray( (n+1)*(n+1), J, M0 );

    for ( i = 0; i < nt; ++i )
    {
	int ryny = 0;
	//printf("t: %d ru: %d ry: %d\n", t[i], ru, ry);

	Jexp = gsl_matrix_calloc( n+1, n+1 );

	// input dependent changes in Jacobian
	// -------------------------------------------------------------------
	if ( t[i] == su[ru] )
	{
	    Ucol1 = gsl_matrix_column( &U.matrix, t[i] );
	    arr_cpyArray( (n+1)*(n+1), J, M0 );
	    for ( j = 0; j < m; ++j )
	    {
		double* tmpJ = (double*) malloc( sizeof(double) * (n+1)*(n+1) );
		arr_mulScalar( (n+1)*(n+1), tmpJ, gsl_vector_get( &Ucol1.vector, j ), &M1[(n+1)*(n+1)*j] );
		arr_addArray( (n+1)*(n+1), J, tmpJ, 1.0 );
		free( tmpJ );
	    }



	    ru++;
	}

	// output sampled
	// ---------------------------------------------------------------
	arr_i_cpyArray( Model->ny, toSample, nullVect );
	ryny = ry*Model->ny;
	nSample = 0;

	for ( j = 0; j < Model->ny; ++j )
	{
	    if ( t[i] == sy[ryny+j] )
	    {
		toSample[nSample++] = j;
	    }
	}

	if ( nSample )
	{
	    Bold->x = &x[1];
	    Bold->ny = Model->ny;
	    Bold->epsilon = Model->epsilon;
	    Bold->twostate = Model->twostate;

	    status += gx_fmri( Bold, q );

	    for ( j = 0; j < nSample; ++j )
	    {
		int offs = toSample[j];
		Y[ryny+offs] = q[offs];
	    }

	    nodeCnt += nSample;
	    if ( nodeCnt == Model->ny )
	    {
		ry++;
		nodeCnt = 0;
	    }
	}

	// compute updated states x = expm(J*dt)*x
	// ---------------------------------------------------------------
	arr_mulScalar( (n+1)*(n+1), Jtmp, tDiff[i], J );
	gsl_linalg_exponential_ss( &Jview.matrix, Jexp, GSL_PREC_DOUBLE );
	gsl_blas_dgemv( CblasNoTrans, 1.0, Jexp, &xview.vector, 0.0, xtmp );
	arr_cpyArray( n+1, x, xtmp->data );

	if ( nSample )
	{
	    //arr_cpyArray( n, &X[(ry-1)*(n)], &x[1] );
	    for ( j = 0; j < n+1; ++j )
	    {
		X[(ry-1)*n + j] = gsl_vector_get( &JcolOrig.vector, j );
	    }
	}


	// check for convergence
	// ---------------------------------------------------------------
	if ( arr_norm1( &x[1], states, Model->ny ) > 1e6 )
	{
	    break;
	}

	gsl_matrix_free( Jexp );

    }






    free( Bold );
    free( M0 );
    free( M1 );
    //free( L1 );
    //free( L2 );
    free( J );
    free( Jtmp );
    free( x );
    free( D );
    free( su );
    free( sy );
    free( t );
    free( tDiff );
    free( nullVect );
    free( q );
    free( toSample );
    //gsl_matrix_free( Jexp );
    gsl_vector_free( xtmp );

    return status;

}





int int_dfdp( int_par* ModelInt, gsl_matrix* J )
{
    //int i, j;
    int status = 0;

    int f0_size = ModelInt->ns * ModelInt->Model->ny;
    //int J_size = f0_size * ModelInt->Vc;

    gsl_vector* fx = gsl_vector_alloc( f0_size );
    gsl_vector_view x = gsl_vector_view_array( ModelInt->Pr, ModelInt->Vc );
    gsl_multifit_function_fdf f;
;

    f.f = &_int_det_eval;
    f.df = NULL;
    f.fdf = NULL;
    f.n = f0_size;
    f.p = ModelInt->Vc;
    f.params = ModelInt;

    //status = gsl_multifit_fdfsolver_dif_fdf( &x.vector, &f, fx, J );
    status = gsl_multifit_fdfsolver_dif_fdf_simple( &x.vector, &f, fx, J );


    gsl_vector_free( fx );

    return status;
}



int _int_det_eval( const gsl_vector* x, void* params, gsl_vector* f )
{

    int i, j;
    int status = 0;

    int_par* ModelInt = (int_par*) params;
    red_par* Model = ModelInt->Model;

    int f0_size = ModelInt->ns * ModelInt->Model->ny;
    double* f0 = (double*) malloc( sizeof(double) * f0_size );

    int nu = Model->nu;
    int ny = Model->ny;
    int nyny = Model->ny * Model->ny;
    int postD = (!Model->nonlin) ? ((!Model->backward) ? nyny+nu*nyny+nu*ny : 2*nyny+nu*nyny+nu*ny) :
				((!Model->backward) ? nyny+nu*nyny+nu*ny+ny*nyny : 2*nyny+nu*nyny+nu*ny+ny*nyny );
    int A_offs = (!Model->backward) ? 0 : nyny;
    int preC = (!Model->backward) ? nu+1 : nu+2;


    // transform back to full parameter space
    double* P = (double*) malloc( sizeof(double) * ModelInt->Vr );
    gsl_vector_view tmpP = gsl_vector_view_array( P, ModelInt->Vr );
    gsl_matrix_view tmpV = gsl_matrix_view_array( ModelInt->V, ModelInt->Vr, ModelInt->Vc );
    gsl_matrix_view tmpMx;
    double* C = (double*) malloc( sizeof(double) * nu*ny );

    gsl_blas_dgemv( CblasNoTrans, 1.0, &tmpV.matrix, x, 0.0, &tmpP.vector);

    // transpose A and B
    for ( i = 0; i < preC; ++i )
    {
	tmpMx = gsl_matrix_view_array( &P[i*nyny], ny, ny );
	gsl_matrix_transpose( &tmpMx.matrix );
    }
    // transpose C
    tmpMx = gsl_matrix_view_array( &P[nyny+A_offs + nu*nyny], nu, ny );
    for ( i = 0; i < ny; ++i )
    {
	for ( j = 0; j < nu; ++j )
	{
	    C[i*nu + j] = gsl_matrix_get( &tmpMx.matrix, j, i );
	}
    }
    // transpose D
    if ( Model->nonlin )
    {
	for ( i = 0; i < ny; ++i )
	{
	    tmpMx = gsl_matrix_view_array( &P[i*nyny + nyny+A_offs + nu*nyny + nu*ny], ny, ny );
	    gsl_matrix_transpose( &tmpMx.matrix );
	}
    }

    Model->A = &P[0];							// linear parameters
    Model->B = &P[nyny+A_offs];						// bilinear parameters
    Model->C = C;							// exogenous parameters
    Model->D = (!Model->nonlin) ? NULL : &P[nyny+A_offs+nu*nyny+nu*ny];	// nonlinear parameters
    Model->transit = &P[postD];						// transit time coefficient
    Model->decay = &P[postD+ny];					// signal decay coefficient
    Model->epsilon = P[postD+ny+ny];

/*
printf("\n");
for ( i = 0; i < 4; ++i )
{
    for ( j = 0; j < 3; ++j )
	printf("%f  ", Model->C[i*3+j] );
    printf("\n");
}
for ( j = 0; j < 4; ++j )
    printf("%f  ", Model->decay[j] );
printf("\n%f\n", Model->epsilon);

printf("\n");
for ( i = 0; i < 85; ++i )
{
    for ( j = 0; j < 23; ++j )
	printf("%f  ", gsl_matrix_get( &tmpV.matrix, i, j ) );
    printf("\n");
}
printf("\n");
*/


    if ( !Model->nonlin )
    {
	status = int_det( ModelInt, f0 );
    }
    else
    {
	status = int_det_D( ModelInt, f0 );
    }

/*
printf("\n");
for ( j = 0; j < 12; ++j )
    printf("%f  ", f0[j] );
printf("\n");
*/

    for ( i = 0; i < f0_size; ++i )
    {
	gsl_vector_set( f, i, f0[i] );
    }


    free ( C );
    free ( P );
    free ( f0 );

    return status;

}


int _compare_i ( const void* a, const void* b )
{
    return ( *(int*)a - *(int*)b );
}



