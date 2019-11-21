#include "fx_fmri.h"


int _dfdxeval( const gsl_vector* x, void* params, gsl_vector* f );
int _feval( const gsl_vector* x, void* params, gsl_vector* f );




int fx_fmri ( fx_par* data, double* f0 )
{

    // declarations
    /* ---------------------------------------------------------------- */
    
    int status = 0;
    unsigned int i, j;
    int states = data->twostate ? 6 : 5;
    int ny = data->ny;
    int nu = data->nu;
    const double C_scale = 0.0625;
    const int mx_len = ny * ny;
    //double H[6] = {0.64, 0.32, 2.0, 0.32, 0.4, 0.0};
    double* H = data->H;

    double* x = (double*) malloc( sizeof(double) * ny * states );
    double* u = (double*) malloc( sizeof(double) * nu );
    double* A = (double*) malloc( sizeof(double) * mx_len );
    double* C = (double*) malloc( sizeof(double) * ny * nu );

    int c1 = 0;
    int c2 = ny*1;
    int c3 = ny*2;
    int c4 = ny*3;
    int c5 = ny*4;
    int c6 = ny*5;

    double *sd = (double*) malloc( sizeof(double) * ny );
    double *tt = (double*) malloc( sizeof(double) * ny );
    double *fv = (double*) malloc( sizeof(double) * ny );
    double *ff = (double*) malloc( sizeof(double) * ny );
    
    gsl_matrix_view EE = gsl_matrix_view_array( A, ny, ny );
    gsl_matrix_view Cm = gsl_matrix_view_array( C, ny, nu );
    gsl_vector* Cu     = gsl_vector_calloc( ny );

    gsl_matrix* IE = gsl_matrix_calloc( ny, ny );
    gsl_matrix* EI = gsl_matrix_alloc( ny, ny );
    gsl_matrix* SE = gsl_matrix_alloc( ny, ny );
    gsl_matrix* SI = gsl_matrix_alloc( ny, ny );
    gsl_vector_view EEdiag = gsl_matrix_diagonal( &EE.matrix );
    gsl_vector_view IEdiag = gsl_matrix_diagonal( IE );
    
    gsl_vector_view uv;
    gsl_vector_view x1;
    gsl_vector_view x6;

    // initialization
    /* ---------------------------------------------------------------- */
    //f0 = (double*) realloc( f0, sizeof(double) * ny*states );
    arr_cpyArray( ny*states, f0, data->x );
    arr_cpyArray( ny*states, x, data->x );
    arr_cpyArray( nu, u, data->u );
    arr_cpyArray( mx_len, A, data->A );
    arr_cpyArray( ny*nu, C, data->C );

    x1 = gsl_vector_view_array( data->x, ny );
    x6 = gsl_vector_view_array( &data->x[c6], ny );
    uv = gsl_vector_view_array( u, nu );

    data->backward = 0;
    data->A_eigen = 0;
    data->symmetry = 0;

    /* ================================================================ */
    /* (1) - Neuronal motion */
    /* ================================================================ */
    gsl_matrix_scale( &Cm.matrix, C_scale );
    
    // implement differential state equation y = dx/dt (neuronal)
    /* ---------------------------------------------------------------- */
    
    // Single-state DCM - Hidden states: 5 (4 hemodynamic + 1 neuronal)
    /* ================================================================ */
    if ( states == 5 )
    {

	// average connections are encoded explicitly
	/* ============================================================ */
	if ( !data->A_eigen )
	{
	
	    // input dependent modulation
	    /* -------------------------------------------------------- */
	    for ( i = 0; i < nu; ++i )
	    {
		arr_addArray( mx_len, A, &data->B[i*mx_len], u[i] );
	    }
	
	    // and nonlinear (state) terms
	    /* -------------------------------------------------------- */
	    if ( data->nonlin )
	    {
		for ( i = 0; i < ny; ++i )
		{
		    arr_addArray( mx_len, A, &data->D[i*mx_len], x[c1+i] );
		}
	    }
	
	    // one neuronal state per region: diag(A) is a self-inhibition
	    /* -------------------------------------------------------- */
	    for ( i = 0; i < ny; ++i )
	    {
		double diagAi = A[i*ny + i];
		A[i*ny + i] = diagAi - ( exp(diagAi) * 0.5 + diagAi );
	    }
	
	    // symmetry constraints for demonstration purposes
	    /* -------------------------------------------------------- */
	    if ( data->symmetry )
	    {
		gsl_matrix* tEE = gsl_matrix_alloc( ny, ny );
		gsl_matrix_transpose_memcpy( tEE, &EE.matrix );
		gsl_matrix_add( &EE.matrix, tEE );
		gsl_matrix_scale( &EE.matrix, 0.5 );
		gsl_matrix_free( tEE );
	    }
	
	}
	// otherwise A encodes the eigenvalues of the conn. matrix
	/* ============================================================ */
	else
	{
	    ;
	}
	
	// flow
	/* ------------------------------------------------------------ */
	// R: f[,1] == EE %*% x[,1] + P@C %*% u
	gsl_blas_dgemv( CblasNoTrans, 1.0, &Cm.matrix, &uv.vector, 0.0, Cu );
	gsl_blas_dgemv( CblasNoTrans, 1.0, &EE.matrix, &x1.vector, 1.0, Cu );
	//gsl_vector_add( Cu, &x1.vector );
	arr_cpyArray( ny, f0, Cu->data );
	gsl_vector_free( Cu );

    }
    // Two-state DCM - Hidden states: 6 (4 hemodynamic + 2 neuronal)
    /* ================================================================ */
    else
    {
	// input dependent modulation
	/* ------------------------------------------------------------ */
	for ( i = 0; i < nu; ++i )
	{
	    arr_addArray( mx_len, A, &data->B[i*mx_len], u[i] );
	}
	
	// and nonlinear (state) terms
	/* ------------------------------------------------------------ */
	if ( data->nonlin )
	{
	    for ( i = 0; i < ny; ++i )
	    {
		arr_addArray( mx_len, A, &data->D[i*mx_len], x[c1+i] );
	    }
	    // !! derivatives difference between spm_fx_fmri and spm_diff
	}
	
	// extrinsic - two neuronal state per region
	/* ------------------------------------------------------------ */
	arr_evalUnary( exp, mx_len, A, A, 1.0 );
	arr_mulScalar( mx_len, A, 0.125, A );		// 1/8
	gsl_vector_memcpy( &IEdiag.vector, &EEdiag.vector );
	arr_addArray( mx_len, A, IE->data, -1.0 );
	gsl_matrix_set_identity( EI );
	gsl_matrix_set_identity( SI );
	gsl_matrix_set_identity( SE );
	gsl_matrix_scale( SE, 0.5 );
	
	// excitatory proportion
	/* ------------------------------------------------------------ */
	arr_addArray( mx_len, A, SE->data, -1.0 );

	// motion - excitatory and inhibitory: f = dx/dt
	/* ------------------------------------------------------------ */
	// M.: f(:,1) = EE*x(:,1) - IE*x(:,6) + P.C*u(:);
	gsl_blas_dgemv( CblasNoTrans, 1.0, &Cm.matrix, &uv.vector, 0.0, Cu );
	gsl_blas_dgemv( CblasNoTrans, 1.0, &EE.matrix, &x1.vector, 1.0, Cu );
	gsl_blas_dgemv( CblasNoTrans, -1.0, IE, &x6.vector, 1.0, Cu );
	arr_cpyArray( ny, f0, Cu->data );

	// M.: f(:,6) = EI*x(:,1) - SI*x(:,6);
	gsl_blas_dgemv( CblasNoTrans, 1.0, EI, &x1.vector, 0.0, Cu );
	gsl_blas_dgemv( CblasNoTrans, -1.0, SI, &x6.vector, 1.0, Cu );
	arr_cpyArray( ny, &f0[c6], Cu->data );
	gsl_vector_free( Cu );


	gsl_matrix_free( IE );
	gsl_matrix_free( EI );
	gsl_matrix_free( SE );
	gsl_matrix_free( SI );
    }



    /* ======================================================================= */
    /* (2) - Hemodynamic motion */
    /* ======================================================================= */

    // hemodynamic parameters
    /* ----------------------------------------------------------------------- */
    // H(1) - signal decay                                   d(ds/dt)/ds)
    // H(2) - autoregulation                                 d(ds/dt)/df)
    // H(3) - transit time                                   (t0)
    // H(4) - exponent for Fout(v)                           (alpha)
    // H(5) - resting oxygen extraction                      (E0)
    // H(6) - ratio of intra- to extra-vascular components   (epsilon)
    //        of the gradient echo signal
    /* ----------------------------------------------------------------------- */
    
    // exponentiation of hemodynamic state variables
    /* ----------------------------------------------------------------------- */
    // R: x[,3:5] = exp( x[,3:5] )
    arr_expArray( 3*ny, &x[c3], 1.0, &x[c3] );
    
    // signal decay
    /* ----------------------------------------------------------------------- */
    // R: sd = H[1] * exp(P@decay)
    arr_expArray( ny, sd, H[0], data->decay );
    
    // transit time
    /* ----------------------------------------------------------------------- */
    // R: tt = H[3]* exp(P@transit)
    arr_expArray( ny, tt, H[2], data->transit );
    
    // Fout = f(v) - outflow
    /* ----------------------------------------------------------------------- */
    // R: fv = x[,4] ^ (1/H[4])
    arr_powScalar( ny, fv, 1/H[3], &x[c4] );
    
    // e = f(f) - oxygen extraction
    /* ----------------------------------------------------------------------- */
    // R: ff = ( 1 - (1-H[5])^(1/x[,3]) ) / H[5]	// signal decay
    for ( i = 0; i < ny; ++i )
    {
	ff[i] = (1 - pow(1.0-H[4], 1/x[c3+i])) / H[4];
    }

    // implement differential state equation y = dx/dt (hemodynamic)
    /* ----------------------------------------------------------------------- */
    for ( i = 0; i < ny; ++i )
    {
	// y(:,2) = x(:,1) - sd.*x(:,2) - H(2)*(x(:,3) - 1);
	f0[c2+i] = x[c1+i] - sd[i] * x[c2+i] - H[1] * (x[c3+i] - 1.0);
	// y(:,3) = x(:,2)./x(:,3);
	f0[c3+i] = x[c2+i] / x[c3+i];
	// y(:,4) = (x(:,3) - fv)./(tt.*x(:,4));
	f0[c4+i] = (x[c3+i] - fv[i]) / (tt[i] * x[c4+i]);
	// y(:,5) = (ff.*x(:,3) - fv.*x(:,5)./x(:,4))./(tt.*x(:,5));
	f0[c5+i] = (ff[i] * x[c3+i] - fv[i] * x[c5+i] / x[c4+i]) / (tt[i] * x[c5+i]);

	//printf("sd: %f, tt: %f, fv: %f, ff: %f\n", sd[i], tt[i], fv[i], ff[i]);
	//printf("c2: %f, c3: %f, c4: %f, c5: %f\n", f0[c2+i], f0[c3+i], f0[c4+i], f0[c5+i]);


    }


    free(sd);
    free(tt);
    free(fv);
    free(ff);

    free(A);
    free(C);
    free(x);
    free(u);


    return status;

}




int dfdx( fx_par* Model, gsl_matrix* J )
{
    //int i, j;
    int status = 0;

    int states = Model->twostate ? 6 : 5;
    int f0_size = Model->ny * states;

    gsl_vector* fx = gsl_vector_alloc( f0_size );
    gsl_vector_view x = gsl_vector_view_array( Model->x, f0_size );
    gsl_multifit_function_fdf f;
;

    f.f = &_feval;
    f.df = NULL;
    f.fdf = NULL;
    f.n = f0_size;
    f.p = f0_size;
    f.params = Model;

    //status = gsl_multifit_fdfsolver_dif_fdf( &x.vector, &f, fx, J );
    status = gsl_multifit_fdfsolver_dif_fdf_simple( &x.vector, &f, fx, J );


    gsl_vector_free( fx );

    return status;
}


int dfdxdu( fx_par* Model, gsl_matrix* J )
{
    //int i, j;
    int status = 0;

    int states = Model->twostate ? 6 : 5;
    int f0_size = Model->ny * states;

    gsl_vector* fx = gsl_vector_alloc( f0_size * f0_size );
    gsl_vector_view x = gsl_vector_view_array( Model->u, Model->nu );
    gsl_multifit_function_fdf f;


    f.f = &_dfdxeval;
    f.df = NULL;
    f.fdf = NULL;
    f.n = f0_size * f0_size;
    f.p = Model->nu;
    f.params = Model;

    //status = gsl_multifit_fdfsolver_dif_fdf( &x.vector, &f, fx, J );
    status = gsl_multifit_fdfsolver_dif_fdf_simple( &x.vector, &f, fx, J );


    gsl_vector_free( fx );

    return status;
}


int dfdxdx( fx_par* Model, gsl_matrix* J )
{
    //int i, j;
    int status = 0;

    int states = Model->twostate ? 6 : 5;
    int f0_size = Model->ny * states;

    gsl_vector* fx = gsl_vector_alloc( f0_size * f0_size );
    gsl_vector_view x = gsl_vector_view_array( Model->x, f0_size );
    gsl_multifit_function_fdf f;


    f.f = &_dfdxeval;
    f.df = NULL;
    f.fdf = NULL;
    f.n = f0_size * f0_size;
    f.p = f0_size;
    f.params = Model;

    //status = gsl_multifit_fdfsolver_dif_fdf( &x.vector, &f, fx, J );
    status = gsl_multifit_fdfsolver_dif_fdf_simple( &x.vector, &f, fx, J );


    gsl_vector_free( fx );

    return status;
}


int dfdu( fx_par* Model, gsl_matrix* J )
{

    //int i, j;
    int status = 0;

    int states = Model->twostate ? 6 : 5;
    int f0_size = Model->ny * states;

    gsl_vector* fx = gsl_vector_alloc( f0_size );
    gsl_vector_view x = gsl_vector_view_array( Model->u, Model->nu );
    gsl_multifit_function_fdf f;


    f.f = &_feval;
    f.df = NULL;
    f.fdf = NULL;
    f.n = f0_size;
    f.p = Model->nu;
    f.params = Model;


    //status = gsl_multifit_fdfsolver_dif_fdf( &x.vector, &f, fx, J );
    status = gsl_multifit_fdfsolver_dif_fdf_simple( &x.vector, &f, fx, J );


    gsl_vector_free( fx );

    return status;
}



int _feval( const gsl_vector* x, void* params, gsl_vector* f )
{

    int i;
    int status = 0;

    fx_par* Model = (fx_par*) params;

    int states = Model->twostate ? 6 : 5;
    int f0_size = Model->ny * states;

    double* f0 = (double*) malloc( sizeof(double) * f0_size );


    status = fx_fmri( Model, f0 );

    for ( i = 0; i < f0_size; ++i )
    {
	gsl_vector_set( f, i, f0[i] );
    }


    free ( f0 );

    return status;

}


int _dfdxeval( const gsl_vector* x, void* params, gsl_vector* f )
{

    int i, j;
    int status = 0;

    fx_par* Model = (fx_par*) params;

    int states = Model->twostate ? 6 : 5;
    int f0_size = Model->ny * states;

    gsl_matrix* J = gsl_matrix_calloc( f0_size, f0_size );


    status = dfdx( Model, J );

    for ( i = 0; i < f0_size; ++i )
    {
	for ( j = 0; j < f0_size; ++j )
	{
	    gsl_vector_set( f, i*f0_size + j, gsl_matrix_get( J, i, j ) );
	}
    }


    gsl_matrix_free( J );

    return status;

}

