

# INCLUDE
#################################################################################################
#source('class/ReDCM_generics.R')
#source('class/ReDCM_data.R')
#source('class/ReDCM_external.R')
#source('class/ReDCM_parameters.R')
#source('class/ReDCM_model.R')
#source('class/ReDCM_estimate_struct.R')

#source('ReDCM_nlsi_GN.r')
#source('ReDCM_fx_fmri.r')
#source('ReDCM_gx_fmri.r')
#source('ReDCM_bireduce.r')
#source('ReDCM_int.r')
#source('ReDCM_c.r')
#source('ReDCM_kernels.r')
#source('ReDCM_bilinear_condition.r')
#source('ReDCM_generate.r')
#source('ReDCM_evidence.r')
#source('ReDCM_BMS.r')



# Creates basis functions for Discrete Cosine Transform.
#-----------------------------------------------------------------------------

ReDCM_dctmtx = function( N, K=0 )
{

  if ( K == 0 )
  {
    K = N
  }

  n = matrix( seq(0, N-1, 1) )
  
  C = matrix( nrow=dim(n)[1], ncol=K, data=0)
  
  
  C[,1] = rep( 1, dim(n)[1] ) / sqrt(N)
  for (k in 2:K)
  {
    C[,k] = sqrt(2/N) * cos( pi * (2*n + 1) * (k - 1) / (2*N) )
  }
  
  return (C)

}




# return error covariance constraints (for serially correlated data)
#-----------------------------------------------------------------------------

ReDCM_Ce = function ( v, a=NULL )
{
  
  C = list()
  l = length(v)
  n = sum(v)
  k = 0
  
  if ( l > 1 )
  {
    for (i in 1:l)
    {
      dCda = ReDCM_Ce( v[i] )
      for (j in 1:length(dCda))
      {        
        s = ReDCM_find( dCda[[j]] )
        
        s$r = s$r + k
        s$c = s$c + k
        
        C = append( C, list( sparseMatrix( i=s$r, j=s$c, x=s$v, dims=c(n, n) ) ) )
      }
      k = v[i] + k
    }
  }
  else
  {
  C = list( Matrix(diag(rep(1,v))) )
  }
  
  
  return (C)
  
}



# returns non-zero elements of a matrix
#-----------------------------------------------------------------------------

ReDCM_find = function( mx, rows.first=TRUE )
{
  
  s = NULL
  
  if ( is.null(dim(mx)) )
  {
    mx = as.matrix(mx)
  }
  
  if ( rows.first )
  {
    for ( i in 1:dim(mx)[2] )
    {
      if ( is.null(s) )
      {
        s = data.frame( r=which( mx[,i] != 0 ), c=rep( i, length(which( mx[,i] != 0 )) ), v=mx[which( mx[,i] != 0 ),i] )
      }
      else
      {
        s = rbind( s, data.frame( r=which( mx[,i] != 0 ), c=rep( i, length(which( mx[,i] != 0 )) ), v=mx[which( mx[,i] != 0 ),i] ) )
      }
    }
  }
  else
  {
#     for ( i in 1:dim(mx)[1] )
#     {
#       if ( is.null(s) )
#       {
#         s = data.frame( r=which( mx[,i] != 0 ), c=rep( i, length(which( mx[,i] != 0 )) ), v=mx[which( mx[,i] != 0 ),i] )
#       }
#       else
#       {
#         s = rbind( s, data.frame( r=which( mx[,i] != 0 ), c=rep( i, length(which( mx[,i] != 0 )) ), v=mx[which( mx[,i] != 0 ),i] ) )
#       }
#     }
  }
  
  return(s)
}


ReDCM_find2 = function( mx, cols.first=TRUE )
{
  
  # if ( is.null(dim(mx)) )
  # {
  #   mx = as.matrix(mx)
  # }
  mx = as.matrix(mx)
  m = dim(mx)[1]
  n = dim(mx)[2]
  N = length( which(mx != 0) )
  s = data.frame( r=rep(0, N), c=rep(0, N), v=rep(0, N) )
  
  if ( cols.first )
  {
    k = 1
    j = which( mx != 0 )
    for ( i in j )
    {
      s[k, 1] = i %% m
      if ( s[k, 1] == 0 )
      {
        s[k, 1] = m
      }
      s[k, 2] = ceiling( i / m )
      s[k, 3] = mx[s[k, 1], s[k, 2]]
      k = k + 1
    }
  }
  else
  {
    mx = t(mx)
    m = dim(mx)[1]
    n = dim(mx)[2]
    
    k = 1
    j = which( mx != 0 )
    for ( i in j )
    {
      s[k, 1] = i - floor( i / m ) * m
      if ( s[k, 1] == 0 )
      {
        s[k, 1] = m
      }
      s[k, 2] = ceiling( i / m )
      s[k, 3] = mx[s[k, 1], s[k, 2]]
      k = k + 1
    }
    tmp = s[,1]
    s[,1] = s[,2]
    s[,2] = tmp
  }
  
  return(s)
}


ReDCM_find3 = function(mx)
{
  mx = as.matrix(mx)
  s = which( as.matrix(mx)!=0, arr.ind=TRUE )
  
  return(s)
}


## usage
# 
# f2= function( x, y, a, b ) { a*x^2 + b*y }
# A=2
# B=-2
# x=0:5
# y=10:15
# f2(x,y,A,B)
# pgrad( f2, x, y, 2, A, B )
# pgrad( f2, x, y, 1, A, B )
# pgrad2( f2, x, y,  A, B )


# df/dx1 (n=1) or df/dx2 (n=2) 
pgrad = function( fun, x1, x2, n, ... )
{
  if ( n ==1 )
  {
    f = function(y) { fun(y,x2, ... ) }
    x = x1
  }
  else
  {
    f = function(y) { fun(x1,y, ... ) }
    x = x2
  }
  
  return( numDeriv::grad( func = f, x=x, method="Richardson", side=NULL, method.args=list() ) )
}


# df/dx1 (n=1) or df/dx2 (n=2) 
pjacob = function( fun, x1, x2, n, ... )
{
  if ( n ==1 )
  {
    f = function(y) { fun(y,x2, ... ) }
    x = x1
  }
  else
  {
    f = function(y) { fun(x1,y, ... ) }
    x = x2
  }
  
  return( numDeriv::jacobian( func = f, x=x, method="Richardson", side=NULL, method.args=list() ) )
  #return( numDeriv::jacobian( func = f, x=x, method="simple", side=NULL, method.args=list() ) )
  #return( numDeriv::jacobian( func = f, x=x, method="complex", side=NULL, method.args=list() ) )
}



ReDCM_cat = function ( ..., ncol=1 )
{
  
  c = NULL
  l = list(...)
  nrow = length( l ) / ncol
  
  if ( ncol == 1 && nrow == 1 )
  {
    c = l[[1]]
  }
  else if ( ncol < 1 || nrow < 1 )
  {
    return( c )
  }
  else
  {
    
    if ( length( l ) %% ncol != 0 )
    {
      cat("Error (ReDCM_cat): Inconsistent number of inputs\n")
      return( c )
    }
    else
    {
      for( col in 1:ncol )
      {
        rtmp = NULL
        for ( row in 1:nrow )
        {
          idx = (col-1)*nrow + row
          if ( !is.null(dim(l[[idx]])) || length(l[[idx]]) == 1 )
          {
            rtmp = rbind( rtmp, l[[idx]] )
          }
          else 
          {
            rtmp = rbind( rtmp, matrix(l[[idx]]) )
          }
        }
        
        if ( length(rtmp) == 1 )
        {
          c = cbind( c, rtmp[1] )
        }
        else
        {
          c = cbind( c, rtmp )
        }
      }
    }
  }

  return( c )

}


ReDCM_trace = function(A, B=NULL)
{

  if ( is.null(B) )
  {
    return( sum(A) )
  }
  else
  {
    return( sum( t(A) %*% B ) )
  }

}


ReDCM_dx = function(dfdx, f, t=Inf)
{
  
  # TODO: only works if dfdx is matrix, and f is vector
  
  f = c(f)

  # t is a regularizer
  #------------------------------------------------------------------
  if ( is.finite( sum(t) ) )
  {
    t = Re( exp( t - log( complex(real=diag( -dfdx )) ) ) )
  }
  
  # use a [pseudo]inverse if all t > TOL
  #==================================================================
  if ( min(t) > exp(16) )
  {
  
    dx = -matrixcalc::svd.inverse(dfdx) %*% f
  
  }
  else
  {
    # ensure t is a scalar or matrix
    #----------------------------------------------------------------
    if ( is.vector(t) )
    {
      if ( length(t) == 1 )
      {
        t = matrix(t)
      }
      else
      {
        t = diag(t)
      }
    }
    
    # augment Jacobian and take matrix exponential
    #==================================================================
    J = ReDCM_cat(0, t %*% f, t(rep(0,dim(dfdx)[2])), t %*% dfdx, ncol=2)
    J = as( Matrix(J), "dgeMatrix" )
    dx = Matrix::expm(J)[-1, 1]
  
  }
  
  return( Re(dx) )
  
}


ReDCM_ginv = function(X, tol = sqrt(.Machine$double.eps))
{
  ## Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
    dimnames = dnx[2:1])
}


ReDCM_spdiags = function( B, d, m, n, extract=0 )
{

  if ( extract )
  {
    ;
  }
  else
  {
    moda = 0
    p = length(d)
    A = array(0, dim=c(m, n))
    
    if ( moda )
    {
      ;
    }
    else
    {
      len = array(0, dim=c(p+1))
      for ( k in 1:p )
      {
        len[k+1] = len[k] + length( c( max(1, 1-d[k]):min(m, n-d[k]) ) )
      }
      a = array(0, dim=c(2*max(m,n)-1, 3))
      for ( k in 1:p )
      {
        i = t( max(1, 1-d[k]):min(m, n-d[k]) )
        a[(len[k]+1):(len[k+1]), ] = c( i, i+d[k], B[i+(m>=n)*d[k], k] )
      }
    }
    
    Aout = array(0, dim=c(m, n))
    Aout[a[,1:2]] = a[,3]
    
    res = Aout
  }

  return ( res )

}


ReDCM_Ncdf = function( x, u=0, v=1 )
{
  
  # if ( any( v <= 0 ) )
  # {
  #   cat('ReDCM_Ncdf: variance is restricted to positive values\n')
  #   return(-1)
  # }
  if ( is.null(dim(x)) | length(dim(x)) < 2 )
  {
    x = array(x, dim=c(length(x), 1))
  }
  if ( is.null(dim(u)) | is.null(dim(v)) | length(dim(u)) < 2 | length(dim(v)) < 2 )
  {
    u = array(u, dim=c(length(u), 1))
    v = array(v, dim=c(length(v), 1))
  }
  if ( (dim(x)[1] != dim(u)[1] | dim(x)[1] != dim(v)[1]) & (dim(u)[1] != 1 | dim(v)[1] != 1) )
  {
    cat('ReDCM_Ncdf: arguments don\'t match in size\n ')
    return(-1)
  }
  
  ad = c( max(2, length(dim(x))), max(2, length(dim(u))), max(2, length(dim(v))) )
  rd = max(ad)
  as = array( c( dim(x), rep(1, rd-ad[1]), dim(u), rep(1, rd-ad[2]), dim(v), rep(1, rd-ad[3]) ), dim=c(3, rd) )
  rs = as[which.max(as), ]
  xa = (apply(as, 1, prod) > 1)*1
  
  # init results
  # -------------------------------------------------
  F0 = array(0, dim=rs)
  Q = seq( 1, dim(x)[1], 1 )
  
  if ( dim(x)[1] == dim(u)[1] )
  {
    md = (array(1, dim=dim(x)) & array(1, dim=dim(u)) & v > 0) * 1
    Q = ReDCM_find(md)[,1]
  }
  # else
  #   md = apply( array(1, dim=dim(x)), 1, function(x){x * 1} )
  
  if ( is.null(Q) )
  {
    cat('ReDCM_Ncdf: Q is empty\n')
    return(-1)
  }
  if ( xa[1] ) Qx = Q
  else Qx = 1
  if ( xa[2] ) Qu = Q
  else Qu = 1
  if ( xa[3] ) Qv = Q
  else Qv = 1
  
  # compute
  # -------------------------------------------------
  F0[Q] = 0.5 + 0.5*ReDCM_erf( (x[Qx]-u[Qu]) / sqrt(2*v[Qv]) )
  F0[-Q] = NaN
  
  return ( F0 )
  
}


ReDCM_inv = function( A, TOL=NULL )
{
  m = dim(A)[1]
  n = dim(A)[2]
  
  if ( is.null(TOL) )
  {
    TOL = max( pracma::eps( Matrix::norm(A, 'i') ) * max(m,n), exp(-32) )
  }
  
  eye = matrix(0, m, n)
  diag(eye) = 1
  
  return( base::solve( A + eye * TOL ) )
}


ReDCM_matrix2binary = function( x, m )
{

  l = m*m - m
  A = array( 0, l )
  M = diag(m)
  
  for ( i in 1:l )
  {
    A[i] = x %% 2; x = floor( x / 2 );
  }
  A = rev(A)
  
  M[which(M!=1)] = A
  
  return( M )
  
}

ReDCM_dec2binary = function( x )
{

  a = NULL
  
  while ( x>0 )
  {
    a = c( x %% 2, a); x = floor( x / 2 );
  }
  
  return( a )

}

ReDCM_binary2matrix = function( x )
{

  l = length(x)
  M = diag(sqrt(l))
  x = as.logical(x)*1
  x = rev( x )
  x = x[-which(M==1)]
  k=0
  
  for ( i in 1:length(x) )
  {
    k = k + x[i] * 2^(i-1)
  }
  
  return( k )
  
}


## error functions
ReDCM_erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
ReDCM_erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
ReDCM_erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)
ReDCM_erfcinv <- function (x) qnorm(x/2, lower = FALSE)/sqrt(2)


## estimation test
ReDCM_test = function( )
{
  return( ReDCM_estimate('data/DCM_att_3final.mat') )
}
