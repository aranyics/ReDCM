

ReDCM_kernels = function( M0, M1, L1, L2, N, dt, o2=FALSE )
{

  
  # parameters
  #------------------------------------------------------------------------
  N = trunc(N)                        # kernel depth
  n = dim(M0)[1]                      # state variables
  m = dim(M1)[3]                      # inputs
  l = dim(L2)[3]                      # outputs
  X0 = array(0, dim=c(n,1))
  H1 = array(0, dim=c(N,n,m))
  K1 = array(0, dim=c(N,l,m))
  K2 = array(0, dim=c(N,N,l,m,m))


  # pre-compute matrix exponentials
  #------------------------------------------------------------------------
  e1 = Matrix::expm( dt * M0 )
  e2 = Matrix::expm(-dt * M0 )
  M = array( 0, dim=c(n,n,m,N) )
  for ( p in 1:m )
  {
    M[,,p,1] = array( e1 %*% M1[,,p] %*% e2 )
  }
  
  ei = e1
  for ( i in 2:N )
  {
    ei = ei %*% e1
    for ( p in 1:m )
    {
      M[,,p,i] = array( e1 %*% M[,,p,i-1] %*% e2 )
    }
  }
  
  
  # check for convergence and apply a more robust scheme if necessary
  #------------------------------------------------------------------------
  q = FALSE
  for ( p in 1:m )
  {
    q = q | base::norm( M[,,p,N], 'I' ) > base::norm( M[,,p,1], 'I' )
    q = q | base::norm( M[,,p,N], 'I' ) > exp(16)
    q = q | is.nan( base::norm( M[,,p,N], "2" ) ) | is.na( base::norm( M[,,p,N], "2" ) )
  }
  if ( q )
  {
    M0 = ReDCM_bilinear_condition(M0, N, dt)
    e1 = Matrix::expm( dt * M0 )
    
    ei = diag(n)
    for( i in 1:N )
    {
      ei = ei %*% e1
      ie = matrixcalc::svd.inverse( as.matrix(ei) )
      for ( p in 1:m )
      {
        M[,,p,i] = array( ei %*% M1[,,p] %*% ie )
      }
    }
  }
  
  
  # 0th order kernel
  #------------------------------------------------------------------------
  X0[1,1] = 1
  H0 = ei %*% X0
  K0 = L1 %*% H0
  
  
  # 1st order kernel
  #------------------------------------------------------------------------
  for ( p in 1:m )
  {
    for ( i in 1:N )
    {
      H1[i,,p] = M[,,p,i] %*% as.matrix(H0)
      K1[i,,p] = H1[i,,p] %*% t(L1)
    }
  }
  
  return ( list(H0, H1, K0, K1) )
  
  
  if ( o2 )
  {
    # 2nd order kernels
    #------------------------------------------------------------------------
    for ( p in 1:m )
    {
      for ( q in 1:m )
      {
        for ( j in 1:N )
        {
          H = array( apply( H1, 1, function(x){ L1 %*% M[,,q,j] %*% as.matrix(x[,p]) } ), dim=c(l,N) )
          for ( i in j:N )
          {
            K2[j,i,,q,p] = t(H)[i,]
            K2[i,j,,p,q] = t(H)[i,]
          }
        }
      }
    }
    
    # add output nonlinearity
    #------------------------------------------------------------------------
    for ( i in 1:m )
    {
      for ( j in 1:m )
      {
        for ( p in 1:l )
        {
          K2[,,p,i,j] = K2[,,p,i,j] + H1[,,i] %*% L2[,,p] %*% t(H1[,,j])
        }
      }
    }
    
  }
  
  return ( list(H0, H1, K0, K1, K2) )
  
}
