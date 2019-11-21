

ReDCM_generate = function ( Y, P, M, U, SNR=1, DCM=NULL )
{
  
  if ( !is.null(DCM) & class(DCM) == "Estimates" )
  {
    if ( is.na( DCM@Ep@A[1] ) )
    {
      cat('ReDCM_generate: DCM not estimated yet... ignoring argument DCM\n')
    }
    else
    {
      Y = DCM@Y
      P = DCM@Ep
      M = DCM@M[[1]]
      U = DCM@U
    }
  }

  v = M@ns
  n = M@l
  m = M@m
  
  # integrate and compute hemodynamic response at v sample points
  y = c_int_det( P, M, U )
  
  # compute standard deviation of additive noise for all areas
  r = diag( apply( Y@y, 2, sd ) / SNR )
  
  # add noise
  p = 1
  a = 1/16
  a = c( 1, -a )
  B = array(1, dim=c(v)) %*% t(a)
  K = base::solve( ReDCM_spdiags( B, -c(0:p), v, v ) )
  K = K * sqrt( v/sum(Matrix::diag(K %*% t(K))) )
  z = array( rnorm(v*n), dim=c(v, n) )
  e = K %*% z
  Y@Q = ReDCM_Ce( array(v, dim=c(n)) )
  
  Y@y = y[,1:n] + e %*% r

  return ( Y )

}
