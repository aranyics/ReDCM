

ReDCM_bilinear_condition = function(M0, N, dt)
{
  
  # regulariser (1/8 of kernel support)
  #--------------------------------------------------------------------------
  t = 16/(N*dt);
  
  # remove unstable modes from Jacobian
  #--------------------------------------------------------------------------
  dfdx  = M0[2:dim(M0)[1], 2:dim(M0)[2]]
  us    = eigen( dfdx )
  s     = us$values
  u     = us$vectors
  i     = ReDCM_find( s )
  i     = i$r[which( Re(s) > -t )]
  s[i]  = sqrt( as.complex(-1) ) * Im(s[i]) - log( exp(t) + exp( Re(-s[i]) ) )
  
  # replace in bilinear operator
  #--------------------------------------------------------------------------
  M0[2:dim(M0)[1], 2:dim(M0)[2]] = Re( u %*% diag(s) %*% ReDCM_ginv(u) )
  
  return ( M0 )

}