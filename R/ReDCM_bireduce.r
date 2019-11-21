

#################################################################################################
# source:       ReDCM_bireduce.r
#------------------------------------------------------------------------------------------------
# type:         Function
#------------------------------------------------------------------------------------------------
# title:        Reduction of a fully nonlinear MIMO system to Bilinear form
#------------------------------------------------------------------------------------------------
# description:  Returns Matrix operators for the Bilinear approximation to the MIMO
#               system described by
#
#                   dx/dt = f(x,u,P)
#                   y(t)  = g(x,u,P)
#
#               A Bilinear approximation is returned, where the states are
#
#                   q(t)  = [1; x(t) - x(0)]
#
#               Input:    M   - model specification object
#                         P   - model parameters object
#
#               Output:   M0  - Bilinear operator - M0
#                         M1  - Bilinear operator - M1 = dM0/du
#                         [L1 - Output operators  - L1]
#                         [L2 - Output operators  - L2 = dL1/dx]
#------------------------------------------------------------------------------------------------
# based on:     spm_bireduce.m
#               DCM12 (http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
#------------------------------------------------------------------------------------------------
# reference:    ---
#------------------------------------------------------------------------------------------------
# TODO:         ---
#################################################################################################


ReDCM_bireduce = function( M, P, out1.op=FALSE, out2.op=FALSE )
{
  
  # set up
  #=========================================================================
  
  # expansion point
  #-------------------------------------------------------------------------
  x = M@x
  xv = c(x)
  u = rep(0, M@m)
  m = M@m
  
  
  # Partial derivatives for 1st order Bilinear operators
  #=========================================================================
  
  #f(x[0], 0) and derivatives
  #-------------------------------------------------------------------------
  # try wether M@dfdxu, M@dfdx, M@dfdu, M@f0 exists
  
#   f0 = ReDCM_fx_fmri(x, u, P, M)
#   n = length(f0)
#   
#   dfdx.fn = function(x, u, P, M, ...) {pjacob( ReDCM_fx_fmri, x, u, P, M, n=1 )}
#   dfdu.fn = function(x, u, P, M, ...) {pjacob( ReDCM_fx_fmri, x, u, P, M, n=2 )}
#   dfdxdu.fn = function(x, u, P, M, ...) {pjacob( dfdx.fn, x, u, P, M, n=2 )}
#   
#   dfdx = dfdx.fn( x, u, P, M )
#   dfdu = dfdu.fn( x, u, P, M )
#   dfdxdu = array( dfdxdu.fn( x, u, P, M ), dim=c(n,n,m) )

  f0 = c_fx_fmri(x, u, P, M)
  n = length(f0)

  dfdx = c_dfdx(x, u, P, M)
  dfdu = c_dfdu(x, u, P, M)
  dfdxdu = c_dfdxdu(x, u, P, M)

  # delay operator
  #-------------------------------------------------------------------------
  # try wether M@D exists
  
  
  # Bilinear operators
  #=========================================================================
  
  # Bilinear operator M0
  #-------------------------------------------------------------------------
  M0 = ReDCM_cat( 0, f0 - dfdx%*%xv, t(rep(0, dim(dfdx)[2])), dfdx, ncol=2 )
  
  # Bilinear operators - M1 = dM0/du
  #-------------------------------------------------------------------------
  M1 = NULL
  for ( i in 1:m )
  {
    M1 = c( M1, c(ReDCM_cat( 0, dfdu[,i] - dfdxdu[,,i]%*%xv, t(rep(0, dim(dfdxdu[,,i])[2])), dfdxdu[,,i], ncol=2 )) )
  }
  M1 = array( M1, dim=c(1+dim(dfdu)[1], 1+dim(dfdxdu)[2], m) )

  
  if ( !out1.op )
  {
    return( list(M0, M1) )
  }
  
  
  # Output operators ## UNTESTED ##
  #=========================================================================
  
  # g(x(0), 0)
  #-------------------------------------------------------------------------
  
  g0 = ReDCM_gx_fmri(x, u, P, M)
  
  dgdx.fn = function(x, u, P, M, ...) {pjacob( ReDCM_gx_fmri, x, u, P, M, n=1 )}
  
  dgdx = dgdx.fn(x, u, P, M)
  
  l = length(g0)
  
  # Output matrices - L1
  #-------------------------------------------------------------------------
  L1 = ReDCM_cat( g0 - dgdx%*%xv, dgdx, ncol=2 )
  
  
  if ( !out2.op )
  {
    return( list(M0, M1, L1) )
  }
  
  
  # Output matrices - L2
  #-------------------------------------------------------------------------
  
  dgdxdx.fn = function(x, u, P, M, ...) {pjacob( dgdx.fn, x, u, P, M, n=1 )}
  
  dgdxdx = dgdxdx.fn(x, u, P, M)
  dgdxdx = array( dgdxdx, dim=c( dim(dgdx), length(x)) )
  
  D = array(0, dim=c(n,length(dgdxdx[1,,1]),l))
  L2 = NULL
  for ( i in 1:l )
  {
    for ( j in 1:n )
    {
      D[j,,i] = dgdxdx[i,,j]
    }
    
    L2 = c( L2, ReDCM_cat( 0, rep(0, dim(D[,,i])[1]), t(rep(0, dim(D[,,i])[2])), D[,,i], ncol=2 ) )
  }
  L2 = array( L2, dim=c(1+dim(D)[1], 1+dim(D)[2], l) )
  
  
  return( list(M0, M1, L1, L2) )
  
  
  
}