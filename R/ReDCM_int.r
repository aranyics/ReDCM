

#################################################################################################
# source:       ReDCM_int.r
#------------------------------------------------------------------------------------------------
# type:         Function
#------------------------------------------------------------------------------------------------
# title:        Integrates a MIMO bilinear system dx/dt = f(x,u) = A*x + B*x*u + Cu + D
#------------------------------------------------------------------------------------------------
# description:  Integrates the bilinear approximation to the MIMO system described by
#
#                   dx/dt = f(x,u,P) = A*x + u*B*x + C*u + D
#                   y     = g(x,u,P) = L*x
#
#               Fast integrator that uses a bilinear approximation to the Jacobian  evaluated
#               using ReDCM_bireduce. This routine will also allow for sparse sampling of the 
#               solution and delays in observing outputs, specified by M@delays. 
#               It is used primarily for integrating fMRI models (see also spm_int_D)
#
#               Input:    P   - model parameters object
#                         M   - model specification object
#                         U   - external input object
#
#               Output:   y   - response
#------------------------------------------------------------------------------------------------
# based on:     spm_int.m
#               DCM12 (http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
#------------------------------------------------------------------------------------------------
# reference:    ---
#------------------------------------------------------------------------------------------------
# TODO:         - optimize index matrices
#################################################################################################


ReDCM_int = function(P, M, U, reduc=FALSE, V=NULL)
{
  
  if( reduc )
  {
    if ( is.null(V) || !is.vector(P) )
    {
      if ( class(P)[1] == "Params" )
      {
        reduc = FALSE
      }
      else
      {
        cat( "Error (ReDCM_int): unrecognised P argument\n" )
        return( rep(0, M@l) )
      }
    }
    else      # if V is transform matrix, and P is parameter vetor
    {
      if( dim(V)[2] != length(P) )
      {
        cat( "Error (ReDCM_int): transform matrix dimensions don't match with parameter vetrox length\n")
        return( rep(0, M@l) )
      }
      else
      {
        pv = V %*% P
        P = .Params()
        P = pUnvect( P, pv, M@m, M@l )
      }
    }
  }

  # get number of time samples (v) and microtime bins (u)
  #-------------------------------------------------------------------------
  dt = U@dt
  
  u = dim(U@u)[1]
  v = M@ns
  
  
  # get expansion point
  #-------------------------------------------------------------------------
  x = c(1, M@x)
  
  # specify fx and gx functions
  #-------------------------------------------------------------------------
  #
  
  # Bilinear approximation (1st order)
  #-------------------------------------------------------------------------
  M0.M1 = c_bireduce( M, P )
  #M0.M1 = ReDCM_bireduce(M, P)
  M0 = M0.M1[[1]]
  M1 = M0.M1[[2]]
  
  m = dim(M1)[3]                  # m inputs
  
  # delays        # only if M@delays > U@dt !!!
  #-------------------------------------------------------------------------
  D = round(M@delays/U@dt)      # round(0.5) = 0 !!!
  
  
  # Evaluation times (t) and indicator array for inputs (su) and output (sy)
  #=========================================================================
  
  # get times that the input changes
  #-------------------------------------------------------------------------
  tu = c( 1, 1+ReDCM_find(diff(U@u))[,1] )
  su = rep(0, u)
  su[tu] = 1
  
  # get times that the response is sampled
  #-------------------------------------------------------------------------
  s = ceiling( seq(0, v-1, 1) * u / v )
  sy = NULL
  #return(su)
  for ( j in 1:M@l )
  {
    i = s + D[j]
    tmp = rep(0, u)
    tmp[i] = seq(1,v,1)
    sy = rbind( sy, tmp )
  }
  
  # time in seconds
  #-------------------------------------------------------------------------
  t = sort(unique( c( tu, ceiling( which(sy != 0) / M@l ) ) ));   # TODO: optimize
  su = su[t]
  sy = sy[,t]
  dt = c(diff(t), 0) * U@dt
  
  
  # Integrates
  #-------------------------------------------------------------------------
  y = matrix( nrow=M@l, ncol=v )
  J = M0
  U.u = as.matrix(U@u)
  
  for ( i in 1:length(t) )
  {
    
    # input dependent changes in Jacobian
    #-----------------------------------------------------------------------
    if( su[i] )
    {
      u = U.u[t[i], ]
      J = M0
      for ( j in 1:m )
      {
        J = J + u[j] * M1[,,j]
      }
    }
    
    # output sampled
    #-----------------------------------------------------------------------
    if( any(sy[, i]) )
    {
      #q = ReDCM_gx_fmri(x[-1], u, P, M)
      q = c_gx_fmri(x[-1], u, P, M)
      j = ReDCM_find( sy[,i] )[,1]
      s = sy[j[1], i]
      y[j, s] = q[j]
    }
    
    # compute updated states x = expm(J*dt)*x
    #-----------------------------------------------------------------------
    x = expm::expm( J*dt[i] ) %*% x
    
    # check for convergence
    #-----------------------------------------------------------------------
    if ( Matrix::norm(x) > 1e6 )
    {
      break
    }
    
  }
  
  y = Re(y)
  
  
  return( t(y) )


}




ReDCM_hemodynamics1 = function( P, M, U, freq=0.20125, len=40.25 )
{

    bins = len/freq
  
    U@dt = freq
    U@u = Matrix(0, nrow=bins, ncol=dim(U@u)[2])
    U@u[2,1] = 1
    
    P@A = diag( diag( P@A ) )
    P@B = P@B * 0
    P@C = P@C * 0
    P@C[,1] = 1
    
    M@delays = rep(0, M@l)
    M@Ydt = U@dt
    M@ns = bins
    
    res = c_int_det_hemodyn( P, M, U )
    
#     f0 = array(0, dim=c(M@n, bins+1))
#     for ( i in 2:bins )
#     {
#       cat(i, " ", f0[,i], "\n")
#       f0[9:M@n,i] = exp( f0[9:M@n, i] )
#       f0[,i+1] = diag( c_dfdx(f0[,i], U@u[i,], P, M) )
#       #f0[,i+1] = colSums( c_dfdx(f0[,i], U@u[i,], P, M) )
#       #cat( diag( c_dfdx(f0[,i], U@u[i,], P, M) ) )
#     }
#     
#     return ( list( res[[1]], f0) )
    
    return ( res[[1]] )

}




ReDCM_hemodynamics = function( freq=0.20125, len=40.25, H=NULL, pTransit=0, pDecay=0, pEpsilon=0 )
{

  bins = len/freq
  
  U = .External()
  U@name = 'spike'
  U@dt = freq
  U@u = Matrix(0, nrow=bins, ncol=1)
  U@u[1,1] = 1
  
  P = .Params()
  P@A = array(0)
  P@B = array(0)
  P@C = array(1)
  P@D = array(dim=0)
  P@transit = array(pTransit)
  P@decay = array(pDecay)
  P@epsilon = array(pEpsilon)
  
  M = .Model()
  M@l = 1
  M@m = 1
  M@x = array(0, dim=c(1,5))
  M@n = 5
  M@delays = rep(0, M@l)
  M@Ydt = U@dt
  M@ns = bins
  if ( !is.null(H) && is.numeric(H) && length(H) == 5 )
  {
    M@H = H
  }
  # return ( list( P, M, U ) )
  
  res = c_int_det_hemodyn( P, M, U )
  
  cat("\n")
  cat("H1(",M@H[1],"), signal decay: H[1]*exp(decay) =                  ", M@H[1]*exp(pDecay), "\n")
  cat("H2(",M@H[2],"), autoregulation: H[2] =                           ", M@H[2], "\n")
  cat("H3(",M@H[3],"   ), transit time: H[3]*exp(transit) =                ", M@H[3]*exp(pTransit), "\n")
  cat("H4(",M@H[4],"), Grubb's exponent for outflow: H[4] =             ", M@H[4], "\n")
  cat("H5(",M@H[5]," ), resting oxygene extraction: H[5] =               ", M@H[5], "\n")
  cat("H6(",1,"   ), intra-extravascular signal ratio: exp(epsilon) = ", exp(pEpsilon), "\t(ratio of the intra- to extra-vascular components of the gradient echo signal)\n")
  cat("\n")
  
  return ( res[[1]] )
  
}
