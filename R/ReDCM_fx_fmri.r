

#################################################################################################
# source:       ReDCM_fx_fmri.r
#------------------------------------------------------------------------------------------------
# type:         Function
#------------------------------------------------------------------------------------------------
# title:        State equation for a dynamic model of fMRI responses
#------------------------------------------------------------------------------------------------
# description:  Solves state equations at time point u for a dynamic [bilinear/nonlinear/Balloon]
#               model of fMRI responses
#
#               Input:    x     - state variables
#                                   x[,1]   - excitatory neuronal activity      ue
#                                   x[,2]   - vascular signal                   s
#                                   x[,3]   - rCBF                              ln(f)
#                                   x[,4]   - venous volume                     ln(v)
#                                   x[,5]   - deoxyHb                           ln(q)
#                                   [x[,6]  - inhibitory neuronal activity      ui]
#                         u     - external input time point
#                         P     - model parameters object
#                         M     - model specifiactions object
#
#               Output:   f     - dx/dt
#                         [dfdx - df/dx]
#                         [dfdu - df/du]
#                         [D    - delays]
#------------------------------------------------------------------------------------------------
# based on:     spm_fx_fmri.m
#               DCM12 (http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
#------------------------------------------------------------------------------------------------
# reference:    hemodynamic and neuronal state equation
#                   1.  Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
#                       changes during brain activation: The Balloon model. MRM 39:855-864,
#                       1998.
#                   2.  Friston KJ, Mechelli A, Turner R, Price CJ. Nonlinear responses in
#                       fMRI: the Balloon model, Volterra kernels, and other hemodynamics.
#                       Neuroimage 12:466-477, 2000.
#                   3.  Stephan KE, Kasper L, Harrison LM, Daunizeau J, den Ouden HE,
#                       Breakspear M, Friston KJ. Nonlinear dynamic causal models for fMRI.
#                       Neuroimage 42:649-662, 2008.
#                   4.  Marreiros AC, Kiebel SJ, Friston KJ. Dynamic causal modelling for
#                       fMRI: a two-state model.
#                       Neuroimage. 2008 Jan 1;39(1):269-78.
#------------------------------------------------------------------------------------------------
# TODO:         - implement 2-state branch
#               - implement nonlinear branch
#               - return derivatives
#################################################################################################


ReDCM_fx_fmri = function (x, u, P, M=NULL, f.only=TRUE, ...)
{

  M = NULL
  symmetry = 0
  
  if ( length(list(...)) == 0 )
  {
    symmetry = 0
  }
  else
  {
    # unargin = NULL
    for ( i in 1:length(list(...)) )
    {
      if ( class( list(...)[[i]] )[1] == "Model" )
      {
        M = list(...)[[i]]
      }
      # else
      # {
      #   unargin = c(unargin, (i+3))
      # }
    }
    # if ( !is.null(unargin) )
    # {
    #   cat('Unused arguments (ReDCM_fx_fmri): ', unargin, '\n')
    # }
  }
  
  if ( class(P)[1] != "Params" )
  {
    cat('Error: P is not object of class "Params"\n')
    return(0)
  }
  
  
  #==================================================================
  # (1) - Neuronal motion
  #==================================================================
  #P@A                                # linear parameters
  #P@B                                # bilinear parameters
  P@C     = P@C / 16                  # exogenous parameters
  #P@D                                # nonlinear parameters
  
  # implement differential state equation y = dx/dt (neuronal)
  #------------------------------------------------------------------
  x       = matrix( x, ncol=5 )
  f       = x
  
  # Single-state DCM - Hidden states: 5 (4 hemodynamic + 1 neuronal)
  #==================================================================
  if ( dim(x)[2] == 5 )
  {
    
    SE = NULL               # self-inhibition (self excitatory)
    EE = NULL               # external excitation
    
    # average connections are encoded explicitly
    #================================================================
    if ( length(dim(P@A)) == 2  &&  is.matrix(P@A) )
    {
      
      # input dependent modulation
      #--------------------------------------------------------------
      for ( i in 1:dim(P@B)[3] )
      {
        P@A = P@A + u[i] * P@B[,,i]
      }
      
      # and nonlinear (state) terms
      #--------------------------------------------------------------
      if ( dim(P@D)[3] > 0 )
      {
        for ( i in 1:dim(P@D)[3] )
        {
          P@A = P@A + x[i,1] * P@D[,,i]
        }
      }
      
      # combine forward and backward connections (if necessary)
      #--------------------------------------------------------------
      if ( !is.na(dim(P@A)[3]) && dim(P@A)[3] > 1 )
      {
        P@A = exp(P@A[,,1]) - exp(P@A[,,2])
      }
      
      # one neuronal state per region: diag(A) is a self-inhibition
      #--------------------------------------------------------------
      SE = diag(P@A)
      EE = P@A - diag( exp(SE)/2 + SE )
      
      # symmetry constraints for demonstration purposes
      #--------------------------------------------------------------
      if ( symmetry )
      {
        EE = (EE + t(EE)) / 2
      }
      
      #cat('Bingo\n')
    }
    # otherwise P@A encodes the eigenvalues of the connectivity matrix
    #================================================================
    else
    {
      # TODO
    }
    
    # flow
    #----------------------------------------------------------------
    f[,1] = EE %*% x[,1] + P@C %*% u
    
  }
  # Two-state DCM - Hidden states: 5 (4 hemodynamic + 2 neuronal)
  #==================================================================
  else
  {
    # TODO
  }
  
  
  #==========================================================================
  # (2) - Hemodynamic motion
  #==========================================================================
    
  # hemodynamic parameters
  #--------------------------------------------------------------------------
  # H(1) - signal decay                                   d(ds/dt)/ds)
  # H(2) - autoregulation                                 d(ds/dt)/df)
  # H(3) - transit time                                   (t0)
  # H(4) - exponent for Fout(v)                           (alpha)
  # H(5) - resting oxygen extraction                      (E0)
  # H(6) - ratio of intra- to extra-vascular components   (epsilon)
  #        of the gradient echo signal
  #--------------------------------------------------------------------------
  
  H         = c(0.64, 0.32, 2.0, 0.32, 0.4)
  if (!is.null(M))
  {
    H = M@H
  }
  
  # exponentiation of hemodynamic state variables
  #--------------------------------------------------------------------------
  x[,3:5]   = exp( x[,3:5] )
  
  # signal decay
  #--------------------------------------------------------------------------
  sd        = H[1] * exp(P@decay)
  
  # transit time
  #--------------------------------------------------------------------------
  tt        = H[3] * exp(P@transit)
  
  # Fout = f(v) - outflow
  #--------------------------------------------------------------------------
  fv        = x[,4] ^ (1/H[4])
  
  # e = f(f) - oxygen extraction
  #--------------------------------------------------------------------------
  ff        = ( 1 - (1 - H[5]) ^ (1 / x[,3]) ) / H[5]
  
  # implement differencial state equation f = dx/dt (hemodynamic)
  #--------------------------------------------------------------------------
  f[,2]     = x[,1] - sd * x[,2] - H[2] * (x[,3] - 1)
  f[,3]     = x[,2] / x[,3]
  f[,4]     = (x[,3] - fv) / (tt * x[,4])
  f[,5]     = (ff * x[,3] - fv * x[,5] / x[,4]) / (tt * x[,5])
  f         = c( f )
  
  if ( f.only )
  {
    return(f)
  }
  
  
  #==========================================================================
  # (3) - Compute Jacobians
  #==========================================================================
  
  # TODO
  
  
  
  return(EE)

}