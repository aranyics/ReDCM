

#################################################################################################
# source:       ReDCM_gx_fmri.r
#------------------------------------------------------------------------------------------------
# type:         Function
#------------------------------------------------------------------------------------------------
# title:        Simulated BOLD response to input
#------------------------------------------------------------------------------------------------
# description:  Returns a predicted response from state variables passed through the
#               output nonlinearity:
#
#                   g   = V0 * (k1(1 - q) + k2(1 - q/v) + k3(1 - v))
#
#               Input:    x     - state variables
#                         P     - model parameters object
#
#               Output:   g     - BOLD response
#                         [dgdx - dg/dx]
#------------------------------------------------------------------------------------------------
# based on:     spm_gx_fmri.m
#               DCM12 (http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
#------------------------------------------------------------------------------------------------
# reference:    BOLD signal model
#                   1.  Stephan KE, Weiskopf N, Drysdale PM, Robinson PA, Friston KJ (2007)
#                       Comparing hemodynamic models with DCM. NeuroImage 38: 387-401.
#------------------------------------------------------------------------------------------------
# TODO:         - return derivatives
#################################################################################################


ReDCM_gx_fmri = function (x, u, P, M=NULL, g.only=TRUE, ...)
{
  
  M = NULL
  
  if ( length(list(...)) > 0 )
  {
    for ( i in 1:length(list(...)) )
    {
      if ( class( list(...)[[i]] )[1] == "Model" )
      {
        M = list(...)[[i]]
      }
    }
  }
  
  if ( class(P)[1] != "Params" )
  {
    cat('Error: P is not object of class "Params"\n')
    return(0)
  }
  
  
  x       = matrix( x, ncol=5 )
  
  #==========================================================================
  # (1) - Biophysical constants for 1.5T
  #==========================================================================
    
  # time to echo (TE) (default 0.04 sec)
  #--------------------------------------------------------------------------
  TE  = 0.04
  
  # resting venous volume (%)
  #--------------------------------------------------------------------------
  V0  = 4

  # estimated region-specific ratios of intra- to extra-vascular signal 
  #--------------------------------------------------------------------------
  ep  = exp(P@epsilon)

  # slope r0 of intravascular relaxation rate R_iv as a function of oxygen 
  # saturation S:  R_iv = r0*[(1 - S)-(1 - S0)] (Hz)
  #--------------------------------------------------------------------------
  r0  = 25

  # frequency offset at the outer surface of magnetized vessels (Hz)
  #--------------------------------------------------------------------------
  nu0 = 40.3

  # resting oxygen extraction fraction
  #--------------------------------------------------------------------------
  E0  = 0.4

  #==========================================================================
  # (2) - Coefficients in BOLD signal model
  #==========================================================================
  k1  = 4.3 * nu0 * E0 * TE
  k2  = ep * r0 * E0 * TE
  k3  = 1 - ep

  #==========================================================================  
  # (3) - Output equation of BOLD signal model
  #==========================================================================
  v   = exp( x[,4] )
  q   = exp( x[,5] )
  g   = V0 * (k1 - k1*q + k2 - k2*q/v + k3 - k3*v);
  
  g   = c(g)
  
  if ( g.only )
  {
    return(g)
  }


  #==========================================================================
  # (4) - derivative dgdx
  #==========================================================================
  
  # TODO
  
  
  
  return(g)
  
}
