

ReDCM_nlsi_GN = function(M, U, Y)
{
  
  # TODO: try FS (feature selection) after integrator
  # TODO: check existence of required fields
  # TODO: spm_svd
  # TODO: vectorization
  # TODO: try revert
  
  # Read objects
  #=========================================================================  
    
  # response information
  #-------------------------------------------------------------------------
  y = Y@y                                 # response
  ns = M@ns                               # number of scans
  ny = length(y)                          # total number of response variables
  nr = ny / ns                            # number of response components
  
  dt = Y@dt                               # time-step (TR)
  
  Q = Y@Q                                 # precision components Q (list)
  nh = length(Q)                          # number of precision components
  nq = ny / dim(Q[[1]])[2]                #
  
  nb = dim(Y@X0)[1]                       # number of bins (confounds)
  nx = nr * ns / nb                       # number of blocks
  dfdu = Matrix::kronecker( diag(1, nx), Y@X0 )
  
  # input
  #-------------------------------------------------------------------------
  #u = U@u
  ne = dim(U@u)[2]                        # number of external inputs
  
  # save initial model
  #-------------------------------------------------------------------------
  M.init = M
  
  pE = M@pE                               # prior moments (parameter exceptations)
  pC = M@pC                               # prior moments (parameter covariance)
  
  hE = M@hE#[1:nh]                        # hyperpriors (expectations)
  h = hE
  #ihC = base::solve(M@hC)                 # hyperpriors (covariance)
  ihC = ReDCM_inv(M@hC)
  
  
  # dimension reduction of parameter space
  #-------------------------------------------------------------------------
  nu = dim(dfdu)[2]                       # number of parameters (confounds)
  np = length( which(colMeans(pC) != 0) ) # number of parameters (effective)
  iu = t(c(1:nu)) + np
  ip = t(c(1:np))
  V = svd(pC)$u[,1:np]         # spm_svd??? 
  
  # second order moments (in reduced space)
  #-------------------------------------------------------------------------
  pC = t(V) %*% pC %*% V
  uC = diag( rep(1, nu) * 1e8  )
  ipC = base::solve( ReDCM_cat( pC, matrix(rep(0,nu*np), ncol=np), matrix(rep(0,np*nu), ncol=nu), uC, ncol=2 ) )
  
  # initialize conditional density
  #-------------------------------------------------------------------------
  Eu = matrixcalc::svd.inverse(dfdu) %*% c(y)
  pEv = pVect(pE)
  p = array( c( t(V) %*% (pVect(M@pE) - pEv), Eu ), dim=c(np+nu, 1) )
  Ep = .Params()
  Ep = pUnvect( Ep, pEv + V %*% p[ip], ne, nr )
  
  # unlist Q
  #-------------------------------------------------------------------------
  # Q = array( unlist( lapply(Q, function(x) { Matrix::as.matrix(x) }) ), dim=c(ny,ny,nr) )
  
  
  # EM
  #=========================================================================
  criterion = c(FALSE, FALSE, FALSE, FALSE)
  
  F0 = NULL
  C.F = -Inf                                    # free energy
  C.p = NULL
  C.h = NULL
  C.Cp = NULL
  
  v = -4                                        # log ascent rate
  dFdh = array(0, dim=c(nh,1))
  dFdhh = array(0, dim=c(nh,nh))
  
  for (k in 1:128)    # or M@N
  {
    #time1 = system.time({ 
    # TODO: calculate time spent in the loop
    # tStart
    
    revert = FALSE
    
    # E-Step: prediction f, and gradients; dfdp
    #=======================================================================
    
    #try
    
      # gradients
      #-----------------------------------------------------------------------
      Ep.rv = c( t(V) %*% pVect(Ep) )
      
      #f = ReDCM_int( Ep, M, U )
      f = c_int_det( Ep, M, U )
    
      #dfdp.fn = function(pv, M, U, ...) {pjacob(ReDCM_int, pv, M, U, reduc=TRUE, V=V, n=1)}
      #dfdp = array( dfdp.fn( Ep.rv, M, U ), dim=c(ny, np) )
    
      dfdp = c_int_dfdp( Ep.rv, M, U, V )
      
      # check for stability
      #-----------------------------------------------------------------------
      normdfdp = norm( dfdp, 'i' )
      #revert = is.nan(normdfdp) || normdfdp > exp(32)
      revert = is.na(normdfdp) || normdfdp > exp(32)
      cat("revert: ", revert, "\n")
    
    #catch revert = true

    
    # if revert then
    # - increase regularization and try again
    if ( revert && k > 1 )
    {
      for ( i in 1:4 )
      {
        # reset expansion point and increase regularization
        #---------------------------------------------------------------------
        v = min(v-2, -4)
        
        # E-Step: update
        #---------------------------------------------------------------------
        p = C.p + ReDCM_dx( dFdpp, dFdp, v )
        Ep = pUnvect( Ep, pEv + V %*% p[ip], ne, nr )
        Ep.rv = c( t(V) %*% pVect(Ep) )
        
        f = c_int_det( Ep, M, U )
        dfdp = c_int_dfdp( Ep.rv, M, U, V )
        
        # check for stability
        #-----------------------------------------------------------------------
        normdfdp = norm( dfdp, 'i' )
        revert = is.na(normdfdp) || normdfdp > exp(32)
        cat("revert: ", revert, "\n")
        
        if ( revert == FALSE )
        {
          break;
        }
      }
    }
    # if still revert then
    # - convergence failure
    if ( revert )
    {
      cat( "nlsi_GN: convergence failure" )
      return ( 0 )
    }
    
    # prediction error and full gradients
    #-----------------------------------------------------------------------
    
    e = Matrix( c(y) - c(f) - dfdu %*% p[iu] )
    J = Matrix( - cbind(dfdp, dfdu) )
    

    #Rprof()
    # M-step; Fischer scoring scheme to find h = max( F(p,h) )
    #=======================================================================
    #P = simplify2array(Q)

    for ( m in 1:8 )
    {
      # TODO: care for complex numbers
      
      # precision and conditional covariance
      #---------------------------------------------------------------------
      iS = Q[[1]] * (exp(-32) + exp(h[1]))
      for ( i in 2:nh )
      {
        iS = iS + Q[[i]] * (exp(-32) + exp(h[i]))
      }
      S = Matrix::chol2inv( Matrix::chol( iS ) )
      #iS = Matrix::as.matrix(iS)
      #S = base::solve( iS )
      #S = Matrix( S )
      iS = Matrix::kronecker( Matrix::diag(1, nq), iS )

      Pp = t(J) %*% iS %*% J
      # Pp = base::Re( t(J) %*% iS %*% J )
      # Pp = base::Re( c_blas_dgemm( c_blas_dgemm( t(J), iS ), J ) )
      # Cp = base::solve( Pp + ipC )
      Cp = ReDCM_inv( Pp + ipC )
    
      # precision operators for M-step
      #---------------------------------------------------------------------
      P = Q
      PS = list()
      JPJ = list()
      for ( i in 1:nh )
      {
        P[[i]] = P[[i]] * exp(h[i])
        PS[[i]] = P[[i]] %*% S
        #P[[i]] = Matrix::as.matrix(P[[i]])
        # PS[[i]] = c_blas_dgemm( P[[i]], S )
        P[[i]] = Matrix::kronecker( Matrix::diag(1,nq), P[[i]] )
        JPJ[[i]] = t(J) %*% P[[i]] %*% J
        # JPJ[[i]] = base::Re( t(J) %*% P[[i]] %*% J )
        # JPJ[[i]] = base::Re( c_blas_dgemm( c_blas_dgemm( t(J), P[[i]] ), J ) )
      }
      
      # derivatives: dLdh = dL/dh, ...
      #---------------------------------------------------------------------
      for( i in 1:nh )
      {
        #dFdh[i] = ReDCM_trace( PS[[i]] ) * nq / 2 -
        #          base::Re( t(e) %*% Matrix::as.matrix(P[[i]]) %*% e ) / 2 -
        #          ReDCM_trace( Cp %*% JPJ[[i]] ) / 2
        #dFdh[i] = matrixcalc::matrix.trace( PS[[i]] ) * nq / 2 -
        #          ( t(e) %*% P[[i]] %*% e ) / 2 -
        #          matrixcalc::matrix.trace( t(Cp) %*% JPJ[[i]] ) / 2
        dFdh[i] = sum( Matrix::diag( PS[[i]] ) ) * nq / 2 -
                  as.numeric( t(e) %*% P[[i]] %*% e ) / 2 -
                  sum( Matrix::diag( t(Cp) %*% JPJ[[i]] ) ) / 2
        #dFdh[i] = matrixcalc::matrix.trace( PS[[i]] ) * nq / 2 -
        #          base::Re( c_blas_dgemm( c_blas_dgemm( t(e), Matrix::as.matrix(P[[i]]) ), e ) ) / 2 -
        #          matrixcalc::matrix.trace( c_blas_dgemm( t(Cp), JPJ[[i]] ) ) / 2
        for ( j in i:nh )
        {
          #dFdhh[i, j] = - ReDCM_trace( PS[[i]] %*% PS[[j]] ) * nq / 2
          #dFdhh[i, j] = - matrixcalc::matrix.trace( PS[[i]] %*% PS[[j]] ) * nq / 2
          dFdhh[i, j] = - sum( Matrix::diag( PS[[i]] %*% PS[[j]] ) ) * nq / 2
          #dFdhh[i, j] = - matrixcalc::matrix.trace( c_blas_dgemm( t( PS[[i]] ), PS[[j]] ) ) * nq / 2
          dFdhh[j, i] = dFdhh[i, j]
        }
      }
      
      # add hyperpriors
      #---------------------------------------------------------------------
      d = array( h - hE, dim=c(nh, 1) )
      dFdh = dFdh - ihC %*% d
      # dFdh = dFdh - c_blas_dgemm( ihC, d )
      dFdhh = dFdhh - ihC
      #Ch = base::solve( -dFdhh )
      Ch = ReDCM_inv( -dFdhh )
      
      # update ReML estimate
      #---------------------------------------------------------------------
      dh = ReDCM_dx( dFdhh, dFdh, 4 )
      dh = aaply( aaply( dh, 1, function(x){max(x,-1)}), 1, function(x){min(x,1)} )
      h = h + dh
      
      # convergence
      #---------------------------------------------------------------------
      dF = t(dFdh) %*% dh
      # dF = c_blas_dgemm( t(dFdh), dh )
      if ( dF < 1e-2 )
      {
        break
      }
      
    }
    cat('\n')

    # Rprof(NULL)
    # return (summaryRprof()$by.total)

    # E-step with Levensberg-Marquardt regularization
    #=======================================================================
    
    # objective function: F(p) = log_evidence - divergence
    #-----------------------------------------------------------------------
    Fe = - base::Re( (t(e) %*% iS %*% e)[1] ) / 2 -
          t(p) %*% ipC %*% p / 2 -
          t(d) %*% ihC %*% d / 2 -
          ny * base::log(8 * base::atan(1)) / 2 -
          Matrix::determinant( S, logarithm=TRUE )$modulus * nq / 2 +
          Matrix::determinant( ipC %*% Cp, logarithm=TRUE )$modulus / 2 +
          Matrix::determinant( ihC %*% Ch, logarithm=TRUE )$modulus / 2
    
    # record increases and reference log-evidence for reporting
    #-----------------------------------------------------------------------
    if ( k == 1 )
    {
      F0 = Fe
    }
    else
    {
      # TODO print Free-energy difference
    }
    
    # if F has increased, update gradients and curvatures for E-step
    #-----------------------------------------------------------------------
    if ( k < 4 || Fe > C.F )
    {
      
      # accept current estimates
      #---------------------------------------------------------------------
      C.p = p
      C.h = h
      C.F = Fe
      C.Cp = Cp
      
      # E-step: Conditional update of gradients and curvature
      #---------------------------------------------------------------------
      dFdp = - base::Re( Matrix::as.matrix( t(J) %*% iS %*% e ) ) - ipC %*% p
      dFdpp = - base::Re( Matrix::as.matrix( t(J) %*% iS %*% J ) ) - ipC
      
      # decrease regularization
      #---------------------------------------------------------------------
      v = min( v + 0.5, 4 )
      cat( "+\n" )
    }
    # if F has decreased, reset changes
    #-----------------------------------------------------------------------
    else
    {
      
      # reset expansion point
      #---------------------------------------------------------------------
      p = C.p
      h = C.h
      Cp = C.Cp
      
      # increase reguarization
      #---------------------------------------------------------------------
      v = min( v - 2, -4 )
      cat( "-\n" )
    }
    
    # E-step: update
    #=======================================================================
    dp = ReDCM_dx( dFdpp, dFdp, v )
    p = p + dp
    Ep = pUnvect( Ep, pEv + V %*% p[ip], ne, nr )
    
    # graphics
    #=======================================================================
    # TODO
    
    # convergence
    #=======================================================================
    dF = t(dFdp) %*% dp
    # TODO print iteration information

    cat(k, '. F: ', C.F, ' dF: ', dF, '\n', sep='')
    
    criterion = c( (dF < 1e-1), criterion[1:length(criterion)-1] )
    if ( all( criterion ) )
    {
      cat(' convergence\n')
      break
    }
  
    #})
    #cat( "running time: ", time1, "\n" )
  }
  
  # figures, stats...
  
  
  # outputs
  #=========================================================================
  Ep = pUnvect( Ep, pEv + V %*% C.p[ip], ne, nr )
  Cp = V %*% C.Cp[ip, ip] %*% t(V)
  Eh = C.h
  Fe = C.F

  output = list( Ep, Cp, Eh, Fe, k )
  
  
  return(output)

}
