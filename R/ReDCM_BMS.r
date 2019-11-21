

ReDCM_run_BMS = function( Fe, method, model.names=NULL )
{
  
  N = NULL                                # number of subjects
  M = NULL                                # number of models
  sumF = 0
  if ( length( dim(Fe) ) < 2 )
  {
    # one subject -> fixed effect
    # method = 'ffx'
    N = 1
    M = length(Fe)
    sumF = Fe
  }
  else
  {
    N = dim(Fe)[1]
    M = dim(Fe)[2]
    sumF = colMeans(Fe)
  }
  
  if ( M < 2 )
  {
    cat('ReDCM_BMS: select more than one model\n')
    return ( 0 )
  }
  
  if ( method != 'ffx' && method != 'rfx' )
  {
    cat('ReDCM_BMS: unknown method (use \'ffx\' or \'rfx\')\n')
    return ( 0 )
  }
  
  
  # Bayesian Model Selection
  # ====================================================================================
  if ( method == 'ffx' )                # single subject, or fixed effects group BMS
  {
    names( sumF ) = model.names
    sumF = sumF - min( sumF )
    i = sumF < ( max( sumF ) - 32 )
    P = sumF
    P[which(i)] = max( sumF ) - 32
    P = P - min( P )
    P = exp( P )
    P = P / sum( P )
    
    barplot( sumF, ylab='Log-evidences', main=str_c('Bayesian Model Selection: ', toupper(method)) )
    barplot( P, ylim=c(0,1), ylab='Model Posterior Probability', main=str_c('Bayesian Model Selection: ', toupper(method)) )
  }
  else                                  # random effects BMS
  {
    bms     = ReDCM_BMS_gibbs( Fe )
    exp_r   = bms[[1]]
    xp      = bms[[2]]
    r_samp  = bms[[3]]
    g_post  = bms[[4]]
    
    # compute protected xp's
    # bms     = spm_BMS( Fe )
    # pxp     = bms[[1]]
    # bor     = bms[[2]]
    
    names( exp_r ) = model.names
    names( xp ) = model.names
    
    barplot( exp_r, ylab='Model Expected Probability', main=str_c('Bayesian Model Selection: ', toupper(method)) )
    barplot( xp, ylim=c(0,1), ylab='Model Exceedance Probability', main=str_c('Bayesian Model Selection: ', toupper(method)) )
  }
  
  return( 0 )
  
}


ReDCM_BMS_gibbs = function( lme, Nsamp=1e4, alpha0=NULL )
{
  
  if ( is.null(dim(lme)) && length(lme) > 1 )
  {
    lme = t( matrix(lme) )
  }
  N = dim(lme)[1]                                # number of subjects
  M = dim(lme)[2]                                # number of models
  if ( is.null(alpha0) )
  {
    # prior observations
    alpha0 = rep( 1, M )
  }
  
  # initialize: sample r from prior
  r = t( matrix( rep( 0, M ) ) )
  for ( k in 1:M )
  {
    r[,k] = rgamma( 1, shape=alpha0[k], scale=1 )
  }
  sr = rowSums( r )
  r = t( matrix( apply( r, 2, function(x){ x = x / sr } ) ) )
  
  # subtract subject means
  lme = lme - rowMeans(lme) %*% t(rep(1, M))
  
  # ensure all log model evidence differences are now within machine range
  max_val = log( .Machine$double.xmax )
  max_val = max_val / M
  for ( i in 1:N )
  {
    for ( k in 1:M )
    {
      lme[i,k] = sign( lme[i,k] ) * min( max_val, abs( lme[i,k] ) )
    }
  }
  
  # Gibbs sampling
  r_samp = array( 0, dim=c(Nsamp, M) )
  g_post = array( 0, dim=c(N, M) )
  gmat = array( 0, dim=c(N, M) )
  
  for ( samp in 1:(2*Nsamp) )
  {
  
    mod_vec = array( 0, dim=c(N, M) )
    
    # sample m's given y, r
    for ( i in 1:N )
    {
      # pick a model for this subject
      u = exp( lme[i,] + log(r) ) + .Machine$double.eps
      g = u / sum(u)
      gmat[i,] = g
      modnum = ReDCM_find3( rmultinom( 1, size=1, prob=g ) )[,1]
      mod_vec[i, modnum] = 1
    }
    
    # sample r's given y, m
    beta = colSums( mod_vec )
    alpha = alpha0 + beta
    for ( k in 1:M )
    {
      r[,k] = rgamma( 1, shape=alpha[k], scale=1 )
    }
    sr = rowSums( r )
    r = t( matrix( apply( r, 2, function(x){ x = x / sr } ) ) )
    
    # only keep last Nsamp samples
    if ( samp > Nsamp )
    {
      r_samp[samp-Nsamp,] = c(r)
      g_post = g_post + gmat
    }
    
    if ( samp %% (Nsamp/2) == 0 )
    {
      cat(samp, ' samples out of ', 2*Nsamp, '\n')
    }
  
  }
  
  g_post = g_post / Nsamp
  
  # Posterior Mean
  exp_r = colMeans( r_samp )
  
  # Exceedence Probabilities
  xp = rep( 0, M )
  j = max.col( r_samp )
  tmp = hist( j, breaks=0:M, include.lowest=TRUE, plot=FALSE )$counts
  xp = tmp / Nsamp
  
  
  return ( list(exp_r, xp, r_samp, g_post) )
  
}

