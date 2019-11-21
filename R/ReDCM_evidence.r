

ReDCM_evidence = function( DCM )
{

  M = DCM@M[[1]]
  v = M@ns
  n = M@l
  
  # only look at those parameters with non-zero prior covariance
  # -------------------------------------------------------------------------
  wsel = ReDCM_find( diag(M@pC) )[,1]
  
  # look at costs of coding prediction errors by region
  # -------------------------------------------------------------------------
  Ce = array( DCM@Ce )
  if ( is.null(Ce) )
  {
    cat('ReDCM_evidence: hyperparameters not estimated\n')
    return( -1 )
  }
  
  region_cost = rep(0, n)
  if ( length(dim(Ce)) > 1 )
  {
    # Ce is error covariance
    for ( i in 1:n )
    {
      lambda_i = Ce[i*v, i*v]
      region_cost[i] = - 0.5*v*log(lambda_i) - (0.5 * t(DCM@R[,i])) %*% ((1/lambda_i) * diag(v)) %*% DCM@R[,i]
    }
  }
  else if ( dim(array(Ce)) > 1 )
  {
    # Ce is a regional hyperparameter
    for ( i in 1:n )
    {
      lambda_i = Ce[i]
      region_cost[i] = - 0.5*v*log(lambda_i) - (0.5 * t(DCM@R[,i])) %*% ((1/lambda_i) * diag(v)) %*% DCM@R[,i]
    }
  }
  else
  {
    # Ce is the hyperparameter
    for ( i in 1:n )
    {
      lambda_i = Ce
      region_cost[i] = - 0.5*v*log(lambda_i) - (0.5 * t(DCM@R[,i])) %*% ((1/lambda_i) * diag(v)) %*% DCM@R[,i]
    }
  }
  
  # results
  # -------------------------------------------------------------------------
  aic_penalty = length(wsel)
  bic_penalty = 0.5 * length(wsel) * log(v)
  aic_overall = sum(region_cost) - aic_penalty
  bic_overall = sum(region_cost) - bic_penalty
  

  return( list(aic_overall, bic_overall, region_cost, aic_penalty, bic_penalty) )

}
