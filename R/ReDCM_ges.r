
ReDCM_GES = function( DCM.mat, search, out.dir=NULL, A=NULL, B=NULL, C=NULL, missB=NULL, iter=0, backward=TRUE )
{

  DCM = readMat(DCM.mat)
  DCM = DCM$DCM[,,1]
  DCM.n = length(DCM$Y[,,1]$name)
  DCM.u = length(DCM$U[,,1]$name)
  
  #iter = 0
  DCMe.traj = NULL
  ges.dir=NULL
  
  if ( backward )
  {
    ges.dir = str_c(out.dir, '/ges/', search, '_bw', sep='')
  }
  else
  {
    ges.dir = str_c(out.dir, '/ges/', search, sep='')
  }
  
  #A.len = DCM.n * DCM.n - DCM.n
  A.en = c(0, which( (matrix(1,DCM.n,DCM.n) - diag(DCM.n)) == 1 ))
  #B.len = (length( which( DCM$a != 0 ) ) - DCM.n) * DCM.u
  B.en = c(0, maskeach.mx( which((DCM$a-diag(DCM.n)) != 0), DCM.u, DCM.n*DCM.n, missB ))
  #C.len = DCM.n * DCM.u
  C.en = seq(0, DCM.n*DCM.u)
  
  if ( !is.null(A) )
  {
    DCM$a = array( A, dim=c(DCM.n,DCM.n) )
    B.en = c(0, maskeach.mx( which((DCM$a-diag(DCM.n)) != 0), DCM.u, DCM.n*DCM.n, missB ))
  }
  if ( !is.null(B) )
  {
    DCM$b = array( B, dim=c(DCM.n,DCM.n,DCM.u) )
  }
  if ( !is.null(C) )
  {
    DCM$c = array( C, dim=c(DCM.n,DCM.u) )
  }
  
  
  # search methods by matrix type
  # -------------------------------------------------------------
  
  if ( search == 'A' )
  {
    
    DCM$b = DCM$b * 0
    if ( backward )
    {
      DCM$a = matrix(1,DCM.n,DCM.n)
    }
    else
    {
      DCM$a = diag( DCM.n )
    }
    
    peak.F = -Inf
    best.F = -Inf
    best.idx = NULL
    repeat
    {
      
      A.len = length( A.en )
      if ( A.len == 0 )
      {
        return( DCM$a )
      }
      
      for ( i in A.en )
      {
        if ( i != 0 )
        {
          if ( backward )
          {
            DCM$a[i] = 0
          }
          else
          {
            DCM$a[i] = 1
          }
        }
        
        cat( str_c(DCM$a, sep='', collapse=''), '\n' )
        
        DCMe = ReDCM_GES_setmodel( DCM )
        
        ptm = proc.time()
        DCMe = ReDCM_estimate( DCMe )
        ptm = proc.time() - ptm
        DCMe = setupEstimates( DCMe, tt=ptm[3] )
        
        if ( DCMe@Fe > best.F )
        {
          best.F = DCMe@Fe
          best.idx = i
          DCMe.traj = DCMe
        }
        
        if ( backward )
        {
          DCM$a[i] = 1
        }
        else
        {
          DCM$a[i] = 0
        }
      }
      
      if ( best.F > peak.F )
      {
        peak.F = best.F
        iter = iter + 1
        dir.create( ges.dir, showWarnings=FALSE, recursive=TRUE )
        save( DCMe.traj, file=str_c(ges.dir, '/DCM_', iter, sep='') )
      }
      else
      {
        return( DCM$a )
      }
      
      if ( best.idx == 0 )
      {
        return( DCM$a )
      }
      
      if ( backward )
      {
        DCM$a[best.idx] = 0
      }
      else
      {
        DCM$a[best.idx] = 1
      }
      
      A.en = A.en[-which(A.en == best.idx)]
      if ( A.en[1] == 0 )
      {
        A.en = A.en[-1]
      }
      
    }
    
    return( DCM$a )
  
  }
  else if ( search == 'B' )
  {
  
    if ( backward )
    {
      DCM$b = DCM$b * 0
      DCM$b[B.en] = DCM$b[B.en] + 1
    }
    else
    {
      DCM$b = DCM$b * 0
    }
    
    peak.F = -Inf
    best.F = -Inf
    best.idx = NULL
    repeat
    {
      
      best.F = -Inf
      
      B.len = length( B.en )
      if ( B.len == 0 )
      {
        return( DCM$b )
      }
      
      for ( i in B.en )
      {
        if ( i != 0 )
        {
          if ( backward )
          {
            DCM$b[i] = 0
          }
          else
          {
            DCM$b[i] = 1
          }
        }
        
        cat( str_c(DCM$b, sep='', collapse=''), '\n' )
        
        DCMe = ReDCM_GES_setmodel( DCM )
        
        ptm = proc.time()
        DCMe = ReDCM_estimate( DCMe )
        ptm = proc.time() - ptm
        DCMe = setupEstimates( DCMe, tt=ptm[3] )
        
        if ( DCMe@Fe > best.F )
        {
          best.F = DCMe@Fe
          best.idx = i
          DCMe.traj = DCMe
        }
        
        if ( backward )
        {
          DCM$b[i] = 1
        }
        else
        {
          DCM$b[i] = 0
        }
      }
      
      if ( best.F >= peak.F )
      {
        peak.F = best.F
        iter = iter + 1
        dir.create( ges.dir, showWarnings=FALSE, recursive=TRUE )
        save( DCMe.traj, file=str_c(ges.dir, '/DCM_', iter, sep='') )
      }
      else
      {
        return( DCM$b )
      }
      
      if ( best.idx == 0 )
      {
        return( DCM$b )
      }
      
      if ( backward )
      {
        DCM$b[best.idx] = 0
      }
      else
      {
        DCM$b[best.idx] = 1
      }
      
      B.en = B.en[-which(B.en == best.idx)]
      if ( length(B.en) > 0 )
      {
        if ( B.en[1] == 0 )
        {
          B.en = B.en[-1]
        }
      }
      
    }
    
    return( DCM$b )
    
  }
  else if ( search == 'C' )
  {
  
  }
  else if ( search == 'AB' )
  {

    if ( backward )
    {
      if ( is.null(A) )
      {
        DCM$a = matrix(1,DCM.n,DCM.n)
      }
      if ( is.null(B) )
      {
        B.en = c(0, maskeach.mx( which((DCM$a-diag(DCM.n)) != 0), DCM.u, DCM.n*DCM.n, missB ))
        DCM$b = DCM$b * 0
        DCM$b[B.en] = DCM$b[B.en] + 1
      }
    }
    else
    {
      DCM$b = DCM$b * 0
      DCM$a = diag( DCM.n )
    }
    #DCM$a = array(c(1,1,1,1,1,1,0,1,1), dim=c(3,3))
    #DCM$b = array(c(0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,0,0,0,1,1,1,0,1,0,1,0), dim=c(3,3,3))
    AB = c(DCM$a, DCM$b)
    AB.bckp=NULL
    
    peak.F = -Inf
    best.F = -Inf
    best.idx = NULL
    repeat
    {

      if ( backward )
      {
        A.en = c(which( AB[1:(DCM.n*DCM.n)] - diag(DCM.n) == 1 ))
        B.en = c(maskeach.mx( which( (AB[1:(DCM.n*DCM.n)]-c(diag(DCM.n)) ) != 0), DCM.u, DCM.n*DCM.n, missB ))
        AB.en = c(A.en, B.en+DCM.n*DCM.n) # null modell kihagyva itt
        AB.len = length(AB.en)
      }
      else
      {
        A.en = c(which( AB[1:(DCM.n*DCM.n)] == 0 ))
        B.en = c(maskeach.mx( which( (AB[1:(DCM.n*DCM.n)]-c(diag(DCM.n)) ) != 0), DCM.u, DCM.n*DCM.n, missB ))
        AB.en = c(A.en, B.en+DCM.n*DCM.n) # null modell kihagyva itt
        AB.len = length(AB.en)
      }

      
      if ( AB.len != 0 )
      {
        occupied = NULL
        for( i in 1:AB.len )
        {
          if ( backward )
          {
            if ( AB[AB.en[i]] == 0 && AB.en[i] > DCM.n*DCM.n )
            {
              occupied = c(occupied, i)
            }
          }
          else
          {
            if ( AB[AB.en[i]] != 0 && AB.en[i] > DCM.n*DCM.n )
            {
              occupied = c(occupied, i)
            }
          }
        }
        if ( !is.null(occupied) )
        {
          AB.en = AB.en[-occupied]
        }
      }        
      
      AB.len = length( AB.en )
      if ( AB.len == 0 )
      {
        return( AB )
      }
      
      for ( i in AB.en )
      {
        if ( i != 0 )
        {
          if ( backward )
          {
            AB.bckp = AB
            AB[i] = 0
            if ( i <= DCM.n*DCM.n )
            {
              B.en2 = c(maskeach.mx( which( (AB[1:(DCM.n*DCM.n)]-c(diag(DCM.n)) ) != 0), DCM.u, DCM.n*DCM.n, missB ))
              AB.tmp = AB[-(1:(DCM.n*DCM.n))]*0
              AB.tmp[B.en2] = AB.tmp[B.en2] + 1
              AB = c(AB[1:(DCM.n*DCM.n)], AB[-(1:(DCM.n*DCM.n))] * AB.tmp )
            }
          }
          else
          {
            AB[i] = 1
          }
          DCM$a = array( AB[1:(DCM.n*DCM.n)], dim=c(DCM.n,DCM.n))
          DCM$b = array( AB[-(1:(DCM.n*DCM.n))], dim=c(DCM.n,DCM.n, DCM.u))
        }
        
        cat( str_c(AB, sep='', collapse=''), '\n' )
        
        DCMe = ReDCM_GES_setmodel( DCM )
        
        ptm = proc.time()
        DCMe = ReDCM_estimate( DCMe )
        ptm = proc.time() - ptm
        DCMe = setupEstimates( DCMe, tt=ptm[3] )
        
        if ( DCMe@Fe > best.F )
        {
          best.F = DCMe@Fe
          best.idx = i
          DCMe.traj = DCMe
        }
        
        if ( backward )
        {
          AB = AB.bckp
        }
        else
        {
          AB[i] = 0
        }

      }
      
      if ( best.F > peak.F )
      {
        peak.F = best.F
        iter = iter + 1
        dir.create( ges.dir, showWarnings=FALSE, recursive=TRUE )
        save( DCMe.traj, file=str_c(ges.dir, '/DCM_', iter, sep='') )
      }
      else
      {
        return( AB )
      }
      
      if ( best.idx == 0 )
      {
        return( AB )
      }
      
      if ( backward )
      {
        AB[best.idx] = 0
      }
      else
      {
        AB[best.idx] = 1
      }

      AB.en = AB.en[-which(AB.en == best.idx)]
      if ( AB.en[1] == 0 )
      {
        AB.en = AB.en[-1]
      }
      
    }
    
    return( AB )
  }

}



ReDCM_GES_setmodel = function( DCM )
{
  ## unpack DCM structure
  # --------------------------------------------------------------
  
  # inputs
  # cat('Reading external inputs - U\n------------------------------------------------\n\n')
  U = .External()
  U = setData(U, DCM)
  # ----------
  
  # data
  # cat('Reading data - Y\n------------------------------------------------\n\n')
  Y = .Data()
  Y = setData(Y, DCM)
  # ----------
  
  # complete prior model specification
  # cat('\nCreating model specifications - M\n------------------------------------------------\n\n')
  M = .Model()
  M = setupModel(M, DCM, U)
  # ----------
  
  # correct external input length
  if ( (M@Ydt / U@dt) * M@ns != dim(U@u)[1] )
  {
    U = extCorr(U, M)
  }
  
  ## estimate DCM
  # --------------------------------------------------------------
  
  DCMe = .Estimates()
  DCMe = setupEstimates( DCMe, 
                         M=list(M), 
                         Y=Y, 
                         U=U, 
                         T0=0,
                         options=DCM$options[,,1] )
  
  return( DCMe )
}


maskeach.mx = function( v, u, l, miss=NULL )
{
  
  k = NULL
  j = 1:u
  if ( !is.null(miss) )
  {
    for( i in miss )
    {
      j = j[-which(j==i)]
    }
  }
  
  for ( i in j-1 )
  {
    k = c( k, v + i*l )
  }
  
  return( k )
  
}
