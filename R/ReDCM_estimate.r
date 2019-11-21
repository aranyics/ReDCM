
#script.dir = dirname(sys.frame(1)$ofile)
#cat('script.dir: ', script.dir, '\n')
#setwd(script.dir)

#suppressPackageStartupMessages( require('R.matlab') )
#suppressPackageStartupMessages( require('Matrix') )
#suppressPackageStartupMessages( require('matrixcalc') )
# suppressPackageStartupMessages( require('expm') )
#suppressPackageStartupMessages( require('numDeriv') )
#suppressPackageStartupMessages( require('plyr') )
#suppressPackageStartupMessages( require('stringr') )
#suppressPackageStartupMessages( require('pracma') )

#source('ReDCM_utils.r')


#DCM.mat = 'data/DCM_natt2mot_d1l.mat'

#################################################################################################
# source:       ReDCM_estimate.r
#------------------------------------------------------------------------------------------------
# type:         Function
#------------------------------------------------------------------------------------------------
# title:        Estimates parameters of a DCM for fMRI data
#------------------------------------------------------------------------------------------------
# description:  ---
#------------------------------------------------------------------------------------------------
# based on:     spm_dcm_estimate.m
#               DCM12 (http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
#------------------------------------------------------------------------------------------------
# reference:    ---
#------------------------------------------------------------------------------------------------
# TODO:         - implement DCM variations (DCM.options)
#################################################################################################


ReDCM_estimate = function(DCM.mat, A=NULL, B=NULL, C=NULL) {
  
  cat(' TEST function \n\n')  
  
  if ( class(DCM.mat) == 'Estimates' )
  {
    M = DCM.mat@M[[1]]
    U = DCM.mat@U
    Y = DCM.mat@Y
    DCM.options = DCM.mat@options
  }
  else
  {
  
    ## load DCM structure
    # --------------------------------------------------------------
    
    DCM = readMat(DCM.mat)
    DCM = DCM$DCM[,,1]
    DCM.options = DCM$options[,,1]
    
    ## swap matrices
    # --------------------------------------------------------------
    
    if ( !is.null(A) )
    {
      # DCM$a = array( A, dim=c(DCM.n,DCM.n) )
      DCM$a = A
    }
    if ( !is.null(B) )
    {
      # DCM$b = array( B, dim=c(DCM.n,DCM.n,DCM.u) )
      DCM$b = B
    }
    if ( !is.null(C) )
    {
      # DCM$c = array( C, dim=c(DCM.n,DCM.u) )
      DCM$c = C
    }
    
    ## initialize estimate struct
    # --------------------------------------------------------------
    
    DCMe = .Estimates()
    DCMe = initEstimates( DCMe, DCM )
    
    M = DCMe@M[[1]]
    U = DCMe@U
    Y = DCMe@Y
  
  }
  
  
  
  ## test functions
  # --------------------------------------------------------------
  #cat('\n________________________________________________\n')
  #cat('Test - return input for ReDCM_kernels (ReDCM_bireduce output)\n')
  #result = ReDCM_bireduce(M, M@pE, TRUE, TRUE)
  #result = append(result, c(M@N, M@dt), after=4)
  #names(result) = c('M0', 'M1', 'L1', 'L2', 'M.N', 'M.dt')
  #cat( '    list of ', length(result), ': M0, M1, L1, L2, M@N, M@dt\n')
  #return( result )
  
  ## estimation
  # --------------------------------------------------------------
  
  #cat(DCM.options$stochastic)
  
  #if ( !DCM.options$stochastic[1] )
  if ( TRUE )
  {
    
    # deterministic DCM
    ####################
    
    cat( "Running deterministic DCM estimation\n------------------------------------------------\n\n" )
    
    # nonlinear system identification (Variational EM) - deterministic DCM
    # ---------------------
    nlsiOut = ReDCM_nlsi_GN(M, U, Y)
    cat( "\nF: ", nlsiOut[[4]], "\n\n")
    
    # predicted responses (y) and residuals (R)
    # ---------------------
    y = c_int_det(nlsiOut[[1]], M, U)
    R = Y@y - y
    R = R - Y@X0 %*% Matrix::solve( t(Y@X0) %*% Y@X0 ) %*% ( t(Y@X0) %*% R )
    Ce = exp( - nlsiOut[[3]] )
    

    # return( list( nlsiOut[[1]], nlsiOut[[2]], nlsiOut[[3]], nlsiOut[[4]], y, M, U, Y ) )
    
  }
  else
  {
    
    # stochastic DCM
    ####################
    # TODO
    
    cat( "Running stochastic DCM estimation\n------------------------------------------------\n\n" )
    cat( "stochastic DCM not implemented yet... exiting\n\n" )
    return( list(M,U,Y) )
    
    # Decimate U.u from micro-time
    # ---------------------
    
    u = U@u
    y = Y@y
    Dy = ReDCM_dctmtx( dim(y)[1], dim(y)[1] )
    Du = ReDCM_dctmtx( dim(u)[1], dim(y)[1] )
    Dy = Dy * sqrt( dim(y)[1] / dim(u)[1] )
    u = matrix( nrow=dim(Dy)[1], ncol=dim(u)[2], Dy %*% ( t(Du) %*% u ) )
    
    # DEM structure
    # ---------------------
    
    M = extendModel(M)
    
    #DEM = .Dem()
    #DEM = setupDEM(M, u, Y)
  
  }
  
  
  # bilinear representation (swapped?)
  # first order hemodynamic (H) and neuronal (K) kernels
  # ---------------------
  if ( ! M@options$nonlinear )
    reduc = c_bireduce( M, nlsiOut[[1]], out=1 )
  else
    reduc = c_soreduce( M, nlsiOut[[1]], out=1 )
  kern = ReDCM_kernels( reduc[[1]], reduc[[2]], reduc[[3]], reduc[[4]], M@N, M@dt )
  H0 = array( kern[[1]] )
  H1 = kern[[2]]
  K0 = kern[[3]]
  K1 = kern[[4]]
  # K2 = kern[[5]]
  
  # bayesian inference and variance
  # ---------------------
  Tp = 1 - ReDCM_Ncdf( pVect( M@pE ), abs(pVect(nlsiOut[[1]])), diag(nlsiOut[[2]]) )
  Pp = .Params()
  Pp = pUnvect( Pp, Tp, M@m, M@l )
  Vp = .Params()
  Vp = pUnvect( Vp, diag(nlsiOut[[2]]), M@m, M@l )
  
  # store parameter estimates
  # ---------------------
  DCMe = .Estimates()
  DCMe = setupEstimates( DCMe, 
                         M=list(M), 
                         Y=Y, 
                         U=U, 
                         Ce=Ce, 
                         Ep=nlsiOut[[1]], 
                         Cp=nlsiOut[[2]], 
                         Pp=Pp, 
                         Vp=Vp, 
                         H1=H1, 
                         K1=K1, 
                         R=R, 
                         y=y, 
                         T0=0,
                         k=nlsiOut[[5]],
                         options=DCM.options )
  
  # compute approximations to model evidence
  # ---------------------
  # evid = ReDCM_evidence( DCMe )
  DCMe = evidence( DCMe )
  DCMe = setupEstimates( DCMe, 
                         Fe=nlsiOut[[4]] )
    
  
  
  #return( c_bireduce( M, nlsiOut[[1]], out=1 ) )
  #return( list( nlsiOut[[1]], nlsiOut[[2]], nlsiOut[[3]], nlsiOut[[4]], y, M, U, Y, H0, H1, K0, K1, Pp ) )
  return(DCMe)
  
  
  cat("\nkesz\n")
  
  #return(M)

}



ReDCM_estimate_loop = function( DCM.mat, output.dir, mListFile=NULL )
{
  
  ## load DCM structure
  # --------------------------------------------------------------
  
  DCM = readMat(DCM.mat)
  DCM = DCM$DCM[,,1]
  DCM.options = DCM$options[,,1]
  DCM.n = length(DCM$Y[,,1]$name)
  DCM.uN = length(DCM$U[,,1]$name)
  
  ## load model list file
  # --------------------------------------------------------------
  
  if ( !is.null(mListFile) )
    modelList = read.csv( mListFile, header=FALSE )
  else
  {
    #modelList = 0 : (2^(DCM.n*DCM.n-DCM.n)-1) # A
    modelList = 0 : (2^(DCM.n*DCM.n*DCM.uN - DCM.uN*(DCM.n + length(which(DCM$a==0)))) - 1) # B
  }
  
  
  for (modelCnt in modelList)
  {
    
    ## set prior parameter matrices
    # --------------------------------------------------------------
    
    B = array( 0, DCM.uN*DCM.n*DCM.n )
    
    B[26] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
    B[25] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
    B[24] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
    B[21] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
    B[17] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
    B[16] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
    B[15] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
    B[12] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
    B[8] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
    B[7] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
    B[6] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
    B[3] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
    
    
#     A = as.array( diag(DCM.n) )
#     
#     A[15] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
#     A[14] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
#     A[13] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
#     A[12] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
#     A[10] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
#     A[9] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
#     A[8] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
#     A[7] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
#     A[5] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
#     A[4] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
#     A[3] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
#     A[2] = modelCnt %% 2; modelCnt = floor( modelCnt / 2 );
    
    #DCM$a = A
    DCM$b = B
    #modelID = str_c(A, sep='',collapse='')
    modelID = str_c(B, sep='',collapse='')
    
    #DCM$b[4,2,3] = 0; DCM$b[2,4,3] = 1; DCM$b[3,4,2] = 1;
    
    dir.create( output.dir, recursive=TRUE )
    write.table(str_c(DCM$a, sep='',collapse=''), str_c(output.dir, '/A3.csv'), row.names=FALSE, col.names=FALSE, append=TRUE)
    write.table(str_c(DCM$b, sep='',collapse=''), str_c(output.dir, '/B3.csv'), row.names=FALSE, col.names=FALSE, append=TRUE)
    write.table(str_c(DCM$c, sep='',collapse=''), str_c(output.dir, '/C3.csv'), row.names=FALSE, col.names=FALSE, append=TRUE)
    write.table(str_c(DCM$d, sep='',collapse=''), str_c(output.dir, '/D3.csv'), row.names=FALSE, col.names=FALSE, append=TRUE)
    next;
    
    
    cat('modelID: ', modelID, ',  density: ', (sum(A)-DCM.n) / (DCM.n*DCM.n-DCM.n), '\n' )
    
    ## unpack DCM structure
    # --------------------------------------------------------------
    
    # inputs
    cat('Reading external inputs - U\n------------------------------------------------\n\n')
    U = .External()
    U = setData(U, DCM)
    # ----------
    
    # data
    cat('Reading data - Y\n------------------------------------------------\n\n')
    Y = .Data()
    Y = setData(Y, DCM)
    # ----------
    
    # complete prior model specification
    cat('\nCreating model specifications - M\n------------------------------------------------\n\n')
    M = .Model()
    M = setupModel(M, DCM, U)
    #return(M)
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
                           options=DCM.options )
    
    ptm = proc.time()
    DCMe = ReDCM_estimate( DCMe )
    ptm = proc.time() - ptm
    cat(ptm, '\n')
    
    ## print results
    # --------------------------------------------------------------
    
    dir.create(output.dir, recursive=TRUE)
    
    # A matrix
    A = format( DCMe@Ep@A, digits=5 )
    colnames(A) = DCMe@Y@name
    write.csv(A, str_c(output.dir, '/', modelID, '_Mr_A.csv'), row.names=FALSE)
    A = format( DCMe@Pp@A, digits=5 )
    colnames(A) = DCMe@Y@name
    write.csv(A, str_c(output.dir, '/', modelID, '_Mr_pA.csv'), row.names=FALSE)
    
    # B matrices
    for ( i in 1:M@m )
    {
      B = format( DCMe@Ep@B[,,i], digits=5 )
      colnames(B) = DCMe@Y@name
      write.csv(B, str_c(output.dir, '/', modelID, '_Mr_B', i, '.csv'), row.names=FALSE)
      B = format( DCMe@Pp@B[,,i], digits=5 )
      colnames(B) = DCMe@Y@name
      write.csv(B, str_c(output.dir, '/', modelID, '_Mr_pB', i, '.csv'), row.names=FALSE)
    }
    
    # C matrix
    C = format( t(DCMe@Ep@C), digits=5 )
    colnames(C) = DCMe@Y@name
    write.csv(C, str_c(output.dir, '/', modelID, '_Mr_C.csv'), row.names=FALSE)
    C = format( t(DCMe@Pp@C), digits=5 )
    colnames(C) = DCMe@Y@name
    write.csv(C, str_c(output.dir, '/', modelID, '_Mr_pC.csv'), row.names=FALSE)
    
    # hemodyn
    Tra = format( t(DCMe@Ep@transit), digits=5 )
    Trap = format( t(DCMe@Pp@transit), digits=5 )
    Dec = format( t(DCMe@Ep@decay), digits=5 )
    Decp = format( t(DCMe@Pp@decay), digits=5 )
    colnames(Tra) = DCMe@Y@name
    colnames(Trap) = DCMe@Y@name
    colnames(Dec) = DCMe@Y@name
    colnames(Decp) = DCMe@Y@name
    
    write.csv(Tra, str_c(output.dir, '/', modelID, '_Mr_transit.csv'), row.names=FALSE)
    write.csv(Trap, str_c(output.dir, '/', modelID, '_Mr_ptransit.csv'), row.names=FALSE)
    write.csv(Dec, str_c(output.dir, '/', modelID, '_Mr_decay.csv'), row.names=FALSE)
    write.csv(Decp, str_c(output.dir, '/', modelID, '_Mr_pdecay.csv'), row.names=FALSE)
    write.table(format( DCMe@Ep@epsilon, digits=5 ), str_c(output.dir, '/', modelID, '_Mr_epsilon.csv'), row.names=FALSE, col.names=FALSE, sep=',')
    write.table(format( DCMe@Pp@epsilon, digits=5 ), str_c(output.dir, '/', modelID, '_Mr_pepsilon.csv'), row.names=FALSE, col.names=FALSE, sep=',')
    
    # evidence and running time
    write.table(format( DCMe@Fe, digits=5 ), str_c(output.dir, '/', modelID, '_Mr_F.csv'), row.names=FALSE, col.names=FALSE, sep=',')
    write.table(format( ptm[3], digits=5 ), str_c(output.dir, '/', modelID, '_Mr_T.csv'), row.names=FALSE, col.names=FALSE, sep=',')
    write.table(DCMe@k, str_c(output.dir, '/', modelID, '_Mr_k.csv'), row.names=FALSE, col.names=FALSE, sep=',')
    
    #return(DCMe)
    
  }
  
  return( DCMe )
  
}


ReDCM_prepare_dcm = function( DCM.mat )
{
  
  ## load DCM structure from file if necessary
  # --------------------------------------------------------------
  if ( is.character(DCM.mat) )
  {
    DCM = readMat(DCM.mat)
    DCM = DCM$DCM[,,1]
  }
  else
  {
    DCM = DCM.mat
  }
  DCM.options = DCM$options[,,1]
  DCM.n = length(DCM$Y[,,1]$name)
  
  
  ## unpack DCM structure
  # --------------------------------------------------------------
  # inputs
  cat('Reading external inputs - U\n------------------------------------------------\n\n')
  U = .External()
  U = setData(U, DCM)
  # ----------
  
  # data
  cat('Reading data - Y\n------------------------------------------------\n\n')
  Y = .Data()
  Y = setData(Y, DCM)
  # ----------
  
  # complete prior model specification
  cat('\nCreating model specifications - M\n------------------------------------------------\n\n')
  M = .Model()
  M = setupModel(M, DCM, U)
  # ----------
  
  # correct external input length
  if ( (M@Ydt / U@dt) * M@ns != dim(U@u)[1] )
  {
    U = extCorr(U, M)
  }
  
  
  ## setup DCM
  # --------------------------------------------------------------
  DCMe = .Estimates()
  DCMe = setupEstimates( DCMe, 
                         M=list(M), 
                         Y=Y, 
                         U=U, 
                         T0=0,
                         options=DCM.options )
  
  return( DCMe )
}
