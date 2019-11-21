

#################################################################################################
# source:       ReDCM_model.R
#------------------------------------------------------------------------------------------------
# type:         Class definition
#------------------------------------------------------------------------------------------------
# title:        Model specifications
#------------------------------------------------------------------------------------------------
# description:  Contains model specifications for DCM.
#
#               Slots:  delays      - TR (repetition time)
#                       TE          - TE (echo time)
#                       x           - initial condition (states)
#                       pE          - prior expectation (parameters)
#                       pC          - prior covariance  (parameters)
#                       hE          - prior expectation (precisions)
#                       hC          - prior covariance  (precisions)
#                       m           - number of inputs
#                       n           - number of state variables
#                       l           - number of nodes
#                       N           - number of maximum EM iterations
#                       dt          - 
#                       Ydt         - TR
#                       ns          - number of scans
#                       IS          - fn handler: integration of MIMO bilinear system dx/dt
#                       f           - fn handler: equations of motion
#                       g           - fn handler: observation equation
#                       H           - hemodynamic constants for Balloon-model
#                       -----------------------------------------------------------
#                       E           - inversion parameters
#                         - form    - form of random fluctuations
#                         - s       - smoothness of fluctuations
#                         - d       - embedding dimension
#                         - n       - embedding dimension
#                         - nN      - maximum number of iterations
#                       Q           - precision components
#                       xP          - precision (hidden-sates)
#                       W           - precision (hidden-motion)
#                       V           - precision (hidden-cause)
#------------------------------------------------------------------------------------------------
# based on:     M structure
#               DCM12 (http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
#------------------------------------------------------------------------------------------------
# TODO:         - use appropriate R matrix representation
#################################################################################################



#source('ReDCM_parameters.R')

.Model <- setClass("Model",
                   
                   slots = list(delays  = "numeric",
                                TE      = "numeric",
                                x       = "matrix",
                                pE      = "Params",
                                pC      = "matrix",
                                hE      = "matrix",
                                hC      = "matrix",
                                m       = "numeric",
                                n       = "numeric",
                                l       = "numeric",
                                N       = "numeric",
                                dt      = "numeric",
                                Ydt     = "numeric",
                                ns      = "numeric",
                                IS      = "character",
                                f       = "character",
                                g       = "character",
                                H       = "numeric",
                                options = "list",
                                E       = "list",
                                Q       = "list",
                                xP      = "numeric",
                                W       = "matrix",
                                V       = "numeric"),
                   
                   prototype = list(delays  = double(),
                                    TE      = double(),
                                    x       = matrix(),
                                    pE      = new("Params"),
                                    pC      = matrix(),
                                    hE      = matrix(),
                                    hC      = matrix(),
                                    m       = integer(),
                                    n       = integer(),
                                    l       = integer(),
                                    N       = integer(),
                                    Ydt     = double(),
                                    dt      = double(),
                                    ns      = integer(),
                                    IS      = 'ReDCM_int',
                                    f       = 'ReDCM_fx_fmri',
                                    g       = 'ReDCM_gx_fmri',
                                    H       = c(0.64, 0.32, 2.0, 0.32, 0.4),
                                    options = list(),
                                    E       = list(),
                                    Q       = list(),
                                    xP      = double(),
                                    W       = matrix(),
                                    V       = double())
)


setMethod(f           = "setupModel",
          signature   = "Model",
          definition  = function(.Object, DCM, Ext, H=NULL){
            
            DCM.n = length(DCM$Y[,,1]$name)
            DCM.v = dim(DCM$Y[,,1]$y)[1]
            DCM.options = DCM$options[,,1]
            
            pA = 64
            dA = 1
            pE = .Params()
            pE = setPriors(pE, DCM)
            pC = setPriorCov(pE, DCM, pA)
            
            .Object@delays  = as.vector(DCM$delays)                               # TR (repetition time)
            .Object@TE      = DCM$TE[1]                                           # TE (echo time)
            if ( DCM.options$two.state )
              .Object@x       = matrix(nrow=DCM.n, ncol=6, data=0)                # initial condition (states)
            else
              .Object@x       = matrix(nrow=DCM.n, ncol=5, data=0)                # initial condition (states)
            .Object@pE      = pE                                                  # prior expectation (parameters)
            .Object@pC      = pC                                                  # prior covariance  (parameters)
            #.Object@hE      = sparseMatrix(i=c(1:DCM.n), j=rep(1,DCM.n), x=6)     # prior expectation (precisions)
            #.Object@hC      = sparseMatrix(i=c(1:DCM.n), j=c(1:DCM.n), x=1/128)   # prior covariance  (precisions)
            .Object@hE      = matrix(data=6, nrow=DCM.n, ncol=1)                  # prior expectation (precisions)
            .Object@hC      = diag( rep(x=1/128, times=DCM.n) )                   # prior covariance  (precisions)
            .Object@m       = dim(Ext@u)[2]                                       # number of inputs
            .Object@n       = length(.Object@x)                                   # number of state variables
            .Object@l       = length(.Object@x[,1])                               # number of nodes
            .Object@N       = maxIter                                             # number of maximum EM iterations
            .Object@dt      = 16/.Object@N
            .Object@Ydt     = DCM$Y[,,1]$dt[1]                                    # TR
            .Object@ns      = DCM.v                                               # number of scans
            #if ( DCM.options$nonlinear )
            #  .Object@IS      = '_r_int_det_D'                                    # integration of MIMO bilinear system dx/dt
            #else
              .Object@IS      = '_r_int_det'                                      # integration of MIMO bilinear system dx/dt
            .Object@f       = '_r_fx_fmri'                                        # equations of motion
            .Object@g       = '_r_gx_fmri'                                        # observation equation
            .Object@options = DCM.options
            
            if ( !is.null(H) )
            {
              if ( length(H) == 5 )
              {
                .Object@H = H
              }
              else
              {
                cat("Warning: Build model: wrong size of H (remains default)")
              }
            }
            
            remove(pE)
            
            return(.Object)
          }
)


setMethod(f = "extendModel",
          signature = "Model",
          #definition = function(.Object, Y=NULL, V=exp(16)){
          definition = function(.Object, ...){
            
            Y = NULL
            V = exp(16)
            
            if ( class(.Object)[1] != "Model" )
            {
              cat('Error: Wrong type of object.')
              return(0)
            }
            
            argin = list(...)
            if ( any( c(1,2) == length(argin) ) )
            {
              for ( i in 1:length(argin) )
              switch( class(argin[[i]])[1],
                      "Data"      = { Y = argin[[i]] },
                      "numeric"   = { V = argin[[i]] },
                                    { cat('Error: Unrecognised argument!')
                                      return(.Object) }
              )
            }
            else if( length(argin) != 0 )
            {
              cat('Error: Wrong number of arguments!')
              return(.Object)
            }            
            
            
            if ( !is.null(Y) )
            {
              n = dim(Y@y)[2]
              
              # set inversion parameters
              # --------------------------------------------------------
              E = list(form     = 'Gaussian',
                       s        = 1/2,
                       d        = 2,
                       n        = 4,
                       nN       = maxIter)
              
              .Object@E         = E
              
              # DEM works in time bins, not seconds
              # --------------------------------------------------------
              .Object@delays = .Object@delays / Y@dt
            
              # Specify hyper-priors on (log-precision of) observation noise
              # --------------------------------------------------------
              .Object@Q         = ReDCM_Ce( rep(1, n) )
              .Object@hE        = .Object@hE
              .Object@hC        = .Object@hC
              
              # allow (only) neuronal [x, s, f, q, v] hidden states to fluctuate
              # --------------------------------------------------------
              W = rep(1, n) %*% t(matrix( c(12, 16, 16, 16, 16) ))
              .Object@xP        = exp(6)
              .Object@W         = matrix( diag(exp( W )) )
            }
            else
            {
              .Object@V         = V
            }
            
            return(.Object)
            
          }
)


setMethod(f = "setHemodyn",
          signature = "Model",
          #definition = function(.Object, Y=NULL, V=exp(16)){
          definition = function(.Object, H){
            
            if ( length(H == 5) )
            {
              .Object@H = H
            }
            else
            {
              cat("Warning: Model set hemodynamic constants: wrong size of H")
            }
            

            return(.Object)

          }
)
