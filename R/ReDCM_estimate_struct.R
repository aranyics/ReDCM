

#################################################################################################
# source:       ReDCM_estimate_struct.R
#------------------------------------------------------------------------------------------------
# type:         Class definition
#------------------------------------------------------------------------------------------------
# title:        Model estimates for DCM
#------------------------------------------------------------------------------------------------
# description:  Contains priors and results for model estimation.
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
# based on:     DCM structure
#               DCM12 (http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
#------------------------------------------------------------------------------------------------
# TODO:         - 
#################################################################################################




.Estimates <- setClass("Estimates",
                   
                   slots = list(M       = "list",
                                Y       = "Data",
                                U       = "External",
                                Ce      = "matrix",
                                Ep      = "Params",
                                Cp      = "matrix",
                                Pp      = "Params",
                                Vp      = "Params",
                                H1      = "array",
                                K1      = "array",
                                R       = "matrix",
                                y       = "matrix",
                                T0      = "numeric",
                                qU      = "array",
                                qP      = "array",
                                qH      = "array",
                                Fe      = "numeric",
                                AIC     = "numeric",
                                BIC     = "numeric",
                                ID      = "numeric",
                                ver     = "numeric",
                                options = "list",
                                xY      = "list",
                                k       = "numeric",
                                tt       = "numeric"),
                   
                   prototype = list(M       = list(),
                                    Y       = new("Data"),
                                    U       = new("External"),
                                    Ce      = matrix(),
                                    Ep      = new("Params"),
                                    Cp      = matrix(),
                                    Pp      = new("Params"),
                                    Vp      = new("Params"),
                                    H1      = array(),
                                    K1      = array(),
                                    R       = matrix(),
                                    y       = matrix(),
                                    T0      = numeric(),
                                    qU      = array(),
                                    qP      = array(),
                                    qH      = array(),
                                    Fe      = numeric(),
                                    AIC     = numeric(),
                                    BIC     = numeric(),
                                    ID      = numeric(),
                                    ver     = numeric(),
                                    options = list(),
                                    xY      = list(),
                                    k       = numeric(),
                                    tt      = numeric())
)



setMethod(f           = "setupEstimates",
          signature   = "Estimates",
          definition  = function(.Object,
                                M=NULL,
                                Y=NULL,
                                U=NULL,
                                Ce=NULL,
                                Ep=NULL,
                                Cp=NULL,
                                Pp=NULL,
                                Vp=NULL,
                                H1=NULL,
                                K1=NULL,
                                R=NULL,
                                y=NULL,
                                T0=NULL,
                                qU=NULL,
                                qP=NULL,
                                qH=NULL,
                                Fe=NULL,
                                AIC=NULL,
                                BIC=NULL,
                                ID=NULL,
                                ver=NULL,
                                options=NULL,
                                xY=NULL,
                                k=NULL,
                                tt=NULL){
            
            if ( !is.null(M) & class(M) == "list" & class(M[[1]]) == "Model" )
              .Object@M = M
            if ( !is.null(Y) & class(Y) == "Data" )
              .Object@Y = Y
            if ( !is.null(U) & class(U) == "External" )
              .Object@U = U
            if ( !is.null(Ce) )
              .Object@Ce = Ce
            if ( !is.null(Ep) & class(Ep) == "Params" )
              .Object@Ep = Ep
            if ( !is.null(Cp) )
              .Object@Cp = Cp
            if ( !is.null(Pp) & class(Pp) == "Params" )
              .Object@Pp = Pp
            if ( !is.null(Vp) & class(Vp) == "Params" )
              .Object@Vp = Vp
            if ( !is.null(H1) )
              .Object@H1 = H1
            if ( !is.null(K1) )
              .Object@K1 = K1
            if ( !is.null(R) )
              .Object@R = R
            if ( !is.null(y) )
              .Object@y = y
            if ( !is.null(T0) )
              .Object@T0 = T0
            if ( !is.null(qU) )
              .Object@qU = qU
            if ( !is.null(qP) )
              .Object@qP = qP
            if ( !is.null(qH) )
              .Object@qH = qH
            if ( !is.null(Fe) )
              .Object@Fe = as.numeric(Fe)
            if ( !is.null(AIC) )
              .Object@AIC = AIC
            if ( !is.null(BIC) )
              .Object@BIC = BIC
            if ( !is.null(ID) )
              .Object@ID = ID
            if ( !is.null(ver) )
              .Object@ver = ver
            if ( !is.null(options) & class(options) == "list" )
              .Object@options = options
            if ( !is.null(xY) & class(xY) == "list" )
              .Object@xY = xY
            if ( !is.null(k) )
              .Object@k = k
            if ( !is.null(tt) )
              .Object@tt = tt
            
            return ( .Object )
          }
)



setMethod(f           = "initEstimates",
          signature   = "Estimates",
          definition  = function(.Object, DCM){
            
            
            # input
            # cat('Reading external input - U\n------------------------------------------------\n\n')
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
            
            ## init DCM
            # --------------------------------------------------------------
            
            .Object = setupEstimates( .Object, 
                                      M=list(M), 
                                      Y=Y, 
                                      U=U, 
                                      T0=0,
                                      options=DCM$options[,,1] )
            
            return ( .Object )
          }
)


setMethod(f           = "evidence",
          signature   = "Estimates",
          definition  = function(.Object){
            
            evid = ReDCM_evidence( .Object )
            
            .Object@AIC = evid[[1]]
            .Object@BIC = evid[[2]]
            
            return ( .Object )
          }
)

