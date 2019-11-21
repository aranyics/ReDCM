

#################################################################################################
# source:       ReDCM_parameters.R
#------------------------------------------------------------------------------------------------
# type:         Class definition
#------------------------------------------------------------------------------------------------
# title:        Connection and hemodynamic parameter expectations
#------------------------------------------------------------------------------------------------
# description:  Contains expectations of connection and hemodynamic parameter distribution.
#
#               Slots:  A       - intrinsic (endogenous) connections
#                       B       - (bilinear) modulation on connections
#                       C       - exogenous connections
#                       D       - nonlinear modulations
#                       transit
#                       decay
#                       epsilon
#------------------------------------------------------------------------------------------------
# based on:     DCM.Ep structure
#               DCM12 (http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
#------------------------------------------------------------------------------------------------
# TODO:         - use appropriate R matrix representation
#################################################################################################



.Params <- setClass("Params",
                    
                    slots = list(A            = "array",
                                 B            = "array",
                                 C            = "array",
                                 D            = "array",
                                 transit      = "array",
                                 decay        = "array",
                                 epsilon      = "array"),
                    
                    prototype = list(A        = array(),
                                     B        = array(),
                                     C        = array(),
                                     D        = array(),
                                     transit  = array(),
                                     decay    = array(),
                                     epsilon  = array() )
)



setMethod( f          = "setPriors",
           signature  = "Params",
           definition = function(.Object, DCM){
             
             DCM.options = DCM$options[,,1]
             DCM.n = length(DCM$a[,1])
             
             if ( !DCM.options$two.state )
             {
               DCM.a = (DCM$a != 0) * 1                         # A = logical(A)
               DCM.a = DCM.a - diag(diag( DCM.a ))
               .Object@A          = (DCM.a / 128)
               .Object@B          = (DCM$b*0)
               .Object@C          = (DCM$c*0)
               .Object@D          = (DCM$d*0)
               .Object@transit    = array(dim=DCM.n, data=0)
               .Object@decay      = array(dim=DCM.n, data=0)
               .Object@epsilon    = array(dim=1, data=0)
             }
             else  # two_state
             {
               DCM.a = (DCM$a != 0) * 1                         # A = logical(A)
               DCM.a = DCM.a - diag(diag( DCM.a ))
               .Object@A          = (DCM.a + diag(DCM.n)) * 32 - 32
               .Object@B          = (DCM$b*0)
               .Object@C          = (DCM$c*0)
               .Object@D          = (DCM$d*0)
               .Object@transit    = array(dim=DCM.n, data=0)
               .Object@decay      = array(dim=DCM.n, data=0)
               .Object@epsilon    = array(dim=1, data=0)
             }
             
             return(.Object)
             
           }
)


setMethod( f          = "setPriorCov",
           signature  = "Params",
           definition = function(.Object, DCM, pA){
             
             DCM.options = DCM$options[,,1]
             DCM.n = length(DCM$a[,1])
             
             if ( !DCM.options$two.state )
             {
               if ( length(dim(.Object@A)) == 2 )
               {
                 pC.A = ( ((.Object@A != 0) * 1) + diag(rep(1, DCM.n)) ) / pA
               }
               pC.B         = (DCM$b)
               pC.C         = (DCM$c)
               pC.D         = (DCM$d)
               pC.transit   = array(dim=DCM.n, data=0) + exp(-6)
               pC.decay     = array(dim=DCM.n, data=0) + exp(-6)
               pC.epsilon   = array(dim=1, data=0) + exp(-6)
               pC = diag( c(as.vector(pC.A), 
                            as.vector(pC.B), 
                            as.vector(pC.C),
                            as.vector(pC.D),
                            as.vector(pC.transit),
                            as.vector(pC.decay),
                            as.vector(pC.epsilon)) )
             }
             else     # two state
             {
               pA = 16
               if ( length(dim(.Object@A)) == 2 )
               {
                 pC.A = ( ((.Object@A+32) != 0) * 1 ) / pA
               }
               pC.B         = (DCM$b) / 4
               pC.C         = (DCM$c) * 4
               pC.D         = (DCM$d) / 4
               pC.transit   = array(dim=DCM.n, data=0) + exp(-6)
               pC.decay     = array(dim=DCM.n, data=0) + exp(-6)
               pC.epsilon   = array(dim=1, data=0) + exp(-6)
               pC = diag( c(as.vector(pC.A), 
                            as.vector(pC.B), 
                            as.vector(pC.C),
                            as.vector(pC.D),
                            as.vector(pC.transit),
                            as.vector(pC.decay),
                            as.vector(pC.epsilon)) )
             }
             
             return(pC)
           }
)


setMethod( f          = "readParams",
           signature  = "Params",
           definition = function(.Object, Ep.mat){
             
             Ep = readMat( Ep.mat )
             Ep = Ep$Ep[,,1]
             
             DCM.n = dim(Ep$A)[1]
             
             .Object@A          = as.array( Ep$A )
             .Object@B          = as.array( Ep$B )
             .Object@C          = as.array( Ep$C )
             .Object@D          = as.array( Ep$D )
             .Object@transit    = array( Ep$transit, dim=DCM.n )
             .Object@decay      = array( Ep$decay, dim=DCM.n )
             .Object@epsilon    = array( Ep$epsilon, dim=1 )
             
             return (.Object)
           }
)


setMethod( f          = "pVect",
           signature  = "Params",
           definition = function(.Object, transpose=FALSE){
             
             if( transpose )
             {
               n = dim(.Object@A)[1]
               u = length(.Object@C) / n
               
               tB = NULL
               tD = NULL
               tA = NULL
               
               if ( length( dim(.Object@B) ) == 2 )
               {
                 tB = t(.Object@B)
               }
               else
               {
                 for ( i in 1:u )
                 {
                   tB = c( tB, t(.Object@B[,,i]) )              
                 }
               }
               
               if ( length(.Object@D) > 0 )
               {
                 for ( i in 1:n )
                 {
                   tD = c( tD, t(.Object@D[,,i]) )              
                 }
               }
               
               if ( length(dim(.Object@A)) > 2 )
               {
                 for ( i in 1:dim(.Object@A)[3] )
                 {
                   tA = c( tA, t(.Object@A[,,i]) )              
                 }
               }
               else
               {
                tA = t(.Object@A)
               }
               
               return( c(tA, tB, t(.Object@C), tD, .Object@transit, .Object@decay, .Object@epsilon) )
             }
             
             return( c(.Object@A, .Object@B, .Object@C, .Object@D, .Object@transit, .Object@decay, .Object@epsilon) )
             
           }
)


setMethod( f = "pUnvect",
           signature = "Params",
           definition = function(.Object, v, nU, nY){
             
             sB = nY * nY + 1
             sC = (nU + 1) * (nY * nY) + 1
             sD = sC + nY * nU
             sH = length(v) - 2*nY
             
             .Object@A        = array( v[1:(sB-1)], dim=c(nY,nY) )
             .Object@B        = array( v[sB:(sC-1)], dim=c(nY,nY,nU) )
             .Object@C        = array( v[sC:(sD-1)], dim=c(nY,nU))
             if ( (sH - sD) < nY * nY )
             {
              .Object@D       = array( 0, dim=c(nY,nY,0))
             }
             else
             {
              .Object@D       = array( v[sD:(sH-1)], dim=c(nY,nY,nY) )
              #.Object@D       = array( 0, dim=c(nY,nY,0))
             }
             .Object@transit  = array( v[sH:(sH+nY-1)] )
             .Object@decay    = array( v[(sH+nY):(sH+2*nY-1)] )
             .Object@epsilon  = array( v[(sH+2*nY)] )
             
             return(.Object)
             
           }
)
