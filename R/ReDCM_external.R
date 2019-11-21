

#################################################################################################
# source:       ReDCM_external.R
#------------------------------------------------------------------------------------------------
# type:         Class definition
#------------------------------------------------------------------------------------------------
# title:        External input class definition
#------------------------------------------------------------------------------------------------
# description:  Contains information of external input and design.
#
#               Slots:  u       -   exogenous inputs
#                       name    -   input names
#                       dt      -   sampling interval
#------------------------------------------------------------------------------------------------
# based on:     DCM.U structure
#               DCM12 (http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
#------------------------------------------------------------------------------------------------
# TODO:         - use appropriate R matrix representation
#               - rename 'External'
#################################################################################################



.External <- setClass("External",
                  
                  slots = list(#u       = "dgCMatrix",
                               u       = "Matrix",
                               name    = "character",
                               dt      = "numeric"),
                  
                  prototype = list(#u       = new("dgCMatrix"),
                                   u       = Matrix(),
                                   name    = character(),
                                   dt      = double() )
                  
)


setMethod( f          = "setData",
           signature  = "External",
           definition = function(.Object, DCM){
             
             DCM.n = length(DCM$U[,,1]$name)
             DCM.options = DCM$options[,,1]
             
             for (i in 1:DCM.n)
             {
               .Object@name[i] = DCM$U[,,1]$name[[i]][[1]][1]    # region name
             }
             .Object@u     = Matrix(DCM$U[,,1]$u)                # responses / output time series #2018.09.13
             .Object@dt    = DCM$U[,,1]$dt[1]                    # sampling interval
             
             
             cat("microtime bins: ", dim(.Object@u)[1], "\t(", dim(.Object@u)[1]/16, " scans by default)", '\n\n')
             
             # centre input
             if ( FALSE )
             #if ( DCM.options$centre )
             {
               .Object = detrend(.Object)
               cat("centering input")
             }
             
             return(.Object)
             
           }
)


setMethod( f          = "detrend",
           signature  = "External",
           definition = function(.Object){
             
             if ( is.na(.Object@u)[1] )
             {
               cat("Err: Uninitialized data")
               return(0)
             }
             
             r = length(.Object@u[,1])
             
             .Object@u = .Object@u - rep(1, r) %*% t(colMeans(.Object@u))
             .Object@u = as( Matrix(.Object@u), "dgCMatrix" )
             
             return(.Object)
           }
)


setMethod( f          = "extSpike",
           signature  = "External",
           definition = function(.Object, freq=0.20125){
             
             len = 20.125 #sec
             
             .Object@name = "spike"
             .Object@dt = freq
             .Object@u = Matrix(0, nrow=len/freq, ncol=1)
             .Object@u = .Object@u[2] = 1
             
             return(.Object)
           }
)


setMethod( f          = "extCorr",
           signature  = "External",
           definition = function(.Object, M){
             
             u = dim(.Object@u)[1]
             period = M@Ydt / .Object@dt
             k = dim(.Object@u)[1] - M@ns * period
             #k = M@ns * period
             
             if ( k > 0 )
             #if ( (u-k) > 0 )
             {
                .Object@u = as( Matrix(.Object@u[-(1:k),]), "dgCMatrix" )
               #.Object@u = as( Matrix(.Object@u[-((k+1):u),]), "dgCMatrix" )
             }
             
             return(.Object)
           }
)
