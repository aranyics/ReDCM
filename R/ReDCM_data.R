

#################################################################################################
# source:       ReDCM_data.R
#------------------------------------------------------------------------------------------------
# type:         Class definition
#------------------------------------------------------------------------------------------------
# title:        BOLD-signal class definition for fMRI
#------------------------------------------------------------------------------------------------
# description:  Contains information of measured responses for fMRI.
#
#               Slots:  y       -   responses / output timeseries
#                       name    -   region names
#                       dt      -   sampling interval (TR)
#                       X0      -   confounds
#                       Q       -   array of precision components
#                       scale   -   BOLD-signal scaling
#
#               The measured BOLD-signal is stored in the y matrix. Scaling is applied to y to
#               enforce a maximum change of 4%.
#
#               X0 contains a constant term, and drift components given by a discrete cosine set
#               to simulate low frequency drift.
#------------------------------------------------------------------------------------------------
# based on:     DCM.Y structure
#               DCM12 (http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
#------------------------------------------------------------------------------------------------
# TODO:         - try wavelet filter in place of discrete cosine set
#               - use appropriate R matrix representation
#################################################################################################



.Data <- setClass("Data",
                  
                  slots = list(y       = "matrix",
                               name    = "character",
                               dt      = "numeric",
                               Q       = "list",
                               X0      = "matrix",
                               scale   = "numeric"),
                  
                  prototype = list(y       = matrix(),
                                   name    = character(),
                                   dt      = double(),
                                   Q       = list(),
                                   X0      = matrix(),
                                   scale   = double() )
                  
)



setMethod( f          = "setData",
           signature  = "Data",
           definition = function(.Object, DCM){
            
            DCM.n = length(DCM$Y[,,1]$name)
            DCM.v = dim(DCM$Y[,,1]$y)[1]
            
            Q = list()
            for ( i in 1:length(DCM$Y[,,1]$Q) )
            {
              Q = append(Q, list(DCM$Y[,,1]$Q[[i]][[1]]))
            }
            
            for ( i in 1:DCM.n )
            {
              .Object@name[i] = DCM$Y[,,1]$name[[i]][[1]][1]    # region name
            }
            .Object@y     =	(DCM$Y[,,1]$y)                      # responses / output time series
            if ( is.na(.Object@X0)[1] )
            {
              .Object@X0  = as.matrix(rep(1,DCM.v))
            }
            if ( is.null(DCM$Y[,,1]$X0) )                       #2018.09.13
            {
              .Object@X0  = as.matrix(rep(1,DCM.v))
            }
            else
            {
              .Object@X0    = (DCM$Y[,,1]$X0)                   # confounds
            }
            .Object@dt    = DCM$Y[,,1]$dt[1]                    # TR
            .Object@Q     = Q                                   # array of precision components
            
            cat( '- measured bold data ----------------\n' )
            cat( '----original:   ', .Object@y[1:10], '\n')
            .Object = detrend(.Object)
            
            cat( '----detrended:  ', .Object@y[1:10], '\n')
            .Object = setScale(.Object)
            
            cat( '----scaled:     ', .Object@y[1:10], '\n')
            return(.Object)
            
          }
)


setMethod( f          = "detrend",
           signature  = "Data",
           definition = function(.Object){
             
             if ( is.na(.Object@y)[1] )
             {
               cat("Err: Uninitialized data")
               return(0)
             }
             
             r = length(.Object@y[,1])
             
             .Object@y = .Object@y - rep(1, r) %*% t(colMeans(.Object@y))
             
             return(.Object)
           }
)


setMethod( f          = "setScale",
           signature  = "Data",
           definition = function(.Object){
             
             
             scale = max(.Object@y) - min(.Object@y)
             scale = 4 / max(scale, 4)
             
             
             .Object@y = .Object@y * scale
             .Object@scale = scale
             
             return(.Object)
           }
)
