

#################################################################################################
# source:       ReDCM_generics.R
#------------------------------------------------------------------------------------------------
# type:         -
#------------------------------------------------------------------------------------------------
# title:        General constant and name definitions, and source including
#------------------------------------------------------------------------------------------------
# description:  S4 generics are listed here, and constant (prior) values
#------------------------------------------------------------------------------------------------
# based on:     -
#------------------------------------------------------------------------------------------------
# TODO:         -
#################################################################################################


# GENERICS
#################################################################################################

setGeneric("setData",
           
           def = function(.Object, DCM){
             standardGeneric("setData")
           }
)

setGeneric("extSpike",
           
           def = function(.Object, freq=0.20125){
             standardGeneric("extSpike")
           }
)

setGeneric("extCorr",
           
           def = function(.Object, M){
             standardGeneric("extCorr")
           }
)

setGeneric("detrend",
           
           def = function(.Object){
             standardGeneric("detrend")
           }
)

setGeneric("setScale",
           
           def = function(.Object){
             standardGeneric("setScale")
           }
)

setGeneric("setupModel",
           
           def = function(.Object, DCM, Ext, H=NULL){
             standardGeneric("setupModel")
           }
)

setGeneric("extendModel",
           
           def = function(.Object, ...){
             standardGeneric("extendModel")
           }
)

setGeneric("setPriors",
           
           def = function(.Object, DCM){
             standardGeneric("setPriors")
           }
)

setGeneric("setPriorCov",
           
           def = function(.Object, DCM, pA){
             standardGeneric("setPriorCov")
           }
)

setGeneric("setHemodyn",
           
           def = function(.Object, H){
             standardGeneric("setHemodyn")
           }
)

setGeneric("readParams",
           
           def = function(.Object, Ep.mat){
             standardGeneric("readParams")
           }
)

setGeneric("pVect", 
           
           def = function(.Object, transpose=FALSE){
             standardGeneric("pVect")
           }
)

setGeneric("pUnvect",
           
           def = function(.Object, v, nU, nY){
             standardGeneric("pUnvect")
           }
)

setGeneric("setupEstimates",
           def = function(.Object,
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
             standardGeneric("setupEstimates")
           }
)

setGeneric("initEstimates",
           
           def = function(.Object, DCM){
             standardGeneric("initEstimates")
           }
)

setGeneric("evidence",
           
           def = function(.Object){
             standardGeneric("evidence")
           }
)


# CONSTANTS
#################################################################################################

maxIter = 32

