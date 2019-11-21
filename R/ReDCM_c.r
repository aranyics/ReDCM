
# C interface for utility functions =========================================#
# ===========================================================================#

c_blas_dgemm = function ( A, B )
{

  mx_dims = c(dim(as.matrix(A)), dim(as.matrix(B)))
  
  Y = array(0, dim=mx_dims[1]*mx_dims[4])
  
  res = .C( "_r_cblas_dgemm",
            A = as.double( t(A) ),
            B = as.double( t(B) ),
            d = as.integer(mx_dims),
            Y = as.double(Y),
            PACKAGE = 'libredcmc' )
  
  res = t( array( res$Y, dim=c(mx_dims[4], mx_dims[1]) ) )
  
  return(res)

}







# C interface for model integration =========================================#
# ===========================================================================#

c_call_integrate = function ( func, x.par, u.par, P.par, M.par, U.par, H.par, Y.out, X.out=NULL )
{

  if ( is.null(X.out) )
  {
    res = .C( func, x.par = as.double( x.par ),
              u.par = as.double( u.par ),
              P.par = as.double( P.par ),
              M.par = as.integer( M.par ),
              U.par = as.double( U.par ),
              H.par = as.double( H.par ),
              Y.out = as.double( Y.out ),
              PACKAGE = 'libredcmc' )
  
    return ( res$Y.out )
  }
  else
  {
    res = .C( func, x.par = as.double( x.par ),
              u.par = as.double( u.par ),
              P.par = as.double( P.par ),
              M.par = as.integer( M.par ),
              U.par = as.double( U.par ),
              H.par = as.double( H.par ),
              Y.out = as.double( Y.out ),
              X.out = as.double( X.out ),
              PACKAGE = 'libredcmc' )
    
    return ( list( res$Y.out, res$X.out ) )
  }

}


c_int_det = function( P, M, U )
{
  
  bins = M@ns * ( M@Ydt / U@dt )
  
  P.par = pVect( P, transpose=TRUE)
  M.par = c(M@m, M@l, M@options$nonlinear, M@options$two.state, 0, 0, 0)
  x.par = c(M@x)
  #Ext.par = c(t( as.matrix(U@u) ))
  Ext.par = c( as.matrix(U@u) )
  u.par = c( M@Ydt, U@dt, bins, M@ns, M@delays )
  H.par = c( M@H )
  
  Y.out = array(0, dim=M@l*M@ns)
  
  res_c = c_call_integrate( M@IS, x.par, u.par, P.par, M.par, Ext.par, H.par, Y.out )
  res_c = t( array( res_c, dim=c(M@l, M@ns) ) )
  
  return ( res_c )
  
}


c_int_det_hemodyn = function( P, M, U )
{
  
  if ( M@l > 1 )
    P.par = pVect( P, transpose=TRUE)
  else
    P.par = pVect( P )
  M.par = c(M@m, M@l, 0, 0, 0, 0, 0)
  x.par = c(M@x)
  #Ext.par = c(t( as.matrix(U@u) ))
  Ext.par = c( as.matrix(U@u) )
  u.par = c( M@Ydt, U@dt, dim(U@u)[1], M@ns, M@delays )
  H.par = c( M@H )
  
  Y.out = array(0, dim=M@l*M@ns)
  X.out = array(0, dim=M@n*M@ns)
  
  res_c = c_call_integrate( '_r_int_det', x.par, u.par, P.par, M.par, Ext.par, H.par, Y.out, X.out )
  Y.out = t( array( res_c[[1]], dim=c(M@l, M@ns) ) )
  X.out = t( array( res_c[[2]], dim=c(M@n, M@ns) ) )
  
  return ( list( Y.out, X.out ) )
  
}


c_call_intdiff = function ( func, x.par, u.par, Pr.par, M.par, U.par, H.par, V.par, dY.out )
{

  res = .C( func, x.par = as.double( x.par ),
            u.par = as.double( u.par ),
            Pr.par = as.double( Pr.par ), 
            M.par = as.integer( M.par ),
            U.par = as.double( U.par ),
            H.par = as.double( H.par ),
            V.par = as.double( V.par ),
            dY.out = as.double( dY.out ),
            PACKAGE = 'libredcmc' )
  
  return ( res$dY.out )

}


c_int_dfdp = function( Pr, M, U, V )
{
  
  bins = M@ns * ( M@Ydt / U@dt )
  
  Pr.par = c(Pr)
  M.par = c(M@m, M@l, M@options$nonlinear, M@options$two.state, 0, 0, 0, dim(V))
  x.par = c(M@x)
  #Ext.par = c(t( as.matrix(U@u) ))
  Ext.par = c( as.matrix(U@u) )
  u.par = c( M@Ydt, U@dt, bins, M@ns, M@delays )
  H.par = c( M@H )
  V.par = c( t( V ) )
  
  Y.out = array(0, dim=M@l*M@ns*length(Pr))
  
  res_c = c_call_intdiff( '_r_int_dfdp', x.par, u.par, Pr.par, M.par, Ext.par, H.par, V.par, Y.out )
  #res_c = array( res_c, dim=c(M@l*M@ns, length(Pr)) )
  res_c = matrix( res_c, M@l*M@ns, length(Pr), byrow=TRUE )
  res_c = t( aaply(res_c, 2, function(x) { c(t(array(x, dim=c( M@l, M@ns )))) }) )
  
  return ( res_c )
  
}



# C interface for bireduce ==================================================#
# ===========================================================================#

c_call_reduce = function( func, x.par, H.par, M.par, P.par, M0.par, M1.par, L1.par, L2.par, out.par, M2.par=NULL )
{

  if ( M.par[3] )
  {
    res = .C( func, x.par = as.double( x.par ),
              H.par = as.double( H.par ),
              M.par = as.integer( M.par ),
              P.par = as.double( P.par ),
              M0.par = as.double( M0.par ),
              M1.par = as.double( M1.par ),
              M2.par = as.double( M2.par ),
              L1.par = as.double( L1.par ),
              L2.par = as.double( L2.par ),
              out.par = as.integer( out.par ),
              PACKAGE = 'libredcmc' )
    
    return ( list( res$M0.par, res$M1.par, res$M2.par, res$L1.par, res$L2.par ) )
  }
  else
  {
    res = .C( func, x.par = as.double( x.par ),
              H.par = as.double( H.par ),
              M.par = as.integer( M.par ),
              P.par = as.double( P.par ),
              M0.par = as.double( M0.par ),
              M1.par = as.double( M1.par ),
              L1.par = as.double( L1.par ),
              L2.par = as.double( L2.par ),
              out.par = as.integer( out.par ),
              PACKAGE = 'libredcmc' )
    
    return ( list( res$M0.par, res$M1.par, res$L1.par, res$L2.par ) )
  }
  

}


c_bireduce = function( M, P, out=0 )
{

  P.par = pVect( P, transpose=TRUE )
  M.par = c(M@m, M@l, M@options$nonlinear, M@options$two.state, 0, 0, 0)
  x.par = c(M@x)
  H.par = c(M@H)
  out.par = out
  
  mx.size = (M@n+1) * (M@n+1)
  M0 = array(0, dim=mx.size)
  M1 = array(0, dim=mx.size * M@m)
  L1 = array(0, dim=mx.size)
  L2 = array(0, dim=mx.size * M@l)
  
  res_c = c_call_reduce( '_r_bireduce', x.par, H.par, M.par, P.par, M0, M1, L1, L2, out.par )
  
  M0 = array( res_c[[1]], dim=c(M@n+1, M@n+1) )
  M1 = array( res_c[[2]], dim=c(M@n+1, M@n+1, M@m) )
  L1 = array( res_c[[3]], dim=c(M@n+1, M@n+1) )[1:M@l,]
  L2 = array( res_c[[4]], dim=c(M@n+1, M@n+1, M@l) )
  
  return ( list( M0, M1, L1, L2 ) )

}


c_soreduce = function( M, P, out=0 )
{
  
  P.par = pVect( P, transpose=TRUE )
  M.par = c(M@m, M@l, M@options$nonlinear, M@options$two.state, 0, 0, 0)
  x.par = c(M@x)
  H.par = c(M@H)
  out.par = out
  
  mx.size = (M@n+1) * (M@n+1)
  M0 = array(0, dim=mx.size)
  M1 = array(0, dim=mx.size * M@m)
  M2 = array(0, dim=mx.size * M@n)
  L1 = array(0, dim=mx.size)
  L2 = array(0, dim=mx.size * M@l)
  
  res_c = c_call_reduce( '_r_soreduce', x.par, H.par, M.par, P.par, M0, M1, M2, L1, L2, out.par )
  
  M0 = array( res_c[[1]], dim=c(M@n+1, M@n+1) )
  M1 = array( res_c[[2]], dim=c(M@n+1, M@n+1, M@m) )
  M2 = array( res_c[[3]], dim=c(M@n+1, M@n+1, M@n) )
  L1 = array( res_c[[4]], dim=c(M@n+1, M@n+1) )[1:M@l,]
  L2 = array( res_c[[5]], dim=c(M@n+1, M@n+1, M@l) )
  
  return ( list( M0, M1, L1, L2, M2 ) )
  
}



# C interface for fx_fmri and gx_fmri calls =================================#
# ===========================================================================#

c_call_fx = function( func, x.par, u.par, P.par, M.par, H.par, f.out )
{
  
  res = .C( func, x.par=as.double(x.par),
              u.par=as.double(u.par),
              P.par=as.double(P.par),
              M.par=as.integer(M.par),
              H.par=as.double(H.par),
              f.out=as.double(f.out),
              PACKAGE = 'libredcmc' )
  
  return ( res$f.out )
  
}

c_call_gx = function( func, x.par, u.par, P.par, M.par, f.out )
{
  
  res = .C( func, x.par=as.double(x.par),
            u.par=as.double(u.par),
            P.par=as.double(P.par),
            M.par=as.integer(M.par),
            f.out=as.double(f.out),
            PACKAGE = 'libredcmc' )
  
  return ( res$f.out )
  
}


c_fx_fmri = function( x, u, P, M )
{
  
  P.par = pVect( P, transpose=TRUE )
  M.par = c(M@m, M@l, M@options$nonlinear, M@options$two.state, 0, 0, 0)
  H.par = c( M@H )
  
  f.out = array( 0, dim=M@n )
  res_c = c_call_fx( '_r_fx_fmri', x, u, P.par, M.par, H.par, f.out )
  
  return ( res_c )
  
}


c_dfdx = function( x, u, P, M )
{
  
  P.par = pVect( P, transpose=TRUE )
  M.par = c(M@m, M@l, M@options$nonlinear, M@options$two.state, 0, 0, 0)
  H.par = c( M@H )
  
  f.out = array( 0, dim=M@n*M@n )
  res_c = c_call_fx( '_r_dfdx', x, u, P.par, M.par, H.par, f.out )
  res_c = t(array( res_c, dim=c(M@n, M@n) ))
  
  return ( res_c )
  
}


c_dfdu = function( x, u, P, M )
{
  
  P.par = pVect( P, transpose=TRUE )
  M.par = c(M@m, M@l, M@options$nonlinear, M@options$two.state, 0, 0, 0)
  H.par = c( M@H )
  
  f.out = array( 0, dim=M@n*M@m )
  res_c = c_call_fx( '_r_dfdu', x, u, P.par, M.par, H.par, f.out )
  res_c = t(array( res_c, dim=c(M@m, M@n) ))
  
  return ( res_c )
  
}


c_dfdxdu = function( x, u, P, M )
{
  
  P.par = pVect( P, transpose=TRUE )
  M.par = c(M@m, M@l, M@options$nonlinear, M@options$two.state, 0, 0, 0)
  H.par = c( M@H )
  
  f.out = array( 0, dim=M@n*M@n*M@m )
  res_c = c_call_fx( '_r_dfdxdu', x, u, P.par, M.par, H.par, f.out )
  res_c = t(array( res_c, dim=c(M@m, M@n*M@n) ))
  res_c = array( res_c, dim=c(M@n, M@n, M@m) )
  res_c = array( apply(res_c, 3, t), dim=c(M@n, M@n, M@m) )
  
  return ( res_c )
  
}


c_dfdxdx = function( x, u, P, M )
{
  
  P.par = pVect( P, transpose=TRUE )
  M.par = c(M@m, M@l, M@options$nonlinear, M@options$two.state, 0, 0, 0)
  H.par = c( M@H )
  
  f.out = array( 0, dim=M@n*M@n*M@n )
  res_c = c_call_fx( '_r_dfdxdx', x, u, P.par, M.par, H.par, f.out )
  res_c = t(array( res_c, dim=c(M@n, M@n*M@n) ))
  res_c = array( res_c, dim=c(M@n, M@n, M@n) )
  res_c = array( apply(res_c, 3, t), dim=c(M@n, M@n, M@n) )
  
  return ( res_c )
  
}


c_gx_fmri = function( x, u, P, M )
{
  
  P.par = P@epsilon
  M.par = c(M@l, M@options$two.state, M@n)
  
  f.out = array( 0, dim=M@l )
  res_c = c_call_gx( '_r_gx_fmri', x, u, P.par, M.par, f.out )
  
  return ( res_c )
  
}


c_dgdx = function( x, u, P, M )
{
  
  P.par = P@epsilon
  M.par = c(M@l, M@options$two.state, M@n)
  
  f.out = array( 0, dim=M@l*M@n )
  res_c = c_call_gx( '_r_dgdx', x, u, P.par, M.par, f.out )
  res_c = t(array( res_c, dim=c(M@n, M@l) ))
  
  return ( res_c )
  
}


c_dgdxdx = function( x, u, P, M )
{
  
  P.par = P@epsilon
  M.par = c(M@l, M@options$two.state, M@n)
  
  f.out = array( 0, dim=M@l*M@n*M@n )
  res_c = c_call_gx( '_r_dgdxdx', x, u, P.par, M.par, f.out )
  res_c = t(array( res_c, dim=c(M@n, M@l*M@n) ))
  res_c = array( res_c, dim=c( M@l, M@n, M@n ) )
  
  return ( res_c )
  
}
