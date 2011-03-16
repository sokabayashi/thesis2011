require( rcdd )

#######################################
## kappa
##
#######################################
kappa <- function( theta, t.y, nu )
{
	sum( exp( t.y %*% theta * nu ) )
}

#######################################
## kappa.a
##
#######################################
kappa.a <- function( theta, t.y, nu )
{
	a <- max( t.y %*% theta )
	sum( exp( t.y %*% theta - a ) * nu )
}

#######################################
## loglike
##
#######################################
loglike <- function( theta, y.obs, Y.all.mat )
{
	t.y <- Y.all.mat[ ,c(1,2)]
	nu  <- Y.all.mat[ ,3]
	a 	<- max( t.y %*% theta )
	kappa_a <- kappa.a(theta, t.y, nu )
	ll <- y.obs %*% theta - a - log( kappa_a )
	ll
}


#######################################
## loglike.alpha
##
#######################################
loglike.alpha <- function( alph, theta, p, y.obs, Y.all.mat )
{
	t.y <- Y.all.mat[ ,c(1,2)]
	nu  <- Y.all.mat[ ,3]
	theta.new <- theta + alph*p
	a 	<- max( t.y %*% theta.new )
	kappa_a <- kappa.a(theta.new, t.y, nu )
	ll <- y.obs %*% theta.new - a - log( kappa_a )
	ll
}


#######################################
## loglike.constrained
##
#######################################
loglike.constrained <- function( theta, y.obs, Y.all.mat )
{
	theta[1] <- 0
	loglike( theta, y.obs, Y.all.mat )
}

#######################################
## logl
##	for use with trust.
#######################################
logl <- function( theta, y.obs, Y.all.mat )
{
	t.y <- Y.all.mat[,c(1,2)]
	nu <- Y.all.mat[,3]
	a <- max( t.y %*% theta )
	e.ttheta.a <- as.vector( exp( t.y %*% theta - a ) )
	kappa.a <- sum(  e.ttheta.a * nu )
	value <- y.obs %*% theta - a - log( kappa.a )
	value <- as.double( value )

	mu <- 1 / kappa.a * colSums( t.y * e.ttheta.a * nu ) 
	grad <- as.matrix( y.obs - mu )

	t.mu <- t( t(t.y) - mu )
	hess <- t( apply( t.mu, 1, function(x) x %*% t(x) ) ) * e.ttheta.a * nu
	hess <- -matrix( 1/kappa.a * colSums( hess ), ncol = 2 )
	
	return( list( value = value, gradient = grad, hessian = hess ) )
}

#######################################
## logl.constrained
##	for use with trust.
#######################################
logl.constrained <- function( theta, y.obs, Y.all.mat )
{
	theta[1] <- 0
	t.y <- Y.all.mat[,c(1,2)]
	nu <- Y.all.mat[,3]
	a <- max( t.y %*% theta )
	e.ttheta.a <- as.vector( exp( t.y %*% theta - a ) )
	kappa.a <- sum(  e.ttheta.a * nu )
	value <- y.obs %*% theta - a - log( kappa.a )
	value <- as.double( value )

	mu <- 1 / kappa.a * colSums( t.y * e.ttheta.a * nu ) 
	grad <- as.matrix( y.obs - mu )

	t.mu <- t( t(t.y) - mu )
	hess <- t( apply( t.mu, 1, function(x) x %*% t(x) ) ) * e.ttheta.a * nu
	hess <- -matrix( 1/kappa.a * colSums( hess ), ncol = 2 )
	
	return( list( value = value, gradient = grad, hessian = hess ) )
}

#######################################
## pair2g
##
#######################################
pair2g <- function( pair, Y.all.mat )  # must be the full Y.all.mat
{
#	which( Y.all.mat$edges==pair[1] & Y.all.mat$triangles==pair[2] )
	which( Y.all.mat[,1]==pair[1] & Y.all.mat[,2]==pair[2] )
}

#######################################
## t.likelihood
##
#######################################
t.likelihood <- function( theta, pair, Y.all.mat )
{
	t.y <- Y.all.mat[ ,c(1,2)]
	nu  <- Y.all.mat[ ,3]
	freq <- nu[ pair2g( pair, Y.all.mat ) ]
	a <- max( t.y %*% theta )
	kappa_a <- kappa.a( theta, t.y, nu )

	value <- 1/kappa_a * exp( pair %*% theta - a ) * freq
	value
}


#######################################
## t.likelihood.tset
##  apply t.likelihood to a set of t values
#######################################
t.likelihood.tset <- function( theta, tset, Y.all.mat, pt.flag = FALSE )
{
	if( is.matrix( tset ) )
		probs <- apply( tset, 1, function(ty) t.likelihood( theta, ty, Y.all.mat ) )
	else
		probs <- t.likelihood( theta, tset, Y.all.mat )

	if( pt.flag )	
		return( probs )
	else
		return( sum( probs ) )
}

#######################################
## theta2mu
##
#######################################
theta2mu <- function( theta, Y.all.mat )
{
	t.y <- Y.all.mat[ ,c(1,2), drop=FALSE]
	if( is.matrix( t.y ) )
	{
		prob <- apply( t.y, 1, function(pair) t.likelihood( theta, pair, Y.all.mat ) )
 		return( colSums( t.y * prob ) )
	} else
	{
		prob <- t.likelihood( theta, t.y, Y.all.mat )
 		return( t.y * prob )
	}
}


#######################################
## Hrep
##
#######################################
Hrep <- function( Y.all.mat )
{
	t.y <- as.matrix( Y.all.mat[ ,c(1,2), drop=FALSE ] )
	m <- cbind( 0, 1, t.y )
	
	scdd( d2q(m), representation = "V" )	
}

#######################################
## ty2Hrep
##
#######################################
ty2Hrep <- function( t.y )
{
	m <- cbind( 0, 1, t.y )	
	scdd( d2q(m), representation = "V" )	
}

#######################################
## get.tcone
## 	return H-rep of tangent cone at t.pt
#######################################
get.tcone <- function( t.pt, t.y )
{
	t.y <- rbind( t.y, t.pt )
	out <- ty2Hrep( t.y )
	b <- out$output[, 2]
	v <- out$output[, -c(1,2) ]
	a <- qneg(v) 

	axb <- qmatmult(a, as.matrix(t.pt) )
	axb <- sweep( axb, 1, b, FUN = qmq ) # subtract
	active <- axb==0
	if( all( active==FALSE) ) return( out$output ) 

	tcone <- out$output[ active, ,drop=FALSE]
	tcone[ ,2] <- 0

	return( tcone )
}



#######################################
## get.face.pts
## given a pt that is on face, return other pts from Y.all.freq that are on that same face
#######################################
get.face.pts <- function( t.pt, Y.all.freq )
{
	# find the Hrep inequalities that are active for t.pt
	# 	put the coefs in a.active, b.active
	out <- Hrep( Y.all.freq )
	b <- out$output[, 2]
	v <- out$output[, -c(1,2), drop=FALSE ]
	a <- qneg(v) 

	axb <- qmatmult(a, as.matrix(t.pt) )
	axb <- sweep( axb, 1, b, FUN = qmq ) # subtract
	if( max(qsign(axb)) > 0 ) return( NULL ) # pt is outside hull

	active <- axb==0
	if( all( active==FALSE) ) return( Y.all.freq ) # pt is on interior of hull
	a.active <- a[active, , drop=FALSE] 
	b.active <- b[active]
	
	# find pts in space that meet this new subset of inequalities
	axb.active <- qmatmult(a.active, t(as.matrix(Y.all.freq[,c(1,2)])) )
	axb.active <- sweep( axb.active, 1, b.active, FUN = qmq ) # subtract
	active.pts <- apply( axb.active, 2, function(foo) max(qsign(foo)) )+1

	# construct rays for all pts relative to t.pt
#	rays <-  sweep( Y.all.freq[ active.pts==1, c(1,2) ], 2, t.pt, FUN = qmq ) 
	rays <-t( Y.all.freq[ active.pts==1, c(1,2), drop=FALSE ] ) - t.pt
	rays <- t( rays )
	
	# find the rays that satisfy linearity
	Vnew <- as.matrix( cbind( 0, 0, rays ) ) # 2nd column must be zero--these are RAYS!
	Vnew.dat <- as.data.frame( Vnew )		 # need the data frame to preserve the original rownames
	lin.active 	<- linearity( Vnew, rep = "V" )
	indices 	<- rownames( Vnew.dat[ lin.active, ] )

	return( Y.all.freq[indices,] )
}


#######################################
## get.ncone
##		
#######################################
get.ncone <- function( t.pt, t.y )
{
	tcone <- get.tcone( t.pt, t.y )
	if( is.null(tcone) ) return( NULL )
	
	ncone <- qneg( tcone )
	ncone[,1] <- tcone[ ,1]

	if( nrow( ncone ) > 1 )		
	{
		ncone <- redundant( ncone, rep = "V" )
		ncone <- ncone$output
	}

	return( ncone )
}


#######################################
## get.network.sample
##  returns a sample of INDICES to Y.all.mat
##  these indices correspond to the graphs
#######################################
get.network.sample <- function( size, theta, Y.all.mat )
{
	n <- nrow( Y.all.mat )
	t.y <- Y.all.mat[ ,c(1,2)]
	prob <- apply( t.y, 1, function(pair) t.likelihood( theta, pair, Y.all.mat ) )
	
	sample.int( n, size, replace = TRUE, prob = prob )
}


#######################################
## count.pair
##   
#######################################
count.pair <- function( pair, sample.ty )
{
	sum( sample.ty[,1] == pair[1] & sample.ty[,2] == pair[2] )
}


#######################################
## nabla.loglike.t.p
##   a non-MCMC version.
#######################################
nabla.loglike.t.p <- function( alpha, theta, p, y.obs, Y.all.mat, adj=0 )
{
	nabla.loglike <- y.obs - theta2mu( theta + alpha*p, Y.all.mat )
	t(nabla.loglike) %*% p - adj
}


#######################################
## closest.ascent.dir
##   
#######################################
closest.ascent.dir <- function( p.k, del.ll )
{
	p.dim <- length( p.k )
	del.ll.perp <- matrix( rep( 1 / del.ll, p.dim ), nrow = p.dim, byrow = TRUE )

	diag( del.ll.perp ) <- -0.99*diag( del.ll.perp )
	dotprods <- del.ll.perp %*% p.k

	if ( dotprods[1] > dotprods[2] )
	{ 
		p.current <- del.ll.perp[1,]
	} else p.current <- del.ll.perp[2,]
	
	return( p.current )
}


#######################################
## in.exterior
##    use strongly separating hyperplane
##    to see if a pt q is outside the hull defined by ty.hull
#######################################
in.exterior <- function( q, ty.hull )
{
	hrep <- cbind( 0, 0, 1, qneg(ty.hull) )	# s.t. 1st ineq
	hrep <- rbind( hrep, c(0, 1, 1, -q ) )	# 2nd (artifical) ineq
	f <- lpcdd( hrep, d2q( c(-1, q) ), minimize = FALSE )$optimal.value
	if( f > 0 ) return( TRUE )
	else return( FALSE )
}



