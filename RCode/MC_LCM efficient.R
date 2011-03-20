##########################################################################
## last touched 3/16/11
library( rcdd    )
library( trust   )
library( ggplot2 )
library( lattice )
source( "~/Tako/THESIS/RCode/LIB_edgetri_fns.R" )
##########################################################################
# graph9.stats.RData
# 	stats.total
# 	Y.all.freq
load( "~/Tako/THESIS/RCode/graph9.stats.RData" )
Y.all.mat <- as.matrix( Y.all.freq )
t.y 	<- Y.all.mat[ ,c(1,2) ]
nu 		<- Y.all.mat[ ,3 ]
##########################################################################
y.obs 	  <- c( 21, 4 ) # NOT on boundary - problematic. (21,4), (23,12) also bad but (24,16) ok.
#y.obs 	  <- c( 31, 50 ) # boundary face
#y.obs 	  <- c( 27,27 ) # extreme pt
#y.obs 	  <- c( 29,47 ) # interior - no problem
#y.obs 	  <- c( 19, 0 ) # 
verbose.flag 	<- TRUE
MC.flag			<- TRUE
c 				<- 0.2
del.len.cutoff	<- 1e-2
samp.size 		<- 10000
seed 			<- 2
set.seed( seed )
max.k			<- 25
face.cutoff 	<- 0.60
##########################################################################
# ANSWERS - use knowledge of hull to look up answers which we can compare against
##########################################################################
face.freq 		<- get.face.pts( y.obs, Y.all.freq )  # face pts for y.obs, complete with frequencies.
face.freq.mat 	<- as.matrix( face.freq )
face.ty 		<- face.freq.mat[,c(1,2), drop=FALSE] # just the pts 
# Extract MLE if it exists
foo <- Y.all.freq[ pair2g( y.obs, Y.all.mat ), ]
theta.MLE <- c( foo$trust$theta.1, foo$trust$theta.2 )
# Find MLE in LCM
BN.optim.out <- optim( c(0,0), loglike.constrained, method = c( "BFGS" ),  control=list( maxit=200,fnscale=-1), y.obs = y.obs, Y.all.mat = face.freq.mat )
BN.optim.out
( theta.LCM <- BN.optim.out$par )
# what prob does this MLE give to each pt in face?
face.probs <- t.likelihood.tset( theta.LCM, face.ty, face.freq.mat, pt.flag = TRUE )
print.table( cbind( face.ty, face.probs ) )
# Find normal cone
ncone <- get.ncone( y.obs, t.y )
ncone <- q2d( ncone )
ncone <- ncone[ , -c(1,2), drop=FALSE ]
gdor.true 	<- colMeans( ncone ) # any vector off boundaries
doc			<- gdor.true

# just to make sure it's a doc
theta.new <- theta.LCM + doc
# Check that these give the same ll value
loglike( theta.LCM, y.obs, face.freq.mat )
loglike( theta.new, y.obs, face.freq.mat )
theta2mu( theta.LCM, face.freq.mat )

# from previous work for (31,50)
s.alpha <- 1.517749 # CI bound to 5% on these pts
( theta.alpha <- theta.LCM  + t( s.alpha %*% t(gdor.true) )  )
face.alpha.probs <- t.likelihood.tset( theta.alpha, face.ty, Y.all.mat, pt.flag=T ) 
cbind( face.ty, face.probs, face.alpha.probs )

##########################################################################
# PROBLEM SETUP
##########################################################################
theta.hist  	<- matrix( nrow = max.k+10, ncol = 2 )
theta.current 	<- c(0,0)
d 				<- length( theta.current )
theta.hist[1,]  <- theta.current
# initial sample
sample.g 		<- get.network.sample( samp.size, theta.current, Y.all.mat )
sample.mat 		<- Y.all.mat[ sample.g, ]
sample.ty  		<- sample.mat[ , c(1,2) ]
sample.Vrep 		<- unique( cbind( 0, 1, sample.ty ) )
sample.reduced 		<- redundant( d2q(sample.Vrep), rep = "V" )
sample.reduced.ty 	<- q2d( sample.reduced$output[,-c(1:2)] )
ty.hull 		<- sample.reduced.ty
# plot them!
plot( sample.ty, xlim = c(0,38), ylim= c(0,85) )
points( t(y.obs), pch = 19, col = "blue" )
polygon( sample.reduced.ty[ chull( sample.reduced.ty ), ], border = "blue", lty = "dotted" )

################################
# Gradient quantities 
if( MC.flag )
{	mu.current		<- colMeans( sample.ty )
} else {
	mu.current		<- theta2mu( theta.current, Y.all.mat )
}
p.current 		<- del.ll <- y.obs - mu.current
del.ll.len 		<- sqrt( del.ll %*% del.ll )
del.ll.tp.current 	<- del.ll.tp.old <- t( del.ll ) %*% p.current
	
k 		<- 2				
LCM.k 	<- 1
LCM.flag	<- FALSE
on.boundary <- FALSE
on.interior <- FALSE
gdor <- c("0","0")
####  main loop ##########################################
while( k <= max.k && del.ll.len > del.len.cutoff )
{
  ####  Step size search ##########################################
	# Find alpha that works - exact as possible, cheat and use uniroot (requires knowing exact distribution)
	if( LCM.flag )
	{	uni.mat <- temp.face.mat
	} else uni.mat <- Y.all.mat			# exact distribution, to pass on to uniroot

	####  error check for MC sample size  ##########################################
	del.ll.exact <- y.obs - theta2mu( theta.current, uni.mat )  
	if( del.ll.exact %*% p.current < 0 )
	{
		cat( "p: ", p.current, 
			"\n del.ll.exact: ", del.ll.exact, 
			"\n inner prod: ", p.current %*% del.ll.exact, "\n" )
		cat( "\n", k, ": The MC estimate and exact value for del.ll are in opposite direction\n, making
		it impossible to meet curvature condition.  This likely occurs because the MC sample size \n
		is too small.  On the flip side, this only happens when del.ll is very small as well,\n  
		so we are close to the maximum.  Exiting loop!\n" )
		break
	}	
	####  don't use uniroot when p is a GDOR: no root!  ##########################################
	if( all( p.current == q2d(gdor) ) ) # GDOR direction
	{	alpha.current <- 20 # 20 seems a large enough step size to go in GDOR direction 
	} else{
	####  uniroot - at last  ##########################################
		uni.exact.out <- uniroot( nabla.loglike.t.p, lower = .Machine$double.eps, upper = 10000000, tol = .Machine$double.eps, theta=theta.current, p=p.current, y.obs=y.obs, Y.all.mat=uni.mat ) 
		alpha.current <- uni.exact.out$root
	}
	alpha.current <- min( alpha.current, 20 )
	nabla.loglike.t.p( alpha.current, theta.current, p.current, y.obs, uni.mat ) # better be small!
	#############################################
	# Save previous iteration vars
	theta.old  		<- theta.current
	p.old			<- p.current
	del.ll.old		<- del.ll
	del.ll.len.old  <- del.ll.len
	del.ll.tp.old   <- del.ll.tp.current	
	#### update theta ##########################################
	theta.current 	<- theta.old + alpha.current * p.current
	theta.hist[k,]	<- theta.current
	cat( sprintf("\n\n %s: New step length: %11.6f del.ll.len old: %11.6f\n", k, alpha.current, del.ll.len  ) )
	cat( "theta current (updated)", theta.current, "\n" )
	
	#### New sample at this new theta, from ORIGINAL model ####################################
	sample.g.new 			<- get.network.sample( samp.size, theta.current, Y.all.mat )
	sample.ty.new  			<- Y.all.mat[ sample.g.new, c(1,2), drop = FALSE ]
	
	#### Importance weights #########################################################
	ty.diff <- sweep( sample.ty.new, 2, y.obs, FUN = "-" )
	importance.weights <- exp( ty.diff %*% (theta.current - theta.old ) )
	cat( stem(importance.weights ) )
	cat( summary( importance.weights ), "\n" )
	
	#### IF LCM, then we want the subsample restricted to face ######################
	if( LCM.flag == TRUE ) 
	{
		sample.ty.rays  <- sweep( sample.ty.new, 2, y.obs, FUN = "-" )
		sample.dp		<- sample.ty.rays %*% q2d(gdor)
		missed.pts		<- sample.ty.new[ sample.dp > 0, , drop = FALSE ]
		if( nrow( missed.pts ) > 0 )			
		{
			# Our sampler discovered a pt that has positive inner product
			# with GDOR.  means it is a point outside our empirical boundary.
			# need to reevaluate everything.  Restore ty.hull to old
			cat( "NEW POINT: ", missed.pts, "\n" )
			LCM.flag    <- FALSE
			on.boundary <- FALSE
			ty.hull <- ty.hull.old
		} else {
			sample.ty.old	<- sample.ty.new # keep a copy around -- if just for plotting
			temp.face.ty 	<- sample.ty.new[ sample.dp == 0, , drop = FALSE ]
			# assign to sample.ty.new
			sample.ty.new <- temp.face.ty # may have different size than default samp.size
		}
	}	
	
	#### Find nonredundant points in sample ##########################################
	sample.Vrep.new 		<- unique( cbind( 0, 1, sample.ty.new ) ) # leave sample.ty.new intact
	if( nrow( sample.Vrep.new ) > 1 )
	{
		sample.reduced.new 		<- redundant( d2q(sample.Vrep.new), rep = "V" )
		sample.reduced.ty.new 	<- q2d( sample.reduced.new$output[,-c(1:2), drop =FALSE] )
	} else sample.reduced.ty.new <- sample.Vrep.new[,-c(1:2), drop =FALSE]
	#### Construct CUMULATIVE convex hull ##########################################
	temp.hull 				<- unique( rbind( ty.hull, sample.reduced.ty.new ) )
	temp.hull.Vrep 		<- cbind( 0, 1, temp.hull )
	if( nrow( temp.hull.Vrep ) > 1 )
	{
		temp.hull.reduced 	<- redundant( d2q(temp.hull.Vrep), rep = "V" )
		temp.hull.ty			<- q2d( temp.hull.reduced$output[,-c(1:2), drop=FALSE] )
	} else temp.hull.ty 		<- temp.hull
	
 	ty.hull 				<- temp.hull.ty[ chull( temp.hull.ty ), , drop=FALSE]
 	ty.hull.extreme.count  	<- nrow( ty.hull )
#	temp.ncone <- get.ncone( y.obs, ty.hull )

	#### Plot hulls ##########################################
	par( mfrow=c(1,2) )
	if( LCM.flag == FALSE )
	{
		plot( sample.ty.new, xlim = c(0,38), ylim= c(0,85) ); points( t(y.obs), pch = 19, col = "blue" )
		polygon( ty.hull, border = "orange", lty = "dashed" )
	} else {
		plot( sample.ty.old, xlim = c(0,38), ylim= c(0,85) ); points( t(y.obs), pch = 19, col = "blue" ) 
		polygon( ty.hull.old, border = "orange", lty = "dashed" )
	}
	hullplot.pts <- sample.reduced.ty.new[ chull( sample.reduced.ty.new ), , drop=FALSE]
	polygon( hullplot.pts, border = "blue", lty = "dotted" )
	
	#### Plot thetas ##########################################
	plot( theta.hist, pch=20 )
	arrows( theta.old[1], theta.old[2], theta.current[1], theta.current[2], length = 0.05 )
	if( gdor.true[1] != 0 ) {
		abline( a = 0, b = gdor.true[2]/gdor.true[1], col = "orange", lty = "dotted" )
	}	
	##### y.obs in exterior? #########################################
	# lpcdd approach described in rcdd vignette Section 4.4 using strongly separating hyperplane
#	is.exterior <- in.exterior(y.obs, sample.reduced.ty.new )  # calls lpcdd
	hrep <- cbind( 0, 0, 1, qneg(ty.hull) )	# s.t. 1st ineq
	hrep <- rbind( hrep, c(0, 1, 1, -y.obs ) )	# 2nd (artifical) ineq
	f <- q2d( lpcdd( hrep, d2q( c(-1, y.obs) ), minimize = FALSE )$optimal.value )
	is.exterior <- ( f > 0 )

	if( is.exterior ) # y.obs in exterior
	{
		##### exterior case #########################################
		cat( k, ": y.obs is outside the convex hull of the sample.  Keep sampling!\n" )
		on.boundary <- FALSE
	} else{
		##### NOT exterior case: boundary or interior #########################################
		cat( k, ": y.obs is in interior or boundary.\n" )
		##### Find the empirical face at y.obs #########################################
		# apply linearity to CUMULATIVE HULL to find face, but only if ty.hull has more than 1 point
		if( nrow( ty.hull ) > 1 )
		{	# linearity breaks if there is only one point
			W <- sweep( ty.hull, 2, y.obs, FUN = "-" )
			W <- d2q( cbind( 0, 0, W ) )
			lin.active 		<- linearity( W, rep = "V" )
		} else lin.active <- 1 # if ty.hull has only 1 point, then that's the face!

		face.hull.ty 	<- ty.hull[ lin.active, , drop = FALSE ]
		face.extreme.count 		<- nrow( face.hull.ty )
		if( face.extreme.count == 0 )
		{	cat( k, ":??  How did we get here??  0 dim face is not possible when not exterior.\n" )
		} else if( (face.extreme.count == ty.hull.extreme.count) ) 		{
			##### Face is full dim, so y.obs in interior #########################################
			# full dimensional face.  
			cat( k, ": y.obs is in interior of the convex hull of samples.\n")
			on.boundary <- FALSE
			on.interior <- TRUE
		} else
		{
			##### y.obs on boundary, there's a new lower dim empirical face #######################
			on.boundary <- TRUE
			on.interior <- FALSE
			cat( k, ": y.obs is on the boundary of the convex hull of samples.\n" )
			cat( "Empirical face has dimension: ", length(lin.active)-1, "\n" )

			##### Find a GDOR to this face #######################
			objv <- d2q( c( rep(0,d), 1) )
			gdor.hrep <- cbind( qneg(W), 0 ) 
			gdor.hrep[  lin.active, 1 ] <- 1
			gdor.hrep[ -lin.active, ncol(gdor.hrep) ] <- "-1"
			gdor.hrep <- rbind( gdor.hrep, c( 0, 1, rep(0,d), -1 ) )
			lpout <- lpcdd( gdor.hrep, objv, minimize = FALSE )
			gdor  <- lpout$primal.solution[ -(d+1) ]
			cat( "GDOR", gdor, "\n" )
	
			##### From CURRENT sample, find pts on face #######################
			sample.ty.rays  <- sweep( sample.ty.new, 2, y.obs, FUN = "-" )
			# Find sample.ty.new restricted to face
			temp.face.ty 	<- sample.ty.new[ sample.ty.rays %*% q2d(gdor) == 0, , drop = FALSE ]
			# unique face pts and the count
			temp.face.ty.unique <- unique( temp.face.ty )
			temp.face.counts <- apply( temp.face.ty.unique, 1, function( pair ) count.pair( pair, temp.face.ty ) )
			temp.face.total <- sum( temp.face.counts )
			# Cheat and get exact probabilities to compare to.
			temp.face.theta.probs <- t.likelihood.tset( theta.current, temp.face.ty.unique, Y.all.mat, pt.flag = TRUE )
			temp.face.theta.probsum <- sum( temp.face.theta.probs )
			cat( "The empirical face points are: \n" )
			print( cbind( temp.face.ty.unique, temp.face.counts/temp.face.total,temp.face.theta.probs ) )
			cat( temp.face.theta.probsum, "\n" )
			samp.size.effective <- nrow( sample.ty.new )
			temp.face.prop  <- temp.face.total / samp.size.effective

			# stuff to cheat with for uniroot.  I need face pts and actual freq on these pts
			temp.face.g  	<- sort( apply( temp.face.ty.unique, 1, function( pair ) pair2g( pair, Y.all.mat ) ) )
			temp.face.freq 	<- Y.all.freq[ temp.face.g, , drop = FALSE]
			temp.face.mat 	<- as.matrix( temp.face.freq ) # to pass to theta2mu
		
			# the truth - how far are we off?
#			face.counts 	<- apply( face.ty, 1, function( pair ) count.pair( pair, sample.ty.new ) )
#			if( verbose.flag )
#			{ 	cat( "The true face points are: \n" )
#				print.table( cbind( face.ty, face.probs, face.counts ) )
#			}
			##### Enough pts on face? #########################################
			# Do we have more than X% on the face?
			if( temp.face.prop < face.cutoff ) # keep sampling
			{
				cat( sprintf( "Only %s point(s) out of %s occurred on face.  Keep sampling!\n\n", temp.face.total, samp.size.effective ) )
			} else {
				##### It's a real FACE.  Let's go #########################################
				cat( sprintf( "We've got %s points out of %s on this face.  MOVE ON TO LCM!\n\n", temp.face.total, samp.size ) )
				LCM.flag <- TRUE
				ty.hull.old				<- ty.hull
			 	ty.hull 				<- temp.face.ty.unique
			 	ty.hull.extreme.count  	<- nrow( ty.hull )

				num.face.pts <- nrow( temp.face.ty.unique )
				if( num.face.pts == 1) # single face pt
				{
				# This means the constancy space is ANY vector in R^2.  So, at this point we should
				# try a variety of different DOCs and see what happens.
				}
			} # else.  temp.face.prop >= face.cutoff
		}  # else.  on boundary
	} # else.  in interior or boundary
	
	#### ##################################################
	# Propose a new search direction
	if( MC.flag )
	{	
		mu.current		<- colMeans( sample.ty.new )
	} else {
		if( LCM.flag == FALSE )
		{	mu.current	<- theta2mu( theta.current, Y.all.mat )
		} else {
			mu.current	<- theta2mu( theta.current, temp.face.mat )
		}
	}
	if( verbose.flag ) cat( "mu.current: ", mu.current, "\n" )
	
	#### Update del.ll ##################################################
	del.ll 		<- y.obs - mu.current
	del.ll.len 	<- sqrt( del.ll %*% del.ll )

	if( LCM.flag == FALSE && on.boundary == TRUE && k %% 5 == 0 )
	{
		cat( "Use GDOR:\n" )
		p.current <- q2d( gdor )
	} else if( LCM.flag == FALSE && k>=9 && k%%3==0 )
	{
		cat( "Original ML: use regression direction:\n" )
		temp.data <- theta.hist[ (k-7):k, ]
		m.temp <- lm( temp.data[,2] ~ temp.data[,1] )
		p.current <- gdor.prop <- c(1, coef(m.temp)[2])
		# if this is *not* an ascent direction, use the closet allowable ascent direction
		if( t(p.current) %*% del.ll <= 0 )
			p.current <- -p.current
#			p.current <- closest.ascent.dir( p.current, del.ll )
	} else {
		cat( "Use steepest ascent direction:\n" )
		p.current <- del.ll
	}

	cat( "Next search direction: " )
	print.table( p.current )
	del.ll.tp.current <- t(del.ll) %*% p.old
	face.current.probs <- t.likelihood.tset( theta.current, face.ty, Y.all.mat, pt.flag = TRUE )
	if( LCM.flag ) cat( "Maximizing LCM ..." )

#	cat( sprintf("Step length: %11.6f del.ll.len: %11.6f\n", alpha.current, del.ll.len  ) )
#	if( verbose.flag )
#		cat( "Current theta face probs: ", face.current.probs, ":", sum( face.current.probs ), "\n\n\n" )
	
	if( LCM.flag ) LCM.k <- LCM.k+1

	k <- k+1
	cat( "--------------------------\n\n\n\n" )
	
}
cat( "Exited main loop!  del.ll.len: ", del.ll.len, "theta", theta.current, "\n" )

cbind( theta.current, theta.MLE, theta.LCM )
cat( gdor )
m1 <- lm( V2 ~ V1, data = as.data.frame(theta.hist[-1,]) )
#summary( m1 )
abline( m1 )
abline( a = coef(m1)[1], b = gdor.true[2]/gdor.true[1], col = "orange", lty = "dotted" )
##########################################################################
# What does the current distribution look like according to theta.current?
prob <- t.likelihood.tset( theta.current, t.y, Y.all.mat, pt.flag = TRUE )
prob.data <- data.frame( t.y, prob )
#qplot( edges, triangles, data=prob.data[ prob.data$prob > 0.001,], colour=prob ) + scale_colour_gradient(low="steelblue", high="yellow")
#print.table( as.matrix( prob.data[ prob > 0.0001, ] ) )
#if( is.matrix( face.ty ) )
#{	print.table( cbind( face.ty, face.probs,face.current.probs ) )
#} else cat( face.ty, face.probs,face.current.probs )
theta.current
(theta2mu( theta.current, Y.all.mat ) )  # using full space for support
(theta2mu( theta.current, face.freq.mat ) )  # using "true" support, based on our full knowledge






#########################################################################################
# One-sided CI
#########################################################################################

if( LCM.flag )
{
(theta2mu( theta.current, temp.face.mat ) )
#Y.all.freq[ pair2g( y.obs, Y.all.mat ), ]

temp.face.mat
temp.face.ty.unique
#delta <- colMeans( ncone.mat )
delta <- q2d( gdor )

#########################################################################################
prob.on.face <- function( s, theta.current, delta, temp.face.ty.unique, Y.all.mat, samp.size )
{
	theta.CI <- theta.current + s*delta
	sample.CI <- get.network.sample( samp.size, theta.CI, Y.all.mat )
	sample.CI.ty <- Y.all.mat[ sample.CI, c(1,2) ]

	face.ty.list <- apply( temp.face.ty.unique, 1, paste, collapse = "-" )
	sample.CI.list <- apply( sample.CI.ty, 1, paste, collapse = "-" )
	sample.table <- table( sample.CI.list )
	face.ty.sums <- sapply( face.ty.list, function( pair ) sum( sample.table[names(sample.table) == pair ] ) )
	face.ty.sums
	sum( face.ty.sums ) / samp.size
}

#########################################################################################
prob.on.face.minus <- function( s, face.prob, theta.current, delta, temp.face.ty.unique, Y.all.mat, samp.size )
{
	prob.on.face( s, theta.current, delta, temp.face.ty.unique, Y.all.mat, samp.size ) - face.prob
}

s.val <- uniroot( prob.on.face.minus, c(-20,20), 0.05, theta.current, delta, temp.face.ty.unique, Y.all.mat, samp.size )
s.val
prob.on.face( s.val$root, theta.current, delta, temp.face.ty.unique, Y.all.mat, samp.size )
( theta.endpt <- theta.current + s.val$root*delta )
cat( "CI for theta: ", theta.endpt )

# in mean value parameterization?
#sample.CI <- get.network.sample( samp.size, theta.endpt, Y.all.mat ) 
#sample.CI.mat 	<- Y.all.mat[ sample.CI, ]
#sample.CI.ty  		<- sample.CI.mat[ , c(1,2) ]
#colMeans( sample.CI.ty )

#theta2mu( theta.endpt, Y.all.mat )
#theta2mu( theta.endpt + 100*delta, Y.all.mat )
#theta2mu( theta.MLE, Y.all.mat )
}
