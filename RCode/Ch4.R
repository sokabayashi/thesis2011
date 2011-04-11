library( rcdd )
library( ggplot2 )
sample.gy <- t(matrix( c(
19, 11,
21, 14,
19, 11,
21, 13,
21, 13,
23, 21,
14, 3,
21, 16,
15, 4,
17, 8), ncol = 10) )

sample.gy.unique <- unique( sample.gy ) 
( sample.gy.Vrep 	 <- cbind( 0, 1, sample.gy.unique ) )
( sample.gy.reduced <- redundant( d2q(sample.gy.Vrep), rep = "V" ) )
( sample.gy.Hrep <- scdd( sample.gy.reduced$output, rep = "V" ) )

sample.gy2 <- t(matrix( c(
20, 14,
11, 1,
19, 10,
12, 3,
15, 5,
18, 7,
15, 6,
20, 13,
17, 9,
18, 5), ncol=10) ) 
sample.gy2.unique 	<- unique( sample.gy2 ) 
( sample.gy2.Vrep 	 <- cbind( 0, 1, sample.gy2.unique ) )
( sample.gy2.reduced <- redundant( d2q(sample.gy2.Vrep), rep = "V" ) )
( sample.gy2.Hrep <- scdd( sample.gy2.reduced$output, rep = "V" ) )


sample.gy.combined <- unique( rbind( sample.gy.reduced$output, sample.gy2.reduced$output ) )
( sample.gy.combined <- redundant( sample.gy.combined, rep = "V" ) )
( sample.gy.Hrep <- scdd( sample.gy.combined$output, rep = "V" ) )

#########################################
# Get the figures of the hulls
mydat <- data.frame( rbind( cbind(sample.gy,1), cbind(sample.gy2, 2) ) )
mydat.polygon <- data.frame( rbind( cbind(sample.gy[chull(sample.gy),],1), cbind(sample.gy2[chull(sample.gy2),], 2) ) )
combined.polygon <- data.frame( unique( rbind(sample.gy, sample.gy2 ) ) )
combined.polygon <- combined.polygon[chull( combined.polygon), ]
p <- ggplot( mydat, aes(x=X1,y=X2), aspect.ratio = 1 ) 
figureplot <- p + geom_polygon( data=combined.polygon,fill=0,colour="grey" )  + geom_polygon(data=mydat.polygon[mydat.polygon$X3==1,],fill=0,colour="blue", size = 1, linetype = 2 ) + geom_polygon(data=mydat.polygon[mydat.polygon$X3==2,],fill=0,colour="red",  size = 1, linetype = 4 ) + geom_point(aes(shape=X3) )+ theme_bw()+ coord_equal(ratio = 1)
#ggsave( "~/Tako/THESIS/Figures/combined-hull-norm.pdf" )


# very messy!
sample.gy.combined <- rbind( sample.gy.Hrep$output, sample.gy2.Hrep$output )
( sample.gy.combined <- redundant( sample.gy.combined, rep = "H" ) )
combined.gy.scdd <- q2d( scdd( sample.gy.combined$output, rep = "H" )$output )
combined.gy.scdd <- combined.gy.scdd[ ,-c(1:2)]
points( combined.gy.scdd, pch=19 )

mydat <- data.frame( rbind( cbind(sample.gy,"1"), cbind(sample.gy2, "2") ) )
mydat.unique <- data.frame( rbind( cbind(sample.gy.unique,1), cbind(sample.gy2.unique, 2) ) )
mydat.polygon <- data.frame( rbind( cbind(sample.gy[chull(sample.gy),],1), cbind(sample.gy2[chull(sample.gy2),], 2) ) )
combined.gy <- sample.gy.combined$output
combined.polygon <- data.frame( unique( rbind(sample.gy, sample.gy2 ) ) )

scdd( d2q(as.matrix( cbind(0,1,combined.polygon) ) ), rep = "V" )

p <- ggplot( mydat, aes(x=X1,y=X2) ) 
p + geom_point(aes(shape=X3))

p + geom_polygon(data=mydat[mydat$X3==1,] )
p <- ggplot( mydat, aes(x=X1,y=X2) ) 
p  + geom_polygon(data=mydat.polygon[mydat.polygon$X3==1,],fill=0,colour="blue", linetype = 2 ) + geom_polygon(data=mydat.polygon[mydat.polygon$X3==2,],fill=0,colour="red",  linetype = 2 ) + geom_polygon( data=combined.polygon,fill=0,colour="grey" ) + geom_point(aes(shape=X3),size=2.5, fill=1 )
#ggsave( "~/Tako/THESIS/Figures/combined-hull-norm.pdf" )

# the old way was so much easier ...
plot( rbind(sample.gy, sample.gy2) )
polygon( sample.gy[chull(sample.gy),] )
polygon( sample.gy2[chull(sample.gy2),] )
polygon( combined.polygon )
#########################################

# use lpcdd to test for exterior point
# exterior point
q <- c(20,10) 
( hull.pts <- sample.gy.combined$output[ , -c(1:2) ] ) # get points in hull
hrep <- cbind( 0, 0, 1, qneg(hull.pts) )	# s.t. 1st ineq
( hrep <- rbind( hrep, c(0, 1, 1, -q) ) )	# 2nd (artifical) ineq
lpcdd( hrep, d2q( c(-1, q) ), minimize = FALSE )$optimal.value

# interior point
q <- c(20,12) 
 hull.pts <- sample.gy.combined$output[ , -c(1:2) ]  # get points in hull
hrep <- cbind( 0, 0, 1, qneg(hull.pts) )	# s.t. 1st ineq
( hrep <- rbind( hrep, c(0, 1, 1, -q) ) )	# 2nd (artifical) ineq
lpcdd( hrep, d2q( c(-1, q) ), minimize = FALSE )$optimal.value

# boundary point
q <- c(21,13) 
 hull.pts <- sample.gy.combined$output[ , -c(1:2) ]  # get points in hull
hrep <- cbind( 0, 0, 1, qneg(hull.pts) )	# s.t. 1st ineq
( hrep <- rbind( hrep, c(0, 1, 1, -q) ) )	# 2nd (artifical) ineq
lpcdd( hrep, d2q( c(-1, q) ), minimize = FALSE )$optimal.value

# interior/exterior the boneheaded way
x <- rbind( c(20,10), c(20,12), c(21,13) )
l <- sample.gy.Hrep$output[ ,1]
b <- sample.gy.Hrep$output[ ,2]
a <- sample.gy.Hrep$output[ ,-c(1,2)]
a <- qneg( a )
axb <- qmatmult( a, t(x) )
axb <- sweep( axb, 1, b, FUN=qmq ) # subtract b
apply( axb, 2, function(pair) max(qsign(pair)) )
# pos = exterior, neg = interior, 0 = boundary
#########################################

# linearity
# interior point
q <- c(20,12)
( W <- sweep( hull.pts, 2, q, FUN = qmq ) )
W <- cbind( 0, 0, W )
( lin.active <- linearity( W, rep = "V" ) )
( L <- W[ lin.active, , drop=FALSE] )

# boundary point - 1 dim
q <- c(21,13)
( W <- sweep( hull.pts, 2, q, FUN = qmq ) )
W <- cbind( 0, 0, W )
( lin.active <- linearity( W, rep = "V" ) )
( L <- W[ lin.active, , drop=FALSE] )

# boundary point - 2 dim
# need to find a point
sample.gy.Hrep
b <- sample.gy.Hrep$output[,2]
A <- qneg( sample.gy.Hrep$output[,-c(1:2)] )
xcord <- as.matrix( rep("20",5) )
y <- qdq( qmq( b , qxq( A[,1], xcord ) ), A[,2] )
# got one
q <- c("20","31/3")
( W <- sweep( hull.pts, 2, q, FUN = qmq ) )
W <- cbind( 0, 0, W )
( lin.active <- linearity( W, rep = "V" ) )
( L <- W[ lin.active, , drop=FALSE] )
hull.pts[ lin.active, , drop=FALSE ]
sweep( L[, -c(1:2), drop = FALSE ], 2, q, FUN = qpq )

#########################################
# tangent cone at q <- c(21,13)
q <- as.matrix( c("21","13") )
sample.gy.Hrep
l <- sample.gy.Hrep$output[ ,1]
b <- sample.gy.Hrep$output[ ,2]
a <- sample.gy.Hrep$output[ ,-c(1,2)]
( active <- qpq( qmatmult( a, q ), b ) == "0" )
tcone.Hrep <- sample.gy.Hrep$output[ active, , drop=FALSE]
tcone.Hrep[ ,2 ] <- "0"
tcone.Hrep
# normal cone
( ncone.Vrep <- qneg( tcone.Hrep ) )
q.d <- q2d(q )
ncone.dir <- q2d( ncone.Vrep[ ,-c(1:2), drop = FALSE] )
arrowplot <- geom_segment( aes( x=q.d[1], y = q.d[2], xend=q.d[1] + ncone.dir[,1], yend = q.d[2] + ncone.dir[,2] ), colour="blue", arrow= arrow(length=unit(0.2,"cm")) ) 
####
q <- as.matrix( c("20","31/3") )
sample.gy.Hrep
l <- sample.gy.Hrep$output[ ,1]
b <- sample.gy.Hrep$output[ ,2]
a <- sample.gy.Hrep$output[ ,-c(1,2)]
( active <- qpq( qmatmult( a, q ), b ) == "0" )
tcone.Hrep <- sample.gy.Hrep$output[ active, , drop=FALSE]
tcone.Hrep[ ,2 ] <- "0"
tcone.Hrep
# normal cone
( ncone.Vrep <- qneg( tcone.Hrep ) )
#########
# add normals to figureplot
q.d2 <- q2d(q )
ncone.dir2 <- q2d( ncone.Vrep[ ,-c(1:2), drop = FALSE] )
arrowplot2 <- geom_segment( aes( x=q.d2[1], y = q.d2[2], xend=q.d2[1] + ncone.dir2[,1], yend = q.d2[2] + ncone.dir2[,2] ), colour="blue", arrow= arrow(length=unit(0.2,"cm")) ) 

figureplot + arrowplot + arrowplot2 + coord_equal(ratio = 1) + theme_bw(base_size = 14 )  + opts( legend.position = "none" ) 
#ggsave( "~/Tako/THESIS/Figures/combined-hull-norm.pdf" )



#########################################
# GDOR using lpcdd
q <-  c("20","31/3")
#q <-  c("21","13")
d <- length( q )
( hull.pts <- sample.gy.combined$output[ , -c(1:2) ] ) # get points in hull
( W <- sweep( hull.pts, 2, q, FUN = qmq ) )
W <- cbind( 0, 0, W )
( lin.active <- linearity( W, rep = "V" ) )
( objv <- d2q( c( rep(0,d), 1) ) )
hrep <- cbind( qneg(W), 0 ) 
hrep[ lin.active, 1 ] <- 1
hrep[ -lin.active, ncol(hrep) ] <- "-1"
hrep <- rbind( hrep, c( 0, 1, rep(0,d), -1 ) )
hrep
lpout <- lpcdd( hrep, objv, minimize = FALSE )
( gdor  <- lpout$primal.solution[ -(d+1) ] )

#########################################
gdor <- as.matrix(c( "1", "-3/8" ) )
q <- c( "20", "31/3" )
sample.gy.3 <- t(matrix( c(
"17", "7/3",
"18", "5",
"19", "23/3",
"20", "31/3",
"21", "13",
"22", "47/3",
"20", "10",
"20", "32/3",
"19", "8",
"19", "25/3"
), ncol = 10) )
W <- sweep( sample.gy.3, 2, q, FUN = qmq )
qmatmult( W, gdor )
