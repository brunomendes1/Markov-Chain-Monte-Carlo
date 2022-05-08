####################################################################################################
# This is an implementation of a Bayesian statistical method that uses the Reversible Jump Markov
# Chain Monte Carlo (RJMCMC) algorithm [Green P (1995). Reversible jump Markov chain Monte
# Carlo computation and Bayesian model determination. Biometrika, 82, 711–732.] to estimate the
# number and duration of regimes that were present in an tracer experiment performed at the 
# Waste Isolation Pilot Plant [Gonzales D, Bentley C (1984). Field test for effective
# porosity and dispersivity in fractured dolomite: the WIPP, Southeastern New Mexico. In Groundwater
# Hydraulics, Rosenshein JS, Bennett GD (editors), Washington DC: American Geophysical Union, 207–
# 221.]
# RJMCMC is bascially a Markov Chain Monte Carlo simulation where the number of parameters of the model
# is allowed to change during the simulation.
# 
# The main idea is to investigate if the experimental data curve can be explained by a single analytical
# solution (single regime) or if splitting the curve in sub-regimes (i.e. a different analytical solution
# per regime) improves the fit to the data.
# More details in poster-isba-08-v1-no-seal.pdf
# rjmcmc parameter vector theta
# theta = ( l.nu, l.D, l.sigma2 )
#   ( theta ) ~ flat
#   ( y_i | theta ) ~indep N( mu( time_i, exp( l.nu ), exp( l.D ) ), 
#     exp( l.sigma^2 ) ), i = 1, ..., n, mu( time, exp( l.nu ), 
#     exp( l.D ) ) = the complicated function involving A2( x, time )
#
# Inputs:
#
#   y = vector of log concentrations, of length n = sample size
#   time = vector of time points, also of length n
#   kappa = scaling factor for Metropolis proposal distribution
#     (affects the acceptance rate R; to increase R, decrease 
#     kappa) 
#   Sigma = proposal distribution covariance matrix
#   theta.0 = initial value for ( l.nu, l.D, l.sigma2 ) in 
#     Metropolis iterations
#   n.burnin = length of burn-in period
#   n.monitor = length of monitoring period
#   n.thin = thinning constant (only every n.thin-th iteration
#     in the monitoring period will be written to disk)
#   seed = random number seed (for generating repeatable 
#     sequences of Metropolis iterations); must be an integer 
#     from 0 to 1000
#   output.file.prefix = character string naming where you want 
#     the MCMC data set to go; for example, output.file.prefix 
#     = "wipp1" would write the MCMC data set to the file 
#     "wipp1.d"
# 
# Outputs:
# 
#   Acceptance rate R returned when iterations are finished
#   A file called paste( output.file.prefix, ".d", sep = "" ) is
#     written (in the same directory where R has been called) 
#     containing one row for each monitored iteration and six 
#     columns: the monitored iteration number (from 1 to 
#     n.monitor/n.thin), the simulated draws from the posterior 
#     for theta = ( l.nu, l.D, l.sigma2 ) for that iteration, and 
#     the corresponding simulated draws from the posterior for 
#     ( nu, D, sigma ). If the output file exists before the 
#     function is invoked, it will be over-written




# Loading experimental data

time.points <- scan("wipp-t.d")

y <- scan("wipp-y.d")

y <- exp(y)

y <- 0.001*y

y <- log(y)

plot(time.points,y, xlab='time',ylab='Log(Concentration)',main='WIPP experimental data')



 
PD.sim <- function( theta, p, L.kappa ) {
# Proposal distribution simulation  

  Z <- matrix( rnorm( p ), p, 1 )
   
  Mu <- matrix( theta, p, 1 )
  theta.star <- c( Mu + ( L.kappa %*% Z ) )

  return( theta.star )

}


alpha <- function( theta.old, theta.new, y ) {
# Acceptance probability calculation

  return( min( 1, exp( log.post( theta.new, y ) - log.post( 
    theta.old, y ) ) ) )

}

log.post <- function( theta, y ) {
# log( posterior ) calculation

  return( log.prior( theta ) + log.lik( theta, y,time.points ) )

}

log.prior <- function( theta ) {
# flat prior
  return( 0 )

}


log.lik <- function( theta, y, time.points ) {
# log( likelihood ) calculation

  l.nu <- theta[ 1 ]
  l.D <- theta[ 2 ]
  l.sigma2 <- theta[ 3 ]
  n <- length( y )

# Compensating the fact there is one too many contributions for the log-lik at the change-points          
  if(time.points[1]!=75){
  time.points <- time.points[-1]
  y <- y[-1]
  n <- n-1}

  ll <- - n * l.sigma2 / 2 - sum( ( y - lnC.new( time.points, l.nu, 
    l.D ) )^2 ) / ( 2 * exp( l.sigma2 ) ) 

  return( ll )

}

#######################################################################################
# Definition of the analytical solution of a advection-diffusion differential equation
#

#error function
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower=FALSE)

# There is an extra parameter for the analyical solution - retention-,
# This is a physical parameter that represents the level the tracer interacts with the rock.
# a value of 1 means there is not interaction between the tracer and the rock
# this keeps the problem bi-dimensional (velocity,dispersivity) and therefore easier to visualize. 
# Previous exploratory work also showed this physical parameter is not very important
# it's something we could include on future iterarions of this code
 
Ret <- 1.0

#x is distance, in meters, from the source where the tracer is measured
x <- 22.858


#t stands for time, in days, l.nu for log-velocity of the water, l.D stands for log-dispersivity
A2.new <- function( t, l.nu, l.D ) {

a1 <- exp( exp( l.nu ) * x / exp( l.D ) )
a2 <- erfc( ( Ret * x +  exp( l.nu )* t ) / ( 2 * sqrt( exp( l.D ) * Ret * t ) ) )
i <- a1*a2
i[is.infinite(i)] <- 0

if(is.infinite(a1)){ i <- 0}

a <- erfc( ( Ret * x - exp( l.nu ) * t ) / ( 2 * sqrt( exp( l.D ) * Ret * 
    t ) ) ) / 2 + sqrt( exp( l.nu )^2 * t / ( pi * exp( l.D ) * Ret ) ) * 
    exp( - ( Ret * x - exp( l.nu ) * t )^2 / ( 4 * exp( l.D ) * Ret * t ) ) - 
    ( 1 + exp( l.nu ) * x / exp( l.D ) + exp( l.nu )^2 * t / ( exp( l.D ) * 
    Ret ) ) * i / 2
    
return( a )
}

# Time at which we start measuring the amount of tracer 22.8 meters downstream from the source
t0 <- 28

# C0 is the amount of tracer introduced in the systerm at x=0
C0 <- 0.721

lnC.new <- function( t, l.nu, l.D ) {
# calculation of the tracer profile

  diff <- A2.new( t, l.nu, l.D ) - A2.new( t - t0, l.nu, l.D )

  diff[ diff<=0 ] <- 10^(-32)
  
  diff[ is.nan(diff) ] <- 10^(-32)

  return( log( C0 ) + log( diff ) )
}

######################################
#exploring the log-likelihood function
# it's a 2-d surface

nr.points <- 100

matrix.ll <- array(dim=c(nr.points,nr.points))

lv <- seq(log(0.001),log(10),length.out=nr.points )
ld <- seq(log(0.001),log(50),length.out=nr.points )

# Producing matrix with all mle coefficients
for ( k in 1: nr.points ) {
 for ( i in 1:nr.points ) {
   matrix.ll[k,i] <- -(- n * log( sum( ( y - lnC.new(time.points,lv[k],ld[i]) )^2 ) / n )- n / 2)
 }
 print(k)
}

alarm()

contour(lv,ld,matrix.ll,xlab='Log(V)',ylab='LogD', nlevels=18,
main='Log-likelihood \n Fitting 1-regime model to \n artificial data with 2 change-points')

min(matrix.ll,na.rm=T)
[1] -19.78463

matrix.ll[54,62]
[1] -19.78463

plot(time.points,y, xlab='Time',ylab='Log(C)')
lines(time.points,lnC.new(time.points,lv[54],ld[62]))

 c(lv[54],ld[62])
[1] -1.9769670 -0.2410232
> 
>  exp(c(lv[54],ld[62]))
[1] 0.1384886 0.7858234


#######################################################################
# Defining initial values for MCMC and covariance matrix


ssq <- function (z){
#sum of squares function

l.nu <- z[1]
l.D <- z[2]

return(sum( ( y - lnC.new(  time.points, l.nu, l.D ) )^2 ))
}

lsigma <- ssq(c( lv[54],ld[62]) )/n

 0.5170822
 
loglik.new <- function( theta ) {
l.nu <- theta[1]
l.D <- theta[2]
sigma1 <- theta[3]

return(-(- n * sigma1/2 - sum( ( y - lnC.new(  time.points, l.nu, l.D ) )^2 ) / ( 2 * exp( sigma1 ))))
}


#####################################################
# Using the hessian of the log-lik function in order
# to optimize the proposal distribution
# this makes is more likely the mcmc will stay around the maximum of
# the log-lik function during the simulation


nlm (loglik.new,c(lv[54],ld[62],lsigma),hessian=TRUE)

$minimum
[1] 21.09038

$estimate
[1] -1.9841378 -0.2284500 -0.6598327

$gradient
[1]  6.815598e-06  3.990124e-06 -2.785292e-06

$hessian
             [,1]         [,2]        [,3]
[1,] 745.89914831 405.66688071 -0.03729923
[2,] 405.66688071 438.06558594 -0.02191243
[3,]  -0.03729923  -0.02191243 61.99380351

$code
[1] 1

$iterations
[1] 9

 
exp(c(-1.9841378, -0.2284500, -0.6598327))
[1] 0.1374991 0.7957661 0.5169378

nlm.art <- nlm (loglik.new,c(lv[54],ld[62],lsigma),hessian=TRUE)

solve(nlm.art$hessian)

              [,1]          [,2]         [,3]
[1,]  2.700989e-03 -2.501227e-03 7.409907e-07
[2,] -2.501227e-03  4.599003e-03 1.206810e-07
[3,]  7.409907e-07  1.206810e-07 1.613064e-02

# Load numDeriv first

hessian.art <- hessian(loglik.new, c(lv[54],ld[62],lsigma), method.args = list(d = 0.6),method="Richardson")

solve(hessian.art)
            [,1]          [,2]          [,3]
[1,]  0.0087136486 -0.0080637878  0.0003725726
[2,] -0.0080637878  0.0147222301 -0.0006523497
[3,]  0.0003725726 -0.0006523497  0.0523425693


# Load package nlme first

fdhess.art <- fdHess(c(lv[54],ld[62],lsigma), loglik.new)

solve(fdhess.art$Hessian)

             [,1]          [,2]          [,3]
[1,]  0.008711961 -0.0080613350  0.0003741910
[2,] -0.008061335  0.0147185971 -0.0006532587
[3,]  0.000374191 -0.0006532587  0.0523385517
> 

# 1 regime only
Sigma <- matrix( c( 0.008711961, -0.0080613350, 0,
                   -0.008061335, 0.0147185971, 0, 
                   0, 0, 0.0523385517 ), 3, 3 )
                   
theta.0 <- c(-1.9841378, -0.2284500, -0.6598327)


#####################################################
# Main rjmcm function
# In each iteration of the simulation one of four moves is chosen randomly
# there are four move types: new analytical curves, change position of current changepoints
# birth of new changepoint, deletion of a current changepoint
# Include the first and last time points in the list of changepoints
# kappa.simple,
# kappa.birth1,kappa.birth2,kappa.birth3,output.file.prefix,p.stay,p.change,
# p.birth,p.delete


main.rjmcmc <- function( n.monitor,n.burnin,n.thin,theta.0, Sigma, kappa.simple,
kappa.birth1,kappa.birth2,kappa.birth3,output.file.prefix,p.stay,p.change,
p.birth,p.delete,seed  ) {
    
  n <- 124
  
  list.current.chgpts <- c(1,n)

  list.open.pts <- 2:(n-1)
   
  set.seed( seed )

  L <- t( chol( Sigma ) )

  L.kappa <- sqrt( kappa.simple ) * L

  theta = theta.0
  
  ll <- log.lik( theta, y, time.points )
     
  sd1 <- kappa.birth1
  sd2 <- kappa.birth2
  sd3 <- kappa.birth3
     
  R <- 0
  N.simple = 0
  N.birth = 0
  N.delete = 0
  N.move = 0
  R.simple = 0
  R.birth = 0
  R.delete = 0
  R.move = 0
  A.simple = 0
  A.birth = 0
  A.delete = 0
  A.move = 0
   
  for ( i in 1:( n.burnin + n.monitor ) ) {
    u <- runif( 1 )
    max <- length(list.current.chgpts)
    if( u<=p.stay )
    {
# No change in number of change.points or location of current ones
#      print('simple')
      N.simple = N.simple+1
      z <- simple.metropolis(theta,ll,list.current.chgpts,list.open.pts,A.simple,L.kappa,p.stay)
      theta <- z$theta
      ll <- z$ll
      list.current.chgpts <- z$list.current.chgpts
      list.open.pts <- z$list.open.pts
      A.simple <- z$A.simple
      R.simple <- A.simple/N.simple
    }else
    {
      if( u > p.stay && u <= (p.change+p.stay) && max>2 ) 
      {
# Move the location of a current change point
#        print('move')
        N.move = N.move + 1
        z <- change.position(theta,ll,list.current.chgpts,list.open.pts,A.move,p.move)
        theta <- z$theta
        ll <- z$ll
        list.current.chgpts <- z$list.current.chgpts
        list.open.pts <- z$list.open.pts
        A.move <- z$A.move
        R.move <- A.move/N.move
      }else
      {
        if( u>(p.change+p.stay) && 
        u<=(p.change+p.stay+p.birth) && max<n-2)  
        {
# introduce a new change point
#            print('birth')
            N.birth = N.birth + 1
            z <- birth.changepoint(theta,ll,list.current.chgpts,
                   list.open.pts,A.birth,sd1,sd2,sd3,p.birth,p.delete)
            theta <- z$theta
            ll <- z$ll
            list.current.chgpts <- z$list.current.chgpts
            list.open.pts <- z$list.open.pts
            A.birth <- z$A.birth
            R.birth <- A.birth/N.birth
        }else
          if( u>(p.change+p.stay+p.birth)&&max>2) 
          { 
# Remove a current change point
#            print('delete')
            N.delete = N.delete + 1
            z<- delete.breakpoint(theta,ll,list.current.chgpts,
                  list.open.pts,A.delete,sd1,sd2,sd3,p.delete,p.birth) 
            theta <- z$theta
            ll <- z$ll
            list.current.chgpts <- z$list.current.chgpts
            list.open.pts <- z$list.open.pts
            A.delete <- z$A.delete
            R.delete <- A.delete/N.delete
          }
      }
    }
    
# Taking point '1' and '124' out of the list of change points
    max <- length(list.current.chgpts)
    current.chgpts <- list.current.chgpts[-max]
    current.chgpts <- current.chgpts[-1]

    if ( ( i > n.burnin ) && ( ( i - n.burnin ) %% n.thin == 0 ) )
      write( c(( i - n.burnin ) / n.thin, max-2, theta,
      current.chgpts,R.simple,R.birth,R.delete,R.move ),
      paste( output.file.prefix, ".d", sep = "" ),
      ncolumns = length(theta)+ 8 + max,
      append = ( i > n.burnin + n.thin ) )      
    
   }
   
   alarm()
   return( theta,list.current.chgpts,ll )
}





#############################################
# Maintain nr breakpoints and their positions
#############################################

simple.metropolis <- function(theta,ll,list.current.chgpts,list.open.pts,A.simple,L.kappa,p.stay)
{

 #backup current theta and log-lik
 theta.old <- theta
 ll.current <- ll
 
 #backup current lists
 old.list.chgpts <- list.current.chgpts 
 old.list.openpts <- list.open.pts
 
 max <- length(list.current.chgpts)
 # Choose a random branch
 index <- ceiling(runif(1)*(max-1))

 # range of components of vector theta    
 sub.theta <- (1+(index-1)*3):((1+(index-1)*3)+2)
 branch <- list.current.chgpts[index]:list.current.chgpts[index+1]
 p <- length( theta.old[sub.theta] )
 
 # Take out the value of log-lik associated with this branch from total log-lik

 ll <- ll - log.lik( theta[sub.theta], y[branch], time.points[branch] )

 # sample new values for the choosen branch
 theta[sub.theta] <- PD.sim( theta.old[sub.theta], p, L.kappa ) 

 # Calculate likelihood with the new values for parameters in choosen branch

 ll <- ll + log.lik( theta[sub.theta], y[branch], time.points[branch] )
 
 #alpha         
 if (ll==-Inf && ll.current==-Inf) {alpha <- 0
 }else{alpha <- min( 1, exp( ll - ll.current ))}
 
 u <- runif( 1 )
 
 if( u > alpha )
 {         
  # If new position is NOT accepted, restore old list of chgpts
  # Restore old value of theta
  list.current.chgpts <- old.list.chgpts
  list.open.pts <- old.list.openpts
  ll <- ll.current
  theta <-theta.old
 }else{
   A.simple = A.simple + 1
 }
   return(theta,ll,list.current.chgpts,list.open.pts,A.simple)

}



######################
# Delete a breakpoint
######################

delete.breakpoint <- function(theta,ll,list.current.chgpts,list.open.pts,A.delete,sd1,sd2,
sd3,p.delete,p.birth)
{
  ll.current <- ll
  theta.old <- theta
           
  # backup current lists
  old.list.chgpts <- list.current.chgpts 
  old.list.openpts <- list.open.pts
  
  max <- length(list.current.chgpts)
  
  # choose randomly one of the current breakpoints to remove
  # Choose randomly a position in the vector that contains the positions 
  # of current chgpts 
  # This variable contains the POSITION in the list of break points
  # that is to be deleted          
  index.t.b.deleted <- ceiling(runif(1)*(max-2) + 1)
  point.t.b.deleted <- list.current.chgpts[index.t.b.deleted]
  
  # Take the amount of log-likelihood to the right and left of the point to be deleted
  # left
  index <- index.t.b.deleted
  parameter1.1 = (3*(index-1)-2)
  parameter2.1 = (3*(index-1)-1)
  parameter3.1 = (3*(index-1))
  
  sub.theta1 <- parameter1.1:parameter3.1
  branch1 <- list.current.chgpts[(index-1)]:list.current.chgpts[index]
  
  ll <- ll - log.lik( theta[sub.theta1], y[branch1], time.points[branch1] )

  #right
  parameter1.2= (3*index-2) # parameter 1 in branch 2
  parameter2.2= (3*index-1) # parameter 2 in branch 2
  parameter3.2= (3*index) # parameter 3 in branch 2
  sub.theta2 <- parameter1.2:parameter3.2 
  branch2 <- list.current.chgpts[index]: list.current.chgpts[index+1]
  
  ll <- ll - log.lik( theta[sub.theta2], y[branch2], time.points[branch2] )

  # New values for parameter vector. Average values of parameters to the LEFT and 
  # RIGHT of the the chosen one
  mean1 = (theta.old[ parameter1.1 ] + theta.old[ parameter1.2 ])/2
  mean2 = (theta.old[ parameter2.1 ] + theta.old[ parameter2.2 ])/2
  mean3 = (theta.old[ parameter3.1 ] + theta.old[ parameter3.2 ])/2

  theta.proposed <-theta.old
  theta.proposed[sub.theta1] = c( mean1, mean2, mean3 ) 

  theta.proposed = theta.proposed[-(sub.theta2)]

  branch.new <-  list.current.chgpts[(index-1)]:list.current.chgpts[(index+1)]
   
  # Part of the log-likelihood to be added to total log-likelihood
  ll <- ll + log.lik( theta[sub.theta1], y[branch.new], time.points[branch.new] )

  # remove the deleted point from current lists of breakpoints
  list.current.chgpts <- list.current.chgpts[ -index.t.b.deleted ]
               
  # Update lists of open points
  count <-1
  max <- length(list.open.pts)
  if(max==0){ list.open.pts <- point.t.b.deleted
  }else{ 
    while(list.open.pts[count]<point.t.b.deleted&&count<max)
    { 
    count=count+1
    }
    list.open.pts <- append(list.open.pts,point.t.b.deleted,after=count-1)}
     
  #alpha
  u1= (theta.old[ parameter1.1 ] - theta.old[ parameter1.2 ])/2
  u2= (theta.old[ parameter2.1 ] - theta.old[ parameter2.2 ])/2
  u3= (theta.old[ parameter3.1 ] - theta.old[ parameter3.2 ])/2
  
  #  if there are much more chgpts than open points the proposal ratio below will increase
  #  making alpha bigger and the delete move more likely to be accepted
  
  ratio.loglik = exp( ll - ll.current )
  
  proposal.undelete = p.birth*( dnorm(u1,mean=mean1,sd=sd1)*dnorm(u2,mean=mean2,sd=sd2)*
            dnorm(u3,mean=mean3,sd=sd3) ) * 1/ (length(old.list.openpts)+1)
  proposal.delete = p.delete / (length(old.list.chgpts)-2)
  
  if (ll==-Inf && ll.current ==-Inf) { alpha <-0
  }else{ if(ratio.loglik == Inf){ alpha <- 1
         }else{
            alpha <- min( 1, 
            ratio.loglik*proposal.undelete/proposal.delete/8 ) }
  }
   
  u <- runif( 1 )
                             
  # If new position is NOT accepted, restore old list of chgpts
  if( u > alpha )
  {  
   list.current.chgpts <- old.list.chgpts
   list.open.pts <- old.list.openpts
   ll <- ll.current
   theta <-theta.old
  }else{
   theta <-theta.proposed
   A.delete = A.delete + 1
 }
  
  return(theta,ll,list.current.chgpts,list.open.pts,A.delete)
}




##################
# Add a breakpoint
###################

birth.changepoint <- function(theta,ll,list.current.chgpts,list.open.pts,A.birth,sd1,sd2,
sd3,p.birth,p.delete)
{
 theta.old <-theta
 ll.current <- ll
 
 #backup current lists
 old.list.chgpts <- list.current.chgpts 
 old.list.openpts <- list.open.pts
         
 # choose new location for breakpoint randomly
 new.index <- ceiling(runif(1)*(length(list.open.pts)))
 new.point <- list.open.pts[new.index]
  
 #update lists of breakpoints
 count <-1
 max <- length(list.current.chgpts)
 while(list.current.chgpts[count]<new.point&&count<max)
 {count=count+1
 }

 parameter1 = (3*(count-2)+1)
 parameter2 = (3*(count-2)+2)
 parameter3 = (3*(count-2)+3)

 sub.theta <- parameter1:parameter3

 branch <- list.current.chgpts[(count-1)]:list.current.chgpts[count]
 
 # Extract the amount of log-likelihood associated with the branch that is about to be deleted
 ll <- ll - log.lik( theta[sub.theta], y[branch], time.points[branch] )
          
 # New values for parameter vector.
 vu <- rnorm(3,mean=c( theta.old[parameter1],theta.old[parameter2],theta.old[parameter3] ), 
 sd=c(sd1,sd2,sd3))
 
 u1 <- vu[1]
 u2 <- vu[2]
 u3 <- vu[3]
 
 theta[parameter1:parameter3] = 
   c(theta.old[parameter1]+u1,theta.old[parameter2]+u2,
    theta.old[parameter3]+u3) 
 theta = append(theta, c(theta.old[parameter1]-u1,
   theta.old[parameter2]-u2,theta.old[parameter3]-u3),after=3*(count-1))

 # Updating lists of current and open points
 list.current.chgpts <- append(list.current.chgpts,new.point,after=count-1)
 
 list.open.pts <- list.open.pts[-new.index]
  
 # Updating log-likelihood
 #left 
 sub.theta.left <- parameter1:parameter3
 branch.left <- list.current.chgpts[(count-1)]:list.current.chgpts[count]

 ll <- ll + log.lik( theta[sub.theta.left], y[branch.left], time.points[branch.left] )

 #right 
 index <- count + 1
 sub.theta.right <- (3*(index-2)+1):(3*(index-2)+3)
 branch.right <-  list.current.chgpts[(index-1)]:list.current.chgpts[index] 
 
 ll <- ll + log.lik( theta[sub.theta.right], y[branch.right], time.points[branch.right] )
       
 #alpha   
 # if there are many more open points than change-points the proposal ratio below will be high
 # which will make alpha higher and the likelihood of the birth move to be more easily accepted         
 # alpha is incomplete it needs the prob of the inverse transformation (deleting the new chgpt)
 p1 = dnorm(u1,mean=theta.old[parameter1],sd=sd1)
 p2 = dnorm(u2,mean=theta.old[parameter2],sd=sd2)
 p3 = dnorm(u3,mean=theta.old[parameter3],sd=sd3)
  
 ratio.loglik = exp( ll - ll.current )

 proposal.unbirth = p.delete / (length(old.list.chgpts)-1)
 proposal.birth = p.birth*p1*p2*p3 / length(old.list.openpts)

 if (ll==-Inf && ll.current ==-Inf) { alpha <-0
 }else{ if(ratio.loglik == Inf || p1==0 || p2==0 || p3==0){ alpha <- 1
         }else{ alpha <- min( 1, 
            ratio.loglik*proposal.unbirth*8/proposal.birth )
 }
 } 
 u <- runif( 1 )

 # If new position is NOT accepted, restore old list of chgpts
 if( u > alpha )
 {  
   list.current.chgpts <- old.list.chgpts
   list.open.pts <- old.list.openpts
   ll <- ll.current
   theta <-theta.old
 }else{
   A.birth = A.birth +1
 }

  return(theta,ll,list.current.chgpts,list.open.pts,A.birth)  
}



###########################################
# Maintain nr chgpts, change one position
###########################################

change.position <- function(theta,ll,list.current.chgpts,list.open.pts,A.move,p.move)
{

  ll.current <- ll
  theta.tmp <- theta
        
  #backup current lists
  old.list.chgpts <- list.current.chgpts 
  old.list.openpts <- list.open.pts
           
  max <- length(list.current.chgpts)
  # Choose randomly a position in the vector that contains the positions 
  # of current chgpts 
  # This variable contains the POSITION in the list of break points
  # that is to be changed          
  index.t.b.changed <- ceiling(runif(1)*(max-2) + 1)
  point.t.b.changed <- list.current.chgpts[index.t.b.changed]

  # Remove ammount of log-likelihood associated with this point
  # left
  index <- index.t.b.changed
  ll <- ll - log.lik( theta[(3*(index-2)+1):(3*(index-2)+3)], 
  y[ list.current.chgpts[(index-1)]:list.current.chgpts[index] ],
  time.points[ list.current.chgpts[(index-1)]:list.current.chgpts[index] ] )
  # right
  index <- index + 1
  ll <- ll - log.lik( theta[(3*(index-2)+1):(3*(index-2)+3)], 
  y[ list.current.chgpts[(index-1)]:list.current.chgpts[index] ], 
  time.points[ list.current.chgpts[(index-1)]:list.current.chgpts[index] ] )

  #remove that point from current lists of breakpoints
  list.current.chgpts <- list.current.chgpts[-index.t.b.changed]
                  
  #Choose new position for that point
  new.index <- ceiling(runif(1)*(length(list.open.pts)))
  new.point <- list.open.pts[new.index]
  
  #update lists of breakpoints
  max <- length(list.current.chgpts)

  count <-1
  while(list.current.chgpts[count]<new.point&&count<max)
  {count=count+1}
 
  list.current.chgpts <- append(list.current.chgpts,new.point,after=count-1)
 
  list.open.pts <- list.open.pts[-new.index]
         
  count2 <-1
  max <- length(list.open.pts)
  if( max==0){ list.open.pts <- point.t.b.changed
  }else{ while(list.open.pts[count2]<point.t.b.changed&&count2<max)
  {count2=count2+1}
  list.open.pts <- append(list.open.pts,point.t.b.changed,after=count2-1)}

# Add log-likelihood associated with new point
# left
  ll <- ll + log.lik( theta[(3*(count-2)+1):(3*(count-2)+3)], 
  y[ list.current.chgpts[(count-1)]:list.current.chgpts[count] ], 
  time.points[ list.current.chgpts[(count-1)]:list.current.chgpts[count] ] )
# right
  count.plus <- count + 1
  ll <- ll + log.lik( theta[(3*(count.plus-2)+1):(3*(count.plus-2)+3)], 
  y[ list.current.chgpts[(count.plus-1)]:list.current.chgpts[count.plus] ], 
  time.points[ list.current.chgpts[(count.plus-1)]:list.current.chgpts[count.plus] ] )
   
  #alpha         
  if (ll!=-Inf && ll.current !=-Inf) {alpha <- min( 1, exp( ll - ll.current ) ) 
  }else{alpha <-0}
  
  u <- runif( 1 )
       
  if( u > alpha )
  {         
   # If new position is NOT accepted, restore old list of chgpts
   list.current.chgpts <- old.list.chgpts
   list.open.pts <- old.list.openpts
   ll <- ll.current
   theta <-theta.tmp
  }else
  {
      A.move = A.move + 1
  }
  
  return(theta,ll,list.current.chgpts,list.open.pts,A.move)
}

############################################################
#different simulations with slightly different parameters

main.rjmcmc(100000,10000,1,theta.0,Sigma,0.1,0.01,0.01,0.1,
"art.2chgpts.0.1-0.01-0.01-0.1.100000.0.5-0.25-0.125-0.125",
0.5,0.25,0.125,0.125,32612)

main.rjmcmc(100000,10000,1,theta.0,Sigma,0.1,0.01,0.01,0.01,
"art.2chgpts.0.1-0.01-0.01-0.01.100000.0.5-0.25-0.125-0.125",
0.5,0.25,0.125,0.125,32612)

main.rjmcmc(100000,10000,1,theta.0,Sigma,0.1,0.01,0.01,0.001,
"art.2chgpts.0.1-0.01-0.01-0.001.100000.0.5-0.25-0.125-0.125",
0.5,0.25,0.125,0.125,32612)
  
main.rjmcmc(100000,10000,1,theta.0,Sigma,0.1,0.01,0.01,0.0001,
"art.2chgpts.0.1-0.01-0.01-0.0001.100000.0.5-0.25-0.125-0.125",
0.5,0.25,0.125,0.125,32612)

main.rjmcmc(100000,10000,1,theta.0,Sigma,0.1,0.1,0.1,0.1,
"art.2chgpts.0.1-0.1-0.1-0.1.100000.0.5-0.25-0.125-0.125",
0.5,0.25,0.125,0.125,32612)

main.rjmcmc(100000,10000,1,theta.0,Sigma,0.1,0.1,0.1,0.01,
"art.2chgpts.0.1-0.1-0.1-0.01.100000.0.5-0.25-0.125-0.125",
0.5,0.25,0.125,0.125,32612)

main.rjmcmc(100000,10000,1,theta.0,Sigma,0.1,0.1,0.1,0.001,
"art.2chgpts.0.1-0.1-0.1-0.001.100000.0.5-0.25-0.125-0.125",
0.5,0.25,0.125,0.125,32612)

main.rjmcmc(100000,10000,1,theta.0,Sigma,0.01,0.1,0.1,0.1,
"art.2chgpts.0.01-0.01-0.01-0.001.100000.0.5-0.25-0.125-0.125",
0.5,0.25,0.125,0.125,32612)

main.rjmcmc(100000,10000,1,theta.0,Sigma,0.01,0.1,0.1,0.01,
"art.2chgpts.0.5-0.01-0.01-0.0001.100000.0.5-0.25-0.125-0.125",
0.5,0.25,0.125,0.125,32612)

main.rjmcmc(100000,10000,1,theta.0,Sigma,0.01,0.1,0.1,0.001,
"art.2chgpts.0.01-0.1-0.1-0.001.100000.0.5-0.25-0.125-0.125",
0.5,0.25,0.125,0.125,32612)

main.rjmcmc(100000,10000,1,theta.0,Sigma,0.01,0.01,0.01,0.1,
"art.2chgpts.0.01-0.01-0.01-0.1.100000.0.5-0.25-0.125-0.125",
0.5,0.25,0.125,0.125,32612)



main.rjmcmc(100000,10000,1,theta.0,Sigma,0.1,0.1,0.1,0.1,
"art.2chgpts.0.01-0.01-0.01-0.1.100000.0.5-0.25-0.125-0.125",
0.5,0.25,0.125,0.125,32612)

###
# Debugging
###



######
main.rjmcmc(1000000,1,1,theta.0,Sigma,0.8,1,1,1,
"1million-real.data-0.8-1-1-1-0.5-0.25-0.125-0.125",0.50,0.25,0.125,0.125,32612)

main.rjmcmc(1000000,1,1,theta.0,Sigma,1.5,0.9,0.9,0.9,
"1million-real.data-1.5-0.9-0.9-0.9-0.5-0.25-0.125-0.125",0.50,0.25,0.125,0.125,32612)


main.rjmcmc(1000000,1,1,theta.0,Sigma,0.8,1,1,1,
"1million-real.data-0.8-1-1-1-0.6-0.2-0.1-0.1",0.60,0.20,0.1,0.1,32612)


##################################
# post-processing of the output of the simulation
##################################

memory.limit(size = 3000)

chdata <- scan("1million-real.data-0.8-1-1-1-0.6-0.2-0.1-0.1.d",sep="\n",what="",nlines=1000000)

run.nr <- NULL
nr.chg.pts <- NULL
params.list <- list()
pos.chg.pts.list <- list()
acceptance <- NULL
pos.chg.pts.vector <- NULL
index <-1
count <- 1

for(i in seq(from=10000, to=length(chdata),by=100)) 
{
  tmp <- strsplit(chdata[i]," ")
  run.nr[count] <- as.numeric(tmp[[1]][1])
  nr.chg.pts[count] <- as.numeric(tmp[[1]][2])
  # Doesn't work if nr. change-points =0
  params.list[[count]] <- as.numeric(tmp[[1]][3:( 2+(nr.chg.pts[count]+1)*3 )])
  if(nr.chg.pts[count]!=0) {pos.chg.pts.list[[count]] <- as.integer(tmp[[1]][(3+(nr.chg.pts[count]+1)*3):(3+(nr.chg.pts[count]+1)*3+nr.chg.pts[count]-1)]) }
  else{ pos.chg.pts.list[[count]] <- NA}
  pos.chg.pts.vector[index:(index+nr.chg.pts[count]-1)] <- pos.chg.pts.list[[count]]
  index <- index+nr.chg.pts[count]
  acceptance[count] <- as.numeric( tmp[[1]][(4+(nr.chg.pts[count]+1)*3+nr.chg.pts[count]-1):(4+(nr.chg.pts[count]+1)*3+nr.chg.pts[count])] )
  if(i%%1000==0) print(count)
  count <- count + 1
}

alarm()

hist(pos.chg.pts.vector,breaks=1.5:123.5,
 main = paste ("WIPP data with 2 change points \n Frequency histogram of position of change points"), xlab="Position of change point i")

  [1]      0      0      0      0      0      0      0      0      0      4      0
 [12]      0      0      0     24   1052   8418  24114  43970  82928   8930   1336
 [23]     86      0      0     10      0      0      0      0      0      0      0
 [34]      0      0      0      0      0     26      0      0      0      0      0
 [45]      0      0     10      0      0      0      0      0      0      0      0
 [56]      0      0      0      0      0      0      0      0      0      0      0
 [67]      0      0      0      0      0      0      0      0      0      0      0
 [78]   2278    500     74      0      0      0      8      0      0   3669    254
 [89] 124380  40350  12948   9201   5529   2256    954    668    447    172    229
[100]    240     60     26     72    142    164    194     66     82    168    126
[111]     86   1178    576    662    550    204    456    848    636    268    586
[122]    242
 

hist(nr.chg.pts,breaks=0.5:5.5,
 main = paste ("WIPP data with 2 change points \n Frequency histogram of number of change points"),
 xlab="Number of change points")

ts.plot(nr.chg.pts,main = paste("Time series of number of change points"),ylab="Number of change points",xlab="MCMC step")

##
# Plotting some randomly sampled curves that were produced by a simulation 
# against the original curve
##
#
n.monitor <- 100000

tmp.pos.chgpt <- pos.chg.pts.list

for( i in seq(from=1, to=n.monitor,by=1000) )
{
 if(nr.chg.pts[i]!=0)
 {
   count <-0
   total.params <- params.list[[i]]
   for ( j in seq(from=i, to=1,by=-1) )
   {
     if( identical(pos.chg.pts.list[[j]],pos.chg.pts.list[[i]]) ){
       total.params <- total.params + params.list[[j]]
       count <- count+1
     }
   }
   average.params <- total.params/count
   
   curve <- list()  
   time1<-1
   tmp.pos.chgpt[[i]][nr.chg.pts[i]+1] <- 124
   plot(time,y,col='red')
   for( k in 1:(nr.chg.pts[i]+1) )
   {
     time2 <- tmp.pos.chgpt[[i]][k] 
     curve[[k]] <- lnC.new(time[time1:time2],average.params[2*k-1],average.params[(2*k)])
     points(time1:time2,curve[[k]])
     time1 <- pos.chg.pts.list[[i]][k]+1
   }
   par(ask=T)
 }
}


##############
# Inference on the model's parameters, conditional on two change points
# at 21 and 90
##############

n.monitor <-100000
tmp.pos.chgpt <- pos.chg.pts.list
vector.param1 <- NULL
vector.param2 <- NULL
vector.param3 <- NULL
vector.param4 <- NULL
vector.param5 <- NULL
vector.param6 <- NULL
vector.param7 <- NULL
vector.param8 <- NULL
vector.param9 <- NULL

count <-1
for( i in 1:length(nr.chg.pts) )
{
 if(nr.chg.pts[i]==2)
 {
     if( identical(as.integer(pos.chg.pts.list[[i]]),as.integer(c(21,90)) ) ){
       write(params.list[[i]],"params.real-data.d",append=T,ncol=6)
       vector.param1[count:(count+1)] <- exp(params.list[[i]][1])
       vector.param2[count:(count+1)] <- exp(params.list[[i]][2])
       vector.param3[count:(count+1)] <- sqrt(exp(params.list[[i]][3]))
       vector.param4[count:(count+1)] <- exp(params.list[[i]][4])
       vector.param5[count:(count+1)] <- exp(params.list[[i]][5])
       vector.param6[count:(count+1)] <- sqrt(exp(params.list[[i]][6]))
       vector.param7[count:(count+1)] <- exp(params.list[[i]][7])
       vector.param8[count:(count+1)] <- exp(params.list[[i]][8])
       vector.param9[count:(count+1)] <- sqrt(exp(params.list[[i]][9]))           
       count <-count+2
     }
 }
}

alarm()

# V1
hist(vector.param1,xlab="V", main="WIPP data. First branch
Conditional on Change Points Located at 100 and 209 days")

postscript("real-data-branch1-v.ps")

density.v1 <- density(vector.param1,bw=0.0005)
plot(density.v1,xlab="V", main="WIPP data. First branch
Conditional on Change Points Located at 100 and 209 days")

graphics.off()

# D1
hist(vector.param2,xlab="D", main="WIPP data. First branch
Conditional on Change Points Located at 100 and 209 days ")

postscript("real-data-branch1-D.ps")

density.d1 <- density(vector.param2,bw=0.0005)
plot(density.d1,xlab="D", main="WIPP data. First branch
Conditional on Change Points Located at 100 and 209 days ")

graphics.off()

# sigma1
hist(vector.param3,xlab="sigma", main="WIPP data. First branch 
Conditional on Change Points Located at 100 and 209 days")

postscript("real-data-branch1-sigma.ps")

density.sigma1 <- density(vector.param3,bw=0.007)
plot(density.sigma1,xlab="sigma", main="WIPP data. First branch 
Conditional on Change Points Located at 100 and 209 days")

graphics.off()

# V2
hist(vector.param4,xlab="V", main="WIPP data. Second branch 
Conditional on Change Points Located at 100 and 209 days")

postscript("real-data-branch2-v.ps")

density.v2 <- density(vector.param4,bw=0.0001)
plot(density.v2,xlab="V", main="WIPP data. Second branch 
Conditional on Change Points Located at 100 and 209 days")

graphics.off()

# D2
hist(vector.param5,xlab="D", main="WIPP data. Second branch 
Conditional on Change Points Located at 100 and 209 days")

postscript("real-data-branch2-D.ps")

density.d2 <- density(vector.param5,bw=0.003)
plot(density.d2,xlab="D", main="WIPP data. Second branch 
Conditional on Change Points Located at 100 and 209 days")

graphics.off()

# sigma2
hist(vector.param6,xlab="sigma", main="WIPP data. Second branch 
Conditional on Change Points Located at 100 and 209 days")

postscript("real-data-branch2-sigma.ps")

density.sigma2 <- density(vector.param6,bw=0.002)
plot(density.sigma2,xlab="sigma", main="WIPP data. Second branch 
Conditional on Change Points Located at 100 and 209 days")

graphics.off()

# V3
hist(vector.param7,xlab="V", main="WIPP data. Third branch 
Conditional on Change Points Located at 100 and 209 days")

postscript("real-data-branch3-v.ps")

density.v3 <- density(vector.param7,bw=0.05)
plot(density.v3,xlab="V", main="WIPP data. Third branch 
Conditional on Change Points Located at 100 and 209 days")

graphics.off()

# D3
hist(vector.param8,xlab="D", main="WIPP data. Third branch 
Conditional on Change Points Located at 100 and 209 days")

postscript("real-data-branch3-D.ps")

density.d3 <- density(vector.param8,bw=0.5)
plot(density.d3,xlab="D", main="WIPP data. Third branch 
Conditional on Change Points Located at 100 and 209 days")

graphics.off()

# sigma3
hist(vector.param9,xlab="sigma", main="WIPP data. Third branch 
Conditional on Change Points Located at 100 and 209 days")

postscript("real-data-branch3-sigma.ps")

density.sigma3 <- density(vector.param9,bw=0.002)
plot(density.sigma3,xlab="sigma", main="WIPP data. Third branch 
Conditional on Change Points Located at 100 and 209 days")

graphics.off()

##############################
# Comparing original WIPP curve with a curve obtinaed with using the mean
# of the posterior distributions for each of the parameters
##############################

plot(time.points,y,col="red",
main="WIPP data. Change Point Located at 100 and 209 days
Black line - curve obtained with posterior mean values
Red circles - experimental data")
lines(time.points[1:21],lnC.new(time.points[1:21],mean(log(vector.param1)),mean(log(vector.param2))))
lines(time.points[21:90],lnC.new(time.points[21:90],mean(log(vector.param4)),mean(log(vector.param5))))
lines(time.points[90:124],lnC.new(time.points[90:124],mean(log(vector.param7)),mean(log(vector.param8))))

######
# at 21
#######

n.monitor <-100000
tmp.pos.chgpt <- pos.chg.pts.list
vector1.param1 <- NULL
vector1.param2 <- NULL
vector1.param3 <- NULL
vector1.param4 <- NULL
vector1.param5 <- NULL
vector1.param6 <- NULL

count <-1
for( i in 1:length(nr.chg.pts) )
{
 if(nr.chg.pts[i]==1)
 {
     if( identical(pos.chg.pts.list[[i]],c(21) ) ){
       write(params.list[[i]],"params.real-data.d",append=T,ncol=6)
       vector1.param1[count:(count+1)] <- params.list[[i]][1]
       vector1.param2[count:(count+1)] <- params.list[[i]][2]
       vector1.param3[count:(count+1)] <- params.list[[i]][3]
       vector1.param4[count:(count+1)] <- params.list[[i]][4]
       vector1.param5[count:(count+1)] <- params.list[[i]][5]
       vector1.param6[count:(count+1)] <- params.list[[i]][6]          
       count <-count+2
     }
 }
}



plot(time.points,y,col="red",main="WIPP data. Change Point Located at 16 
Black circles - curve obtained with posterior mean values
Red circles - experimental data")
points(time[1:16],lnC.new(time.points[1:16],mean(vector1.param1),mean(vector1.param2)))
points(time[16:124],lnC.new(time.points[16:124],mean(vector1.param4),mean(vector1.param5)))
