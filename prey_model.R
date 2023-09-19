###prey density model
a <- 0.2
b <- 0.15
c <- 0.1
d <- 0.15
K <- 15
h <- 0.5

N.start <- d/(c*b)*1.1
P.start <- (a*(1-(d/(c*b))/K)/b)*1.1

prey.density <- function(a,b,c,d,K,h,N.start,P.start){
  
  #equilibrium points
  N.star1 <- 0
  P.star1 <- 0
  N.star2 <- K
  P.star2 <- 0
  N.star3 <- d/(c*b)
  P.star3 <- (a*(1-(d/(c*b))/K)/b)
  
  #isoclines
  N.points <- seq(0, 1.5*K, by=0.1)
  P.points <- seq(0, 1.5*K, by=0.1)
  N.isocline.2 <- (a*(1-(N.points/K)))/b
  
  plot(N.points,N.isocline.2, type= 'l', col='blue', lwd=2,
       xlab='Prey (N)', ylab='Predator (P)', 
       main='Prey Density Dependance Model \n Zoomed in on equilibrium point',
       #xlim=c(0,K),ylim=c(0,3))
       xlim=c(9,11.5), ylim=c(0.35,0.6)) #0,K, 0,3
  #legend('topleft', inset=0.05, legend=c('Prey(N)', 'Predator(P)', 'Trajectory'), 
  #       col=c('blue','red', 'black'), lty=1)
  
  #red is P isoclines
  abline(h=0, col="red", lwd=2)
  abline(v=d/(c*b), col='red', lwd=2)
  #blue is N isoclines
  abline(v=0, col='blue', lwd=2)

  
  points(N.star1,P.star1, pch = 21, bg = "white", cex = 1.5)
  points(N.star2,P.star2, pch = 21, bg = "white", cex = 1.5)
  points(N.star3,P.star3, pch = 21, bg = "white", cex = 1.5)
  
  
  #trajectory
  points(N.start,P.start,pch=21,bg="black",cex=1.5)
  
  time <- 300
  timelist<-seq(1,time)
  
  N.output<-matrix(NA,nrow=time,ncol=3)
  N.output[1,1:3]<-c(1,N.start,P.start)
  
  for (i in 2:time){
    N<-N.output[i-1,2]
    P<-N.output[i-1,3]
    N.output[i,2]<- N+a*N*(1-N/K)-b*N*P
    N.output[i,3]<- P+c*b*N*P-d*P
  }
  
  lines(N.output[,2],N.output[,3],col="black",lwd=2)
  
  #time population data
  #2 is N(blue) and 3 is P (red)
  plot(seq(1,time,1),N.output[,2], col='blue', ylim=c(0,15), type='l', 
       xlab='time (years)', ylab='Population', 
       main='Predator and Prey populations over time \n (Prey density dependence model)')
  lines(seq(1,time,1),N.output[,3],col='red')
  legend('topleft',legend=c('Prey(N)', 'Predator(P)'), 
         col=c('blue','red'), lty=1)
  
  
  #calculate eigenvalues
  jacob <- matrix(0,2,2)
  jacob[1,1] <- 1+a-2*a*N.star3/K-b*P.star3 #df/dn
  jacob[1,2] <- -b*N.star3 #df/dp
  jacob[2,1] <- c*b*P.star3 #dg/dn
  jacob[2,2] <- 1+c*b*N.star3-d #dg/dp
  
  print(jacob)
  ev<-eigen(jacob)
  print(ev$values)
  
  Eigen.Magnitude=(Re(ev$values[1])^2+Im(ev$values[1])^2)^0.5
  print(paste("Eigenvalue magnitude = ", round(Eigen.Magnitude, 4)))

}

prey.density(a,b,c,d,K,h,N.start,P.start)