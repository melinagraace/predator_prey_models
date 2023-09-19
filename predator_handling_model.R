###predator density model
a <- 0.2
b <- 0.15
c <- 0.3
d <- 0.15
K <- 10
h <- 0.5

N.start <- d/(c*b-d*b*h)*1.1
P.start <- (a*c/(c*b-d*b*h))*(1-d/(c*b*K-d*b*h*K))*1.1

predator <- function(a,b,c,d,K,h,N.start,P.start){
  
  #equilibrium points
  N.star1 <- 0
  P.star1 <- 0
  N.star2 <- K
  P.star2 <- 0
  N.star4 <- d/(c*b-d*b*h)
  P.star4 <- (a*c/(c*b-d*b*h))*(1-d/(c*b*K-d*b*h*K))
  
  #isoclines
  N.points <- seq(0, 1.5*K, by=0.1)
  P.points <- seq(0, 1.5*K, by=0.1)
  N.isocline <- (a/b)*(1-N.points/K)*(1+b*h*N.points)
  
  plot(N.points,N.isocline, type= 'l', col='blue', lwd=2,
       xlab='Prey (N)', ylab='Predator (P)', 
       main=paste('Predator Handling Time Model, K=',K),
       xlim=c(0,K), ylim=c(0,3))
  legend('topright', inset=0.05, legend=c('Prey(N)', 'Predator(P)', 'Trajectory'), 
         col=c('blue','red','black'), lty=1)
  
  #red is P isoclines
  abline(h=0, col="red", lwd=2)
  abline(v=d/(c*b-d*b*h), col='red', lwd=2)
  #blue is N isoclines
  abline(v=0, col='blue', lwd=2)
  
  
  points(N.star1,P.star1, pch = 21, bg = "white", cex = 1.5)
  points(N.star2,P.star2, pch = 21, bg = "white", cex = 1.5)
  points(N.star4,P.star4, pch = 21, bg = "white", cex = 1.5)  
  
  #trajectory
  points(N.start,P.start,pch=21,bg="black",cex=1.5)
  
  time <- 100
  timelist<-seq(1,time)
  
  N.output<-matrix(NA,nrow=time,ncol=3)
  N.output[1,1:3]<-c(1,N.start,P.start)
  
  for (i in 2:time){
    N<-N.output[i-1,2]
    P<-N.output[i-1,3]
    N.output[i,2]<- N+a*N*(1-N/K)-(b*N*P/(1+b*h*N))
    N.output[i,3]<- P+(c*b*N*P/(1+b*h*N))-d*P
  }
  
  lines(N.output[,2],N.output[,3],col="black",lwd=2)
  
  
  #time population data
  #2 is N(blue) and 3 is P (red)
  plot(seq(1,time,1),N.output[,2], col='blue', ylim=c(0,10), type='l', 
       xlab='time (years)', ylab='Population', 
       main=paste('Predator and Prey populations over time \n (Predator Handling Model)
        K=',K
       ))
  lines(seq(1,time,1),N.output[,3],col='red')
  legend('topleft',legend=c('Prey(N)', 'Predator(P)'), 
         col=c('blue','red'), lty=1)
  
  
  #calculate eigenvalues
  jacob <- matrix(0,2,2)
  jacob[1,1] <- 1+a-2*a*N.star4/K-b*P.star4/((1+b*h*N.star4)^2) #df/dn
  jacob[1,2] <- -b*N.star4/(1+b*h*N.star4) #df/dp
  jacob[2,1] <- c*b*P.star4/((1+b*h*N.star4)^2)  #dg/dn
  jacob[2,2] <- 1+c*b*N.star4/(1+b*h*N.star4)-d #dg/dp
  
  print(jacob)
  ev<-eigen(jacob)
  print(ev$values)
  
  Eigen.Magnitude=(Re(ev$values[1])^2+Im(ev$values[1])^2)^0.5
  print(paste("Eigenvalue magnitude = ", round(Eigen.Magnitude, 4)))
  
}

predator(a,b,c,d,10,h,N.start,P.start)

predator(a,b,c,d,25,h,N.start,P.start)


