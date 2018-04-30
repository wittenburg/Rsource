### general mapping functions assuming NCI (-> 0<=theta<=0.5)
tiff('map-functions.tiff',width=20,height=20,units='cm',res=200,pointsize=12); par(mfrow=c(2,2))


haldane=function(x){
  theta=c()
  for (i in 1:length(x)){
    theta[i]=0.5*(1-exp(-2*x[i]))
  }
  theta
}


curve(haldane,0,2,main="Haldane map function",xlab="Morgan units",ylab="recombination rate")



kosambi=function(x){
  theta=c()
  for (i in 1:length(x)){
    theta[i]=0.5*tanh(2*x[i])
  }
  theta
}

curve(kosambi,0,2,main="Kosambi map function",xlab="Morgan units",ylab="recombination rate")



# Rao 1977

w=function(p,x){
  y=c()
  for(i in 1:length(x)){
    theta=x[i]
    y[i]=(p*(2*p-1)*(1-4*p)*log(1-2*theta)+16*p*(p-1)*(2*p-1)*atan(2*theta)+
    2*p*(1-p)*(8*p+2)*atanh(2*theta)+6*(1-p)*(1-2*p)*(1-4*p)*theta)/6
  }
  y
}


r=seq(0,0.5,1e-3)
count=0
for (p in c(0.5,1,1.5)){
  count=count+1
  if(count==1) plot(w(p,r),r,type='l',xlim=c(0,2),main='Rao map function',xlab="Morgan units",ylab="recombination rate")
  if(count>1) lines(w(p,r),r,col=count)
}




# Felsenstein 1979

f=function(K,x){
  r=c()
  for(i in 1:length(x)){
    d=x[i]
    r[i]=(1-exp(2*(K-2)*d))/2/(1-(K-1)*exp(2*(K-2)*d))
  }
  r
}


x=seq(0,2,1e-3)
count=0
for (K in c(-1,0,1,1.5)){
  count=count+1
  if(count==1) plot(x,f(K,x),type='l',main='Felsenstein map function',xlab="Morgan units",ylab="recombination rate")
  if(count>1) lines(x,f(K,x),col=count)
}

dev.off()
