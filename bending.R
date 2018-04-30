# unweighted bending algorithm according to Jorjani et al. 2003
bend=function(K){
  K.neu=K
  epsilon=1e-4; counter=0
  repeat{
    counter=counter+1
    if(counter>1000) break
    if(counter%%100==1) cat('iteration',counter,'\n')
    e=eigen(K.neu,symmetric=T)
    d=e$values
    if(min(d)>0) break
    d[d<epsilon]=epsilon
    K.neu=e$vectors%*%diag(d)%*%t(e$vectors)
  }
  cat(counter,'iterations\n')
  K.neu
}
