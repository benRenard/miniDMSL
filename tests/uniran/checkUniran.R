u=read.table('uniranDeviates.txt')[,1]

# function implementing the test for uniformity (KS)
unifTest <- function(u){
  w=ks.test(x=u,y='punif')
  return(w$p.value)
}

# function to split a vector of deviates into 1, 2, ..., 1024 blocks of identical size
getBlocks <- function(u,nBlocks=2^(0:10)){
  n=length(u)
  out=vector(mode='list',length=length(nBlocks))
  for(i in 1:length(out)){
    out[[i]]=matrix(u,nrow=trunc(n/nBlocks[i]),ncol=nBlocks[i])
  }
  return(out)
}

# Pack u into blocks of various sizes
Bu=getBlocks(u)

# Open pdf result file
h=7;w=h*16/9
pdf('checkUniran.pdf',height=h,width=w,useDingbats=F)

# Test uniformity within each block
par(mfrow=c(3,4),mai=rep(0.5,4),omi=rep(0.5,4))
for(i in 1:length(Bu)){
  pv=apply(Bu[[i]],2,unifTest) # p-values of the test for each block
  tit=paste('Block length:',NROW(Bu[[i]]))
  plot(pv,ylim=c(0,1),xlab='block',ylab='Unif. test pvalue',main=tit,pch=19)
}

# Test non-autocorrelation
ci.level=0.95
par(mfrow=c(3,4),mai=rep(0.5,4),omi=rep(0.5,4))
for(i in 1:length(Bu)){
  M=Bu[[i]]
  col=sample(1:NCOL(M),1) # pick one block at random
  z=acf(M[,col],lag.max=min(1000,NROW(M)),plot=F) # autocorrelation function
  cor=z$acf[2:length(z$acf)] # remove lag 0
  CI=c(1,-1)*qnorm(0.5*(1-ci.level),sd=sqrt(1/z$n.used)) # compute CI around 0 under an iid assumption
  percentOut=100*sum(abs(cor)>max(CI))/length(cor) # Count autocor. coeff. outside of the CI
  tit=paste('Block length:',NROW(M),'- Outside CI:',percentOut,'%')
  plot(z$acf[2:length(z$acf)],type='h',xlab='lag',ylab='cor.',main=tit)
  abline(h=CI,col='red')
}

# Visualize Spectrums
par(mfrow=c(3,4),mai=rep(0.5,4),omi=rep(0.5,4))
for(i in 1:length(Bu)){
  M=Bu[[i]]
  tit=paste('Block length:',NROW(M))
  col=sample(1:NCOL(M),1) # pick up one block at random
  spans=NROW(M)/10^3;if(spans<2) {spans=NULL} # Smoothing of the periodogram ()
  z=spectrum(M[,col],spans=spans,main=tit)
}

# Visualize bitmaps ----
toBitmap <- function(M,main=''){
  n=dim(M)
  image(M,col=hcl.colors(6,'Blues'),asp=1,axes=F,
        main=paste(main,n[1],'x',n[2]))
}
fourBitmaps <- function(M){
  toBitmap(M,'u -') # raw deviates
  toBitmap(M>0.1,'u>0.1 -') # transform to 0-1
  toBitmap(M>0.5,'u>0.5 -') # transform to 0-1
  toBitmap(M>0.9,'u>0.9 -') # transform to 0-1
}
n=length(u)
k=sqrt(n) # e.g. 1024*1024
M=matrix(u,k,k) # pack to matrix
C=array(u,dim=c(k/2,k/2,4)) # pack to 4 smaller matrices
par(mfrow=c(3,8),mai=rep(0.1,4),omi=rep(0.1,4))
fourBitmaps(M)
for(i in 1:4){fourBitmaps(C[,,i])}

# Close pdf result file
dev.off() 
