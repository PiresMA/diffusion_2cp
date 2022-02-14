
##******************************************************************##
#  Coarse-grained Allee-Diffusion model (multi-groups)
##******************************************************************##


tmax       = 201  # number of time steps to simulate over
nsubp      = L = 10
sizesubpop = 10^5
N          = sizesubpop*nsubp
alp        = 0.9  
lam        = 1.0
r          = 0.86
fao        = 0.48 # fraction of ocupation a t0  
fbo        = 1-fao

###***************** import  ********************#####
vecn   = round(rep(N/nsubp,nsubp))
###***************** import  ********************#####

simu    = 4

nsourc  = 1
fatorIo = 1


vecD    = c( 0.1 )
#vecD    = sort( c( seq(0.0, 0.3, 0.1)) )


                   

for(indexD in 1:length(vecD) )
{


  D = vecD[indexD]
  
  
 prepatt =  sprintf("2cp_simu%d_nS%d_r%1.3f_alp%1.3f_N%d_L%d_D%1.2f",simu,nsourc,r,alp,N,L,D)

  
 #exportFigure   <- FALSE
 exportFigure   <- TRUE
 if(exportFigure==TRUE){
    nameOUT=gsub("_", " ",prepatt) 
    nameOUT=gsub("[[:punct:]]", "", nameOUT)
    nameOUT=gsub(" ", "-", nameOUT)
    nameOUT=paste(nameOUT,".png",sep="") 
    png(nameOUT) 
 }
 

  
  if(D==0) ylim=c(0, 1*sizesubpop) 
  else ylim=c(0, 1*sizesubpop)



  xlim=c(0,tmax-1)
  
  
  
  plot(1,1,type="n",axes=FALSE,xlab=" ",ylab=" ",xlim=xlim,ylim=ylim)
  axis(side=1, cex.axis=1.5, at=seq(xlim[1],xlim[2],len=5) )
  

  axis(side=2, cex.axis=1.5, at=seq(ylim[1],ylim[2],len=5) ) 
        
  box()
  mtext(side=1,cex=1.5,line=2.5,"time (mcs)")
  mtext(side=2,cex=1.5,line=2.5,"Total number of individuals")


  
  
  
  
  
  
  
  
  
  
  #********** pre: Numerical Solution (EULER) ********************#####
  
  
  nn      = 2
  Dbar    = 1-D;  #mobility parameter
  maxIo   = 1/L

  vecAo = vecBo = c(NULL)
  vecAo = c( rep( fao*(N/L)*(1/nsourc), nsourc) , rep(0, L-nsourc) )
  vecBo = c( rep( fbo*(N/L)*(1/nsourc), nsourc) , rep(0, L-nsourc) )
  
  vecVo = vecn  - vecAo- vecBo
  soma  = vecVo + vecAo + vecAo
  print(data.frame(vecAo,vecVo,vecBo,soma))
  #print(data.frame(vecAo/vecn,vecVo/vecn,vecBo/vecn,soma/vecn))

  Aoglobal = sum(vecAo)/N
  Boglobal = sum(vecBo)/N
  cat("maxIo",maxIo,"fator",fatorIo*maxIo,"Aoglobal:",Aoglobal,"Boglobal:",Boglobal)

  
  # Homogeous landscape
  vecalpcommu = rep(alp,L) 
  veclamcommu = rep(lam,L) 
  
  
  neighbors_mat = matrix(NA,nrow=L,ncol=nn+1)
  vec1 = c(L,1:L,1) 
  label=1
  for(i  in 2:(L+1)){
    neighbors_mat[label,] = c(vec1[i],vec1[i-1], vec1[i+1])
    label=label+1;
  }
  #print(neighbors_mat)
  
  
  
  
  #********** now: Numerical Solution (EULER) ********************#####  
  vecnt = lA = lB = lV = list(NULL)
  for (pop in 1:L) {

    lV[[pop]] <- c(vecVo[pop])     
    lA[[pop]] <- c(vecAo[pop])     
    lB[[pop]] <- c(vecBo[pop])     
  }
  
  for( t in 1:tmax){
    for (u in 1:L){
      nutA = lV[[u]][t]+lA[[u]][t]+lB[[u]][t]  
      nutB = lV[[u]][t]+lB[[u]][t]+lA[[u]][t]  
      
      lameff = Dbar*veclamcommu[u] 
      alpeff = Dbar*vecalpcommu[u] 
      
      left  = neighbors_mat[u,2]
      right = neighbors_mat[u,3]
      
      diffuA = -D*lA[[u]][t] + (1/2)*D*(lA[[left]][t]+lA[[right]][t])
      diffuB = -D*lB[[u]][t] + (1/2)*D*(lB[[left]][t]+lB[[right]][t])
      
      lV[[u]][t+1] = lV[[u]][t] - lameff*lV[[u]][t]*( lA[[u]][t]/nutA + r*lB[[u]][t]/nutB ) + alpeff*( lA[[u]][t] + r*lB[[u]][t] )
      lA[[u]][t+1] = lA[[u]][t] +   lameff*lV[[u]][t]*lA[[u]][t]/nutA -   alpeff*lA[[u]][t] + diffuA
      lB[[u]][t+1] = lB[[u]][t] + r*lameff*lV[[u]][t]*lB[[u]][t]/nutB - r*alpeff*lB[[u]][t] + diffuB
    }
  }
  

  matIeulerA = matrix(unlist(lA),  ncol=length(lA)  )
  matIeulerB = matrix(unlist(lB),  ncol=length(lB)  )
  
  lines(apply(matIeulerA, 1, sum),type="l",col="red", lwd=3 )
  lines(apply(matIeulerB, 1, sum),type="l",col="blue",lwd=3 )
  








  vecindexA = seq(2,2+3*(L-1),by=3)





  print(prepatt)
  namefiles = dir(pattern=prepatt)

  #print(namefiles)

  for(ii in 1:length(namefiles) ){

    mat = data.matrix(read.table( namefiles[ii] ) )
    print(  namefiles[ii] )

    time = mat[,1]
    matA = mat[,vecindexA]
    matV = mat[,vecindexA+1]
    matB = mat[,vecindexA+2]
    
    lines(time, apply(matA, 1, sum), type="p", pch=0,col="red",  cex=0.69 )
    lines(time, apply(matB, 1, sum), type="p", pch=1,col="blue", cex=0.69 )  
  }


  mtext(side=3,cex=0.9,line=2,sprintf("L:%d  fao:%1.2f  fbo:%1.2f  r:%1.2f lam:%1.2f  alp:%1.2f",L,fao,fbo,r,lam,alp))  
  mtext(side=3,cex=1.5,line=0.1,sprintf("D:%1.2f",D))

#    mtext(side=3,cex=1.2,line=1.9,bquote( fa[o]==.(fao)~~fb[o]==.(fbo)~~r==.(r)  ) )  
#    mtext(side=3,cex=1.3,line=0.1,bquote( D==.(D) ) )

    if(D==0) positionLeg="topright"
    else positionLeg="bottomright"
 
    lamA <- format(lam,digits=3L);
     
    legend( positionLeg, pt.cex=1, cex=1.3, bty="n", horiz=F,
       #lty = c(NA), lwd = c(NA),
       pch = c(15,16),
       col = c("red","blue"),
       legend = c(  parse(text=sprintf('lambda[A] == %s ~ alpha[A] == %s ~ (fast)' ,lam,alp  )), 
                    parse(text=sprintf('lambda[B] == %s ~ alpha[B] == %s ~ (slow)',r*lam,r*alp)) ) ) 

if(exportFigure==TRUE){
    dev.off()
}
#####*********************** Figure: end ***************************###

}



