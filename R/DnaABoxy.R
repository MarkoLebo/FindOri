DnaABoxy<-function(oristart,faza){
  DnaAfaza<-c(-0.7853982, -0.7853982,  0.7853982, -0.7853982, -2.3561945, -2.3561945,  0.7853982, -2.3561945,  0.7853982)
  end<-oristart+299
  prahfazahore<-0.75
  prahfazadole<--0.75
  repeat{
    part<-faza[oristart:end]
    r<-c()
    z<-1
    for (i in 1:length(part)){
      pom<-cor(part[z:(z+8)], DnaAfaza)
      r<-c(r,pom)
      z<-z+1
    }
    r<-r[-which(is.na(r)==TRUE)] 
    x<-seq(from=oristart, to=oristart+length(r)-1,by=1)
    plot(x,r,type='l',main='DnaA-box correlation', xlab='position in genome [bp]',ylab='r [-]')
    poz<-which(r>=prahfazahore)+oristart
    poz<-c(poz,which(r<=prahfazadole)+oristart)
    poz<-sort(poz)
    if (length(poz)<=2){
      prahfazahore<-prahfazahore-0.02
      prahfazadole<-prahfazadole+0.02
      next
    }
    rozdiely<-c()
    for (i in 2:length(poz)){
      pom<-poz[i]-poz[i-1]
      rozdiely<-c(rozdiely,pom)
    }
    minpocet<-4
    for (i in 1:length(rozdiely)){
      if (rozdiely[i]>100 && i<minpocet){
        poz<-poz[-c(1:i)]
      }
      if (rozdiely[i]>100 && i>minpocet){
        poz<-poz[-c((i+1):length(poz))]
      }
    }
    p<-length(part)-(poz[which.max(poz)]-oristart)
    if (p<100){end<-end+(100-p)} else {break}
    if (end>1000){break}
  }
  
  oriend<-poz[which.max(poz)]
  return(list(c(poz),oriend))
}