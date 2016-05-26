OriStart<-function(min,sekchar){
  prirastok<-0
  end<-1000
  oristart<-0
  repeat{
    if (end>2000){
      peaks<-findPeaks(AT)
      peakhight<-AT[peaks]
      peaks<-which(peakhight>=max(AT)*0.9)
      valleys<-findValleys(AT[peaks[length(peaks)]:length(AT)])
      valleys<-valleys[1]
      break}
    A<-0
    C<-0
    G<-0
    T<-0
    zacatek<-24
    konec<-24
    AT<-c()
    while (zacatek<=end){
      okno<-sekchar[(min-zacatek):(min+konec)]
      for (i in 1:length(okno)){
        if (okno[i]=='A'){
          A<-A+1}
        if (okno[i]=='C'){
          C<-C+1    }
        if (okno[i]=='G'){
          G<-G+1    }
        if (okno[i]=='T'){
          T<-T+1    }
      }
      ATpom<-(A+T)/(G+C)
      AT<-c(AT,ATpom)
      konec<-konec-1
      zacatek<-zacatek+1
    }
    if (max(AT)<1){min<-min+100
    prirastok<-prirastok+100}
    else if (min(AT)>1){end<-end+1000
    } else {break}
  }
  # x<-seq(from=min, to=min+length(AT)-1,by=1)
  # plot(x,AT,type='l',main='AT',ylab='AT [-]',xlab='position in genome [bp]')
  if (end<2000){
    prirastok<-prirastok-which.max(AT)
    AT<-AT[which.max(AT):length(AT)]
    hodnoty<-which(AT>=1)
    prirastok<-oristart+prirastok-hodnoty[length(hodnoty)]
    oristart<-oristart+min+prirastok
  } else {oristart<-min+valleys}
  return(oristart)
}
