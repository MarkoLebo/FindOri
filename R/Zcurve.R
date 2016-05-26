#' @title Compute and display the GC-disparity curve
#'
#' @description A function, that computes the GC-disparity curve from Z-curve of the input genome and displays it consecutively.
#' @param sek a vector with nucleotide sequence of genome, saved in DNAString
#' @author Marko Lebo <markolebo48@gmai.com>
#' @export
#' @return returns the value of minimum of the GC-disparity curve and displays the plot of the same curve
#' @examples \dontrun{
#' sek<-readDNAStringSet('sequence.fasta')
#' sek<-sek[[1]]
#' Min <- Zcurve(sek)
#'}
Zcurve<-function(sek){

  c<-length(sek)

  # vytvorenie vektorov četností nukleotidov
  A<-vector(mode='integer', length=c)
  C<-vector(mode='integer', length=c)
  G<-vector(mode='integer', length=c)
  T<-vector(mode='integer', length=c)
  sek <- unlist(strsplit(as.character(sek), split=''))
  A[which(sek=='A')]<-1
  C[which(sek=='C')]<-1
  G[which(sek=='G')]<-1
  T[which(sek=='T')]<-1

  # kumulovane sumy
  cumA<-cumsum(A)
  cumC<-cumsum(C)
  cumG<-cumsum(G)
  cumT<-cumsum(T)

  # vypocet Z-krivky
  Z<-matrix(0,length(sek),3)
  Z[,1]=(cumA+cumG)-(cumC+cumT)
  Z[,2]=(cumA+cumC)-(cumG+cumT)
  Z[,3]=(cumA+cumT)-(cumG+cumC)
  # projekcia Z-krivky
  GC<-(Z[,1]-Z[,2])/2
  # AT<-(Z[,1]+Z[,2])/2
  # doplnit osu x pre spravne rozlisenie
  x<-seq(from=1,by=10000,to=length(GC))
  plot(x,GC[seq(from=1,to=length(GC),by=10000)],type='l',main='GC-disparity',xlab='position in genome [bp]',ylab='GC [-]')
  # plot(x,AT[seq(from=1,to=length(AT),by=10000)],type='l',main='AT-disparity',xlab='pozícia v genóme [bp]',ylab='AT [-]')
  min<-which.min(GC)
  return(min)
}
