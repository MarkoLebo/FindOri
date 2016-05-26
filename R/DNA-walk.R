#' @title Compute and display the DNA-walk of the input genome
#'
#' @description A function, that is used for computation and displaying the DNA-walk of the input genome
#' @param sek a vector with nucleotide sequence of genome, saved in DNAString.
#' @author Marko Lebo <markolebo48@gmai.com>
#' @export
#' @return Displays the DNA-walk of the input genome.
#' @examples \dontrun{
#' sek<-readDNAStringSet('sequence.fasta')
#' sek<-sek[[1]]
#' DNAWalk <- Dnawalk(sek)
#' }
Dnawalk<-function(sek){
  A<-c(0,1,1)
  T<-c(0,-1,1)
  G<-c(1,0,1)
  C<-c(-1,0,1)
  walk<-matrix(0,length(sek),3)
  okno<-1000
  zaciatok<-1
  for (i in 2:length(sek)){
    if (sek[i]=='A'){
      walk[i,c(1:ncol(walk))]=walk[i,c(1:ncol(walk))]+A
    }
    if (sek[i]=='C'){
      walk[i,c(1:ncol(walk))]=walk[i,c(1:ncol(walk))]+C
    }
    if (sek[i]=='G'){
      walk[i,c(1:ncol(walk))]=walk[i,c(1:ncol(walk))]+G
    }
    if (sek[i]=='T'){
      walk[i,c(1:ncol(walk))]=walk[i,c(1:ncol(walk))]+T
    }}
x<-matrix(c(walk[c(1:nrow(walk))]),1,ncol=1)
y<-matrix(c(walk[c(1:nrow(walk))]),2,ncol=1)
z<-matrix(c(walk[c(1:nrow(walk))]),3,ncol=1)
plot(z,x,type='l',main='DNA-walk',xlab='position in genome[bp]',ylab='DNA-walk[-]')
}
