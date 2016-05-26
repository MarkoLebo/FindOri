#' @title Compute and display the cumulted phase signal of the input genome
#'
#' @description A function, that is used for computation and displaying the cumulted phase signal of the input genome
#' @param sek a vector with nucleotide sequence of genome, saved in DNAString.
#' @author Marko Lebo <markolebo48@gmai.com>
#' @export
#' @return Displays the cumulted phase signal of the input genome.
#' @examples \dontrun{
#' sek<-readDNAStringSet('sequence.fasta')
#' sek<-sek[[1]]
#' cumulphase <- cumphase(sek)
#' }
cumphase<-function(sek){
  sekchar <- unlist(strsplit(as.character(sek), split=''))
  c<-length(sekchar)
  faza<-vector(mode='integer', length=c)
  faza[which(sekchar=='A')]<-pi/4
  faza[which(sekchar=='C')]<--3*pi/4
  faza[which(sekchar=='G')]<-3*pi/4
  faza[which(sekchar=='T')]<--pi/4
  kumfaza<-cumsum(faza)
  x<-seq(from=1,to=length(sekchar))
  plot(x,kumfaza[seq(from=1, to=length(kumfaza),by=10000)],type='l',main='Cumulated Phase',xlab='position in genome [bp]',ylab='phase[rad]')
}
