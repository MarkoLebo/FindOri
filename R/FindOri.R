#' @title Detect the oriC region of a bacterial genome
#'
#' @description A function, that is used for oriC region detection in bacterial genomes
#' @param sek a vector with nucleotide sequence of genome, saved in DNAString.
#' @author Marko Lebo <markolebo48@gmai.com>
#' @export
#' @return returns a number of values in a \code{list}. In the first element of the \code{list} are parameters of the detected oriC region. The order of the parameteres are:
#'
#'      \code{OriStart} - position genome, where the oriC region starts,
#'
#'      \code{OriEnd} - position genome, where the oriC region starts,
#'
#'      \code{OriLenght} - length of the detected oriC,
#'
#'      \code{OriGC} - percentual GC-content of the oriC region,
#'
#'      \code{OriAT} - percentual AT-content of the oriC region,
#'
#'      \code{Min} - position of minimum of  GC-disparity curve,
#'
#'      \code{DnaAcount} - number of detected DnaA-boxes in the detected oriC region.
#'
#'      The second element of the returning \code{list} is a vector called \code{DistDnaA} in which the distances of the beginning of each
#'      detected DnaA-box from the beginning of oriC is returned.
#'
#'      The last element of the returned \code{list} is the whole nucleotide sequence of the detected oriC region.
#'
#'      This function also displays the GC-disparity curve of the input genome.
#' @examples \dontrun{
#' sek<-readDNAStringSet('sequence.fasta')
#' sek<-sek[[1]]
#' oriC <- FindOri(sek)
#' }

FindOri<-function(sek){
  sekchar <- unlist(strsplit(as.character(sek), split=''))
  done<-'Genome already starts in oriC region'
  Min<-Zcurve(sek)
  if (Min<(length(sek)-4000) && Min>4000){
    Start<-OriStart(Min,sekchar)
    Faza<-fazac(sekchar)
    DnaA<-DnaABoxy(Start,Faza)
    Poz<-DnaA[[1]]
    OriEnd<-DnaA[[2]]
    Poz<-Poz-Start
    Ori<-as.character(sek[Start:OriEnd])
    Ori<-toupper(Ori)
    OriString<-sek[Start:OriEnd]
    pom<-alphabetFrequency(OriString)
    GCcontent<-(pom[2]+pom[3])/(pom[1]+pom[4]+pom[3]+pom[2])
    ATcontent<-(pom[1]+pom[4])/(pom[1]+pom[4]+pom[3]+pom[2])
    OriLength<-OriEnd-Start
    vysledky<-matrix(0,nrow=1,ncol=7)
    colnames(vysledky)<-c('OriStart',' OriEnd',' OriLength',' OriGC',' oriAT',' Min',' DnaAcount')
    vysledky[1,1]<-Start
    vysledky[1,2]<-OriEnd
    vysledky[1,3]<-OriLength
    vysledky[1,4]<-GCcontent
    vysledky[1,5]<-ATcontent
    vysledky[1,6]<-Min
    vysledky[1,7]<-length(Poz)
    poz<-matrix(0,nrow=1,ncol=length(Poz))
    for ( i in 1:length(Poz)){
      poz[1,i]<-Poz[i]
    }
    rownames(poz)<-'DistDnaA'
    return(list(vysledky,poz,Ori[[1]]))
  } else {return(done)}
}
