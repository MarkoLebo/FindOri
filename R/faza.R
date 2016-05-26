fazac<-function(sekchar){

  c<-length(sekchar)
  faza<-vector(mode='integer', length=c)
  faza[which(sekchar=='A')]<-pi/4
  faza[which(sekchar=='C')]<--3*pi/4
  faza[which(sekchar=='G')]<-3*pi/4
  faza[which(sekchar=='T')]<--pi/4
  return(faza)
}