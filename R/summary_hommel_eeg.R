summary_hommel_eeg <- function(hommel,ix,alpha, clusters){
  Total=length(hommel@p[ix])
  False_Null=hommel::discoveries(hommel, alpha=alpha, ix=ix)
  True_Null=Total-False_Null
  Active_Proportion= tdp1(hommel, ix=ix, alpha = alpha)
  clustermass <- eval(parse(text=paste0("model$multiple_comparison$", effect, "$clustermass$cluster$clustermass[clusters]")))
  pvalue <- eval(parse(text=paste0("model$multiple_comparison$", effect, "$clustermass$cluster$pvalue[clusters]")))
  out=c(clusters,Total, clustermass,pvalue,False_Null,True_Null,Active_Proportion)
  names(out)=c("ID","Total", "clustermass", "pvalue", "False Null", "True Null", "Active Proportion" )
  out
}


tdp1 <- function (hommel, ix, alpha) 
{
  m <- length(hommel@p)
  if (missing(ix)) {
    d <- discoveries(hommel, alpha = alpha)
    k <- m
  }
  else {
    p <- hommel@p[ix]
    k <- length(p)
    d <- discoveries(hommel, ix, alpha = alpha)
  }
  d/k
}