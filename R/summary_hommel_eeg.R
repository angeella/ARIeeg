#' @title summary cluster EEG
#' @description gives summary for each cluster EEG from clusterMass method
#' @usage summary_hommel_eeg(hommel,ix,alpha, clusters)
#' @param hommel hommel object from hommel package
#' @param ix set of interest
#' @param alpha alpha level
#' @param clusters id clusters
#' @author Angela Andreella
#' @return Returns a list with the following objects: discoveries number of discoveries in the set selected, cluster id, p-value
#' @export
#' @importFrom hommel discoveries



summary_hommel_eeg <- function(hommel,ix,alpha, clusters,eff){
  Total=length(hommel@p[ix])
  False_Null= discoveries(hommel, alpha=alpha, ix=ix)
  True_Null=Total-False_Null
  Active_Proportion= False_Null/Total
  clustermass <- eval(parse(text=paste0("model$multiple_comparison$", eff, "$clustermass$cluster$clustermass[clusters]")))
  pvalue <- eval(parse(text=paste0("model$multiple_comparison$", eff, "$clustermass$cluster$pvalue[clusters]")))
  out=c(clusters,Total, clustermass,pvalue,False_Null,True_Null,Active_Proportion)
  names(out)=c("ID","Total", "clustermass", "pvalue", "False Null", "True Null", "Active Proportion" )
  out
}


