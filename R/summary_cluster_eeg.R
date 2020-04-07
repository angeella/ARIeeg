#' @title summary cluster EEG
#' @description gives summary for each cluster EEG from clusterMass method
#' @usage summary_cluster_eeg(clusters,model,cv,ix,pvalues)
#' @param clusters id cluster
#' @param model clusterMass model from permuco4brain package
#' @param cv critical vector
#' @param ix set of interest
#' @param pvalues pvalues matrix where rows indicate the permutations
#' @author Angela Andreella
#' @return Returns a list with the following objects: discoveries number of discoveries in the set selected, cluster id, p-value
#' @export
#' @importFrom ARIpermutation dI




summary_cluster_eeg <- function(clusters,model,cv,ix,pvalues){
  p <- pvalues[1,ix]
  Total = length(p)
  False_Null= dI(ix = ix,cv = cv,praw = pvalues[1,])
  True_Null=Total - False_Null
  Active_Proportion= False_Null / Total
  #info <- eval(parse(text=paste0("model$multiple_comparison$", effect, "$clustermass$cluster")))
  clustermass <- eval(parse(text=paste0("model$multiple_comparison$", effect, "$clustermass$cluster$clustermass[clusters]")))
  pvalue <- eval(parse(text=paste0("model$multiple_comparison$", effect, "$clustermass$cluster$pvalue[clusters]")))

  #csize <- model$multiple_comparison$stimuli$clustermass$cluster$csize[clusters]
  #membership <- toString(unique(model$multiple_comparison$condition$clustermass$data$channel[model$multiple_comparison$condition$clustermass$data$cluster_id==1]))
  #out=c(clustermass)
  #names(out)=c("clustermass")
  out=c(clusters,Total, clustermass,pvalue,False_Null,True_Null,Active_Proportion)
  names(out)=c("ID","Total", "clustermass", "pvalue", "False Null", "True Null", "Active Proportion" )
  out
}
