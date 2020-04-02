#' @title ARI Permutation-based for EEG data
#' @description Performs ARI using permutation local test for EEG data
#' @usage ARIpermCT(data, alpha, summary_stat, silent, family, dist, ct, time, ...)
#' @param data data
#' @param alpha alpha level
#' @param summary_stat Choose among "max", "center-of-mass"
#' @param silent FALSE by default.
#' @param family which family for the confidence envelope? simes, finner, beta or higher.criticism. default is simes
#' @param delta do you want to consider at least delta size set?
#' @param B number of permutation, default 1000
#' @param ct set of thresholds
#' @param time time signal to select, default NULL, i.e. no selection
#' @param dist an double defining the maximal distance for adjacency of two channels
#' @param formula formula object defining the design of the model.
#' @param var variables to use in the right part of the formula
#' @param B number of permutation
#' @author Angela Andreella
#' @return Returns a list with the following objects: discoveries number of discoveries in the set selected, cluster id, maximum test statistic and relative coordinates
#' @export
#' @importFrom plyr laply
#' @importFrom tidyr nest
#' @importFrom abind abind
#' @importFrom dplyr group_by
#' @importFrom signal_tbl eeguana
#' @importFrom dplyr mutate
#' @importFrom purrr  invoke
#' @importFrom eeguana segments_tbl
#' 
ARIeeg <- function(data, alpha = 0.1, family = "Simes", delta = 0, ct = c(0,1), alternative = "two.sided",time = NULL,dist = 50,formula,var, B = 5000,...){
  
  if(!is_eeg_lst(data)){
    data <- utilsTOlst(data, ...)
  }
  
  #Data from normal or segmented?? <-
  signal <- 
    data%>%
    signal_tbl()%>%
    group_by(.id)%>%
    nest()%>%
    mutate(data = map(data,~as.matrix(.x[-1])))%>%
    pull(data)%>%
    invoke(abind,.,along = 3)%>%
    aperm(c(3,1,2))
  
  if(!is.null(time)){signal <- signal[,time,]}
  
  design <- 
    segments_tbl(data)%>%
    select(var)

  graph <- position_to_graph(channels_tbl(data), name = .channel, delta = dist,
                             x = .x, y = .y, z = .z)
  
  #formula <- signal ~ condition + Error(.subj/(condition))
  
  model <- permuco4brain::brainperm(formula = formula,
                                    data = design,
                                    graph = graph,
                                    np = B,
                                    method = NULL,
                                    type = "permutation",
                                    test = "fisher",
                                    aggr_FUN = NULL,
                                    threshold = NULL,
                                    multcomp = "clustermass",
                                    effect = NULL,
                                    return_distribution = TRUE)
  
  
  Test <- model$multiple_comparison$stimuli$uncorrected$distribution
  dim(Test) <- c(dim(Test)[1], dim(Test)[2]*dim(Test)[3])
  pvalues <- switch(alternative, 
                    "two.sided" = 2*(pnorm(abs(Test), lower.tail=FALSE)),
                    "greater" = pnorm(Test, lower.tail=FALSE),
                    "less" = 1-pnorm(Test, lower.tail=FALSE))
  test_obs <- model$multiple_comparison$stimuli$clustermass$data$statistic
  pvalues_obs <- switch(alternative, 
                    "two.sided" = 2*(pnorm(abs(test_obs), lower.tail=FALSE)),
                    "greater" = pnorm(test_obs, lower.tail=FALSE),
                    "less" = 1-pnorm(test_obs, lower.tail=FALSE))

  pvalues <-rbind(pvalues_obs,pvalues)
  #pvalues <- pvalues[,which(model$multiple_comparison[[1]]$clustermass$data$cluster_id!=0)]
  pvalues_ord <- rowSortC(pvalues)
  lambda <- lambdaOpt(pvalues = pvalues_ord, family = family, ct = ct, alpha = alpha, delta = delta) 
  cvOpt = cv(pvalues = pvalues_ord, family = family, alpha = alpha, lambda= lambda, delta = delta)
  
 #clstr_id <- model$multiple_comparison[[1]]$clustermass$data$cluster_id[which(model$multiple_comparison[[1]]$clustermass$data$cluster_id!=0)]
  clstr_id <- model$multiple_comparison$stimuli$clustermass$data$cluster_id
  #clusters <- c(1:model$multiple_comparison[[1]]$clustermass$cluster$no)[which(model$multiple_comparison[[1]]$clustermass$cluster$pvalue<=0.1)]
  #clusters <- which(model$multiple_comparison[[1]]$clustermass$cluster$pvalue<=0.1)
  clusters <- c(1:model$multiple_comparison$stimuli$clustermass$cluster$no)
  out=lapply(clusters,function(i){
    ix= which(clstr_id == i)
    
    #cluster_ids=which(ix,arr.ind = TRUE)
    #cluster_ids=cbind(cluster_ids,Stat=StatFun(ix))
    #Error if I put pvalues[,mask] instead of pvalues in SingleStepCT
    #perm <- SingleStepCT(pvalues = pvalues,ct =ct, ix =as.vector(which(ix[mask])), alpha = alpha, shift = shift, family = 'Simes', lambda = lambda)
    #perm <- discoveriesPerm(praw = praw, ix = ix[mask], cvh = cvh)
    summary_cluster_eeg(clusters = i,model = model, 
                        cv = cvOpt,ix=ix,pvalues = pvalues)
  
  })
  
  out_d <- t(as.data.frame(out))
  rownames(out_d) = NULL
  out_d
}

