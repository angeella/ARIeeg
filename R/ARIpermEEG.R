#' @title ARI Permutation-based for EEG data
#' @description Performs ARI using permutation local test for EEG data
#' @usage ARIpermEEG(data, alpha, family, delta, alternative, timeS, dist, formula, variable, B, effect, model, ...)
#' @param data data, default is NULL. You need to insert data or model
#' @param alpha alpha level
#' @param family which family for the confidence envelope? simes, finner, beta or higher.criticism. default is simes
#' @param delta do you want to consider at least delta size set?
#' @param alternative alternative hypothesis
#' @param timeS time signal to select, default NULL, i.e. no selection
#' @param dist an double defining the maximal distance for adjacency of two channels
#' @param formula formula object defining the design of the model.
#' @param variable variables to use in the right part of the formula
#' @param B number of permutation, default 5000
#' @param eff effect of interest where apply ARI
#' @model permuco4brain model, default is NULL
#' @author Angela Andreella
#' @return Returns a data.frame with number of true discoveries for each cluster
#' @export
#' @importFrom plyr laply
#' @importFrom tidyr nest
#' @importFrom abind abind
#' @importFrom dplyr group_by
#' @importFrom eeguana signal_tbl
#' @importFrom dplyr mutate
#' @importFrom purrr invoke
#' @importFrom purrr map
#' @importFrom eeguana segments_tbl
#' @importFrom pARI lambdaOpt
#' @importFrom pARI criticalVector
#' @importFrom  permuco4brain brainperm
#' @importFrom permuco4brain position_to_graph
#' 
ARIpermEEG <- function(data=NULL, alpha = 0.1, family = "Simes", delta = 0, alternative = "two.sided",timeS = NULL,dist = 50,formula=NULL,variable=NULL, B = 5000,eff = "condition", model = NULL,...){
  
  if(!is.null(data) & !is.null(model)){stop("Please insert data or model object")}
  if(!is.null(data)){
    
  
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
  
  if(!is.null(timeS)){signal <- signal[,timeS,]}
  
  design <- 
    segments_tbl(data)%>%
    select(variable)

  graph <- position_to_graph(channels_tbl(data), name = .channel, delta = dist,
                             x = .x, y = .y, z = .z)
  
  #formula <- signal ~ condition + Error(.subj/(condition))
  
  model <- brainperm(formula = formula,
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
  }
  
  
  Test <- eval(parse(text=paste0("model$multiple_comparison$", eff, "$uncorrected$distribution")))
  dim(Test) <- c(dim(Test)[1], dim(Test)[2]*dim(Test)[3])

  #test_obs <- eval(parse(text=paste0("model$multiple_comparison$", eff, "$clustermass$data$statistic")))
  #TT <- rbind(test_obs,Test)
  pvalues <- switch(alternative, 
               "two.sided" =  matrixStats::colRanks(-abs(Test)) / nrow(Test),
               "greater" = matrixStats::colRanks(-Test) / nrow(Test),
               "less" = matrixStats::colRanks(Test) / nrow(Test))
  
  
  pvalues <-t(pvalues)
  #pvalues <- pvalues[,which(model$multiple_comparison[[1]]$clustermass$data$cluster_id!=0)]
  pvalues_ord <- rowSortC(pvalues)
  lambda <- lambdaOpt(pvalues = t(pvalues_ord), family = family, alpha = alpha, delta = delta) 
  cvOpt = criticalVector(pvalues = t(pvalues_ord), family = family, alpha = alpha, lambda= lambda, delta = delta)
  if(is.unsorted(cvOpt)){
    idS = which(sapply(c(1:length(cvOpt)), function(x) is.unsorted(cvOpt[1:x])))[1]
    cvOpt = c(cvOpt[1:(idS-1)], rep(length(cvOpt)-1, max(cvOpt[1:(idS-1)])))  
    
  }
 #clstr_id <- model$multiple_comparison[[1]]$clustermass$data$cluster_id[which(model$multiple_comparison[[1]]$clustermass$data$cluster_id!=0)]
  clstr_id <- eval(parse(text=paste0("model$multiple_comparison$", eff, "$clustermass$data$cluster_id")))
  
  #clusters <- c(1:model$multiple_comparison[[1]]$clustermass$cluster$no)[which(model$multiple_comparison[[1]]$clustermass$cluster$pvalue<=0.1)]
  #clusters <- which(model$multiple_comparison[[1]]$clustermass$cluster$pvalue<=0.1)
  clusters <- c(1:eval(parse(text=paste0("model$multiple_comparison$", eff, "$clustermass$cluster$no"))))
  #hom <-hommel(pvalues[1,])
  
  out=lapply(clusters,function(i){
    ix= which(clstr_id == i)
    
    #cluster_ids=which(ix,arr.ind = TRUE)
    #cluster_ids=cbind(cluster_ids,Stat=StatFun(ix))
    #Error if I put pvalues[,mask] instead of pvalues in SingleStepCT
    #perm <- SingleStepCT(pvalues = pvalues,ct =ct, ix =as.vector(which(ix[mask])), alpha = alpha, shift = shift, family = 'Simes', lambda = lambda)
    #perm <- discoveriesPerm(praw = praw, ix = ix[mask], cvh = cvh)
    summary_cluster_eeg(clusters = i,model = model, 
                        cv = cvOpt,ix=ix,pvalues = pvalues, eff = eff)
    #summary_hommel_eeg(hommel = hom,ix=ix, alpha = 0.1, clusters = i)
    
  
  })
  
  out_d <- t(as.data.frame(out))
  rownames(out_d) = NULL
  out_d
}


