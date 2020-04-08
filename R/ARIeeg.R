#' @title Parametric ARI for EEG data
#' @description Performs ARI using permutation local test for EEG data
#' @usage ARIeeg(data, alpha,alternative,timeS,dist,formula,variable, B, effect,...)
#' @param data data
#' @param alpha alpha level
#' @param alternative alternative hypothesis
#' @param timeS time signal to select, default NULL, i.e. no selection
#' @param dist an double defining the maximal distance for adjacency of two channels
#' @param formula formula object defining the design of the model.
#' @param variable variables to use in the right part of the formula
#' @param B number of permutation, default 5000
#' @param effect effect of interest where apply ARI
#' @author Angela Andreella
#' @return Returns a data.frame with number of true discoveries for each cluster
#' @export
#' @importFrom plyr laply
#' @importFrom tidyr nest
#' @importFrom abind abind
#' @importFrom dplyr group_by
#' @importFrom signal_tbl eeguana
#' @importFrom dplyr mutate
#' @importFrom purrr  invoke
#' @importFrom eeguana segments_tbl
#' @importFrom  permuco4brain brainperm
#' @importFrom permuco4brain position_to_graph
#' @importFrom hommel hommel
#' 
ARIeeg <- function(data, alpha = 0.1,alternative ="two.sided", timeS = NULL,dist = 50,formula,variable, B = 5000, effect = "condition",...){
  
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
    eeguana::segments_tbl(data)%>%
    select(variable)
  
  graph <- position_to_graph(channels_tbl(data), name = .channel, delta = dist,
                             x = .x, y = .y, z = .z)
  
  #formula <- signal ~ condition + Error(.subj/(condition))
  
  model <- brainperm(formula = formula,
                     data = design,
                     graph = graph,
                     np = B,
                     method = NULL,
                     type = "signflip",
                     test = "fisher",
                     aggr_FUN = NULL,
                     threshold = NULL,
                     multcomp = "clustermass",
                     effect = NULL,
                     return_distribution = TRUE)
  
  Test <- eval(parse(text=paste0("model$multiple_comparison$", eff, "$uncorrected$distribution")))
  dim(Test) <- c(dim(Test)[1], dim(Test)[2]*dim(Test)[3])

  test_obs <- eval(parse(text=paste0("model$multiple_comparison$",eff, "$clustermass$data$statistic")))

  TT <- rbind(test_obs,Test)
  pvalues <- switch(alternative, 
                    "two.sided" =  matrixStats::colRanks(-abs(TT)) / (nrow(TT)),
                    "greater" = matrixStats::colRanks(-TT) / (nrow(TT)),
                    "less" = matrixStats::colRanks(TT) / (nrow(TT)))

  clstr_id <- eval(parse(text=paste0("model$multiple_comparison$", eff, "$clustermass$data$cluster_id")))
  
  clusters <- c(1:eval(parse(text=paste0("model$multiple_comparison$", eff, "$clustermass$cluster$no"))))
  hom <-hommel(pvalues[,1])
  
  out=lapply(clusters,function(i){
    ix= which(clstr_id == i)
    summary_hommel_eeg(hommel = hom,ix=ix, alpha = alpha, clusters = i,eff = eff)
    
    
  })
  
  out_d <- t(as.data.frame(out))
  rownames(out_d) = NULL
  out_d
}


