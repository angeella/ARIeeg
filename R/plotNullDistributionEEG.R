#' @title plot null distribution EEG 
#' @description plot null distribution EEG
#' @usage plotNullDistributionEEG(P,family,alpha, ct, path, name, delta)
#' @param P pvalues matrix permutation
#' @param alpha alpha level
#' @param family which family for the confidence envelope? simes, finner, beta or higher.criticism. default is simes
#' @param delta do you want to consider at least delta size set?
#' @param ct set of thresholds
#' @param path path
#' @param name plot name
#' @author Angela Andreella
#' @return Returns plot null distribution
#' @export
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom grDevices rainbow
#' @importFrom graphics legend
#' 


plotNullDistributionEEG <- function(P,family="simes",alpha = 0.1, ct = c(0,1), path = getwd(), name = "plot", delta = NULL){
  
  family_set <- c("simes", "finner", "beta", "higher.criticism")
  fam_match <- function(x) {match.arg(tolower(x), family_set)}
  if(!is.null(family)){family <- unlist(lapply(family, fam_match))}
  if(is.null(P)){stop('Please insert pvalues matrix')}
  
  if(!is.null(P) & is.unsorted(P[1,])){pvalues_ord <- rowSortC(P)}
  
  if(is.null(family)){
    png(paste0(path,"/", name, ".png")) 
    plot(pvalues_ord[1,], type = 'l', col = ' red', xlab = expression(i), ylab = expression(p[(i)]),ylim=c(0,1))
    for(i in 2:nrow(pvalues_ord)){
      
      lines(pvalues_ord[i,],col='black',type="l")
      
    }
    lines(pvalues_ord[1,], lwd =2, col= 'red')
    dev.off()
  }else{
    lcv <- function(family,delta=NULL, cols = "blue"){
      lambdaO <- lambdaOpt(pvalues = pvalues_ord,family=family,ct=ct,alpha=alpha, delta = delta)
      cvO<- cv(pvalues = pvalues_ord, family = family, alpha = alpha, lambda = lambdaO, delta = delta)
      if(is.unsorted(cvO)){
        idS = which(sapply(c(1:length(cvO)), function(x) is.unsorted(cvO[1:x])))[1]
        cvO = c(cvO[1:(idS-1)], rep(length(cvO)-1, max(cvO[1:(idS-1)])))  
        
      }
      lines(cvO, lwd =2, col= cols)
      
    }
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
    lcvV <- Vectorize(lcv,vectorize.args = c("family", "delta", "cols"))
    cols = rainbow(length(family))
    cols[3] = "#FF0099FF"
    png(paste0(path,"/", name, ".png")) 
    plot(pvalues_ord[1,], type = 'l', col = ' green', xlab = expression(i), ylab = expression(p[(i)]))
    for(i in 2:nrow(pvalues_ord)){
      
      lines(pvalues_ord[i,],col='black',type="l")
      
    }
    lines(pvalues_ord[1,], lwd =2, col= 'green')
    #lines(cvO, col= 'blue', lwd =2)
    mapply(lcv, family, delta, cols)
    family <- firstup(family)
    fn = sum(family!="beta" | family!="higher.criticism")
    legend('top',legend=c(sapply(c(1:fn), 
                                 function(x) if(family[x]=="Beta" | family[x]=="Higher.criticism") {paste0(" ",family[x])}else{as.expression(bquote(~ .(family[x]) ~ delta == .(delta[x]) ))} ), 
                          " Observed Pvalues"), col= c(cols, "green"),lwd =2)
    
    dev.off()
  }
  
  
}