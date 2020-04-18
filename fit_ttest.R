#' Perform t test 
#' 
#' @param features A data frame of features. 
#' Must have the exact same columns (samples) as the rows in metadata. 
#' @param metadata A data frame of metadata to be associated.
#' Must have the exact same rows (samples) as the columns in features.
#' @param MultTestCorrection Multiple testing correction method. Default is 'BH'.
#' @export
#' 
fit.ttest <- function(countdata,
                      metadata,
                      MultTestCorrection = 'BH') {
  name_metadata <- names(metadata)
    idx <- seq_len(nrow(countdata))
    ttest_p <- sapply(idx, function(i) {
      t.test(as.matrix(countdata)[i, ] ~ metadata[,name_metadata])$p.value
    })

    paras <- data.frame(
      coef = NA,
      pval = ttest_p,
      metadata = name_metadata,
      feature = row.names(countdata)
    )
    row.names(paras) <- NULL
    
    paras$qval<-as.numeric(p.adjust(paras$pval, method = MultTestCorrection))
    paras<-paras[order(paras$qval, decreasing=FALSE),]
    paras<-dplyr::select(paras, c('feature', 'metadata'), everything())
    rownames(paras)<-NULL
    return(paras)   
}