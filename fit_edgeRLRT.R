#' Fit a EdgeRLRT model to the feature and metadata tables.
#' 
#' @param features A data frame of features. 
#' Must have the exact same columns (samples) as the rows in metadata. 
#' @param metadata A data frame of metadata to be associated.
#' Must have the exact same rows (samples) as the columns in features.
#' @param MultTestCorrection Multiple testing correction method. Default is 'BH'.
#' @export
#' 
fit.edgeRLRT <- function(countdata,
                         metadata,
                         MultTestCorrection = 'BH') {
  name_metadata <- names(metadata)
    dge <- DGEList(countdata, group = metadata[,name_metadata])
    dge <- calcNormFactors(dge)
    design <- model.matrix(~metadata[,name_metadata])
    dge <- estimateDisp(dge, design = design)
    fit <- glmFit(dge, design = design)
    lrt <- glmLRT(fit)
    
    paras <- data.frame(
      coef = lrt$coefficients[,-1],
      pval = lrt$table$PValue,
      metadata = name_metadata,
      feature = row.names(lrt$coefficients)
    )
    row.names(paras) <- NULL
    
    paras$qval<-as.numeric(p.adjust(paras$pval, method = MultTestCorrection))
    paras<-paras[order(paras$qval, decreasing=FALSE),]
    paras<-dplyr::select(paras, c('feature', 'metadata'), everything())
    rownames(paras)<-NULL
    return(paras)   
}
