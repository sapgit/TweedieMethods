#' Fit a ZICP model to the feature and metadata tables.
#' 
#' @param features A data frame of features. 
#' Must have the exact same rows (samples) as metadata. 
#' @param metadata A data frame of metadata to be associated.
#' Must have the exact same rows (samples) as features.
#' @param libSize A vector of library size across samples to be included as an offset in the model.
#' @param MultTestCorrection Multiple testing correction method. Default is 'BH'.
#' @export
#' 
fit.ZICP <- function(features, 
                     metadata, 
                     libSize,
                     MultTestCorrection = "BH"){
  paras <- pbapply::pbsapply(1:ncol(features), simplify=FALSE, function(x){
    
    featuresVector <- features[, x] 
    
    # Fit Model
    dat_sub <- data.frame(expr = as.numeric(featuresVector), metadata, libSize)
    formula<-as.formula(paste("expr ~ ", paste(colnames(metadata), collapse= "+")))
    
    # Library size adjustment
    formula<-update(formula, . ~ . - offset(log(libSize)))
    # formula<-update(formula, . ~ log(libSize))
    
    # Fit model
    fit <- tryCatch({
      fit1 <- cplm::zcpglm(formula, data = dat_sub)
    }, error=function(err){
      fit1 <- try({cplm::zcpglm(formula, data = dat_sub)}) 
      return(fit1)
    })
    
    if (class(fit) != "try-error"){
      sink("/dev/null")
      para<-as.data.frame(summary(fit)$coefficients$tweedie)[-1,-c(2:3)]
      sink()
      colnames(para)<-c('coef', 'pval')
      para$metadata<-colnames(metadata)
      para$feature<-colnames(features)[x]
      rownames(para)<-NULL
    }
    else{
      print(paste("Fitting problem for feature", x, "returning NA"))
      para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=2))
      colnames(para)<-c('coef', 'pval')
      para$metadata<-colnames(metadata)
      para$feature<-colnames(features)[x]
      rownames(para)<-NULL
    }
    return(para)
  })

  paras<-do.call(rbind, paras)
  paras$qval<-as.numeric(p.adjust(paras$pval, method = MultTestCorrection))
  paras<-paras[order(paras$qval, decreasing=FALSE),]
  paras<-dplyr::select(paras, c('feature', 'metadata'), everything())
  rownames(paras)<-NULL
  return(paras)   
}
