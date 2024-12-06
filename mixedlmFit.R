
mixedlmFit = function(design.x, data, contrast, keep){
  library(lme4)
  library(multcomp)
  
  library(parallel)
  # An mc-version of the sapply function.
  mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
    FUN <- match.fun(FUN)
    answer <- parallel::mclapply(X = X, FUN = FUN, ...)
    if (USE.NAMES && is.character(X) && is.null(names(answer))) 
      names(answer) <- X
    if (!isFALSE(simplify) && length(answer)) 
      simplify2array(answer, higher = (simplify == "array"))
    else answer
  }
  
  do_one_gene_analysis = function(gene_name){
    geneidx=which(keep==gene_name)
    cat(gene_name, ": Gene ",geneidx,"/",length(keep),round(geneidx*100/length(keep), digits=3)," % \n")
    data$gene = log(data[,gene_name]+1)
    
    
    form = as.formula(paste( "gene ~", design.x))
    
    fit11 = lmer(form, data=data,
                 REML = F,
                 control = lmerControl(calc.derivs=F))
    # for(contrast in contrasts){
    smr = summary(glht(fit11, contrast))
    #smr
    one.result= c(AIC(fit11), gene_name,smr$test$coefficients, smr$test$sigma, smr$test$tstat,smr$test$pvalues, contrast)
    names(one.result) = c("AIC","gene","coefficients","sigma","tstat","pvalues","contrast")
    # res[res[,"gene"]==gene_name & res[,"contrast"]==contrast,names(one.result) ] <- one.result
    #}
    one.result
  }
  mods = mcsapply(keep,FUN=do_one_gene_analysis)
  
  mods=as.data.frame(t(mods))
  num_vars=setdiff(colnames(mods),c("gene","contrast"))
  mods[,num_vars] = apply(mods[,num_vars],2, function(x){as.numeric(x)})
  
  mods=mods[order(abs(mods$coefficients), decreasing=T),]
  
  mods$pval_adj = p.adjust(round(mods$pvalues, digits=10),"BH")
  mods=mods[order(mods$pval_adj),]
  mods$sig = (mods$pval_adj<.01)
  mods$sig5 = (mods$pval_adj<.05)
  mods$sigeff = (mods$pval_adj<.05)*(abs(mods$coefficients)>=0.5)
  mods
}