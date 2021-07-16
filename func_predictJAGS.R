
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
##*** Function to predict relationship with covariates
##*** in linear models fitted in JAGS using jagsUI
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
#* Ismael V. Brack
#* Jun/2021

# fm = fitted model (jagsUI output)
# params = parameters of the linear equation. Names must be equal to fm object
# newdata = data.frame with covariate values to predict
# link = link functions for the response allowed: "identity", "log" and "logit"

predictJAGS <- function(fm, params, newdata, quants=c(0.025,0.5,0.975),link="identity",appendData=TRUE){

  if(is(out1) != "jagsUI"){
    stop("fm is not a 'jagsUI' object")
  }
  if(is.null(fm$sims.list)){
    stop("Provide a fm object containing a 'sims.list'")
  }
  if(any((params %in% names(fm$sims.list))==FALSE)){
    stop("Names in params must be equal to fm object")
  }
  if(ncol(newdata) != (length(params)-1)){
    stop("ncol(newdata) must be equal to the number of params - 1")
  }
  
  # get MCMC posterior values for parameters
  pars.post <- list()
  for(p in 1:length(params)){
    pars.post[[p]] <- as.vector(fm$sims.list[[params[p]]])
  }
  pars.post <- do.call(cbind,pars.post)
  
  # Include an intercept to newdata
  newdata.m <- as.matrix(cbind(1,newdata))
  
  # matrix to receive predictions
  fit.post <- matrix(NA, nrow=nrow(newdata.m), ncol=nrow(pars.post))
  
  # predict...
  for(i in 1:nrow(pars.post)){ # each line correspond to a covariate value
    fit.post[,i] <- newdata.m %*% pars.post[i,] # each column correpond to sample of the 
  }                                               # posteriors distributions (MCMC values)
  
  if(link=="log"){
    fit.post <- exp(fit.post)
  }
  if(link=="logit"){
    fit.post <- plogis(fit.post)
  }
  
  # arrange table with summary results for prediction
  preds <- as.data.frame(t(apply(fit.post,1,
                                 function(x) c(mean(x),quantile(x,probs=quants)) ))
  )
  colnames(preds)[1] <- "mean"
  
  if(appendData==TRUE){
    preds <- cbind(preds, newdata)
  }
  
  return(preds)
}


