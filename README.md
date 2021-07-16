# predictJAGS_func
Function to predict relationships with covariates from outputs of Bayesian linear models fitted with jagsUI
\.
It was written to predict "regression-type" relationships, so far.. but it can be easily used for "ANCOVA-type" linear models by separately predicting relations for each factor level.

#### Arguments
  - fm = fitted model (jagsUI output)
  - params = parameters of the linear equation. Names must be equal to fm object
  - newdata = data.frame with covariate values to predict
  - link = link functions for the response allowed: "identity", "log" and "logit"
