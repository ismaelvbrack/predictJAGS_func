# predictJAGS_func
Function to predict relationships with covariates from outputs of Bayesian linear models fitted with jagsUI
\.
It was built to predict "regression-type" and "ANCOVA-type" relationships, so far.. but it can be adapted to other applications...

#### Arguments
  - fm = fitted model (jagsUI output)
  - params = parameters of the linear equation. Names must be equal to fm object
  - newdata = data.frame with covariate values to predict
  - link = link functions for the response allowed: "identity", "log" and "logit"
  - appendData = logical. Include covariates values in the output?

\.
See the example with Swiss Jays from the AHM1 book

