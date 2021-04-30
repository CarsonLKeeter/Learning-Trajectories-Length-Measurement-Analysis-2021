#######################################
## Run "Model_construction.R" First! ##
#######################################

## Traceplot for Chain Convergence 
model_trace <- summary_Model$Traceplot 

## Pareto K/PSIS Diagnostic Plot 
plot(loo_ex)

####################################################################################
