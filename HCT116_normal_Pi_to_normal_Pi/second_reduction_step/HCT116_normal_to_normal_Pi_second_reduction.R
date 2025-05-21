#' 2025-05-15
#' Script to analyse the data for the condition HCT116 normal to normal Pi in the 
#' the article https://doi.org/10.1101/2025.03.27.645712
#' Second reduction model
#' author: Jacques Hermes <jacques.hermes@fdm.uni-freiburg.de>

# Sets current working directory to the location of this script
try(setwd(dirname(rstudioapi::getSourceEditorContext()$path)))

## Loading Libraries --------------------
library(deSolve)
library(trust)
library(parallel)
library(ggplot2)
library(ggthemes)
library(cOde)
library(dMod)
library(cowplot)
library(magrittr)
library(tidyr)
library(fBasics)
library(data.table)
library(purrr)

## Model Definition - Equations --------------------

tolerances <- 1e-10
set.seed(42)

model_name <- "HCT116_normal_to_normal_Pi"
flist_O <- NULL %>% 

  addReaction("InsP6", "InsP75", "k2*ratio_0O*InsP6") %>%
  addReaction("InsP6","InsP75_1O", "k2*ratio_1O*InsP6") %>%
  addReaction("InsP6","InsP75_2O", "k2*ratio_2O*InsP6") %>%
  addReaction("InsP6","InsP75_3O", "k2*ratio_3O*InsP6") %>%

  addReaction("InsP6", "InsP71", "k3*ratio_0O") %>%
  addReaction("InsP6","InsP71_1O", "k3*ratio_1O") %>%
  addReaction("InsP6","InsP71_2O", "k3*ratio_2O") %>%
  addReaction("InsP6","InsP71_3O", "k3*ratio_3O") %>%

  addReaction( "InsP75","InsP6", "k4*InsP75") %>%
  addReaction("InsP75_1O","InsP6", "k4*InsP75_1O") %>%
  addReaction("InsP75_2O","InsP6", "k4*InsP75_2O") %>%
  addReaction("InsP75_3O","InsP6", "k4*InsP75_3O") %>%

  addReaction("InsP71"   ,"InsP6", "k5*InsP71") %>%
  addReaction("InsP71_1O","InsP6", "k5*InsP71_1O") %>%
  addReaction("InsP71_2O","InsP6", "k5*InsP71_2O") %>%
  addReaction("InsP71_3O","InsP6", "k5*InsP71_3O") %>%

  addReaction("InsP71","InsP8"          , "k6*ratio_0O") %>%
  addReaction("InsP71","InsP8_1O_10O51O", "k6*ratio_1O") %>%
  addReaction("InsP71","InsP8_2O_10O52O", "k6*ratio_2O") %>%
  addReaction("InsP71","InsP8_3O_10O53O", "k6*ratio_3O") %>%
  
  addReaction("InsP71_1O","InsP8_1O_11O50O", "k6*ratio_0O") %>%
  addReaction("InsP71_1O","InsP8_2O_11O51O", "k6*ratio_1O") %>%
  addReaction("InsP71_1O","InsP8_3O_11O52O", "k6*ratio_2O") %>%
  addReaction("InsP71_1O","InsP8_4O_11O53O", "k6*ratio_3O") %>%

  addReaction("InsP71_2O","InsP8_2O_12O50O", "k6*ratio_0O") %>%
  addReaction("InsP71_2O","InsP8_3O_12O51O", "k6*ratio_1O") %>%
  addReaction("InsP71_2O","InsP8_4O_12O52O", "k6*ratio_2O") %>%
  addReaction("InsP71_2O","InsP8_5O_12O53O", "k6*ratio_3O") %>%
  
  addReaction("InsP71_3O","InsP8_3O_13O50O", "k6*ratio_0O") %>%
  addReaction("InsP71_3O","InsP8_4O_13O51O", "k6*ratio_1O") %>%
  addReaction("InsP71_3O","InsP8_5O_13O52O", "k6*ratio_2O") %>%
  addReaction("InsP71_3O","InsP8_6O"       , "k6*ratio_3O") %>%

  addReaction("InsP8"          ,"InsP71"   , "k8*InsP8") %>%
  addReaction("InsP8_1O_10O51O","InsP71"   , "k8*InsP8_1O_10O51O") %>%
  addReaction("InsP8_2O_10O52O","InsP71"   , "k8*InsP8_2O_10O52O") %>%
  addReaction("InsP8_3O_10O53O","InsP71"   , "k8*InsP8_3O_10O53O") %>%
  
  addReaction("InsP8_1O_11O50O","InsP71_1O", "k8*InsP8_1O_11O50O") %>%
  addReaction("InsP8_2O_11O51O","InsP71_1O", "k8*InsP8_2O_11O51O") %>%
  addReaction("InsP8_3O_11O52O","InsP71_1O", "k8*InsP8_3O_11O52O") %>%
  addReaction("InsP8_4O_11O53O","InsP71_1O", "k8*InsP8_4O_11O53O") %>%
  
  addReaction("InsP8_2O_12O50O","InsP71_2O", "k8*InsP8_2O_12O50O") %>%
  addReaction("InsP8_3O_12O51O","InsP71_2O", "k8*InsP8_3O_12O51O") %>%
  addReaction("InsP8_4O_12O52O","InsP71_2O", "k8*InsP8_4O_12O52O") %>%
  addReaction("InsP8_5O_12O53O","InsP71_2O", "k8*InsP8_5O_12O53O") %>%
  
  addReaction("InsP8_3O_13O50O","InsP71_3O", "k8*InsP8_3O_13O50O") %>%
  addReaction("InsP8_4O_13O51O","InsP71_3O", "k8*InsP8_4O_13O51O") %>%
  addReaction("InsP8_5O_13O52O","InsP71_3O", "k8*InsP8_5O_13O52O") %>%
  addReaction("InsP8_6O"       ,"InsP71_3O", "k8*InsP8_6O") %>%
  
  addReaction("InsP8"          ,"InsP75"   , "k9*InsP8") %>%
  addReaction("InsP8_1O_11O50O","InsP75"   , "k9*InsP8_1O_11O50O") %>%
  addReaction("InsP8_2O_12O50O","InsP75"   , "k9*InsP8_2O_12O50O") %>%
  addReaction("InsP8_3O_13O50O","InsP75"   , "k9*InsP8_3O_13O50O") %>%
  
  addReaction("InsP8_1O_10O51O","InsP75_1O", "k9*InsP8_1O_10O51O") %>%
  addReaction("InsP8_2O_11O51O","InsP75_1O", "k9*InsP8_2O_11O51O") %>%
  addReaction("InsP8_3O_12O51O","InsP75_1O", "k9*InsP8_3O_12O51O") %>%
  addReaction("InsP8_4O_13O51O","InsP75_1O", "k9*InsP8_4O_13O51O") %>%
  
  addReaction("InsP8_2O_10O52O","InsP75_2O", "k9*InsP8_2O_10O52O") %>%
  addReaction("InsP8_3O_11O52O","InsP75_2O", "k9*InsP8_3O_11O52O") %>%
  addReaction("InsP8_4O_12O52O","InsP75_2O", "k9*InsP8_4O_12O52O") %>%
  addReaction("InsP8_5O_13O52O","InsP75_2O", "k9*InsP8_5O_13O52O") %>%
  
  addReaction("InsP8_3O_10O53O","InsP75_3O", "k9*InsP8_3O_10O53O") %>%
  addReaction("InsP8_4O_11O53O","InsP75_3O", "k9*InsP8_4O_11O53O") %>%
  addReaction("InsP8_5O_12O53O","InsP75_3O", "k9*InsP8_5O_12O53O") %>%
  addReaction("InsP8_6O"       ,"InsP75_3O", "k9*InsP8_6O") %>%
  
  addReaction("ratio_0O","ratio_1O", "phospho0*ratio_0O") %>%
  addReaction("ratio_0O","ratio_2O", "phospho1*ratio_0O") %>%
  addReaction("ratio_0O","ratio_3O", "phospho2*ratio_0O") %>%
  addReaction("ratio_1O","ratio_0O", "phospho_0*ratio_1O") %>%
  addReaction("ratio_2O","ratio_0O", "phospho_1*ratio_2O") %>%
  addReaction("ratio_3O","ratio_0O", "phospho_2*ratio_3O ") %>%
  

  addReaction("InsP6","", "deg_InsP6*InsP6") %>%
  addReaction("","InsP6", "prod_InsP6")

## Model Definition - ODE model --------------------

modelIPD_O <- odemodel(flist_O, modelname = paste0("odemodel_", model_name),
                       jacobian = "inz.lsodes", compile = TRUE)

## Model Definition - Observables --------------------

observables <-  eqnvec( 
  
  Obs_InsP6 = "InsP6",
  
  Obs_ratio_0O = "scale0*ratio_0O",

  Obs_ratio_1O = "scale1*ratio_1O",
  Obs_ratio_2O = "scale2*ratio_2O",
  Obs_ratio_3O = "scale3*ratio_3O",
  
  
  
  Obs_InsP75 = " InsP75",
  Obs_InsP75_1O = "InsP75_1O",
  Obs_InsP75_2O = "InsP75_2O",
  
  Obs_InsP71 = " InsP71",
  Obs_InsP71_1O = "InsP71_1O",
  Obs_InsP71_2O = "InsP71_2O",
  
  Obs_InsP8 = "InsP8",
  Obs_InsP8_1O = "InsP8_1O_10O51O + InsP8_1O_11O50O ",
  Obs_InsP8_2O = "InsP8_2O_10O52O + InsP8_2O_12O50O + InsP8_2O_11O51O"
)

## Model Definition - Error model --------------------

obsError <- c(names(observables))
errors <- as.eqnvec(
  c(paste0("sigma_",obsError,"_abs")
  ), names = c(obsError)
)

errPars <- getSymbols(errors)
errPars <- errPars[which(grepl("sigma", errPars))]

## Data processing --------------------

data <- as.datalist(wide2long(read.csv(
  file ="Data_Labeling_HCT116_normal_to_normal_Pi.csv",
  sep = ","),keep = c(1,16)))

data <- subset(data, !is.na(value))

data_w <- data[1]

covtable <- data.frame(
  condition = names(data_w),
  row.names =names(data_w))

## Model Definition - Parameter transformations --------------------

constraints <- resolveRecurrence(c(
  
  InsP75_2O = "0",
  InsP75_3O = "0",
  
  InsP71_2O = "0",
  InsP71_3O = "0",
  
  InsP8_1O_11O50O="0",
  
  InsP8_2O_12O50O="0",
  InsP8_2O_10O52O="0",
  InsP8_2O_11O51O="0",
  
  InsP8_3O_13O50O="0",
  InsP8_3O_10O53O="0",
  InsP8_3O_12O51O="0",
  InsP8_3O_11O52O="0",
  
  InsP8_4O_13O51O="0",
  InsP8_4O_11O53O="0",
  InsP8_4O_12O52O="0",
  
  InsP8_5O_13O52O="0",
  InsP8_5O_12O53O="0",
  
  scale0="1",
  
  InsP8_6O="0",
  prod_InsP6 = "0"
  
))

innerpars <- unique(c(getParameters(modelIPD_O), getSymbols(observables), errPars)) #parameters seen by model
names(innerpars) <- innerpars

trafo <- replaceSymbols(names(observables), observables, innerpars)
trafo <- replaceSymbols(names(constraints), constraints, trafo)
trafoL <- branch(trafo, table=covtable) %>%
  insert(x~y,x="ratio_0O", y="(1-ratio_1O-ratio_2O-ratio_3O)") %>% # Ensures sum of ratios always = 1
  
  insert(x~"0.01", x = c("ratio_1O")) %>%
  
  insert(x~"0", x = c("ratio_2O")) %>%
  insert(x~"0", x = c("ratio_3O")) %>%
  
  insert(x~y, x = c("k5"), y= c("rat_k5_k3 * k3 "), conditionMatch = "HCT116") %>%
  
  insert(x~(x_H), x=getParameters(flist_O)[which(!getParameters(flist_O) %in% names(constraints))], conditionMatch = "HCT116") %>%
  insert(x~(x_H), x=errPars, conditionMatch = "HCT116") %>%
  
  insert("x~10^(x)", x = .currentSymbols)

trafoL <- trafoL[names(data_w)]

## Observation and Prediction functions --------------------

p0 <- x_o <- NULL
for (C in names(trafoL)) {
  p0 <- p0 + P(trafoL[[C]], condition = C)
  x_o <- x_o + Xs(modelIPD_O, optionsOde = list(method = "lsoda", rtol = tolerances, atol = tolerances, maxsteps = 7000),
                  optionsSens = list(method = "lsodes", lrw=200000, rtol = tolerances, atol = tolerances),
                  condition = C)
}

g <- Y(observables, f=x_o, compile=TRUE, modelname="Obsfn")

e <- Y(errors, observables, attach.input=F, compile=TRUE, modelname="Err")

outerpars <- getParameters(p0) 
pouter <-  structure(rep(-1,length(outerpars)), names = outerpars)

## Objective function --------------------

obj_o <- normL2(data_w, g * x_o * p0,errmodel = e) + constraintL2(pouter, sigma = 5)
times <- sort(unique(c(seq(0**0.5, max(21)**0.5, len=100)**2,do.call(rbind, data)$time)))

## Multistart fit --------------------

out_mstrust_o <- mstrust(obj_o, pouter, studyname = "HCT116_normal_to_normal_pi_second_reduction_fits",rinit = 1, rmax = 10, iterlim = 2000,
                         sd = 4, cores = 14, fits = 200, blather=TRUE)

print("Fits Model O done")

## Fit processing and plotting --------------------

myframe_o <- as.parframe(out_mstrust_o)
save(myframe_o, file = "Results_HCT116_normal_to_normal_pi_second_reduction.RData")

p1 <-plotValues(myframe_o, tol = 0.01, value<=2000)
ggsave("Waterfall_HCT116_normal_to_normal_pi_second_reduction.pdf",p1, height = 8,width = 10)

bestfit_o <- as.parvec(myframe_o)
save(bestfit_o, file = "bestfit_HCT116_normal_to_normal_pi_second_reduction.RData")

prediction_o <- (g * x_o * p0)(seq(1,130,1), bestfit_o, fixed=NULL, deriv=FALSE)
pred_o <- as.data.frame(prediction_o, data = data, errfn = e)

p2 <-plotCombined(prediction_o, data_w)+
  geom_ribbon(data = pred_o,aes(ymin = value - sigma, ymax = value + sigma), alpha = 0.25, color=NA)+ theme(legend.position = "bottom")

ggsave("Bestfit_HCT116_normal_to_normal_pi_second_reduction.pdf",p2, height = 8,width = 10)

pred_only_fitted <- subset(pred_o, name %in% c(as.character(unique(data_w$HCT116_18O_Water$name))))

p3<-plotCombined(prediction_o, data_w, name %in%  c(as.character(unique(data_w$HCT116_18O_Water$name))))+
  geom_ribbon(data = pred_only_fitted,aes(ymin = value - sigma, ymax = value + sigma), alpha = 0.25, color=NA)+ theme(legend.position = "bottom")

ggsave("Bestfit_HCT116_normal_to_normal_pi_second_reduction_only_fitted.pdf",p3, height = 8,width = 10)

BIC<- length(bestfit_o)*log(length(data_w$HCT116_18O_Water$value))+myframe_o$value[1]

## Profile Likelihood calculation --------------------

profiles_o <- profile(
  obj = obj_o,
  pars =  bestfit_o,
  whichPar = names(bestfit_o),
  cores = 12,
  alpha=0.01,
  limits=c(-7, 7),
  method = "optimize"
)

save(profiles_o, file = "Profiles_HCT116_normal_to_normal_pi_second_reduction.RData")
load(file = "Profiles_HCT116_normal_to_normal_pi_second_reduction.RData")
print("Profiles done")

## Plotting profiles --------------------

p8<-plotProfile(profiles_o, mode == "data")

ggsave("Profiles_HCT116_normal_to_normal_pi_second_reduction.pdf",p8, height = 10,width = 12)

