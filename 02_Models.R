



# Created 13.10.2021



# This code runs linear models for asthma. With correct specification of the parameters
# adj, lag, dec and groupsid, the code examined adjusted/unadjusted models, different lags (lag5 is the max)
# percentiles for flexible temperature fits and different age and sex groups.



#########################################################################


library(INLA)
library(dplyr)
library(tidyr)
library(sf)
library(spdep)

options(scipen = 999)


main_analysis <- TRUE
assess_RW <- TRUE
spaceEM <- FALSE
timeEM <- FALSE


# data working directory
wd_data <- "data/"
# results to store working directory
wd_res <- "results/"


asthma_dat <- readRDS(paste0(wd_data, "dat"))

head(asthma_dat)

asthma_dat$sex <- ifelse(asthma_dat$sex == "male", 1, 2)

# names to store
names2store <- expand.grid(age = c("0_4", "5_14", "15-65", "65plus", "Total"), sex = c(1,2,3))
GROUPS2RUN <- expand.grid(age = c("0-4", "5-14", "15-65", "65+", "Total"), sex = c(1,2,3))

# global INLA options
# priors
hyper.iid <- list(theta = list(prior="pc.prec", param=c(0.1, 0.01)))
hyper.bym <- list(theta1 = list("PCprior", c(0.1, 0.01)), theta2 = list("PCprior", c(0.5, 0.5)))

# default settings
control.family=inla.set.control.family.default()


RunINLACA <- function(groupsid, adj, metric, rw2 = FALSE, thet = NULL, rerun = TRUE, store = "all", year_start = 2002, year_end = 2019, efmod = "no", region = NULL){
  
  #--------------------------------------------------
  # groupsid                     as per GROUPS2RUN
  # adj                          if adj FALSE model is unadjusted for covariates, if adj != FALSE then you need to input the name of the covariates to be adjusted
  # rw2                          If rw2 TRUE temperature is included as a Gaussian process
  # thet                         the mode of the hyperparameters for INLA
  # rerun                        if true the model reruns with better initial values
  # store                        "all" or a vector of which INLA elements to be stored
  # year_start                   the first year of the analysis
  # year_end                     the last year of the analysis
  # efmod                        takes the values: no, space, time and spacetime (includes separable space time effects for the temperature in a sinlge model) and examines space and time efect modification
  #--------------------------------------------------

  asthma_dat %>% filter(year >= year_start, year <= year_end) -> asthma_dat
  asthma_dat$ID_year <- (as.numeric(asthma_dat$year) - year_start) + 1

  if(!is.null(region)){
  asthma_dat %>% filter(regid %in% region) -> asthma_dat
}
  
  if(groupsid <= 10){
    if(GROUPS2RUN[,1][groupsid] %in% "Total"){
      asthma_dat %>% filter(sex %in% GROUPS2RUN[,2][groupsid]) ->
        asthma_dat
    }else{
      asthma_dat %>% filter(age.group %in% GROUPS2RUN[,1][groupsid]) %>% filter(sex %in% GROUPS2RUN[,2][groupsid]) ->
        asthma_dat
    }
  }
  
  if((groupsid >= 11) & (groupsid < 15)){
    asthma_dat %>% filter(age.group %in% GROUPS2RUN[,1][groupsid]) -> asthma_dat
  }
  
  if(groupsid == 15){
    asthma_dat -> asthma_dat
  }
  
 
  # fix the id
  asthma_dat$ID_2 <- as.numeric(as.factor(asthma_dat$ID))
  asthma_dat$ID.year <- as.numeric(asthma_dat$year) - 2001

  asthma_dat$tmp.ID <- inla.group(asthma_dat$temperature, n = 80, method = "cut", idx.only = FALSE)

  if(efmod == "no"){
    efmodnam = ""
    if(adj[1] == FALSE){
      if(rw2){
        formula = 
          O ~ 1 + 
          f(tmp.ID, model='rw2', hyper=hyper.iid, constr = TRUE, scale.model = TRUE) +
          f(ID_2, model='iid', hyper = list(prec = list(prior = "gamma", initial = 1/100, fixed = TRUE))) 
        nam.rw2 <- "_rw2"
      }else{
        formula = 
          O ~ 1 + 
          temperature +
          f(ID_2, model='iid', hyper = list(prec = list(prior = "gamma", initial = 1/100, fixed = TRUE)))
        nam.rw2 <- "" 
      }
      mod.cov <- "unadjusted"
     }else{
      if(rw2){
        formula = paste0("O~", paste0(c(1, "hol", adj), collapse  = "+"), 
                         "+f(tmp.ID, model='rw2', hyper=hyper.iid, constr = TRUE, scale.model = TRUE)+f(ID_2, model='iid', hyper = list(prec = list(prior = 'gamma', initial = 1/100, fixed = TRUE)))+f(ID.rec, model='iid', hyper=hyper.iid)")
        nam.rw2 <- "_rw2"
      }else{
        formula = paste0("O~", paste0(c(1, "hol", "temperature", adj), collapse  = "+"), 
                         "+f(ID_2, model='iid', hyper = list(prec = list(prior = 'gamma', initial = 1/100, fixed = TRUE)))+f(ID.rec, model='iid', hyper=hyper.iid)")
        nam.rw2 <- ""
      }
      mod.cov <- paste0(c("hol", "recur", adj), collapse  = "_")
    }
  }
  
  
  ##
  ## separate models for space and time interactions with temperature:

   if(efmod == "time"){ 
      nam.rw2 <- ""
      efmodnam = "_emtime"
      mod.cov <- paste0(c("hol", "recur", adj), collapse  = "_") 
      formula = paste0("O~", paste0(c(1, "temperature", "hol", adj), collapse  = "+"), 
                "+f(ID_2, model='iid', hyper = list(prec = list(prior = 'gamma', initial = 1/100, fixed = TRUE)))+f(ID.rec, model='iid', hyper=hyper.iid)")
  }
    

   if(efmod == "space"){ 
     mod.cov <- paste0(c("hol", "recur", adj), collapse  = "_") 

     nam.rw2 <- ""
     formula = paste0("O~", paste0(c(1, "temperature", "hol", adj), collapse  = "+"), 
                       "+f(ID_2, model='iid', hyper = list(prec = list(prior = 'gamma', initial = 1/100, fixed = TRUE)))+f(ID.rec, model='iid', hyper=hyper.iid)")
     efmodnam = "_emspace_"

  }

 
  # INLA call
  res = inla(as.formula(formula),
             data=asthma_dat,
             family="Poisson",
             
             # priors for the fixed effects
             control.fixed = list(
               
               
               mean.intercept = 0, 
               prec.intercept = 1,
               mean = list(
                 hol = 0, rh = 0, temperature = 0, ws = 0, rainfall = 1
               ), 
               
               prec = list(
                 hol = 1, rh = 1, temperature = 1, ws = 1, rainfall = 1
               )
             ),
             
             control.family=control.family,
             verbose = TRUE,
             num.threads = 10,
             control.compute=list(config = TRUE, dic = TRUE, waic = TRUE), 
             control.predictor=list(link = 1), 
             control.mode=list(theta = thet, restart = TRUE))

  
  if(rerun){
  # sometime it crashes on the rerun, so make sure you store it before
  if(store[1] == "all"){

    saveRDS(res, file = paste0(wd_res, "AsthmaLinear_", names2store[,1][groupsid],"_", 
                               names2store[,2][groupsid], "_", mod.cov, nam.rw2,  
                               "_", year_start, efmodnam, region))
  }else{
    saveRDS(res[store], file = paste0(wd_res, "AsthmaLinear_", names2store[,1][groupsid],"_", 
                                      names2store[,2][groupsid], "_", mod.cov, nam.rw2,  
                                      "_", year_start, efmodnam, region))
  }  
    res <- inla.rerun(res)
  }
  
  
  if(store[1] == "all"){
   saveRDS(res, file = paste0(wd_res, "AsthmaLinear_", names2store[,1][groupsid],"_", 
                               names2store[,2][groupsid], "_", mod.cov, nam.rw2,  
                               "_", year_start, efmodnam, region))
  }else{

    saveRDS(res[store], file = paste0(wd_res, "AsthmaLinear_", names2store[,1][groupsid],"_", 
                                      names2store[,2][groupsid], "_", mod.cov, nam.rw2,  
                                      "_", year_start, efmodnam, region))
  }  
}


t_0_all <- Sys.time() 


# 1. Main analysis

list_loop <- list(
  list(adj = FALSE, year_start = 2002, year_end = 2019),
  list(adj = c("rainfall", "ws", "rh"), year_start = 2007, year_end = 2019)
)

t_0 <- Sys.time()
if(main_analysis == TRUE){
  print("main_analysis")
  for(i in 1:length(list_loop)){
    print(paste0("i=", i))
    for(k in 1:15){
      print(paste0("k=", k))
      RunINLACA(groupsid = k, 
        adj = list_loop[[i]]$adj, 
        rw2 = FALSE, 
        thet = NULL, 
        rerun = FALSE, 
        store = c("summary.fixed")
        )
      gc()
    }
  }
}
t_1 <- Sys.time()
t_1 - t_0


# 2. Linearity assumption

t_0 <- Sys.time()
if(assess_RW == TRUE){
  print("assess_RW")
  for(i in 1:length(list_loop)){
    print(paste0("i=", i))
    for(k in 1:15){
      print(paste0("k=", k))
      RunINLACA(groupsid = k, 
         adj = list_loop[[i]]$adj, 
         rw2 = TRUE, 
         thet = NULL, 
         rerun = FALSE, 
         store = c("summary.random", "summary.fixed", "marginals.random", "mode"))
      gc()
    }
  }
}
t_1 <- Sys.time()
t_1 - t_0


# 3. Space time effect modification

# time
t_0 <- Sys.time()
if(timeEM == TRUE){
  print("timeEM")
    for(k in 1:15){
      print(paste0("k=", k))
      RunINLACA(groupsid = k, 
        adj = list_loop[[2]]$adj, 
        rw2 = FALSE, 
        thet = NULL, 
        rerun = FALSE, 
        store = c("summary.fixed"), 
        year_start = 2002, 
        year_end = 2007, 
        efmod = "time")

      RunINLACA(groupsid = k, 
        adj = list_loop[[2]]$adj, 
        rw2 = FALSE, 
        thet = NULL, 
        rerun = FALSE, 
        store = c("summary.fixed"), 
        year_start = 2008, 
        year_end = 2013, 
        efmod = "time")

      RunINLACA(groupsid = k, 
        adj = list_loop[[2]]$adj, 
        rw2 = FALSE, 
        thet = NULL, 
        rerun = FALSE, 
        store = c("summary.fixed"), 
        year_start = 2014, 
        year_end = 2019, 
        efmod = "time")
      gc()
  }
}
t_1 <- Sys.time()
t_1 - t_0


# space

# to run for y = c(5, 10, 15)

reg <- unique(asthma_dat$regid)

t_0 <- Sys.time()
if(spaceEM == TRUE){
  for(j in 1:3){
    for(i in 1:length(reg)){
      RunINLACA(groupsid = c(5, 10, 15)[j], 
                adj = list_loop[[2]]$adj, 
                rw2 = FALSE, 
                thet = NULL, 
                rerun = FALSE, 
                store = c("summary.fixed"), 
                region = reg[i], 
                efmod = "space")
    }    
  }
}
t_1 <- Sys.time()
t_1 - t_0



t_1_all <- Sys.time() 
t_1_all - t_0_all




print("code has run")


##########################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################


