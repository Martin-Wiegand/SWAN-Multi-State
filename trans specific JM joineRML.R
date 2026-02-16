library("joineRML")
library("JMbayes2")

clinical_symptoms = readRDS(rds_file("clinical_symptoms"))
combined_data = readRDS(rds_file("msm_data"))
# mp_class_m5 = readRDS(rds_file("mp_class_m5"))
lcmm_m4_all = readRDS(rds_file("lcmm_m4_all"))

# Agreed no clin. symptoms after all
# lcmm_m5_all %>%
#   group_by(SWANID) %>%
#   mutate(VISIT = cumsum(replace_na(round(AGE_mod - lag(AGE_mod),digits = 0),replace = 0))) %>%
#   left_join(clinical_symptoms %>% dplyr::select(SWANID,VISIT,VMS,S),by = "SWANID")

surv_data_frame_class = readRDS(rds_file("surv_data_frame_class_mult5")) %>%
  mutate(trans = factor(trans))

minimal_surv_data = combined_data %>%
  group_by(SWANID) %>%
  dplyr::select(SWANID,VISIT,AGE_mod,STATUS_num,FSH,E2AVE) %>%
  mutate(nextstate = lead(STATUS_num)) %>%
  mutate(Tstop = lead(AGE_mod)) %>%
  mutate(Tstart = AGE_mod) %>%
  rename(from = STATUS_num)

transmissions = tibble(from = c(1,1,1,
                                2,2,
                                3),
                       to = c(2,3,4,
                             3,4,
                             4),
                       trans = factor(1:6))

full_surv_frame = minimal_surv_data %>%
  right_join(transmissions,by = "from",relationship = "many-to-many") %>%
  group_by(SWANID,VISIT,from) %>%
  arrange(SWANID,VISIT,from,to) %>%
  mutate(status = (to == nextstate)) %>%
  dplyr::select(SWANID,VISIT,from,to,status,nextstate,Tstart,Tstop,AGE_mod,everything()) %>%
  filter(from != 4) %>%
  group_by(SWANID) %>%
  filter(!is.na(nextstate)) 


full_surv_class = full_surv_frame %>%
  left_join(lcmm_m4_all %>% dplyr::select(SWANID,AGE_mod,prob1:prob4),by = c("SWANID","Tstart" = "AGE_mod"))

saveRDS(full_surv_class,rds_file("full_surv_class"))


## Create production MS frames


by_long = lcmm_m4_all %>% filter(!is.na(FSH),!is.na(E2AVE),!is.na(prob1)) %>%
  mutate(AGE_mod = AGE_mod - 42)

by_surv = full_surv_class %>%
  filter(SWANID %in% lcmm_m4_all$SWANID) %>%
  mutate(Tstart = Tstart - 42,
         Tstop = Tstop -42) %>%
  group_by(SWANID) %>%
  arrange(AGE_mod) %>%
  fill(prob1,prob2,prob3,prob4,.direction = "down") %>%
  dplyr::select(SWANID,trans,Tstart,Tstop,status,prob1:prob4,FSH,E2AVE) %>%
  drop_na()

by_long = by_long %>% filter(SWANID %in% by_surv$SWANID)
by_surv = by_surv %>% filter(SWANID %in% by_long$SWANID)

last_time = by_surv %>%
  group_by(SWANID) %>%
  reframe(last_time  =max(Tstop))
  
by_long %>% left_join(last_time,by = "SWANID") %>% filter(AGE_mod <= last_time)

by_long %>% group_by(SWANID) %>% n_groups()
by_surv %>% group_by(SWANID) %>% n_groups()

by_long %>% summary
by_surv %>% summary

by_surv = by_surv %>% arrange(SWANID,Tstart,Tstop,trans)
by_long = by_long %>% arrange(SWANID,AGE_mod)

saveRDS(by_surv,rds_file("by_surv_full"))
saveRDS(by_long,rds_file("by_long_full"))

by_surv = readRDS(rds_file("by_surv_full"))
by_long = readRDS(rds_file("by_long_full"))

# Long
FSHFit <- lme(fixed = FSH ~ ns(AGE_mod,3),
              data = by_long,
              random = ~ AGE_mod | SWANID,
              method = "REML",
              control = list(opt = "optim",
                             niter = 10000))

E2Fit <- lme(fixed = E2AVE ~ ns(AGE_mod,3),
              data = by_long,
              random = ~ AGE_mod | SWANID,
              method = "REML",
              control = list(opt = "optim"))

summary(FSHFit)
summary(E2Fit)

# Cox model
TransCox <- coxph(Surv(Tstart, Tstop, status) ~ prob1 + prob3 + prob4 + strata(trans) + cluster(SWANID), 
                  data=by_surv,
                  x = T,
                  model = T)

TransCoxVar <- coxph(Surv(Tstart, Tstop, status) ~ prob1 + prob3 + prob4 + scale(FSH) + scale(E2AVE) + strata(trans) + cluster(SWANID),
                  data=by_surv,
                  x = T,
                  model = T)

summary(TransCox)

# Define both slope & value
# Tried with slopes
fForms <- list("FSH" = ~ value(FSH) + slope(FSH),
               "E2AVE" = ~ value(E2AVE) + slope(E2AVE)) 

# Better without
fForms <- list("FSH" = ~ value(FSH),
               "E2AVE" = ~ value(E2AVE)) 

fit <- jm(Mixed_objects = list(FSHFit,E2Fit),
          Surv_object = TransCox,
          functional_forms = fForms,
          id_var = "SWANID",
          time_var = "AGE_mod") 

summary(fit)

# saveRDS(fit,rds_file("JM_fit_full_no_slope"))
# saveRDS(fit,rds_file("JM_fit_full_slope"))

fit = readRDS(rds_file("JM_fit_full_no_slope"))



# Get predictions out of model
# pred_fit = predict(fit,newdata = list(newdataL = by_long %>% filter(SWANID == 10056),
#                                       newdataE = by_surv %>% filter(SWANID == 10056)),process = "event")

saveRDS(pred_fit,rds_file("pred_fit"))


##### Pred funciton
surv_by_row = by_surv %>%
  dplyr::select(-trans,-status) %>%
  distinct()

epsilon = 0.0000001

 wrapper =  function(x){
    
    ldat = by_long %>%
      filter(SWANID %in% x[1],
             AGE_mod <= x[2])
    
    sdat = tibble(as_tibble(matrix(as.numeric(x),nrow = 6,ncol = length(x),byrow = T)),status = rep(F,6),trans = 1:6)
    names(sdat) = c("SWANID","Tstart","Tstop","prob1","prob2","prob3","prob4","prob5","SLEEP","VMS","FSH","E2AVE","status","trans")
    
    sdat = sdat %>%
      mutate(Tstop = Tstart + epsilon)
    
    temp = predict(fit, 
                   newdata = list(newdataL = ldat,
                                  newdataE = sdat),
                   process = "event", 
                   times = x[3],
                   type_pred = "link",
                   type = "subject_specific",
                   return_newdata = F)
    
    return(tibble(SWANID = temp$id,Tstop = temp$times,Prob = temp$pred,trans = temp$`_strata`)[2*(1:6),])
  }
 
 ## Parallel computation
library("parallel") 
detectCores()
cl <- makeCluster(11)
clusterExport(cl = cl,varlist = ls(envir = .GlobalEnv))
clusterEvalQ(cl = cl,library("tidyverse"))
clusterEvalQ(cl = cl,library("JMbayes2"))

stopwatch <- Sys.time()
res = parApply(cl = cl,FUN = wrapper,X = surv_by_row,MARGIN = 1)
# res = apply(FUN = wrapper,X = surv_by_row,MARGIN = 1)
difftime(Sys.time(),stopwatch)

saveRDS(res,rds_file("res"))
res = readRDS(rds_file("res"))

res_tibble = do.call(rbind,res)

res_tibble = res_tibble %>% as_tibble()

res_tibble = res_tibble %>%
  mutate(trans = as.factor(substr(trans,nchar(trans),nchar(trans))))

by_prob = by_surv %>%
  left_join(res_tibble,by = c("SWANID","Tstop","trans"))

summary(by_prob)

##
library("pROC")
library("MLmetrics")

auc_res = auc(roc(by_prob$status,by_prob$Prob))
ci.auc(auc_res)

PRAUC(y_pred = by_prob$Prob,y_true = by_prob$status)

for(i in 1:6){
  data = by_prob %>%
    filter(trans == i)
  
  print(c(i,mean(data$status),as.numeric(auc(roc(data$status,data$Prob))),
          as.numeric(PRAUC(y_pred = data$Prob,y_true = data$status)),
          as.numeric(data %>% group_by(trans) %>% reframe(mean((status - Prob)^2)))[2]))
}


ggplot(by_prob) +
  geom_point(aes(x = Prob,y = status)) +
  facet_wrap(. ~ trans)









### Cox model only
bh = basehaz(TransCoxVar,centered = F)
pred_trans = predict(TransCoxVar)

cox_surv = by_surv %>% ungroup %>% mutate(transprob = pred_trans)



ggplot(bh) + geom_line(aes(x = time + 42,y = hazard,colour = strata)) +
  theme_minimal() +
  labs(x = "Age in years",
       y = "Hazard",
       colour = "Transition") +
  theme(legend.position = "bottom") 


lp = predict(TransCoxVar,newdata = by_surv,type = "lp")
bh = as_tibble(bh)
by_surv_lp = by_surv %>% ungroup %>% mutate(lp = lp)

dat1 = bh %>% filter(strata == 1) 
dat2 = bh %>% filter(strata == 2) 
dat3 = bh %>% filter(strata == 3) 
dat4 = bh %>% filter(strata == 4) 
dat5 = bh %>% filter(strata == 5) 
dat6 = bh %>% filter(strata == 6) 

af1 = approxfun(x = dat1$time,y = dat1$hazard,method = "linear",yleft = 0)
af2 = approxfun(x = dat2$time,y = dat2$hazard,method = "linear",yleft = 0)
af3 = approxfun(x = dat3$time,y = dat3$hazard,method = "linear",yleft = 0)
af4 = approxfun(x = dat4$time,y = dat4$hazard,method = "linear",yleft = 0)
af5 = approxfun(x = dat5$time,y = dat5$hazard,method = "linear",yleft = 0)
af6 = approxfun(x = dat6$time,y = dat6$hazard,method = "linear",yleft = 0)

get_hazard <- function(timex,stratax){
  # dat = bh %>% filter(strata == stratax) 
  # af = approxfun(x = dat$time,y = dat$hazard,method = "linear",yleft = 0)
  
  if(stratax == 1) return(af1(timex))
  if(stratax == 2) return(af2(timex))
  if(stratax == 3) return(af3(timex))
  if(stratax == 4) return(af4(timex))
  if(stratax == 5) return(af5(timex))
  if(stratax == 6) return(af6(timex))
}

# get_hazard(9 ,1)

by_surv_lp = by_surv_lp %>%
  rowwise() %>%
  mutate(haz_delta = get_hazard(Tstop,trans)-get_hazard(Tstart,trans))

by_surv_lp = by_surv_lp %>%
  mutate(exp_lp = exp(lp)) %>%
  mutate(Hik = exp_lp*haz_delta)

bh$H <- bh$hazard * exp(lp)

probtrans(TransCoxVar)



by_surv_lp = by_surv_lp %>%
  mutate(from = case_when(trans %in% c(1,2,3) ~ 1,
                          trans %in% c(4,5) ~ 2,
                          trans == 6 ~ 3))%>%
  mutate(to = case_when(trans %in% c(3,5,6) ~ 4,
                        trans %in% c(2,4) ~ 3,
                          trans == 1 ~ 2))

by_surv_lp = by_surv_lp %>%
  group_by(SWANID,from) %>%
  mutate(Htot = sum(Hik)) %>%
  mutate(prob = (Hik/Htot)*(1-exp(-Htot)))

# Finished COX model!
saveRDS(by_surv_lp,rds_file("by_surv_lp"))

library(pROC)
library(DescTools)
library(MLmetrics)

by_surv_lp = by_surv_lp %>% filter(!is.na(prob))

auc(roc(by_surv_lp$status,by_surv_lp$prob))
BrierScore(by_surv_lp$status,by_surv_lp$prob)
PRAUC(by_surv_lp$prob,by_surv_lp$status) 

par(mfrow = c(2,3))

for(i in 1:6){
  data = by_surv_lp %>% filter(trans == i)
  auc(roc(data$status,data$prob)) %>% print
  BrierScore(data$status,data$prob) %>% print
  PRAUC(data$prob,data$status) %>% print
  mean(data$status) %>% print
  
  # val.prob.ci.2(data$prob,data$status)
}

install.packages("CalibrationCurves")
library("CalibrationCurves")

val.prob.ci.2(by_surv_lp$prob,by_surv_lp$status)
