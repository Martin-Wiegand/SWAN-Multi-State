library("JMbayes2")

by_long = readRDS(rds_file("by_long"))
by_surv = readRDS(rds_file("by_surv"))


trans_jm_func <- function(sv,by_long = by_long){
  last_time = sv %>% group_by(SWANID) %>%  reframe(last_time =max(Tstop))
  
  data_long = by_long %>% left_join(last_time,by = "SWANID") %>% filter(AGE_mod <= last_time)
  
  sv = sv %>% arrange(SWANID,Tstart,Tstop,trans)
  data_long = data_long %>% arrange(SWANID,AGE_mod)
  
  # Longmod
  FSHFit <- lme(fixed = FSH ~ ns(AGE_mod,3),
                data = data_long,
                random = ~ AGE_mod | SWANID,
                method = "REML",
                control = list(opt = "optim",
                               niter = 10000))
  
  E2Fit <- lme(fixed = E2AVE ~ ns(AGE_mod,3),
               data = data_long,
               random = ~ AGE_mod | SWANID,
               method = "REML",
               control = list(opt = "optim"))
  
  # Cox model
  TransCox <- coxph(Surv(Tstart, Tstop, status) ~ prob2 + prob3 + prob4 + prob5 + SLEEP + VMS + strata(trans) + cluster(SWANID), 
                    data=sv,
                    x = T,
                    model = T)
  
  fForms <- list("FSH" = ~ value(FSH) + slope(FSH),
                 "E2AVE" = ~ value(E2AVE) + slope(E2AVE)) 
  
  innerfit <- jm(Mixed_objects = list(FSHFit,E2Fit),
            Surv_object = TransCox,
            functional_forms = fForms,
            id_var = "SWANID",
            time_var = "AGE_mod") 
  
  return(list(FSHFit = FSHFit,
              E2Fit = E2Fit,
              TransCox = TransCox,
              JMmodel = innerfit))
}

res_state1 = trans_jm_func(sv = by_surv %>% filter(trans == 1),by_long = by_long)
res_state2 = trans_jm_func(sv = by_surv %>% filter(trans == 2),by_long = by_long)
res_state3 = trans_jm_func(sv = by_surv %>% filter(trans == 3),by_long = by_long)
res_state4 = trans_jm_func(sv = by_surv %>% filter(trans == 4),by_long = by_long)
res_state5 = trans_jm_func(sv = by_surv %>% filter(trans == 5),by_long = by_long)
res_state6 = trans_jm_func(sv = by_surv %>% filter(trans == 6),by_long = by_long)

saveRDS(res_state1,rds_file("res_state1"))
saveRDS(res_state2,rds_file("res_state2"))
saveRDS(res_state3,rds_file("res_state3"))
saveRDS(res_state4,rds_file("res_state4"))
saveRDS(res_state5,rds_file("res_state5"))
saveRDS(res_state6,rds_file("res_state6"))

res_state1 = readRDS(rds_file("res_state1"))
res_state2 = readRDS(rds_file("res_state2"))
res_state3 = readRDS(rds_file("res_state3"))
res_state4 = readRDS(rds_file("res_state4"))
res_state5 = readRDS(rds_file("res_state5"))
res_state6 = readRDS(rds_file("res_state6"))

by_long = by_long %>% group_by(SWANID) %>% mutate(newID = cur_group_id())
by_surv = by_surv %>% group_by(SWANID) %>% mutate(newID = cur_group_id())

surv_by_row1 = by_surv %>% filter(trans == 1) %>%  distinct()
surv_by_row2 = by_surv %>% filter(trans == 2) %>%  distinct()
surv_by_row3 = by_surv %>% filter(trans == 3) %>%  distinct()
surv_by_row4 = by_surv %>% filter(trans == 4) %>%  distinct()
surv_by_row5 = by_surv %>% filter(trans == 5) %>%  distinct()
surv_by_row6 = by_surv %>% filter(trans == 6) %>% distinct()

epsilon = 0.000000001


# wrapper =  function(leckmichdoch,SWANID,trans,Tstart,Tstop,FSH,E2AVE,prob1,prob2,prob3,prob4,prob5,SLEEP,VMS,fit){
wrapper =  function(surv,fit){
  # sdat = tibble(ID = leckmichdoch,trans = trans,Tstart = Tstart,Tstop = Tstop,FSH = FSH,E2AVE = E2AVE,SWANID = SWANID,
  #               prob1 = prob1,prob2 = prob2,prob3 = prob3,prob4 = prob4,prob5 = prob5,SLEEP = SLEEP,VMS = VMS,status = F)%>%
  #   mutate(Tstop = Tstart + epsilon,
  #          status = 0)
  
  sdat = surv %>%
    mutate(Tpred = Tstop,
           Tstop = Tstart + epsilon,
           status = F)
  
  ldat = by_long  %>%
    filter(SWANID == (surv$SWANID),
           as.numeric(AGE_mod) <= as.numeric(surv$Tstart))
    

  temp = predict(fit,
                 newdata = list(newdataL = ldat,
                                newdataE = sdat),
                 process = "event",
                 times = sdat$Tpred,
                 type_pred = "link",
                 type = "subject_specific",
                 return_newdata = F)
  
  return(tibble(SWANID = surv$SWANID,Tstop = temp$times,Prob = temp$pred,trans = surv$trans)[2,])
}

## Trans 1
results1 <- list()

for(i in 1:nrow(surv_by_row1)){
  print(i)
  results1[[i]] = wrapper(surv_by_row1[i,],fit = res_state1$JMmodel)
}

saveRDS(results1,rds_file("results1"))

# Trans 2
results2 <- list()

for(i in 1:nrow(surv_by_row2)){
  print(i)
  results2[[i]] = wrapper(surv_by_row2[i,],fit = res_state2$JMmodel)
}

saveRDS(results2,rds_file("results2"))

# Trans 3
results3 <- list()

for(i in 1:nrow(surv_by_row3)){
  print(i)
  results3[[i]] = wrapper(surv_by_row3[i,],fit = res_state3$JMmodel)
}

saveRDS(results3,rds_file("results3"))


# Trans 4
results4 <- list()

for(i in 1:nrow(surv_by_row4)){
  print(i)
  results4[[i]] = wrapper(surv_by_row4[i,],fit = res_state4$JMmodel)
}

saveRDS(results4,rds_file("results4"))


# Trans 5
results5 <- list()

for(i in 1:nrow(surv_by_row5)){
  print(i)
  results5[[i]] = wrapper(surv_by_row5[i,],fit = res_state5$JMmodel)
}

saveRDS(results5,rds_file("results5"))



# Trans 6
results6 <- list()

for(i in 1:nrow(surv_by_row6)){
  print(i)
  results6[[i]] = wrapper(surv_by_row6[i,],fit = res_state6$JMmodel)
}

saveRDS(results6,rds_file("results6"))



res_tibble1 = do.call(rbind,results1) %>% as_tibble()
by_prob1 = surv_by_row1 %>%
  left_join(res_tibble1,by = c("SWANID","Tstop","trans"))

res_tibble2 = do.call(rbind,results2) %>% as_tibble()
by_prob2 = surv_by_row2 %>%
  left_join(res_tibble2,by = c("SWANID","Tstop","trans"))

res_tibble3 = do.call(rbind,results3) %>% as_tibble()
by_prob3 = surv_by_row3 %>%
  left_join(res_tibble3,by = c("SWANID","Tstop","trans"))

res_tibble4 = do.call(rbind,results4) %>% as_tibble()
by_prob4 = surv_by_row4 %>%
  left_join(res_tibble4,by = c("SWANID","Tstop","trans"))

res_tibble6 = do.call(rbind,results6) %>% as_tibble()
by_prob6 = surv_by_row6 %>%
  left_join(res_tibble6,by = c("SWANID","Tstop","trans"))



##
library("pROC")
library("MLmetrics")

auc(roc(by_prob1$status,by_prob1$Prob))
auc(roc(by_prob2$status,by_prob2$Prob))
auc(roc(by_prob3$status,by_prob3$Prob))
auc(roc(by_prob4$status,by_prob4$Prob))
auc(roc(by_prob5$status,by_prob5$Prob))
auc(roc(by_prob6$status,by_prob6$Prob))



for(i in 1:6){
  data = by_prob %>%
    filter(trans == i)
  
  print(c(as.numeric(auc(roc(data$status,data$Prob))),as.numeric(PRAUC(y_pred = data$Prob,y_true = data$status)),as.numeric(data %>% group_by(trans) %>% reframe(mean((status - Prob)^2)))[2]))
}


ggplot(by_prob) +
  geom_point(aes(x = Prob,y = status)) +
  facet_wrap(. ~ trans)



# prevalence
fit = readRDS(rds_file("JM_fit_ms_class"))

by_long$AGE_mod %>% range

times = 0:20

new_surv <- expand.grid(
  SWANID = 10056,
  trans  = 1:6,
  Tstop  = times
) |>
  mutate(
    Tstart = 0,
    status = FALSE
  )

CIF <- predict(fit,
               newdata = list(newdataL = by_long, newdataE = new_surv),
               process = "event",
               type = "mean_subject")
