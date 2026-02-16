
combined_data = readRDS(rds_file("msm_data")) %>%
  as.data.frame()

mp_class_m5 = readRDS(rds_file("mp_class_m5"))
surv_data_frame_class = readRDS(rds_file("surv_data_frame_class_mult5")) %>%
  mutate(trans = factor(trans))

validation_combined <- combined_data %>%
  filter(!SWANID %in% mp_class_m5$SWANID) %>%
  as_tibble()

library(mstate) # Please use the version 0.2.7
library(survival)
library(JM)


n = validation_combined %>% group_by(SWANID) %>% n_groups()

from_to <- tibble(from = c(1,1,1,
                           2,2,
                           3),
                  to = c(2,3,4,
                         3,4,
                         4))

trans_frame <- data.frame(SWANID = rep(unique(validation_combined$SWANID),4),
                          STATUS = c(rep(1,n),rep(2,n),rep(3,n),rep(4,n))) %>%
  arrange(SWANID,STATUS)

# trans times
trans_times <-  validation_combined %>% 
  dplyr::select(SWANID,AGE_mod,STATUS_num) %>%
  group_by(SWANID) %>%
  arrange(SWANID,AGE_mod) %>%
  group_by(SWANID,STATUS_num) %>%
  reframe(Start = min(AGE_mod),
          Stop = max(AGE_mod)) %>%
  rename(STATUS = STATUS_num)


# Add transition times to raw frame long format
trans_int <- trans_frame %>%
  left_join(trans_times,by = c("SWANID","STATUS")) %>%
  group_by(SWANID) %>%
  mutate(Reached = ifelse(is.na(Start),0,1)) %>%
  fill(Start,.direction = "up") %>%
  fill(Stop,.direction = "down") %>%
  mutate(Start = case_when(is.na(Start) ~ Stop,
                           TRUE ~ Start))


# wide format
trans_wide <- trans_int %>%
  dplyr::select(SWANID,STATUS,Start,Reached) %>%
  pivot_wider(id_cols = SWANID,names_from = STATUS,values_from = c(Start,Reached))

saveRDS(trans_wide,rds_file("trans_wide_val"))
trans_wide = readRDS(rds_file("trans_wide_val")) 

wide_frame = trans_wide %>% 
  dplyr::select(-Start_1,-Reached_1) %>%
  rename(t_State_2 = Start_2,
         t_State_3 = Start_3,
         t_State_4 = Start_4,
         State_2 = Reached_2,
         State_3 = Reached_3,
         State_4 = Reached_4) %>%
  as.data.frame() %>%
  mutate(t_State_2 = t_State_2 - 42,
         t_State_3 = t_State_3 - 42,
         t_State_4 = t_State_4 - 42)

tmatmp <- matrix(NA, 4, 4)
tmatmp[1, 2:4] <- 1:3
tmatmp[2, 3:4] <- 4:5
tmatmp[3, 4] <- 6
dimnames(tmatmp) <- list(from = c("State_1", "State_2", "State_3","State_4"),
                         to = c("State_1", "State_2", "State_3","State_4"))

mp_mstate_val <- msprep(time = c(NA, "t_State_2","t_State_3","t_State_4"),
                    status = c(NA,"State_2","State_3","State_4"),
                    data = wide_frame,
                    trans = tmatmp,
                    id = "SWANID") 

mp_mstate_val = mp_mstate_val %>% filter(Tstart < Tstop)

saveRDS(mp_mstate_val,rds_file("mp_mstate_val"))
mp_mstate_val = readRDS(rds_file("mp_mstate_val"))

#### Longitudinal data
mp_long_val = validation_combined %>%
  filter(STATUS_num < 4) %>%
  dplyr::select(SWANID,AGE_mod,FSH,E2AVE,STATUS_num) %>%
  # mutate(AGE_mod = AGE_mod - 42) %>%
  filter(!is.na(FSH) & !is.na(E2AVE)) %>%
  ungroup %>%
  mutate(SWANID = factor(SWANID))

mp_mstate_val = mp_mstate_val %>%
  filter(SWANID %in% as.character(mp_long_val$SWANID))


mp_mstate_val %>% group_by(SWANID) %>% n_groups()
mp_long_val %>% group_by(SWANID) %>% n_groups()

## Add LCMM group probabilities
lcmm_data_val = mp_long_val %>%
  dplyr::select(-STATUS_num) %>%
  group_by(SWANID) %>%
  mutate(ID = cur_group_id())

truncator <- function(x,data = lcmm_data_val){
  data %>% 
    as_tibble() %>%
    group_by(SWANID) %>%
    filter(n() >= x) %>%
    filter(1:n() <= x) %>%
    mutate(trunc = x) %>%
    return()
}

trunc_long_val = 
  bind_rows(truncator(1),
            truncator(2),
            truncator(3),
            truncator(4),
            truncator(5),
            truncator(6),
            truncator(7),
            truncator(8),
            truncator(9),
            truncator(10),
            truncator(11)) %>%
  arrange(SWANID,trunc,AGE_mod) %>%
  mutate(truncid = as.numeric(paste0(SWANID,trunc)))

trunc_long_val = trunc_long_val %>% arrange(SWANID,AGE_mod,trunc)
lcmm_m5 = readRDS(rds_file("lcmm_m5"))

trunc_pivot_val = trunc_long_val %>% 
  pivot_longer(cols = c("FSH","E2AVE"),names_to = "marker",values_to = "value")

m5 = readRDS(rds_file("m5"))

library("lcmm")
pred_m5_val = predictClass(model = m5,
                       newdata = trunc_pivot_val,
                       subject = "truncid")

lcmm_m5_val = lcmm_data_val %>%
  group_by(SWANID) %>%
  mutate(truncid = paste0(SWANID,1:n())) %>%
  left_join(pred_m5_val,
            by = "truncid")


## Join

mp_long_val = mp_long_val %>% left_join(lcmm_m5_val %>% dplyr::select(SWANID,AGE_mod,prob1:prob5),by = c("SWANID","AGE_mod"))


clinical_symptoms = readRDS(rds_file("clinical_symptoms")) 
mp_long_val = mp_long_val %>%
  group_by(SWANID) %>%
  mutate(VISIT = 0:(n()-1)) %>%
  left_join(clinical_symptoms %>% mutate(SWANID = factor(SWANID)) %>% dplyr::select(SWANID,VISIT,VMS,SLEEP),by = c("SWANID","VISIT"))

saveRDS(mp_long_val,rds_file("mp_long_val_class"))
mp_long_val = readRDS(rds_file("mp_long_val_class"))

# Surv frame for evaluation
sv = validation_combined %>%
  dplyr::select(SWANID,AGE_mod,STATUS_num)

surv_data_frame_val = sv %>%
  group_by(SWANID) %>%
  mutate(Tstart = AGE_mod) %>%
  mutate(Tstop = lead(AGE_mod)) %>%
  mutate(nexstate = lead(STATUS_num)) %>%
  mutate(repetitions = case_when(STATUS_num == 1 ~ 3,
                                 STATUS_num == 2 ~ 2,
                                 T ~ 1)) %>%
  uncount(repetitions) %>%
  group_by(SWANID,AGE_mod) %>%
  mutate(to = 1:n()+STATUS_num) %>%
  rename(from = STATUS_num)  %>%
  filter(from != 4) %>%
  mutate(status = (nexstate == to)) %>%
  dplyr::select(SWANID,AGE_mod,from,to,Tstart,Tstop,status)

surv_data_frame_val %>% print(n = 100)
surv_data_frame_val = surv_data_frame_val %>%
  filter(!is.na(status)) %>%
  mutate(SWANID = factor(SWANID))

surv_data_frame_val = surv_data_frame_val%>% left_join(mp_long_val %>% dplyr::select(SWANID,AGE_mod,prob1:SLEEP),by = c("SWANID","AGE_mod"))
surv_data_frame_val = surv_data_frame_val %>% filter(!is.na(prob1) & !is.na(VMS) & !is.na(SLEEP)) 

### Evaluation
library("JMbayes2")
library("JM")
library("joineRML")

fit = readRDS(rds_file("JM_fit_ms_class"))
mp_mstate_val = as_tibble(mp_mstate_val)


surv_data_frame_val %>% group_by(SWANID)
mp_long_val%>% group_by(SWANID)

surv_data_frame_val = surv_data_frame_val %>%
  mutate(trans = case_when(from == 1 & to == 2 ~ 1,
                           from == 1 & to == 3 ~ 2,
                           from == 1 & to == 4 ~ 3,
                           from == 2 & to == 3 ~ 4,
                           from == 2 & to == 4 ~ 5,
                           from == 3 & to == 4 ~ 6,))

surv_data_frame_val = surv_data_frame_val %>%
  left_join(mp_long_val %>% dplyr::select(SWANID,AGE_mod,FSH,E2AVE),by = c("SWANID","AGE_mod"))


saveRDS(surv_data_frame_val,rds_file("sdat_val"))
saveRDS(mp_long_val,rds_file("ldat_val"))

surv_data_frame_val = readRDS(rds_file("sdat_val"))
mp_long_val = readRDS(rds_file("ldat_val"))

mp_long_val = mp_long_val %>% left_join(surv_data_frame_val %>%
                            group_by(SWANID) %>%
                            reframe(maxtime = max(Tstop)),
                          by = "SWANID") %>%
  filter(AGE_mod < maxtime)

temp = predict(fit, 
               newdata = list(newdataL = mp_long_val,
                              newdataE = surv_data_frame_val),
               process = "event", 
               times = surv_data_frame_val$Tstop,
               type_pred = "link",
               type = "subject_specific",
               return_newdata = F)


saveRDS(temp,rds_file("val_probs_predicted"))

ids = mp_long_val$SWANID %>% unique 


for(i in 1:100){
  print(i)
  
  temptest = predict(fit, 
                     newdata = list(newdataL = mp_long_val %>% filter(SWANID %in% ids[i]),
                                    newdataE = surv_data_frame_val %>% filter(SWANID %in% ids[i]) ),
                     process = "event", 
                     # times = surv_data_frame_val %>% filter(SWANID %in% ids[i])  %>% pull(Tstop)%>% slice(1),
                     times = c(46,46,46,47,47,48,48) +0.5,
                     type_pred = "link",
                     type = "subject_specific",
                     return_newdata = F)
}






surv_by_row = surv_data_frame_val %>%
  dplyr::select(-trans,-status) %>%
  distinct()

epsilon = 0.0000001

wrapper =  function(x){
  
  ldat = mp_long_val %>%
    filter(SWANID %in% x$SWANID,
           AGE_mod <= x$AGE_mod) %>%
    mutate(SWANID = 10056)
  
  # sdat = tibble(as_tibble(matrix(x,nrow = 6,ncol = length(x),byrow = T)),status = rep(F,6),trans = 1:6)
  sdat = tibble(x %>% slice(rep(1,6)),status = rep(F,6),trans = 1:6)
  # names(sdat) = c("SWANID","AGE_mod","from","to","Tstart","Tstop","prob1","prob2","prob3","prob4","prob5","VISIT","VMS","SLEEP","FSH","E2AVE","status","trans")
  
  sdat = sdat %>%
    mutate(Tstop = Tstart + epsilon) %>%
    mutate(SWANID = 10056)
  
  temp = predict(fit, 
                 newdata = list(newdataL = ldat,
                                newdataE = sdat),
                 process = "event", 
                 times = x$Tstop,
                 type_pred = "link",
                 type = "subject_specific",
                 return_newdata = F)
  
  return(tibble(SWANID = temp$id,Tstop = temp$times,Prob = temp$pred,trans = temp$`_strata`)[2*(1:6),])
}

## Parallel computation
res = apply(FUN = wrapper,X = surv_by_row,MARGIN = 1)


saveRDS(res,rds_file("res"))
res = readRDS(rds_file("res"))

x = surv_by_row[7,]
