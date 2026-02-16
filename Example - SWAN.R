load(data_file("data.RData"))
ls()

combined_data = readRDS(rds_file("msm_data")) %>%
  as.data.frame()

library(mstate) # Please use the version 0.2.7
library(survival)
library(JM)

ids = combined_data$SWANID %>% unique
ids = ids[1:2400]
combined_data = combined_data %>% filter(SWANID %in% ids)

n = length(ids)

from_to <- tibble(from = c(1,1,1,
                2,2,
                3),
       to = c(2,3,4,
              3,4,
              4))

trans_frame <- data.frame(SWANID = rep(unique(combined_data$SWANID),4),
       STATUS = c(rep(1,n),rep(2,n),rep(3,n),rep(4,n))) %>%
  arrange(SWANID,STATUS)

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

saveRDS(trans_wide,rds_file("trans_wide"))
trans_wide = readRDS(rds_file("trans_wide")) 
  
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
  
mp_mstate <- msprep(time = c(NA, "t_State_2","t_State_3","t_State_4"),
                        status = c(NA,"State_2","State_3","State_4"),
                        data = wide_frame,
                        trans = tmatmp,
                        id = "SWANID") 

mp_mstate = mp_mstate %>% filter(Tstart < Tstop)

saveRDS(mp_mstate,rds_file("mp_mstate"))
mp_mstate = readRDS(rds_file("mp_mstate"))

#### Longitudinal data
mp_long = combined_data %>%
  filter(SWANID %in% ids) %>%
  filter(STATUS_num < 4) %>%
  dplyr::select(SWANID,AGE_mod,FSH,E2AVE,STATUS_num) %>%
  mutate(AGE_mod = AGE_mod - 42) %>%
  filter(!is.na(FSH) & !is.na(E2AVE)) %>%
  ungroup %>%
  mutate(FSH = scale(FSH),
         E2AVE = scale(E2AVE)) %>%
  mutate(SWANID = factor(SWANID))

# mp_long_class = readRDS(rds_file("mp_long_class")) %>%
#   dplyr::select(SWANID,AGE_mod,FSH,E2AVE,class:prob3) 


mp_mstate = mp_mstate %>%
  filter(SWANID %in% as.character(mp_long$SWANID))

# Weird rounding issues in these, fix later
exclusions <- c(10464, 14334, 14392, 17209, 20411, 25090, 25207, 25400, 27811, 
                31543, 31630, 36875, 51785, 53185, 60611, 64002, 66203, 66793,
                68316, 68936, 71184, 74606, 75490, 78109,17626, 18355, 20124, 23334)

mp_mstate = mp_mstate %>% filter(!SWANID %in% exclusions)
mp_long = mp_long %>% filter(!SWANID %in% exclusions)


mp_mstate %>% group_by(SWANID) %>% n_groups()
mp_long %>% group_by(SWANID) %>% n_groups()

# Add class predictions
mp_long = mp_long %>%
  group_by(SWANID) %>%
  mutate(trunc = 1:n()) %>%
  left_join(mp_long_class %>% dplyr::select(SWANID,trunc,class:prob3) %>% mutate(SWANID = factor(SWANID)),by = c("SWANID","trunc"))

mp_long %>% summary

# saveRDS(mp_long,rds_file("mp_long_predclass"))
saveRDS(mp_long,rds_file("mp_long"))
mp_long = readRDS(rds_file("mp_long"))

# Check class membership
# mp_mstate_class = mp_mstate %>%
#   








# Load data
# mp_surv = readRDS(rds_file("mp_surv"))
mp_long = readRDS(rds_file("mp_long"))

# get list of R files
# This is all loading JointModel stuffs from the Ferrer paper.
runR = list.files(path = main_path)
runR = runR[grepl(".R",runR)]

sapply(X = runR,FUN = function(x){source(code_file(x))})


### Linear model
FSHFit <- lme(
  # fixed = FSH ~ AGE_mod,
              fixed = FSH ~ poly(AGE_mod,2),
              data = mp_long,
              random = ~ AGE_mod| SWANID,
              method = "REML",
              control = list(opt = "optim",
                             niter = 10000))

# E2Fit <- lme(fixed = E2AVE ~ ns(AGE_mod,3),
#               data = mp_long,
#               random = ~ AGE_mod| SWANID,
#               method = "REML",
#               control = list(opt = "optim"))


# Survival model
coxMP <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans),
                data = mp_mstate, 
                method = "breslow", 
                x = TRUE, 
                model = TRUE)


## Fixed and random parts that go from long. to survival model
dForm <- list(fixed = ~ c(AGE_mod,slope(AGE_mod)),
              # indFixed = c(2:3 ,5:6),
              random = ~ c(AGE_mod,slope(AGE_mod)))

dForm <- list(fixed = ~ c(AGE_mod),
              # indFixed = c(2:3 ,5:6),
              random = ~ c(AGE_mod))


# Fix joint model
fit0 <- JMstateModel(lmeObject = FSHFit,
                     survObject = coxMP,
                     timeVar = "AGE_mod",
                     parameterization = "value",
                     method = "spline-PH-GH",
                     derivForm = dForm,
                     Mstate = TRUE,
                     data.Mstate = mp_mstate,
                     ID.Mstate = "SWANID",
                     control=list(GHk = 1,lng.in.kn = 0), 
                     verbose=TRUE)

fit0 <- JMstateModel(lmeObject = FSHFit,
                     survObject = coxMP,
                     timeVar = "AGE_mod",
                     parameterization = "value",
                     method = "spline-PH-GH",
                     derivForm = dForm,
                     Mstate = TRUE,
                     data.Mstate = mp_mstate,
                     ID.Mstate = "SWANID",
                     control = list(
                       GHk = 1,
                       lng.in.kn = 0,  # no internal knots
                       tol1 = 1e-3,
                       tol2 = 1e-3,
                       itermax = 50
                     ), 
                     verbose=TRUE)




summary(fit0)

jm_example <- JMstateModel(lmeObject = FSHFit,
                           survObject = coxMP,
                           timeVar = "AGE_mod",
                           parameterization = "value",
                           method = "spline-PH-GH",
                           derivForm = dForm,
                           Mstate = TRUE,
                           interFact = list(value = ~strata(trans),  # This was - 1 in example, why?
                                            data = mp_mstate),
                           # init = fit0,
                           data.Mstate = mp_mstate,
                           ID.Mstate = "SWANID",
                           control = list(GHk = 7, lng.in.kn = 3),
                           verbose = TRUE)

summary(jm_example)

saveRDS(jm_example,rds_file("jm_example_strat_trans"))
jm_example = readRDS(rds_file("jm_example"))

# JM_pred = predict(jm_example,newdata = mp_long)




# ## Evaluation
# t0 <- 0
# tmax <- max(mp_mstate$Tstop, na.rm=TRUE)
# times <- seq(t0, tmax, length.out = 200) 
# 
# subjects <- unique(mp_mstate$SWANID)
# new_long <- expand.grid(SWANID = factor(subjects), 
#                         AGE_mod = times, 
#                         KEEP.OUT.ATTRS = FALSE, 
#                         stringsAsFactors = FALSE) %>%
#   arrange(SWANID, AGE_mod)



# mhat_vec <- predict(jm_example, newdata = new_long, type = "Subject")  # length = nrow(new_long)
# new_long$mhat <- mhat_vec

# survpred <- survfitJM(jm_example,
#           idVar = "SWANID",
#           newdata = tibble(SWANID = factor(1),trans = 1:6,AGE_mod = c(3),Tstart = c(2),Tstop = c(4),status = c(0),FSH = c(2)))
# 
# plot(survpred)

surv_data_frame_class = readRDS(rds_file("surv_data_frame_class"))
surv_data_frame_class = surv_data_frame_class %>%
  filter(from != to) %>%
  mutate(trans = case_when(from == 1 & to == 2 ~ 1,
                           from == 1 & to == 3 ~ 2,
                           from == 1 & to == 4 ~ 3,
                           from == 2 & to == 3 ~ 4,
                           from == 2 & to == 4 ~ 5,
                           from == 3 & to == 4 ~ 6))

surv_data_frame_class = surv_data_frame_class %>% 
  mutate(FSH = as.numeric(FSH),E2AVE = as.numeric(E2AVE))

survpred <- survfitJM(jm_example,
                      idVar = "SWANID",
                      newdata = surv_data_frame_class %>%
                        filter(SWANID == 30747))


plot(survpred,include.y = T)



survfitJM()
pred = predict(jm_example,
        newdata = surv_data_frame_class %>%
          filter(SWANID == 30747),
        type = "Subject",
        process = "event",
        # FtTimes = 0:5,
        # transition = 6,
        idVar = "SWANID",
        timeVar = "AGE_mod")


survpred4 = survfitJM(jm_example,
        newdata = surv_data_frame_class %>%
          filter(SWANID == 30747) %>% slice(1:4),
        survTimes = c(0:10),
        transition = 4,
        idVar = "SWANID",
        timeVar = "AGE_mod")

survpred5 = survfitJM(jm_example,
                     newdata = surv_data_frame_class %>%
                       filter(SWANID == 10514) %>% slice(1:9),
                     survTimes = c(6:10),
                     transition = 4,
                     idVar = "SWANID",
                     timeVar = "AGE_mod")
