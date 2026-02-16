library(mstate) # Please use the version 0.2.7
library(survival)
library(JM)


mp_mstate = readRDS(rds_file("mp_mstate"))
mp_long = readRDS(rds_file("mp_long"))

# Weird rounding issues in these, fix later
exclusions <- c(10464, 14334, 14392, 17209, 20411, 25090, 25207, 25400, 27811, 
                31543, 31630, 36875, 51785, 53185, 60611, 64002, 66203, 66793,
                68316, 68936, 71184, 74606, 75490, 78109,17626, 18355, 20124, 23334)

mp_mstate = mp_mstate %>% filter(!SWANID %in% exclusions)
mp_long = mp_long %>% filter(!SWANID %in% exclusions)


mp_mstate %>% group_by(SWANID) %>% n_groups()
mp_long %>% group_by(SWANID) %>% n_groups()

mp_long = mp_long %>%
  mutate(FSH = as.numeric(FSH))%>%
  mutate(E2AVE = as.numeric(E2AVE))

mp_mstate_1 = mp_mstate %>% filter(trans == 1)
mp_mstate_2 = mp_mstate %>% filter(trans == 2)
mp_mstate_3 = mp_mstate %>% filter(trans == 3)
mp_mstate_4 = mp_mstate %>% filter(trans == 4)
mp_mstate_5 = mp_mstate %>% filter(trans == 5)
mp_mstate_6 = mp_mstate %>% filter(trans == 6)

mp_long_1 = mp_long %>% 
  filter(SWANID %in% mp_mstate_1$SWANID)

# get list of R files
# This is all loading JointModel stuffs from the Ferrer paper.
# runR = list.files(path = main_path)
# runR = runR[grepl(".R",runR)]
# 
# sapply(X = runR,FUN = function(x){source(code_file(x))})




### Linear model
FSHFit1 <- lme(
  # fixed = FSH ~ AGE_mod,
  fixed = FSH ~ poly(AGE_mod,2),
  data = mp_long_1,
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
coxM1P <- coxph(Surv(Tstart, Tstop, status) ~ cluster(SWANID),
               data = mp_mstate_1, 
               method = "breslow", 
               x = TRUE, 
               model = TRUE)

coxM2P <- coxph(Surv(Tstart, Tstop, status) ~ 1,
               data = mp_mstate_2, 
               method = "breslow", 
               x = TRUE, 
               model = TRUE)



# Fix joint model
dForm <- list(fixed = ~ 1 + FSH, 
              # indFixed = c(2, 4), 
              random = ~ 1)

fit1 <- JM::jointModel(lmeObject = FSHFit1,
                       survObject = coxM1P,
                       timeVar = "AGE_mod",
                       parameterization = "value")

summary(fit1)
