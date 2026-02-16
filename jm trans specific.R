library("JM")
library("lme4")
library("joineRML")


# load
mp_class_m5 = readRDS(rds_file("mp_class_m5"))
surv_data_frame_class = readRDS(rds_file("surv_data_frame_class_mult5")) %>%
  mutate(trans = factor(trans))

AGE_spline_basis <- ns(mp_class_m5$AGE_mod, df = 3)
mp_class_m5$AGE_spline = AGE_spline_basis

# id_test = (surv_data_frame_class$SWANID %>% unique)[1:100]

sf1 = surv_data_frame_class %>% filter(trans == 1) %>%
  dplyr::select(SWANID,Tstart,Tstop,from,to,trans,status,FSH,E2AVE,prob1:prob5) %>% drop_na
sf2 = surv_data_frame_class %>% filter(trans == 2)%>%
  dplyr::select(SWANID,Tstart,Tstop,from,to,trans,status,FSH,E2AVE,prob1:prob5) %>% drop_na
sf3 = surv_data_frame_class %>% filter(trans == 3)
sf4 = surv_data_frame_class %>% filter(trans == 4)
sf5 = surv_data_frame_class %>% filter(trans == 5)
sf6 = surv_data_frame_class %>% filter(trans == 6)

FSH_data_1 =  mp_class_m5 %>% filter(!is.na(FSH)) %>% filter(SWANID %in% sf1$SWANID)
FSH_data_2 =  mp_class_m5 %>% filter(!is.na(FSH)) %>% filter(SWANID %in% sf2$SWANID)
FSH_data_3 =  mp_class_m5 %>% filter(!is.na(FSH)) %>% filter(SWANID %in% sf3$SWANID)
FSH_data_4 =  mp_class_m5 %>% filter(!is.na(FSH)) %>% filter(SWANID %in% sf4$SWANID)
FSH_data_5 =  mp_class_m5 %>% filter(!is.na(FSH)) %>% filter(SWANID %in% sf5$SWANID)
FSH_data_6 =  mp_class_m5 %>% filter(!is.na(FSH)) %>% filter(SWANID %in% sf6$SWANID)


ids_1 = intersect(sf1$SWANID,FSH_data_1$SWANID)
ids_2 = intersect(sf2$SWANID,FSH_data_2$SWANID)
ids_3 = intersect(sf3$SWANID,FSH_data_3$SWANID)
ids_4 = intersect(sf4$SWANID,FSH_data_4$SWANID)
ids_5 = intersect(sf5$SWANID,FSH_data_5$SWANID)
ids_6 = intersect(sf6$SWANID,FSH_data_6$SWANID)

sf1 = sf1 %>% filter(SWANID %in% ids_1)
sf2 = sf2 %>% filter(SWANID %in% ids_2)
sf3 = sf3 %>% filter(SWANID %in% ids_3)
sf4 = sf4 %>% filter(SWANID %in% ids_4)
sf5 = sf5 %>% filter(SWANID %in% ids_5)
sf6 = sf6 %>% filter(SWANID %in% ids_6)

FSH_data_1 = FSH_data_1 %>% filter(SWANID %in% ids_1)
FSH_data_2 = FSH_data_2 %>% filter(SWANID %in% ids_2)
FSH_data_3 = FSH_data_3 %>% filter(SWANID %in% ids_3)
FSH_data_4 = FSH_data_4 %>% filter(SWANID %in% ids_4)
FSH_data_5 = FSH_data_5 %>% filter(SWANID %in% ids_5)
FSH_data_6 = FSH_data_6 %>% filter(SWANID %in% ids_6)

FSH_data_2 = surv_data_frame_class %>% filter(SWANID %in% sf2$SWANID) %>%
  group_by(SWANID) %>%
  dplyr::select(SWANID,AGE_mod,FSH,E2AVE) %>% distinct() %>% filter(!is.na(FSH),!is.na(E2AVE))
  reframe()

#### Survival


# Trans-specific
cox_class_1 <- coxph(Surv(Tstart, Tstop, status) ~ cluster(SWANID),
                     data = sf1 %>% group_by(SWANID) %>% reframe(Tstart = min(Tstart),
                                                                 Tstop = max(Tstop),
                                                                 status = max(status)), 
                     method = "breslow", 
                     x = TRUE, 
                     model = TRUE)

cox_class_2 <- coxph(Surv(Tstart, Tstop, status) ~ prob1 + prob2 + prob4 + prob5 + E2AVE+ cluster(SWANID),
                     data = sf2, 
                     method = "breslow", 
                     x = TRUE, 
                     model = TRUE)

cox_class_3 <- coxph(Surv(Tstart, Tstop, status) ~ prob1 + prob2 + prob4 + prob5 + E2AVE + cluster(SWANID),
                     data = sf3, 
                     method = "breslow", 
                     x = TRUE, 
                     model = TRUE)

cox_class_4 <- coxph(Surv(Tstart, Tstop, status) ~ prob1 + prob2 + prob4 + prob5 + E2AVE + VMS + SLEEP+ cluster(SWANID),
                     data = sf4, 
                     method = "breslow", 
                     x = TRUE, 
                     model = TRUE)

cox_class_5 <- coxph(Surv(Tstart, Tstop, status) ~ prob1 + prob2 + prob4 + prob5 + E2AVE + VMS + SLEEP+ cluster(SWANID),
                     data = sf5, 
                     method = "breslow", 
                     x = TRUE, 
                     model = TRUE)

cox_class_6 <- coxph(Surv(Tstart, Tstop, status) ~ prob1 + prob2 + prob4 + prob5 + E2AVE + VMS + SLEEP + cluster(SWANID),
                     data = sf6, 
                     method = "breslow", 
                     x = TRUE, 
                     model = TRUE)


#### Longitudinal

FSH_lm1 = lme(fixed = FSH ~ ns(AGE_mod,3),
             random = ~ AGE_mod | SWANID,method = "ML",
             data = FSH_data_1)

FSH_lm2 = lme(fixed = FSH ~ ns(AGE_mod,3),
              random = ~ AGE_mod | SWANID,control = lmeControl(opt='optim'),
              data = FSH_data_2)

FSH_lm3 = lme(fixed = FSH ~ ns(AGE_mod,3),
              random = ~ AGE_mod | SWANID,
              data = FSH_data_3)

FSH_lm4 = lme(fixed = FSH ~ ns(AGE_mod,3),
              random = ~ AGE_mod | SWANID,
              data = FSH_data_4)

FSH_lm5 = lme(fixed = FSH ~ ns(AGE_mod,3),
              random = ~ AGE_mod | SWANID,
              data = FSH_data_5)

FSH_lm6 = lme(fixed = FSH ~ ns(AGE_mod,3),
              random = ~ AGE_mod | SWANID,
              data = FSH_data_6)


### Joint model
jm_1 <- jointModel(FSH_lm1,
                   cox_class_1,
                   timeVar = "AGE_mod",
                   parameterization = "value",
                   CompRisk = F)

summary(jm_1)

jm_2 <- jointModel(FSH_lm2,
                   cox_class_2,
                   timeVar = "AGE_mod",
                   parameterization = "value",
                   CompRisk = F)

jm_3 <- jointModel(FSH_lm3,
                   cox_class_3,
                   timeVar = "AGE_mod",
                   parameterization = "value",
                   CompRisk = F)



FSH_lm = lme(fixed = FSH ~ ns(AGE_mod,3),
             random = ~ AGE_mod | SWANID,
             data = mp_class_m5 %>% filter(!is.na(FSH)))

cox_class <- coxph(Surv(Tstart, Tstop, status) ~ cluster(SWANID),
                   data = surv_data_frame_class, 
                   method = "breslow", 
                   x = TRUE, 
                   model = TRUE)

jm_class <- jointModel(FSH_lm,
                       cox_class,
                   timeVar = "AGE_mod",
                   parameterization = "value",
                   CompRisk = F)
