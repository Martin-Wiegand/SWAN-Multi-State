library("lcmm")
library("ggplot2")

combined_data = readRDS(rds_file("msm_data"))
mp_long = readRDS(rds_file("mp_long"))
lcmm_data = combined_data %>%
  # filter(SWANID %in% mp_long$SWANID) %>%
  dplyr::select(SWANID,AGE_mod,FSH,E2AVE) %>%
  group_by(SWANID) %>%
  mutate(ID = cur_group_id())

library(lcmm)
library(mstate)
library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(splines)

truncator <- function(x,data = lcmm_data){
  data %>% 
    as_tibble() %>%
    group_by(SWANID) %>%
    filter(n() >= x) %>%
    filter(1:n() <= x) %>%
    mutate(trunc = x) %>%
    return()
}

trunc_long = 
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

trunc_long = trunc_long %>% arrange(SWANID,AGE_mod,trunc)

# Init
fit_FSH_1 <- hlme(
  fixed = FSH ~ ns(AGE_mod,3),     # population quadratic; adapt as needed
  random = ~ AGE_mod,                     # random intercept + slope
  subject = "ID",
  ng = 1,
  data = lcmm_data                           # allow class-specific residual variance (optional)
)

saveRDS(fit_FSH_1,rds_file("fit_FSH_1"))

# Groups 2-5
fit_FSH_2<-hlme(FSH ~ ns(AGE_mod,3),
          mixture=~ ns(AGE_mod,3),
          random=~AGE_mod,
          subject='ID', 
          ng=2,
          data=lcmm_data,
          B=fit_FSH_1) 

saveRDS(fit_FSH_2,rds_file("fit_FSH_2"))

fit_FSH_3<-hlme(FSH ~ ns(AGE_mod,3),
                mixture = ~ ns(AGE_mod,3),
                random = ~AGE_mod,
                subject='ID', 
                ng=3,
                data=lcmm_data,
                B=fit_FSH_1) 

saveRDS(fit_FSH_3,rds_file("fit_FSH_3"))

fit_FSH_4<-hlme(FSH ~ ns(AGE_mod,3),
                mixture=~ ns(AGE_mod,3),
                random=~ AGE_mod,
                subject='ID', 
                ng=4,
                data=lcmm_data,
                B=fit_FSH_1) 

saveRDS(fit_FSH_4,rds_file("fit_FSH_4"))

fit_FSH_5<-hlme(FSH ~ ns(AGE_mod,3),
                mixture=~ ns(AGE_mod,3),
                random=~AGE_mod,
                subject='ID',
                ng=5,
                data= lcmm_data,
                B=fit_FSH_1)

saveRDS(fit_FSH_5,rds_file("fit_FSH_5"))

fit_FSH_6<-hlme(FSH ~ ns(AGE_mod,3),
                mixture=~ ns(AGE_mod,3),
                random=~AGE_mod,
                subject='ID',
                ng=6,
                data= lcmm_data,
                B=fit_FSH_1)

saveRDS(fit_FSH_6,rds_file("fit_FSH_6"))

fit_FSH_7<-hlme(FSH ~ ns(AGE_mod,3),
                mixture=~ ns(AGE_mod,3),
                random=~AGE_mod,
                subject='ID',
                ng=7,
                data= lcmm_data,
                B=fit_FSH_1)

saveRDS(fit_FSH_7,rds_file("fit_FSH_7"))


# ESTRADIOL

fit_E2_1 <- hlme(
  fixed = E2AVE ~ ns(AGE_mod,3),     # population quadratic; adapt as needed
  random = ~ AGE_mod,                     # random intercept + slope
  subject = "ID",
  ng = 1,
  data = lcmm_data                           # allow class-specific residual variance (optional)
)

saveRDS(fit_E2_1,rds_file("fit_E2_1"))

fit_E2_2<-hlme(E2AVE ~ ns(AGE_mod,3),
               mixture=~ ns(AGE_mod,3),
               random=~AGE_mod,
               subject='ID', 
               ng=2,
               data=lcmm_data,
               B= random(fit_E2_1)) 

saveRDS(fit_E2_2,rds_file("fit_E2_2"))

fit_E2_3<-hlme(E2AVE ~ ns(AGE_mod,3),
                mixture = ~ ns(AGE_mod,3),
                random = ~AGE_mod,
                subject='ID', 
                ng=3,
                data=lcmm_data,
                B=fit_E2_1) 

saveRDS(fit_E2_3,rds_file("fit_E2_3"))

fit_E2_4<-hlme(E2AVE ~ ns(AGE_mod,3),
                mixture=~ ns(AGE_mod,3),
                random=~ AGE_mod,
                subject='ID', 
                ng=4,
                data=lcmm_data,
                B=fit_E2_1) 

saveRDS(fit_E2_4,rds_file("fit_E2_4"))

fit_E2_5<-hlme(FSH ~ ns(AGE_mod,3),
                mixture=~ ns(AGE_mod,3),
                random=~AGE_mod,
                subject='ID',
                ng=5,
                data= lcmm_data,
                B=fit_E2_1)

saveRDS(fit_E2_5,rds_file("fit_E2_5"))

fit_E2_6<-hlme(FSH ~ ns(AGE_mod,3),
               mixture=~ ns(AGE_mod,3),
               random=~AGE_mod,
               subject='ID',
               ng=6,
               data= lcmm_data,
               B=fit_E2_1)

saveRDS(fit_E2_6,rds_file("fit_E2_6"))

fit_E2_7<-hlme(FSH ~ ns(AGE_mod,3),
               mixture=~ ns(AGE_mod,3),
               random=~AGE_mod,
               subject='ID',
               ng=7,
               data= lcmm_data,
               B=fit_E2_1)

saveRDS(fit_E2_7,rds_file("fit_E2_7"))













# pred_E2_5  = predictY(fit_E2_5,newdata = tibble(AGE_mod = seq(0,18,length.out = 100)))
# pred_E2_5_tibble = bind_cols(tibble(times = as.numeric(unlist(pred_E2_5$times))),as.matrix(pred_E2_5$pred))
# names(pred_E2_5_tibble) = c("times",paste0("C",1:5))
# pred_E2_5_long = pivot_longer(pred_E2_5_tibble,cols = C1:C5,names_to = "class",values_to = "values")
# 
# ggplot(pred_E2_5_long) +
#   geom_line(aes(x = times,y = values,colour = class))+
#   coord_cartesian(ylim = c(-5,30)) +
#   theme_minimal() +
#   scale_colour_manual(values = c("darkgreen","darkblue","orange","firebrick","purple"))


ggplot(tibble(x = c(1,2,3,4,5,6,7),
              E_AIC = c(fit_E2_1$AIC,fit_E2_2$AIC,fit_E2_3$AIC,fit_E2_4$AIC,fit_E2_5$AIC,fit_E2_6$AIC,fit_E2_7$AIC),
              E_BIC = c(fit_E2_1$BIC,fit_E2_2$BIC,fit_E2_3$BIC,fit_E2_4$BIC,fit_E2_5$BIC,fit_E2_6$BIC,fit_E2_7$BIC),
              F_AIC = c(fit_FSH_1$AIC,fit_FSH_2$AIC,fit_FSH_3$AIC,fit_FSH_4$AIC,fit_FSH_5$AIC,fit_FSH_6$AIC,fit_FSH_7$AIC),
              F_BIC = c(fit_FSH_1$BIC,fit_FSH_2$BIC,fit_FSH_3$BIC,fit_FSH_4$BIC,fit_FSH_5$BIC,fit_FSH_6$BIC,fit_FSH_7$BIC))) +
  geom_line(aes(x = x,y = E_AIC),colour = "firebrick")+
  geom_line(aes(x = x,y = E_BIC),colour = "firebrick",linetype = "dashed") +
  geom_line(aes(x = x,y = F_AIC),colour = "darkblue")+
  geom_line(aes(x = x,y = F_BIC),colour = "darkblue",linetype = "dashed") +
  theme_minimal() +
  labs(x = "Number of groups",
       y = "AIC/ BIC") +
  scale_x_continuous(minor_breaks = NULL,
                     breaks = 0:7)+
  scale_y_continuous(minor_breaks = NULL) +
  coord_cartesian(xlim = c(2,7))



### Group trajectory plots
ggplot(mp_long %>%
         left_join(fit_FSH_5$pprob %>% dplyr::select(SWANID,class),by = "SWANID") %>%
         group_by(class) %>%
         mutate(class = paste0(class," n = ",length(unique(SWANID))," (",round(length(unique(SWANID))/(mp_long %>% group_by(SWANID) %>% n_groups())*100,1),"%)")) %>%
         mutate(AGE_mod = AGE_mod + 42)) +
  geom_line(aes(x = AGE_mod,y = FSH,colour = class,group = SWANID),alpha = 0.2) +
  geom_smooth(aes(x = AGE_mod,y = FSH,colour = class,fill = class)) +
  theme_minimal() +
  labs(x = "Age",
       y = "FSH trajectory",
       colour = "",
       fill = "") +
  theme(legend.position = "bottom")
  
ggsave(plot_file("LC_5groups.png"),width = cm(3),height = cm(2))

ggplot(mp_long %>%
         left_join(fit_FSH_4$pprob %>% dplyr::select(SWANID,class),by = "SWANID") %>%
         group_by(class) %>%
         mutate(class = paste0(class," n = ",length(unique(SWANID))," (",round(length(unique(SWANID))/(mp_long %>% group_by(SWANID) %>% n_groups())*100,1),"%)")) %>%
         mutate(AGE_mod = AGE_mod + 42)) +
  geom_line(aes(x = AGE_mod,y = FSH,colour = class,group = SWANID),alpha = 0.2) +
  geom_smooth(aes(x = AGE_mod,y = FSH,colour = class,fill = class)) +
  theme_minimal() +
  labs(x = "Age",
       y = "FSH trajectory",
       colour = "",
       fill = "") +
  theme(legend.position = "bottom")

ggsave(plot_file("LC_4groups.png"),width = cm(3),height = cm(2))

ggplot(mp_long %>%
         left_join(fit_FSH_3$pprob %>% dplyr::select(SWANID,class),by = "SWANID") %>%
         group_by(class) %>%
         mutate(class = paste0(class," n = ",length(unique(SWANID))," (",round(length(unique(SWANID))/(mp_long %>% group_by(SWANID) %>% n_groups())*100,1),"%)")) %>%
         mutate(AGE_mod = AGE_mod + 42)) +
  geom_line(aes(x = AGE_mod,y = FSH,colour = class,group = SWANID),alpha = 0.2) +
  geom_smooth(aes(x = AGE_mod,y = FSH,colour = class,fill = class)) +
  theme_minimal() +
  labs(x = "Age",
       y = "FSH trajectory",
       colour = "",
       fill = "") +
  theme(legend.position = "bottom")

ggsave(plot_file("LC_3groups.png"),width = cm(3),height = cm(2))


## Group trajec straight
pred3 = predictY(fit_FSH_3,newdata = mp_long,var.time = "AGE_mod")

d3 = 
  tibble(time = pred3$times %>% unlist %>% as.numeric,
         C1 = as.numeric(pred3$pred[,1]),
         C2 = as.numeric(pred3$pred[,2]),
         C3 = as.numeric(pred3$pred[,3]))

ggplot(data = d3) +
  geom_line(aes(x = time,y = C1),colour = "red")+
  geom_line(aes(x = time,y = C2),colour = "darkblue")+
  geom_line(aes(x = time,y = C3),colour = "darkgreen")+
  theme_minimal()



#############
mp_long_class <-  mp_long %>%
  group_by(SWANID) %>%
  mutate(trunc = 1:n()) %>%
  left_join(trunc_pred_3 %>%
              mutate(SWANID = as.numeric(substr(truncid,start = 1,stop = 5)),
                     trunc = as.numeric(substr(truncid,start = 6,stop = nchar(truncid)))) %>%
              dplyr::select(-truncid),
            by = c("SWANID","trunc")) 

saveRDS(mp_long_class,rds_file("mp_long_class"))












####### Multivariate FSH + E2
combined_data = readRDS(rds_file("msm_data"))

combined_data = combined_data %>%
  filter(!is.na(FSH) & !is.na(E2AVE)) %>%
  dplyr::select(SWANID,VISIT,AGE_mod,FSH,E2AVE)

combined_data_long = combined_data %>%
  pivot_longer(cols = c("FSH","E2AVE"),
               names_to = "marker",
               values_to = "value") %>%
  mutate(marker = (marker == "FSH")) %>%
  mutate(AGE_mod = AGE_mod - 42)

lcmm_long = lcmm_data %>%
  pivot_longer(cols = c("FSH","E2AVE"),
               names_to = "marker",
               values_to = "value")

m1 <- multlcmm(
  fixed = value ~ marker + ns(AGE_mod,3)*marker,
  random = ~ AGE_mod,
  subject = "ID",
  ng = 1,                # number of latent classes
  idiag = TRUE,          # diagonal random effect covariance
  link = "linear",
  data = lcmm_long
)

# Two classes
m2 <- multlcmm(
  fixed = value ~ marker + ns(AGE_mod,3)*marker,
  random = ~ AGE_mod,
  mixture = ~ ns(AGE_mod,3)*marker,
  subject = "ID",
  ng = 2,                # number of latent classes
  idiag = TRUE,          # diagonal random effect covariance
  link = "linear",
  data = lcmm_long,
  B = m1
)

m3 <- multlcmm(
  fixed = value ~ marker + ns(AGE_mod,3)*marker,
  random = ~ AGE_mod,
  mixture = ~ ns(AGE_mod,3)*marker,
  subject = "ID",
  ng = 3,                # number of latent classes
  idiag = TRUE,          # diagonal random effect covariance
  link = "linear",
  data = lcmm_long,
  B = m1
)

m4 <- multlcmm(
  fixed = value ~ marker + ns(AGE_mod,3)*marker,
  random = ~ AGE_mod,
  mixture = ~ ns(AGE_mod,3)*marker,
  subject = "ID",
  ng = 4,                # number of latent classes
  idiag = TRUE,          # diagonal random effect covariance
  link = "linear",
  data = lcmm_long,
  B = m1
)

m5 <- multlcmm(
  fixed = value ~ marker + ns(AGE_mod,3)*marker,
  random = ~ AGE_mod,
  mixture = ~ ns(AGE_mod,3)*marker,
  subject = "ID",
  ng = 5,                # number of latent classes
  idiag = TRUE,          # diagonal random effect covariance
  link = "linear",
  data = lcmm_long,
  B = m1
)

m6 <- multlcmm(
  fixed = value ~ marker + ns(AGE_mod,3)*marker,
  random = ~ AGE_mod,
  mixture = ~ ns(AGE_mod,3)*marker,
  subject = "ID",
  ng = 6,                # number of latent classes
  idiag = TRUE,          # diagonal random effect covariance
  link = "linear",
  data = lcmm_long,
  B = m1
)

m7 <- multlcmm(
  fixed = value ~ marker + ns(AGE_mod,3)*marker,
  random = ~ AGE_mod,
  mixture = ~ ns(AGE_mod,3)*marker,
  subject = "ID",
  ng = 7,                # number of latent classes
  idiag = TRUE,          # diagonal random effect covariance
  link = "linear",
  data = lcmm_long,
  B = m1
)


saveRDS(m1,rds_file("m1_unscaled"))
saveRDS(m2,rds_file("m2_unscaled"))
saveRDS(m3,rds_file("m3_unscaled"))
saveRDS(m4,rds_file("m4_unscaled"))
saveRDS(m5,rds_file("m5_unscaled"))
saveRDS(m6,rds_file("m6_unscaled"))
saveRDS(m7,rds_file("m7_unscaled"))

m1 = readRDS(rds_file("m1"))
m2 = readRDS(rds_file("m2"))
m3 = readRDS(rds_file("m3"))
m4 = readRDS(rds_file("m4"))
m5 = readRDS(rds_file("m5"))
m6 = readRDS(rds_file("m6"))
m7 = readRDS(rds_file("m7"))
#

ggplot(tibble(x = c(2,3,4,5),
              AIC = c(m2$AIC,m3$AIC,m4$AIC,m5$AIC),
              BIC = c(m2$BIC,m3$BIC,m4$BIC,m5$BIC))) +
  geom_line(aes(x = x,y = AIC),colour = "red")+
  geom_line(aes(x = x,y = BIC),colour = "blue") +
  theme_minimal() +
  labs(x = "Number of groups",
       y = "AIC/ BIC") +
  scale_x_continuous(minor_breaks = NULL,
                     breaks = 0:7)+
  scale_y_continuous(minor_breaks = NULL)

ggsave(plot_file("AIC_multi_all_4groups.svg"),width = cm(3),height = cm(2))

## Group trajec straight
trunc_pivot = trunc_long %>%     
  pivot_longer(cols = c("FSH","E2AVE"),names_to = "marker",values_to = "value")
 
pred_m4 = predictClass(model = m4,
                       newdata = trunc_pivot,
                       subject = "truncid")

lcmm_m4 = lcmm_data %>%
  group_by(SWANID) %>%
  mutate(truncid = paste0(SWANID,1:n())) %>%
  left_join(pred_m4,
            by = "truncid")


lcmm_m4 %>% group_by(SWANID) %>% slice(n()) %>% group_by(class) %>% tally

saveRDS(lcmm_m4,rds_file("lcmm_m4_all"))
# saveRDS(lcmm_m4,rds_file("lcmm_m4"))
# lcmm_m5 = readRDS(rds_file("lcmm_m5"))

ggplot(lcmm_data %>%
         mutate(ID = as.character(ID)) %>%
         left_join(m4$pprob,by = "ID") %>%
         left_join(m4$pprob %>%   
                     group_by(class) %>% 
                     tally %>%
                     mutate(freq = paste0(class,": n =",n," (",round(n/3002,3)*100,"%)")),by = "class") %>%
         filter(!is.na(class))) +
  geom_smooth(aes(x = AGE_mod,y = FSH,group = freq,colour = freq,fill = freq)) +
  geom_smooth(aes(x = AGE_mod,y = E2AVE,group = freq,colour = freq,fill = freq),linetype = "dashed") +
  theme_minimal() +
  theme(legend.position = "NULL") +
  facet_wrap(. ~ freq,ncol = 1) +
  geom_hline(yintercept = 0,linewidth = 0.5) +
  coord_cartesian(ylim = c(0,350)) +
  labs(x = "Age in years",
       y = "FSH (solid) E2 (dashed)",
       fill = "",
       colour = "")

ggsave(plot_file("Multi_FSH_E2_all_4groups.svg"),width = cm(3),height = cm(7))
