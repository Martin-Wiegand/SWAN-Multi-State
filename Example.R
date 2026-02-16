load(data_file("data.RData"))
ls()

library(mstate) # Please use the version 0.2.7
library(survival)
library(JM)

# get list of R files
runR = list.files(path = main_path)
runR = runR[grepl(".R",runR)]

sapply(X = runR,FUN = function(x){source(code_file(x))})

data_long = as_tibble(data_long)
data_surv = as_tibble(data_surv)

#
# library("ggplot2")

# plot_long <- (ggplot(data_long) +
#                 geom_line(aes(x = times, y = Y, group = id), color = "grey30", alpha = 0.8) +
#                 stat_smooth(aes(x = times, y = Y), method = "loess", size = 0.75) +
#                 theme_bw() +
#                 xlab("Time") +
#                 ylab("Marker value"))
# plot_long


### Linear model
lmeFit <- lme(fixed = Y ~ times,
              data = data_long,
              random = ~ times| id,
              method = "REML",
              control = list(opt = "optim"))


## Multistate subpart
tmat <- matrix(NA, 3, 3)
tmat[1, 2:3] <- 1:2
tmat[2, 3] <- 3
dimnames(tmat) <- list(from = c("State_0", "State_1", "State_2"),
                       to = c("State_0", "State_1", "State_2"))
tmat

# Describe covariates for ms subpart
covs <- "X"

data_mstate <- msprep(time = c(NA, "t_State_1", "t_State_2"),
                      status = c(NA, "State_1", "State_2"),
                      data = data_surv,
                      trans = tmat,
                      keep = covs,
                      id = "id") 

data_mstate <- expand.covs(data_mstate, covs,
                           append = TRUE, longnames = FALSE)
head(data_mstate)

## Perc transitions
events(data_mstate)

# Survival model
coxFit <- coxph(Surv(Tstart, Tstop, status) ~ X + strata(trans),
                data = data_mstate, 
                method = "breslow", 
                x = TRUE, 
                model = TRUE)


## Fixed and random parts that go from long. to survival model
dForm <- list(fixed = ~ times,
              # indFixed = c(2:3 ,5:6),
              random = ~ times)

jm_example <- JMstateModel(lmeObject = lmeFit,
             survObject = coxFit,
             timeVar = "times",
             parameterization = "value",
             method = "spline-PH-GH",
             derivForm = dForm,
             Mstate = TRUE,
             data.Mstate = data_mstate,
             ID.Mstate = "id",
             control = list(GHk = 3, lng.in.kn = 1),
             verbose = TRUE)

summary(jm_example)
