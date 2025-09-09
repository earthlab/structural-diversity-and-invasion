# adam's analysis
library(tidyverse)
library(ranger)
library(randomForest)
library(pdp)
d <- read_csv("data_quality_control/forest_div_quality_control_v1.csv") %>%
  dplyr::mutate(latitude = latitude.x,
                cosign_aspect = cos(aspect),
                slope_aspect = slope*cosign_aspect)
glimpse(d)
names(d)
summary(d)
preds <- c("mean.max.canopy.ht.aop","rumple.aop",
           "deepgap.fraction.aop", "latitude",
          "vertCV.aop","entropy.aop","GFP.AOP.aop","VAI.AOP.aop",
           "slope_aspect","minElev ") 
# higly correlated (r>80) vars: ,"q25.aop ","VCI.AOP.aop ""max.canopy.ht.aop", "cover.fraction.aop","vert.sd.aop","top.rugosity.aop",

# update methods in google doc
# add domain, site -- good random effect for LMM -  spatial regime model for lmm?
# maybe do an NMDS as well
# get climate normals from PRISM, aridity index
#"cover.fraction.aop","vert.sd.aop","top.rugosity.aop",
# super basic random forest analysis

mods <- list()
mods[[1]]  <- formula(paste("rel_cover_exotic ~", paste(preds, collapse = " + "))) %>%
  randomForest(data =d)
mods[[2]]  <- formula(paste("turnover ~", paste(preds, collapse = " + "))) %>%
  randomForest(data =d)
mods[[3]]  <- formula(paste("nestedness~", paste(preds, collapse = " + "))) %>%
  randomForest(data =d)
mods[[4]]  <- formula(paste("cover_exotic~", paste(preds, collapse = " + "))) %>%
  randomForest(data =d)
mods[[5]]  <- formula(paste("shannon_exotic~", paste(preds, collapse = " + "))) %>%
  randomForest(data =d)
mods[[6]]  <- formula(paste("shannon_notexotic~", paste(preds, collapse = " + "))) %>%
  randomForest(data =d)

# basically, elevation and lat are really good predictors of general diversity 
# (see plots for turnover and nestedness)

lapply(mods, varImpPlot)

cor(d %>% dplyr::select_if(is.numeric))

impvars<- importance(mods[[1]]) %>%
  as_tibble(rownames="variable") %>%
  arrange(desc(IncNodePurity)) %>%
  slice(1:6) %>%
  pull(variable)

pdf<-list()
for(i in impvars) pdf[[i]] <- partial(mods[[1]], i) %>%
  mutate(variable = i) %>%
  dplyr::rename(value = i)

bind_rows(pdf) %>%
  ggplot(aes(x=value, y=yhat)) +
  geom_line() +
  facet_wrap(~variable, scales = "free")


# glm analysis
library(lmerTest)
library(performance)
library(brms)

# might rescale variables, change to a non-gaussian distribution
mod <- formula(paste("I(rel_cover_exotic/100) ~", 
                     paste(preds, collapse = " + "),
                     "+ (1|site)")) %>%
  lmer(data =d,)

summary(mod)
performance(mod)
check_model(mod)
