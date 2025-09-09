# Random forest analsysi: In this scripts we investigate the ability of forst structural diversity to detect the invasion and provide a detailed analsysis on best metrics to detect the invasion at NEON sites.

library(tidyverse)
library(ranger)
# library(randomForest)
library(sf)
library(terra)
library(tidymodels)
library(pdp)
library(ggpubr)
library(ggrepel)
library(topomicro)
library(purrr)


lut_vars <- c("mean.max.canopy.ht"="Structure",
              "rumple"="Structure",
              "map"="Climate", "mat"="Climate",
              "deepgap.fraction"="Structure", 
              "latitude"="Site", "site"="Site", 
              'days_since_change' = 'Site',
              "tmax"="Climate", 
              "tmin"="Climate",
              "vertcv"="Structure",
              "entropy"="Structure",
              "gfp"="Structure",
              "vai"="Structure",
              'nspp_native' = "Community",
              'shannon_native' = 'Community',
              "folded_aspect"="Site",
              "slope" = "Site",
              "nlcd_plot_des_main" = "Site",
              "minelev"="Site")

lut_vars1 <- c("tmax"="C. Tmax", 
               "tmin"="D. Tmin",
               'days_since_change' = "TSC",
               "map"="A. MAP", 
               "mat"="B. MAT",
               "latitude"="E. Latitude", 
               "folded_aspect"="F. Folded Aspect",
               "slope" = "G. Slope",
               "nlcd_plot_des_main" = "H. NLCD",
               "minelev"="I. Elevation",
               "vertcv"="J. VertCV",
               "entropy"="K. Entropy",
               "gfp"="L. GFP",
               "vai"="M. VAI",
               "deepgap.fraction"="N. DeepGap Fraction", 
               "mean.max.canopy.ht"="O. Avg Max Canopy Ht",
               "rumple"="P. Rumple",
               'nspp_native' = "Q. Native Richness",
               'shannon_native' = 'R. Native Diversity')

lut_vars2 <- c("tmax"="C. Tmax", 
               "map"="A. MAP", 
               "mat"="B. MAT",
               "latitude"="E. Latitude", 
               "slope" = "D. Slope",
               "entropy"="F. Entropy",
               "mean.max.canopy.ht"="G. Avg Max Canopy Ht",
               "rumple"="H. Rumple",
               'nspp_native' = "I. Native Richness")

lut_vars3 <- c("entropy"="Entropy",
               "mean.max.canopy.ht"="Canopy Height",
               "rumple"="Rumple",
               'nspp_native' = "Native Richness")

d<- read_csv('data/data_with_climate_norms.csv') |>
  mutate(i_cat = as.factor(i_cat)) |>
  dplyr::select(-days_since_change)


# get rid of some sites


# listing out the predictors
predlist <- c("mean.max.canopy.ht","rumple",  "MAP", "MAT", #"days_since_change",
                  "deepgap.fraction", "latitude", "tmax",
                  "tmin", "vertCV","entropy","GFP","VAI", "nspp_native", "shannon_native",
                  "folded_aspect","minElev", "slope") |> str_to_lower() 


resp <- c("rel_cover_exotic","rer", 'i_cat')
# higly correlated (r>80) vars: ,"q25.aop ","VCI.AOP.aop ""max.canopy.ht.aop", "cover.fraction.aop","vert.sd.aop","top.rugosity.aop",


# big random forest function ===================================================


drf <- d |>
  na.omit()
  # dplyr::select(-days_since_change)


mod_form_pi <- as.formula(paste("i_cat ~", paste(predlist, collapse = " + ")))
mod_form_ir <- as.formula(paste("rer ~", paste(predlist, collapse = " + ")))
mod_form_ii <- as.formula(paste("rel_cover_exotic ~", paste(predlist, collapse = " + ")))

rs_obj <- group_vfold_cv(drf, group = plotID,v=10 )

# rs_obj$splits

holdout_results_pi <- function(splits, mod_form){
  rf_fit <-
    rand_forest(trees = 1500) %>%
    set_engine("ranger") %>%
    set_mode("classification") %>%
    fit(formula = mod_form, data = analysis(splits))
  
  holdout <- assessment(splits)
  res <- broom::augment(rf_fit, new_data = holdout)
  acc <- yardstick::accuracy(res, truth = as.character(mod_form[2]), .pred_class)
  acc
}

holdout_results_r <- function(splits, mod_form){
  rf_fit <-
    rand_forest(trees = 1500) %>%
    set_engine("ranger") %>%
    set_mode("regression") %>%
    fit(formula = mod_form, data = analysis(splits))
  
  holdout <- assessment(splits)
  res <- broom::augment(rf_fit, new_data = holdout)
  rsq <- yardstick::rsq(res, truth = as.character(mod_form)[2], .pred)
  
  
  # Return the assessment data set with the additional columns
  rsq
}

fit_mods_c <- function(splits, mod_form){
  rf_mod <- 
    rand_forest(trees = 1500) %>% 
    set_engine("ranger", importance = 'permutation') %>% 
    set_mode("classification")
  rf_fit <- 
    rf_mod %>% 
    fit(mod_form, data = analysis(splits))
}

fit_mods_r <- function(splits, mod_form){
  rf_mod <- 
    rand_forest(trees = 1500) %>% 
    set_engine("ranger", importance = 'permutation') %>% 
    set_mode("regression")
  rf_fit <- 
    rf_mod %>% 
    fit(mod_form, data = analysis(splits))
}

get_imps <- function(fit) {
  as_tibble(fit$fit$variable.importance, rownames = 'var')
}

rs_obj$fits_pi <- map(rs_obj$splits, fit_mods_c, mod_form_pi)
rs_obj$fits_ii <- map(rs_obj$splits, fit_mods_r, mod_form_ii)
rs_obj$fits_ir <- map(rs_obj$splits, fit_mods_r, mod_form_ir)


rs_obj$imps_pi <- map(rs_obj$fits_pi, get_imps)
rs_obj$imps_ii <- map(rs_obj$fits_ii, get_imps)
rs_obj$imps_ir <- map(rs_obj$fits_ir, get_imps)

rs_obj$results_pi <- map(rs_obj$splits, holdout_results_pi, mod_form_pi)
rs_obj$results_ii <- map(rs_obj$splits, holdout_results_r, mod_form_ii)
rs_obj$results_ir <- map(rs_obj$splits, holdout_results_r, mod_form_ir)

# save(rs_obj, file = 'data/rs_obj.rda')

# load("data/rs_obj.rda")

p_r2 <- bind_rows(
  rs_obj$results_ii |> bind_rows() |> dplyr::select(est = .estimate) |> mutate(respo = "C. RANN", metric = 'R2'),
  rs_obj$results_ir |> bind_rows() |> dplyr::select(est = .estimate) |> mutate(respo = "B. RRNN", metric = 'R2'),
  # rs_obj$results_nr |> bind_rows() |> dplyr::select(est = .estimate) |> mutate(respo = "D. NR", metric = 'R2'),
  # rs_obj$results_sn |> bind_rows() |> dplyr::select(est = .estimate) |> mutate(respo = "E. ND", metric = 'R2'),
  rs_obj$results_pi |> bind_rows() |> dplyr::select(est = .estimate) |> mutate(respo = "A. P(I)", metric = 'Acc.')) |>
  mutate(panel = 'D. Accuracy') |>
  ggplot(aes(x=respo, y=est, fill = metric)) +
  geom_boxplot() +
  facet_wrap(~panel) +
  ggtitle("_") +
  theme_bw() + 
  ylab("Estimate") +
  scale_y_continuous(limits = c(0,1)) +
  theme(axis.title = element_blank(), title = element_text(color = 'white'),
        legend.background = element_rect(fill=NA),
        legend.position = c(0,0),
        legend.justification = c(0,0)); p_r2
# ggsave('output/r2_boxplot.png', width = 4, height = 4, bg = 'white')

df_imp <- rs_obj$imps_pi |>  bind_rows() |> mutate(respo = 'A. P(Invasion)') |>
  bind_rows(rs_obj$imps_ii |> bind_rows() |> mutate(respo = 'C. Rel. Abundance Non-Natives')) |>
  bind_rows(rs_obj$imps_ir |> bind_rows() |> mutate(respo = 'B. Rel. Richness Non-Natives')) |>
  # bind_rows(rs_obj$imps_sn |> bind_rows() |> mutate(respo = 'E. Native Diversity')) |>
  # bind_rows(rs_obj$imps_nr |> bind_rows() |> mutate(respo = 'D. Native Richness')) |>
  group_by(respo, var) |>
  summarise(mean = mean(value),
            sd = sd(value)) |>
  ungroup() |>
  mutate(cat = lut_vars[var],
         cat = forcats::fct_relevel(cat, c("Climate", "Site", "Structure", "Community")),
         var = ifelse(var == "folded_aspect", 'aspect', var),
         var = ifelse(var == "mean.max.canopy.ht", 'max_ht', var),
         var = ifelse(var == "deepgap.fraction", 'dgf', var),
         var = ifelse(var == 'nspp_native', "N. Rich", var),
         var = ifelse(var == 'days_since_change', 'TSC', var),
         var = ifelse(var == 'shannon_native', "N. Div", var),
         var = ifelse(nchar(var) == 3, str_to_upper(var), var),
         var = fct_reorder2(var, mean, cat))

pf <- df_imp |>
  ggplot(aes(x=mean, y=var, fill = cat)) +
  geom_bar(stat = 'identity', color = 'black') +
  facet_wrap(~respo, scales = 'free_x') +
  scale_fill_brewer(palette = "Dark2") +
  ggtitle("Among Sites") +
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(.99,0.01),
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA),
        legend.justification = c(1,0))


ggarrange(pf, p_r2, ncol = 2, nrow = 1, widths = c(3,1)) 
ggsave("output/among_sites_imp_fig.png", width =10, height = 5, bg = 'white')


# still need to redo pdps

# pdps - doing a loop

load('data/rs_obj.rda')
rf_fit <- list()
pdplist <- list()
cc <- 1
for(i in 1:nrow(rs_obj)){
  rf_fit$pi <- rs_obj$fits_pi[[i]]
  rf_fit$ir <- rs_obj$fits_ir[[i]]
  rf_fit$ii <- rs_obj$fits_ii[[i]]
  # rf_fit$nr <- rs_obj$fits_nr[[i]]
  # rf_fit$sn <- rs_obj$fits_sn[[i]]
  
  splits <- rs_obj$splits[[i]]
  for(j in 1:length(rf_fit)){
    preds <- rf_fit[[j]]$fit$variable.importance |> names()
    for(k in 1:length(preds)){
      if(preds[k] != "nlcd_plot_des_main"){
        pdplist[[cc]] <- 
          pdp::partial(rf_fit[[j]],
                       approx = TRUE,
                       grid.resolution = 30,
                       train = analysis(splits),
                       pred.var = preds[k]) |>
          mutate(rep = as.character(i),
                 name = preds[k],
                 response = names(rf_fit)[j]) |>
          dplyr::rename(value = 1)
      cc <- cc +1
      }else{print("skip")}
  print(paste(preds[k], names(rf_fit)[j], i))
}}}

lapply(pdplist, function(x)dplyr::mutate(x, value = as.numeric(value))) |>
  bind_rows() |>
  write_csv('data/pdp_df_new.csv')

lapply(pdplist, function(x)dplyr::mutate(x, value = as.numeric(value))) |>
  bind_rows() |>
  # filter(response == 'pi') |>
  # group_by(response, name, value) |>
  # summarise(yhat = mean(yhat)) |>
  # ungroup() |>
  group_by(response) |>
  mutate(yhat= scale(yhat) |> as.numeric()) |>
  ungroup() |>
  ggplot(aes(x=value, y=yhat,# group = rep, 
             color = response)) +
  # geom_line(aes(#group=rep
  #               )) +
  geom_smooth(se=F) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~name, scales = 'free')

pdp_df <- lapply(pdplist, function(x)dplyr::mutate(x, value = as.numeric(value))) |>
  bind_rows() |>
  mutate(response = case_when(
    response == 'ii' ~ 'Invasion Impact',
    response == 'ir' ~ 'Invasion Rate',
    response == 'pi' ~ "P(Invasion)"),
    type = lut_vars[name])

# rumple entropy height nspp native ============================================

pdp_df |>
  filter(name %in% c("rumple", "entropy", 'mean.max.canopy.ht', 'nspp_native')) |>
  mutate(name = str_replace_all(name, "rumple", "Rumple") |>
           str_replace_all("entropy", "Entropy") |>
           str_replace_all("mean.max.canopy.ht", "Canopy Height") |>
           str_replace_all("nspp_native", "Native Richness") |>
           forcats::fct_relevel("Native Richness", "Canopy Height", "Entropy", 'Rumple'),
         response = response |> 
           str_replace_all("Invasion Rate", "Rel. Richness NN (%)") |> 
           str_replace_all("Invasion Impact", "Rel. Cover NN (%)") |>
           forcats::fct_relevel("P(Invasion)", "Rel. Richness NN (%)", "Rel. Cover NN (%)")) |>
  ggplot(aes(x=value, y=yhat)) +
  geom_line(aes(group=rep), color = 'grey') +
  geom_smooth(se=F, color = 'black', span = .5) +
  # scale_color_brewer(palette = "Set1") +
  facet_grid(response~name, scales = 'free') +
  theme_bw() +
  ggtitle("Among Sites")

  ggsave("output/rf_pdps_among_sites_structure.png", width=8, height = 6, bg='white')



# all variables - for supplement -==============================================
plotlist <- list()
cc <- 1
for(i in unique(pdp_df$response)){
  
  plotlist[[cc]] <- pdp_df |>
    filter(response == i, type %in% c("Climate", "Site")) |>
    ggplot(aes(x=value, y=yhat)) +
    geom_line(aes(group=rep), color = 'grey') +
    geom_smooth(se=F, color = 'black', span = .5) +
    # scale_color_brewer(palette = "Set1") +
    facet_wrap(~name, scales = 'free_x') +
    theme_bw() +
    ylab("Standardized Values") +
    ggtitle(i)
  cc <- cc + 1
  
  ggsave(paste0('output/pdp_', i, "_climate", '.png') |> str_replace_all(" ", "_"), width=6, height = 6, bg='white')
  
  plotlist[[cc]] <- pdp_df |>
    filter(response == i, type %in% c("Community", "Structure")) |>
    ggplot(aes(x=value, y=yhat)) +
    geom_line(aes(group=rep), color = 'grey') +
    geom_smooth(se=F, color = 'black', span = .5) +
    # scale_color_brewer(palette = "Set1") +
    facet_wrap(~name, scales = 'free_x') +
    theme_bw() +
    ylab("Standardized Values") +
    ggtitle(i)
  cc <- cc + 1
  
  ggsave(paste0('output/pdp_', i, "_struct", '.png') |> str_replace_all(" ", "_"), width=6, height = 6, bg='white')
}

# standardize response x variables =============================================
pdp_df <- read_csv('data/pdp_df_new.csv')

rep1_vals <- pdp_df |> filter(rep == 1) |> pull(value)

pdp_df1 <- mutate(pdp_df, value = rep1_vals |> rep(10))

# cols <- c("pink", 'firebrick', "red", 'skyblue', 'darkblue')

pdp_df1 |>
  mutate(#response = case_when(
    # response == 'ii' ~ 'Invasion Impact',
    # response == 'ir' ~ 'Invasion Rate',
    # response == 'pi' ~ "P(Invasion)"),
    type = lut_vars[name],
    name = lut_vars1[name],
    response = forcats::fct_relevel( 
      response, c("P(Invasion)",'Invasion Rate',
                  'Invasion Impact')))|>
  dplyr::filter(type %in% c('Climate', "Site")) |>
  group_by(response, value, name, type) |>
  summarise(yhat= mean(yhat)) |>
  ungroup() |>
  group_by(name, type) |>
  mutate(value = scale(value)) |>
  ggplot(aes(x=value, y=yhat,
             color = name)) +
  geom_line(linewidth = 1) +
  ggsci::scale_color_futurama() +
  # scale_color_brewer(palette = "Set1") +
  # scale_color_manual(values = cols) +
  facet_grid(response~type, scales = 'free') +
  theme_bw() +
  ylab("yhat") +
  theme(#legend.position = c(.75, 0.025),
        legend.title = element_blank(),
        legend.justification = c(1,0),
        axis.title.x  = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave('output/pdp_means_climate_site.png', width = 8, height = 6, bg = 'white')

pdp_df1 |>
  mutate(#response = case_when(
    # response == 'ii' ~ 'Invasion Impact',
    # response == 'ir' ~ 'Invasion Rate',
    # response == 'pi' ~ "P(Invasion)"),
    type = lut_vars[name],
    name = lut_vars1[name],
    response = forcats::fct_relevel( 
      response, c("P(Invasion)",'Invasion Rate',
                  'Invasion Impact')))|>
  dplyr::filter(type %in% c('Community', "Structure")) |>
  group_by(response, value, name, type) |>
  summarise(yhat= mean(yhat)) |>
  ungroup() |>
  group_by(name, type) |>
  mutate(value = scale(value)) |>
  ggplot(aes(x=value, y=yhat,
             color = name)) +
  geom_line(linewidth = 1) +
  ggsci::scale_color_simpsons() +
  # scale_color_brewer(palette = "Set1") +
  # scale_color_manual(values = cols) +
  facet_grid(response~type, scales = 'free') +
  theme_bw() +
  ylab("yhat") +
  theme(#legend.position = c(.75, 0.025),
    legend.title = element_blank(),
    legend.justification = c(1,0),
    axis.title.x  = element_blank(),
    # axis.text.y = element_blank(),
    axis.ticks.y = element_blank())
ggsave('output/pdp_means_community_structure.png', width = 8, height = 6, bg = 'white')

# fewer variables ADAMS FAV ====================================================

pdp_df1 |>
  filter(name %in% c('rumple', 'entropy', 'mean.max.canopy.ht', 
                     'nspp_native')) |>
  mutate(
    name = lut_vars3[name],
    response =  response |>
      str_replace_all('Invasion Rate', "B. Rel. Richness NN (%)") |>
      str_replace_all("Invasion Impact", "C. Rel. Cover NN (%)") |>
      str_replace_all("P\\(Invasion\\)", "A. P\\(Invasion\\)")) |> 
  group_by(response, value, name) |>
  summarise(yhat= mean(yhat)) |>
  ungroup() |>
  group_by(name) |>
  mutate(value = scale(value)) |>
  ggplot(aes(x=value, y=yhat,
             color = name)) +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "Set1") +
  # scale_color_manual(values = cols) +
  facet_wrap(~response, scales = 'free') +
  theme_bw() +
  ggtitle("Among Sites") +
  xlab("Scaled values") +
  ylab("Mean Prediciton") +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        # legend.justification = c(1,0),
        # axis.title.x  = element_blank(),
        #axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave('output/pdp_means_reduced_keyvars.png', width = 9, height = 4, bg = 'white')

#### ==================+++++++++++++++++++++++++++++
# 
# do_rf <- function(drf){
#   cell_split <- group_initial_split(drf, group = siteID)
#   cell_train <- training(cell_split)
#   cell_test  <- testing(cell_split)
#   rf_mod <- 
#     rand_forest(trees = 5000) %>% 
#     set_engine("ranger") %>% 
#     set_mode("regression")
#   rf_fit <- 
#     rf_mod %>% 
#     fit(LDMC ~ ., data = cell_train)
#   rf_training_pred <- 
#     predict(rf_fit, cell_train) %>% 
#     # Add the true outcome data back in
#     bind_cols(cell_train %>% 
#                 dplyr::select(LDMC))
#   train_rsq <- rf_training_pred %>%                # training set predictions
#     yardstick::rsq(truth = LDMC, .pred) |>
#     dplyr::rename(train_rsq = 3)
#   
#   rf_testing_pred <- 
#     predict(rf_fit, cell_test) %>% 
#     bind_cols(cell_test %>% dplyr::select(LDMC))
#   test_rsq <- rf_testing_pred %>%                   # test set predictions
#     yardstick::rsq(truth = LDMC, .pred) |>
#     dplyr::rename(test_rsq = 3) |>
#     left_join(train_rsq)
#   return(test_rsq)
# }
# 
# 
# mods <- list()
# counter <- 1
# result_df <- data.frame(resp=NA, r2=NA, mse=NA, predlist=NA)
# 
# for(i in 1:length(resp)) {
#   for(j in 1:length(predlist)){
#     mods[[counter]]  <- formula(paste(resp[i], " ~ ",
#                                       paste(predlist[[j]], 
#                                             collapse = " + "))) %>%
#       randomForest(data =d, ntree = 1000)
#     result_df[counter,1] <- resp[i]
#     result_df[counter,2] <- mods[[counter]]$rsq %>% median() %>% signif(3)
#     result_df[counter,3] <- mods[[counter]]$mse %>% median() %>% signif(3)
#     result_df[counter,4] <- names(predlist)[j]
#     print(resp[i])
#     counter <- counter + 1
#   }
# }
# result_df %>%
# write_csv("output/rf_accuracy.csv")
# 
# names(mods) <- result_df %>%
#   mutate(name = str_c(resp, "_", predlist)) %>%
#   pull(name)
# # basically, elevation and lat are really good predictors of general diversity 
# # (see plots for turnover and nestedness)
# # lapply(mods, plot)
# 
# 
# # multipanel importance plot ====================
# lut_vars <- c("mean.max.canopy.ht.aop"="Structure",
#               "rumple.aop"="Structure",
#               "MAP"="Climate", "MAT"="Climate",
#               "deepgap.fraction.aop"="Structure", 
#               "latitude"="Site", "site"="Site", 
#               "vpdmx"="Climate", "vpdmn"="Climate",
#               "vertCV.aop"="Structure",
#               "entropy.aop"="Structure",
#               "GFP.AOP.aop"="Structure",
#               "VAI.AOP.aop"="Structure",
#               "folded_aspect"="Site",
#               "slope" = "Site",
#               "NLCD_plot_des_main" = "Site",
#               "minElev"="Site")
# 
# dfs<-list()
# for(i in 1:length(mods)){
#   dfs[[i]] <- importance(mods[[i]]) %>%
#     as_tibble(rownames="variable") %>%
#     mutate(response = names(mods)[i])
# }
# 
# 
# # mods_all <- mods[c(1,4,7,10)]
# mods_all <- mods
# impvars_df<- list()
# for(i in 1:length(names(mods_all))){
#   impvars_df[[i]] <- importance(mods_all[[i]]) %>%
#     as_tibble(rownames="variable") %>%
#     arrange(desc(IncNodePurity)) %>%
#     mutate(rank = rank(IncNodePurity),
#            response = names(mods_all)[i])
# }
# impvars<- bind_rows(impvars_df) %>%
#   filter(variable != "NLCD_plot") %>%
#   pull(variable) %>%unique()
# 
# avg_rank <-
#   impvars_df %>%
#   bind_rows() %>%
#   group_by(variable) %>%
#   summarise(mean_rank_imp = mean(rank)) %>%
#   ungroup() %>%
#   mutate(type = lut_vars[variable]) %>%
#   arrange(type, desc(mean_rank_imp)) %>%
#   mutate(rank = letters[1:16],
#          variable = str_remove_all(variable, ".aop"),
#          variable = str_remove_all(variable, ".AOP"))
# 
# topranks <-
#   impvars_df %>%
#   bind_rows() %>%
#   mutate(type = lut_vars[variable],
#          toprank = ifelse(rank > 11, "top5", "bottom11")) %>%
#   dplyr::select(toprank, variable, response, Importance = IncNodePurity) %>%
#   filter(toprank == "top5")%>%
#   mutate(response = ifelse(response == "shannon_native_all", "E. Native Alpha Diversity", response),
#          response = ifelse(response == "nspp_native_all", "D. Native Richness", response),
#          response = ifelse(response == "i_cat_all", "A. P(Invasion)", response),
#          response = ifelse(response == "rer_all", "B. Invasion Rate", response),
#          response = ifelse(response == "rel_cover_exotic_all", "C. Invasion Impact", response),
#          variable = variable %>% str_remove_all(".aop") %>% 
#            str_remove_all(".AOP") |> str_remove_all("_plot_des_main"))
# 
# 
# barplot <- dfs %>%
#   bind_rows() %>%
#   filter(response %in% c("rel_cover_exotic_all", "rer_all", 'i_cat_all',
#                          "nspp_native_all", "shannon_native_all")) |>
#   mutate(response = ifelse(response == "shannon_native_all", "E. Native Alpha Diversity", response),
#          response = ifelse(response == "nspp_native_all", "D. Native Richness", response),
#          response = ifelse(response == "i_cat_all", "A. P(Invasion)", response),
#          response = ifelse(response == "rer_all", "B. Invasion Rate", response),
#          response = ifelse(response == "rel_cover_exotic_all", "C. Invasion Impact", response)) %>%
#   dplyr::rename(Importance = IncNodePurity) %>%
#   mutate(type = lut_vars[variable],
#          variable = variable %>% str_remove_all(".aop") %>%str_remove_all(".AOP") %>%
#            str_remove_all("_plot_des_main") |>
#            fct_reorder2(Importance ,type),
#          type = factor(type)) %>%
#   ggplot() +
#   geom_bar(stat = "Identity",aes(x = Importance, y=variable, fill = type)) +
#   geom_bar(data = topranks, stat = "Identity", color ="black",
#            fill="transparent",aes(x = Importance, y=variable)) +
#   facet_wrap(~response, scales = "free_x", nrow =2) +
#   scale_fill_brewer(palette = "Dark2", name = "Type") +
#   theme_bw() +
#   theme(axis.title.y = element_blank(),
#         legend.background = element_rect(fill=NA),
#         legend.position = c(.95,0.1),
#         legend.justification = c(1,0),
#         legend.title = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_blank())
# barplot
# 
# ggsave(plot = barplot, filename = "output/rf_fig.png", width = 8.5, height = 7, bg="white")
# 
# 
# mods_all <- mods[c(1,4,7,10)]
# if(!file.exists("data/pdp_df.csv")){
#   df_partial <- list()
#   for(i in 1:length(resp)){
#     pdf<-list()
#     for(j in unique(impvars)) {
#       pdf[[j]] <- partial(mods_all[[i]], j) %>%
#       mutate(variable = j) %>%
#       dplyr::rename(value = j)
#     }
#   
#   
#     df_partial[[i]]<- bind_rows(pdf) %>%
#       mutate(response = resp[i])
#   }
#   
#   
#   
#   pdp_df <- bind_rows(df_partial) %>%
#     group_by(response) %>%
#     mutate(yhat = scale(yhat) %>% as.numeric()) %>%
#     ungroup() %>%
#     mutate(Type = lut_vars[variable],
#            variable = variable %>% str_remove_all(".aop") %>% str_remove_all(".AOP")) %>%
#     mutate(response = ifelse(response == "shannon_native", "Native Alpha Diversity", response),
#            response = ifelse(response == "nspp_native", "Native Richness", response),
#            response = ifelse(response == "rer", "Exotic Relative Richness", response),
#            response = ifelse(response == "rel_cover_exotic", "Exotic Relative Cover", response))
# 
# write_csv(pdp_df, "data/pdp_df.csv")
# }else(pdp_df <- read_csv("data/pdp_df.csv"))

# pdp_plot <- pdp_df %>%
#   left_join(avg_rank) %>%
#   mutate(variable = paste0(rank, ". ", variable)) %>%
#   ggplot(aes(x=value, y=yhat, color = response)) +
#   geom_line(key_glyph = "timeseries", lwd=1) +
#   facet_wrap(~variable, scales = "free", nrow = 4) +
#   scale_color_brewer(palette = "Set1") +
#   theme_bw() +
#   ylab("Standardized Values\n                         Structure                                                       Topography                            Climate") +
#   theme(axis.title.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         legend.position = c(.98, .02),
#         legend.justification = c(1,0),
#         legend.background = element_rect(color = "black"),
#         legend.title = element_blank());pdp_plot
#   
# ggsave(plot = pdp_plot,
#        filename = "output/pdps_site.png", 
#        bg="white", width=9.5, height=9)


# figure out hold-out plot numbers
library(dplyr)
for(i in 1:nrow(rs_obj)){
paste(rs_obj$splits[[i]] |> rsample::assessment() |> pull(plotid) |> unique() |> length(), " - ",
rs_obj$splits[[i]] |> rsample::analysis() |> pull(plotid) |> unique() |> length()) |> print()
}
