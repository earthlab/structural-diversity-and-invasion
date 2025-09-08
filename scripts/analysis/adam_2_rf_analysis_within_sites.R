# trying to do a good random forest analysis - within sites
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

# dnew <- read_csv('output/insitu_and_lidar_all_in_one_final_5_31_2023.csv')

d <- read_csv('data/data_with_climate_norms.csv') |>
  mutate(i_cat = as.factor(i_cat)) |>
  dplyr::select(-days_since_change)


# get rid of some sites


# listing out the predictors
predlist <- c("mean.max.canopy.ht","rumple",  "map", "mat",# "days_since_change",
                  "deepgap.fraction", "latitude", "tmax",
                  "tmin", "vertCV","entropy","GFP","VAI", "nspp_native", "shannon_native",
                  "folded_aspect","minElev", "slope")  |> str_to_lower()


resp <- c("rel_cover_exotic","rer", 'i_cat')
# higly correlated (r>80) vars: ,"q25.aop ","VCI.AOP.aop ""max.canopy.ht.aop", "cover.fraction.aop","vert.sd.aop","top.rugosity.aop",


# big random forest function ===================================================


drf <- d |>
  na.omit()
  # dplyr::select(-days_since_change)


mod_form_pi <- as.formula(paste("i_cat ~", paste(predlist, collapse = " + ")))
mod_form_ir <- as.formula(paste("rer ~", paste(predlist, collapse = " + ")))
mod_form_ii <- as.formula(paste("rel_cover_exotic ~", paste(predlist, collapse = " + ")))


# functions ====================================================================
get_impz <- function(rs_obj){
  df_imp <- rs_obj$imps_pi |>  bind_rows() |> mutate(respo = 'A. P(Invasion)') |>
    bind_rows(rs_obj$imps_ii |> bind_rows() |> mutate(respo = 'C. Rel. Cover Non-Natives')) |>
    bind_rows(rs_obj$imps_ir |> bind_rows() |> mutate(respo = 'B. Rel. Richness Non-Natives'))
}

holdout_results_pi <- function(splits, mod_form){
  rf_fit <-
    rand_forest(trees = 5000) %>%
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
    rand_forest(trees = 5000) %>%
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
    rand_forest(trees = 5000) %>% 
    set_engine("ranger", importance = 'permutation') %>% 
    set_mode("classification")
  rf_fit <- 
    rf_mod %>% 
    fit(mod_form, data = analysis(splits))
}

fit_mods_r <- function(splits, mod_form){
  rf_mod <- 
    rand_forest(trees = 5000) %>% 
    set_engine("ranger", importance = 'permutation') %>% 
    set_mode("regression")
  rf_fit <- 
    rf_mod %>% 
    fit(mod_form, data = analysis(splits))
}

get_imps <- function(fit) {
  as_tibble(fit$fit$variable.importance, rownames = 'var')
}

get_r2z <- function(rs_obj){
  bind_rows(
    rs_obj$results_ii |> bind_rows() |> dplyr::select(est = .estimate) |> mutate(respo = "C. RANN", metric = 'R2'),
    rs_obj$results_ir |> bind_rows() |> dplyr::select(est = .estimate) |> mutate(respo = "B. RRNN", metric = 'R2'),
    rs_obj$results_pi |> bind_rows() |> dplyr::select(est = .estimate) |> mutate(respo = "A. P(I)", metric = 'Acc.'))
}
# doing the stuff ==============================================================

by_site_list <- list()

sytes <- unique(drf$site)
# sytes <- unique(drf$domainID)

c <- 1
for(i in 1:length(sytes)){
  
  syte <- sytes[i]
  dddd <- drf |> dplyr::filter(site == syte)
  print(nrow(dddd))
  if(nrow(dddd) > 30){
    rs_obj <- group_vfold_cv(dddd, group = plotID, v=10)
    
    rs_obj$fits_pi <- map(rs_obj$splits, fit_mods_c, mod_form_pi)
    rs_obj$fits_ii <- map(rs_obj$splits, fit_mods_r, mod_form_ii)
    rs_obj$fits_ir <- map(rs_obj$splits, fit_mods_r, mod_form_ir)
    
    
    rs_obj$imps_pi <- map(rs_obj$fits_pi, get_imps)
    rs_obj$imps_ii <- map(rs_obj$fits_ii, get_imps)
    rs_obj$imps_ir <- map(rs_obj$fits_ir, get_imps)
    
    rs_obj$results_pi <- map(rs_obj$splits, holdout_results_pi, mod_form_pi)
    rs_obj$results_ii <- map(rs_obj$splits, holdout_results_r, mod_form_ii)
    rs_obj$results_ir <- map(rs_obj$splits, holdout_results_r, mod_form_ir)
    
    by_site_list[[c]] <- rs_obj
    c <- c+1
  }else{print('skipped')}
}

save(rs_obj, file = "within_site_rf.rda")

# importance plot ==============================================================

pfw <- lapply(by_site_list, get_impz) |>
  bind_rows()|>
  group_by(respo, var) |>
  summarise(mean = mean(value),
            sd = sd(value)) |>
  ungroup() |>
  mutate(mean = ifelse(mean <0, 0, mean)) |>
  mutate(cat = lut_vars[var],
         cat = forcats::fct_relevel(cat, c("Climate", "Site", "Structure", "Community")),
         var = ifelse(var == "folded_aspect", 'aspect', var),
         var = ifelse(var == "mean.max.canopy.ht", 'max_ht', var),
         var = ifelse(var == "deepgap.fraction", 'dgf', var),
         var = ifelse(var == 'nspp_native', "N. Rich", var),
         var = ifelse(var == 'shannon_native', "N. Div", var),
         var = ifelse(var == 'days_since_change', "TSC", var),
         var = ifelse(nchar(var) == 3, str_to_upper(var), var),
         var = fct_reorder2(var, mean, cat))|>
  dplyr::mutate(respo = str_replace_all(respo, "A\\.", "E."),
                respo = str_replace_all(respo, "B\\.", "F."),
                respo = str_replace_all(respo, "C\\.", "G.") |>
                  str_replace_all("Abundance", "Cover")) |>
  ggplot(aes(y=var)) +
  geom_bar(aes(x=mean, fill = cat), stat= 'identity' , color = 'black') +
  facet_wrap(~respo, scales = 'free_x', nrow =1, ncol=3) +
  ggtitle("Within Sites") +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.background = element_rect(fill = NA),
        legend.justification = c(1,0))

ggsave("output/within_site_importance.png", width = 7, height = 7)

# r2 plot ======================================================================

p_r2w <- lapply(by_site_list, get_r2z) |>
  bind_rows() |>
  na.omit()|>
  mutate(panel = 'H. Accuracy') |>
  ggplot(aes(x=respo, y=est, fill = metric)) +
  geom_boxplot() +
  facet_wrap(~panel) +
  ggtitle("_") +
  theme_bw() + 
  ylab("Estimate") +
  scale_y_continuous(limits = c(0,1)) +
  theme(axis.title = element_blank(),
        title = element_text(color = 'white'),
        legend.background = element_rect(fill=NA),
        legend.position = 'none',
        legend.justification = c(0,0)); p_r2w


ggarrange(pfw, p_r2w, ncol = 2, nrow = 1, widths = c(3,1)) 
ggsave("output/new_imp_fig_bottom.png", width =9, height = 4, bg = 'white')

# getting the pdps =============================================================

rf_fit <- list()
pdplist <- list()
cc <- 1

for(bs in 1:length(by_site_list)){
  
  oooo <- by_site_list[[bs]]
  for(i in 1:nrow(oooo)){
    rf_fit$pi <- oooo$fits_pi[[i]]
    rf_fit$ir <- oooo$fits_ir[[i]]
    rf_fit$ii <- oooo$fits_ii[[i]]
    
    splits <- oooo$splits[[i]]
    ssss <- analysis(splits) |> pull(site) |> unique()
    for(j in 1:length(rf_fit)){
      preds <- predlist
      for(k in 1:length(preds)){
          pdplist[[cc]] <- 
            pdp::partial(rf_fit[[j]],
                         approx = FALSE,
                         grid.resolution = 30,
                         train = analysis(splits),
                         pred.var = preds[k]) |>
            mutate(rep = as.character(i),
                   name = preds[k],
                   response = names(rf_fit)[j],
                   site = ssss) |>
            dplyr::rename(value = 1)
        cc <- cc +1
    print(paste(preds[k], names(rf_fit)[j], i, ssss))
      }
    }
  }
}

# turning the result into a df =================================================
lapply(pdplist, function(x)dplyr::mutate(x, value = as.numeric(value))) |>
  bind_rows() |>
  write_csv('data/pdp_df_new_bysite.csv')


pdp_df <- read_csv('data/pdp_df_new_bysite.csv') |> 
  mutate(response = case_when(
    response == 'ii' ~ 'Rel. Cover NN',
    response == 'ir' ~ 'Rel. Richness NN',
    response == 'pi' ~ "P(Invasion)"),
    type = lut_vars[name])

# plotting the pdps ============================================================

pdp_df |>
  group_by(site, rep, name, response) |>
  mutate(value = scale(value),
         yhat = scale (yhat)) |>
  na.omit() |>
  filter(name %in% c('rumple', 'entropy', 'nspp_native', 'mean.max.canopy.ht')) |>
  mutate(name = case_when(
    name == 'rumple' ~ 'Rumple',
    name == 'entropy' ~ 'Entropy',
    name == 'nspp_native' ~ "Native Richness",
    name == 'mean.max.canopy.ht' ~ "Canopy Height")) |>
  ggplot(aes(x=value, y=yhat)) +
  geom_line(aes(color = name, 
                group = paste(name,rep)),
            alpha = 0.75) +
  # geom_smooth(#method = 'lm',
  #   aes(color = name)) +
  facet_grid(response~site, scales = 'free') +
  scale_color_brewer(palette = "Set1") +
  ylab("Prediction (Scaled Values)") +
  xlab("Scaled Values") +
  theme_bw() +
  theme(axis.text = element_blank(),
        legend.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'bottom',#' c(.95,0.05),
        legend.justification = c(.95,0.05))
ggsave(filename = 'output/within_site_pdps.png', width = 10, height = 5, bg = 'white')

