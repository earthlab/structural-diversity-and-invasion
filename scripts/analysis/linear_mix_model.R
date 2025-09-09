# glmm analysis 

sapply(c('tidyverse', 'glmmTMB', 'performance', 'splines', 'ggeffects', 'ggrepel', 'broom.mixed', 'ggpubr'),
       library, character.only=T) |> invisible()

# data input

d <- read_csv('data/data_with_climate_norms.csv') |>
  mutate(nlcd = as.factor(NLCD_code_plot)) |>
  dplyr::select(-sd.sd, -vertcv) |>
  dplyr::filter(!is.na(days_since_change))

# # doing a princomp ===========================================================
pca <- d |>
  dplyr::select(vai, vci, vert.sd, q25, q50, gfp, entropy, rumple, vertcv,
                top.rugosity, deepgap.fraction, mean.max.canopy.ht, max.canopy.ht,
                cover.fraction) |>
  na.omit() |>
  prcomp(scale=TRUE)
# 
pca$rotation |>
  as_tibble(rownames = "var") |>
  ggplot(aes(x=PC1, y=PC2)) +
  geom_text_repel(aes(label = var)) +
  xlab("PCA Axis 1: order to disorder") +
  ylab("PCA Axis 2: high 3d surface area?") +
  theme_bw()
ggsave("output/pca_plot.png", bg='white')
d_pca <- pca$x |>
  as_tibble() |>
  dplyr::select(pc1_order = PC1, pc2_sa_vol = PC2)

d1 <- cbind(d, d_pca)

# p(inv) =======================================================================
mpi <- glmmTMB(i_cat ~ nspp_native + rumple + entropy + mean.max.canopy.ht + #days_since_change +#vai +# deepgap.fraction + 
                 #tmax * map + 
                 (1|site/plotID), 
               data = d, family = 'binomial')

summary(mpi); car::Anova(mpi);performance::r2(mpi)
performance::check_model(mpi)

row_pi <- ggpubr::ggarrange(plotlist = list(plot(ggeffects::predict_response(mpi, terms = c('nspp_native [all]'))) +
                                    ylab("P(Invasion)") + theme(title = element_blank()),
                                  plot(ggeffects::predict_response(mpi, terms = c('entropy [all]'))) +
                                    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), title = element_blank()),
                                  plot(ggeffects::predict_response(mpi, terms = c('rumple [all]')))+
                                    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), title = element_blank()),
                                  plot(ggeffects::predict_response(mpi, terms = c('mean.max.canopy.ht [all]')))+
                                    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), title = element_blank())),
                  nrow =1, ncol = 4, widths = c(1, .8,.8,.8))

# Invasion Rate ================================================================
mir <- glmmTMB(rer ~ nspp_native  +  rumple + entropy + mean.max.canopy.ht  + #vai +#deepgap.fraction +# mat+map +#nlcd +
                 (1|site/plotID), 
               data = d,
               ziformula = ~ 1 +
                 (1|site/plotID),
               family = beta_family())

summary(mir);performance::r2(mir); car::Anova(mir)
performance::check_model(mir)
performance::check_overdispersion(mir)
DHARMa::testDispersion(mir)
DHARMa::simulateResiduals(mir, plot = TRUE)


row_ir <- ggpubr::ggarrange(plotlist = list(plot(ggeffects::predict_response(mir, terms = c('nspp_native [all]'))) +
                                              ylab("Rel. Rich. Non-Native") + theme(title = element_blank()),
                                            plot(ggeffects::predict_response(mir, terms = c('entropy [all]'))) +
                                              theme(axis.title.y = element_blank(), axis.text.y = element_blank(), title = element_blank()),
                                            plot(ggeffects::predict_response(mir, terms = c('rumple [all]')))+
                                              theme(axis.title.y = element_blank(), axis.text.y = element_blank(), title = element_blank()),
                                            plot(ggeffects::predict_response(mir, terms = c('mean.max.canopy.ht [all]')))+
                                              theme(axis.title.y = element_blank(), axis.text.y = element_blank(), title = element_blank())),
nrow =1, ncol = 4, widths = c(1, .8,.8,.8))


# Invasion impact ================================================================
mii <- glmmTMB(rel_cover_exotic ~ rumple + entropy + mean.max.canopy.ht +  #vai +#nlcd +
                 nspp_native  + 
                 (1|site/plotID), 
               data = d,
               ziformula = ~ 1 +
                 (1|site/plotID),
               family = beta_family())

summary(mii);performance::r2(mii); car::Anova(mii)
performance::check_model(mii)
performance::check_overdispersion(mii)
DHARMa::testDispersion(mii)
DHARMa::simulateResiduals(mii, plot = TRUE)


row_ii <- ggpubr::ggarrange(plotlist = list(plot(ggeffects::predict_response(mii, terms = c('nspp_native [all]'))) +
                                              ylab("Rel. Ab. Non-Native") + theme(title = element_blank()),
                                            plot(ggeffects::predict_response(mii, terms = c('entropy [all]'))) +
                                              theme(axis.title.y = element_blank(), axis.text.y = element_blank(), title = element_blank()),
                                            plot(ggeffects::predict_response(mii, terms = c('rumple [all]')))+
                                              theme(axis.title.y = element_blank(), axis.text.y = element_blank(), title = element_blank()),
                                            plot(ggeffects::predict_response(mii, terms = c('mean.max.canopy.ht [all]')))+
                                              theme(axis.title.y = element_blank(), axis.text.y = element_blank(), title = element_blank())),
nrow =1, ncol = 4, widths = c(1, .8,.8,.8))


ggarrange(row_pi, row_ir, row_ii, nrow=3, ncol =1)
ggsave(filename = 'output/glmmtmb_plots.png', width =8, height = 6, bg = 'white')


tidy(car::Anova(mii)) |> 
  mutate(resp = "RANN") |>
  bind_rows(tidy(car::Anova(mir)) |> mutate(resp = "RRNN")) |>
  bind_rows(tidy(car::Anova(mpi)) |> mutate(resp = "P(I)")) |>
  mutate(star = ifelse(p.value < 0.01, ".", ""),
         star = ifelse(p.value < 0.05, '*', star),
         star = ifelse(p.value < 0.01, '**', star),
         star = ifelse(p.value < 0.001, '***', star)) |>
  mutate_if(is.numeric, round, 2) |>
  dplyr::select(resp, term, statistic, df, p.value, star) |>
  write_csv('output/glmm_tables.csv')

# old stuff ====================================================================

# library(brms)
# bir <- brm(I((rer/100)) ~ latitude + mat + pc1_order +
#                  vpdmn + nspp_native + (1|siteid/plotid)
#            hu ~ latitude + mat + vpdmn + folded_aspect + shannon_native + nlcd_plot_des_main + 
#              (1|siteid/plotid), 
#                data = d1,
#                family = 'zero_inflated_beta')

# not too bad diagnostics

tidy(mir, conf.int = T) |>
  filter(effect != 'ran_pars', term !='(Intercept)') |>
  mutate(sig = ifelse(conf.low * conf.high > 0, "Significant", 'Not\nSignificant'),
         component = ifelse(component == 'cond', "Invasion Rate, Conditional on Invasion",
         'Effect on Invasion')) |>
  ggplot(aes(y = term))  +
  geom_vline(xintercept = 0, linetype = 2, color = 'grey30') +
  geom_point(aes(x=estimate#, color = sig
                 )) +
  geom_segment(aes(x=conf.low, xend=conf.high#, color = sig
                   ), linewidth = 1) +
  facet_wrap(~component, scales = 'free') +
  theme_bw() +
  ggtitle("Invasion Rate, Hurdle Model")
ggsave("output/hurdle_model_invasion_rate.png", width = 10, height =5, bg = 'white')


emmeans::emmeans(mir, specs = 'nlcd_plot_des_main', )
predict(mir, newdata = d1)

ggeffects::ggpredict(mir, terms = c('vai')) |> plot()
ggeffects::predict_response(mir, terms = c('nlcd_plot_des_main', 'entropy')) |> plot()

# ii
mii <- glmmTMB(rel_cover_exotic ~ pc1_order + latitude + vpdmx + (1|siteid/plotid), 
               data = d1, 
               ziformula = ~ latitude + mat + vpdmx +  nspp_native + pc2_sa_vol + egf +
                 (1|siteid/plotid), 
               family = beta_family())
summary(mii);performance::r2(mii);car::Anova(mii)

performance::check_model(mii)

tidy(mii, conf.int = T) |>
  filter(effect != 'ran_pars', term !='(Intercept)') |>
  mutate(sig = ifelse(conf.low * conf.high > 0, "Significant", 'Not\nSignificant'),
         component = ifelse(component == 'cond', "Invasion Impact, Conditional on Invasion",
                            'Effect on Invasion')) |>
  ggplot(aes(y = term))  +
  geom_vline(xintercept = 0, linetype = 2, color = 'grey30')+
  geom_point(aes(x=estimate#, color = sig
                 )) +
  geom_segment(aes(x=conf.low, xend=conf.high#, color = sig
                   ), linewidth = 1) +
  facet_wrap(~component, scales = 'free') +
  theme_bw() +
  ggtitle("Invasion Impact, Hurdle Model")


# make an effects plot =========================================================
p_df0 <- bind_rows(
  lapply(ggeffects::ggpredict(mpi), as.data.frame) |> 
    bind_rows() |> mutate(response = "P(Invasion)"),
  lapply(ggeffects::ggpredict(mir), as.data.frame) |> 
    bind_rows() |> mutate(response = "Invasion Rate"),
  lapply(ggeffects::ggpredict(mpi), as.data.frame) |> 
    bind_rows() |> mutate(response = "Invasion Impact")
)

p_df <- p_df0 |>
  mutate(type = case_when(
    str_sub(group,1,3) %in% c('mat', 'vpd', 'map') ~ 'Climate',
    group %in% c('minelev', 'slo', 'asp') ~ "Site",
    group %in% c('nspp_native') ~ "Community",
    str_sub(group,1,4) %in% c(c('rump', 'entr', 'mean')) ~ "Structure"),
    # group = case_when(
    #   group == 'cwd_z' ~ "CWD Z: Sample Year",
    #   group == 'def_norm' ~ "CWD: 30-year Normal",
    #   group == 'hli' ~ "Heat Load Index",
    #   group == 'total_cover' ~ "Total Veg Cover",
    #   group == 'tpa' ~ "Trees/ha",
    #   group == 'fwd' ~ "Fine Woody Debris",
    #   group == 'quadratic_mean_diameter' ~ "QMD",
    #   group == 'slope' ~ "Slope",
    #   group == 'ba_m2perha' ~ 'Basal Area (m2/ha)',
    #   group == 'twi' ~ "TWI",
    # #   group == 'def_z_trt' ~ "CWD Z: Year of Treatment"),
    predicted = ifelse(str_sub(response,1,8) %in% c("Invasion"), predicted * 100, predicted), # readjusting the beta 0/1 to percent
    conf.high = ifelse(str_sub(response,1,8) %in% c("Invasion"), conf.high  * 100, conf.high ),
    conf.low = ifelse(str_sub(response,1,8) %in% c("Invasion"), conf.low  * 100, conf.low )
  )


plot_row <- function(p_df, legend = F){
  p <- ggplot(p_df, aes(x=x, y = predicted)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = type)) +
    geom_line(lwd=1) +
    geom_line(aes(y=conf.low), lty = 3, lwd=.5) +
    geom_line(aes(y=conf.high), lty = 3, lwd=.5) +
    facet_wrap(~group, scales = 'free_x', nrow = 1) +
    scale_fill_brewer(palette = "Set2")+
    ylab(p_df$response) +
    theme_bw() +
    theme(axis.title.x = element_blank())
  if(legend){return(p + theme(legend.title = element_blank()))}else(return(p + theme(legend.position = 'none')))
}
ggarrange(
  plot_row(p_df |> filter(response == "P(Invasion)"), T),
  plot_row(p_df |> filter(response == "Invasion Rate"), T),
  plot_row(p_df |> filter(response == "Invasion Impact"), T),
  nrow = 3, common.legend = T, legend = "bottom")

ggsave('output/glmmTMB_plots.png', width = 8, height = 6, bg = 'white')


# side quest =================

ggplot(d, aes(x=year, y = entropy, color = siteid, group=plotid)) +
  geom_line(show.legend = T) 

ggplot(d, aes(x=rer, y=rel_cover_exotic,color = siteid)) +
  geom_point()

