library(here)
library(posterior)
library(cmdstanr)
library(tidyverse)
library(cowplot)
library(bayestestR)

build_density_slice_df = 
  function(post, param, g){
    result = 
      data.frame(
        x = 
          post %>% 
          filter(parameter==param&group_id==g) %>%
          pull(value) %>%
          density() %>%
          .$x,
        y =
          post %>%
          filter(parameter==param & group_id==g) %>%
          pull(value) %>%
          density() %>%
          .$y,
        xmin = 
          post %>%
          filter(parameter == param & group_id == g) %>%
          pull(value) %>%
          #quantile(., 0.025),
          bayestestR::ci(x = ., ci = 0.89, method = 'SPI') %>%
          .$CI_low,
        xmax = 
          post %>%
          filter(parameter == param & group_id == g) %>%
          pull(value) %>%
          #quantile(., 0.975),
          bayestestR::ci(x = ., ci = 0.89, method = 'SPI') %>%
          .$CI_high,
        mort_class = as.character(g)) %>%
      filter(x >= xmin & x <= xmax)
    return(result)
  }

#### litter ################################################################

litter_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'litter_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'alpha[1]', 'alpha[2]', 'alpha[3]',
                                     'rho[1]', 'rho[2]', 'rho[3]',
                                     'sigmaPlot[1]', 'sigmaPlot[2]', 'sigmaPlot[3]',
                                     'kappa[1]', 'kappa[2]', 'kappa[3]')
                     ) %>%
  as_draws_df()

litter_posterior_long = 
  litter_posterior %>%
  select(contains(c('alpha', 'beta', 'rho', 'sigmaPlot', 'kappa'))) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter_full',
               values_to = 'value') %>%
  mutate(group_id = 
           as.integer(gsub(x = parameter_full, pattern = '^.*\\[|\\]$',replacement='')),
         parameter = 
           gsub(x = parameter_full, pattern = '\\[.*$', replacement = '')) %>%
  left_join(data.frame(group_id = 1:3,
                        mort_class = factor(c('Low', 'Med.', 'High'),
                                            levels = c('Low', 'Med.', 'High'))))



litter_beta_plot = 
  ggplot()+
  geom_density(
    data = litter_posterior_long %>% filter(parameter == 'beta'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = litter_posterior_long,
                                  param = 'beta',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = litter_posterior_long,
                                  param = 'beta',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = litter_posterior_long,
                                  param = 'beta',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(litter_posterior_long %>% filter(parameter=='beta') %>% pull(value) %>% min(),
                         litter_posterior_long %>% filter(parameter=='beta') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = dnorm(x = x, mean = 0, sd = 2)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'a) Intercept')

litter_beta_plot

litter_alpha_plot = 
  ggplot()+
  geom_density(
    data = litter_posterior_long %>% filter(parameter == 'alpha'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = litter_posterior_long,
                                  param = 'alpha',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = litter_posterior_long,
                                  param = 'alpha',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = litter_posterior_long,
                                  param = 'alpha',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(litter_posterior_long %>% filter(parameter=='alpha') %>% pull(value) %>% min(),
                         litter_posterior_long %>% filter(parameter=='alpha') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 2, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'b) GP Magnitude')

litter_alpha_plot

litter_rho_plot = 
  ggplot()+
  geom_density(
    data = litter_posterior_long %>% filter(parameter == 'rho'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = litter_posterior_long,
                                  param = 'rho',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = litter_posterior_long,
                                  param = 'rho',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = litter_posterior_long,
                                  param = 'rho',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(litter_posterior_long %>% filter(parameter=='rho') %>% pull(value) %>% min(),
                         litter_posterior_long %>% filter(parameter=='rho') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = invgamma::dinvgamma(x = x, shape = 5, rate = 40)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'c) GP Length Scale')

litter_rho_plot

litter_sigmaPlot_plot = 
  ggplot()+
  geom_density(
    data = litter_posterior_long %>% filter(parameter == 'sigmaPlot'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = litter_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = litter_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = litter_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(litter_posterior_long %>% filter(parameter=='sigmaPlot') %>% pull(value) %>% min(),
                         litter_posterior_long %>% filter(parameter=='sigmaPlot') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 2, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'd) Plot effect SD')

litter_sigmaPlot_plot

litter_kappa_plot = 
  ggplot()+
  geom_density(
    data = litter_posterior_long %>% filter(parameter == 'kappa'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = litter_posterior_long,
                                  param = 'kappa',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = litter_posterior_long,
                                  param = 'kappa',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = litter_posterior_long,
                                  param = 'kappa',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(litter_posterior_long %>% filter(parameter=='kappa') %>% pull(value) %>% min(),
                         litter_posterior_long %>% filter(parameter=='kappa') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 2, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  labs(title = 'e) NB dispersion',
       fill = 'Mortality\nClass',
       color = 'Mortality\nClass')

litter_kappa_plot

litter_plot = 
  plot_grid(litter_beta_plot,
            litter_alpha_plot,
            litter_rho_plot,
            litter_sigmaPlot_plot,
            litter_kappa_plot,
            ncol = 1)

litter_plot

ggsave(litter_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'litter_plot.png'),
       height = 7.5, width = 6.5, units = 'in')

litter_summary = 
  litter_posterior %>%
  select(`beta[1]`, `beta[2]`, `beta[3]`, `alpha[1]`, `alpha[2]`, `alpha[3]`,
         `rho[1]`, `rho[2]`, `rho[3]`, `sigmaPlot[1]`, `sigmaPlot[2]`, 
         `sigmaPlot[3]`, `kappa[1]`, `kappa[2]`, `kappa[3]`, draw = `.draw`) %>%
  pivot_longer(cols = c(-draw),
               names_to = 'param',
               values_to = 'value') %>%
  mutate(group_id = gsub(x = param, replacement = '', pattern = '^.*\\[|\\]'),
         parameter = gsub(x = param, replacement = '', pattern = '\\[.*$'),
         mort_class = 
           factor(ifelse(group_id==1,
                  'low',
                  ifelse(group_id==2,
                         'medium',
                         'high')),
                  levels = c('low', 'medium', 'high'))) %>%
  group_by(parameter, mort_class) %>%
  summarise(value = paste0(round(median(value), 2),' (',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_low, 2),
                           ', ',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_high, 2),
                           ')')) %>%
  ungroup() %>%
  pivot_wider(id_cols = c('parameter'),
              names_from = 'mort_class',
              values_from = 'value') %>%
  left_join(
    data.frame(
      parameter = c('beta', 'alpha', 'rho', 'tau', 'sigmaPlot', 'kappa'),
      p.order = 1:6)
    ) %>%
  arrange(p.order)

litter_summary  

write.csv(litter_summary,
          here::here('04-communication',
                     'tables',
                     'litter_summary.csv'),
          row.names = FALSE)

#### duff ################################################################

duff_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'duff_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'alpha[1]', 'alpha[2]', 'alpha[3]',
                                     'rho[1]', 'rho[2]', 'rho[3]',
                                     'sigmaPlot[1]', 'sigmaPlot[2]', 'sigmaPlot[3]',
                                     'kappa[1]', 'kappa[2]', 'kappa[3]')
                     ) %>%
  as_draws_df()


duff_posterior_long = 
  duff_posterior %>%
  select(contains(c('alpha', 'beta', 'rho', 'sigmaPlot', 'kappa'))) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter_full',
               values_to = 'value') %>%
  mutate(group_id = 
           as.integer(gsub(x = parameter_full, pattern = '^.*\\[|\\]$',replacement='')),
         parameter = 
           gsub(x = parameter_full, pattern = '\\[.*$', replacement = '')) %>%
  left_join(data.frame(group_id = 1:3,
                        mort_class = factor(c('Low', 'Med.', 'High'),
                                            levels = c('Low', 'Med.', 'High'))))



duff_beta_plot = 
  ggplot()+
  geom_density(
    data = duff_posterior_long %>% filter(parameter == 'beta'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = duff_posterior_long,
                                  param = 'beta',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = duff_posterior_long,
                                  param = 'beta',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = duff_posterior_long,
                                  param = 'beta',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(duff_posterior_long %>% filter(parameter=='beta') %>% pull(value) %>% min(),
                         duff_posterior_long %>% filter(parameter=='beta') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = dnorm(x = x, mean = 0, sd = 2)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'a) Intercept')

duff_beta_plot

duff_alpha_plot = 
  ggplot()+
  geom_density(
    data = duff_posterior_long %>% filter(parameter == 'alpha'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = duff_posterior_long,
                                  param = 'alpha',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = duff_posterior_long,
                                  param = 'alpha',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = duff_posterior_long,
                                  param = 'alpha',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(duff_posterior_long %>% filter(parameter=='alpha') %>% pull(value) %>% min(),
                         duff_posterior_long %>% filter(parameter=='alpha') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 2, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'b) GP Magnitude')

duff_alpha_plot

duff_rho_plot = 
  ggplot()+
  geom_density(
    data = duff_posterior_long %>% filter(parameter == 'rho'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = duff_posterior_long,
                                  param = 'rho',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = duff_posterior_long,
                                  param = 'rho',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = duff_posterior_long,
                                  param = 'rho',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(duff_posterior_long %>% filter(parameter=='rho') %>% pull(value) %>% min(),
                         duff_posterior_long %>% filter(parameter=='rho') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = invgamma::dinvgamma(x = x, shape = 5, rate = 40)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'c) GP Length Scale')

duff_rho_plot

duff_sigmaPlot_plot = 
  ggplot()+
  geom_density(
    data = duff_posterior_long %>% filter(parameter == 'sigmaPlot'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = duff_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = duff_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = duff_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(duff_posterior_long %>% filter(parameter=='sigmaPlot') %>% pull(value) %>% min(),
                         duff_posterior_long %>% filter(parameter=='sigmaPlot') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 2, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'd) Plot effect SD')

duff_sigmaPlot_plot

duff_kappa_plot = 
  ggplot()+
  geom_density(
    data = duff_posterior_long %>% filter(parameter == 'kappa'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = duff_posterior_long,
                                  param = 'kappa',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = duff_posterior_long,
                                  param = 'kappa',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = duff_posterior_long,
                                  param = 'kappa',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(duff_posterior_long %>% filter(parameter=='kappa') %>% pull(value) %>% min(),
                         duff_posterior_long %>% filter(parameter=='kappa') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 2, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  labs(title = 'e) NB dispersion',
       fill = 'Mortality\nClass',
       color = 'Mortality\nClass')

duff_kappa_plot

duff_plot = 
  plot_grid(duff_beta_plot,
            duff_alpha_plot,
            duff_rho_plot,
            duff_sigmaPlot_plot,
            duff_kappa_plot,
            ncol = 1)

duff_plot

ggsave(duff_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'duff_plot.png'),
       height = 7.5, width = 6.5, units = 'in')

duff_summary = 
  duff_posterior %>%
  select(`beta[1]`, `beta[2]`, `beta[3]`, `alpha[1]`, `alpha[2]`, `alpha[3]`,
         `rho[1]`, `rho[2]`, `rho[3]`, `sigmaPlot[1]`, `sigmaPlot[2]`, 
         `sigmaPlot[3]`, `kappa[1]`, `kappa[2]`, `kappa[3]`, draw = `.draw`) %>%
  pivot_longer(cols = c(-draw),
               names_to = 'param',
               values_to = 'value') %>%
  mutate(group_id = gsub(x = param, replacement = '', pattern = '^.*\\[|\\]'),
         parameter = gsub(x = param, replacement = '', pattern = '\\[.*$'),
         mort_class = 
           factor(ifelse(group_id==1,
                  'low',
                  ifelse(group_id==2,
                         'medium',
                         'high')),
                  levels = c('low', 'medium', 'high'))) %>%
  group_by(parameter, mort_class) %>%
  summarise(value = paste0(round(median(value), 2),' (',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_low, 2),
                           ', ',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_high, 2),
                           ')')) %>%
  ungroup() %>%
  pivot_wider(id_cols = c('parameter'),
              names_from = 'mort_class',
              values_from = 'value')%>%
  left_join(
    data.frame(
      parameter = c('beta', 'alpha', 'rho', 'tau', 'sigmaPlot', 'kappa'),
      p.order = 1:6)
    ) %>%
  arrange(p.order)

duff_summary  

write.csv(duff_summary,
          here::here('04-communication',
                     'tables',
                     'duff_summary.csv'),
          row.names = FALSE)


#### fwd tallies ###############################################################

fwdtallies_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'fwdtallies_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'alpha[1]', 'alpha[2]', 'alpha[3]',
                                     'rho[1]', 'rho[2]', 'rho[3]',
                                     'sigmaPlot[1]', 'sigmaPlot[2]', 'sigmaPlot[3]',
                                     'tau[1]', 'tau[2]', 'tau[3]',
                                     'kappa[1]', 'kappa[2]', 'kappa[3]')
                     ) %>%
  as_draws_df()


fwdtallies_posterior_long = 
  fwdtallies_posterior %>%
  select(contains(c('alpha', 'beta', 'rho', 'sigmaPlot', 'tau', 'kappa'))) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter_full',
               values_to = 'value') %>%
  mutate(group_id = 
           as.integer(gsub(x = parameter_full, pattern = '^.*\\[|\\]$',replacement='')),
         parameter = 
           gsub(x = parameter_full, pattern = '\\[.*$', replacement = '')) %>%
  left_join(data.frame(group_id = 1:3,
                        mort_class = factor(c('Low', 'Med.', 'High'),
                                            levels = c('Low', 'Med.', 'High'))))



fwdtallies_beta_plot = 
  ggplot()+
  geom_density(
    data = fwdtallies_posterior_long %>% filter(parameter == 'beta'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'beta',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'beta',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'beta',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(fwdtallies_posterior_long %>% filter(parameter=='beta') %>% pull(value) %>% min(),
                         fwdtallies_posterior_long %>% filter(parameter=='beta') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = dnorm(x = x, mean = 0, sd = 2)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'a) Intercept')

fwdtallies_beta_plot

fwdtallies_alpha_plot = 
  ggplot()+
  geom_density(
    data = fwdtallies_posterior_long %>% filter(parameter == 'alpha'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'alpha',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'alpha',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'alpha',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(fwdtallies_posterior_long %>% filter(parameter=='alpha') %>% pull(value) %>% min(),
                         fwdtallies_posterior_long %>% filter(parameter=='alpha') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 2, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'b) GP Magnitude')

fwdtallies_alpha_plot

fwdtallies_rho_plot = 
  ggplot()+
  geom_density(
    data = fwdtallies_posterior_long %>% filter(parameter == 'rho'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'rho',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'rho',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'rho',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(fwdtallies_posterior_long %>% filter(parameter=='rho') %>% pull(value) %>% min(),
                         fwdtallies_posterior_long %>% filter(parameter=='rho') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = invgamma::dinvgamma(x = x, shape = 5, rate = 40)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'c) GP Length Scale')

fwdtallies_rho_plot

fwdtallies_tau_plot = 
  ggplot()+
  geom_density(
    data = fwdtallies_posterior_long %>% filter(parameter == 'tau'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'tau',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'tau',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'tau',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(fwdtallies_posterior_long %>% filter(parameter=='tau') %>% pull(value) %>% min(),
                         fwdtallies_posterior_long %>% filter(parameter=='tau') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 2, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'd) GP Nugget')

fwdtallies_tau_plot

fwdtallies_sigmaPlot_plot = 
  ggplot()+
  geom_density(
    data = fwdtallies_posterior_long %>% filter(parameter == 'sigmaPlot'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(fwdtallies_posterior_long %>% filter(parameter=='sigmaPlot') %>% pull(value) %>% min(),
                         fwdtallies_posterior_long %>% filter(parameter=='sigmaPlot') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 2, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'e) Plot effect SD')

fwdtallies_sigmaPlot_plot

fwdtallies_kappa_plot = 
  ggplot()+
  geom_density(
    data = fwdtallies_posterior_long %>% filter(parameter == 'kappa'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'kappa',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'kappa',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwdtallies_posterior_long,
                                  param = 'kappa',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(fwdtallies_posterior_long %>% filter(parameter=='kappa') %>% pull(value) %>% min(),
                         fwdtallies_posterior_long %>% filter(parameter=='kappa') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 2, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  labs(title = 'f) NB dispersion',
       fill = 'Mortality\nClass',
       color = 'Mortality\nClass')

fwdtallies_kappa_plot

fwdtallies_plot = 
  plot_grid(fwdtallies_beta_plot,
            fwdtallies_alpha_plot,
            fwdtallies_rho_plot,
            fwdtallies_tau_plot,
            fwdtallies_sigmaPlot_plot,
            fwdtallies_kappa_plot,
            ncol = 1)

fwdtallies_plot

ggsave(fwdtallies_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwdtallies_plot.png'),
       height = 7.5, width = 6.5, units = 'in')

fwdtallies_summary = 
  fwdtallies_posterior %>%
  select(`beta[1]`, `beta[2]`, `beta[3]`, `alpha[1]`, `alpha[2]`, `alpha[3]`,
         `rho[1]`, `rho[2]`, `rho[3]`, `tau[1]`, `tau[2]`, `tau[3]`,
         `sigmaPlot[1]`, `sigmaPlot[2]`, 
         `sigmaPlot[3]`, `kappa[1]`, `kappa[2]`, `kappa[3]`, draw = `.draw`) %>%
  pivot_longer(cols = c(-draw),
               names_to = 'param',
               values_to = 'value') %>%
  mutate(group_id = gsub(x = param, replacement = '', pattern = '^.*\\[|\\]'),
         parameter = gsub(x = param, replacement = '', pattern = '\\[.*$'),
         mort_class = 
           factor(ifelse(group_id==1,
                  'low',
                  ifelse(group_id==2,
                         'medium',
                         'high')),
                  levels = c('low', 'medium', 'high'))) %>%
  group_by(parameter, mort_class) %>%
  summarise(value = paste0(round(median(value), 2),' (',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_low, 2),
                           ', ',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_high, 2),
                           ')')) %>%
  ungroup() %>%
  pivot_wider(id_cols = c('parameter'),
              names_from = 'mort_class',
              values_from = 'value')%>%
  left_join(
    data.frame(
      parameter = c('beta', 'alpha', 'rho', 'tau', 'sigmaPlot', 'kappa'),
      p.order = 1:6)
    ) %>%
  arrange(p.order)


fwdtallies_summary  

write.csv(fwdtallies_summary,
          here::here('04-communication',
                     'tables',
                     'fwdtallies_summary.csv'),
          row.names = FALSE)


#### fwd diameters #############################################################

fwddiams_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'fwddiams_fit.rds'))$draws(
                       variables = c('intercept[1]', 'intercept[2]', 'intercept[3]',
                                     'phi[1]', 'phi[2]', 'phi[3]',
                                     'SD_plot[1]', 'SD_plot[2]', 'SD_plot[3]',
                                     'SD_trans[1]', 'SD_trans[2]', 'SD_trans[3]')
                     ) %>%
  as_draws_df()


fwddiams_posterior_long = 
  fwddiams_posterior %>%
  select(contains(c('intercept', 'phi', 'SD_plot', 'SD_trans'))) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter_full',
               values_to = 'value') %>%
  mutate(group_id = 
           as.integer(gsub(x = parameter_full, pattern = '^.*\\[|\\]$',replacement='')),
         parameter = 
           gsub(x = parameter_full, pattern = '\\[.*$', replacement = '')) %>%
  left_join(data.frame(group_id = 1:3,
                        mort_class = factor(c('Low', 'Med.', 'High'),
                                            levels = c('Low', 'Med.', 'High'))))

fwddiams_intercept_plot = 
  ggplot()+
  geom_density(
    data = fwddiams_posterior_long %>% filter(parameter == 'intercept'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = fwddiams_posterior_long,
                                  param = 'intercept',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwddiams_posterior_long,
                                  param = 'intercept',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwddiams_posterior_long,
                                  param = 'intercept',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(fwddiams_posterior_long %>% filter(parameter=='intercept') %>% pull(value) %>% min(),
                         fwddiams_posterior_long %>% filter(parameter=='intercept') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = dnorm(x = x, mean = 0, sd = 1)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'a) Intercept')

fwddiams_intercept_plot

fwddiams_SD_plot_plot = 
  ggplot()+
  geom_density(
    data = fwddiams_posterior_long %>% filter(parameter == 'SD_plot'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = fwddiams_posterior_long,
                                  param = 'SD_plot',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwddiams_posterior_long,
                                  param = 'SD_plot',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwddiams_posterior_long,
                                  param = 'SD_plot',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(fwddiams_posterior_long %>% filter(parameter=='SD_plot') %>% pull(value) %>% min(),
                         fwddiams_posterior_long %>% filter(parameter=='SD_plot') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 0.25, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'b) Plot effect SD')

fwddiams_SD_plot_plot

fwddiams_SD_trans_plot = 
  ggplot()+
  geom_density(
    data = fwddiams_posterior_long %>% filter(parameter == 'SD_trans'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = fwddiams_posterior_long,
                                  param = 'SD_trans',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwddiams_posterior_long,
                                  param = 'SD_trans',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwddiams_posterior_long,
                                  param = 'SD_trans',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(fwddiams_posterior_long %>% filter(parameter=='SD_trans') %>% pull(value) %>% min(),
                         fwddiams_posterior_long %>% filter(parameter=='SD_trans') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 0.25, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'c) Subtransect effect SD')

fwddiams_SD_trans_plot

fwddiams_phi_plot = 
  ggplot()+
  geom_density(
    data = fwddiams_posterior_long %>% filter(parameter == 'phi'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = fwddiams_posterior_long,
                                  param = 'phi',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwddiams_posterior_long,
                                  param = 'phi',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = fwddiams_posterior_long,
                                  param = 'phi',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(fwddiams_posterior_long %>% filter(parameter=='phi') %>% pull(value) %>% min(),
                         fwddiams_posterior_long %>% filter(parameter=='phi') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 5, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  labs(title = 'd) Precision', fill = 'Mortality\nClass',color='Mortality\nClass')

fwddiams_phi_plot

fwddiams_plot = 
  plot_grid(fwddiams_intercept_plot,
            fwddiams_SD_plot_plot,
            fwddiams_SD_trans_plot,
            fwddiams_phi_plot,
            ncol = 1)

fwddiams_plot

ggsave(fwddiams_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'fwddiams_plot.png'),
       height = 6.5, width = 6.5, units = 'in')

fwddiams_summary = 
  fwddiams_posterior %>%
  select(`intercept[1]`, `intercept[2]`, `intercept[3]`,
         `SD_plot[1]`, `SD_plot[2]`, `SD_plot[3]`,
         `SD_trans[1]`, `SD_trans[2]`, `SD_trans[3]`,
         `phi[1]`, `phi[2]`, `phi[3]`, draw = `.draw`) %>%
  pivot_longer(cols = c(-draw),
               names_to = 'param',
               values_to = 'value') %>%
  mutate(group_id = gsub(x = param, replacement = '', pattern = '^.*\\[|\\]'),
         parameter = gsub(x = param, replacement = '', pattern = '\\[.*$'),
         mort_class = 
           factor(ifelse(group_id==1,
                  'low',
                  ifelse(group_id==2,
                         'medium',
                         'high')),
                  levels = c('low', 'medium', 'high'))) %>%
  group_by(parameter, mort_class) %>%
  summarise(value = paste0(round(median(value), 2),' (',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_low, 2),
                           ', ',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_high, 2),
                           ')')) %>%
  ungroup() %>%
  pivot_wider(id_cols = c('parameter'),
              names_from = 'mort_class',
              values_from = 'value')%>%
  left_join(
    data.frame(
      parameter = c('intercept', 'SD_plot', 'SD_trans', 'phi'),
      p.order = 1:4)
    ) %>%
  arrange(p.order)


fwddiams_summary  

write.csv(fwddiams_summary,
          here::here('04-communication',
                     'tables',
                     'fwddiams_summary.csv'),
          row.names = FALSE)


#### cwd tallies ###############################################################


cwdtallies_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'cwdtallies_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'alpha[1]', 'alpha[2]', 'alpha[3]',
                                     'rho[1]', 'rho[2]', 'rho[3]',
                                     'sigmaPlot[1]', 'sigmaPlot[2]', 'sigmaPlot[3]',
                                     'sigmaTrans[1]', 'sigmaTrans[2]', 'sigmaTrans[3]')
                     ) %>%
  as_draws_df()


cwdtallies_posterior_long = 
  cwdtallies_posterior %>%
  select(contains(c('alpha', 'beta', 'rho', 'sigmaPlot', 'sigmaTrans'))) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter_full',
               values_to = 'value') %>%
  mutate(group_id = 
           as.integer(gsub(x = parameter_full, pattern = '^.*\\[|\\]$',replacement='')),
         parameter = 
           gsub(x = parameter_full, pattern = '\\[.*$', replacement = '')) %>%
  left_join(data.frame(group_id = 1:3,
                        mort_class = factor(c('Low', 'Med.', 'High'),
                                            levels = c('Low', 'Med.', 'High'))))



cwdtallies_beta_plot = 
  ggplot()+
  geom_density(
    data = cwdtallies_posterior_long %>% filter(parameter == 'beta'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = cwdtallies_posterior_long,
                                  param = 'beta',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = cwdtallies_posterior_long,
                                  param = 'beta',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = cwdtallies_posterior_long,
                                  param = 'beta',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(cwdtallies_posterior_long %>% filter(parameter=='beta') %>% pull(value) %>% min(),
                         cwdtallies_posterior_long %>% filter(parameter=='beta') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = dnorm(x = x, mean = 0, sd = 0.5)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'a) Intercept')

cwdtallies_beta_plot

cwdtallies_alpha_plot = 
  ggplot()+
  geom_density(
    data = cwdtallies_posterior_long %>% filter(parameter == 'alpha'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = cwdtallies_posterior_long,
                                  param = 'alpha',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = cwdtallies_posterior_long,
                                  param = 'alpha',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = cwdtallies_posterior_long,
                                  param = 'alpha',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(cwdtallies_posterior_long %>% filter(parameter=='alpha') %>% pull(value) %>% min(),
                         cwdtallies_posterior_long %>% filter(parameter=='alpha') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 0.5, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'b) GP Magnitude')

cwdtallies_alpha_plot

cwdtallies_rho_plot = 
  ggplot()+
  geom_density(
    data = cwdtallies_posterior_long %>% filter(parameter == 'rho'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = cwdtallies_posterior_long,
                                  param = 'rho',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = cwdtallies_posterior_long,
                                  param = 'rho',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = cwdtallies_posterior_long,
                                  param = 'rho',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(cwdtallies_posterior_long %>% filter(parameter=='rho') %>% pull(value) %>% min(),
                         cwdtallies_posterior_long %>% filter(parameter=='rho') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = invgamma::dinvgamma(x = x, shape = 5, rate = 40)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'c) GP Length Scale')

cwdtallies_rho_plot

cwdtallies_sigmaPlot_plot = 
  ggplot()+
  geom_density(
    data = cwdtallies_posterior_long %>% filter(parameter == 'sigmaPlot'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = cwdtallies_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = cwdtallies_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = cwdtallies_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(cwdtallies_posterior_long %>% filter(parameter=='sigmaPlot') %>% pull(value) %>% min(),
                         cwdtallies_posterior_long %>% filter(parameter=='sigmaPlot') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 0.5, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'd) Plot effect SD')

cwdtallies_sigmaPlot_plot

cwdtallies_sigmaTrans_plot = 
  ggplot()+
  geom_density(
    data = cwdtallies_posterior_long %>% filter(parameter == 'sigmaTrans'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = cwdtallies_posterior_long,
                                  param = 'sigmaTrans',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = cwdtallies_posterior_long,
                                  param = 'sigmaTrans',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = cwdtallies_posterior_long,
                                  param = 'sigmaTrans',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(cwdtallies_posterior_long %>% filter(parameter=='sigmaTrans') %>% pull(value) %>% min(),
                         cwdtallies_posterior_long %>% filter(parameter=='sigmaTrans') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 0.5, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  labs(title = 'e) Azimuth effect SD',
       fill = 'Mortality\nClass',
       color = 'Mortality\nClass')

cwdtallies_sigmaTrans_plot



cwdtallies_plot = 
  plot_grid(cwdtallies_beta_plot,
            cwdtallies_alpha_plot,
            cwdtallies_rho_plot,
            cwdtallies_sigmaPlot_plot,
            cwdtallies_sigmaTrans_plot,
            ncol = 1)

cwdtallies_plot

ggsave(cwdtallies_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'cwdtallies_plot.png'),
       height = 7.5, width = 6.5, units = 'in')



cwdtallies_summary = 
  cwdtallies_posterior %>%
  select(`beta[1]`, `beta[2]`, `beta[3]`, `alpha[1]`, `alpha[2]`, `alpha[3]`,
         `rho[1]`, `rho[2]`, `rho[3]`, `sigmaTrans[1]`, `sigmaTrans[2]`, `sigmaTrans[3]`,
         `sigmaPlot[1]`, `sigmaPlot[2]`, `sigmaPlot[3]`, draw = `.draw`) %>%
  pivot_longer(cols = c(-draw),
               names_to = 'param',
               values_to = 'value') %>%
  mutate(group_id = gsub(x = param, replacement = '', pattern = '^.*\\[|\\]'),
         parameter = gsub(x = param, replacement = '', pattern = '\\[.*$'),
         mort_class = 
           factor(ifelse(group_id==1,
                  'low',
                  ifelse(group_id==2,
                         'medium',
                         'high')),
                  levels = c('low', 'medium', 'high'))) %>%
  group_by(parameter, mort_class) %>%
  summarise(value = paste0(round(median(value), 2),' (',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_low, 2),
                           ', ',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_high, 2),
                           ')')) %>%
  ungroup() %>%
  pivot_wider(id_cols = c('parameter'),
              names_from = 'mort_class',
              values_from = 'value')%>%
  left_join(
    data.frame(
      parameter = c('beta', 'alpha', 'rho', 'sigmaPlot', 'sigmaTrans'),
      p.order = 1:5)
    ) %>%
  arrange(p.order)


cwdtallies_summary  

write.csv(cwdtallies_summary,
          here::here('04-communication',
                     'tables',
                     'cwdtallies_summary.csv'),
          row.names = FALSE)

#### cwd diameters #############################################################

cwddiams_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'cwddiams_fit.rds'))$draws(
                       variables = c('mu[1]', 'mu[2]', 'mu[3]',
                                     'phi[1]', 'phi[2]', 'phi[3]')) %>%
  as_draws_df()


cwddiams_posterior_long = 
  cwddiams_posterior %>%
  select(contains(c('mu', 'phi'))) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter_full',
               values_to = 'value') %>%
  mutate(group_id = 
           as.integer(gsub(x = parameter_full, pattern = '^.*\\[|\\]$',replacement='')),
         parameter = 
           gsub(x = parameter_full, pattern = '\\[.*$', replacement = '')) %>%
  left_join(data.frame(group_id = 1:3,
                        mort_class = factor(c('Low', 'Med.', 'High'),
                                            levels = c('Low', 'Med.', 'High'))))

cwddiams_mu_plot = 
  ggplot()+
  geom_density(
    data = cwddiams_posterior_long %>% filter(parameter == 'mu'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = cwddiams_posterior_long,
                                  param = 'mu',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = cwddiams_posterior_long,
                                  param = 'mu',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = cwddiams_posterior_long,
                                  param = 'mu',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(cwddiams_posterior_long %>% filter(parameter=='mu') %>% pull(value) %>% min(),
                         cwddiams_posterior_long %>% filter(parameter=='mu') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = dunif(x = x, 
                       min = 
                         readRDS(here::here('02-data', '05-for_analysis', 
                                            'cwddiams_training_data.rds'))$lb,
                       max = 
                         readRDS(here::here('02-data', '05-for_analysis', 
                                            'cwddiams_training_data.rds'))$ub)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'a) Mean')

cwddiams_mu_plot


cwddiams_phi_plot = 
  ggplot()+
  geom_density(
    data = cwddiams_posterior_long %>% filter(parameter == 'phi'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = cwddiams_posterior_long,
                                  param = 'phi',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = cwddiams_posterior_long,
                                  param = 'phi',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = cwddiams_posterior_long,
                                  param = 'phi',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(cwddiams_posterior_long %>% filter(parameter=='phi') %>% pull(value) %>% min(),
                         cwddiams_posterior_long %>% filter(parameter=='phi') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 10, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  labs(title = 'b) Precision', fill = 'Mortality\nClass',color='Mortality\nClass')

cwddiams_phi_plot

cwddiams_plot = 
  plot_grid(cwddiams_mu_plot,
            cwddiams_phi_plot,
            ncol = 1)

cwddiams_plot

ggsave(cwddiams_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'cwddiams_plot.png'),
       height = 4.5, width = 6.5, units = 'in')

cwddiams_summary = 
  cwddiams_posterior %>%
  select(`mu[1]`, `mu[2]`, `mu[3]`,
         `phi[1]`, `phi[2]`, `phi[3]`, draw = `.draw`) %>%
  pivot_longer(cols = c(-draw),
               names_to = 'param',
               values_to = 'value') %>%
  mutate(group_id = gsub(x = param, replacement = '', pattern = '^.*\\[|\\]'),
         parameter = gsub(x = param, replacement = '', pattern = '\\[.*$'),
         mort_class = 
           factor(ifelse(group_id==1,
                  'low',
                  ifelse(group_id==2,
                         'medium',
                         'high')),
                  levels = c('low', 'medium', 'high'))) %>%
  group_by(parameter, mort_class) %>%
  summarise(value = paste0(round(median(value), 2),' (',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_low, 2),
                           ', ',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_high, 2),
                           ')')) %>%
  ungroup() %>%
  pivot_wider(id_cols = c('parameter'),
              names_from = 'mort_class',
              values_from = 'value')%>%
  left_join(
    data.frame(
      parameter = c('mu', 'phi'),
      p.order = 1:2)
    ) %>%
  arrange(p.order)


cwddiams_summary  

write.csv(cwddiams_summary,
          here::here('04-communication',
                     'tables',
                     'cwddiams_summary.csv'),
          row.names = FALSE)


#### veg pres/abs ##############################################################


vegpa_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'vegpa_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'alpha[1]', 'alpha[2]', 'alpha[3]',
                                     'rho[1]', 'rho[2]', 'rho[3]',
                                     'sigmaPlot[1]', 'sigmaPlot[2]', 'sigmaPlot[3]')
                     ) %>%
  as_draws_df()


vegpa_posterior_long = 
  vegpa_posterior %>%
  select(contains(c('alpha', 'beta', 'rho', 'sigmaPlot'))) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter_full',
               values_to = 'value') %>%
  mutate(group_id = 
           as.integer(gsub(x = parameter_full, pattern = '^.*\\[|\\]$',replacement='')),
         parameter = 
           gsub(x = parameter_full, pattern = '\\[.*$', replacement = '')) %>%
  left_join(data.frame(group_id = 1:3,
                        mort_class = factor(c('Low', 'Med.', 'High'),
                                            levels = c('Low', 'Med.', 'High'))))



vegpa_beta_plot = 
  ggplot()+
  geom_density(
    data = vegpa_posterior_long %>% filter(parameter == 'beta'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = vegpa_posterior_long,
                                  param = 'beta',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = vegpa_posterior_long,
                                  param = 'beta',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = vegpa_posterior_long,
                                  param = 'beta',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(vegpa_posterior_long %>% filter(parameter=='beta') %>% pull(value) %>% min(),
                         vegpa_posterior_long %>% filter(parameter=='beta') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = dnorm(x = x, mean = 0, sd = 1.5)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'a) Intercept')

vegpa_beta_plot

vegpa_alpha_plot = 
  ggplot()+
  geom_density(
    data = vegpa_posterior_long %>% filter(parameter == 'alpha'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = vegpa_posterior_long,
                                  param = 'alpha',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = vegpa_posterior_long,
                                  param = 'alpha',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = vegpa_posterior_long,
                                  param = 'alpha',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(vegpa_posterior_long %>% filter(parameter=='alpha') %>% pull(value) %>% min(),
                         vegpa_posterior_long %>% filter(parameter=='alpha') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 1.5, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'b) GP Magnitude')

vegpa_alpha_plot

vegpa_rho_plot = 
  ggplot()+
  geom_density(
    data = vegpa_posterior_long %>% filter(parameter == 'rho'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = vegpa_posterior_long,
                                  param = 'rho',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = vegpa_posterior_long,
                                  param = 'rho',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = vegpa_posterior_long,
                                  param = 'rho',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(vegpa_posterior_long %>% filter(parameter=='rho') %>% pull(value) %>% min(),
                         vegpa_posterior_long %>% filter(parameter=='rho') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = invgamma::dinvgamma(x = x, shape = 5, rate = 40)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'c) GP Length Scale')

vegpa_rho_plot

vegpa_sigmaPlot_plot = 
  ggplot()+
  geom_density(
    data = vegpa_posterior_long %>% filter(parameter == 'sigmaPlot'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = vegpa_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = vegpa_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = vegpa_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(vegpa_posterior_long %>% filter(parameter=='sigmaPlot') %>% pull(value) %>% min(),
                         vegpa_posterior_long %>% filter(parameter=='sigmaPlot') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 1.5, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  labs(title = 'd) Plot effect SD',
       fill = 'Mortality\nClass',
       color = 'Mortality\nClass')

vegpa_sigmaPlot_plot


vegpa_plot = 
  plot_grid(vegpa_beta_plot,
            vegpa_alpha_plot,
            vegpa_rho_plot,
            vegpa_sigmaPlot_plot,
            ncol = 1)

vegpa_plot

ggsave(vegpa_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'vegpa_plot.png'),
       height = 6.5, width = 6.5, units = 'in')

vegpa_summary = 
  vegpa_posterior %>%
  select(`beta[1]`, `beta[2]`, `beta[3]`, `alpha[1]`, `alpha[2]`, `alpha[3]`,
         `rho[1]`, `rho[2]`, `rho[3]`, `sigmaPlot[1]`, `sigmaPlot[2]`, 
         `sigmaPlot[3]`, draw = `.draw`) %>%
  pivot_longer(cols = c(-draw),
               names_to = 'param',
               values_to = 'value') %>%
  mutate(group_id = gsub(x = param, replacement = '', pattern = '^.*\\[|\\]'),
         parameter = gsub(x = param, replacement = '', pattern = '\\[.*$'),
         mort_class = 
           factor(ifelse(group_id==1,
                         'low',
                         ifelse(group_id==2,
                                'medium',
                                'high')),
                  levels = c('low', 'medium', 'high'))) %>%
  group_by(parameter, mort_class) %>%
  summarise(value = paste0(round(median(value), 2),' (',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_low, 2),
                           ', ',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_high, 2),
                           ')')) %>%
  ungroup() %>%
  pivot_wider(id_cols = c('parameter'),
              names_from = 'mort_class',
              values_from = 'value') %>%
  left_join(
    data.frame(
      parameter = c('beta', 'alpha', 'rho', 'sigmaPlot'),
      p.order = 1:4)
  ) %>%
  arrange(p.order)

vegpa_summary  

write.csv(vegpa_summary,
          here::here('04-communication',
                     'tables',
                     'vegpa_summary.csv'),
          row.names = FALSE)


#### trees #####################################################################

trees_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'trees_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'alpha[1]', 'alpha[2]', 'alpha[3]',
                                     'rho[1]', 'rho[2]', 'rho[3]',
                                     'sigmaPlot[1]', 'sigmaPlot[2]', 'sigmaPlot[3]')
                     ) %>%
  as_draws_df()


trees_posterior_long = 
  trees_posterior %>%
  select(contains(c('alpha', 'beta', 'rho', 'sigmaPlot'))) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter_full',
               values_to = 'value') %>%
  mutate(group_id = 
           as.integer(gsub(x = parameter_full, pattern = '^.*\\[|\\]$',replacement='')),
         parameter = 
           gsub(x = parameter_full, pattern = '\\[.*$', replacement = '')) %>%
  left_join(data.frame(group_id = 1:3,
                        mort_class = factor(c('Low', 'Med.', 'High'),
                                            levels = c('Low', 'Med.', 'High'))))



trees_beta_plot = 
  ggplot()+
  geom_density(
    data = trees_posterior_long %>% filter(parameter == 'beta'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = trees_posterior_long,
                                  param = 'beta',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = trees_posterior_long,
                                  param = 'beta',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = trees_posterior_long,
                                  param = 'beta',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(trees_posterior_long %>% filter(parameter=='beta') %>% pull(value) %>% min(),
                         trees_posterior_long %>% filter(parameter=='beta') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = dnorm(x = x, mean = -1, sd = 0.5)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'a) Intercept')

trees_beta_plot

trees_alpha_plot = 
  ggplot()+
  geom_density(
    data = trees_posterior_long %>% filter(parameter == 'alpha'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = trees_posterior_long,
                                  param = 'alpha',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = trees_posterior_long,
                                  param = 'alpha',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = trees_posterior_long,
                                  param = 'alpha',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(trees_posterior_long %>% filter(parameter=='alpha') %>% pull(value) %>% min(),
                         trees_posterior_long %>% filter(parameter=='alpha') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 0.5, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'b) GP Magnitude')

trees_alpha_plot

trees_rho_plot = 
  ggplot()+
  geom_density(
    data = trees_posterior_long %>% filter(parameter == 'rho'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = trees_posterior_long,
                                  param = 'rho',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = trees_posterior_long,
                                  param = 'rho',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = trees_posterior_long,
                                  param = 'rho',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(trees_posterior_long %>% filter(parameter=='rho') %>% pull(value) %>% min(),
                         trees_posterior_long %>% filter(parameter=='rho') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = invgamma::dinvgamma(x = x, shape = 5, rate = 40)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'c) GP Length Scale')

trees_rho_plot

trees_sigmaPlot_plot = 
  ggplot()+
  geom_density(
    data = trees_posterior_long %>% filter(parameter == 'sigmaPlot'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = trees_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = trees_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = trees_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(trees_posterior_long %>% filter(parameter=='sigmaPlot') %>% pull(value) %>% min(),
                         trees_posterior_long %>% filter(parameter=='sigmaPlot') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 0.5, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  labs(title = 'd) Plot effect SD',
       fill = 'Mortality\nClass',
       color = 'Mortality\nClass')

trees_sigmaPlot_plot


trees_plot = 
  plot_grid(trees_beta_plot,
            trees_alpha_plot,
            trees_rho_plot,
            trees_sigmaPlot_plot,
            ncol = 1)

trees_plot

ggsave(trees_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'trees_plot.png'),
       height = 6.5, width = 6.5, units = 'in')

trees_summary = 
  trees_posterior %>%
  select(`beta[1]`, `beta[2]`, `beta[3]`, `alpha[1]`, `alpha[2]`, `alpha[3]`,
         `rho[1]`, `rho[2]`, `rho[3]`, `sigmaPlot[1]`, `sigmaPlot[2]`, 
         `sigmaPlot[3]`, draw = `.draw`) %>%
  pivot_longer(cols = c(-draw),
               names_to = 'param',
               values_to = 'value') %>%
  mutate(group_id = gsub(x = param, replacement = '', pattern = '^.*\\[|\\]'),
         parameter = gsub(x = param, replacement = '', pattern = '\\[.*$'),
         mort_class = 
           factor(ifelse(group_id==1,
                  'low',
                  ifelse(group_id==2,
                         'medium',
                         'high')),
                  levels = c('low', 'medium', 'high'))) %>%
  group_by(parameter, mort_class) %>%
  summarise(value = paste0(round(median(value), 2),' (',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_low, 2),
                           ', ',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_high, 2),
                           ')')) %>%
  ungroup() %>%
  pivot_wider(id_cols = c('parameter'),
              names_from = 'mort_class',
              values_from = 'value') %>%
  left_join(
    data.frame(
      parameter = c('beta', 'alpha', 'rho', 'tau', 'sigmaPlot', 'kappa'),
      p.order = 1:6)
    ) %>%
  arrange(p.order)

trees_summary  

write.csv(trees_summary,
          here::here('04-communication',
                     'tables',
                     'trees_summary.csv'),
          row.names = FALSE)


#### saplings #####################################################################

saplings_posterior = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'saplings_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'alpha[1]', 'alpha[2]', 'alpha[3]',
                                     'rho[1]', 'rho[2]', 'rho[3]',
                                     'sigmaPlot[1]', 'sigmaPlot[2]', 'sigmaPlot[3]')
                     ) %>%
  as_draws_df()

saplings_posterior_long = 
  saplings_posterior %>%
  select(contains(c('alpha', 'beta', 'rho', 'sigmaPlot'))) %>%
  rowid_to_column('draw') %>%
  pivot_longer(cols = c(-draw),
               names_to = 'parameter_full',
               values_to = 'value') %>%
  mutate(group_id = 
           as.integer(gsub(x = parameter_full, pattern = '^.*\\[|\\]$',replacement='')),
         parameter = 
           gsub(x = parameter_full, pattern = '\\[.*$', replacement = '')) %>%
  left_join(data.frame(group_id = 1:3,
                       mort_class = factor(c('Low', 'Med.', 'High'),
                                           levels = c('Low', 'Med.', 'High'))))



saplings_beta_plot = 
  ggplot()+
  geom_density(
    data = saplings_posterior_long %>% filter(parameter == 'beta'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = saplings_posterior_long,
                                  param = 'beta',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = saplings_posterior_long,
                                  param = 'beta',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = saplings_posterior_long,
                                  param = 'beta',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(saplings_posterior_long %>% filter(parameter=='beta') %>% pull(value) %>% min(),
                         saplings_posterior_long %>% filter(parameter=='beta') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = dnorm(x = x, mean = -1, sd = 0.5)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'a) Intercept')

saplings_beta_plot

saplings_alpha_plot = 
  ggplot()+
  geom_density(
    data = saplings_posterior_long %>% filter(parameter == 'alpha'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = saplings_posterior_long,
                                  param = 'alpha',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = saplings_posterior_long,
                                  param = 'alpha',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = saplings_posterior_long,
                                  param = 'alpha',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(saplings_posterior_long %>% filter(parameter=='alpha') %>% pull(value) %>% min(),
                         saplings_posterior_long %>% filter(parameter=='alpha') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 0.5, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'b) GP Magnitude')

saplings_alpha_plot

saplings_rho_plot = 
  ggplot()+
  geom_density(
    data = saplings_posterior_long %>% filter(parameter == 'rho'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = saplings_posterior_long,
                                  param = 'rho',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = saplings_posterior_long,
                                  param = 'rho',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = saplings_posterior_long,
                                  param = 'rho',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(saplings_posterior_long %>% filter(parameter=='rho') %>% pull(value) %>% min(),
                         saplings_posterior_long %>% filter(parameter=='rho') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = invgamma::dinvgamma(x = x, shape = 5, rate = 40)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = 'none')+
  labs(title = 'c) GP Length Scale')

saplings_rho_plot

saplings_sigmaPlot_plot = 
  ggplot()+
  geom_density(
    data = saplings_posterior_long %>% filter(parameter == 'sigmaPlot'),
    aes(x = value, color = mort_class)
  )+
  geom_ribbon(
    data = build_density_slice_df(post = saplings_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 1),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = saplings_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 2),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_ribbon(
    data = build_density_slice_df(post = saplings_posterior_long,
                                  param = 'sigmaPlot',
                                  g = 3),
    aes(ymax = y, ymin = 0, x = x, fill = mort_class),
    alpha = 0.5)+
  geom_line(
    data = 
      data.frame(x = seq(saplings_posterior_long %>% filter(parameter=='sigmaPlot') %>% pull(value) %>% min(),
                         saplings_posterior_long %>% filter(parameter=='sigmaPlot') %>% pull(value) %>% max(),
                         length.out = 1000)) %>%
      mutate(y = truncnorm::dtruncnorm(x = x, mean = 0, sd = 0.5, a = 0, b = Inf)),
    aes(x = x, y = y),
    lty = 2, color = 'red'
  )+
  theme_minimal()+
  scale_color_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C',
                       breaks = c('1', '2', '3'),
                       labels = c('Low', 'Med.', 'High'))+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  labs(title = 'd) Plot effect SD',
       fill = 'Mortality\nClass',
       color = 'Mortality\nClass')

saplings_sigmaPlot_plot


saplings_plot = 
  plot_grid(saplings_beta_plot,
            saplings_alpha_plot,
            saplings_rho_plot,
            saplings_sigmaPlot_plot,
            ncol = 1)

saplings_plot

ggsave(saplings_plot,
       filename = here::here('04-communication',
                             'figures',
                             'manuscript',
                             'saplings_plot.png'),
       height = 6.5, width = 6.5, units = 'in')

saplings_summary = 
  saplings_posterior %>%
  select(`beta[1]`, `beta[2]`, `beta[3]`, `alpha[1]`, `alpha[2]`, `alpha[3]`,
         `rho[1]`, `rho[2]`, `rho[3]`, `sigmaPlot[1]`, `sigmaPlot[2]`, 
         `sigmaPlot[3]`, draw = `.draw`) %>%
  pivot_longer(cols = c(-draw),
               names_to = 'param',
               values_to = 'value') %>%
  mutate(group_id = gsub(x = param, replacement = '', pattern = '^.*\\[|\\]'),
         parameter = gsub(x = param, replacement = '', pattern = '\\[.*$'),
         mort_class = 
           factor(ifelse(group_id==1,
                  'low',
                  ifelse(group_id==2,
                         'medium',
                         'high')),
                  levels = c('low', 'medium', 'high'))) %>%
  group_by(parameter, mort_class) %>%
  summarise(value = paste0(round(median(value), 2),' (',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_low, 2),
                           ', ',
                           round(bayestestR::ci(x = value, ci = 0.89, method = 'SPI')$CI_high, 2),
                           ')')) %>%
  ungroup() %>%
  pivot_wider(id_cols = c('parameter'),
              names_from = 'mort_class',
              values_from = 'value') %>%
  left_join(
    data.frame(
      parameter = c('beta', 'alpha', 'rho', 'tau', 'sigmaPlot', 'kappa'),
      p.order = 1:6)
    ) %>%
  arrange(p.order)

saplings_summary  

write.csv(saplings_summary,
          here::here('04-communication',
                     'tables',
                     'saplings_summary.csv'),
          row.names = FALSE)




#### scratch ###################################################################


# load fitted models and extract median parameter values,
litterduff_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'litterduff_fit.rds'))$draws(
                       variables = 
                         c('beta[1]', 'beta[2]', 'beta[3]',
                           'alpha[1]', 'alpha[2]', 'alpha[3]',
                           'rho[1]', 'rho[2]', 'rho[3]',
                           'sigmaPlot', 'kappa')
                     ) %>%
  as_draws_df()


fwddiams_params = 
  readRDS(here::here('02-data', 
                     '06-results', 
                     'real_fits',
                     'fwddiams_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'alpha[1]', 'alpha[2]', 'alpha[3]')
                     ) %>%
  as_draws_df()

fwd1h_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'fwd1h_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'alpha[1]', 'alpha[2]', 'alpha[3]',
                                     'rho[1]', 'rho[2]', 'rho[3]',
                                     'tau[1]', 'tau[2]', 'tau[3]',
                                     'sigmaPlot', 'kappa')
                     ) %>%
  as_draws_df()


fwd10h_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'fwd10h_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'alpha[1]', 'alpha[2]', 'alpha[3]',
                                     'rho[1]', 'rho[2]', 'rho[3]',
                                     'tau[1]', 'tau[2]', 'tau[3]',
                                     'sigmaPlot', 'kappa')
                     ) %>%
  as_draws_df()

fwd100h_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'fwd100h_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'alpha[1]', 'alpha[2]', 'alpha[3]',
                                     'rho[1]', 'rho[2]', 'rho[3]',
                                     'tau[1]', 'tau[2]', 'tau[3]',
                                     'sigmaPlot', 'kappa')
                     ) %>%
  as_draws_df()

veg_params = 
  readRDS(here::here('02-data',
                     '06-results',
                     'real_fits',
                     'veg_fit.rds'))$draws(
                       variables = c('beta[1]', 'beta[2]', 'beta[3]',
                                     'alpha[1]', 'alpha[2]', 'alpha[3]',
                                     'rho[1]', 'rho[2]', 'rho[3]',
                                     'sigmaPlot')
                     ) %>%
  as_draws_df()

#### compare across mortality classes ##########################################

head(litterduff_params)

litterduff_params %>%
  select(`beta[1]`, `beta[2]`, `beta[3]`,
         `alpha[1]`, `alpha[2]`, `alpha[3]`,
         `rho[1]`, `rho[2]`, `rho[3]`) %>%
  pivot_longer(cols = c('beta[1]', 'beta[2]', 'beta[3]',
                        'alpha[1]', 'alpha[2]', 'alpha[3]',
                        'rho[1]', 'rho[2]', 'rho[3]'),
               names_to = c('parameter', 'group_id'),
               names_sep= '\\[') %>%
  mutate(group_id = gsub(x = group_id, pattern = '\\]', replacement = '')) %>%
  bind_rows(data.frame(parameter = 
                         c(rep('alpha', 4000),
                           rep('beta', 4000),
                           rep('rho', 4000)),
                       group_id = rep('prior', 12000),
                       value = 
                         c(truncnorm::rtruncnorm(n = 4000, a = 0, mean = 0, sd = 0.5),
                           rnorm(n = 4000, mean = 0, sd = 2),
                           invgamma::rinvgamma(n = 4000, shape = 5, rate = 40)))) %>%
  ggplot(aes(x = value, fill = group_id))+
  geom_density(alpha = 0.5)+
  facet_wrap(~parameter,ncol = 1, scales = 'free')+
  theme_minimal()+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C')

# zoom in on rho
litterduff_params %>%
  select(`rho[1]`, `rho[2]`, `rho[3]`) %>%
  pivot_longer(cols = c('rho[1]', 'rho[2]', 'rho[3]'),
               names_to = c('parameter', 'group_id'),
               names_sep= '\\[') %>%
  mutate(group_id = gsub(x = group_id, pattern = '\\]', replacement = '')) %>%
  bind_rows(data.frame(parameter = 
                         c(rep('rho', 4000)),
                       group_id = rep('prior', 4000),
                       value = 
                         c(invgamma::rinvgamma(n = 4000, shape = 5, rate = 40)))) %>%
  ggplot(aes(x = value, fill = group_id))+
  geom_density(alpha = 0.5)+
  facet_wrap(~parameter,ncol = 1, scales = 'free')+
  theme_minimal()+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
  scale_x_continuous(limits = c(0, 20))


# litterduff posteriors well informed by data; looks like:
#   - mortality class has no effect on loading (beta[2] and beta[3] overlap 0)
#   - mid mortality has a higher range than low/high
#   - high mortality maybe a little lower alpha than low/mid, but lots of overlap

fwd1h_params %>%
  select(`beta[1]`, `beta[2]`, `beta[3]`,
         `alpha[1]`, `alpha[2]`, `alpha[3]`,
         `rho[1]`, `rho[2]`, `rho[3]`,
         `tau[1]`, `tau[2]`, `tau[3]`) %>%
  pivot_longer(cols = c('beta[1]', 'beta[2]', 'beta[3]',
                        'alpha[1]', 'alpha[2]', 'alpha[3]',
                        'rho[1]', 'rho[2]', 'rho[3]',
                        'tau[1]', 'tau[2]', 'tau[3]'),
               names_to = c('parameter', 'group_id'),
               names_sep= '\\[') %>%
  mutate(group_id = gsub(x = group_id, pattern = '\\]', replacement = '')) %>%
  bind_rows(data.frame(parameter = 
                         c(rep('alpha', 4000),
                           rep('beta', 4000),
                           rep('rho', 4000),
                           rep('tau', 4000)),
                       group_id = rep('prior', 16000),
                       value = 
                         c(truncnorm::rtruncnorm(n = 4000, a = 0, mean = 0, sd = 0.5),
                           rnorm(n = 4000, mean = 0, sd = 2),
                           invgamma::rinvgamma(n = 4000, shape = 5, rate = 40),
                           truncnorm::rtruncnorm(n = 4000, a = 0, mean = 0, sd = 0.5)))) %>%
  ggplot(aes(x = value, fill = group_id))+
  geom_density(alpha = 0.5)+
  facet_wrap(~parameter,ncol = 1, scales = 'free')+
  theme_minimal()+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C')

# - FWD 1h params well informed by data
# - no effect of mortality on loading, rho, tau
# - mild positive effect of high mortality on alpha


fwd10h_params %>%
  select(`beta[1]`, `beta[2]`, `beta[3]`,
         `alpha[1]`, `alpha[2]`, `alpha[3]`,
         `rho[1]`, `rho[2]`, `rho[3]`,
         `tau[1]`, `tau[2]`, `tau[3]`) %>%
  pivot_longer(cols = c('beta[1]', 'beta[2]', 'beta[3]',
                        'alpha[1]', 'alpha[2]', 'alpha[3]',
                        'rho[1]', 'rho[2]', 'rho[3]',
                        'tau[1]', 'tau[2]', 'tau[3]'),
               names_to = c('parameter', 'group_id'),
               names_sep= '\\[') %>%
  mutate(group_id = gsub(x = group_id, pattern = '\\]', replacement = '')) %>%
  bind_rows(data.frame(parameter = 
                         c(rep('alpha', 4000),
                           rep('beta', 4000),
                           rep('rho', 4000),
                           rep('tau', 4000)),
                       group_id = rep('prior', 16000),
                       value = 
                         c(truncnorm::rtruncnorm(n = 4000, a = 0, mean = 0, sd = 0.5),
                           rnorm(n = 4000, mean = 0, sd = 2),
                           invgamma::rinvgamma(n = 4000, shape = 5, rate = 40),
                           truncnorm::rtruncnorm(n = 4000, a = 0, mean = 0, sd = 0.5)))) %>%
  ggplot(aes(x = value, fill = group_id))+
  geom_density(alpha = 0.5)+
  facet_wrap(~parameter,ncol = 1, scales = 'free')+
  theme_minimal()+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C')

# - other than rho, posteriors well informed by the data. posteriors for 
# rho[1] and rho[2] pretty similar to the prior
# - no effect of mortality on loading
# - high mortality has a higher alpha than low/mid
# - high mortality has a lower rho than low/mid
# - high mortality has a lower tau than mid, which has a slightly lower tau than 
#  low


fwd100h_params %>%
  select(`beta[1]`, `beta[2]`, `beta[3]`,
         `alpha[1]`, `alpha[2]`, `alpha[3]`,
         `rho[1]`, `rho[2]`, `rho[3]`,
         `tau[1]`, `tau[2]`, `tau[3]`) %>%
  pivot_longer(cols = c('beta[1]', 'beta[2]', 'beta[3]',
                        'alpha[1]', 'alpha[2]', 'alpha[3]',
                        'rho[1]', 'rho[2]', 'rho[3]',
                        'tau[1]', 'tau[2]', 'tau[3]'),
               names_to = c('parameter', 'group_id'),
               names_sep= '\\[') %>%
  mutate(group_id = gsub(x = group_id, pattern = '\\]', replacement = '')) %>%
  bind_rows(data.frame(parameter = 
                         c(rep('alpha', 4000),
                           rep('beta', 4000),
                           rep('rho', 4000),
                           rep('tau', 4000)),
                       group_id = rep('prior', 16000),
                       value = 
                         c(truncnorm::rtruncnorm(n = 4000, a = 0, mean = 0, sd = 0.5),
                           rnorm(n = 4000, mean = 0, sd = 2),
                           invgamma::rinvgamma(n = 4000, shape = 5, rate = 40),
                           truncnorm::rtruncnorm(n = 4000, a = 0, mean = 0, sd = 0.5)))) %>%
  ggplot(aes(x = value, fill = group_id))+
  geom_density(alpha = 0.5)+
  facet_wrap(~parameter,ncol = 1, scales = 'free')+
  theme_minimal()+
  scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C')

# most parameters (other than beta) not well informed by the data; esp rho

#### compare across FWD size classes ###########################################

litterduff_params %>%
  select(contains('beta')) %>%
  mutate(component = 'litterduff') %>%
  bind_rows(
    fwd1h_params %>%
      select(contains('beta')) %>%
      mutate(component = 'fwd1h')
  ) %>%
  bind_rows(
    fwd10h_params %>%
      select(contains('beta')) %>%
      mutate(component = 'fwd10h')
  ) %>%
  bind_rows(
    fwd100h_params %>%
      select(contains('beta')) %>%
      mutate(component = 'fwd100h')
  ) %>%
  pivot_longer(cols = c('beta[1]', 'beta[2]', 'beta[3]'),
               names_to = c('parameter', 'group_id'),
               names_sep = '\\[') %>%
  mutate(group_id = gsub(x = group_id, pattern = '\\]', replacement = ''))%>%
  
  ggplot(aes(x = value, fill = group_id))+
    geom_density(alpha = 0.5)+
    facet_grid(component~., scales = 'free_y') +
    scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
    theme_minimal()

# beta varies by component

litterduff_params %>%
  select(contains('alpha')) %>%
  mutate(component = 'litterduff') %>%
  bind_rows(
    fwd1h_params %>%
      select(contains('alpha')) %>%
      mutate(component = 'fwd1h')
  ) %>%
  bind_rows(
    fwd10h_params %>%
      select(contains('alpha')) %>%
      mutate(component = 'fwd10h')
  ) %>%
  bind_rows(
    fwd100h_params %>%
      select(contains('alpha')) %>%
      mutate(component = 'fwd100h')
  ) %>%
  pivot_longer(cols = c('alpha[1]', 'alpha[2]', 'alpha[3]'),
               names_to = c('parameter', 'group_id'),
               names_sep = '\\[') %>%
  mutate(group_id = gsub(x = group_id, pattern = '\\]', replacement = ''))%>%
  
  ggplot(aes(x = value, fill = group_id))+
    geom_density(alpha = 0.5)+
    facet_grid(component~., scales = 'free_y') +
    scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
    theme_minimal()

# alpha seems to vary more by mortality class than by component, except that 
# L/D is different from fwd

litterduff_params %>%
  select(contains('rho')) %>%
  mutate(component = 'litterduff') %>%
  bind_rows(
    fwd1h_params %>%
      select(contains('rho')) %>%
      mutate(component = 'fwd1h')
  ) %>%
  bind_rows(
    fwd10h_params %>%
      select(contains('rho')) %>%
      mutate(component = 'fwd10h')
  ) %>%
  bind_rows(
    fwd100h_params %>%
      select(contains('rho')) %>%
      mutate(component = 'fwd100h')
  ) %>%
  pivot_longer(cols = c('rho[1]', 'rho[2]', 'rho[3]'),
               names_to = c('parameter', 'group_id'),
               names_sep = '\\[') %>%
  mutate(group_id = gsub(x = group_id, pattern = '\\]', replacement = ''))%>%
  
  ggplot(aes(x = value, fill = group_id))+
    geom_density(alpha = 0.5)+
    facet_grid(component~., scales = 'free_y') +
    scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
    theme_minimal()+
  scale_x_continuous(limits = c(0, 20))

# rho is different between L/D and FWD1h, and the posteriors for fwd10h and fwd100h 
# are mostly just the prior (as seen above)

fwd1h_params %>%
      select(contains('tau')) %>%
      mutate(component = 'fwd1h') %>%
  bind_rows(
    fwd10h_params %>%
      select(contains('tau')) %>%
      mutate(component = 'fwd10h')
  ) %>%
  bind_rows(
    fwd100h_params %>%
      select(contains('tau')) %>%
      mutate(component = 'fwd100h')
  ) %>%
  pivot_longer(cols = c('tau[1]', 'tau[2]', 'tau[3]'),
               names_to = c('parameter', 'group_id'),
               names_sep = '\\[') %>%
  mutate(group_id = gsub(x = group_id, pattern = '\\]', replacement = ''))%>%
  
  ggplot(aes(x = value, fill = group_id))+
    geom_density(alpha = 0.5)+
    facet_grid(component~., scales = 'free_y') +
    scale_fill_viridis_d(begin = 0.05, end = 0.85, option = 'C')+
    theme_minimal()

# tau very similar across all mortality classes for 1h fuels, classes diverge 
# for 10h fuels, tau seems generally higher but less certain fore 100h fuels


# takeaway: 
#   - realistic simulations clearly require correlation amongst the FWD classes, 
#   which is clear in the raw data. 
#   - 100h fwd, and to a lesser extent 10h fwd, are rare enough that the posterior 
#   distributions for the GP parameters are mostly informed by the prior, rather than 
#   the data
#   - there is mixed evidence about the GP parameters varying by size class. 
#   Some parameters vary more by mortality class than size class (alpha), and 
#  the data are not informative about some parameters for 10h and 100h fuels (rho), 
#  and the data are suggestive about differences among size classes for tau.
#  - given the 3 above points, I think the most realistic simulations will be 
# for a generic FWD model, lumping the tallies to include 1-100h and then assigning 
# particle diameters from the FWD diameters model.
