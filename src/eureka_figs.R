## title: figures for eureka seminar
## author: trond
## date: 2.3.2020

## house keeping
library(rgdal)
library(ggplot2)
library(data.table)
library(rAzureBatch)
library(doAzureParallel)
library(sf)
library(xtable)
library(gridExtra)

## login
system('az login')

## after filled in set credentials
setCredentials("credentials.json")

## useful functions
source('handy_functions.R')

## dependent (1) and independent (-1) variables
model_vars <- c('V1', 'AANT_INW', 'dist_rd', 'nb', 'cent_rd')

## download data

## input data
load('../output/sp_dt.Rdata')

## testing on a subset (comment out for full)
sp_dt <- subset(sp_dt, dist <= 35)
##sp_dt <- subset(sp_dt, code_orig == unique(sp_dt@data$code_orig)[c(1,2)])
##sp_dt <- subset(sp_dt, prov_nm %in% unique(sp_dt@data$prov_nm)[c(1:3)])

## test data
test_dt <- readRDS('../output/test_dt.rds')
test_dt <- subset(test_dt, dist <= 35  & code_orig %in% sp_dt@data$code_orig & code_orig != code_dest)

## download results from gwr estimation
gwr_poisson_out <- getJobResult("eureka-seminar-19")[[2]]
##gwr_poisson_out <- readRDS('../output/dummy_output.rds')[[2]]
    
## download results from ols run
gwr_ols_out <- getJobResult("eureka-seminar-19")[[1]]
##gwr_ols_out <- readRDS('../output/dummy_output.rds')[[1]]

## download output from calculation of prediction
cv_out <- getJobResult("calculate-errors-eureka-12")[[1]]

## download geographical data
epsg <- data.table(make_EPSG())
reg_data_dir <- '/home/trond/Documents/CLAiR-City_data/RegionalData/'
gem_nl <- list()
gem_nl[[2016]] <-read_sf(paste0(reg_data_dir, 'buurt_en_wijk_kaarten/2016/Uitvoer_shape/'), 'gem_2016')
gem_nl[[2016]] <- st_transform(gem_nl[[2016]], epsg[code == 7415, prj4]) # set RDS planar projection
gem_nl[[2016]]$code <- as.numeric(gsub('GM', '', gem_nl[[2016]]$GM_CODE))
gem_nl[[2016]] <- subset(gem_nl[[2016]], code %in% sp_dt@data$code_orig & WATER == 'NEE')
       
## out of sample predictions
oos_pred <- rbindlist(list(pred_fn(test_dt,
                                   data.table(code_orig = sp_dt@data$code_orig,
                                              code_dest = sp_dt@data$code_dest,
                                              gwr_ols_out[, 1:4]
                                              )
                                       )[,
                                         mod := 'ols'
                                         ],
                              pred_fn(test_dt,
                                       data.table(code_orig = sp_dt@data$code_orig,
                                                  code_dest = sp_dt@data$code_dest,
                                                  gwr_poisson_out[, 1:4]
                                                  )
                                       )[,
                                         mod := 'poisson'
                                         ]
                           ))
oos_pred <- merge(oos_pred,
                  sp_dt@data[,c('code_orig', 'code_dest', 'V1')],
                  by = c('code_orig', 'code_dest'), all.x = T)


## out of sample prediction errors
oos_error <- merge(
    oos_pred[,
           .(rmse = rmse(actual, prediction),
             mae = mae(actual, prediction)),
       by = mod
       ],
    oos_pred[,
             .(rmse(actual, prediction)/sum(actual),
               mean(abs((actual - prediction)/ mean(abs(V1 - mean(V1)))))),
           by = c('code_orig', 'mod')
           ][,
             .(srmse = mean(V1),
               mase = mean(V2)),
             by = mod
             ],
    by = 'mod'
)[,
  mod := c('OLS (bisquare, 340)', 'Poisson (boxcar, 136)')
  ]

## print latex table with errors
xtable(oos_error, digits = 4)

## some descriptive stats of explanatory variables
t(sapply(model_vars[-1],
       function(x) cbind(round(mean(sp_dt@data[, x]), 3), round(sd(sp_dt@data[, x]), 3)
                         )
       ))

## figures

## figures from cross validation
cv_out_long <- melt(cv_out, id.vars = c(1:3)
                    )[,
                      c('var', 'error', 'type') := tstrsplit(as.character(variable),
                                                             split = '_', fixed = T
                                                             )
                      ][,
                        type := gsub('ws', 'training error',
                                     gsub('oos', 'test error',
                                          type)
                                     )
                        ]

ggplot(cv_out_long[var == 'mean' & error == 'rmse']) +
    geom_errorbar(data = cv_out_long[error == 'rmse',
                                     .(value[var == 'mean'] + value[var == 'se'],
                                       value[var == 'mean'] - value[var == 'se']),
                                     by = c('idx', 'kernel', 'bw', 'type')],
                  aes(x = bw, ymin = V1, ymax = V2), width = 3, col = 'grey') +
    geom_line(aes(bw, value, group = 1)) +
    facet_grid(type ~ kernel, scales = 'free') +
    theme_bw() +
    xlab('bandwidth') +
    ylab('RMSE') +
    theme(strip.background = element_blank(), strip.placement = "outside")
ggsave('../figs/cv_oos_ws_comparison.png', width = 10, scale = 0.8)    

ggplot(cv_out_long[var == 'mean' & error == 'mae']) +
    geom_errorbar(data = cv_out_long[error == 'mae',
                                     .(value[var == 'mean'] + value[var == 'se'],
                                       value[var == 'mean'] - value[var == 'se']),
                                     by = c('idx', 'kernel', 'bw', 'type')],
                  aes(x = bw, ymin = V1, ymax = V2), width = 3, col = 'grey') +
    geom_line(aes(bw, value, group = 1)) +
    facet_grid(type ~ kernel, scales = 'free') +
    theme_bw() +
    xlab('bandwidth') +
    ylab('MAE') +
    theme(strip.background = element_blank(), strip.placement = "outside")
ggsave('../figs/cv_oos_ws_comparison2.png', width = 10, scale = 0.8)    

ggplot(cv_out, aes(bw, mean_rmse_oos, group = kernel)) +
    geom_line() +
    geom_errorbar(aes(ymin = mean_rmse_oos - se_rmse_oos, ymax = mean_rmse_oos + se_rmse_oos), width = 3, col = 'grey') +
    geom_hline(yintercept = min(cv_out$mean_rmse_oos), col = 'red') +
    #geom_line(data = cv_out, aes(bw, mean_rmse_ws, group = kernel), col = 'grey') +
    #geom_errorbar(aes(ymin = mean_rmse_ws - se_rmse_ws, ymax = mean_rmse_ws + se_rmse_ws), col = 'grey', width = 0.8) +
    facet_wrap(~kernel, ncol = 2, scales = 'free') + 
    theme_bw() +
    xlab('Bandwidth') +
    ylab('RMSE')
ggsave('../figs/cv_result.png', scale = 0.8)

ggplot(cv_out, aes(bw, mean_mae_oos, group = kernel)) +
    geom_line() +
    geom_errorbar(aes(ymin = mean_mae_oos - se_mae_oos, ymax = mean_mae_oos + se_mae_oos), width = 3, col = 'grey') +
    geom_hline(yintercept = min(cv_out$mean_mae_oos), col = 'red') +
    #geom_line(data = cv_out, aes(bw, mean_mae_ws, group = kernel), col = 'grey') +
    #geom_errorbar(aes(ymin = mean_mae_ws - se_mae_ws, ymax = mean_mae_ws + se_mae_ws), col = 'grey', width = 0.8) +
    facet_wrap(~kernel, ncol = 2, scales = 'free') + 
    theme_bw() +
    xlab('Bandwidth') +
    ylab('MAE')
ggsave('../figs/cv_result2.png', scale = 0.8)    

## map with coefficients
coeff_map <- merge(gem_nl[[2016]],
                   melt(unique(data.table(sp_dt@data$code_orig,
                                          gwr_poisson_out[, 1:4])
                               ),
                               id.vars = 'V1'),
                   by.x = 'code', by.y = 'V1', all.y = T
                   )

coeff_dict <- data.table(sort(unique(coeff_map$variable)),
                         c('population', 'distance', 'dwellings', 'centrality')
                         )


coeff_map <- merge(coeff_map,
                   coeff_dict,
                   by.x = 'variable',
                   by.y = 'V1',
                   all.x = T
                   )

coeff_plot <- lapply(unique(coeff_map$V2),
                     function(x) {
                         df <- subset(coeff_map, V2 == x)
                         mid_point <- ifelse(min(df$value) > 0,
                                             min(df$value),
                                      ifelse(max(df$value) < 0,
                                             max(df$value),
                                             0))
                         ggplot(df) +
                             geom_sf(aes(fill = value), size = 0) +
                             theme_void() +
                             coord_sf(datum=NA) +
                             scale_fill_gradient2(low = "red", midpoint = mid_point,
                                                  high = "forestgreen",
                                                  name = '') +
                             ggtitle(x)
                     })



## ggplot(coeff_map) +
##     geom_sf(aes(fill = value_cats), size = 0) +
##     facet_wrap(~V2) +
##     theme_void() +
##     coord_sf(datum=NA) +
##     ##scale_fill_distiller(palette = "RdBu", direction = 1, name = 'Coefficient') +
##     ##scale_fill_gradient2(low = "red", high = "forestgreen", name = '') +
##     scale_fill_brewer(palette = "RdBu", drop = FALSE, name = 'Parameter \n value')
##     ##theme(legend.position = 'bottom')
ggsave('../figs/coeffs.png',
       grid.arrange(coeff_plot[[1]], coeff_plot[[2]],
             coeff_plot[[3]], coeff_plot[[4]],
             ncol = 2)
       )

## density plot of t value
tvals <- data.table(unique(melt(gwr_poisson_out[, c(paste0(coeff_dict$V1, '_TV'))]))
                    )[,
                      variable := gsub('_TV', '', variable)
                      ]

tvals <- merge(tvals, coeff_dict, by.x = 'variable', by.y = 'V1', all.x = T)


## ggplot(tvals, aes(abs(value))) +
##     stat_ecdf(geom = 'step') +
##     facet_wrap(~V2, scales = 'free', ncol = 2) +
##     theme_bw() +
##     ylab('cumulative density') +
##     xlab('t-value (absolute)')

ggplot(tvals[,
      .(mean_val = mean(abs(value)),
        max_val = mean(abs(value)) + sd(value),
        min_val = mean(abs(value)) - sd(value)),
      by = V2]) +
    geom_errorbar(aes(V2, ymin = min_val, ymax = max_val), col = 'grey', width = 0.2) +     geom_point(aes(V2, mean_val)) +
    coord_flip() +
    xlab('parameter') +
    ylab('t-value (absolute)') + 
    theme_bw()
ggsave('../figs/tvals.png', scale = 0.5)


## predicted versus actual
p1 <- ggplot(oos_pred, aes(actual, prediction)) +
    geom_point(aes(col = mod)) +
    geom_abline(slope = 1) +
    scale_colour_brewer(palette = 'Set1', name = '', labels = c('OLS (bisquare, 340)', 'Poisson (boxcar, 136)')) +
    #facet_wrap(~mod) + 
    theme_bw() +
    ylab('predicted') + 
    theme(legend.position = c(0.3, 0.85),
          legend.background = element_rect(fill=NA))

p2 <- ggplot(rbindlist(list(oos_pred[,
                                     .(error = 'summed over destinations',
                                       actual = sum(actual),
                                       prediction = sum(prediction)),
                                     by = c('code_orig', 'mod')
                                     ]
                            ## oos_pred[,
                            ##          .(error = 'summed over origins',
                            ##            actual = sum(actual),
                            ##            prediction = sum(prediction)),
                            ##          by = c('code_dest', 'mod')
                            ##          ]
                            )
                       ), aes(actual, prediction)) +
    geom_point(aes(col = mod)) +
    geom_abline(slope = 1) +
    ##facet_wrap(~error, ncol =1) +
    theme_bw() +
    ylab('predicted') + 
    scale_colour_brewer(palette = 'Set1', name = '') +
    guides(col = FALSE) +
    ggtitle('Moves summed over destinations')
    ##theme(legend.position = c(0.2, 0.85), legend.background = element_rect(fill=NA))


p3 <- ggplot(merge(oos_pred[code_orig == 363][V1 > quantile(V1, probs = seq(0, 1, 0.1))[10], ],
      data.table(st_drop_geometry((gem_nl[[2016]][,c('code', 'GM_NAAM')]))),
      by.x = 'code_dest', by.y = 'code'),
      aes(GM_NAAM)) +
    geom_point(aes(y = actual), col = 'grey') +
    geom_point(aes(y = prediction, col = mod)) +
    scale_colour_brewer(palette = 'Set1', name = '') +
    guides(col = FALSE) +
    ylab('') +
    xlab('') +
    ggtitle('Moves from Amsterdam') + 
    theme_bw() +
    #coord_flip() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_comb <- grid.arrange(p1, grid.arrange(p2, p3), ncol = 2)

ggsave('../figs/actual_predicted.png', p_comb,  width = 10, height = 7, scale = 0.8)

## plot mase on map
ggplot(merge(gem_nl[[2016]],
             oos_pred[mod == 'poisson',
                      mean(abs((actual - prediction)/ mean(abs(V1 - mean(V1))))),
                      by = code_orig
                      ],
             by.x = 'code',
             by.y = 'code_orig')) +
    geom_sf(aes(fill = V1), size = 0) +
    scale_fill_gradient2(low = "forestgreen", midpoint = 1, high = "red", name = 'MASE') +
    theme_void() +
    coord_sf(datum=NA) +
    theme(legend.position = 'bottom') +
    ggtitle('Poisson (boxcar, 136)')
ggsave('../figs/mase.png', scale = 0.8)    
    

## plot distance decay function
ggplot(unique(gwr_poisson_out[, c(2L, 4L)])) +
    geom_point(aes() +
    geom_smooth(method = 'lm') + 
    theme_bw()

## plot of variables
var_dt <- merge(gem_nl[[2016]],
                melt(unique(data.table(sp_dt@data[, c('code_dest', model_vars[-1])]
                                )[,
                                  c('population', 'distance', 'dwellings', 'centrality') := lapply(.SD, mean),
                                  by = code_dest,
                                  .SDcols = c(model_vars[-1])
                                  ][,
                                    .(code_dest, population, distance, dwellings, centrality)
                                    ]),
                     id.vars = 'code_dest'),
                by.x = 'code', by.y = 'code_dest'
                )

var_plot <- lapply(unique(var_dt$variable),
                     function(x) {
                         df <- subset(var_dt, variable == x)
                         mid_point <- ifelse(x == 'dwellings', 1, mean(df$value))
                         ggplot(df) +
                             geom_sf(aes(fill = value), size = 0) +
                             theme_void() +
                             coord_sf(datum=NA) +
                             scale_fill_gradient2(low = "red", midpoint = mid_point,
                                                  high = "forestgreen",
                                                  name = '') +
                             ggtitle(x)
                     })

      
ggsave('../figs/vars.png',
       grid.arrange(var_plot[[1]], var_plot[[2]],
             var_plot[[3]], var_plot[[4]],
             ncol = 2)
       )

## plot arrivals and departures
arrivals <- data.table(sp_dt@data[,c('code_dest', 'V1')]
                         )[,
                           .(arrivals = sum(V1)),
                           by=code_dest
                           ]
departures <- data.table(sp_dt@data[,c('code_orig', 'V1')]
                       )[,
                         .(departures = sum(V1)),
                         by=code_orig
                         ]

setkey(arrivals, code_dest)
setkey(departures, code_orig)
tmp1 <- arrivals[departures, ]

## subset of destinations from Amsterdam
adam_arrivals <- data.table(subset(sp_dt@data[,c('code_orig', 'code_dest', 'V1')],
                                   code_orig == 363)
                            )[,
                              .(`arrivals from Amsterdam` = sum(V1)),
                              by=code_dest
                              ]

adam_departures <- data.table(subset(sp_dt@data[,c('code_orig', 'code_dest', 'V1')],
                                   code_dest == 363)
                            )[,
                              .(`departures to Amsterdam` = sum(V1)),
                              by=code_orig
                              ]

setkey(adam_arrivals, code_dest)
setkey(adam_departures, code_orig)
tmp2 <- adam_arrivals[adam_departures, ]


p1 <- ggplot(merge(gem_nl[[2016]],
      melt(tmp1, id.vars = 'code_dest'),
      by.x = 'code', by.y = 'code_dest')) +
    geom_sf(aes(fill = value), colour = NA) +
    facet_wrap(~variable, ncol = 1) +
    coord_sf(datum=NA) +
    scale_fill_distiller(palette = 'Greens', direction = 1, name = '') +
    theme_void() +
    theme(legend.position = 'bottom') +
    guides(fill = guide_colourbar(barwidth = 10))

p2 <- ggplot(merge(gem_nl[[2016]],
                   melt(rbindlist(list(tmp2, data.frame(363, NA, NA))),
                        id.vars = 'code_dest'
                        )[,
                          adam := ifelse(code_dest == 363,
                                         'yes',
                                         NA)
                          ],
      by.x = 'code', by.y = 'code_dest')) +
    geom_sf(aes(fill = value, colour = adam)) +
    facet_wrap(~variable, ncol = 1) +
    coord_sf(datum=NA) +
    scale_fill_distiller(palette = 'Greens', direction = 1, na.value = 'white', name = '') +
    scale_colour_manual(values = 'red') + 
    theme_void() +
    theme(legend.position = 'bottom') +
    guides(colour = F)

grid.arrange(p1, p2, ncol = 2)

ggsave('../figs/arrivals_and_departures.png',
       grid.arrange(p1, p2, ncol = 2)
       )

## distance decay versus centrality
dist_decay_dt <- data.table(sp_dt@data[, c('code_dest', '')])

p1 <- ggplot(data.table(sp_dt@data$cent_rd, gwr_poisson_out$`log(dist_rd)`),
       aes(V1, abs(V2))) +
    geom_point(col = 'grey') +
    stat_smooth(method = 'lm') +
    xlab('centrality of destination') +
    ylab(expression(paste('|', gamma, '|')))+
    theme_bw() +
    ggtitle('distance')

p2 <- ggplot(data.table(sp_dt@data$cent_rd, gwr_poisson_out$`log(AANT_INW)`),
       aes(V1, V2)) +
    geom_point(col = 'grey') +
    stat_smooth(method = 'lm') +
    xlab('centrality of destination') +
    ylab(expression(beta['population'])) +
    theme_bw() +
    ggtitle('population')

ggsave('../figs/distance_decay.png', grid.arrange(p1, p2, ncol = 2), height = 7, width = 12, scale = 0.5)
           






