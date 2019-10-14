## title: estimation of short-distance migration 2008
## author: trond
## date: 03.12.2018

library(readxl)
library(raster)
library(data.table)
library(rgdal)
library(maptools)
library(rgeos)
library(spdep)
library(ggplot2)
library(cbsodataR)
library(parallel)
library(doParallel)
library(rgdal)
library(viridis)
library(gridExtra)
library(devtools)
library(sf)
library(httr)
library(jsonlite)
library(GWmodel)
library(corrplot)
library(caret)

## some nice refs
# https://cran.r-project.org/web/packages/spgwr/index.html
# https://rpubs.com/chrisbrunsdon/101305
# https://rpubs.com/adam_dennett/376877
# http://rspatial.org/analysis/rst/6-local_regression.html
# http://gwr.maynoothuniversity.ie/
# https://rstudio-pubs-static.s3.amazonaws.com/44975_0342ec49f925426fa16ebcdc28210118.hTm
# https://www.researchgate.net/publication/261960122_Geographically_weighted_regression_with_a_non-Euclidean_distance_metric_A_case_study_using_hedonic_house_price_data
# https://www.tandfonline.com/doi/full/10.1080/13658816.2013.865739?src=recsys

## functions
source('handy_functions.R')

## read 2006 - 2008 file with data
#source('short_distance_estimation_08.R')

## read data on house prices and costs
cbs_dt <- data.table(cbs_get_toc(Language="nl"))

cbs_dt[grepl('Voorraad', Title), .(Identifier, Title, ShortTitle, Modified)]

huur_dt <- data.table(cbs_get_data('83368NED'))
woz_dt <- data.table(cbs_get_data('37610'))
sales_dt <- data.table(cbs_get_data('83625NED'))
wtype_dt <- data.table(cbs_get_data('82900NED'))


## fix municipality code and set numeric
sales_dt[grepl('GM', RegioS),
         code := as.numeric(gsub('GM', '', as.character(RegioS)))
         ][,
           GemiddeldeVerkoopprijs_1 := as.numeric(GemiddeldeVerkoopprijs_1)
           ]

woz_dt[grepl('GM', RegioS),
         code := as.numeric(gsub('GM', '', as.character(RegioS)))
         ][,
           GemiddeldeWoningwaarde_4 := as.numeric(GemiddeldeWoningwaarde_4)
           ]

changecols <- names(wtype_dt)[4:ncol(wtype_dt)]

wtype_dt[grepl('GM', RegioS),
         code := as.numeric(gsub('GM', '', as.character(RegioS)))
         ][,
           (changecols) := lapply(.SD, as.numeric),
           .SDcols = changecols
           ]

## Berg en Dal (GM1945) in 2015: GM0241
## Gooise Meren (GM1942) in 2015: Ontstaan per 01-01-2016, ontvangen van:
##- Bussum (GM0381),
##... 808 hectare met 15681 woningen en 33122 inwoners
##- Muiden (GM0424),
##... 1207 hectare met 2781 woningen en 6206 inwoners
##- Naarden (GM0425),
##... 2142 hectare met 7675 woningen en 17406 inwoners
##- Weesp (GM0457),
##... 4 hectare met 0 woningen en 0 inwoners

## municipality corrections
sales_dt[Perioden == '2015JJ00',
         ':='(verf_1945 = GemiddeldeVerkoopprijs_1[RegioS == 'GM0241'],
              verf_1942 = (15681*GemiddeldeVerkoopprijs_1[RegioS == 'GM0381'] + 2781*GemiddeldeVerkoopprijs_1[RegioS == 'GM0424'] + 7675*GemiddeldeVerkoopprijs_1[RegioS == 'GM0425'])/(15681 + 2781 + 7675)
              )
         ][Perioden == '2015JJ00' & RegioS == 'GM1945',
           GemiddeldeVerkoopprijs_1 := verf_1945
           ][Perioden == '2015JJ00' & RegioS == 'GM1942',
             GemiddeldeVerkoopprijs_1 := verf_1942
             ][,
               c('verf_1942', 'verf_1945') := NULL
               ][,
                 price_index := GemiddeldeVerkoopprijs_1 / GemiddeldeVerkoopprijs_1[Perioden == '2015JJ00'],
                 by = RegioS
                 ]

## municipality corrections
woz_dt[Perioden == '2015JJ00',
       ':=' (verf_1945 = GemiddeldeWoningwaarde_4[RegioS =='GM0241'],
             verf_1940 = GemiddeldeWoningwaarde_4[RegioS == 'GM1921'],
             verf_1942 = (15681*GemiddeldeWoningwaarde_4[RegioS == 'GM0381'] + 2781*GemiddeldeWoningwaarde_4[RegioS == 'GM0424'] + 7675*GemiddeldeWoningwaarde_4[RegioS == 'GM0425'])/(15681 + 2781 + 7675)
             )
       ][Perioden == '2015JJ00' & RegioS == 'GM1945',
         GemiddeldeWoningwaarde_4 := verf_1945
           ][Perioden == '2015JJ00' & RegioS == 'GM1942',
             GemiddeldeWoningwaarde_4 := verf_1942
             ][Perioden == '2015JJ00' & RegioS == 'GM1940',
               GemiddeldeWoningwaarde_4 := verf_1940
               ][,
                 c('verf_1940', 'verf_1942', 'verf_1945') := NULL
                 ][,
                   woz_index := GemiddeldeWoningwaarde_4 / GemiddeldeWoningwaarde_4[Perioden == '2015JJ00'],
                   by = RegioS
                   ]

## plot of sales and woz
ggplot(sales_dt[Perioden == '2016JJ00'],
       aes(price_index)) +
    geom_density() +
    geom_density(data = woz_dt[Perioden == '2016JJ00'], aes(woz_index), col = 'blue') + 
    theme_bw()

## migration data
load(paste0(getwd(), '/../data/reg_bila_dt_2016_code.Rdata'))

## building stock
#wn_dt <- data.table(cbs_get_data('81955NED'))
#save(wn_dt, file = '../data/wn_dt.RData')
load('../data/wn_dt.RData')

## municipality list
gem_list_2016 <- data.table(cbs_get_data('83287NED'))[, code := as.numeric(gsub('GM', '', as.character(RegioS)))]

gem_list_2015 <- data.table(cbs_get_data('82949NED'))[, code := as.numeric(gsub('GM', '', as.character(RegioS)))]

## spatial data
epsg <- data.table(make_EPSG())
epsg[grep('Amersfoort', note, ignore.case = T), ]

reg_data_dir <- '/home/trond/Documents/CLAiR-City_data/RegionalData/'

gem_nl <- list()
vierkant <- list()
for (i in c(2016)) {
    gem_nl[[i]] <- read_sf(paste0(reg_data_dir, 'buurt_en_wijk_kaarten/2016/Uitvoer_shape/'), 'gem_2016')
    vierkant[[i]] <-  read_sf(paste0(reg_data_dir, 'CBS_vierkant/2016/'), 'CBSvierkant100m_2016_v1')
    gem_nl[[i]] <- st_transform(gem_nl[[i]], epsg[code == 7415, prj4]) # set RDS planar projection
    vierkant[[i]] <- st_transform(vierkant[[i]], epsg[code == 7415, prj4])  
}

## clean up mun code
gem_nl[[2016]]$code <- as.numeric(gsub('GM', '', as.character(gem_nl[[2016]]$GM_CODE)))

## create population-weighted centroids

## find centroids in each square and overlap with municipality
vk_cents <- data.table(st_join(st_centroid(vierkant[[2016]][,1:2]),
                               subset(gem_nl[[2016]], WATER == 'NEE')
                               )
                       )

## clean up a bit
vk_cents[,
         ':=' (x = as.numeric(geometry[[1]][1]),
               y = as.numeric(geometry[[1]][2])),
         by = C28992R100
         ]

## create population weighted centroids
pop_w_cents <- vk_cents[INWONER >0,
                        lapply(.SD, weighted.mean, w = INWONER),
                        by = GM_CODE,
                        .SDcols = c('x', 'y')
                        ][,
                          GM_CODE := as.numeric(gsub('GM', '', GM_CODE))
                          ][GM_CODE != 9999,
                            ]

## check whether the points are inside municipality limits
pop_w_cents_sf <- st_as_sf(pop_w_cents,
                           coords = c('x', 'y'),
                           crs = epsg[code == 28992, prj4],
                           agr = "constant")

cent_gem_int <- data.table(mun_id = as.numeric(st_intersects(pop_w_cents_sf, gem_nl[[2016]])))[, row.id := 1:.N]

missing_cents <- pop_w_cents[which(!cent_gem_int$mun_id %in% seq(1, nrow(pop_w_cents))), GM_CODE]

## allocate randomly those that are missing
if (length(missing_cents) > 0 ) {
  set.seed(123)
  missing_cents <- cbind(missing_cents,
                         data.frame(spsample(as(gem_nl[[2016]][as.numeric(gsub('GM', '', gem_nl[[2016]]$GM_CODE)) %in% missing_cents, ],
                                                'Spatial'),
                                             1, 'random')
                                    ))
  pop_w_cents <- rbindlist(list(pop_w_cents[!GM_CODE %in% missing_cents[,1]], missing_cents))
} 

## calculate distance between municipalities
cent_grid <- data.table(expand.grid(pop_w_cents$GM_CODE, pop_w_cents$GM_CODE))
setkey(cent_grid, Var1)
setkey(pop_w_cents, GM_CODE)
cent_grid[pop_w_cents, ':='(o_x = x, o_y = y)]
setkey(cent_grid, Var2)
cent_grid[pop_w_cents, ':='(d_x = x, d_y = y)]
cent_grid[, dist := pointDistance
(data.frame(o_x, o_y), data.frame(d_x, d_y), lonlat = F)/1000][, dist_test := sqrt((d_x - o_x)^2 + (d_y - o_y)^2) / 1000]


## distance over the road
if(!file.exists('../data/road_dist.RData')) {
    source('calculate_distance_osrm.R')
} else {
    load('../data/road_dist.RData')
}

## merge 
setkey(cent_grid, Var1, Var2)
setkey(results, o, d)
cent_grid[results, dist_rd := i.dist/1000]
setkey(cent_grid, Var2, Var1)
cent_grid[results, dist_rd := i.dist/1000]
cent_grid[Var1 == Var2, dist_rd := 0]

## create data table with needed variables: choose average 2014 to 2016 or only 2016
inp_dt <- merge(reg_bila_dt_2016_code[,
                                      sum(value),
                                      by = c(names(reg_bila_dt_2016_code)[c(1:2, 4)])
                                      ][, ## comment if not average over the years
                                        mean(V1),
                                        by = c(names(reg_bila_dt_2016_code)[1:2])
                                        ],
                cent_grid[,.(Var1, Var2, dist, dist_rd)],
                by.x = c('code_orig', 'code_dest'),
                by.y = c('Var1', 'Var2'),
                all.y = T)[is.na(V1),
                           V1 := 0
                           ][is.na(dist_rd),
                             dist_rd := 0
                             #][year == 2016, ## comment if average 
                               ][,
                                 year := NULL 
                                 ]


    
## add inhabitants at destination
inp_dt <- merge(inp_dt,
                data.table(subset(gem_nl[[2016]], WATER == 'NEE')[,c('code', 'AANT_INW', 'OPP_LAND')])[, -'geometry', with = F],
                by.x = 'code_dest', by.y = 'code', all.x = T)


## add percentage owner occupied and social rent
inp_dt <- merge(inp_dt,
                wtype_dt[grepl('2016', Perioden) & StatusVanBewoning == 'T001235' & code %in% gem_nl[[2016]]$code,
                         .(p_koop = Koopwoningen_2/ TotaleWoningvoorraad_1,
                           p_shuur = EigendomWoningcorporatie_4 / TotaleWoningvoorraad_1, code)
                         ],
                by.x = 'code_dest', by.y = 'code', all.x = T
                )
                
## add delta woningbouw voorrad/(voorraad - (nieuwbouw - sloop)) at destination
inp_dt <- merge(inp_dt,
                wn_dt[Gebruiksfunctie == 'A045364' & grepl('GM', RegioS) & Perioden %in% paste0(2014:2016, 'JJ00'), ## change period to include more years
                      ][RegioS %in% wn_dt[Perioden == '2016JJ00' & !is.na(BeginstandVoorraad_1), unique(RegioS)],
                        BeginstandVoorraad_1[Perioden == '2016JJ00']/(BeginstandVoorraad_1[Perioden == '2016JJ00'] - sum((Nieuwbouw_2 - Sloop_4), na.rm = T)), 
                        by = RegioS
                        ][,
                          .(code = as.numeric(as.character(gsub('GM', '', RegioS))), nb = V1)
                          ],
                by.x = 'code_dest', by.y = 'code', all.x = T)


## corop-level totale wonlasten
woonlasten <- merge(gem_list_2016[,.(code, trimws(as.character(Code_8)))],
                    huur_dt[Perioden == '2015JJ00' & grepl('CR', RegioS) & grepl('T001096', EigenaarOfHuurder) & grepl('MW00000', Marges) & grepl('T009004', KenmerkenWoningen),
                            .(trimws(as.character(RegioS)), TotaalWoonlasten_1)
                            ],
                    by.x = 'V2',
                    by.y = 'V1',
                    all.x = T
)

setkey(woonlasten, code)
setkey(inp_dt, code_dest)
inp_dt[woonlasten, h_costs := i.TotaalWoonlasten_1]

## house prices
setkey(sales_dt, code)
setkey(inp_dt, code_dest)
inp_dt[sales_dt[!is.na(code) & grepl('2016JJ00', Perioden)],
       ':=' (hp_dest = i.GemiddeldeVerkoopprijs_1, hp_index_dest = i.price_index)
       ]
setkey(inp_dt, code_orig)
inp_dt[sales_dt[!is.na(code) & grepl('2016JJ00', Perioden)],
       ':=' (hp_orig = i.GemiddeldeVerkoopprijs_1,
             hp_index_orig = i.price_index)
       ][!is.na(hp_index_dest),
         ':='(
             hp_diff = hp_orig / hp_dest,
             hp_index_diff = hp_index_orig / hp_index_dest)
         ]

## woz value
setkey(woz_dt, code)
setkey(inp_dt, code_dest)
inp_dt[woz_dt[!is.na(code) & grepl('2016JJ00', Perioden)], ':=' (woz = i.GemiddeldeWoningwaarde_4, woz_index = i.woz_index)]

## plot of woningbouw
##ggplot(merge(subset(gem_nl[[2016]], WATER == 'NEE'), inp_dt[, unique(nb), by = code_dest], by.x = 'code', by.y = 'code_dest', all.x = T, all.y = F)) +
##    geom_sf(aes(fill = V1), size = 0) +
##    theme_void()

## increase population of diemen with 20%
#inp_dt[code_dest == 384, AANT_INW := AANT_INW * 1.2]

## centrality measure
inp_dt[dist > 0, ':=' (cent = AANT_INW/dist, cent_rd = AANT_INW/dist_rd)]

if(!file.exists('../output/centrality.RData')) {
    ## calculate centrality measure
    no_cores <- detectCores() - 1
    ## Initiate cluster
    cl <- makeCluster(no_cores, type="FORK")
    registerDoParallel(cl)
    setkey(inp_dt, code_orig)
    ## calculate centrality
    system.time(
        centrality <- foreach(x=inp_dt[,unique(code_orig)],
                              .packages = 'data.table') %dopar% {
                                  dt <- rbindlist(lapply(inp_dt[code_orig == x & code_dest != x , unique(code_dest)],
                                                         function(y)
                                                             inp_dt[code_orig == y  & !code_dest %in% c(x,y),
                                                                    .(x,y,
                                                                      sum(cent, na.rm = T),
                                                                      sum(cent_rd, na.rm = T)
                                                                      )
                                                                    ]
                                                         )
                                                  )
                              }
    )
    centrality <- rbindlist(centrality)
    setnames(centrality, 3:4, c('centrality', 'centrality_rd'))
    save(centrality, file = '../output/centrality.RData')
} else {
    load('../output/centrality.RData')
}

centrality[x == 363 & y %in% c(362, 384), ]
#centrality2[x == 363 & y %in% c(362, 384), ]

setkey(inp_dt, code_orig, code_dest)
setkey(centrality, x, y)
inp_dt[centrality, ':=' (cent = i.centrality, cent_rd = i.centrality_rd)]

#plot_centrality(363)    

# add province name 
setkey(inp_dt, code_orig)
setkey(gem_list_2016, code)
inp_dt[gem_list_2016,
       ':=' (prov_code = trimws(Code_34), prov_nm = trimws(Naam_35))
       ]

## plot centroids South Holland
#zh_mun <- inp_dt[grepl('Zuid', prov_nm), unique(code_orig)]

#ggplot(subset(gem_nl[[2016]], WATER == 'NEE' & code %in% zh_mun) ) +
#    geom_sf(col = 'grey', alpha = 0) +
#    geom_point(data = data.frame(st_coordinates(st_centroid(subset(gem_nl[[2016]], WATER == 'NEE' & code %in% zh_mun)))), aes(X, Y), col = 'red') +
#    geom_point(data = pop_w_cents[GM_CODE %in% zh_mun], aes(x,y), col = 'blue') + 
#    theme_void() +
#    coord_sf(datum=NA) +
#    annotate('text', x = 5.05, y = 52.25, label = 'weighted centroid', col = 'blue') +
#    annotate('text', x = 5, y = 52.3, label = 'centroid', col = 'red')
#ggsave('../figs/centroids_zh.pdf')

## create independent variables
model_vars <- c('V1', 'AANT_INW', 'dist_rd', 'nb', 'cent', 'cent_rd'
#                'hp_diff', 'hp_index_dest', 'p_koop', 'OPP_LAND'
                )

## alterations to dependent variable: integerise or replace zero flows with 1
short_distance  <- TRUE
if (short_distance) {
    inp_dt[code_orig != code_dest,
           V1_int := int_trs(V1),
           by = code_orig
           ][code_orig != code_dest & V1 == 0,
             V1 := 0.1
             ]
}

## calculate ind and dep var: log(x) - mean(log(x))
inp_dt[code_orig != code_dest & dist <= 35,
       paste0('V_', model_vars) := lapply(.SD,
                                          function(x) log(x) - mean(log(x))
                                          ),
       .SDcols = model_vars,
       by = code_orig
       ]

## plot correlation between ind vars
h_vars <- c('nb', 'h_costs', 'hp_dest',
            'hp_index_dest',  'woz',
            'woz_index', 'hp_diff', 'hp_index_diff',
            'p_koop', 'p_shuur')

corrplot(inp_dt[code_orig != code_dest,
                lapply(.SD, unique),
                by = c('code_dest'),
                .SDcols = h_vars
                ][,
                  round(cor(.SD), 4),
                  .SDcols = h_vars
                  ],
         method = 'circle')


## out-of-sample srmse

## specify resampling
fitControl <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 10
    )

lapply(unique(inp_dt$prov_nm), run_prov_regression)

## run gravity regression: lm
out_dt_2016  <- gravity_reg(inp_dt[code_orig != code_dest & dist <= 35],
                            paste0('V_', model_vars[-1]),
                            'V_V1')

## run gravity regression: poisson
out_dt_glm <- gravity_glm_reg(inp_dt[code_orig != code_dest & dist <= 35],
                              model_vars[-1],
                              'V1_int')

srmse_glm <- out_dt_glm[[1]][,
                             rmse(actual, predicted) / sum(actual),
                             by = prov_nm
                             ]


ggplot(out_dt_glm[[1]],
       aes(actual, predicted)) +
    geom_point() +
    facet_wrap(~prov_nm, scales = 'free') + 
    theme_bw()

## make consistent with 2016 parameters: make sure these are correct!!
greek_par_nm <- c('population(alpha)', 'centrality(gamma)', 'distance(beta)',
                  'pricerelative(zeta)','pricegrowth(eta)', 'dwellings(delta)')
names(greek_par_nm) <- sort(out_dt_2016[, unique(parameter)])

#lapply(seq(length(greek_par_nm)),
#       function(x) {
#           out_dt_2016[parameter == names(greek_par_nm)[x],
#                      parameter := gsub(names(greek_par_nm)[x],
#                                         greek_par_nm[x],
#                                         parameter)
#                       ]           
#       }
#       )

## compare 2008 and 2016
#ggplot(rbindlist(list(out_dt_2008[variable == 'Estimate'],
#               out_dt_2016[variable == 'Estimate'])
#          )[,
#            year := rep(c(2008, 2016), each = nrow(out_dt_2016[variable == 'Estimate']))
#            ], aes(prov_nm, V1)) +
# geom_ribbon(aes(ymin = est_lo, ymax = est_up, col = factor(year))) +

p1 <- ggplot(out_dt_2016[variable == 'Estimate'], aes(prov_nm, V1)) +
    geom_ribbon(aes(ymin = est_lo, ymax = est_up), col = 'grey70') +
    geom_point() +
    coord_flip() +
    facet_wrap(~parameter, scales = 'free', labeller = label_parsed) +
    theme_bw() +
 #   scale_colour_brewer(name = 'year', palette = 'Set1') + 
    ylab('value') +
    xlab('province')
#ggsave('../figs/est_08_16.pdf')

p2 <- ggplot(out_dt_glm[[2]][!grepl('factor', parameter) & variable == 'Estimate'], aes(prov_nm, value)) +
    geom_ribbon(aes(ymin = est_lo, ymax = est_up), col = 'grey70') +
    geom_point() +
    coord_flip() +
    facet_wrap(~parameter, scales = 'free', labeller = label_parsed) +
    theme_bw() +
 #   scale_colour_brewer(name = 'year', palette = 'Set1') + 
    ylab('value') +
    xlab('province')
#ggsave('../figs/est_08_16.pdf')

grid.arrange(p1, p2)

## gather province-level parameters: ols
out_dt_pars <- dcast(out_dt_2016[variable == 'Estimate',
                            .(prov_nm, V1, parameter)
                            ],
                     'prov_nm~parameter',
                     value.var = 'V1')

## mle
out_dt_pars <- dcast(out_dt_glm[[2]][!grepl('factor', parameter) & variable == 'Estimate',
                            .(prov_nm, value, parameter)
                            ],
                     'prov_nm~parameter',
                     value.var = 'value')

setnames(out_dt_pars, 2:ncol(out_dt_pars), paste0('V_', sort(model_vars[-1])))

## make predictions and calculate srmse
setkey(inp_dt, prov_nm)
setkey(out_dt_pars, prov_nm)

pred_prov_dt <- copy(inp_dt[out_dt_pars,
                            ][code_orig != code_dest & dist <= 35 ,
                              pred_denom := sum(eval(parse(text = gravity_vars(model_vars[-1],mod='ols'))), na.rm = T),
                              by = code_orig,
                              ][code_orig != code_dest & dist <= 35,
                                prediction := sum(V1, na.rm = T)*eval(parse(text = gravity_vars(model_vars[-1],mod='ols')))/pred_denom,
                                by = code_orig
                                ][code_orig != code_dest & dist <= 35,
                                  srmse := (1/sum(V1))*rmse(prediction, V1),
                                  by = prov_nm])

save(out_dt_pars, pred_prov_dt, file = '../data/pred_province_sd_mle.Rdata')

## srmse province level
ggplot(rbindlist(list(pred_prov_dt[!is.na(prediction),
                            .(unique(srmse), model = 'OLS'),
                            by=prov_nm],
               cbind(srmse_glm, rep('MLE', nrow(srmse_glm))))
               ), aes(prov_nm, V1)) +
    geom_col(aes(fill = model), position = 'dodge2') +
    coord_flip() +
    ylab('SRMSE') +
    xlab('Province') + 
    theme_bw()

#ggsave('../figs/srmse_models.png')    
             
## SRMSE
p2 <- ggplot(pred_prov_dt[!is.na(prediction),unique(srmse), by=prov_nm], aes(prov_nm, V1)) +
    geom_col() + 
    coord_flip() +
    ylab('SRMSE') +
    xlab('Province') + 
    theme_bw()

## geographixally weighted regression: save spatial point df
sp_dt <- merge(as.data.frame(inp_dt[code_orig != code_dest]),
               pop_w_cents,
               by.x = 'code_orig',
               by.y = 'GM_CODE',
               )

coordinates(sp_dt)  <- ~x+y
sp_dt@proj4string <- CRS(epsg[code == 7415, prj4])

save(sp_dt, file = '../output/sp_dt.Rdata')

## download results
#load('../data/gwr_results.Rdata')


## end script
