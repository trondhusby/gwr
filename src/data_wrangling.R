## title: data wrangling
## author: trond
## date: 03.12.2018

##
## house keeping
##

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
library(gridExtra)
library(sf)

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

##
## read data
##

## read data on house prices and costs
cbs_dt <- data.table(cbs_get_toc(Language="nl"))
##cbs_dt[grepl('Voorraad', Title), .(Identifier, Title, ShortTitle, Modified)]
##cbs_dt[grepl('Gebieden', Title), .(Identifier, Title, ShortTitle, Modified)]

## migration data
load('../data/reg_bila_dt_2016_code.Rdata')
if(!file.exists('../data/mig_dt_statline.rds')){
    mig_dt <- data.table(cbs_get_data('81734NED',
                       RegioVanVestiging = has_substring('GM'),
                       RegioVanVertrek = has_substring('GM')
                       ))
    saveRDS(mig_dt[!is.na(TussenGemeentenVerhuisdePersonen_1)],
            '../data/mig_dt_statline.rds')
} else {
    mig_dt <- readRDS('../data/mig_dt_statline.rds')
}

## fix orig and destination variable names
mig_dt[,
       c('code_orig', 'code_dest') := lapply(.SD,
                                             function(x) {
                                                 as.numeric(gsub('GM', '', x))
                                             }),
       .SDcols = c(names(mig_dt)[2:1])
       ][,
         year := as.numeric(substr(Perioden, 1, 4))
         ]

## population data
bev_dt <- data.table(cbs_get_data('37259ned',
                                  Geslacht="T001038",
                                  Perioden=c(paste0(2012:2018, 'JJ00'))
                                  )
                     )[grepl('GM', RegioS),
                       .(RegioS, Perioden,
                         BevolkingOp1Januari_1 = as.numeric(BevolkingOp1Januari_1),
                         BevolkingOp31December_25 = as.numeric(BevolkingOp31December_25),
                         code = as.numeric(gsub('GM', '', as.character(RegioS))),
                         year = as.numeric(substr(Perioden, 1, 4))
                         )
                       ]

bev_dt[year < 2018, sum(BevolkingOp1Januari_1, na.rm = T), by= year]

## building stock: NB measured at the end of the year
wn_dt <- data.table(cbs_get_data('81955NED',
                                 Gebruiksfunctie='A045364',
                                 Perioden=c(paste0(2012:2018, 'JJ00'))
                                 )
                    )[grepl('GM', RegioS),
                      .(RegioS,
                        Perioden,
                        BeginstandVoorraad_1 = as.numeric(BeginstandVoorraad_1),
                        EindstandVoorraad_8 = as.numeric(EindstandVoorraad_8),
                        Nieuwbouw_2 = as.numeric(Nieuwbouw_2),
                        Sloop_4 = as.numeric(Sloop_4),
                        #OverigeOnttrekking_5 = as.numeric(OverigeOnttrekking_5),
                        code = as.numeric(gsub('GM', '', as.character(RegioS))),
                        year = as.numeric(substr(Perioden, 1, 4))
                        )
                      ]

## fix de friese meren
friese_meren <- wn_dt[year == 2015 & code %in% c(1940, 1921),
                      .(code=1921,
                        year=2015,
                        BeginstandVoorraad_1 = BeginstandVoorraad_1[code == 1921],
                        EindstandVoorraad_8 = EindstandVoorraad_8[code == 1940],
                        Nieuwbouw_2 = sum(Nieuwbouw_2),
                        Sloop_4 = sum(Sloop_4)
                        )]

wn_dt <- rbindlist(list(wn_dt[!(year == 2015 & code %in% c(1940, 1921)),
                              c(names(friese_meren)),
                              with = FALSE
                              ],
                        friese_meren)
                   )
setkey(wn_dt, code, year)

## spatial data
epsg <- data.table(make_EPSG())
##epsg[grep('Amersfoort', note, ignore.case = T), ]

reg_data_dir <- '/home/trond/Documents/CLAiR-City_data/RegionalData/'

gem_nl <- read_sf(paste0(reg_data_dir, 'buurt_en_wijk_kaarten/2018/'), 'gemeente_2018_v2')
gem_nl <- st_transform(gem_nl, epsg[code == 7415, prj4]) # set RDS planar projection
## clean up mun code
gem_nl$code <- as.numeric(gsub('GM', '', as.character(gem_nl$GM_CODE)))
vierkant <- read_sf(paste0(reg_data_dir, 'CBS_vierkant/2018/'), 'CBS_VK100_2018_v1')
vierkant <- st_transform(vierkant, epsg[code == 7415, prj4])  

## road distance
dist_road <- rbindlist(readRDS('../data/road_dist_l2018.rds'))

##
## municipality corrections
##

## read in data on code changes 
code_changes <- data.table(read_excel('../data/Recode Municipalities - Vector 8_11_2018.xlsm', sheet = 3))

## recode migration data
mig_dt_l2018 <- rbindlist(lapply(seq(2012,2018),
                                 function(x) {
                                     if (x != 2018) {              
                                         recode_matrix_recursive(
                                         dat = mig_dt[year == x,
                                                      .(code_orig, code_dest, year,
                                                        TussenGemeentenVerhuisdePersonen_1)
                                                      ],
                                         var_name = "TussenGemeentenVerhuisdePersonen_1",
                                         in_year = x,
                                         out_year = 2018)
                                     } else {
                                         mig_dt[year == x,
                                                .(code_orig, code_dest, year,
                                                        TussenGemeentenVerhuisdePersonen_1)
                                                      ]
                                     }
                                 })
                          )

## recode population data
pop_dt_l2018 <- rbindlist(lapply(seq(2012, 2018),
                                 function(x) {
                                     if (x == 2018) {
                                         bev_dt[year == 2018 & !is.na(BevolkingOp1Januari_1),
                                                .(code, year,
                                                  BevolkingOp1Januari_1,
                                                  BevolkingOp31December_25)
                                                ]
                                     } else {
                                         recode_vector_vars(
                                             dat = bev_dt[year == x],
                                             vars = c("BevolkingOp1Januari_1",
                                                      "BevolkingOp31December_25"),
                                             in_year = x,
                                             out_year = 2018) 
                                     }
                       }))

## recoded building data
wn_vars <- c("BeginstandVoorraad_1", "EindstandVoorraad_8",
             "Nieuwbouw_2", "Sloop_4"
             #"Beginstand_lead", "Eindstand_lead"
             )

wn_dt_l2018 <- rbindlist(lapply(seq(2012, 2018),
                                 function(x) {
                                     if (x == 2018) {
                                         wn_dt[year == 2018 & !is.na(BeginstandVoorraad_1),
                                               c('code', 'year', wn_vars),
                                               with = F
                                               ]
                                     } else {
                                         recode_vector_vars(
                                             dat = wn_dt[year == x],
                                             vars = wn_vars,
                                             in_year = x,
                                             out_year = 2018) 
                                     }
                       }))
    
## housing variable: 
wn_dt_l2018 <- wn_dt_l2018[,
                           .(code,
                             year,
                             supply_growth = ((BeginstandVoorraad_1 + Nieuwbouw_2 - Sloop_4)/BeginstandVoorraad_1)
                             )
                           ]


wn_dt_l2018[, .(mean(supply_growth), sd(supply_growth)), by = year]



##
## create population-weighted centroids
##

## find centroids in each square and overlap with municipality
vk_cents <- data.table(st_join(st_centroid(vierkant[,1:2]),
                               subset(gem_nl, WATER == 'NEE')
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

cent_gem_int <- data.table(mun_id = as.numeric(st_intersects(pop_w_cents_sf, gem_nl)))[, row.id := 1:.N]

missing_cents <- pop_w_cents[which(!cent_gem_int$mun_id %in% seq(1, nrow(pop_w_cents))), GM_CODE]

## allocate randomly those that are missing
if (length(missing_cents) > 0 ) {
  set.seed(123)
  missing_cents <- cbind(missing_cents,
                         data.frame(spsample(as(gem_nl[as.numeric(gsub('GM', '', gem_nl$GM_CODE)) %in% missing_cents, ],
                                                'Spatial'),
                                             1, 'random')
                                    ))
  pop_w_cents <- rbindlist(list(pop_w_cents[!GM_CODE %in% missing_cents[,1]], missing_cents))
} 

## municipality list
id <- cbs_dt[grepl('Gebieden', Title) & grepl('2018', Title), Identifier]
gem_list <- data.table(cbs_get_data(id))[, code := as.numeric(gsub('GM', '', RegioS))]

##
## create variables and training data
##

## calculate distance between municipalities
cent_grid <- data.table(expand.grid(pop_w_cents$GM_CODE, pop_w_cents$GM_CODE))
setkey(cent_grid, Var1)
setkey(pop_w_cents, GM_CODE)
cent_grid[pop_w_cents, ':='(o_x = x, o_y = y)]
setkey(cent_grid, Var2)
cent_grid[pop_w_cents, ':='(d_x = x, d_y = y)]
cent_grid[,
          dist := pointDistance(data.frame(o_x, o_y),
                                data.frame(d_x, d_y), lonlat = F)/1000
          ][,
            dist_test := sqrt((d_x - o_x)^2 + (d_y - o_y)^2) / 1000
            ]


## merge 
setkey(cent_grid, Var1, Var2)
setkey(dist_road, o, d)
cent_grid[dist_road, dist_rd := i.dist/1000]
setkey(cent_grid, Var2, Var1)
cent_grid[dist_road, dist_rd := i.dist/1000]
cent_grid[Var1 == Var2, dist_rd := 0]

## create data table with needed variables: choose average 2014 to 2015
inp_dt <- merge(mig_dt_l2018[year %in% c(2012:2017),
                                      sum(TussenGemeentenVerhuisdePersonen_1),
                                      by = c("code_orig", "code_dest", "year")
                                      ##by = c("code_orig", "code_dest", "year", "age_group")
                                      ][, ## comment if not average over the years
                                        ##mean(V1),
                                        ##by = c("code_orig", "code_dest")
                                        ##by = c("code_orig", "code_dest", "age_group")
                                        ],
                cent_grid[,.(Var1, Var2, dist, dist_rd)],
                by.x = c('code_orig', 'code_dest'),
                by.y = c('Var1', 'Var2'),
                all.y = T)[is.na(V1),
                           V1 := 0
                           ][is.na(dist_rd),
                             dist_rd := 0
                             ]
    
## add inhabitants at destination
inp_dt <- merge(inp_dt[!is.na(year)],
                pop_dt_l2018[year %in% c(2012:2017),
                             mean(BevolkingOp1Januari_1),
                             by = c('code', 'year')],
                by.x = c('code_dest', 'year'), by.y = c('code', 'year'),
                all.x = T)
                
## add delta woningbouw voorrad/(voorraad - (nieuwbouw - sloop)) at destination
inp_dt <- merge(inp_dt,
                wn_dt_l2018[year %in% c(2012:2017),
                            .(supply_growth = mean(supply_growth)
                              ),
                            by = c('code', 'year')],
                by.x = c('code_dest', 'year'), by.y = c('code', 'year'),
                all.x = T
                )

## centrality measure
inp_dt[dist > 0, ':=' (cent = V1.y/dist, cent_rd = V1.y/dist_rd)]

if(!file.exists('../output/centrality_l2018.rds')) {
    ## calculate centrality measure
    no_cores <- detectCores() - 1
    ## Initiate cluster
    ## cl <- makeCluster(no_cores, type="FORK")
    ## registerDoParallel(cl)
    setkey(inp_dt, code_orig)
    ## calculate centrality
    centrality <- list()
    system.time(
        for (yr in unique(inp_dt$year)) {
            print(paste("Calculating centrality for year", yr))
            centrality[[yr]] <- lapply(unique(inp_dt$code_orig),
                                       function(x) {
                                           d_codes <- unique(inp_dt[code_orig == x & code_dest != x & year == yr, code_dest])
                                           tmp1 <- inp_dt[year == yr, ]
                                           tmp2 <- lapply(d_codes,
                                                         function(y) {
                                                             tmp1[code_orig == y  & !code_dest %in% c(x,y),
                                                                    .(x,y,
                                                                      sum(cent, na.rm = T),
                                                                      sum(cent_rd, na.rm = T)
                                                                      )
                                                                    ]
                                                         })
                                           return(rbindlist(tmp2)[, year := yr])
                               }
                               ##mc.cores = no_cores
                               )
        }
    )
    centrality <- rbindlist(lapply(centrality, rbindlist))
    setnames(centrality, 3:4, c('centrality', 'centrality_rd'))
    saveRDS(centrality, file = '../output/centrality_l2018.rds')
} else {
    centrality <- readRDS('../output/centrality_l2018.rds')
}

setkey(inp_dt, code_orig, code_dest, year)
setkey(centrality, x, y, year)
inp_dt[centrality, ':=' (cent = i.centrality, cent_rd = i.centrality_rd)]

# add province name 
setkey(inp_dt, code_orig)
setkey(gem_list, code)
inp_dt[gem_list,
       ':=' (prov_code = trimws(Code_34), prov_nm = trimws(Naam_35))
       ]

## create independent variables
model_vars <- c('V1', 'AANT_INW', 'dist_rd', 'nb', 'cent', 'cent_rd'
#                'hp_diff', 'hp_index_dest', 'p_koop', 'OPP_LAND'
                )

setnames(inp_dt, c(4, 7, 8), c('V1', 'AANT_INW', 'nb'))

## alterations to dependent variable: integerise or replace zero flows with 0.1
inp_dt[code_orig != code_dest,
       V1_int := ifelse(any(V1 %% 1 != 0),
                        .(int_trs(V1)),
                        .(V1)),
       by = c('code_orig', 'year')
       ][code_orig != code_dest & V1 == 0,
         V1 := 0.1
         ]

## calculate ind and dep var: log(x) - mean(log(x))
inp_dt[code_orig != code_dest & dist <= 35,
       paste0('V_', model_vars) := lapply(.SD,
                                          function(x) log(x) - mean(log(x))
                                          ),
       .SDcols = model_vars,
       by = c('code_orig', 'year')
       ]

## for geographixally weighted regression: save spatial point df
sp_dt <- merge(as.data.frame(inp_dt[code_orig != code_dest]),
               pop_w_cents,
               by.x = 'code_orig',
               by.y = 'GM_CODE',
               )

coordinates(sp_dt)  <- ~x+y
sp_dt@proj4string <- CRS(epsg[code == 7415, prj4])

save(sp_dt, file = '../output/sp_dt_l2018.Rdata')

##
## create test data
##

## create data table with needed variables: choose average 2014 to 2015
test_dt <- merge(mig_dt_l2018[year == 2018,
                                       sum(TussenGemeentenVerhuisdePersonen_1),
                                       c("code_orig", "code_dest", "year")
                                      ][, 
                                        mean(V1),
                                        c("code_orig", "code_dest")
                                        ],
                cent_grid[,.(Var1, Var2, dist, dist_rd)],
                by.x = c('code_orig', 'code_dest'),
                by.y = c('Var1', 'Var2'),
                all.y = T
                )[is.na(V1),
                  V1 := 0
                  ][is.na(dist_rd),
                    dist_rd := 0
                    ]
    
## add inhabitants at destination
test_dt <- merge(test_dt,
                pop_dt_l2018[year == 2018, .(code, BevolkingOp1Januari_1)],
                by.x = 'code_dest', by.y = 'code',
                all.x = T)
                
## add delta woningbouw voorrad/(voorraad - (nieuwbouw - sloop)) at destination
test_dt <- merge(test_dt,
                 wn_dt_l2018[year == 2018,
                             .(code, supply_growth)                            
                            ],
                by.x = 'code_dest', by.y = 'code',
                all.x = T)

## centrality measure
test_dt[dist > 0,
        ':='
        (cent = BevolkingOp1Januari_1/dist,
            cent_rd = BevolkingOp1Januari_1/dist_rd
        )]

if(!file.exists('../output/test_centrality.rds')) {
    ## calculate centrality measure
    no_cores <- detectCores() - 1
    ## Initiate cluster
    ## cl <- makeCluster(no_cores, type="FORK")
    ## registerDoParallel(cl)
    setkey(test_dt, code_orig)
    ## calculate centrality
    system.time(
        centrality <- lapply(test_dt[,unique(code_orig)],
                             function(x) {
                                 d_codes <- test_dt[code_orig == x & code_dest != x, unique(code_dest)]
                                 tmp <- lapply(d_codes,
                                                         function(y) {
                                                             test_dt[code_orig == y  & !code_dest %in% c(x,y),
                                                                    .(x,y,
                                                                      sum(cent, na.rm = T),
                                                                      sum(cent_rd, na.rm = T)
                                                                      )
                                                                    ]
                                                         })
                                 return(rbindlist(tmp))                  
                               }
                               #mc.cores = no_cores
                               )
    )
    centrality <- rbindlist(centrality)
    setnames(centrality, 3:4, c('centrality', 'centrality_rd'))
    saveRDS(centrality, file = '../output/test_centrality_l2018.rds')
} else {
    centrality <- readRDS('../output/test_centrality.rds')
}

setkey(test_dt, code_orig, code_dest)
setkey(centrality, x, y)
test_dt[centrality, ':=' (cent = i.centrality, cent_rd = i.centrality_rd)]

# add province name 
setkey(test_dt, code_orig)
setkey(gem_list, code)
test_dt[gem_list,
       ':=' (prov_code = trimws(Code_34), prov_nm = trimws(Naam_35))
       ]

## replace 0 with something small
test_dt[,
        V1_int := V1
        ][V1 == 0,
          V1 := 0.1
          ]

## create independent variables
model_vars <- c('V1', 'AANT_INW', 'dist_rd', 'nb', 'cent_rd'
                )

setnames(test_dt, c(3, 6, 7), c('V1', 'AANT_INW', 'nb'))

## calculate ind and dep var: log(x) - mean(log(x))
test_dt[code_orig != code_dest,
       paste0('V_', model_vars) := lapply(.SD,
                                          function(x) log(x) - mean(log(x))
                                          ),
       .SDcols = model_vars,
       by = code_orig
       ]

saveRDS(test_dt, file = '../output/test_dt_l2018.rds')

##
## checks and tests
##

## check whether matrix recoding works
check_recoding <- lapply(seq(2011, 2017),
                         function(x) {
                             id <- cbs_dt[grepl('Gebieden', Title) & grepl('2018', Title),
                                          Identifier
                                          ]
                             gem_list <- data.table(cbs_get_data(id))
                             gem_list[,
                                      code := as.numeric(gsub('GM', '', as.character(RegioS)))
                                      ]
                             test <- recode_matrix_recursive(
                                 dat = mig_dt[year == x,
                                              .(code_orig, code_dest, year,
                                                TussenGemeentenVerhuisdePersonen_1)
                                              ],
                                 var_name = "TussenGemeentenVerhuisdePersonen_1",
                                 in_year = x,
                                 out_year = 2018)
                             tmp1 <- unique(test$code_orig)[!unique(test$code_orig) %in% gem_list$code]
                             tmp2 <- gem_list$code[!gem_list$code %in% unique(test$code_orig)]
                             tmp3 <- unique(test$code_dest)[!unique(test$code_dest) %in% gem_list$code]
                             tmp4 <- gem_list$code[!gem_list$code %in% unique(test$code_dest)]
                             return(list(tmp1, tmp2, tmp3, tmp4))         
                         })

## check whether matrix recoding works
check_sums1 <- rbindlist(lapply(seq(2011, 2017),
                     function(x) {
                             test <- recode_matrix_year_on_year(
                                 dat = mig_dt[year == x,
                                              .(code_orig, code_dest, year,
                                                TussenGemeentenVerhuisdePersonen_1)
                                              ],
                                 var_name = "TussenGemeentenVerhuisdePersonen_1",
                                 in_year = x,
                                 out_year = x+1)
                             tmp1 <- mig_dt[year == x, sum(TussenGemeentenVerhuisdePersonen_1)]
                             tmp2 <- test[, sum(TussenGemeentenVerhuisdePersonen_1)]
                             return(data.table(actual = tmp1, recoded = tmp2))
                     }))

## check whether matrix works
check_sums2 <- lapply(seq(2011, 2017),
                     function(x) {
                             test <- recode_matrix_recursive(
                                 dat = mig_dt[year == x,
                                              .(code_orig, code_dest, year,
                                                TussenGemeentenVerhuisdePersonen_1)
                                              ],
                                 var_name = "TussenGemeentenVerhuisdePersonen_1",
                                 in_year = x,
                                 out_year = 2018)
                             tmp1 <- mig_dt[year == x, sum(TussenGemeentenVerhuisdePersonen_1)]
                             tmp2 <- test[, sum(TussenGemeentenVerhuisdePersonen_1)]
                             return(data.table(actual = tmp1, recoded = tmp2))
                     })

## check whether vector recoding works
check_sums3 <- rbindlist(lapply(seq(2012, 2017),
                     function(x) {
                             test <- recode_vector_year_on_year(
                                 dat = bev_dt[year == x,
                                              .(code, year, BevolkingOp1Januari_1)
                                              ],
                                 var_name = "BevolkingOp1Januari_1",
                                 in_year = x,
                                 out_year = x+1)
                             tmp1 <- bev_dt[year == x, sum(BevolkingOp1Januari_1, na.rm = T)]
                             tmp2 <- test[, sum(BevolkingOp1Januari_1)]
                             return(data.table(actual = tmp1, recoded = tmp2))
                     }))

## check whether vector recoding works
check_sums4 <- lapply(seq(2012, 2017),
                     function(x) {
                         test <- recode_vector_recursive(
                             dat = bev_dt[year == x & !is.na(BevolkingOp1Januari_1),
                                              .(code, year, BevolkingOp1Januari_1)
                                              ],
                                 var_name = "BevolkingOp1Januari_1",
                                 in_year = x,
                             out_year = 2018)
                         tmp1 <- bev_dt[year == x, sum(BevolkingOp1Januari_1, na.rm = T)]
                         tmp2 <- test[, sum(BevolkingOp1Januari_1)]
                         return(data.table(actual = tmp1, recoded = tmp2))
                     })


## plot against previous migration data
ggplot(merge(reg_bila_dt_2016_code[, sum(value), by = c('code_orig', 'code_dest', 'year')],
      mig_dt_l2018,
      by = c('code_orig', 'code_dest', 'year'),
      all.x=T,all.y=T
      )[code_orig == 363, ]
      ) +
    geom_density(aes(x=V1), col = 'blue', alpha = 0.2) +
    geom_density(aes(x=TussenGemeentenVerhuisdePersonen_1), col = 'red', alpha = 0.2) +
    theme_bw() +
    facet_wrap(~year)

## end script
