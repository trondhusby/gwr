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
library(viridis)
library(gridExtra)
library(devtools)
library(sf)
library(GWmodel)


setwd('/home/trond/gwr')

## functions
source('src/handy_functions.R')

## read 2006 - 2008 file with data
#source('short_distance_estimation_08.R')

## read data on house prices and costs
cbs_dt <- data.table(cbs_get_toc(Language="nl"))
##cbs_dt[grepl('Voorraad', Title), .(Identifier, Title, ShortTitle, Modified)]

bev_dt <- data.table(cbs_get_data('37259ned',
                                  Geslacht="T001038",
                                  Perioden=c("2014JJ00","2015JJ00", "2016JJ00")
                                  )
                     )[grepl('GM', RegioS),
                       .(RegioS, Perioden,
                         BevolkingOp1Januari_1 = as.numeric(BevolkingOp1Januari_1),
                         code = as.numeric(gsub('GM', '', as.character(RegioS))),
                         year = as.numeric(substr(Perioden, 1, 4))
                         )
                       ]

## building stock
wn_dt <- data.table(cbs_get_data('81955NED',
                                 Gebruiksfunctie='A045364',
                                 Perioden=c(paste0(2013:2016, 'JJ00')))
                    )[grepl('GM', RegioS),
                       .(RegioS, Perioden,
                         BeginstandVoorraad_1 = as.numeric(BeginstandVoorraad_1),
                         Nieuwbouw_2 = as.numeric(Nieuwbouw_2),
                         Sloop_4 = as.numeric(Sloop_4),
                         code = as.numeric(gsub('GM', '', as.character(RegioS))),
                         year = as.numeric(substr(Perioden, 1, 4))
                         )
                       ]

## migration data
load('data/reg_bila_dt_2016_code.Rdata')


## municipality list
gem_list_2016 <- data.table(cbs_get_data('83287NED'))[, code := as.numeric(gsub('GM', '', as.character(RegioS)))]

gem_list_2015 <- data.table(cbs_get_data('82949NED'))[, code := as.numeric(gsub('GM', '', as.character(RegioS)))]
  
## spatial data frame
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

## municipality corrections
recode_vector <- data.table(read_excel('data/Recode Municipalities - Vector.xlsm',
                                       sheet = 3))


## recode population data
pop_dt_l2016 <- rbindlist(list(bev_dt[year == 2016 & !is.na(BevolkingOp1Januari_1),
                                      .(code, year, V2 = BevolkingOp1Januari_1)
                                      ],
                               recode_mun(dat = bev_dt, var_name = "BevolkingOp1Januari_1", in_year = 2015, out_year = 2016),
                               recode_mun(dat = bev_dt, var_name = "BevolkingOp1Januari_1", in_year = 2014, out_year = 2016)
                               )
                          )

## recoded building data
wn_dt_l2016 <- rbindlist(list(wn_dt[year == 2016 & !is.na(BeginstandVoorraad_1),
                                    .(code, year, BeginstandVoorraad_1, Nieuwbouw_2, Sloop_4)
                                    ],
                              recode_mun_vars(dat = wn_dt,
                                              vars = c("BeginstandVoorraad_1", "Nieuwbouw_2", "Sloop_4"), in_year = 2015, out_year = 2016),
                              recode_mun_vars(dat = wn_dt,
                                              vars = c("BeginstandVoorraad_1", "Nieuwbouw_2", "Sloop_4"), in_year = 2014, out_year = 2016),
                              recode_mun_vars(dat = wn_dt,
                                              vars = c("BeginstandVoorraad_1", "Nieuwbouw_2", "Sloop_4"), in_year = 2013, out_year = 2016)
                         ))

## create housing index
wn_dt_l2016[,
            c(paste(paste0('V', 2014:2016), sep = ',')) := lapply(2014:2016,
                                                     function(x)
                                                         BeginstandVoorraad_1[year==x]/BeginstandVoorraad_1[year==x-1]
                                                     ),
            by=code
            ]


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
cent_grid[,
          dist := pointDistance(data.frame(o_x, o_y),
                                data.frame(d_x, d_y), lonlat = F)/1000
          ][,
            dist_test := sqrt((d_x - o_x)^2 + (d_y - o_y)^2) / 1000
            ]


## distance over the road
if(!file.exists('data/road_dist.RData')) {
    source('calculate_distance_osrm.R')
} else {
    load('data/road_dist.RData')
}

## merge 
setkey(cent_grid, Var1, Var2)
setkey(results, o, d)
cent_grid[results, dist_rd := i.dist/1000]
setkey(cent_grid, Var2, Var1)
cent_grid[results, dist_rd := i.dist/1000]
cent_grid[Var1 == Var2, dist_rd := 0]

## create training data

## create data table with needed variables: choose average 2014 to 2015
inp_dt <- merge(reg_bila_dt_2016_code[year %in% c(2014, 2015),
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
                             ]
    
## add inhabitants at destination
inp_dt <- merge(inp_dt,
                pop_dt_l2016[year < 2016, mean(V2), code],
                by.x = 'code_dest', by.y = 'code',
                all.x = T)
                
## add delta woningbouw voorrad/(voorraad - (nieuwbouw - sloop)) at destination
inp_dt <- merge(inp_dt,
                unique(wn_dt_l2016[, .(code, (V2014 + V2015)/2)]),
                by.x = 'code_dest', by.y = 'code',
                all.x = T)

## centrality measure
inp_dt[dist > 0, ':=' (cent = V1.y/dist, cent_rd = V1.y/dist_rd)]

if(!file.exists('output/centrality.rds')) {
    ## calculate centrality measure
    no_cores <- detectCores() - 1
    ## Initiate cluster
    ## cl <- makeCluster(no_cores, type="FORK")
    ## registerDoParallel(cl)
    setkey(inp_dt, code_orig)
    ## calculate centrality
    system.time(
        ## centrality <- foreach(x=inp_dt[,unique(code_orig)],
        ##                      .packages = 'data.table') %dopar% {
        centrality <- lapply(inp_dt[,unique(code_orig)],
                               function(x) {
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
                               #mc.cores = no_cores
                               )
    )
    centrality <- rbindlist(centrality)
    setnames(centrality, 3:4, c('centrality', 'centrality_rd'))
    saveRDS(centrality, file = 'output/centrality.rds')
} else {
    centrality <- readRDS('output/centrality.rds')
}

setkey(inp_dt, code_orig, code_dest)
setkey(centrality, x, y)
inp_dt[centrality, ':=' (cent = i.centrality, cent_rd = i.centrality_rd)]

# add province name 
setkey(inp_dt, code_orig)
setkey(gem_list_2016, code)
inp_dt[gem_list_2016,
       ':=' (prov_code = trimws(Code_34), prov_nm = trimws(Naam_35))
       ]

## create independent variables
model_vars <- c('V1', 'AANT_INW', 'dist_rd', 'nb', 'cent', 'cent_rd'
#                'hp_diff', 'hp_index_dest', 'p_koop', 'OPP_LAND'
                )

setnames(inp_dt, c(3, 6, 7), c('V1', 'AANT_INW', 'nb'))

## alterations to dependent variable: integerise or replace zero flows with 1
inp_dt[code_orig != code_dest,
       V1_int := int_trs(V1),
       by = code_orig
       ][code_orig != code_dest & V1 == 0,
         V1 := 0.1
         ]

## calculate ind and dep var: log(x) - mean(log(x))
inp_dt[code_orig != code_dest,
       paste0('V_', model_vars) := lapply(.SD,
                                          function(x) log(x) - mean(log(x))
                                          ),
       .SDcols = model_vars,
       by = code_orig
       ]

## geographixally weighted regression: save spatial point df
sp_dt <- merge(as.data.frame(inp_dt[code_orig != code_dest]),
               pop_w_cents,
               by.x = 'code_orig',
               by.y = 'GM_CODE',
               )

coordinates(sp_dt)  <- ~x+y
sp_dt@proj4string <- CRS(epsg[code == 7415, prj4])

save(sp_dt, file = 'output/sp_dt.Rdata')

## create test data

## create data table with needed variables: choose average 2014 to 2015
test_dt <- merge(reg_bila_dt_2016_code[year == 2016,
                                      sum(value),
                                      by = c(names(reg_bila_dt_2016_code)[c(1:2, 4)])
                                      ][, 
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
                             ]
    
## add inhabitants at destination
test_dt <- merge(test_dt,
                pop_dt_l2016[year == 2016, V2, code],
                by.x = 'code_dest', by.y = 'code',
                all.x = T)
                
## add delta woningbouw voorrad/(voorraad - (nieuwbouw - sloop)) at destination
test_dt <- merge(test_dt,
                unique(wn_dt_l2016[, .(code, V2016)]),
                by.x = 'code_dest', by.y = 'code',
                all.x = T)

## centrality measure
test_dt[dist > 0, ':=' (cent = V2/dist, cent_rd = V2/dist_rd)]

if(!file.exists('output/test_centrality.rds')) {
    ## calculate centrality measure
    no_cores <- detectCores() - 1
    ## Initiate cluster
    ## cl <- makeCluster(no_cores, type="FORK")
    ## registerDoParallel(cl)
    setkey(test_dt, code_orig)
    ## calculate centrality
    system.time(
        ## centrality <- foreach(x=inp_dt[,unique(code_orig)],
        ##                      .packages = 'data.table') %dopar% {
        centrality <- lapply(test_dt[,unique(code_orig)],
                               function(x) {
                                   dt <- rbindlist(lapply(test_dt[code_orig == x & code_dest != x , unique(code_dest)],
                                                         function(y)
                                                             test_dt[code_orig == y  & !code_dest %in% c(x,y),
                                                                    .(x,y,
                                                                      sum(cent, na.rm = T),
                                                                      sum(cent_rd, na.rm = T)
                                                                      )
                                                                    ]
                                                         )
                                                  )
                               }
                               #mc.cores = no_cores
                               )
    )
    centrality <- rbindlist(centrality)
    setnames(centrality, 3:4, c('centrality', 'centrality_rd'))
    saveRDS(centrality, file = 'output/test_centrality.rds')
} else {
    centrality <- readRDS('output/test_centrality.rds')
}

setkey(test_dt, code_orig, code_dest)
setkey(centrality, x, y)
test_dt[centrality, ':=' (cent = i.centrality, cent_rd = i.centrality_rd)]

# add province name 
setkey(test_dt, code_orig)
setkey(gem_list_2016, code)
test_dt[gem_list_2016,
       ':=' (prov_code = trimws(Code_34), prov_nm = trimws(Naam_35))
       ]

## replace 0 with something small
test_dt[,
        V1_int := V1
        ][V1 == 0,
          V1 := 0.1
          ]

## create independent variables
model_vars <- c('V1', 'AANT_INW', 'dist_rd', 'nb', 'cent', 'cent_rd'
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

saveRDS(test_dt, file = 'output/test_dt.rds')

## end script
