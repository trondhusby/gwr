## title: calculate distance between two points with the osrm api
## author: trond
## date: 03.12.2018

## build osrm backend
## https://github.com/Project-OSRM/osrm-backend/wiki/Building-OSRM

## set up backend
## https://hub.docker.com/r/osrm/osrm-backend/
## https://www.tutorialspoint.com/docker/docker_images.htm

## house keeping
library(raster)
library(data.table)
library(doParallel)
library(httr)
library(jsonlite)
library(sf)
library(rgdal)
library(maptools)
library(cbsodataR)
library(rgeos)

## directory for the spatial data
reg_data_dir <- '/home/trond/Documents/CLAiR-City_data/RegionalData/'

## read spatial data
epsg <- data.table(make_EPSG())
## RDS 28992 or 74515
## wgs84 4326 
gem_nl <- read_sf(paste0(reg_data_dir, 'buurt_en_wijk_kaarten/2018/'), 'gemeente_2018_v2')
gem_nl <- st_transform(gem_nl, epsg[code == 4326, prj4]) # set RDS planar projection
## clean up mun code
gem_nl$code <- as.numeric(gsub('GM', '', as.character(gem_nl$GM_CODE)))
vierkant <- read_sf(paste0(reg_data_dir, 'CBS_vierkant/2018/'), 'CBS_VK100_2018_v1')
vierkant <- st_transform(vierkant, epsg[code == 4326, prj4])

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
                           crs = st_crs(gem_nl),
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
cbs_dt <- data.table(cbs_get_toc(Language="nl"))
id <- cbs_dt[grepl('Gebieden', Title) & grepl('2018', Title), Identifier]
gem_list <- data.table(cbs_get_data(id))[, code := as.numeric(gsub('GM', '', RegioS))]

## calculate distance between municipalities
cent_grid <- data.table(expand.grid(pop_w_cents$GM_CODE, pop_w_cents$GM_CODE))
setkey(cent_grid, Var1)
setkey(pop_w_cents, GM_CODE)
cent_grid[pop_w_cents, ':='(o_x = x, o_y = y)]
setkey(cent_grid, Var2)
cent_grid[pop_w_cents, ':='(d_x = x, d_y = y)]
cent_grid[,
          dist := pointDistance(data.frame(o_x, o_y),
                                data.frame(d_x, d_y), lonlat = T)/1000
          ]

## needed variables
## url  <- "http://router.project-osrm.org" #demo server
url_backend <- "http://127.0.0.1:5000" ## local backend
##url  <- "http://127.0.0.1:9966" ## frontend

pckgs <- c('data.table', 'httr', 'jsonlite')

## calculate the number of cores
no_cores <- detectCores() - 2

## initiate cluster
cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)

## api query for each od pair (minus the last region)
system.time(
    results <- foreach(i = unique(cent_grid$Var1)[-nrow(gem_list)],
                       .packages=pckgs) %do%
        {
            ## find destinations not previously in the od list
            j_list <- unique(cent_grid$Var1)[unique(cent_grid$Var1) > i]
            dist_dt <- list()
            ## loop over each destination
            tryCatch({
                for (j in j_list) {                    
                    ## create url with od pairs
                    path <- paste0("/route/v1/driving/",
                                   cent_grid[Var1 == i & Var2 == j, o_x],
                                   "_",
                                   cent_grid[Var1 == i & Var2 == j, o_y],
                                   "-",
                                   cent_grid[Var1 == i & Var2 == j, d_x],
                                   "_",
                                   cent_grid[Var1 == i & Var2 == j, d_y],
                                   "?overview=false")
                    path <- gsub('_', ',', gsub('-', ';', path))
                    ## initialisation
                    attempt <- 1
                    raw.result <- NULL
                    ## api query: try 10 times
                    while( is.null(raw.result) && attempt <= 10) {
                        attempt <- attempt + 1
                        if (attempt == 5) {message("Something funky going on...")}
                        try(
                            raw.result <- fromJSON(rawToChar(GET(url = url_backend,
                                                                 path = path)$content))$routes$distance
                        )
                        ## sleep for some times if the request bounces back
                        if (is.null(raw.result)) {Sys.sleep(2)}
                    }
                    ## gather result into a data table
                    dist_dt[[j]] <- data.table(o = i, d =j, dist = ifelse(is.null(raw.result), NA, raw.result))
                }
                return(rbindlist(dist_dt))
            },
            error = function(e) return(paste0("The variable '", i, "'",
                                              " caused the error: '", e, "'"))
            )
        }
)

stopImplicitCluster()
registerDoSEQ()

saveRDS(results, file = '../data/road_dist_l2018.rds')
