## title: tests with gwrnbinom
## date: 
## author: trond

## house keeping
library(GWmodel)
library(rgdal)
library(data.table)
library(cbsshape)
library(sf)
library(ggplot2)
library(MASS)

## download geographical data
epsg <- data.table(make_EPSG())
reg_data_dir <- '/home/trond/Documents/CLAiR-City_data/RegionalData/'
gem_nl <- list()
gem_nl[[2016]] <-read_sf(paste0(reg_data_dir, 'buurt_en_wijk_kaarten/2016/Uitvoer_shape/'), 'gem_2016')
gem_nl[[2016]] <- st_transform(gem_nl[[2016]],
                               epsg[code == 7415, prj4]) # set RDS planar projection
gem_nl[[2016]]$code <- as.numeric(gsub('GM', '', gem_nl[[2016]]$GM_CODE))
gem_nl[[2016]] <- subset(gem_nl[[2016]], WATER == 'NEE')
centroids <- data.frame(st_coordinates(st_centroid(gem_nl[[2016]])))
st_geometry(gem_nl[[2016]]) <- NULL

load(file = '../output/sp_dt.Rdata')

## global models
m1 <- glm(OAD ~ OPP_LAND + AANT_INW,
          family = 'poisson',
          data = gem_nl[[2016]])
m2 <- glm.nb(OAD ~ OPP_LAND + AANT_INW,
            data = gem_nl[[2016]])

out <- data.frame(model = c('Poisson', 'Negative Binomial'),
                     aic = c(AIC(m1), AIC(m2)),
                     bic = c(BIC(m1), BIC(m2)),
           theta = c('', m2$theta)
           )

st_geometry(gem_nl[[2016]]) <- centroids
    
sp_dt <- lapply(


    gem_nl[[2016]][, c('OAD', 'OPP_LAND', 'AANT_INW')],
                function(x) 

    as(gem_nl[[2016]], 'Spatial')
sp_dt


coordinates(sp_dt)  <- ~x+y
sp_dt@proj4string <- CRS(epsg[code == 7415, prj4])


    



                     
                             




