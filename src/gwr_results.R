## house keeping
library(sp)
library(data.table)
library(ggplot2)
library(MASS)
library(mpath)
library(caret)

## vector of models
mdl  <-  c('poisson.fe', 'nbinomial.fe', 'nbinomial.fe.global', 'nbinomial.fe.fixed')[1:2]

setwd('~/gwr')

## read handy functions
source('src/handy_functions.R')

## load training data
load('output/sp_dt.Rdata')
sp_dt <- subset(sp_dt, dist <= 35)

## load result data
result <- readRDS('output/gwr_poisson_nbinomial_result.rds')
for (i in 1:length(result)) {
    result[[i]]$code_orig <- sp_dt@data$code_orig
}

## subset
sp_dt <- subset(sp_dt, prov_nm %in% unique(sp_dt@data$prov_nm)[c(1:3)])

## load test data
test_dt <- readRDS('output/test_dt.rds')[dist <= 35 & code_orig != code_dest]
## subset
test_dt <- subset(test_dt, prov_nm %in% unique(test_dt$prov_nm)[c(1:3)])

summary(result[[1]][, "residual"])
summary(result[[2]][, "residual"])

## dependent (1) and independent (-1) variables
model_vars <- c('V1', 'AANT_INW', 'dist_rd', 'nb', 'cent_rd'
#                'hp_diff', 'hp_index_dest', 'p_koop', 'OPP_LAND'
                )

## find unique combinations of ind vars
mods <- list()
counter <- 1
for (i in 1:4) {
    combs <- combn(model_vars[-1], i, simplify = F)
    for (j in 1:length(combs)) {
        mods[[counter]] <- combs[[j]]
        counter <- counter + 1
    }
}

## does constraints hold?
lapply(1:2,
       function(x) {
           data.table(result[[x]])[,
                                   .(sum(yhat), sum(y)),
                                   by = code_orig
                                   ]
       })


## create test data and calculate test errors
pred_errors <- lapply(mdl,
                    function(x) {
                        if (grepl('poisson', x)) {
                            lst_index <- which(mdl == x)
                            id_vars <- paste0('log(', mods[[15]], ')')
                        } else {
                            lst_index <- which(mdl == x)
                            id_vars <- paste0('log(', mods[[5]], ')')
                        }
                        frml <- paste0('~', paste(id_vars, collapse = '+'),
                                       '+factor(code_orig)-1')
                        dt1 <- data.table(model.matrix(as.formula(frml),
                                                       data = test_dt
                                                       ))
                        dt2 <- data.table(result[[lst_index]][, 1:ncol(dt1)])
                        out <- data.table(predicted = exp(apply(dt1 * dt2, 1, sum)),
                                          actual = test_dt$V1_int
                                          )[, ':='(
                                             mse = rmse(predicted, actual),
                                             mae = mean(abs(predicted - actual))
                                         )]
                        out$code_orig <- sp_dt@data$code_orig
                        out$code_dest <- sp_dt@data$code_dest
                        return(out)
                    })

out_rmse <- merge(pred_errors[[1]][, rmse(predicted, actual), by = code_orig],
                  pred_errors[[2]][, rmse(predicted, actual), by = code_orig],
                  by = 'code_orig')

in_rmse <- merge(pred_errors[[1]][, rmse(predicted, actual), by = code_dest],
                  pred_errors[[2]][, rmse(predicted, actual), by = code_dest],
                 by = 'code_dest')

ggplot(in_rmse, aes(V1.x, V1.y)) +
    geom_point(aes(col = ifelse(V1.x > V1.y, 'red', 'blue'))) +
    theme_bw()

ggplot(out_rmse, aes(V1.x, V1.y)) +
    geom_point(aes(col = ifelse(V1.x > V1.y, 'red', 'blue'))) +
    theme_bw()

ggplot(melt(out_rmse, id.vars = 'code_orig'), aes(x=value)) +
    geom_density(aes(col = variable)) +
    theme_bw()

ggplot(melt(in_rmse, id.vars = 'code_dest'), aes(x=value)) +
    geom_density(aes(col = variable)) +
    theme_bw()




## calculate test error
