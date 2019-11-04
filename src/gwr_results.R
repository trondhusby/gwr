## house keeping
library(sp)
library(data.table)
library(ggplot2)
library(MASS)

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

## run global models
glm_global <- glm(fmla(model_vars, 'mle'), 'poisson', sp_dt@data) 
glm.nb_global <- glm.nb(fmla(model_vars[c(1,2,3)], 'mle'), sp_dt@data) 

glm.nb_global2 <- glm.nb('V1_int ~ log(AANT_INW) + log(dist_rd) + factor(code_orig) + factor(code_dest)', sp_dt@data)

inp_dt <- data.frame(sp_dt@data)
correct_A <- aggregate(V1_int ~ code_orig,  data = sp_dt@data, sum)
init_A <- data.frame(code_orig = unique(inp_dt$code_orig),
                     A = runif(length(unique(inp_dt$code_orig)))
                     )
A <- list()

inp_dt <- merge(inp_dt, aggregate(AANT_INW ~ code_orig, data = inp_dt, sum),
                    by = 'code_orig')

inp_dt <- merge(inp_dt, aggregate(dist_rd ~ code_orig, data = inp_dt, sum),
                    by = 'code_orig')

for (i in 1:3) {
    print(i)
    if (i == 1) {
        tmp <- glm.nb(V1_int ~ log(AANT_INW.x/AANT_INW.y) + log(dist_rd.x/dist_rd.y) + log(A),
                   ##family = 'poisson',
                   data = merge(inp_dt, init_A, by = 'code_orig'))
        A[[i]] <- aggregate(A ~ code_orig,
                       data = data.frame(A = fitted(tmp), code_orig = inp_dt$code_orig),
                       sum)
        inp_dt$y.adj <- fitted(tmp)
    } else {
        tmp <- glm.nb(y.adj ~ log(AANT_INW.x/AANT_INW.y) + log(dist_rd.x/dist_rd.y) + log(A),
                   ##family = 'poisson',
                   merge(inp_dt, A[[i-1]], by = 'code_orig', all.x = T))
        A[[i]] <- aggregate(A ~ code_orig,
                       data = data.frame(A = fitted(tmp), code_orig = inp_dt$code_orig),
                       sum)
        inp_dt$y.adj <- fitted(tmp)
    }
}

out <- rbindlist(lapply(1:3, function(x) {
    data.frame(x, merge(A[[x]], correct_A, by = 'code_orig'))
}))



poisson_pred <- exp(predict(glm_global, newdata = test_dt))
nbinomial_pred <- exp(predict(glm.nb_global, newdata = test_dt))
rmse(poisson_pred, test_dt$V1_int)
rmse(nbinomial_pred, test_dt$V1_int)

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
