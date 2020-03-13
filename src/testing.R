library(sp)
library(data.table)
library(ggplot2)
library(caret)
library(MASS)
library(numDeriv)

rand_dt <- data.frame(a = runif(100), b = runif(100, 0, 100))
apply(rand_dt, 2, function(x) var(scale(x, center = F)))

var(rand_nr)
var(scale(rand_nr, center = F))


setwd('~/gwr')

## read handy functions
source('src/handy_functions.R')

## load training data
load('output/sp_dt.Rdata')
sp_dt <- subset(sp_dt, dist <= 35)


# Generate some artificial training and test data
x <- 1:100
y <- 5 + .1*x + rnorm(100)
xtrain <- sample(x, size=80)
ytrain <- y[xtrain]
xtest <- x[-xtrain]
ytest <- y[-xtrain]

# Compute forecasts from a linear model
forecast <- predict(lm(ytrain~xtrain), newdata=data.frame(xtrain=xtest))

# Plot training data, test data and forecasts
plot(xtrain, ytrain)
lines(xtest,forecast,col='red',pch=19)
points(xtest,ytest,col='blue',pch=19)

# Compute accuracy statistics
accuracy(forecast,ytest)

# with mase
forecast <- forecast(lm(ytrain~xtrain), newdata=data.frame(xtrain=xtest))
accuracy(forecast,ytest)

forecast <- structure(list(mean=rep(10,20), fitted=rep(10,80),
   x=ytrain), class='forecast')

accuracy(forecast,ytest)



## subset
sp_dt <- subset(sp_dt, prov_nm %in% unique(sp_dt@data$prov_nm)[c(7:9)])
inp_dt <- data.frame(sp_dt@data)

## load test data
test_dt <- readRDS('output/test_dt.rds')[dist <= 35 & code_orig != code_dest]
## subset
test_dt <- subset(test_dt, prov_nm %in% unique(test_dt$prov_nm)[c(7:9)])

## dependent (1) and independent (-1) variables
model_vars <- c('V1', 'AANT_INW', 'dist_rd', 'nb', 'cent_rd'
#                'hp_diff', 'hp_index_dest', 'p_koop', 'OPP_LAND'
                )

## simplified version of gwr to test

bal_obj_fn <- function(adj_i, y_i, x_i, betas_i, loss_fn = 'l2') {
    if (loss_fn == 'l2') {
        sum((y_i - exp(colSums(adj_i * betas_i * t(x_i))))^2)
    } else {
        sum(abs(y_i - exp(colSums(adj_i * betas_i * t(x_i)))))
    }
}

bal_obj_fn2 <- function(y_i, x_i, betas_i, loss_fn = 'l2') {
    if (loss_fn == 'l2') {
        sum((y_i - exp(colSums(betas_i * t(x_i))))^2)
    } else {
        sum(abs(y_i - exp(colSums(betas_i * t(x_i)))))
    }
}

bal_obj_fn2 <- function(adj_i, y_i, x_i, betas_i, loss_fn = 'l2') {
    if (loss_fn == 'l2') {
        sum((y_i - exp(betas_i[1] + colSums(adj_i * betas_i[-1] * t(x_i[,-1]))))^2)
    } else {
        sum(abs(y_i - exp(betas_i[1] + colSums(adj_i * betas_i[-1] * t(x_i[,-1])))))
    }
}

mod <- list()
betas <- array(dim = c(length(unique(inp_dt$code_orig)),
                       length(model_vars))
               )
betas_l <- list()
bal_con <- list()
y <- inp_dt$V1_int
mu <-inp_dt$V1_int + 0.1
nu <- log(mu)
theta <- rep(1, length(mu))
alpha  <- 1/theta
code_dt <- data.frame(code = inp_dt$code_orig, idx = 1:nrow(inp_dt))
x_mat <- model.matrix(eval(parse(text = fmla(model_vars, 'nbinomial'))),
                      data = inp_dt)

theta_g <- glm.nb(fmla(model_vars, 'nbinomial'), data = inp_dt)$theta
l_fn <- 'l2'

llik <- 0.0
iter <- 1
repeat {
    ## calculate elements of Fisher Information matrix
    a <- mu / (1 + alpha * mu) +
        (y - mu)*(alpha*mu) / (1 + 2*alpha*mu + alpha**2 * mu**2)
    y.adj <- nu + (y - mu) / (a * (1 + alpha * mu))
    bal_con[[iter]] <- rep(1, nrow(betas))
    for (idx in seq_along(unique(inp_dt$code_orig))) {
        cod <- unique(inp_dt$code_orig)[idx]
        sub_idx <- subset(code_dt, code == cod)$idx        
        if (iter == 1) {
            mod[[idx]] <- lm(y.adj[sub_idx] ~ x_mat[sub_idx,]-1)
            betas[idx, ] <- coefficients(mod[[idx]])
            nu[sub_idx] <- fitted(mod[[idx]])
            bal_con[[iter]][idx] <- ifelse(l_fn == 'l2',
                                           sum((y[sub_idx] - exp(nu[sub_idx]))^2),
                                           sum(abs(y[sub_idx] - exp(nu[sub_idx])))
                                           )
        } else {
            if (any(is.na(betas[idx, ]))) {
                message(paste('NAs in '), cod)
                na_betas <- which(is.na(betas[idx, ]))
                betas[idx, na_betas] <- 1
           }
            out <- optim(rep(1, length(betas[idx, ])),
                         bal_obj_fn,
                         y_i = y[sub_idx],
                         betas_i = betas[idx, ],
                         x_i= x_mat[sub_idx, ],
                         loss_fn = l_fn,
                         hessian = T
                         )
            ## out_grad <- grad(bal_obj_fn,
            ##                  y_i = y[sub_idx],
            ##                  betas_i = betas[idx, ],
            ##                  x_i= x_mat[sub_idx, ],
            ##                  x = out$par
            ##                  )
            ##     betas[idx, ] <- betas[idx, ] + diag(out_grad / out$hessian)
            betas[idx, ] <- betas[idx, ] * out$par
            nu[sub_idx] <- colSums(betas[idx, ] * t(x_mat[sub_idx, ]))           
            bal_con[[iter]][idx] <- out$value
       }
    }
    betas_l[[iter]] <- betas
    mu <- exp(nu)
    old.llik <- llik
    llik <- sum(dnbinom(y,size=theta_g,mu=mu,log=TRUE)) + 0*sum(bal_con[[iter]])
    print(llik)
    if (abs((old.llik - llik)/llik) < 10^-6) break
    iter <- iter+1
    print(iter)
    if (iter >= 200) break
}

rmse(inp_dt$V1, mu)
mae(inp_dt$V1, mu)

glm_global <- glm(fmla(model_vars, 'poisson'), 'poisson', sp_dt@data)
glm.nb_global <- glm.nb(fmla(model_vars, 'nbinomial.fe'), sp_dt@data)
lm_global <- lm(fmla(model_vars, 'ols'), sp_dt@data)
lm_fitted <- calculate_fitness_ols(dat = test_dt,
                                   pars = data.table(code_orig = unique(inp_dt$code_orig),
                                                     t(coefficients(lm_global))
                                                     )
                                   )

betas_out <- data.table(code_orig = unique(inp_dt$code_orig),
                        betas)
setkey(betas_out, code_orig)
setkey(test_dt, code_orig)
test_dt[betas_out, names(betas_out), with=F]

yhat <- betas_out[test_dt,
                  2:6] * test_dt[,
                                 c(model_vars)[-1],
                                 with=F][,
                                         c(1, lapply(.SD, log))
                                         ]

tmp <- data.table(y = test_dt$V1,
                  poisson = exp(predict(glm_global, newdata = test_dt)),
                  nbinom = exp(predict(glm.nb_global, newdata = test_dt)),
                  lm = lm_fitted,
                  my_own = exp(rowSums(yhat))
                  )

ggplot(melt(tmp, id.vars = 'y'), aes(x=y)) +
    geom_point(aes(col = variable, y = value)) +
    scale_colour_brewer(palette = 'Set1') + 
    theme_bw()


tmp[,
    lapply(.SD,
           function(x) {
               rmse(x, y)
           }),
    .SDcols = 2:5
    ]



rmse(inp_dt$V1, fitted(glm_global))
rmse(inp_dt$V1, fitted(glm.nb_global))
rmse(inp_dt$V1, fitted(lm_global))


ggplot(melt(data.table(poisson = resid(glm_global, 'pearson'),
                       nbinom = resid(glm.nb_global, 'pearson'),
                       ols = resid(lm_global, 'pearson')
           )), aes(value)) +
    geom_density() +
    facet_wrap(~variable, scales = 'free_x') +
    theme_bw()

ggplot(residuals(glm_global))


out_cons <- rbindlist(lapply(1:(iter-1),
                        function(x) {
                            data.table(iter = x,
                                       code = unique(inp_dt$code_orig),
                                       cons = bal_con[[x]])
                        }))

ggplot(out_cons[iter > 10], aes(factor(iter), cons)) +
    geom_boxplot() +
    theme_bw()

out_betas <- rbindlist(lapply(1:iter,
                        function(x) {
                            data.table(iter = x,
                                       code = unique(inp_dt$code_orig),
                                       betas_l[[x]])
                        }))

ggplot(melt(out_betas, id.vars = c('iter', 'code')),
       aes(factor(iter), value)) +
    geom_boxplot() +
    facet_wrap(~variable, scales = 'free') +
    theme_bw()
 

data.table(y[sub_idx],
      exp(colSums(betas[idx,] * t(x_mat[sub_idx, ])))
      )[,
          sum((V1 - V2)^2)
        ]

bal_obj_fn(rep(1, 5), y_i = y[sub_idx],
           betas_i = betas[idx, ],
           x_i= x_mat[sub_idx,]
           )
    
    
test <- optim(rep(1, 5),
              bal_obj_fn,
              y_i = y[sub_idx],
              betas_i = betas[idx, ],
              x_i= x_mat[sub_idx,]
              )




test_dt <- model.matrix(~ y + x_i.1 + x_i.2 + factor(a) -1,
                        data = data.frame(a = rep(c(1:3), each = 3),
                                          y = runif(9),
                                          x_i = matrix(runif(18), ncol = 2)
                                          ),
                        )

test_out <- lm(y~.-1, data=data.frame(test_dt))

coef(test_out)
fitted(test_out)
cbind(data.frame(test_dt)$y, fitted(test_out))


optim(rep(1, 4),
      bal_obj_fn,
      y_i = c(1:4),
      #x0 = rep(1, 4),
      x_i = matrix(c(1:8), ncol = 2),
      betas_i = c(1, 2, 3),
      method = "L-BFGS-B"
      )



## run global models
glm_global <- glm(fmla(model_vars, 'poisson'), 'poisson', sp_dt@data) 
glm.nb_global <- glm.nb(fmla(model_vars, 'nbinomial.fe'), sp_dt@data) 
glmeln_lambda <- cv.glmregNB(eval(parse(text = fmla(model_vars, 'nbinomial.fe'))),
                             data = sp_dt@data,
                             parallel = T,
                             n.cores = 4
                             )

glm.nb.el_global <- glmregNB(eval(parse(text = fmla(model_vars, 'nbinomial.fe'))),
                             data = sp_dt@data,
                             lambda = glmeln_lambda$lambda.optim
                           )

cv_5 = trainControl(method = "cv", number = 5)
glm.nb.el2_global <- train(
    eval(parse(text = fmla(model_vars, 'nbinomial.fe'))),
    data = sp_dt@data,
    method = "glmnet",
    trControl = cv_5
)



init_vars <- data.frame(code_orig = inp_dt$code_orig,
                        code_dest = inp_dt$code_dest,
                        A = runif(nrow(inp_dt)),
                        B = runif(nrow(inp_dt))
                        )

## check constraint
sum((aggregate(V1_int - yhat ~ code_orig,
          data = data.frame(inp_dt[, c('code_orig', 'V1_int')],
                            yhat = fitted(glm_global)),
          sum
          )[, 2])^2)


rmse(inp_dt$V1_int, exp(predict(glm_global, newdata = test_dt)))


sum((aggregate(V1_int - yhat ~ code_orig,
          data = data.frame(inp_dt[, c('code_orig', 'V1_int')],
                            yhat = fitted(glm.nb_global)),
          sum
          )[, 2])^2)

rmse(inp_dt$V1_int, exp(predict(glm.nb_global, newdata = test_dt)))


sum((aggregate(V1_int - yhat ~ code_orig,
          data = data.frame(inp_dt[, c('code_orig', 'V1_int')],
                            yhat = fitted(glm.nb.el_global),
          sum
          )[, 2])^2))

rmse(inp_dt$V1_int, exp(predict(glm.nb.el_global,
                                newx=data.frame(code_orig = test_dt$code_orig,
                                                test_dt[,
                                                        c(model_vars)[-1],
                                                        with=F]
                                                )
                                )
                        )
     )

sum((aggregate(V1_int - yhat ~ code_orig,
          data = data.frame(inp_dt[, c('code_orig', 'V1_int')],
                            yhat = fitted(glm.nb.el2_global)),
          sum
          )[, 2])^2)

test_xmat <- model.matrix(~ . +factor(code_orig)-1,
                          data.frame(code_orig = test_dt$code_orig,
                                     test_dt[,
                                             c(model_vars)[-1],
                                             with=F][,
                                                     lapply(.SD, log)
                                                     ]
                                     )
                          )[,-1]



rmse(inp_dt$V1_int, exp(predict(glm.nb.el2_global,
                                newx = test_xmat)
                        )
     )                                       

coefficients(glm.nb_global)[-c(1:4)]

inp_mat <- model.matrix(~ V1_int + log(AANT_INW) + log(dist_rd) +
                            log(nb)+ log(cent_rd) + factor(code_orig)- 1,
                        data = inp_dt
                        )

glm_fe_fn <- function(pars) {
    yhat <- apply(t(t(inp_mat[, -1]) * c(coefficients(glm.nb_global)[c(1:4)],
                              pars *coefficients(glm.nb_global)[-c(1:4)]
                              )
                   ),
                 1,
                 function(x) exp(sum(x))
                 )
    out <- aggregate(V1_int - yhat ~ code_orig,
                     data = data.frame(inp_dt[, c('code_orig', 'V1_int')],
                                       yhat),
                     sum)
    return(sum(out[, 2]^2))
}


glm_fe_fn(rep(1, (length(coefficients(glm.nb_global))-4)))
glm_fe_fn(runif(length(coefficients(glm.nb_global))-4))

glm_

test <- optim(par = rep(1, (length(coefficients(glm.nb_global))-4)),
              fn = glm_fe_fn,
              method = 'L-BFGS-B',
              lower = 0,
              upper = 5
              )

sum(c(6.2, 6.5, (0.5*9), 0.4, 3, 10, 4.6, 21.7, 0.2, (1.5*11), 3, 3.5))

glm_fn <- function(pars) {
    A <- pars[1:(length(pars)/2)]
    B <- pars[(1+length(pars)/2):length(pars)]
    A <- data.frame(code_orig = unique(inp_dt$code_orig), A=A)
    B <- data.frame(code_orig = unique(inp_dt$code_orig), B=B)
    dat <- merge(inp_dt, A, by = 'code_orig', all.x = T)
    dat <- merge(dat, B, by = 'code_orig', all.x = T)
    res <- glm.nb(V1_int ~ log(A) + log(B) + log(AANT_INW) + log(dist_rd)  - 1,
                  data = dat,       
                  )
    yhat <- aggregate(yhat ~ code_orig,
                      data = data.frame(code_orig = inp_dt$code_orig,
                                        yhat = fitted(res)),
                      sum
                      )
    y <- aggregate(V1_int ~ code_orig,  data = inp_dt, sum) 
    return(sum((y - yhat)^2))
}

glm_fn(pars = rep(runif(length(unique(inp_dt$code_orig))), 2))

test <- optim(par = rep(runif(length(unique(inp_dt$code_orig))), 2),
              fn = glm_fn,
              method = 'L-BFGS-B',
              lower = -10,
              upper = 10
              )



A <- data.frame(code_orig = unique(inp_dt$code_orig),
                A = test$par[1:(length(test$par)/2)]
                )
B <- data.frame(code_orig = unique(inp_dt$code_orig),
               B = test$par[(1+length(test$par)/2):length(test$par)]
               )

dat <- merge(inp_dt, A, by = 'code_orig', all.x = T)
dat <- merge(dat, B, by = 'code_orig', all.x = T)
res <- glm.nb(V1_int ~ log(A) + log(B) + log(AANT_INW) + log(dist_rd)  - 1,
              data = dat
              )

yhat_opt <- aggregate(yhat ~ code_orig,
                      data = data.frame(code_orig = inp_dt$code_orig,
                                        yhat = fitted(res)),
                      sum
                      )

yhat <- list()
for (i in 1:250) {
    if (i == 1) {
        res <- glm.nb(V1_int ~ log(A) + log(B) + log(AANT_INW) + log(dist_rd)  - 1,
                      data = merge(inp_dt, init_vars, by = c('code_orig', 'code_dest'))
                      )
    } else {
        dat <- merge(inp_dt, A, by = 'code_orig', all.x = T)
        dat <- merge(dat, B, by = 'code_orig', all.x = T)
        res <- glm.nb(V1_int ~ log(A) + log(B) + log(AANT_INW) + log(dist_rd)  - 1,
                      data = dat
                      )
    }
    A <- aggregate(A ~ code_orig,
                   data = data.frame(code_orig = inp_dt$code_orig, A = fitted(res)),
                   FUN = function(x) 1/sum(x)
                   )
    B <- aggregate(B ~ code_orig,
                   data = data.frame(code_orig = inp_dt$code_orig,
                                     B = apply(res$coefficients[-c(1,2)] * t(res$model[, -c(1, 2, 3)]), 2, sum)),
                   FUN = function(x) 1/sum(x)
                   )
    yhat[[i]] <- cbind(aggregate(. ~ code_orig,
                             data = data.frame(code_orig = inp_dt$code_orig,
                                               yhat = fitted(res),
                                               y = inp_dt$V1_int),
                             sum
                             ),
                         i)
}

res2 <- glm.nb(V1_int ~ log(AANT_INW) + log(dist_rd) + factor(code_orig) - 1,
                      data = dat
                      )

ggplot(rbindlist(yhat)[, .(code_orig, i, abs(yhat - y))],
       aes(i, V3, group = code_orig)) +
    geom_line(aes(col = code_orig)) +
    theme_bw()


correct_A <- aggregate(V1_int ~ code_orig,  data = sp_dt@data, sum)

out <- rbindlist(lapply(1:3, function(x) {
    data.frame(x, merge(A[[x]], correct_A, by = 'code_orig'))
}))



poisson_pred <- exp(predict(glm_global, newdata = test_dt))
nbinomial_pred <- exp(predict(glm.nb_global, newdata = test_dt))
rmse(poisson_pred, test_dt$V1_int)
rmse(nbinomial_pred, test_dt$V1_int)
