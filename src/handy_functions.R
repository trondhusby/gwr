## fitness measure
rmse <- function(m, o) {
  sqrt(mean((m - o)^2))
}

mse <- function(m, o) {
    mean((m - o)^2)
}

mae <- function(m, o) {
    mean(abs(m - o))
}

ploss <- function(m, o) {
    mean(o - m*log(o))
}

## helper function for calculation of flows
gravity_vars <- function(vars, mod = 'ols') {
    if (mod == 'ols') {
        parms <- paste('i.V', vars, sep = '_')
        out <-paste(vars, parms, sep = '^', collapse = '*') 
    } else if (mod == 'mle') {
        parms <- paste('i.V', vars, sep = '_')
        out <- paste(vars, parms, sep = '^', collapse = '*')
    } else if (mod == 'ols_reg') {
        parms <- paste('i.V', vars, sep = '_')
        vars <- paste('V', vars, sep = '_')
        out <- paste(parms, vars, sep = '*', collapse = '+')
    }
    return(out)
}

## check global collinearity
BKWcn <- function(X) {
    p <- dim(X)[2]
    Xscale <- sweep(X, 2, sqrt(colSums(X^2)), "/")
    Xsvd <- svd(Xscale)$d
    Xsvd[1] / Xsvd[p]
}

## generic formula function linear models
fmla_lm <- function(dep_var, ind_var) {
  as.formula(paste(dep_var,'~' , paste(ind_var, collapse= "+")))
}

## generic formula regression formula ols and glm
fmla <- function(vars, mod) {
    if (mod == 'ols') {
        ind_vars <- paste0('V_', vars[-1])
        dep_var <- paste0('V_', vars[1])
        formula <- paste(paste(ind_vars, collapse = ' + '), '- 1')
    } else {
        ind_vars <- paste0('log(', vars[-1], ')')
        dep_var <- paste0(vars[1], '_int')
        ##formula <- paste(ind_vars, collapse = ' + ')
        formula <- paste(paste(ind_vars, collapse = ' + '), ' + factor(code_orig) - 1')
    }
    return(paste(dep_var, formula, sep = ' ~ '))
}

## plot centrality measure
plot_centrality <- function(mun) {
    print(
        ggplot(merge(subset(gem_nl[[2016]], WATER == 'NEE'), centrality[x == mun], by.x = 'code', by.y = 'y')) +
        geom_sf(aes(fill = centrality), size = 0) +
        geom_sf(data= subset(gem_nl[[2016]], code == mun), fill = 'red', size = 0) +
        theme_void() +
        scale_fill_viridis()
    )
}

run_prov_regression <- function(prov) {
    ## create train/test sample split
    trainIndex <- createDataPartition(
        inp_dt[code_orig != code_dest & dist <= 35 & prov_nm == prov, V1],
        p = .75,
        list = FALSE,
        times = 1)
    ## split sample
    dt_train <-inp_dt[code_orig != code_dest & dist <= 35 & prov_nm == prov][trainIndex , ]
    dt_test <- inp_dt[code_orig != code_dest & dist <= 35 & prov_nm == prov][-trainIndex, ]
    ## run regression
    m_fit <- train(as.formula(fmla(model_vars, 'ols')),
               data = dt_train,
               method = 'lm',
               trControl = fitControl
               )
    ## calculate test statistics
    var_importance <- varImp(m_fit, scale = T)
    dt_test$pred <- predict(m_fit, dt_test)
    srmse <- (1/sum(dt_test$V_V1))*rmse(dt_test$pred, dt_test$V1) 
    return(list(prov, var_importance, srmse))
}

## gravity regression function
gravity_reg <- function(dt, ind_vars, dep_var) {
    formula <- paste(paste(ind_vars, collapse = ' + '), '- 1')
    dep_var <- paste(dep_var, '~ ')
    return(dt[,
       summary(lm(paste0(dep_var, formula), .SD))$coefficients,
       by = prov_nm][,
                   ':='(parameter = rep(ind_vars, 4),
                        variable = rep(c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)'), each =length(ind_vars))
                   ),
                   by = prov_nm
                   ][, 
                     ':=' (est_up = V1[variable == 'Estimate'] + 1.96*V1[variable == 'Std. Error'], 
                           est_lo = V1[variable == 'Estimate'] - 1.96*V1[variable == 'Std. Error']
                     ),
                     by = c('prov_nm', 'parameter')
                     ])    
}


gravity_glm_reg <- function(dt, ind_vars, dep_var) {
    ind_vars <- paste0('log(', ind_vars, ')')
    formula <- paste(paste(ind_vars, collapse = ' + '), '- 1')
    dep_var <- paste(dep_var, '~ factor(code_orig) + ')
    out <- lapply(unique(dt$prov_nm),
                  function(x) {
                      df <- glm(paste0(dep_var, formula),
                                        data = dt[prov_nm == x, ],
                                        family = poisson(link = 'log')
                                        )
                      return(list(data.table(prov_nm = x,                                    
                                             actual = as.numeric(dt[prov_nm == x, V1]),
                                             predicted = fitted(df)),
                                  melt(data.table(prov_nm = x,
                                                  parameter = rownames(summary(df)$coefficients),
                                                  summary(df)$coefficients
                                                  )
                                      ,
                                       id.vars = c('prov_nm', 'parameter')
                                       )
                                  )
                             )
                  }
                  )
    return(list(rbindlist(lapply(out, function(x) x[[1]])),
                rbindlist(lapply(out, function(x) x[[2]]))[,
                                                           ':=' (est_up = value[variable == 'Estimate'] + 1.96*value[variable == 'Std. Error'],
                                                                 est_lo = value[variable == 'Estimate'] - 1.96*value[variable == 'Std. Error']
                                                                 ),
                                                           by = c('prov_nm', 'parameter')
                                                           ]
                )
           )
    }

## integerisation
int_trs <- function(x){
  # For generalisation purpose, x becomes a vector
  xv <- as.vector(x) # allows trs to work on matrices
  xint <- floor(xv) # integer part of the weight
  r <- xv - xint # decimal part of the weight
  def <- round(sum(r)) # the deficit population
  # the weights be 'topped up' (+ 1 applied)
  topup <- sample(length(x), size = def, prob = r)
  xint[topup] <- xint[topup] + 1
  dim(xint) <- dim(x)
  dimnames(xint) <- dimnames(x)
  xint
}
