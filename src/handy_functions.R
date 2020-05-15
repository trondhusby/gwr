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

mape <- function(dat, pred) {
    mean(abs((dat - pred) / dat))
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
    } else if (mod == 'mle_cv') {
        parms <- paste('V', vars, sep = '_')
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
fmla <- function(vars, mod = 'generic', year_dummy = TRUE) {
    if (grepl('ols', mod)) {
        ind_vars <- paste0('V_', vars[-1])
        dep_var <- paste0('V_', vars[1])
        formula <- paste(paste(ind_vars, collapse = ' + '), '- 1')
    } else if (mod %in% c('poisson.fe', 'nbinomial.fe')) {
        ind_vars <- paste0('log(', vars[-1], ')')
        dep_var <- paste0(vars[1], '_int')
        formula <- ifelse(year_dummy,
                          paste(paste(ind_vars, collapse = ' + '), ' + factor(code_orig)*factor(year) - 1'),
                          paste(paste(ind_vars, collapse = ' + '), ' + factor(code_orig) - 1'))                        
    } else if (mod %in% c('poisson', 'nbinomial')) {
        ind_vars <- paste0('log(', vars[-1], ')')
        dep_var <- paste0(vars[1], '_int')
        formula <- paste(paste(ind_vars, collapse = ' + '))
    } else {
        ind_vars <- vars[-1]
        dep_var <- vars[1]
        formula <- paste(paste(ind_vars, collapse = ' + '), ' - 1')
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

## helper function to plug in gravity constraints
gravity_constraints <- function(parm, vars) {
    if (length(parm) != length(vars)) {
        stop('Par and var vector of unequal length')
    } else {
        paste(unlist(lapply(1:length(vars),
                            function(x) {
                                paste0(vars[x],
                                '^',
                                parm[x]
                                )}
                            )),
                     collapse = '*')
    }
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

## recoding municipalities vector year on year
recode_vector_year_on_year <- function(dat, var_name, in_year, out_year) {
    if (out_year - in_year > 1) stop("Not year on year input")
    cat(paste('out year', out_year), '\n')
    ## data
    tmp1 <- dat[!is.na(get(var_name)),
               c('code', 'year', var_name),
               with=F]
    ## subset code changes
    rec_dt <- code_changes[(jaar == in_year & maand > 1) | (jaar == out_year & maand == 1),
                           .(jaar, oud_code, frac, nieuw_code, maand, dag)
                           ]
    ## recode data
    tmp2 <- merge(rec_dt,
                 tmp1,
                 by.x = 'oud_code', by.y = 'code', all.x = T, all.y = F
                 )[,
                   .(oud_code, nieuw_code, year, get(var_name) * frac)
                   ][V4>0,
                     sum(V4),
                     by = c('nieuw_code', 'year')
                     ]
    ## combine with non-modified data  
    old_codes <- rec_dt[!frac == 0 | nieuw_code %in% oud_code, unique(oud_code)]
    ## gather
    out <- unique(rbindlist(list(tmp2,
                                 tmp1[!code %in% old_codes]
                                 )
                            )[,
                              sum(V1),
                              by = c('nieuw_code', 'year')
                              ])
    setnames(out, c('code', 'year', var_name))
    return(out)
}
    
## recoding municipalities vector year on year
recode_vector_recursive <- function(dat, var_name, in_year, out_year) {
    ## initialisation
    in_year_tmp <- in_year
    tmp <- dat
    ## update codes recursively each year
    while(out_year - in_year_tmp > 0) {
        if (out_year - in_year_tmp > 1) {
            tmp <- recode_vector_recursive(tmp,
                                           var_name,
                                           in_year_tmp,
                                           in_year_tmp + 1)
        } else {
            out <- recode_vector_year_on_year(tmp,
                                              var_name,
                                              in_year_tmp,
                                              out_year)
        }
        in_year_tmp <- in_year_tmp + 1
    }
    ## update year
    out[, year := in_year]
    return(out)
}
    

## function for recoding municipalities (vector version)
recode_vec <- function(dat, var_name, in_year, out_year) {
    out <- dat[year==in_year & !is.na(get(var_name)),
               c('code', 'year', var_name),
               with=F]
    for (yr in seq(in_year, out_year-1)) {
        rec_dt <- code_changes[jaar == yr + 1 | (jaar == yr & maand > 1),
                                .(jaar, oud_code, frac, nieuw_code, maand, dag)
                                ]
        tmp <- merge(rec_dt,
                     out,
                     by.x = 'oud_code', by.y = 'code', all.x = T, all.y = F
                     )[,
                       .(oud_code, nieuw_code, get(var_name) * frac)
                       ][V3>0,
                         .(in_year, sum(V3)),
                         by = nieuw_code
                         ]
        out <- unique(rbindlist(list(tmp,
                                     out[!code %in% rec_dt[frac<1, oud_code]]
                                     )
                                )[,
                                  .(in_year, sum(V2)),
                                  by = nieuw_code
                                  ])
        setnames(out, c('code', 'year', var_name))
    }
    return(out)
}

## function for recoding municipalities - several variables
recode_vector_vars <- function(dat,vars, in_year, out_year) {
    tmp <- lapply(seq_along(vars),
                  function(x) {
                      recode_vector_recursive(dat, vars[x], in_year, out_year)
                  })
    for (i in 1:length(tmp)) {
        if (i == 1) {
            out <- tmp[[i]]
        } else {
            out <- merge(out, tmp[[i]], by = c('code', 'year'))
        }
    }
    return(out)
}

## function for recoding municipality codes year on year
recode_matrix_year_on_year <- function(dat, var_name, in_year, out_year) {
    if (out_year - in_year > 1) stop("Not year on year input")
    rec_dt <- code_changes[(jaar == in_year & maand > 1) | (jaar == out_year & maand == 1),
                           .(jaar, oud_code, frac, nieuw_code, maand, dag)
                           ]
    ## replace origins
    tmp1 <- merge(dat[code_orig %in% rec_dt$oud_code | code_dest %in% rec_dt$oud_code],
                 rec_dt,
                 by.x = 'code_orig',
                 by.y='oud_code',
                 all.x = T,
                 allow.cartesian = T)
    ## replace destinations
    tmp2 <- merge(tmp1, rec_dt,
                 by.x = 'code_dest',
                 by.y='oud_code',
                 all.x = T,
                 allow.cartesian = T)
    ## fix na's in origins or destinations
    tmp2[is.na(frac.y),
        ':=' (frac.y = 1, nieuw_code.y = code_dest)
        ][is.na(frac.x),
          ':=' (frac.x = 1, nieuw_code.x = code_orig)
          ][,
            moved := frac.x * frac.y * get(var_name)
            ]
    ## calculate moves over new municipalities
    tmp <- tmp2[frac.x * frac.y > 0,
               sum(moved),
               by = c('nieuw_code.x', 'nieuw_code.y', 'year')
               ]
    ## combine with non-modified data  
    old_codes <- rec_dt[!frac == 0 | nieuw_code %in% oud_code, unique(oud_code)]
    out <- unique(
        rbindlist(list(tmp,
                       dat[!(code_orig %in% old_codes | code_dest %in% old_codes)]
                       )
                  )[,
                    sum(V1),
                    by = c('nieuw_code.x', 'nieuw_code.y', 'year')
                    ])
    setnames(out, c('code_orig', 'code_dest', 'year', var_name))
    ## fix friese meeren...
    if (in_year == 2016) {
        out[code_orig == 1921,
            code_orig := 1940
            ][code_dest == 1921,
              code_dest := 1940
              ]
    }
    return(out)
}

## function for recoding municipality codes recursively
recode_matrix_recursive <- function(dat, var_name, in_year, out_year) {
    ## initialisation
    in_year_tmp <- in_year
    tmp <- dat
    ## update codes recursively each year
    while(out_year - in_year_tmp > 0) {
        if (out_year - in_year_tmp > 1) {
            tmp <- recode_matrix_recursive(tmp,
                                           var_name,
                                           in_year_tmp,
                                           in_year_tmp + 1)
        } else {
            out <- recode_matrix_year_on_year(tmp,
                                              var_name,
                                              in_year_tmp,
                                              out_year)
        }
        in_year_tmp <- in_year_tmp + 1
    }
    ## update year
    out[, year := in_year]
    return(out)
}

## predicting flows using the gravity function
pred_fn <- function(dat, pars) {   
    setkey(dat, code_orig, code_dest)
    setkey(pars, code_orig, code_dest)
    ivars <- unlist(lapply(names(pars), function(x) grep('log', x, value = T)))
    setnames(pars,
             ivars,
             gsub('log\\(', 'V_', gsub('\\)', '', ivars))
             ) 
    out <- dat[pars,
               ][,
                 pred_denom := sum(eval(parse(text = gravity_vars(model_vars[-1]))), na.rm = T),
                 by = code_orig,
                 ][,
                   prediction := sum(V1, na.rm = T)*eval(parse(text = gravity_vars(model_vars[-1])))/pred_denom,
                   by = code_orig
                   ]
    return(out[,.(code_orig, code_dest, actual = V1, prediction)])
}

## list with files from storage
azure_storagefile_list <- function(files, job_name) {
    config <- rjson::fromJSON(file = paste0("credentials.json"))
    ## Get Storage Account Credentials
    storageCredentials <- rAzureBatch::SharedKeyCredentials$new(
                                                                name = config$sharedKey$storageAccount$name,
                                                                key = config$sharedKey$storageAccount$key
                                                            )
    storageAccountName <- storageCredentials$name
    ## Generate Storage Client for SAS Token
    storageClient <- rAzureBatch::StorageServiceClient$new(
                                                           authentication = storageCredentials,
                                                           url = sprintf("https://%s.blob.%s",
                                                                         storageCredentials$name,
                                                                         config$sharedKey$storageAccount$endpointSuffix
                                                                         )
                                                       )
    ## create read only credentials
    readSasToken <- storageClient$generateSasToken(permission = "r", "c",
                                                   path = job_name
                                                   )
    ## create list of storage files
    storageFiles <- lapply(1:nrow(files),
                           function(x) {
                               rAzureBatch::createBlobUrl(storageAccountName,
                                                          job_name,
                                                          files[[1]][x],
                                                          readSasToken)
                           })
    return(storageFiles)
}

## end of code
