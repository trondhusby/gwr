## title: estimation of gravity model via azure batch
## author: trond
## date: nov 14

## house keeping
library(sp)
library(doAzureParallel)
library(rAzureBatch)
library(GWmodel)
library(raster)
library(data.table)
library(MASS)
library(ggplot2)
library(rgenoud)

## login
system('az login')

## Create a resource group
system('az group create --name gwr --location westeurope')
## Create a storage account
system('az storage account create --name gwrstorage --resource-group PEARL --location westeurope --sku Standard_LRS --kind StorageV2')
## Create a batch account
system('az batch account create --name gwr --storage-account gwrstorage --resource-group PEARL --location westeurope')

## generate .json file
generateCredentialsConfig("credentials.json")

## generate json for running cluster
generateClusterConfig("cluster.json")

## show storage and batch account info
system('az storage account show-connection-string --name gwrstorage')
system('az batch account list --resource-group PEARL')
system('az batch account keys list --name gwr --resource-group PEARL')

## after filled in set credentials
setCredentials("credentials.json")

## create your cluster if it does not exist; this takes a few minutes
system.time(
    cluster <- makeCluster("cluster.json")
)

## register your parallel backend 
registerDoAzureParallel(cluster) 

## check that the nodes are running 
getDoParWorkers() 

## load necessary data
load('../output/sp_dt_l2018.Rdata')

## testing on a subset (comment out for full)
sp_dt <- subset(sp_dt, dist <= 35)
##sp_dt <- subset(sp_dt, code_orig == unique(sp_dt@data$code_orig)[c(1,2)])
##sp_dt <- subset(sp_dt, prov_nm %in% unique(sp_dt@data$prov_nm)[c(1)])

## add index and fold
set.seed(123)
sp_dt@data$index <- 1:length(sp_dt)
sp_dt@data$V1_int  <-  as.integer(sp_dt@data$V1_int)
n_folds <- 10
n_reps <- 1
rep_fold_grid <- expand.grid(seq(1,n_reps), seq(1, n_folds))
## create folds: make sure that they are evenly distributed within municipality of orig
for (idx in seq(1,n_reps)) {
    tmp <- unlist(lapply(unique(sp_dt@data$code_orig),
                                   function(x) {
                                       n <- nrow(subset(sp_dt@data, code_orig == x))
                                       return(sample(rep(sample(1:n_folds),
                                                         length.out = n)))
                                   })
                  )
    sp_dt@data[, paste0('fold_', idx)]  <-  tmp
}

## load updated gwr poisson functions
source('gwr.poisson_upd.R')
##source('~/GWmodel/R/GeneralizedGWR.r')

## load road distance data
load('../data/road_dist.RData')
road_distance <- rbind(results, data.frame(o = results$d, d = results$o, dist = results$dist))

## if using subset
road_distance <- subset(road_distance, o %in% sp_dt@data$code_orig & d %in% sp_dt@data$code_orig)

## load useful functions
source('handy_functions.R')

## dependent (1) and independent (-1) variables
model_vars <- c('V1', 'AANT_INW', 'dist_rd', 'nb', 'cent_rd')[-4]

## find unique combinations of ind vars
mods <- list()
counter <- 1
for (i in 1:3) {
    combs <- combn(model_vars[-1], i, simplify = F)
    for (j in 1:length(combs)) {
        mods[[counter]] <- combs[[j]]
        counter <- counter + 1
    }
}

## kernel and bandwidth
##kernel_bw <- data.table(Var1 = 'bisquare', Var2 = seq(300, 1200, by = 10))
kernel_bw <- data.table(Var1 = 'bisquare', Var2 = seq(1050, 1200, by = 1))

## check residuals and obtain parameters with non-gwr model
##glm_global <- glm(fmla(model_vars, 'poisson.fe'), 'poisson', sp_dt@data) 
## glm.nb_global <- glm.nb(fmla(model_vars, 'nbinomial.fe'), sp_dt@data) 
## summary(data.table(cbind(sp_dt@data$V1_int, fitted(glm_global)))[, V1 - V2])
## summary(data.table(cbind(sp_dt@data$V1_int, fitted(glm.nb_global)))[, V1 - V2])
## hat_matrix <- influence(glm_global)$hat

## Azure options create name, prevent R session from being open and waiting for results
my_job_id <- "panel-run-estimate2"
#setAutoDeleteJob(TRUE)
opt <- list(job = my_job_id, wait = FALSE, autoDeleteJob = TRUE,
            enableCloudCombine = FALSE, maxTaskRetryCount = 0,
            chunkSize = 1)
#            chunkSize = floor(nrow(kernel_bw) / 13))

## euclidian or road for weight matrix
euclidian_distance <- TRUE

## 10-fold cross validation or not
cv_run <- FALSE

job_id <- foreach(
    k = seq(1, nrow(kernel_bw))[131],
    .packages = c('sp', 'raster', 'GWmodel', 'data.table', 'rgenoud'),
    .options.azure = opt
) %dopar% {
    ## set model type
    ##i = c('ols', 'mle')
    i = c('ols', 'mle')[2]
    ## set model type for mle
    mdl = c('poisson.fe', 'nbinomial', 'nbinomial.fe', 'nbinomial.fe.fixed')[1]
    ## set kernel type
    #k = seq(1, nrow(kernel_bw))[1]
    j = paste(unlist(kernel_bw[k, 1]))
    ## initialisation
    start_time <- Sys.time()
    ## set seed for reproducibility
    set.seed(123)
    ##message(paste('Loaded packages:', Sys.time()))
    ## calculate distance between regression points
    message(paste('Calculate distance matrix, time:', Sys.time()))
    dist_m <- matrix(0, nrow = length(sp_dt), ncol = length(sp_dt))
    for (m in seq(length(sp_dt))) {
        if (euclidian_distance) {
            dist_m[m,] <- matrix(pointDistance(coordinates(sp_dt)[m,],
                                           coordinates(sp_dt),
                                           longlat = F),
                             ncol = length(sp_dt))
        } else {
            code <- unique(sp_dt@data[m,]$code_orig)
            dist <- merge(data.frame(V1 = sp_dt@data$code_orig),
                          road_distance[road_distance$o == code, c('d', 'dist')],
                          by.x = 'V1',
                          by.y = 'd',
                          all.x = T
                          )
            dist_m[m, !is.na(dist$dist)] <- matrix(dist$dist[!is.na(dist$dist)])
        }
    }
    ## determine best bandwidth
    if (i == 'ols') {
        message(paste('Calculate BW OLS, time:', Sys.time()))
        best_bw <- NULL
        ## best_bw <- tryCatch({
        ##     bw.gwr(fmla(model_vars, mdl),
        ##            data = sp_dt, adaptive = T, dMat = dist_m,
        ##            approach = 'CV', kernel = j)
        ## }
       ## ,
       ##  error = function(e) message(paste0('Error in BW calculation:', e))
       ##  )
    } else {
        message(paste('Calculate BW MLE, time:', Sys.time()))
        best_bw <- unlist(kernel_bw[k,2])
        ## best_bw <- tryCatch({
        ##     bw.ggwr2(fmla(model_vars, mdl),
        ##             data = sp_dt, adaptive = T, dMat = dist_m,
        ##             approach = 'cv', kernel = j, family = mdl,
        ##             f.bw = unlist(kernel_bw[k, 2]))
        ## },
        ## error = function(e) message(paste0('Error in BW calculation:', e))
        ## )
    }
    ## if na result return optimal bw from gaussian
    if(!is.null(best_bw)) {
        #message(paste0('Best Bw found :', best_bw))
    } else {
        best_bw <- c(105, 557, 557, 301, 400)[4]
        ##best_bw <- unlist(kernel_bw[k,2])
    }
    ## run gwr model
    if (i == 'ols') {
        message(paste('Run GWR OLS, time:', Sys.time()))
        mod <- tryCatch({
            gwr.basic(fmla(model_vars, 'ols'),
                      data=sp_dt, bw = 563, dMat = dist_m,
                      kernel = 'bisquare', adaptive = T)
        },
        error = function(e) message(paste0('Error in GWR OLS estimation:', e))
        )
    } else {
        ##message(paste('Run GWR MLE, time:', Sys.time()))
        mod <- tryCatch({
            if (cv_run) {                
                message(paste('Estimating model', fmla(model_vars, mdl), '\n Bandwith', best_bw, '\n Kernel', j))
                ## initialisation
                best_bw <- unlist(kernel_bw[k,2])
                out <- list()
                vars <- list()
                params <- list()
                preds <- list()
                for (idx in seq(1, 5)) {
                    yr <- seq(2012,2016)[idx]
                    ## set bandwidth relative to 2018
                    bw_yr <- (yr-2012+1)/(2018-2012+1)
                    ## training data set: training data everything except this fold
                    fld <- rep_fold_grid[idx, 'Var2']
                    rep_var <- paste0('fold_', rep_fold_grid[idx, 'Var1'])
                    message(paste('Now estimating fold',
                                   idx,                                  
                                  Sys.time(),
                                   '\n'))
                    ## train <- subset(sp_dt, eval(parse(text = rep_var)) != fld)
                    ## test <- subset(sp_dt, eval(parse(text = rep_var)) == fld)
                    train <- subset(sp_dt, year <= yr)
                    test <- subset(sp_dt, year == yr + 1)
                    model_formula <- ifelse(idx == 1,
                                            fmla(model_vars, mdl, year_dummy = FALSE),
                                            fmla(model_vars, mdl)
                                            )
                    print(model_formula)
                    ## estimate gwr model on training data
                    ##source('gwr.poisson_upd.R')
                    out[[idx]] <- bw.ggwr2(model_formula,
                                       data = train, f.bw = round(bw_yr*best_bw),
                                       dMat = dist_m[train@data$index, train@data$index],
                                       kernel = j, adaptive = T,
                                       family = mdl, approach = 'kfold.cv'
                                       )
                    message(paste('Finished estimating model:', Sys.time()))
                    ## gather parameters for each municipality 
                    params[[idx]] <- unique(cbind(code_orig = train@data$code_orig,
                                                  out[[idx]][,-c(seq(ncol(out[[idx]]),
                                                                     ncol(out[[idx]])-2))],
                                                  idx)
                                            )
                    ## selection of output
                    out[[idx]] <- cbind(code_orig = train@data$code_orig,
                                      out[[idx]][, c('y', 'yhat', 'residuals')],
                                      idx
                                      )
                    ## ## calculate prediction errors
                    ivars <- unlist(lapply(colnames(params[[idx]]),
                                           function(x) grep('log', x, value = T))
                                    )
                    colnames(params[[idx]])[colnames(params[[idx]]) %in% ivars] <- gsub('log\\(', 'V_', gsub('\\)', '', ivars))
                    preds[[idx]] <- setDT(merge(test@data[, c(1:10)],
                                                params[[idx]][, c(1:5)],
                                                by = 'code_orig'
                                                ))[,
                                                   pred_denom := sum(eval(parse(text = gravity_vars(model_vars[-1], 'mle_cv'))), na.rm = T),
                                                   by = code_orig,
                                                   ][,
                                                     prediction := sum(V1, na.rm = T)*eval(parse(text = gravity_vars(model_vars[-1], 'mle_cv')))/pred_denom,
                                                     by = code_orig
                                                     ][,
                                                       ':=' (
                                                           fold = idx,
                                                           bw = best_bw,
                                                           kernel = j,
                                                           oos_rmse = rmse(V1, prediction),
                                                           oos_mae = mae(V1, prediction),
                                                           ws_rmse = sqrt(mean(out[[idx]][, 'residuals']^2)),
                                                           ws_mae = mean(abs(out[[idx]][, 'residuals']))
                                                       )]
                }
                pred_errors <- unique(rbindlist(preds)[,
                                                       .(fold,
                                                         rep_var,
                                                         bw,
                                                         kernel,
                                                         oos_rmse,
                                                         oos_mae,
                                                         ws_rmse,
                                                         ws_mae
                                                         )])
                message(paste('Finished CV run:', Sys.time()))
            } else {
                message(paste('Estimating model', fmla(model_vars, mdl)))
                if (mdl == 'poisson.fe') {
                    print(best_bw)
                    ## out <- ggwr.basic2(fmla(model_vars, mdl),
                    ##    data=sp_dt, bw = best_bw, dMat = dist_m, maxiter = 50, 
                    ##    kernel = j, adaptive = T, family = mdl, no.hatmatrix = FALSE,
                    ##    cv = F
                    ##    )
                    out <- ggwr.basic2(fmla(model_vars, 'poisson'),
                       data=sp_dt, bw = best_bw, dMat = dist_m, maxiter = 50, 
                       kernel = j, adaptive = T, family = 'poisson', no.hatmatrix = FALSE,
                       cv = F
                       )
                } else {
                    ggwr.basic2(fmla(model_vars, mdl),
                       data=sp_dt, bw = 340, dMat = dist_m, maxiter = 50,
                       kernel = j, adaptive = T, family = mdl, no.hatmatrix = FALSE,
                       cv = F
                       )
                }
            }
        },
        error = function(e) message(paste0('Error in GWR MLE estimation:', e))
        )
    }
    #message(paste('End:', Sys.time()))
    if (!cv_run) {
        return(list(mod$SDF@data, j, best_bw))
    } else if (cv_run) {
        ##return(list(out, params, k, best_bw, kernel_bw[k, ], preds))
        return(pred_errors)
    }
}

rbindlist(lapply(job_id,
       function(x) {
           x[,
             lapply(.SD, function(x) return(c(mean(x), sd(x)))),
             .SDcols = c('oos_rmse', 'oos_mae', 'ws_rmse', 'ws_mae'),
             by = c('bw', 'kernel')
             ][,
               variable := c('mean', 'sd')
               ]
       }))


## does constraint hold?
lapply(1:length(job_id),
       function(x) {
           data.table(sp_dt@data$code_orig,
                      sp_dt@data$year,
                      ##sp_dt@data$age_group,
                      job_id[[x]][[1]]$y,
                      job_id[[x]][[1]]$yhat
                      )[,
                        sum(V4) - sum(V3),
                        by = c('V1', 'V2')
                        ]
       })


## quick out-of-sample check

## gather results
result <- getJobResult(my_job_id)

ggplot(rbindlist(result)[,
                         .(mean(oos_rmse),
                           mean(oos_rmse) + sd(oos_rmse) / sqrt(5),
                           mean(oos_rmse) - sd(oos_rmse) / sqrt(5)
                           ),
                         by = bw
                         ],
       aes(bw, V1)) +
    geom_ribbon(aes(ymax = V2, ymin = V3), fill = 'grey') + 
    geom_line() +
    theme_bw()


## saveRDS(result, file = '../output/gwr_poisson_nbinomial_result.rds')

ggplot(rbindlist(result)[,mean(oos_rmse),by=bw], aes(bw,V1)) +
    geom_line() +
    theme_bw()

rbindlist(result)[, mean(oos_rmse), by = bw]


pred_errors <- rbindlist(lapply(job_id,
       function(x) {
           x[,
             lapply(.SD, function(x) return(c(mean(x), sd(x)))),
             .SDcols = c('oos_rmse', 'oos_mae', 'ws_rmse', 'ws_mae'),
             by = c('bw', 'kernel')
             ][,
               variable := c('mean', 'sd')
               ]
       }))

ggplot(pred_errors[variable == 'mean'], aes(bw, oos_rmse)) +
    geom_line() +
    theme_bw()


## save to rds
saveRDS(result_osfit, file = '../output/gwr_results_sd_os_bw.rds')

## if download does not work...

# List all of the blobs that start with result in container
files <- listStorageFiles("gwr-poisson-tuning", prefix = "result")

## download files
temp_dir <- tempdir(check=T)

result <- lapply(1:length(files[1:100,1]),
       function(x) {
           getStorageFile("gwr-nbionom-40",
                          files[x,1],
                          downloadPath = paste(temp_dir, x, sep = '/')
                          )
           return(readRDS(paste(temp_dir, x, sep = '/')))
       })

saveRDS(result, file = '../output/gwr_results_sd_os_bw2.rds')
unlink(temp_dir,force=T)

## shut down your pool
stopCluster(cluster)

## reading from local
result_bw <- lapply(list.files('~/Downloads', pattern = '.rds')[c(2,5:7)],
                    function(x) readRDS(paste0('~/Downloads/', x))
                    )

## checking estimated parameters
out <- sapply(1:length(result_bw),
                       function(x)
                           subset(result_bw[[x]][[1]][[6]][[1]], code_orig == 518)
                       )
out <- data.frame(t(cbind(names(result_bw[[1]][[1]][[6]][[1]]), out))[-1,])

## shut down idle nodes
system('az batch account login --name gwr --resource-group PEARL --shared-key-auth')

system('IDLE_NODES=$(az batch node list --pool-id gwr --query [*].id --output tsv --filter "state eq 'idle'")')

system('az batch node delete --pool-id gwr --node-list $IDLE_NODES')

system('az batch pool autoscale disable --pool-id gwr')

## end script
