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
load('../output/sp_dt.Rdata')

## testing on a subset (comment out for full)
sp_dt <- subset(sp_dt, dist <= 35)
##sp_dt <- subset(sp_dt, code_orig == unique(sp_dt@data$code_orig)[c(1,2)])
sp_dt <- subset(sp_dt, prov_nm %in% unique(sp_dt@data$prov_nm)[c(1:3)])

## add index and fold
set.seed(123)
sp_dt@data$index <- 1:length(sp_dt)
sp_dt@data$V1_int  <-  as.integer(sp_dt@data$V1_int)
n_folds <- 10
## create folds: make sure that they are evenly distributed within municipality of orig
sp_dt@data$fold_i <- unlist(lapply(unique(sp_dt@data$code_orig),
                                   function(x) {
                                       n <- nrow(subset(sp_dt@data, code_orig == x))
                                       return(sample(rep(sample(1:n_folds), length.out = n)))
                                   })
                           ) 

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
model_vars <- c('V1', 'AANT_INW', 'dist_rd', 'nb', 'cent_rd')

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

## combination of band width and model type
## mods_bw <- expand.grid(seq(140, 200, by = 1),
##                         c(5, 11, 12, 15)
##                         )

## for bisquare
## mods_bw <- expand.grid(seq(300, 360, by = 1),
##                          c(5, 11, 12, 15)
##                        )

## remove from grid
removal_conditions <- "(Var1 %in% c('bisquare', 'tricube') & Var2 < 200)"
removal_conditions <- paste(removal_conditions,
                            "(Var1 %in% c('boxcar') & Var2 < 65)",
                            sep = "|")
removal_conditions <- paste(removal_conditions,
                            "(Var1 %in% c('gaussian') & Var2 > 200)",
                            sep = "|")
removal_conditions <- paste(removal_conditions,
                            "(Var1 %in% 'boxcar' & Var2 > 300)",
                            sep = "|")

## grid with kernel type and size
kernel_bw <- data.table(expand.grid(c('gaussian', 'bisquare', 'tricube', 'boxcar')[4],
                                    seq(50, 200, by = 1))
                        )[!eval(parse(text=removal_conditions)),
                          ]


## check residuals and obtain parameters with non-gwr model
## glm_global <- glm(fmla(model_vars, 'poisson.fe'), 'poisson', sp_dt@data) 
## glm.nb_global <- glm.nb(fmla(model_vars, 'nbinomial.fe'), sp_dt@data) 
## summary(data.table(cbind(sp_dt@data$V1_int, fitted(glm_global)))[, V1 - V2])
## summary(data.table(cbind(sp_dt@data$V1_int, fitted(glm.nb_global)))[, V1 - V2])

## Azure options create name, prevent R session from being open and waiting for results
my_job_id <- "eureka-seminar-19"
#setAutoDeleteJob(TRUE)
opt <- list(job = my_job_id, wait = FALSE, autoDeleteJob = TRUE, enableCloudCombine = FALSE, maxTaskRetryCount = 1)

## euclidian or road for weight matrix
euclidian_distance <- TRUE

## 10-fold cross validation or not
cv_run <- TRUE

job_id <- foreach(
    i = c('ols', 'mle')[2]
     ## .packages = c('sp', 'raster', 'GWmodel', 'data.table'),
     ## .options.azure = opt
) %do% {
    ## set model type
    ##i = c('ols', 'mle')
    ## set model type for mle
    mdl = c('poisson.fe', 'nbinomial', 'nbinomial.fe', 'nbinomial.fe.fixed')[1]
    ## set kernel type
    k = seq(1, nrow(kernel_bw))[1]
    j = paste(unlist(kernel_bw[k, 1]))
    ## initialisation
    start_time <- Sys.time()
    ## set seed for reproducibility
    set.seed(123)
    message(paste('Loaded packages:', Sys.time()))
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
        #best_bw <- tryCatch({
        #    bw.ggwr2(fmla(model_vars, i),
        #            data = sp_dt, adaptive = T, dMat = dist_m,
        #            approach = 'CV', kernel = j, family ="poisson")
        #},
        #error = function(e) message(paste0('Error in BW calculation:', e))
        #)
    }
    ## if na result return optimal bw from gaussian
    if(!is.null(best_bw)) {
        message(paste0('Best Bw found :', best_bw))
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
        message(paste('Run GWR MLE, time:', Sys.time()))
        mod <- tryCatch({
            if (cv_run) {                
                message(paste('Estimating model', fmla(model_vars, mdl), '\n Bandwith', best_bw, '\n Kernel', j))
                ## initialisation
                best_bw <- unlist(kernel_bw[k,2])
                out <- list()
                vars <- list()
                params <- list()
                preds <- list()
                for (fld in 1:n_folds) {
                    cat(paste0('Now estimating fold ', fld, '\n'))
                    ## training data set: training data everything except this fold
                    train <- subset(sp_dt, fold_i != fld)
                    test <- subset(sp_dt, fold_i == fld)
                    ## estimate gwr model on training data
                    out[[fld]] <- ggwr.basic2(fmla(model_vars, mdl),
                                       data = train, bw = best_bw,
                                       dMat = dist_m[train@data$index, train@data$index],
                                       kernel = j, adaptive = T,
                                       cv = F, family = mdl, no.hatmatrix = TRUE,
                                       theta_g = glm.nb_global$theta, null.dev = glm.nb_global$null.deviance
                                       )$SDF@data                   
                    ## gather parameters for each municipality 
                    params[[fld]] <- unique(cbind(code_orig = train@data$code_orig,
                                                  out[[fld]][,-c(seq(ncol(out[[fld]]),
                                                                     ncol(out[[fld]])-2))],
                                                  fld)
                                            )                   
                    ## selection of output
                    out[[fld]] <- cbind(code_orig = train@data$code_orig,
                                      out[[fld]][, c('y', 'yhat', 'residual')],
                                      fld
                                      )
                    ## ## calculate prediction errors
                    ## ivars <- unlist(lapply(names(params[[fld]]),
                    ##                        function(x) grep('log', x, value = T))
                    ##                 )
                    ## names(params[[fld]])[names(params[[fld]]) %in% ivars] <- gsub('log\\(', 'V_', gsub('\\)', '', ivars))
                    ## preds[[fld]] <- setDT(merge(test@data[, c(1:9)],
                    ##                             params[[fld]][, c(1:5)],
                    ##                             by = 'code_orig'
                    ##                             ))[,
                    ##                                pred_denom := sum(eval(parse(text = gravity_vars(model_vars[-1], 'mle_cv'))), na.rm = T),
                    ##                                by = code_orig,
                    ##                                ][,
                    ##                                  prediction := sum(V1, na.rm = T)*eval(parse(text = gravity_vars(model_vars[-1], 'mle_cv')))/pred_denom,
                    ##                                  by = code_orig
                    ##                                  ][,
                    ##                                    ':=' (
                    ##                                        fold = fld,
                    ##                                        rmse_f = rmse(V1, prediction),
                    ##                                        mae_f = mae(V1, prediction)
                    ##                                    )]
                }
                ## pred_errors <- rbindlist(preds)[, .(unique(rmse_f), unique(mae_f)), by = fold]
            } else {
                #message(paste('Estimating model', fmla(model_vars_sub, i)))
                if (mdl == 'poisson.fe') {
                    ggwr.basic2(fmla(model_vars, mdl),
                       data=sp_dt, bw = 136, dMat = dist_m, maxiter = 50, 
                       kernel = 'boxcar', adaptive = T, family = mdl, no.hatmatrix = FALSE,
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
    message(paste('End:', Sys.time()))
    if (!cv_run) {
        return(mod$SDF@data)        
    } else if (cv_run){
        return(list(out, params, k, best_bw, kernel_bw[k, ], preds))
        ##return(list(k, kernel_bw[k], pred_errors))
   }
}

## does constraint hold?
lapply(1:length(job_id),
       function(x) {
           data.table(sp_dt@data$code_orig,
                      ##sp_dt@data$age_group,
                      job_id[[x]]$y,
                      job_id[[x]]$yhat
                      )[,
                        sum(V3) - sum(V2),
                        by = V1
                        ]
       })


## quick out-of-sample check

## gather results
result <- getJobResult(job_id)

## saveRDS(result, file = '../output/gwr_poisson_nbinomial_result.rds')

## postprocessing of results from tuning
##job_name <- c("gwr-poisson-tuning", "gwr-nbinomial-tuning")[2]
job_name <- my_job_id
files <- subset(listStorageFiles(job_name, prefix = "results"),
                nchar(as.character(ContentLength)) > 4)

storageFiles <- azure_storagefile_list(files, job_name)

## calculate mae and rmse
my_job_id <- "calculate-errors-eureka-13"
#setAutoDeleteJob(TRUE)
opt <- list(job = my_job_id, wait = FALSE, autoDeleteJob = TRUE,
            enableCloudCombine = FALSE, maxTaskRetryCount = 0
            )

## run postprocessing job: NB remove large stuff from session
result <- foreach(
    ## dummy index
    idx = 1,   
    .packages = c('sp', 'raster', 'GWmodel', 'data.table'),
    .options.azure = opt
) %dopar% {
    ## download files
    temp_dir <- tempdir(check=T)
    ## download files from storage to temporary directory    
    for (i in 1:length(storageFiles)) {
        download.file(storageFiles[[i]],
                      destfile = paste0(temp_dir, gsub("results/", "/", files[[1]][i]))
                      )        
    }
    ## cross validation over models or over kernels
    cv_type <- c('kernel', 'model')[idx]
    ## read in files 
    result <- lapply(list.files(temp_dir),
                     function(x) {
                         readRDS(paste(temp_dir, x, sep = '/'))
                     })
    ## remove temporary directory
    unlink(temp_dir,force=T)
    ## create model matrix for all models or kernels
    if (cv_type == 'model') {
        dt2 <- lapply(unique(kernel_bw$Var1),
                  function(x) {
                      vars <-  mods[[x]]
                      frml <- paste0('~', paste(paste0('log(', vars,')'), collapse = '+'),
                                     '+factor(code_orig)-1')
                      out <- data.table(code_orig = sp_dt@data$code_orig,
                                        fold_i = sp_dt@data$fold_i,
                                        model.matrix(as.formula(frml),
                                                     data = sp_dt@data
                                                     )
                                        )
                      setkey(out, code_orig, fold_i)
                      return(out)
                  })
    } else if (cv_type == 'kernel') {
        vars <-  mods[[15]]
        frml <- paste0('~', paste(paste0('log(', vars,')'), collapse = '+'),
                       '+factor(code_orig)-1')
        dt2 <- data.table(code_orig = sp_dt@data$code_orig,
                          fold_i = sp_dt@data$fold_i,
                          model.matrix(as.formula(frml),
                                       data = sp_dt@data
                                       )
                          )
        setkey(dt2, code_orig, fold_i)
    }
    ## calculate rmse and mae
    pred_errors <- lapply(seq_along(result),
                          function(x) {                             
                              ## find row number of input grid
                              ## k <- result[[x]][[3]] ## for local run
                              k <- result[[x]][[1]][[3]]
                              ## kernel name and index
                              kernel_type <- paste(unlist(kernel_bw[x, 1]))
                              dt2_index <- grep(kernel_type,
                                                paste(unlist(unique(kernel_bw[,1]))))
                              ## for out of sample err: gather parms for each fold
                              dt1 <- setDT(merge(sp_dt@data[,
                                                      c('code_orig',
                                                        'fold_i')
                                                      ],
                                           rbindlist(result[[x]][[1]][[2]]),
                                           ## rbindlist(result[[x]][[2]]), ## for local run
                                           by.x = c('code_orig',  'fold_i'),
                                           by.y = c('code_orig',  'fld')
                                           ))
                              setkey(dt1, code_orig, fold_i)                            
                              ## out of sample predictions
                              if (cv_type == 'model') {
                                  os_pred <- exp(apply(dt1[,-c(1,2)] * dt2[[dt2_index]][,-c(1,2)],
                                                   1,
                                                   sum))
                              } else if (cv_type == 'kernel') {
                                  os_pred <- exp(apply(dt1[,-c(1,2)] * dt2[,-c(1,2)],
                                                   1,
                                                   sum))
                              }
                              ## calculate within sample errors
                              ws_err <- rbindlist(result[[x]][[1]][[1]]
                              ##  ws_err <- rbindlist(result[[x]][[1]], ## for local run
                                                  )[,
                                                     .(rmse = sqrt(mean(residual^2)),
                                                       mae = mean(abs(residual))
                                                       ),
                                                     by = fld
                                                     ]
                              return(
                                  data.table(predicted = os_pred,
                                             actual = sp_dt@data$V1_int,
                                             fold = sp_dt@data$fold_i
                                             )[, ## calculate mse and mae
                                               .(sqrt(mean((predicted - actual)^2)),
                                                 mean(abs(predicted - actual))
                                                 ),
                                               by = fold
                                               ][, ## return mean and se of mse and mae
                                                 .(idx = k,
                                                   kernel = unlist(kernel_bw[k, 1]),
                                                   bw = unlist(kernel_bw[k, 2]),
                                                   mean_rmse_oos = mean(V1),
                                                   se_rmse_oos = sd(V1)/sqrt(.N),
                                                   mean_mae_oos = mean(V2),
                                                   se_mae_oos = sd(V2)/sqrt(.N),
                                                   mean_rmse_ws = mean(ws_err$rmse),
                                                   se_rmse_ws = sd(ws_err$rmse)/sqrt(.N),
                                                   mean_mae_ws = mean(ws_err$mae),
                                                   se_mae_ws = sd(ws_err$mae)/sqrt(.N)
                                                   )
                                                 ]
                              )
                          })
    return(rbindlist(pred_errors))
}

pred_errors <- getJobResult("calculate-errors-eureka-13")

pred_errors[[1]][mean_rmse_oos == min(mean_rmse_oos)]
## old


ggplot(rbindlist(pred_errors), aes(bw, mean_rmse_oos)) +
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
