bw.ggwr2 <- function (formula, data, family = "poisson", approach = "CV", 
    kernel = "bisquare", adaptive = FALSE, p = 2, theta = 0, 
    longlat = F, dMat) 
{
    if (is(data, "Spatial")) {
        dp.locat <- coordinates(data)
        data <- as(data, "data.frame")
    }
    else {
        stop("Given regression data must be Spatial*DataFrame")
    }
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    dp.n <- nrow(data)
    if (dp.n > 1500) {
        cat("Take a cup of tea and have a break, it will take a few minutes.\n")
        cat("          -----A kind suggestion from GWmodel development group\n")
    }
    if (missing(dMat)) 
        DM.given <- F
    else {
        DM.given <- T
        dim.dMat <- dim(dMat)
        if (dim.dMat[1] != dp.n || dim.dMat[2] != dp.n) 
            stop("Dimensions of dMat are not correct")
    }
    if (adaptive) {
        upper <- dp.n
        lower <- 20
    }
    else {
        if (DM.given) {
            upper <- range(dMat)[2]
            lower <- upper/5000
        }
        else {
            dMat <- NULL
            if (p == 2) {
                b.box <- bbox(dp.locat)
                upper <- sqrt((b.box[1, 2] - b.box[1, 1])^2 + 
                  (b.box[2, 2] - b.box[2, 1])^2)
                lower <- upper/5000
            }
            else {
                upper <- 0
                for (i in 1:dp.n) {
                  dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, 
                    p = p, theta = theta, longlat = longlat)
                  upper <- max(upper, range(dist.vi)[2])
                }
                lower <- upper/5000
            }
        }
    }
    bw <- NA
    if (approach == "cv" || approach == "CV") 
        bw <- gold(ggwr.cv2, lower, upper, adapt.bw = adaptive, 
            x, y, family = family, kernel, adaptive, dp.locat, 
            p, theta, longlat, dMat)
    else if (approach == "aic" || approach == "AIC" || approach == 
        "AICc") 
        bw <- gold(ggwr.aic, lower, upper, adapt.bw = adaptive, 
            x, y, family = family, kernel, adaptive, dp.locat, 
            p, theta, longlat, dMat)
    bw
}


ggwr.cv2 <- function (bw, X, Y, family = "poisson", kernel = "bisquare", 
    adaptive = F, dp.locat, p = 2, theta = 0, longlat = F, dMat) 
{
    dp.n <- length(dp.locat[, 1])
    if (is.null(dMat)) 
        DM.given <- F
    else {
        DM.given <- T
        dim.dMat <- dim(dMat)
        if (dim.dMat[1] != dp.n || dim.dMat[2] != dp.n) 
            stop("Dimensions of dMat are not correct")
    }
    CV <- numeric(dp.n)
    Wt <- matrix(numeric(dp.n * dp.n), ncol = dp.n)
    for (i in 1:dp.n) {
        if (DM.given) 
            dist.vi <- dMat[, i]
        else {
            dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, 
                p = 2, theta = theta, longlat = longlat)
        }
        W.i <- gw.weight(dist.vi, bw, kernel, adaptive)
        W.i[i] <- 0
        Wt[, i] <- W.i
    }
    wt2 <- rep(1, dp.n)
    if (family == "poisson") {
        res1 <- gwr.poisson.wt2(Y, X, bw, Wt)
        wt2 <- res1[[1]]
        y.adj <- res1[[3]]
    }
    else if (family == "binomial") {
        res1 <- gwr.binomial.wt(Y, X, bw, Wt)
        wt2 <- res1[[1]]
        y.adj <- res1[[3]]
    }
    for (i in 1:dp.n) {
        W.i <- Wt[, i] * wt2
        gwsi <- try(gw_reg(X, y.adj, W.i, FALSE, i))
        if (!inherits(gwsi, "try-error")) {
            yhat.noi <- X[i, ] %*% (gwsi[[1]])
            if (family == "poisson") 
                CV[i] <- Y[i] - exp(yhat.noi)
            else if (family == "binomial") 
                CV[i] <- Y[i] - exp(yhat.noi)/(1 + exp(yhat.noi))
        }
        else {
            CV[i] <- Inf
            break
        }
    }
    if (!any(is.infinite(CV))) 
        CV.score <- t(CV) %*% CV
    else {
        CV.score <- Inf
    }
    if (adaptive) 
        cat("Adaptive bandwidth:", bw, "CV score:", CV.score, 
            "\n")
    else cat("Fixed bandwidth:", bw, "CV score:", CV.score, "\n")
    CV.score
}

gwr.poisson.wt2 <- function (y, x, bw, W.mat, verbose = T) 
{   
    tol <- 1e-05
    maxiter <- 20
    var.n <- ncol(x)
    dp.n <- nrow(x)
    betas <- matrix(nrow = dp.n, ncol = var.n)
    betas1 <- betas
    S <- matrix(nrow = dp.n, ncol = dp.n)
    it.count <- 0
    llik <- 0
    mu <- y + 0.1
    nu <- log(mu)
    if (verbose) 
        cat(" Iteration    Log-Likelihood(With bandwidth: ", 
            bw, ")\n=========================\n")
    wt2 <- rep(1, dp.n)
    repeat {
        y.adj <- nu + (y - mu)/mu
        for (i in 1:dp.n) {
            W.i <- W.mat[, i]
            gwsi <- gw_reg(x, y.adj, W.i * wt2, FALSE, i)
            betas1[i, ] <- gwsi[[1]]
        }
        nu <- gw.fitted(x, betas1)
        mu <- exp(nu)
        old.llik <- llik
        #llik <- sum(y * nu - mu - log(gamma(y + 1)))        
        llik <- sum(dpois(y, mu, log = TRUE))
        if (verbose) 
            cat(paste("   ", formatC(it.count, digits = 4, width = 4), 
                "    ", formatC(llik, digits = 4, width = 7), 
                "\n"))
        if (abs((old.llik - llik)/llik) < tol) 
            break
        wt2 <- mu
        it.count <- it.count + 1
        if (it.count == maxiter) 
            break
    }
    res <- list(wt2, llik, y.adj)
    res
}

ggwr.basic2 <-function(formula, data, regression.points, bw, family ="poisson", kernel="bisquare",
              adaptive=FALSE, cv=F, tol=1.0e-5, maxiter=20, p=2, theta=0, longlat=F, dMat,dMat1, no.hatmatrix=FALSE, theta_g=1, null.dev=1)
{
 ##Record the start time
    timings <- list()
    timings[["start"]] <- Sys.time()
    ###################################macth the variables
    this.call <- match.call()
    p4s <- as.character(NA)
    #####Check the given data frame and regression points
    #####Regression points
    if (missing(regression.points))
    {
        rp.given <- FALSE
        regression.points <- data
        hatmatrix<-T
    }
    else
    {
        rp.given <- TRUE
        hatmatrix<-F
    }   
    ##Data points{
    if (is(data, "Spatial"))
    {
        p4s <- proj4string(data)
        dp.locat<-coordinates(data)
        data <- as(data, "data.frame")
    }
    else
    {
        stop("Given regression data must be Spatial*DataFrame")
    }

    ####################
    ######Extract the data frame
    ####Refer to the function lm
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)

    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    ##test_x <<- x
    ##test_y <<- y 
    ############################################
    var.n<-ncol(x)
    if(is(regression.points, "Spatial"))
        rp.locat<-coordinates(regression.points)
    else if(is.numeric(regression.points)&&dim(regression.points)[2]==2)
    {
        rp.locat <- regression.points
    }
    else
        stop("Please use the correct regression points for model calibration!")

    rp.n<-nrow(rp.locat)
    dp.n<-nrow(data)
    betas <-matrix(nrow=rp.n, ncol=var.n)
    betas1<- betas
    if(hatmatrix)
    {
        betas.SE <-matrix(nrow=rp.n, ncol=var.n)
        betas.TV <-matrix(nrow=rp.n, ncol=var.n)
        ##S: hatmatrix
        S<-matrix(nrow=dp.n,ncol=dp.n)
    }
    #C.M<-matrix(nrow=dp.n,ncol=dp.n)
    idx1 <- match("(Intercept)", colnames(x))
    if(!is.na(idx1))
        colnames(x)[idx1]<-"Intercept"
    colnames(betas) <- colnames(x)
    #colnames(betas)[1]<-"Intercept"
    ####################################################GWR
    #########Distance matrix is given or not

    if (missing(dMat))
    {
        DM.given<-F
        if(dp.n + rp.n <= 10000)
        {
            dMat <- gw.dist(dp.locat=dp.locat, rp.locat=rp.locat, p=p, theta=theta, longlat=longlat)
            DM.given<-T
        }
    }
    else
    {
        DM.given<-T
        dim.dMat<-dim(dMat)
        if (dim.dMat[1]!=dp.n||dim.dMat[2]!=rp.n)
            stop("Dimensions of dMat are not correct")
    }
    if(missing(dMat1))
    {
        DM1.given<-F
        if(hatmatrix&&DM.given)
        {
            dMat1 <- dMat
            DM1.given<-T
        }
        else
        {
            if(dp.n < 8000)
            {
                dMat1 <- gw.dist(dp.locat=dp.locat, rp.locat=dp.locat, p=p, theta=theta, longlat=longlat)
                DM1.given<-T
            }
        }
    }
    else
    {
        DM1.given<-T
        dim.dMat1<-dim(dMat1)
        if (dim.dMat1[1]!=dp.n||dim.dMat1[2]!=dp.n)
            stop("Dimensions of dMat are not correct")
    }
    ####Generate the weighting matrix
    #############Calibration the model
    W1.mat<-matrix(numeric(dp.n*dp.n),ncol=dp.n)
    W2.mat<-matrix(numeric(dp.n*rp.n),ncol=rp.n)
    for (i in 1:dp.n)
    {
        if (DM1.given)
            dist.vi<-dMat1[,i]
        else
        {
            dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
        }
        W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
        W1.mat[,i]<-W.i
    }
    if (rp.given)
   { 
        for (i in 1:rp.n)
        {
            if (DM.given)
                dist.vi<-dMat[,i]
            else
            {
                dist.vi<-gw.dist(dp.locat, rp.locat, focus=i, p, theta, longlat)
            }
            W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
            W2.mat[,i]<-W.i
        }
    }
    else
        W2.mat<-W1.mat
    ##model calibration
    if (no.hatmatrix) hatmatrix <- FALSE
    if(family=="poisson")
        res1<-gwr.poisson2(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol, maxiter)
    else if(family=="poisson.fe")
        res1<-gwr.poisson.fe(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol, maxiter)
    else if(family=="binomial")
        res1<-gwr.binomial(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol, maxiter)
    else if(family=="nbinomial.fe")
        res1<-gwr.nbinomial.fe(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol, maxiter)
    else if(family=="nbinomial.fe.global")
        res1<-gwr.nbinomial.fe(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol, maxiter, theta_spec = 'global', theta_g=theta_g, null.dev=null.dev)
    else if(family=="nbinomial.fe.fixed")
        res1<-gwr.nbinomial.fe(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol, maxiter, theta_spec = 'fixed')
    ####################################
    CV <- numeric(dp.n)
    if(hatmatrix && cv)
    {
        CV <- ggwr.cv.contrib2(bw, x, y,family, kernel,adaptive, dp.locat, p, theta, longlat,dMat)
    }
    ####encapsulate the GWR results
    GW.arguments<-list()
    GW.arguments<-list(formula=formula,rp.given=rp.given,hatmatrix=hatmatrix,bw=bw, family=family,
                       kernel=kernel,adaptive=adaptive, p=p, theta=theta, longlat=longlat,DM.given=DM1.given)

    timings[["stop"]] <- Sys.time()
    ##############
    res<-list(GW.arguments=GW.arguments,GW.diagnostic=res1$GW.diagnostic,glms=res1$glms,SDF=res1$SDF,CV=CV,timings=timings,this.call=this.call)
    class(res) <-"ggwrm"
    invisible(res) 
}

############ Possion GWGLM
gwr.poisson2<-function(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol=1.0e-5, maxiter=500, fe.reg = T, fe.col = 1, fe.cutoff = 1.0e-8)
{
    p4s <- as.character(NA)
    if (is(regression.points, "Spatial"))
    {
      p4s <- proj4string(regression.points)
    }
    ############################################
    ##Generalized linear regression
    glms<-glm.fit(x, y, family = poisson()) 
    null.dev <- glms$null.deviance
    glm.dev <-glms$deviance
    glm.pseudo.r2 <- 1- glm.dev/null.dev 
    glms$pseudo.r2 <- glm.pseudo.r2
    var.n<-ncol(x)
    dp.n<-nrow(x)
    ########change the aic
    glms$aic <- glm.dev + 2*var.n
    glms$aicc <- glm.dev + 2*var.n + 2*var.n*(var.n+1)/(dp.n-var.n-1)
    ############################################
    if(is(regression.points, "Spatial"))
    	 rp.locat<-coordinates(regression.points)
    else
       rp.locat <- regression.points
    rp.n<-nrow(rp.locat)
    betas <- matrix(0, nrow=rp.n, ncol=var.n)
    betas1 <- matrix(0, nrow=dp.n, ncol=var.n)
    betas.SE <-matrix(0, nrow=dp.n, ncol=var.n)
    betas.TV <-matrix(0, nrow=dp.n, ncol=var.n)
    ##S: hatmatrix
    S<-matrix(nrow=dp.n,ncol=dp.n)
    #C.M<-matrix(nrow=dp.n,ncol=dp.n)
    colnames(betas) <- colnames(x)
    ## colnames(betas)[1]<-"Intercept"  
    ## fixed effects regression: only include fe-dummies with weight > tol
    if (fe.reg) {
        fe.col.name <- names(regression.points@data)[fe.col]
        fe.cols <- unique(regression.points@data[, fe.col])
        iv.cols <- which(!grepl('factor', colnames(x)))
        x_index <- 1:dp.n
        W1.sub <- list()
        W2.sub <- list()
        wt2.sub <- list()
        x.sub <- list()
        gwr.sub <- list()
        fe.sub.cols <- list()
        for (i in 1:dp.n) {
            ## find rows where weight > tol
            gwr.sub[[i]] <- which(W1.mat[,i] > fe.cutoff)
            ## find corresponding fixed effects and columns
            fe.sub <- unique(regression.points@data[gwr.sub[[i]], fe.col])
            fe.sub.cols[[i]] <- which(fe.cols %in% fe.sub) + length(iv.cols)
            ## subset weight matrix
            W1.sub[[i]] <- W1.mat[gwr.sub[[i]], i]
            W2.sub[[i]] <- W2.mat[gwr.sub[[i]], i]
            ## subset x
            x.sub[[i]] <- cbind(x[gwr.sub[[i]], iv.cols],
                                model.matrix(~ factor(get(fe.col.name)) -1,
                                             regression.points@data[gwr.sub[[i]], ]
                                             )
                                )
        }
    }
    ####################################
    ##model calibration
    it.count <- 0
    llik <- 0.0
    mu <- y + 0.1
    nu <- log(mu)
    cat(paste("Begin LL iterations", Sys.time(), "\n\n"))
    cat(" Iteration    Log-Likelihood\n=========================\n")
    wt2 <- rep(1,dp.n)
    repeat {
     y.adj <- nu + (y - mu)/mu
     for (i in 1:dp.n)
     {
         W.i<-W1.mat[,i]
         if (fe.reg) {
             y.sub <- y.adj[gwr.sub[[i]]]
             wt2.sub <- wt2[gwr.sub[[i]]]
             gwsi<-gw_reg(x.sub[[i]],y.sub,W1.sub[[i]]*wt2.sub,hatmatrix=F,i)
             betas1[i,c(iv.cols, fe.sub.cols[[i]])]<-gwsi[[1]]
         } else {
             gwsi<-gw_reg(x,y.adj,W.i*wt2,hatmatrix=F,i)
             betas1[i,]<-gwsi[[1]]
         }         
     }
     nu <- gw.fitted(x,betas1)
     mu <- exp(nu)
     ##nu_test <<- nu
     ##mu_test <<- mu
     old.llik <- llik
     #llik <- sum(y*nu - mu - log(gamma(y+1)))
     #print(paste0('original: ', llik))

     ################################ edits #####################

     llik <- sum(dpois(y, mu, log = TRUE))
     #print(paste0('new: ', llik))
     ################################ edits #####################
     
     cat(paste("   ",formatC(it.count,digits=4,width=4),"    ",formatC(llik,digits=4,width=7),"\n"))
     if (abs((old.llik - llik)/llik) < tol) break
     wt2 <- as.numeric(mu)
     it.count <- it.count+1
     if (it.count == maxiter) break
    }
    cat(paste("End of LL iterations", Sys.time(), "\n"))
    GW.diagnostic <- NULL
    gw.dev <- 0
    for(i in 1:dp.n) {
        if(y[i]!=0)
            gw.dev <- gw.dev + 2*(y[i]*(log(y[i]/mu[i])-1)+mu[i])
        else
            gw.dev <- gw.dev + 2* mu[i]
    }

     #gw.dev <- 2*sum(y*log(y/mu)-(y-mu))     
     #local.dev <- numeric(dp.n)     
     #local.null.dev <- numeric(dp.n)
     #local.pseudo.r2 <- numeric(dp.n) 
    if(hatmatrix) {
        for (i in 1:dp.n) {
            if (fe.reg) {
                y.sub <- y.adj[gwr.sub[[i]]]
                wt2.sub <- wt2[gwr.sub[[i]]]
                focus_i <- which(x_index[gwr.sub[[i]]] == i)       
                gwsi<-gw_reg(x.sub[[i]],y.sub,W2.sub[[i]]*wt2.sub,hatmatrix,focus_i)
                betas[i, c(iv.cols, fe.sub.cols[[i]])]<-gwsi[[1]]
                S[i,gwr.sub[[i]]]<-gwsi[[2]]
                Ci<-gwsi[[3]]
                invwt2 <- 1.0 /as.numeric(wt2.sub)
                betas.SE[i, c(iv.cols, fe.sub.cols[[i]])] <- diag((Ci*invwt2) %*% t(Ci))
            } else {
                W.i<-W2.mat[,i]
                gwsi<-gw_reg(x,y.adj,W.i*wt2,hatmatrix,i)
                betas[i,]<-gwsi[[1]]
                ##Add the smoother y.adjust, see equation (30) in Nakaya(2005)
                S[i,]<-gwsi[[2]]
                Ci<-gwsi[[3]]
                                        #betas.SE[i,]<-diag(Ci%*%t(Ci))
                invwt2 <- 1.0 /as.numeric(wt2)
                betas.SE[i,] <- diag((Ci*invwt2) %*% t(Ci))# diag(Ci/wt2%*%t(Ci))  #see Nakaya et al. (2005)
            }
            
        }
        cat(paste("End of hatmatrix loop", Sys.time(), "\n"))
        tr.S<-sum(diag(S))
        ####trace(SWS'W^-1) is used here instead of tr.StS
        tr.StS<-sum(S^2)
        ##tr.StS<- sum(diag(S%*%diag(wt2)%*%t(S)%*% diag(1/wt2))) ## NB:: really slow!!
        ###edf is different from the definition in Chris' code
        #edf<-dp.n-2*tr.S+tr.StS
        yhat<-gw.fitted(x, betas)
        residual<-y-exp(yhat)
        ########rss <- sum((y - gwr.fitted(x,b))^2)
        #rss <- sum((y-exp(yhat))^2)
        #sigma.hat <- rss/edf
        #sigma.aic <- rss/dp.n
        cat(paste("Calculated fit and residuals", Sys.time(),  "\n"))
        for(i in 1:dp.n) {
            ##betas.SE[i,]<-sqrt(sigma.hat*betas.SE[i,])
            if (fe.reg) {
                betas.SE[i,c(iv.cols, fe.sub.cols[[i]])]<-sqrt(betas.SE[i,c(iv.cols, fe.sub.cols[[i]])])
                betas.TV[i,c(iv.cols, fe.sub.cols[[i]])]<-betas[i,c(iv.cols, fe.sub.cols[[i]])]/betas.SE[i,c(iv.cols, fe.sub.cols[[i]])]
            } else {
                betas.SE[i,]<-sqrt(betas.SE[i,])
                betas.TV[i,]<-betas[i,]/betas.SE[i,] 
            }
        }
        #AICc <- -2*llik + 2*tr.S*dp.n/(dp.n-tr.S-2) 
        #AICc <- -2*llik + 2*tr.S + 2*tr.S*(tr.S+1)/(dp.n-tr.S-1)  # This is generic form of AICc (TN)
        AIC <- gw.dev + 2*tr.S
        AICc <- gw.dev + 2*tr.S + 2*tr.S*(tr.S+1)/(dp.n-tr.S-1) 
        #yss.g <- sum((y - mean(y))^2)
        #gw.R2<-1-rss/yss.g; ##R Square valeu
        #gwR2.adj<-1-(1-gw.R2)*(dp.n-1)/(edf-1) #Adjusted R squared valu
        
        pseudo.R2 <- 1- gw.dev/null.dev
        GW.diagnostic<-list(gw.deviance=gw.dev,AICc=AICc,AIC=AIC,pseudo.R2 =pseudo.R2)
        cat(paste("End of hatmatrix calculations", Sys.time(),  "\n"))
     }
     else
     {
        for (i in 1:rp.n)
        {
            if (fe.reg) {
                y.sub <- y.adj[gwr.sub[[i]]]
                wt2.sub <- wt2[gwr.sub[[i]]]
                focus_i <- which(x_index[gwr.sub[[i]]] == i)
                gwsi<-gw_reg(x.sub[[i]],y.sub,W2.sub[[i]]*wt2.sub,hatmatrix,focus_i)
                betas[i, c(iv.cols, fe.sub.cols[[i]])]<-gwsi[[1]]
            } else {
                W.i<-W2.mat[,i]
                gwsi<-gw_reg(x,y.adj,W.i*wt2,hatmatrix,i)
                betas[i,]<-gwsi[[1]] ######See function by IG
            }
        }     
     }
    if (hatmatrix)                                         
    {
      gwres.df<-data.frame(betas,y,exp(yhat),residual,betas.SE,betas.TV)
      colnames(gwres.df)<-c(c(c(colnames(betas),c("y","yhat","residual")),
                              paste(colnames(betas), "SE", sep="_")),
                            paste(colnames(betas), "TV", sep="_"))
    }
    else
    {
      gwres.df<-data.frame(betas)
    }
    rownames(rp.locat)<-rownames(gwres.df)
    griddedObj <- F
    if (is(regression.points, "Spatial"))
    { 
        if (is(regression.points, "SpatialPolygonsDataFrame"))
        {
           polygons<-polygons(regression.points)
           #SpatialPolygons(regression.points)
           #rownames(gwres.df) <- sapply(slot(polygons, "polygons"),
                              #  function(i) slot(i, "ID"))
           SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=gwres.df, match.ID=F)
        }
        else
        {
           griddedObj <- gridded(regression.points)
           SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gwres.df, proj4string=CRS(p4s), match.ID=F)
           gridded(SDF) <- griddedObj 
        }
    }
    else
        SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gwres.df, proj4string=CRS(p4s), match.ID=F)
##############
    if(hatmatrix)
      res<-list(GW.diagnostic=GW.diagnostic,glms=glms,SDF=SDF)
    else
      res <- list(glms=glms,SDF=SDF)
}

ggwr.cv.contrib2 <- function (bw, X, Y, family = "poisson", kernel = "bisquare", 
    adaptive = F, dp.locat, p = 2, theta = 0, longlat = F, dMat) 
{
    dp.n <- length(dp.locat[, 1])
    if (is.null(dMat)) 
        DM.given <- F
    else {
        DM.given <- T
        dim.dMat <- dim(dMat)
        if (dim.dMat[1] != dp.n || dim.dMat[2] != dp.n) 
            stop("Dimensions of dMat are not correct")
    }
    CV <- numeric(dp.n)
    Wt <- matrix(numeric(dp.n * dp.n), ncol = dp.n)
    for (i in 1:dp.n) {
        if (DM.given) 
            dist.vi <- dMat[, i]
        else {
            dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, 
                p = p, theta = theta, longlat = longlat)
        }
        W.i <- gw.weight(dist.vi, bw, kernel, adaptive)
        W.i[i] <- 0
        Wt[, i] <- W.i
    }
    wt2 <- rep(1, dp.n)
    if (family == "poisson") {
        res1 <- gwr.poisson.wt2(Y, X, bw, Wt, verbose = F)
        wt2 <- res1[[1]]
        y.adj <- res1[[3]]
    }
    else if (family == "binomial") {
        res1 <- gwr.binomial.wt(Y, X, bw, Wt, verbose = F)
        wt2 <- res1[[1]]
        y.adj <- res1[[3]]
    }
    for (i in 1:dp.n) {
        W.i <- Wt[, i] * wt2
        gwsi <- try(gw_reg(X, y.adj, W.i, FALSE, i))
        if (!inherits(gwsi, "try-error")) {
            yhat.noi <- X[i, ] %*% gwsi[[1]]
            if (family == "poisson") 
                CV[i] <- Y[i] - exp(yhat.noi)
            else if (family == "binomial") 
                CV[i] <- Y[i] - exp(yhat.noi)/(1 + exp(yhat.noi))
        }
        else {
            CV[i] <- Inf
            break
        }
    }
    CV
}

## llik.nbinom <- function(y_i, mu_i, theta_i) {
##                 sum(-dnbinom(x=y_i,size=theta_i,mu=mu_i,log=TRUE))
## }

llik.nbinom <- function(theta_i, mu_i, y_i) {
    -sum((lgamma(theta_i +
                 y_i) - lgamma(theta_i) - lgamma(y_i + 1) + theta_i * log(theta_i) + y_i *
          log(mu_i + (y_i == 0)) - (theta_i + y_i) * log(theta_i + mu_i)))
}


gwr.nbinomial.fe <- function(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol=1.0e-5, maxiter=500, fe.col = 1, theta_g = 1, null.dev = 1, theta_spec = 'local', fe.cutoff = 0)
{
    p4s <- as.character(NA)
    if (is(regression.points, "Spatial"))
    {
        p4s <- proj4string(regression.points)
    }
    ############################################
    ##Generalized linear regression
    ## glms<-glm.fit(x, y, family = binomial())
    ## null.dev <- glms$null.deviance
    ## glm.dev <-glms$deviance
    ## glm.pseudo.r2 <- 1- glm.dev/null.dev
    ## glms$pseudo.r2 <- glm.pseudo.r2
    var.n<-ncol(x)
    dp.n<-nrow(x)
    ## glms$aic <- glm.dev + 2*var.n
    ## glms$aicc <- glm.dev + 2*var.n + 2*var.n*(var.n+1)/(dp.n-var.n-1)
    ############################################
    ##Generalized linear regression: obtain global theta
    ##glms<-glm.nb(y ~ matrix(x, ncol = ncol(x)) - 1)
    ##theta_g <- glms$theta
    ##null.dev <- glms$null.deviance
    ############################################
    rp.locat<-coordinates(regression.points)
    rp.n<-nrow(rp.locat)
    betas <-matrix(0,nrow=rp.n, ncol=var.n)
    betas1<- matrix(0,nrow=dp.n, ncol=var.n)
    betas.SE <-matrix(0,nrow=rp.n, ncol=var.n)
    betas.TV <-matrix(0,nrow=rp.n, ncol=var.n)
    ##S: hatmatrix
    S<-matrix(nrow=dp.n,ncol=dp.n)
    #C.M<-matrix(nrow=dp.n,ncol=dp.n)
    colnames(betas) <- colnames(x)
    #colnames(betas)[1]<-"Intercept"
    ## initialisation fixed effects
    fe.col.name <- names(regression.points@data)[fe.col]
    fe.cols <- unique(regression.points@data[, fe.col])
    iv.cols <- which(!grepl('factor', colnames(x)))
    x_index <- 1:dp.n
    W1.sub <- list()
    W2.sub <- list()
    x.sub <- list()
    gwr.sub <- list()
    fe.sub.cols <- list()
    ## fixed effects regression: only include fe-dummies with weight > cutoff
    for (i in 1:dp.n) {
        if (i>1 & any(rp.locat[i,] == rp.locat[i-1,])) {
            gwr.sub[[i]] <- gwr.sub[[i-1]]
            fe.sub.cols[[i]] <- fe.sub.cols[[i-1]]
            W1.sub[[i]] <- W1.sub[[i-1]]
            W2.sub[[i]] <- W2.sub[[i-1]]
            x.sub[[i]] <- x.sub[[i-1]]
        } else {
            ## find subset of observations where weight > cutoff
            gwr.sub[[i]] <- which(W1.mat[,i] > fe.cutoff)
            ## find non-zero fixed effects and corresponding columns
            fe.sub <- unique(regression.points@data[gwr.sub[[i]], fe.col])
            fe.sub.cols[[i]] <- which(fe.cols %in% fe.sub) + length(iv.cols)
            ## subset weight matrix
            W1.sub[[i]] <- W1.mat[gwr.sub[[i]], i]
            W2.sub[[i]] <- W2.mat[gwr.sub[[i]], i]
            ## subset x and create dummy matrix
            x.sub[[i]] <- cbind(x[gwr.sub[[i]], iv.cols],
                                model.matrix(~ factor(get(fe.col.name)) -1,
                                             regression.points@data[gwr.sub[[i]], ]
                                             )
                                )            
        }        
    }
    ####################################
    ##model calibration
    theta <- ifelse(theta_spec == 'fixed', rep(10^6, dp.n), rep(theta_g, dp.n))
    alpha  <- 1/theta
    cat(" Iteration    Log-Likelihood\n=========================\n")
    a <- rep(1,dp.n)
    it.count <- 0
    llik <- 0.0
    mu <-y + 0.1
    nu <- log(mu)
    repeat {
        ## calculate elements of Fisher Information matrix
        a <- mu / (1 + alpha * mu) +
            (y - mu)*(alpha*mu) / (1 + 2*alpha*mu + alpha**2 * mu**2)
        y.adj <- nu + (y - mu) / (a * (1 + alpha * mu))                     
        for (i in 1:dp.n) {
            if (i>1 & any(rp.locat[i,] == rp.locat[i-1,])) {
                betas1[i,c(iv.cols, fe.sub.cols[[i]])] <- betas1[i-1,c(iv.cols, fe.sub.cols[[i-1]])]
                if (theta_spec == 'local') {
                    theta[i] <- theta[i-1]
                    alpha[i] <- alpha[i-1]
                }
            } else {
                gwsi<-gw_reg(x.sub[[i]],y.adj[gwr.sub[[i]]],
                             W1.sub[[i]]*a[gwr.sub[[i]]],hatmatrix=F,i)
                ## beta_test[[i]] <<- gwsi[[1]]
                ## x_test[[i]] <<- W1.sub[[i]]*x.sub[[i]]
                ##nu_i <- apply(x.sub[[i]], 1, function(z) sum(z*gwsi[[1]]))
                nu_i <- W1.sub[[i]]*x.sub[[i]] %*% gwsi[[1]]
                mu_i <- exp(nu_i)
                if (theta_spec == 'local') {
                    llik_i <- optim(theta_g, lower = 0.0001, llik.nbinom,
                                    method = "L-BFGS-B", y_i = y[gwr.sub[[i]]],
                                    mu_i=mu_i)
                    ##theta_test[[i]] <<- llik_i$par
                    theta[i] <- llik_i$par
                    alpha[i] <- 1 / theta[i]
                }
                betas1[i,c(iv.cols, fe.sub.cols[[i]])]<-gwsi[[1]]
            }            
        }        
        nu <- gw.fitted(x,betas1)
        mu <- exp(nu)        
        ## nu_test <<- nu
        ## mu_test <<- mu
        ## betas_test <<- betas1
        ## x_test <<- x
        old.llik <- llik
        ##llik <- sum(dnbinom(y,size=theta,mu=mu,log=TRUE))
        llik <- -llik.nbinom(y_i=y, theta_i=theta, mu_i=mu)         
        cat(paste("   ",formatC(it.count,digits=4,width=4),"    ",formatC(llik,digits=4,width=7),"\n"))
        if (abs((old.llik - llik)/llik) < tol) break
        it.count <- it.count+1
        if (it.count == maxiter) break    
    }
    GW.diagnostic <- NULL
    gw.dev <- 0
    if(hatmatrix)
    {
        for (i in 1:rp.n) {
            if (i>1 & any(rp.locat[i,] == rp.locat[i-1,])) {
                betas[i, c(iv.cols, fe.sub.cols[[i]])] <- betas[i-1, c(iv.cols, fe.sub.cols[[i-1]])]
                S[i,gwr.sub[[i]]] <- S[i-1,gwr.sub[[i-1]]]
                betas.SE[i, c(iv.cols, fe.sub.cols[[i]])] <- betas.SE[i-1, c(iv.cols, fe.sub.cols[[i-1]])]
            } else {
                y.sub <- y.adj[gwr.sub[[i]]]
                focus_i <- which(x_index[gwr.sub[[i]]] == i)
                gwsi<-gw_reg(x.sub[[i]],y.sub,W2.sub[[i]]*a[gwr.sub[[i]]],hatmatrix,focus_i)
                betas[i, c(iv.cols, fe.sub.cols[[i]])]<-gwsi[[1]]
                S[i,gwr.sub[[i]]]<-gwsi[[2]]
                Ci<-gwsi[[3]]
                invwt2 <- 1.0 /as.numeric(a[i])
                betas.SE[i, c(iv.cols, fe.sub.cols[[i]])] <- diag((Ci*invwt2) %*% t(Ci))
            }
        }
        tr.S<-sum(diag(S))
        #tr.StS<-sum(S^2)
        #tr.StS<- sum(diag(S%*%diag(wt2)%*%t(S)%*% diag(1/wt2)))
        ###edf is different from the definition in Chris' code
        #edf<-dp.n-2*tr.S+tr.StS        
        yhat<-gw.fitted(x, betas)
        ##residual<-y-exp(yhat)/(1+exp(yhat))
        residual<-y-exp(yhat)
        ########rss <- sum((y - gwr.fitted(x,b))^2)
        rss <- sum(residual^2)
        #sigma.hat <- rss/edf
        #sigma.aic <- rss/dp.n   ### can be omitted? (TN)
        ##gw.dev <- sum(log(1/((y-n+exp(yhat)/(1+exp(yhat))))^2))
        gw.dev <- 2*sum(y*log(y/exp(yhat)) - (y+theta)*log((1+alpha*y)/(1+alpha*exp(yhat))))
        for(i in 1:dp.n) {
            if (i>1 & any(rp.locat[i,] == rp.locat[i-1,])) {
                betas.SE[i,c(iv.cols, fe.sub.cols[[i]])]<-betas.SE[i-1,c(iv.cols, fe.sub.cols[[i-1]])]
                betas.TV[i,c(iv.cols, fe.sub.cols[[i]])] <- betas.TV[i-1,c(iv.cols, fe.sub.cols[[i-1]])]
            } else {
                ##betas.SE[i,]<-sqrt(sigma.hat*betas.SE[i,])
                betas.SE[i,c(iv.cols, fe.sub.cols[[i]])]<-sqrt(betas.SE[i,c(iv.cols, fe.sub.cols[[i]])])
                betas.TV[i,c(iv.cols, fe.sub.cols[[i]])]<-betas[i,c(iv.cols, fe.sub.cols[[i]])]/betas.SE[i,c(iv.cols, fe.sub.cols[[i]])]
            }
        }
        #AICc <- -2*llik + 2*tr.S*dp.n/(dp.n-tr.S-2)
        #AICc <- -2*llik + 2*tr.S + 2*tr.S*(tr.S+1)/(dp.n-tr.S-1)
        AICc <- gw.dev + 2*tr.S + 2*tr.S*(tr.S+1)/(dp.n-tr.S-1)
        AIC <- gw.dev + 2*tr.S
        #yss.g <- sum((y - mean(y))^2)
        #gw.R2<-1-rss/yss.g; ##R Square valeu  ### is R2 needed? (TN)
        #gwR2.adj<-1-(1-gw.R2)*(dp.n-1)/(edf-1) #Adjusted R squared value
        pseudo.R2 <- 1 - gw.dev/null.dev
        #GW.diagnostic<-list(rss=rss,AICc=AICc,edf=edf,gw.R2=gw.R2,gwR2.adj=gwR2.adj)
        GW.diagnostic<-list(gw.deviance=gw.dev,AICc=AICc,AIC=AIC,pseudo.R2 =pseudo.R2)
    }
    else
    {
        for (i in 1:rp.n) {
            if (i>1 & any(rp.locat[i,] == rp.locat[i-1,])) {
                betas[i,] <- betas[i-1,]
            } else {
                W.i<-W2.mat[,i]
                y.sub <- y.adj[gwr.sub[[i]]]
                focus_i <- which(x_index[gwr.sub[[i]]] == i)
                gwsi<-gw_reg(x.sub[[i]],y.sub,W2.sub[[i]]*a[gwr.sub[[i]]],hatmatrix,focus_i)
                betas[i, c(iv.cols, fe.sub.cols[[i]])]<-gwsi[[1]]
            }
        }
        yhat<-gw.fitted(x, betas)
        residual<-y-exp(yhat)
    }
    if (hatmatrix)                                         
    {
      gwres.df<-data.frame(betas,alpha,y,exp(yhat),residual,betas.SE,betas.TV)
      colnames(gwres.df)<-c(c(c(colnames(betas),c("alpha", "y","yhat","residual")),paste(colnames(betas), "SE", sep="_")),paste(colnames(betas), "TV", sep="_"))
    }
    else
    {
        gwres.df<-data.frame(betas, y, exp(yhat),residual)
        colnames(gwres.df)<-c(c(colnames(betas),c("y","yhat","residual")))
    }
    rownames(rp.locat)<-rownames(gwres.df)
    griddedObj <- F   
    if (is(regression.points, "Spatial"))
    { 
        if (is(regression.points, "SpatialPolygonsDataFrame"))
        {
           polygons<-polygons(regression.points)
           #SpatialPolygons(regression.points)
           #rownames(gwres.df) <- sapply(slot(polygons, "polygons"),
                              #  function(i) slot(i, "ID"))
           SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=gwres.df, match.ID=F)
        }
        else
        {
           griddedObj <- gridded(regression.points)
           SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gwres.df, proj4string=CRS(p4s), match.ID=F)
           gridded(SDF) <- griddedObj 
        }
    }
    else
        SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gwres.df, proj4string=CRS(p4s), match.ID=F)
   ##############
    if(hatmatrix)
      res<-list(GW.diagnostic=GW.diagnostic,glms=glms,SDF=SDF)
    else
      res <- list(glms=alpha,SDF=SDF)
}

gwr.poisson.fe <-function(y,x,regression.points,W1.mat,W2.mat,hatmatrix,tol=1.0e-5, maxiter=500, fe.reg = T, fe.col = 1, fe.cutoff = 1.0e-8)
{
    p4s <- as.character(NA)
    if (is(regression.points, "Spatial"))
    {
      p4s <- proj4string(regression.points)
    }
    ############################################
    ##Generalized linear regression
    glms<-glm.fit(x, y, family = poisson()) 
    null.dev <- glms$null.deviance
    glm.dev <-glms$deviance
    glm.pseudo.r2 <- 1- glm.dev/null.dev 
    glms$pseudo.r2 <- glm.pseudo.r2
    var.n<-ncol(x)
    dp.n<-nrow(x)
    ########change the aic
    glms$aic <- glm.dev + 2*var.n
    glms$aicc <- glm.dev + 2*var.n + 2*var.n*(var.n+1)/(dp.n-var.n-1)
    ############################################
    if(is(regression.points, "Spatial"))
    	 rp.locat<-coordinates(regression.points)
    else
       rp.locat <- regression.points
    rp.n<-nrow(rp.locat)
    betas <- matrix(0, nrow=rp.n, ncol=var.n)
    betas1 <- matrix(0, nrow=dp.n, ncol=var.n)
    betas.SE <-matrix(0, nrow=dp.n, ncol=var.n)
    betas.TV <-matrix(0, nrow=dp.n, ncol=var.n)
    ##S: hatmatrix
    S<-matrix(nrow=dp.n,ncol=dp.n)
    #C.M<-matrix(nrow=dp.n,ncol=dp.n)
    colnames(betas) <- colnames(x)
    ## colnames(betas)[1]<-"Intercept"  
    ## fixed effects regression: only include fe-dummies with weight > cutoff
    ## initialisation
    fe.col.name <- names(regression.points@data)[fe.col]
    fe.cols <- unique(regression.points@data[, fe.col])
    iv.cols <- which(!grepl('factor', colnames(x)))
    x_index <- 1:dp.n
    W1.sub <- list()
    W2.sub <- list()
    wt2.sub <- list()
    x.sub <- list()
    gwr.sub <- list()
    fe.sub.cols <- list()
    for (i in 1:dp.n) {
        if (i>1 & any(coordinates(regression.points)[i,] == coordinates(regression.points)[i-1,])) {
            gwr.sub[[i]] <- gwr.sub[[i-1]]
            fe.sub.cols[[i]] <- fe.sub.cols[[i-1]]
            W1.sub[[i]] <- W1.sub[[i-1]]
            W2.sub[[i]] <- W2.sub[[i-1]]
            x.sub[[i]] <- x.sub[[i-1]]
        } else {
            ## find subset of observations where weight > cutoff
            gwr.sub[[i]] <- which(W1.mat[,i] > fe.cutoff)
            ## find non-zero fixed effects and corresponding columns
            fe.sub <- unique(regression.points@data[gwr.sub[[i]], fe.col])
            fe.sub.cols[[i]] <- which(fe.cols %in% fe.sub) + length(iv.cols)
            ## subset weight matrix
            W1.sub[[i]] <- W1.mat[gwr.sub[[i]], i]
            W2.sub[[i]] <- W2.mat[gwr.sub[[i]], i]
            ## subset x and create dummy matrix
            x.sub[[i]] <- cbind(x[gwr.sub[[i]], iv.cols],
                                model.matrix(~ factor(get(fe.col.name)) -1,
                                             regression.points@data[gwr.sub[[i]], ]
                                             )
                                )
        }
    }
    ####################################
    ##model calibration    
    it.count <- 0
    llik <- 0.0
    mu <- y + 0.1
    nu <- log(mu)
    cat(paste("Begin LL iterations", Sys.time(), "\n\n"))
    cat(" Iteration    Log-Likelihood\n=========================\n")
    wt2 <- rep(1,dp.n)
    repeat {
     y.adj <- nu + (y - mu)/mu
     for (i in 1:dp.n)
     {
         if (i>1 & any(coordinates(regression.points)[i,] == coordinates(regression.points)[i-1,])) {
             betas1[i,c(iv.cols, fe.sub.cols[[i]])] <- betas1[i-1,c(iv.cols, fe.sub.cols[[i-1]])]
         } else {            
             W.i<-W1.mat[,i]
             y.sub <- y.adj[gwr.sub[[i]]]
             wt2.sub <- wt2[gwr.sub[[i]]]
             gwsi<-gw_reg(x.sub[[i]],y.sub,W1.sub[[i]]*wt2.sub,hatmatrix=F,i)
             betas1[i,c(iv.cols, fe.sub.cols[[i]])]<-gwsi[[1]]
         }
     }
     nu <- gw.fitted(x,betas1)
     mu <- exp(nu)
     ## nu_test <<- nu
     ## mu_test <<- mu
     ## betas_test <<- betas1
     ## x_test <<- x
     old.llik <- llik
     #llik <- sum(y*nu - mu - log(gamma(y+1)))
     #print(paste0('original: ', llik))
     ################################ edits #####################
     llik <- sum(dpois(y, mu, log = TRUE))
     #print(paste0('new: ', llik))
     ################################ edits #####################
     cat(paste("   ",formatC(it.count,digits=4,width=4),"    ",formatC(llik,digits=4,width=7),"\n"))
     if (abs((old.llik - llik)/llik) < tol) break
     wt2 <- as.numeric(mu)
     it.count <- it.count+1
     if (it.count == maxiter) break
    }
    cat(paste("End of LL iterations", Sys.time(), "\n"))
    GW.diagnostic <- NULL
    gw.dev <- 0
    for(i in 1:dp.n) {
        if(y[i]!=0)
            gw.dev <- gw.dev + 2*(y[i]*(log(y[i]/mu[i])-1)+mu[i])
        else
            gw.dev <- gw.dev + 2* mu[i]
    }
     #gw.dev <- 2*sum(y*log(y/mu)-(y-mu))     
     #local.dev <- numeric(dp.n)     
     #local.null.dev <- numeric(dp.n)
     #local.pseudo.r2 <- numeric(dp.n) 
    if(hatmatrix) {
        for (i in 1:dp.n) {
            if (i>1 & any(coordinates(regression.points)[i,] == coordinates(regression.points)[i-1,])) {
                betas[i, c(iv.cols, fe.sub.cols[[i]])] <- betas[i-1, c(iv.cols, fe.sub.cols[[i-1]])]
                S[i,gwr.sub[[i]]] <- S[i-1,gwr.sub[[i-1]]]
                betas.SE[i, c(iv.cols, fe.sub.cols[[i]])] <- betas.SE[i-1, c(iv.cols, fe.sub.cols[[i-1]])]                
            } else {
                y.sub <- y.adj[gwr.sub[[i]]]
                wt2.sub <- wt2[gwr.sub[[i]]]
                focus_i <- which(x_index[gwr.sub[[i]]] == i)
                gwsi<-gw_reg(x.sub[[i]],y.sub,W2.sub[[i]]*wt2.sub,hatmatrix,focus_i)
                betas[i, c(iv.cols, fe.sub.cols[[i]])]<-gwsi[[1]]
                S[i,gwr.sub[[i]]]<-gwsi[[2]]
                Ci<-gwsi[[3]]
                invwt2 <- 1.0 /as.numeric(wt2.sub)
                betas.SE[i, c(iv.cols, fe.sub.cols[[i]])] <- diag((Ci*invwt2) %*% t(Ci))              
            }
        }
        cat(paste("End of hatmatrix loop", Sys.time(), "\n"))
        ## NB: something in the lines below is VERY time consuming
        tr.S<-sum(diag(S))
        ####trace(SWS'W^-1) is used here instead of tr.StS
        tr.StS<-sum(S^2) ## NB: commented out below because...time consuming 
        #tr.StS<- sum(diag(S%*%diag(wt2)%*%t(S)%*% diag(1/wt2)))
        ###edf is different from the definition in Chris' code
        #edf<-dp.n-2*tr.S+tr.StS
        yhat<-gw.fitted(x, betas)
        residual<-y-exp(yhat)
        ########rss <- sum((y - gwr.fitted(x,b))^2)
        #rss <- sum((y-exp(yhat))^2)
        #sigma.hat <- rss/edf
        #sigma.aic <- rss/dp.n
        for(i in 1:dp.n) {
            if (i>1 & any(coordinates(regression.points)[i,] == coordinates(regression.points)[i-1,])) {
                betas.SE[i,c(iv.cols, fe.sub.cols[[i]])]<-betas.SE[i-1,c(iv.cols, fe.sub.cols[[i-1]])]
                betas.TV[i,c(iv.cols, fe.sub.cols[[i]])] <- betas.TV[i-1,c(iv.cols, fe.sub.cols[[i-1]])]
            } else {
                betas.SE[i,c(iv.cols, fe.sub.cols[[i]])]<-sqrt(betas.SE[i,c(iv.cols, fe.sub.cols[[i]])])
            betas.TV[i,c(iv.cols, fe.sub.cols[[i]])]<-betas[i,c(iv.cols, fe.sub.cols[[i]])]/betas.SE[i,c(iv.cols, fe.sub.cols[[i]])]
            }
        }            
        #AICc <- -2*llik + 2*tr.S*dp.n/(dp.n-tr.S-2) 
        #AICc <- -2*llik + 2*tr.S + 2*tr.S*(tr.S+1)/(dp.n-tr.S-1)  # This is generic form of AICc (TN)
        AIC <- gw.dev + 2*tr.S
        AICc <- gw.dev + 2*tr.S + 2*tr.S*(tr.S+1)/(dp.n-tr.S-1) 
        #yss.g <- sum((y - mean(y))^2)
        #gw.R2<-1-rss/yss.g; ##R Square valeu
        #gwR2.adj<-1-(1-gw.R2)*(dp.n-1)/(edf-1) #Adjusted R squared valu
        pseudo.R2 <- 1- gw.dev/null.dev
        GW.diagnostic<-list(gw.deviance=gw.dev,AICc=AICc,AIC=AIC,pseudo.R2 =pseudo.R2)
        cat(paste("End of hatmatrix calculations", Sys.time(),  "\n"))
     }
     else
     {
        for (i in 1:rp.n)
        {
            if (i>1 & any(coordinates(regression.points)[i,] == coordinates(regression.points)[i-1,])) {
                betas[i,] <- betas[i-1,]
            } else {
                y.sub <- y.adj[gwr.sub[[i]]]
                wt2.sub <- wt2[gwr.sub[[i]]]
                focus_i <- which(x_index[gwr.sub[[i]]] == i)
                gwsi<-gw_reg(x.sub[[i]],y.sub,W2.sub[[i]]*wt2.sub,hatmatrix,focus_i)
                betas[i, c(iv.cols, fe.sub.cols[[i]])]<-gwsi[[1]]
                ## W.i<-W2.mat[,i]
                ## gwsi<-gw_reg(x,y.adj,W.i*wt2,hatmatrix,i)
                ## betas[i,]<-gwsi[[1]] ######See function by IG
            }            
        }
        yhat<-gw.fitted(x, betas)
        residual<-y-exp(yhat)
     }
    if (hatmatrix)                                         
    {
      gwres.df<-data.frame(betas,y,exp(yhat),residual,betas.SE,betas.TV)
      colnames(gwres.df)<-c(c(c(colnames(betas),c("y","yhat","residual")),paste(colnames(betas), "SE", sep="_")),paste(colnames(betas), "TV", sep="_"))
    }
    else
    {
        gwres.df<-data.frame(betas, y, exp(yhat),residual)
        colnames(gwres.df)<-c(c(colnames(betas),c("y","yhat","residual")))
    }
    rownames(rp.locat)<-rownames(gwres.df)
    griddedObj <- F
    if (is(regression.points, "Spatial"))
    { 
        if (is(regression.points, "SpatialPolygonsDataFrame"))
        {
           polygons<-polygons(regression.points)
           #SpatialPolygons(regression.points)
           #rownames(gwres.df) <- sapply(slot(polygons, "polygons"),
                              #  function(i) slot(i, "ID"))
           SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=gwres.df, match.ID=F)
        }
        else
        {
           griddedObj <- gridded(regression.points)
           SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gwres.df, proj4string=CRS(p4s), match.ID=F)
           gridded(SDF) <- griddedObj 
        }
    }
    else
        SDF <- SpatialPointsDataFrame(coords=rp.locat, data=gwres.df, proj4string=CRS(p4s), match.ID=F)
##############
    if(hatmatrix)
      res<-list(GW.diagnostic=GW.diagnostic,glms=glms,SDF=SDF)
    else
      res <- list(glms=glms,SDF=SDF)
}
