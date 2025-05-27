######################################################################
###################################################################### 
# Supplementary data providing requested code for performing analysis depicted in

# Garcia-Garzon, E., Abad, F.J. & Garrido, L.E. (in press). Improving Bi-factor Modelling: 
# Empirical Target Rotation based on Loadings Differences.
# Methodology Journal
#
# Functions based on psych, GPArotation packages and SLi function:
#
# Abad, F.J., Garcia-Garzon, E., Garrido, L.E. & Barrada, J.R. (2017).
# Iteration of Partially Specified Target Matrices: Application to the Bi-factor Case.
# Multivariate Behavioral Research, 52(4), 
# http://doi.org/10.1080/00273171.2017.1301244
#
#
# Revelle, W. (2017) psych: Procedures for Personality and Psychological Research, 
# Northwestern University,Evanston, Illinois, USA, 
# https://CRAN.R-project.org/package=psych Version = 1.7.8.
#
#
# Bernaards, Coen A. and Jennrich, Robert I. (2005) Gradient Projection Algorithms 
# and Software forArbitraryRotation Criteria in Factor Analysis, 
# Educational and Psychological Measurement: 65, 676-696.
# <http://www.stat.ucla.edu/research/gpa>
#
#
# Installation of GPArotation and Psych packages is recommended

######################################################################
######################################################################

######################################################################
# Support functions
######################################################################

######################################################################

vgQ.bifactor_geomin <- function (L, epsilon = 0.01) {
  
  # Bi-geomin criterion
  #
  # Args:
  #
  #   L: Unrotated factor solution
  #   epsilon = epsilon parameter
  #
  # Returns:
  #
  #   f  : value of the criterion at L
  #   Gq : Gradient at L
  #   
  #  
  # Error handling  
  
  D <- function(L) {
    L2 <- L^2 + epsilon  
    k <- ncol(L)
    p <- nrow(L)
    pro <- exp(rowSums(log(L2))/k)
    list(f = sum(pro),
         Gq = (2/k) * (L/L2) * matrix(rep(pro, k), p))
  }
  
  lvg <- D(L[, -1, drop = FALSE])
  G <- lvg$Gq
  G <- cbind(G[, 1], G)
  G[, 1] <- 0
  list(f = lvg$f, Gq = G)
  
}
######################################################################

######################################################################
perform_rotation <- function(L, method = NULL, targ = NULL, reps = 50) {
  # Compute factor rotation for geomin and target rotations
  #
  # Args:
  #
  #   L: Unrotated factor solution
  #   method: Factor rotation method, as in GPArotation:
  #           targetT: target
  #           geominQ: oblique geomin
  #           Default is NULL
  #   Targ: Target matrix required for target rotation. 
  #         Default is NULL.
  #   reps: Number of random starts. 
  #         Default is 50.
  #
  # Returns:
  #
  #   Loadings: Rotated factor loading matrix
  #   Phi     : Factor correlation matrix (only when geominQ is specified)
  #   Rotating matrix: Rotation matrix
  #  
  # Error handling  
  
  if (is.null(method)) {
    stop("A rotation method must be specified")
    cat ("\n")
  }
  
  # Compute the factor rotation using GPArotation package  
  results   <- rep( list(list()), reps) 
  criterion <- rep(NA,reps)
  
  for (i in 1:reps) {
    
    if (method == "targetT") {
      x <- GPArotation:::targetT (L, Tmat =GPArotation::: Random.Start(ncol(L)), maxit = 5000, Target = targ)
    }  else if (method == "geominQ") {
      x <- GPArotation:::geominQ (L, Tmat =GPArotation::: Random.Start(ncol(L)), maxit = 5000)
    }  else if (method == "bifactorT"){
      x  <- GPArotation:::bifactorT(L,Tmat = Random.Start(ncol(L)),maxit=5000)
    }  else if (method == "bifactorQ"){
      x  <- GPArotation:::bifactorQ(L,Tmat = Random.Start(ncol(L)),maxit=5000)
    } else if (method == "bigeominT"){
      x  <- GPArotation:::GPForth(L,method="bifactor_geomin",Tmat = Random.Start(ncol(L)),maxit=5000)
    } 
    
    # Selecting the best random start
    if (x$convergence == TRUE) {
      criterion[i] <- min(x$Table[,2])
      results[[i]] <- x
      
    } else { 
      criterion[i] <- NA
      results[[i]] <- NA
      cat("Convergence problem in factor rotation for random start ",i, "\n")
      cat ("\n")
    }
  }
  
  #return the best random start solution  
  
  j <- order(criterion)[1]
  return(results[[j]])
}
######################################################################

######################################################################
get_target_from <- function (L, cutpoint = "dif") {
  # Computes the target matrix 
  #
  # Args:
  #
  #   L: Factor loading matrix
  #   cutpoint: Cut-off for factor loading substantivity. 
  #             Default is "dif":  performs a promin-based estimation of the cut-off point
  #             Fixed cut-offs estimated by providing a value between 0 and 1
  #
  # Returns:
  #
  #   Targ = Partially Specified Target rotation
  #
  
  if (is.null(cutpoint)){
    stop("A cut-off point must be specified")
  }
  
  if (is.character(cutpoint)){
    if (cutpoint != "dif") {
      stop("for applying the SLiD, please write cutpoint = dif")
    }
  }
  
  Targ <- L
  
  if (cutpoint == "dif") {
    
    c2 <- c()
    c2 <- (Targ[, -1]/sqrt(apply(Targ[, -1]^2, 1, sum)))
    c2 <- c2^2
    
    i <- 1
    b <- list()
    dat <- list()
    cuts_mix <- c()
    
    for (i in 1:ncol(c2)) {
      
      I <- sort(c2[, i])
      a <- c()
      j <- 1
      for (j in 1:length(I)) {
        a[j] <- I[j] - I[j - 1]
      }
      
      b[[i]] <- a
      dat[[i]] <- round(data.frame(sort(abs(Targ[, i + 1])), sort(c2[, i]), b[[i]]), 3)
      names(dat[[i]]) <- c("sl.loadings","norm.loadings","diff.loadings")
      cuts_mix[i] <- dat[[i]][which(dat[[i]]$diff.loadings > mean(dat[[i]]$diff.loadings, na.rm = T)) - 1, ][1, ][[1]]
      
    }
    
    c2vec <- c()
    c2ref <- c()
    c2vec <- t(matrix(cuts_mix, ncol(Targ) - 1, nrow(Targ)))
    c2ref <- c2 - c2vec
    
    #FLOATING POINT PROBLEM
    Targ[, -1][c2ref >= 0.001] <- NA
    Targ[, -1][c2ref < 0.001] <- 0
    Targ[, 1] <- NA
    
    ### Identification conditions check
    j <- 1
    Targ2 <- Targ[, -1]
    Targ2[is.na(Targ2)] <- 1
    
    Targ2 <-Targ2*L[,-1] #New
    
    # prepare matrix for rank computation
    multiplier <- L
    a <- matrix(ncol = ncol(Targ2), nrow = 1)
    
    for (j in 1:ncol(Targ[, -1])) {
      
      if (mean(L[, j + 1]) < 0) {
        multiplier[, j + 1] <- L[, j + 1] * (-1)
      }
      
      # submatrix rank assessment  
      m <- Targ2[which(Targ2[, j] == 0), -j]
      if (length(m) <= 2) {
        a[1, j] <- qr(t(m))$rank
      } else {
        a[1, j] <- qr(m)$rank
      }
    }
    
    if (all(a == (ncol(Targ2) - 1))) {
      
      Targ2[Targ2 != 0] <- NA    #changed NEW REMOVED
      Targ[, 2:ncol(Targ)] <- Targ2
      
      return(list (Targ = Targ))
      
    } else {
      
      print(paste("Solution might not be identified"))
      cat ("\n")
      
      c <- which(a != (ncol(Targ2) - 1))
      h <- 1
      Targ2[Targ2 == 0] <- NA
      
      for (h in c) {
        Targ2[which.min(as.matrix(multiplier[, h + 1]) * Targ2[, h]),h] <- NA
        m <- Targ2[which(is.na(Targ2[, h])), -h]
        m[which(is.na(m))] <- 0
        if (length(m) <= 2) {
          a[1, h] <- qr(t(m))$rank
        } else {
          a[1, h] <- qr(m)$rank
        }
      }
      
      Targ2[is.na(Targ2)] <- 0
      Targ2[Targ2 != 0] <- NA ## CHANGED from == 1 to diff than zero
      Targ[, 2:ncol(Targ)] <- Targ2
      
      return(list (Targ = Targ))
    }
    
  } else {
    
    Targ[abs(Targ) >  cutpoint] <- NA
    Targ[abs(Targ) <= cutpoint] <- 0
  }
  
  Targ[,1] <- NA
  return(list(Targ = Targ))
}
######################################################################


######################################################################
target_convergence_check <- function (Targetprev, Targetnew) {
  
  # Convergence check for the SLi rotation
  #
  # Args:
  #
  #   Targetprev: Target matrix from previous iteration
  #   Targetnew : Target matrix from actual iteration
  #
  # Returns:
  #
  #   converg: If TRUE, convergence is achieved
  #            If FALSE, convergence is not achieved
  
  Targetprev[is.na(Targetprev)] <- 1
  Targetnew [is.na(Targetnew) ] <- 1
  
  if (sum((Targetnew - Targetprev)^2) == 0) {
    converg    <- TRUE
  }    else {
    converg    <- FALSE
  }
  return(converg)
}

######################################################################
# Schmid-Leiman algorithm function
######################################################################

SLi <- function (matrix = NULL,
                 specific_factors = NULL,
                 fm = "minres",
                 rotation = "geominQ",
                 cutpoint= "dif",
                 iterations = 20) {
  
  # Compute the iterative target factor rotation
  #
  # Args:
  #
  #   matrix: Correlation or Covariance matrix to be analyzed. Default is NULL.
  #   specific_factors: number of specific factors to be extracted. Default is NULL.
  #   fm: factor estimation method. Default is MINRES. Other alternatives can 
  #       be found in the fa() function documentation (psych) package.
  #   cutpoint: Value for cut-off point criterion (e.g., .20). Default is the Difference loadings criterion (SLiD).
  #   iterations: iterations number. Default is 20. If 0, Schmid-Leiman with target rotation (without iterations) is performed.
  #  
  # Returns:
  #
  #   Loadings: Rotated factor loading matrix
  #   Target: Last target applied
  #  
  # Error handling:  
  
  if (is.null(matrix)){
    stop("A correlation or covariance matrix must be specified")
  }
  
  if (is.null(specific_factors)){
    stop("A number of factors must be specified")
  }
  
  if (is.null(cutpoint)){
    stop("A cut-off point must be specified")
  }
  
  ## Step 1: First order factor analysis with Geomin oblique rotation
  
  lp  <- psych:::fa (r = matrix, 
                     nfactors = specific_factors, 
                     rotate = "none", 
                     fm = fm)$loadings
  
  #lp_rotated     <- perform_rotation (lp, method = "geominQ")
  
  if (rotation == "geominQ") {
    lp_rotated  <- GPArotation:::oblimin(lp,Tmat = Random.Start(ncol(lp)),maxit=5000)
  }  else if (rotation == "oblimin") {
    lp_rotated <- GPArotation:::oblimin(lp, Tmat =GPArotation::: Random.Start(ncol(lp)), maxit = 5000)
  }
  
  convergence_lp <- lp_rotated$convergence
  
  # Convergence check:    
  if (convergence_lp == FALSE) {
    stop("Convergence problems when estimating first order solution for SL solution")
    cat ("\n")
  }
  
  ## Step 2: Second order factor analysis
  lp1 <- psych:::fa (lp_rotated$Phi, nfactors=1, fm = fm)
  
  ## Step 3: Schmid-Leiman transformation
  lpSL1          <- lp_rotated$loadings %*% lp1$loadings
  psl1           <- matrix (0, dim(lp1$loadings), dim(lp1$loadings))
  diag(psl1)     <- sqrt(1- lp1$loadings^2)
  
  if (any(is.na(diag(psl1)))== TRUE) {
    stop(return(list(SL_loadings = NA,loadings = NA)))
  }
  
  lpsl2          <- lp_rotated$loadings %*% psl1
  SL_loadings    <- cbind (lpSL1,lpsl2)
  
  ## Step 4: Calculate an unrotated solution with specific + 1 factors  
  
  unrotated_l     <- psych:::fa (r = matrix, 
                                 nfactors = (specific_factors+1), 
                                 rotate ="none", 
                                 fm = fm)$loadings
  
  # Step 4: Schmid-Leiman iterated target rotation (SLi)
  
  targF           <- get_target_from (SL_loadings, cutpoint = cutpoint)
  SLt_result      <- perform_rotation (L = unrotated_l, 
                                       method = "targetT",
                                       targ = targF$Targ)
  
  # Convergence check
  if ( SLt_result$convergence == FALSE) {
    stop("Convergence problems when estimating SL target rotation")
    cat ("\n")
  }
  # SLt rotation factor loadings
  loadings_targ_prev <- SLt_result$loadings
  prev_target        <- targF
  
  # Step 5: Schmid-Leiman iterated target rotation (SLi)
  
  if (iterations == 0) {
    loadings_targ_new  <- SLt_result$loadings
    return(list(loadings = round(loadings_targ_new,3)))
    
  } else {
    for (it in 1:iterations) {
      
      new_target          <- get_target_from  (loadings_targ_prev, cutpoint = cutpoint)
      new_SLi_result      <- perform_rotation (loadings_targ_prev, 
                                               method = "targetT", 
                                               targ = new_target$Targ)
      
      if ( new_SLi_result$convergence == FALSE) {
        stop("Convergence problems when estimating SLi target rotation")
        cat ("\n")
      }
      
      loadings_targ_new   <- new_SLi_result$loadings
      
      # Check criteria for ending the iterative procedure
      converg <- target_convergence_check(prev_target$Targ, new_target$Targ)
      if (converg == TRUE) {
        cat("Convergence achieved in", it, "iterations \n")
        break()
      }   else {
        loadings_targ_prev   <- loadings_targ_new
        prev_target          <- new_target
      }
    }
    # Convergence check
    if (converg == FALSE) {
      cat("Convergence has not obtained for the target matrix. Please increase the number of iterations.")
      cat ("\n")
    }
  }
  
  # Return rotated factor matrix
  return(list(SL_loadings = SL_loadings,
              loadings = round(loadings_targ_new,3),
              target = new_target$Targ))
}  
######################################################################
