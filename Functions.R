#Notes:
#This has been updated to work with the most recent version of MGCV - 1.8 - 17
#Local version of R uses MGCV 1.8-7 at time of initial creation
#Key changes:
#pls and gdi have additonal parameters in their C functions
#SmoothCon uses a different methodology to scale the penalty matrix - changes reflexted in smoothCon.si


#01
gam.si <- function (
  formula, family = gaussian(), data = list(), weights = NULL, 
  subset = NULL, na.action, offset = NULL, method = "GCV.Cp", 
  optimizer = c("outer", "newton"), control = list(), scale = 0, 
  select = FALSE, knots = NULL, sp = NULL, min.sp = NULL, H = NULL, 
  gamma = 1, fit = TRUE, paraPen = NULL, G = NULL, in.out = NULL, 
  drop.unused.levels = TRUE,
  beta_star, index_data, splineDeg, splineRank, D_Design, D2_Design, S_uncons, ...) {
  XDbeta <- index_data%*%Dbeta(beta_star) 
  D2_beta <- D2beta(beta_star,dbeta = Dbeta(beta_star) )
  
  XD2beta <- list()
  for(i in 1:length(beta_star))
  {
    XD2beta[[i]] <- index_data %*% D2_beta[,i,]
  }
  
  control <- do.call("gam.control", control)
  if (is.null(G)) {
    gp <- interpret.gam(formula)
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    mf$formula <- gp$fake.formula
    mf$family <- mf$control <- mf$scale <- mf$knots <- mf$sp <- mf$min.sp <- mf$H <- mf$select <- mf$gamma <- mf$method <- mf$fit <- mf$paraPen <- mf$G <- mf$optimizer <- mf$in.out <- mf$... <- NULL
    mf$beta_star <- mf$index_data <- mf$splineDeg <- mf$splineRank <- mf$D_Design <- mf$D2_Design <- mf$S_uncons <- NULL
    mf$drop.unused.levels <- drop.unused.levels
    
    mf[[1]] <- as.name("model.frame")
    pmf <- mf
    mf <- eval(mf, parent.frame())
    if (nrow(mf) < 2) 
      stop("Not enough (non-NA) data to do anything meaningful")
    terms <- attr(mf, "terms")
    vars <- all.vars(gp$fake.formula[-2])
    inp <- parse(text = paste("list(", paste(vars, collapse = ","), 
                              ")"))
    if (!is.list(data) && !is.data.frame(data)) 
      data <- as.data.frame(data)
    dl <- eval(inp, data, parent.frame())
    names(dl) <- vars
    var.summary <- mgcv:::variable.summary(gp$pf, dl, nrow(mf))
    rm(dl)
    if (is.list(formula)) {
      environment(formula) <- environment(formula[[1]])
      pterms <- list()
      tlab <- rep("", 0)
      for (i in 1:length(formula)) {
        pmf$formula <- gp[[i]]$pf
        pterms[[i]] <- attr(eval(pmf, parent.frame()), 
                            "terms")
        tlabi <- attr(pterms[[i]], "term.labels")
        if (i > 1 && length(tlabi) > 0) 
          tlabi <- paste(tlabi, i - 1, sep = ".")
        tlab <- c(tlab, tlabi)
      }
      attr(pterms, "term.labels") <- tlab
    }
    else {
      pmf$formula <- gp$pf
      pmf <- eval(pmf, parent.frame())
      pterms <- attr(pmf, "terms")
    }
    if (is.character(family)) 
      family <- eval(parse(text = family))
    if (is.function(family)) 
      family <- family()
    if (is.null(family$family)) 
      stop("family not recognized")
    if (family$family[1] == "gaussian" && family$link == 
        "identity") 
      am <- TRUE
    else am <- FALSE
    if (!control$keepData) 
      rm(data)
    drop.intercept <- if (is.null(family$drop.intercept) || 
                          !family$drop.intercept) 
      FALSE
    else TRUE
    if (is.list(formula)) G <- mgcv:::gam.setup.list(formula = gp, pterms = pterms, 
                                                     data = mf, knots = knots, sp = sp, min.sp = min.sp, 
                                                     H = H, absorb.cons = TRUE, sparse.cons = 0, select = select, 
                                                     idLinksBases = control$idLinksBases, scale.penalty = control$scalePenalty, 
                                                     paraPen = paraPen, drop.intercept = drop.intercept)
    else G <- gam.setup.si(formula = gp, pterms = pterms, 
                           data = mf, knots = knots, sp = sp, min.sp = min.sp, 
                           H = H, absorb.cons = TRUE, sparse.cons = 0, select = select, 
                           idLinksBases = control$idLinksBases, scale.penalty = control$scalePenalty, 
                           paraPen = paraPen, drop.intercept = drop.intercept)
    
    G$var.summary <- var.summary
    G$family <- family
    if ((is.list(formula) && (is.null(family$nlp) || family$nlp != 
                              gp$nlp)) || (!is.list(formula) && !is.null(family$npl) && 
                                           (family$npl > 1))) 
      stop("incorrect number of linear predictors for family")
    if (ncol(G$X) > nrow(G$X)) 
      stop("Model has more coefficients than data")
    G$terms <- terms
    G$mf <- mf
    G$cl <- cl
    G$am <- am
    if (is.null(G$offset)) 
      G$offset <- rep(0, G$n)
    G$min.edf <- G$nsdf
    if (G$m) 
      for (i in 1:G$m) G$min.edf <- G$min.edf + G$smooth[[i]]$null.space.dim
    G$formula <- formula
    G$pred.formula <- gp$pred.formula
    environment(G$formula) <- environment(formula)
  }
  if (!fit) 
    return(G)
  G$conv.tol <- control$mgcv.tol
  G$max.half <- control$mgcv.half
  #Set up first and second derivatives of X
  object <- estimate.gam.si(beta_star, index_data, splineDeg, splineRank, D_Design, D2_Design, S_uncons,
                            G, method, optimizer, control, in.out, 
                            scale, gamma, ...)
  if (!is.null(G$L)) {
    object$full.sp <- as.numeric(exp(G$L %*% log(object$sp) + 
                                       G$lsp0))
    names(object$full.sp) <- names(G$lsp0)
  }
  names(object$sp) <- names(G$sp)
  object$paraPen <- G$pP
  object$formula <- G$formula
  if (is.list(object$formula)) 
    attr(object$formula, "lpi") <- attr(G$X, "lpi")
  object$var.summary <- G$var.summary
  object$cmX <- G$cmX
  object$model <- G$mf
  object$na.action <- attr(G$mf, "na.action")
  object$control <- control
  object$terms <- G$terms
  object$pred.formula <- G$pred.formula
  attr(object$pred.formula, "full") <- reformulate(all.vars(object$terms))
  object$pterms <- G$pterms
  object$assign <- G$assign
  object$contrasts <- G$contrasts
  object$xlevels <- G$xlevels
  object$offset <- G$offset
  object$ZSZ <- G$S
  if (!is.null(G$Xcentre)) 
    object$Xcentre <- G$Xcentre
  if (control$keepData) 
    object$data <- data
  object$df.residual <- nrow(G$X) - sum(object$edf)
  object$min.edf <- G$min.edf
  if (G$am && !(method %in% c("REML", "ML", "P-ML", "P-REML"))) 
    object$optimizer <- "magic"
  else object$optimizer <- optimizer
  object$call <- G$cl
  class(object) <- c("gam", "glm", "lm")
  if (is.null(object$deviance)) 
    object$deviance <- sum(residuals(object, "deviance")^2)
  environment(object$formula) <- environment(object$pred.formula) <- environment(object$terms) <- environment(object$pterms) <- .GlobalEnv
  if (!is.null(object$model)) 
    environment(attr(object$model, "terms")) <- .GlobalEnv
  if (!is.null(attr(object$pred.formula, "full"))) 
    environment(attr(object$pred.formula, "full")) <- .GlobalEnv
  object
}

#01a
gam.setup.si <- function (formula, pterms, data = stop("No data supplied to gam.setup"), 
                          knots = NULL, sp = NULL, min.sp = NULL, H = NULL, absorb.cons = TRUE, 
                          sparse.cons = 0, select = FALSE, idLinksBases = TRUE, scale.penalty = TRUE, 
                          paraPen = NULL, gamm.call = FALSE, drop.intercept = FALSE, 
                          diagonal.penalty = FALSE) {
  if (inherits(formula, "split.gam.formula")) 
    split <- formula
  else if (inherits(formula, "formula")) 
    split <- interpret.gam(formula)
  else stop("First argument is no sort of formula!")
  if (length(split$smooth.spec) == 0) {
    if (split$pfok == 0) 
      stop("You've got no model....")
    m <- 0
  }
  else m <- length(split$smooth.spec)
  G <- list(m = m, min.sp = min.sp, H = H, pearson.extra = 0, 
            dev.extra = 0, n.true = -1, pterms = pterms)
  if (is.null(attr(data, "terms"))) 
    mf <- model.frame(split$pf, data, drop.unused.levels = FALSE)
  else mf <- data
  G$intercept <- attr(attr(mf, "terms"), "intercept") > 0
  G$offset <- model.offset(mf)
  if (!is.null(G$offset)) 
    G$offset <- as.numeric(G$offset)
  if (drop.intercept) 
    attr(pterms, "intercept") <- 1
  X <- model.matrix(pterms, mf)
  if (drop.intercept) {
    xat <- attributes(X)
    ind <- xat$assign > 0
    X <- X[, xat$assign > 0, drop = FALSE]
    xat$assign <- xat$assign[ind]
    xat$dimnames[[2]] <- xat$dimnames[[2]][ind]
    xat$dim[2] <- xat$dim[2] - 1
    attributes(X) <- xat
    G$intercept <- FALSE
  }
  rownames(X) <- NULL
  G$nsdf <- ncol(X)
  G$contrasts <- attr(X, "contrasts")
  G$xlevels <- .getXlevels(pterms, mf)
  G$assign <- attr(X, "assign")
  PP <- mgcv:::parametricPenalty(pterms, G$assign, paraPen, sp)
  if (!is.null(PP)) {
    ind <- 1:length(PP$sp)
    if (!is.null(sp)) 
      sp <- sp[-ind]
    if (!is.null(min.sp)) {
      PP$min.sp <- min.sp[ind]
      min.sp <- min.sp[-ind]
    }
  }
  G$smooth <- list()
  G$S <- list()
  G$S_uncons <- list() 
  if (gamm.call) {
    if (m > 0) 
      for (i in 1:m) attr(split$smooth.spec[[i]], "gamm") <- TRUE
  }
  if (m > 0 && idLinksBases) {
    id.list <- list()
    for (i in 1:m) if (!is.null(split$smooth.spec[[i]]$id)) {
      id <- as.character(split$smooth.spec[[i]]$id)
      if (length(id.list) && id %in% names(id.list)) {
        ni <- length(id.list[[id]]$sm.i)
        id.list[[id]]$sm.i[ni + 1] <- i
        base.i <- id.list[[id]]$sm.i[1]
        split$smooth.spec[[i]] <- clone.smooth.spec(split$smooth.spec[[base.i]], 
                                                    split$smooth.spec[[i]])
        temp.term <- split$smooth.spec[[i]]$term
        for (j in 1:length(temp.term)) id.list[[id]]$data[[j]] <- cbind(id.list[[id]]$data[[j]], 
                                                                        get.var(temp.term[j], data, vecMat = FALSE))
      }
      else {
        id.list[[id]] <- list(sm.i = i)
        id.list[[id]]$data <- list()
        term <- split$smooth.spec[[i]]$term
        for (j in 1:length(term)) id.list[[id]]$data[[j]] <- get.var(term[j], 
                                                                     data, vecMat = FALSE)
      }
    }
  }
  G$off <- array(0, 0)
  first.para <- G$nsdf + 1
  sm <- list()
  G$X_u <- list()
  newm <- 0
  if (m > 0) 
    for (i in 1:m) {
      #We call smooth construct for each smooth S
      id <- split$smooth.spec[[i]]$id
      if (is.null(id) || !idLinksBases) { #This is what we do
        sml <- smoothCon.si(split$smooth.spec[[i]], data, 
                            knots, absorb.cons, scale.penalty = scale.penalty, 
                            null.space.penalty = select, sparse.cons = sparse.cons, 
                            diagonal.penalty = diagonal.penalty)
      }
      else {
        names(id.list[[id]]$data) <- split$smooth.spec[[i]]$term
        sml <- smoothCon.si(split$smooth.spec[[i]], id.list[[id]]$data, 
                            knots, absorb.cons, n = nrow(data), dataX = data, 
                            scale.penalty = scale.penalty, null.space.penalty = select, 
                            sparse.cons = sparse.cons, diagonal.penalty = diagonal.penalty)
      }
      for (j in 1:length(sml)) {
        newm <- newm + 1
        sm[[newm]] <- sml[[j]]
      }
      G$X_u[[i]] <- sml[[1]]$X_uncons
    }
  
  G$m <- m <- newm
  if (m > 0) 
    sm <- gam.side(sm, X, tol = .Machine$double.eps^0.5)
  idx <- list()
  L <- matrix(0, 0, 0)
  lsp.names <- sp.names <- rep("", 0)
  if (m > 0) 
    for (i in 1:m) {
      id <- sm[[i]]$id
      length.S <- length(sm[[i]]$S)
      if (is.null(sm[[i]]$L)) 
        Li <- diag(length.S)
      else Li <- sm[[i]]$L
      if (length.S > 0) {
        if (length.S == 1) 
          spn <- sm[[i]]$label
        else {
          Sname <- names(sm[[i]]$S)
          if (is.null(Sname)) 
            spn <- paste(sm[[i]]$label, 1:length.S, sep = "")
          else spn <- paste(sm[[i]]$label, Sname, sep = "")
        }
      }
      if (is.null(id) || is.null(idx[[id]])) {
        if (!is.null(id)) {
          idx[[id]]$c <- ncol(L) + 1
          idx[[id]]$nc <- ncol(Li)
        }
        L <- rbind(cbind(L, matrix(0, nrow(L), ncol(Li))), 
                   cbind(matrix(0, nrow(Li), ncol(L)), Li))
        if (length.S > 0) {
          sp.names <- c(sp.names, spn)
          lsp.names <- c(lsp.names, spn)
        }
      }
      else {
        L0 <- matrix(0, nrow(Li), ncol(L))
        if (ncol(Li) > idx[[id]]$nc) {
          stop("Later terms sharing an `id' can not have more smoothing parameters than the first such term")
        }
        L0[, idx[[id]]$c:(idx[[id]]$c + ncol(Li) - 1)] <- Li
        L <- rbind(L, L0)
        if (length.S > 0) {
          lsp.names <- c(lsp.names, spn)
        }
      }
    }
  Xp <- NULL
  if (m > 0) 
    for (i in 1:m) { #Creates indexes for first and last parameter per smooth?
      n.para <- ncol(sm[[i]]$X)
      sm[[i]]$first.para <- first.para
      first.para <- first.para + n.para
      sm[[i]]$last.para <- first.para - 1
      Xoff <- attr(sm[[i]]$X, "offset")
      if (!is.null(Xoff)) {
        if (is.null(G$offset)) 
          G$offset <- Xoff
        else G$offset <- G$offset + Xoff
      }
      if (is.null(sm[[i]]$Xp)) {
        if (!is.null(Xp)) 
          Xp <- cbind2(Xp, sm[[i]]$X)
      }
      else {
        if (is.null(Xp)) 
          Xp <- X
        Xp <- cbind2(Xp, sm[[i]]$Xp)
        sm[[i]]$Xp <- NULL
      }
      X <- cbind2(X, sm[[i]]$X)
      sm[[i]]$X <- NULL
      G$smooth[[i]] <- sm[[i]]
    }
  if (is.null(Xp)) {
    G$cmX <- colMeans(X)
  }
  else {
    G$cmX <- colMeans(Xp)
    qrx <- qr(Xp, LAPACK = TRUE)
    R <- qr.R(qrx)
    p <- ncol(R)
    rank <- Rrank(R)
    QtX <- qr.qty(qrx, X)[1:rank, ]
    if (rank < p) {
      R <- R[1:rank, ]
      qrr <- qr(t(R), tol = 0)
      R <- qr.R(qrr)
      G$P <- forwardsolve(t(R), QtX)
    }
    else {
      G$P <- backsolve(R, QtX)
    }
    if (rank < p) {
      G$P <- qr.qy(qrr, rbind(G$P, matrix(0, p - rank, 
                                          p)))
    }
    G$P[qrx$pivot, ] <- G$P
  }
  if (G$nsdf > 0) 
    G$cmX[-(1:G$nsdf)] <- 0
  else G$cmX <- G$cmX * 0
  G$X <- X
  rm(X)
  n.p <- ncol(G$X)
  if (!is.null(sp)) {
    if (length(sp) != ncol(L)) {
      warning("Supplied smoothing parameter vector is too short - ignored.")
    }
    if (sum(is.na(sp))) {
      warning("NA's in supplied smoothing parameter vector - ignoring.")
    }
    G$sp <- sp
  }
  else {
    G$sp <- rep(-1, ncol(L))
  }
  names(G$sp) <- sp.names
  k <- 1
  if (m > 0) 
    for (i in 1:m) {
      id <- sm[[i]]$id
      if (is.null(sm[[i]]$L)) 
        Li <- diag(length(sm[[i]]$S))
      else Li <- sm[[i]]$L
      if (is.null(id)) {
        spi <- sm[[i]]$sp
        if (!is.null(spi)) {
          if (length(spi) != ncol(Li)) 
            stop("incorrect number of smoothing parameters supplied for a smooth term")
          G$sp[k:(k + ncol(Li) - 1)] <- spi
        }
        k <- k + ncol(Li)
      }
      else {
        spi <- sm[[i]]$sp
        if (is.null(idx[[id]]$sp.done)) {
          if (!is.null(spi)) {
            if (length(spi) != ncol(Li)) 
              stop("incorrect number of smoothing parameters supplied for a smooth term")
            G$sp[idx[[id]]$c:(idx[[id]]$c + idx[[id]]$nc - 
                                1)] <- spi
          }
          idx[[id]]$sp.done <- TRUE
          k <- k + idx[[id]]$nc
        }
      }
    }
  k <- 1
  if (length(idx)) 
    for (i in 1:length(idx)) idx[[i]]$sp.done <- FALSE
  if (m > 0) 
    for (i in 1:m) {
      id <- sm[[i]]$id
      if (!is.null(id)) {
        if (idx[[id]]$nc > 0) {
          G$smooth[[i]]$sp <- G$sp[idx[[id]]$c:(idx[[id]]$c + 
                                                  idx[[id]]$nc - 1)]
        }
        if (!idx[[id]]$sp.done) {
          idx[[id]]$sp.done <- TRUE
          k <- k + idx[[id]]$nc
        }
      }
      else {
        if (is.null(sm[[i]]$L)) 
          nc <- length(sm[[i]]$S)
        else nc <- ncol(sm[[i]]$L)
        if (nc > 0) 
          G$smooth[[i]]$sp <- G$sp[k:(k + nc - 1)]
        k <- k + nc
      }
    }
  if (!is.null(min.sp)) {
    if (length(min.sp) != nrow(L)) 
      stop("length of min.sp is wrong.")
    if (sum(is.na(min.sp))) 
      stop("NA's in min.sp.")
    if (sum(min.sp < 0)) 
      stop("elements of min.sp must be non negative.")
  }
  k.sp <- 0
  G$rank <- array(0, 0)
  if (m > 0) 
    for (i in 1:m) {
      sm <- G$smooth[[i]]
      if (length(sm$S) > 0) 
        for (j in 1:length(sm$S)) {
          k.sp <- k.sp + 1
          G$off[k.sp] <- sm$first.para
          G$S[[k.sp]] <- sm$S[[j]]
          G$S_uncons[[k.sp]] <- sm$S_uncons[[j]]
          #G$X_uncons[[k.sp]] <- sm[[j]]$X_uncons
          G$rank[k.sp] <- sm$rank[j]
          if (!is.null(min.sp)) {
            if (is.null(H)) 
              H <- matrix(0, n.p, n.p)
            H[sm$first.para:sm$last.para, sm$first.para:sm$last.para] <- H[sm$first.para:sm$last.para, 
                                                                           sm$first.para:sm$last.para] + min.sp[k.sp] * 
              sm$S[[j]]
          }
        }
    }
  
  if (!is.null(PP)) { #Didn't get here    
    L <- rbind(cbind(L, matrix(0, nrow(L), ncol(PP$L))), 
               cbind(matrix(0, nrow(PP$L), ncol(L)), PP$L))
    G$off <- c(PP$off, G$off)
    G$S <- c(PP$S, G$S)
    G$rank <- c(PP$rank, G$rank)
    G$sp <- c(PP$sp, G$sp)
    lsp.names <- c(PP$full.sp.names, lsp.names)
    G$n.paraPen <- length(PP$off)
    if (!is.null(PP$min.sp)) {
      if (is.null(H)) 
        H <- matrix(0, n.p, n.p)
      for (i in 1:length(PP$S)) {
        ind <- PP$off[i]:(PP$off[i] + ncol(PP$S[[i]]) - 
                            1)
        H[ind, ind] <- H[ind, ind] + PP$min.sp[i] * PP$S[[i]]
      }
    }
  }
  else G$n.paraPen <- 0  
  fix.ind <- G$sp >= 0
  if (sum(fix.ind)) { #Didn't get here
    lsp0 <- G$sp[fix.ind]
    ind <- lsp0 == 0
    ef0 <- indi <- (1:length(ind))[ind]
    if (length(indi) > 0) 
      for (i in 1:length(indi)) {
        ii <- G$off[i]:(G$off[i] + ncol(G$S[[i]]) - 1)
        ef0[i] <- norm(G$X[, ii], type = "F")^2/norm(G$S[[i]], 
                                                     type = "F") * .Machine$double.eps * 0.1
      }
    lsp0[!ind] <- log(lsp0[!ind])
    lsp0[ind] <- log(ef0)
    lsp0 <- as.numeric(L[, fix.ind, drop = FALSE] %*% lsp0)
    L <- L[, !fix.ind, drop = FALSE]
    G$sp <- G$sp[!fix.ind]
  }
  else {
    lsp0 <- rep(0, nrow(L))
  }
  G$H <- H
  if (ncol(L) == nrow(L) && !sum(L != diag(ncol(L)))) 
    L <- NULL
  G$L <- L
  G$lsp0 <- lsp0
  names(G$lsp0) <- lsp.names
  if (absorb.cons == FALSE) { #Didn't get here
    G$C <- matrix(0, 0, n.p)
    if (m > 0) {
      for (i in 1:m) {
        if (is.null(G$smooth[[i]]$C)) 
          n.con <- 0
        else n.con <- nrow(G$smooth[[i]]$C)
        C <- matrix(0, n.con, n.p)
        C[, G$smooth[[i]]$first.para:G$smooth[[i]]$last.para] <- G$smooth[[i]]$C
        G$C <- rbind(G$C, C)
        G$smooth[[i]]$C <- NULL
      }
      rm(C)
    }
  }
  G$y <- data[[split$response]]
  G$n <- nrow(data)
  if (is.null(data$"(weights)")) 
    G$w <- rep(1, G$n)
  else G$w <- data$"(weights)"
  if (G$nsdf > 0) 
    term.names <- colnames(G$X)[1:G$nsdf]
  else term.names <- array("", 0)
  n.smooth <- length(G$smooth)
  if (n.smooth) 
    for (i in 1:n.smooth) {
      k <- 1
      jj <- G$smooth[[i]]$first.para:G$smooth[[i]]$last.para
      if (G$smooth[[i]]$df > 0) 
        for (j in jj) {
          term.names[j] <- paste(G$smooth[[i]]$label, 
                                 ".", as.character(k), sep = "")
          k <- k + 1
        }
    }
  G$term.names <- term.names
  G$pP <- PP
  G$X_uncons <- sm$X_uncons
  G
}

#01b
smoothCon.si <- function (object, data, knots, absorb.cons = FALSE, scale.penalty = TRUE, 
                          n = nrow(data), dataX = NULL, null.space.penalty = FALSE, 
                          sparse.cons = 0, diagonal.penalty = FALSE) {
  sm <- mgcv:::smooth.construct3(object, data, knots)
  X_uncons <- sm$X
  if (!is.null(attr(sm, "qrc"))) 
    warning("smooth objects should not have a qrc attribute.")
  if (is.null(sm$plot.me)) 
    sm$plot.me <- TRUE
  if (is.null(sm$side.constrain)) 
    sm$side.constrain <- TRUE
  drop <- -1
  if (is.null(sm$C)) { # This is where we calculate X?
    if (sparse.cons <= 0) {
      sm$C <- matrix(colMeans(sm$X), 1, ncol(sm$X))
      if (sparse.cons == -1) {
        vcol <- apply(sm$X, 2, var)
        drop <- min((1:length(vcol))[vcol == min(vcol)])
      }
    }
    else if (sparse.cons > 0) {
      if (sum(sm$X == 0) > 0.1 * sum(sm$X != 0)) {
        if (sparse.cons == 1) {
          xsd <- apply(sm$X, 2, FUN = sd)
          if (sum(xsd == 0)) 
            sm$C <- ((1:length(xsd))[xsd == 0])[1]
          else {
            xz <- apply(sm$X, 2, FUN = function(x) {
              sum(x == 0)
            })
            sm$C <- ((1:length(xz))[xz == min(xz)])[1]
          }
        }
        else if (sparse.cons == 2) {
          sm$C = -1
        }
        else {
          stop("unimplemented sparse constraint type requested")
        }
      }
      else {
        sm$C <- matrix(colSums(sm$X), 1, ncol(sm$X))
      }
    }
    else {
      sm$C <- matrix(colSums(sm$X), 1, ncol(sm$X))
    }
    alwaysCon <- FALSE
  }
  else {
    if (is.null(attr(sm$C, "always.apply"))) 
      alwaysCon <- FALSE
    else alwaysCon <- TRUE
  }
  if (is.null(sm$df)) 
    sm$df <- sm$bs.dim
  if (!is.null(object$fixed) && object$fixed) {
    sm$S <- NULL
  }
  sm$S.scale <- rep(1, length(sm$S))
  sm$S_uncons <- list() #This code has been updated so we output the unconstrained, scaled penalty matrices
  
  #This is the method for the version 1.8 - 7
  # if (scale.penalty && length(sm$S) > 0 && is.null(sm$no.rescale)) {
  #   maXX <- mean(abs(t(sm$X) %*% sm$X))
  #   for (i in 1:length(sm$S)) {
  #     maS <- mean(abs(sm$S[[i]]))/maXX
  #     sm$S[[i]] <- sm$S[[i]]/maS
  #     sm$S.scale[i] <- maS
  #     sm$S_uncons[[i]] <- sm$S[[i]]     
  #   }
  # }
  
  #This is the method for version 1.8-17
  if (scale.penalty && length(sm$S) > 0 && is.null(sm$no.rescale)) {
    maXX <- norm(sm$X, type = "I")^2
    for (i in 1:length(sm$S)) {
      maS <- norm(sm$S[[i]])/maXX
      sm$S[[i]] <- sm$S[[i]]/maS
      sm$S.scale[i] <- maS
      sm$S_uncons[[i]] <- sm$S[[i]]     
    }
  }
  
  if (!is.null(dataX)) {
    er <- Predict.matrix3(sm, dataX)
    sm$X <- er$X
    sm$ind <- er$ind
    rm(er)
  }
  if ((is.null(sm$ind) && nrow(sm$X) != n) || (!is.null(sm$ind) && 
                                               length(sm$ind) != n)) {
    matrixArg <- TRUE
    if (is.null(sm$ind)) 
      q <- nrow(sm$X)/n
    else q <- length(sm$ind)/n
    if (!is.null(sm$by.done)) 
      warning("handling `by' variables in smooth constructors may not work with the summation convention ")
  }
  else {
    matrixArg <- FALSE
    if (!is.null(sm$ind)) {
      offs <- attr(sm$X, "offset")
      sm$X <- sm$X[sm$ind, ]
      if (!is.null(offs)) 
        attr(sm$X, "offset") <- offs[sm$ind]
    }
  }
  offs <- NULL
  if (matrixArg || (object$by != "NA" && is.null(sm$by.done))) { #We don't get into this big if statement
    drop <- -1
    if (is.null(dataX)) 
      by <- get.var(object$by, data)
    else by <- get.var(object$by, dataX)
    if (matrixArg && is.null(by)) {
      if (is.null(sm$ind)) 
        by <- rep(1, nrow(sm$X))
      else by <- rep(1, length(sm$ind))
    }
    if (is.null(by)) 
      stop("Can't find by variable")
    offs <- attr(sm$X, "offset")
    if (is.factor(by)) {
      if (matrixArg) 
        stop("factor `by' variables can not be used with matrix arguments.")
      sml <- list()
      lev <- levels(by)
      if (is.ordered(by) && length(lev) > 1) 
        lev <- lev[-1]
      for (j in 1:length(lev)) {
        sml[[j]] <- sm
        by.dum <- as.numeric(lev[j] == by)
        sml[[j]]$X <- by.dum * sm$X
        sml[[j]]$by.level <- lev[j]
        sml[[j]]$label <- paste(sm$label, ":", object$by, 
                                lev[j], sep = "")
        if (!is.null(offs)) {
          attr(sml[[j]]$X, "offset") <- offs * by.dum
        }
      }
    }
    else {
      sml <- list(sm)
      if ((is.null(sm$ind) && length(by) != nrow(sm$X)) || 
          (!is.null(sm$ind) && length(by) != length(sm$ind))) 
        stop("`by' variable must be same dimension as smooth arguments")
      if (matrixArg) {
        if (is.null(sm$ind)) {
          sml[[1]]$X <- as.numeric(by) * sm$X
          ind <- 1:n
          X <- sml[[1]]$X[ind, ]
          for (i in 2:q) {
            ind <- ind + n
            X <- X + sml[[1]]$X[ind, ]
          }
          sml[[1]]$X <- X
          if (!is.null(offs)) {
            offs <- attr(sm$X, "offset") * as.numeric(by)
            ind <- 1:n
            offX <- offs[ind, ]
            for (i in 2:q) {
              ind <- ind + n
              offX <- offX + offs[ind, ]
            }
            attr(sml[[1]]$X, "offset") <- offX
          }
        }
        else {
          ind <- 0:(q - 1) * n
          offs <- attr(sm$X, "offset")
          if (!is.null(offs)) 
            offX <- rep(0, n)
          else offX <- NULL
          sml[[1]]$X <- matrix(0, n, ncol(sm$X))
          for (i in 1:n) {
            ind <- ind + 1
            sml[[1]]$X[i, ] <- colSums(by[ind] * sm$X[sm$ind[ind], 
                                                      ])
            if (!is.null(offs)) {
              offX[i] <- sum(offs[sm$ind[ind]] * by[ind])
            }
          }
          attr(sml[[1]]$X, "offset") <- offX
        }
      }
      else {
        sml[[1]]$X <- as.numeric(by) * sm$X
        if (!is.null(offs)) 
          attr(sml[[1]]$X, "offset") <- offs * as.numeric(by)
      }
      if (object$by == "NA") 
        sml[[1]]$label <- sm$label
      else sml[[1]]$label <- paste(sm$label, ":", object$by, 
                                   sep = "")
      if (!alwaysCon) {
        if (matrixArg) {
          L1 <- as.numeric(matrix(by, n, q) %*% rep(1, 
                                                    q))
          if (sd(L1) > mean(L1) * .Machine$double.eps * 
              1000) {
            sml[[1]]$C <- sm$C <- matrix(0, 0, 1)
            if (!is.null(sm$Cp)) 
              sml[[1]]$Cp <- sm$Cp <- NULL
          }
          else sml[[1]]$meanL1 <- mean(L1)
        }
        else {
          if (sd(by) > mean(by) * .Machine$double.eps * 
              1000) {
            sml[[1]]$C <- sm$C <- matrix(0, 0, 1)
            if (!is.null(sm$Cp)) 
              sml[[1]]$Cp <- sm$Cp <- NULL
          }
        }
      }
    }
  }
  else { #This is where we first set up sml
    sml <- list(sm)
  }
  if (absorb.cons) { #We absorb constraints so this is the if statement we use - Seem to reparam here
    sml[[1]]$X_uncons <- X_uncons #This might need to be readdressed when we have multiple smooths
    k <- ncol(sm$X)
    if (!is.null(sm$Cp) && is.matrix(sm$Cp)) {
      pj <- nrow(sm$Cp)
      qrcp <- qr(t(sm$Cp))
      for (i in 1:length(sml)) {
        sml[[i]]$Xp <- t(qr.qty(qrcp, t(sml[[i]]$X))[(pj + 
                                                        1):k, ])
        sml[[i]]$Cp <- NULL
        if (length(sml[[i]]$S)) {
          sml[[i]]$Sp <- sml[[i]]$S
          for (l in 1:length(sml[[i]]$S)) {
            ZSZ <- qr.qty(qrcp, sml[[i]]$S[[l]])[(pj + 
                                                    1):k, ]
            sml[[i]]$Sp[[l]] <- t(qr.qty(qrcp, t(ZSZ))[(pj + 
                                                          1):k, ])
          }
        }
      }
    }
    else qrcp <- NULL
    if (is.matrix(sm$C)) {
      j <- nrow(sm$C)
      if (j > 0) {
        indi <- (1:ncol(sm$C))[colSums(sm$C) != 0]
        nx <- length(indi)
        if (nx < ncol(sm$C)) {
          nc <- j
          nz <- nx - nc
          qrc <- qr(t(sm$C[, indi, drop = FALSE]))
          for (i in 1:length(sml)) { 
            if (length(sm$S) > 0) 
              for (l in 1:length(sm$S)) {
                ZSZ <- sml[[i]]$S[[l]]
                ZSZ[indi[1:nz], ] <- qr.qty(qrc, sml[[i]]$S[[l]][indi, 
                                                                 , drop = FALSE])[(nc + 1):nx, ]
                ZSZ <- ZSZ[-indi[(nz + 1):nx], ]
                ZSZ[, indi[1:nz]] <- t(qr.qty(qrc, t(ZSZ[, 
                                                         indi, drop = FALSE]))[(nc + 1):nx, 
                                                                               ])
                sml[[i]]$S[[l]] <- ZSZ[, -indi[(nz + 
                                                  1):nx], drop = FALSE]
              }
            sml[[i]]$X[, indi[1:nz]] <- t(qr.qty(qrc, 
                                                 t(sml[[i]]$X[, indi, drop = FALSE]))[(nc + 
                                                                                         1):nx, ])
            sml[[i]]$X <- sml[[i]]$X[, -indi[(nz + 1):nx]]
            attr(sml[[i]], "qrc") <- qrc
            attr(sml[[i]], "nCons") <- j
            attr(sml[[i]], "indi") <- indi
            sml[[i]]$C <- NULL
            sml[[i]]$rank <- pmin(sm$rank, k - j)
            sml[[i]]$df <- sml[[i]]$df - j
            sml[[i]]$null.space.dim <- max(0, sml[[i]]$null.space.dim - 
                                             j)
          }
        }
        else { 
          if (drop > 0) {
            qrc <- c(drop, as.numeric(sm$C)[-drop])
            class(qrc) <- "sweepDrop"
            for (i in 1:length(sml)) {
              sml[[i]]$X <- sml[[i]]$X[, -drop] - matrix(qrc[-1], 
                                                         nrow(sml[[i]]$X), ncol(sml[[i]]$X) - 
                                                           1, byrow = TRUE)
              
              if (length(sm$S) > 0) 
                for (l in 1:length(sm$S)) {
                  sml[[i]]$S[[l]] <- sml[[i]]$S[[l]][-drop, 
                                                     -drop]
                }
            }
          }
          else {
            qrc <- qr(t(sm$C))
            for (i in 1:length(sml)) {
              if (length(sm$S) > 0) 
                for (l in 1:length(sm$S)) {
                  ZSZ <- qr.qty(qrc, sm$S[[l]])[(j + 
                                                   1):k, ]
                  sml[[i]]$S[[l]] <- t(qr.qty(qrc, t(ZSZ))[(j + 
                                                              1):k, ])
                }
              sml[[i]]$X <- t(qr.qty(qrc, t(sml[[i]]$X))[(j + 
                                                            1):k, ])
            }
          }
          for (i in 1:length(sml)) {
            attr(sml[[i]], "qrc") <- qrc
            attr(sml[[i]], "nCons") <- j
            sml[[i]]$C <- NULL
            sml[[i]]$rank <- pmin(sm$rank, k - j)
            sml[[i]]$df <- sml[[i]]$df - j
            sml[[i]]$null.space.dim <- max(0, sml[[i]]$null.space.dim - 
                                             j)
          }
        }
      }
      else {
        for (i in 1:length(sml)) {
          attr(sml[[i]], "qrc") <- "no constraints"
          attr(sml[[i]], "nCons") <- 0
        }
      }
    }
    else if (sm$C > 0) {
      for (i in 1:length(sml)) {
        if (length(sm$S) > 0) 
          for (l in 1:length(sm$S)) {
            sml[[i]]$S[[l]] <- sml[[i]]$S[[l]][-sm$C, 
                                               -sm$C]
          }
        sml[[i]]$X <- sml[[i]]$X[, -sm$C]
        attr(sml[[i]], "qrc") <- sm$C
        attr(sml[[i]], "nCons") <- 1
        sml[[i]]$C <- NULL
        sml[[i]]$rank <- pmin(sm$rank, k - 1)
        sml[[i]]$df <- sml[[i]]$df - 1
        sml[[i]]$null.space.dim <- max(sml[[i]]$null.space.dim - 
                                         1, 0)
      }
    }
    else if (sm$C < 0) {
      for (i in 1:length(sml)) {
        if (length(sm$S) > 0) 
          for (l in 1:length(sm$S)) {
            sml[[i]]$S[[l]] <- diff(t(diff(sml[[i]]$S[[l]])))
          }
        sml[[i]]$X <- t(diff(t(sml[[i]]$X)))
        attr(sml[[i]], "qrc") <- sm$C
        attr(sml[[i]], "nCons") <- 1
        sml[[i]]$C <- NULL
        sml[[i]]$rank <- pmin(sm$rank, k - 1)
        sml[[i]]$df <- sml[[i]]$df - 1
        sml[[i]]$null.space.dim <- max(sml[[i]]$null.space.dim - 
                                         1, 0)
      }
    }
    if (!is.null(qrcp)) {
      for (i in 1:length(sml)) {
        attr(sml[[i]], "qrc") <- qrcp
        if (pj != attr(sml[[i]], "nCons")) 
          stop("Number of prediction and fit constraints must match")
        attr(sml[[i]], "indi") <- NULL
      }
    }
  }
  else for (i in 1:length(sml)) attr(sml[[i]], "qrc") <- NULL
  if (diagonal.penalty && length(sml[[1]]$S) == 1) {
    S11 <- sml[[1]]$S[[1]][1, 1]
    rank <- sml[[1]]$rank
    p <- ncol(sml[[1]]$X)
    if (is.null(rank) || max(abs(sml[[1]]$S[[1]] - diag(c(rep(S11, 
                                                              rank), rep(0, p - rank))))) > abs(S11) * .Machine$double.eps^0.8) {
      np <- nat.param(sml[[1]]$X, sml[[1]]$S[[1]], rank = sml[[1]]$rank, 
                      type = 2, unit.fnorm = FALSE)
      sml[[1]]$X <- np$X
      sml[[1]]$S[[1]] <- diag(p)
      diag(sml[[1]]$S[[1]]) <- c(np$D, rep(0, p - np$rank))
      sml[[1]]$diagRP <- np$P
      if (length(sml) > 1) 
        for (i in 2:length(sml)) {
          sml[[i]]$X <- sml[[i]]$X %*% np$P
          sml[[i]]$S <- sml[[1]]$S
          sml[[i]]$diagRP <- np$P
        }
    }
  }
  if (null.space.penalty) {
    nsm <- length(sml[[1]]$S)
    if (nsm == 1) {
      S11 <- sml[[1]]$S[[1]][1, 1]
      rank <- sml[[1]]$rank
      p <- ncol(sml[[1]]$X)
      if (is.null(rank) || max(abs(sml[[1]]$S[[1]] - diag(c(rep(S11, 
                                                                rank), rep(0, p - rank))))) > abs(S11) * .Machine$double.eps^0.8) 
        need.full <- TRUE
      else {
        need.full <- FALSE
        if (p > rank) 
          for (i in 1:length(sml)) {
            sml[[i]]$S[[2]] <- diag(c(rep(0, rank), rep(1, 
                                                        p - rank)))
            sml[[i]]$rank[2] <- p - rank
            sml[[i]]$S.scale[2] <- 1
            sml[[i]]$null.space.dim <- 0
          }
      }
    }
    else need.full <- if (nsm > 0) 
      TRUE
    else FALSE
    if (need.full) {
      St <- sml[[1]]$S[[1]]
      if (length(sml[[1]]$S) > 1) 
        for (i in 1:length(sml[[1]]$S)) St <- St + sml[[1]]$S[[i]]
        es <- eigen(St, symmetric = TRUE)
        ind <- es$values < max(es$values) * .Machine$double.eps^0.66
        if (sum(ind)) {
          U <- es$vectors[, ind, drop = FALSE]
          Sf <- U %*% t(U)
          M <- length(sm$S)
          for (i in 1:length(sml)) {
            sml[[i]]$S[[M + 1]] <- Sf
            sml[[i]]$rank[M + 1] <- sum(ind)
            sml[[i]]$S.scale[M + 1] <- 1
            sml[[i]]$null.space.dim <- 0
          }
        }
    }
  }
  sml
}

fullRank <-function (x, tol = 1e-07, qrx = qr(x, tol = tol)) 
{
  keep_cols <- drop_cols <- rank <- NULL
  d <- dim(x)
  n <- d[[1L]]
  p <- d[[2L]]
  if (n < p){
    return(t(fullRank(t(x), tol = tol))) #transpose & then perform the routine
  }
  
  rnk <- qrx$rank
  if (rnk == p){
    x_full <- x
  } else {
    keep_cols <- qrx$pivot[seq_len(rnk)]
    drop_cols <- rep(0,q)
    drop_cols <- as.numeric(1:q %in% keep_cols)
    names(drop_cols) <- 1:q
    drop_cols <- names(drop_cols)[drop_cols == 0]
    
    x <- x[, keep_cols, drop = FALSE]
  }
  return(list(x = x, drop_cols = drop_cols, keep_cols = keep_cols, rank = rnk))
}

d_null.coef_fn <- function(beta_star, sp, X){ 
  d_coefstart <- drho_d_coefstart <- d2rho_d_coefstart <- d2_coefstart <- drho_d2_coefstart <- 
    d2rho_d2_coefstart <- d2rho_coefstart <- drho_coefstart <- list()
  
  #We just initialize all coefficients as zero
  zero_vec <- rep(0,times = ncol(X))
  
  for(i in 1:length(sp))
  {
    d2rho_coefstart[[i]] <- list()
    for(j in i:length(sp))
    {
      d2rho_coefstart[[i]][[j]] <- zero_vec
    }
    drho_coefstart[[i]] <- zero_vec
  }
  
  for(d1 in 1:length(beta_star))
  { 
    drho_d_coefstart[[d1]] <- d2rho_d_coefstart[[d1]] <- list()
    d_coefstart[[d1]] <- zero_vec
    for(i in 1:length(sp))
    {
      drho_d_coefstart[[d1]][[i]] <- zero_vec
    }
    
    for(i in 1:length(sp))
    {
      d2rho_d_coefstart[[d1]][[i]] <- list()
      for(j in i:length(sp))
      { 
        d2rho_d_coefstart[[d1]][[i]][[j]] <- zero_vec
      }
    }
  }  
  
  for(d1 in 1:length(beta_star))
  {
    d2_coefstart[[d1]] <- drho_d2_coefstart[[d1]] <- d2rho_d2_coefstart[[d1]] <- list()
    for(d2 in d1:length(beta_star))
    {    
      drho_d2_coefstart[[d1]][[d2]] <- d2rho_d2_coefstart[[d1]][[d2]] <- list()
      d2_coefstart[[d1]][[d2]] <- zero_vec
      
      for(i in 1:length(sp))
      {
        drho_d2_coefstart[[d1]][[d2]][[i]] <- zero_vec
      }
      
      for(i in 1:length(sp))
      {
        d2rho_d2_coefstart[[d1]][[d2]][[i]] <- list()
        for(j in i:length(sp))
        {
          d2rho_d2_coefstart[[d1]][[d2]][[i]][[j]] <- zero_vec
        }
      }
    }    
  }
  
  list(drho_coefstart=drho_coefstart, d2rho_coefstart=d2rho_coefstart, d_coefstart=d_coefstart, d2_coefstart=d2_coefstart, drho_d_coefstart=drho_d_coefstart, 
       d2rho_d_coefstart=d2rho_d_coefstart, drho_d2_coefstart = drho_d2_coefstart, d2rho_d2_coefstart=d2rho_d2_coefstart)
}

#02
estimate.gam.si <- function (beta_star, index_data, splineDeg, splineRank, D_Design, D2_Design, S_uncons,
                             G, method, optimizer, control, in.out, scale, gamma, 
                             ...) {
  
  if (inherits(G$family, "extended.family")) {
    if (!(method %in% c("REML", "ML"))) 
      method <- "REML"
    if (optimizer[1] == "perf") 
      optimizer <- c("outer", "newton")
    if (inherits(G$family, "general.family")) {
      if (!is.null(G$offset)) 
        if (is.list(G$offset)) {
          for (i in 1:length(G$offset)) if (!is.null(G$offset[[i]])) 
            warning("sorry, general families currently ignore offsets")
        }
      else if (sum(G$offset != 0) > 0) 
        warning("sorry, general families currently ignore offsets")
      method <- "REML"
      G$Sl <- Sl.setup(G)
      G$X <- Sl.initial.repara(G$Sl, G$X, both.sides = FALSE)
      if (!is.null(G$family$available.derivs) && G$family$available.derivs == 
          1) 
        optimizer <- c("outer", "bfgs")
    }
  }
  if (!optimizer[1] %in% c("perf", "outer")) 
    stop("unknown optimizer")
  if (!method %in% c("GCV.Cp", "GACV.Cp", "REML", "P-REML", 
                     "ML", "P-ML")) 
    stop("unknown smoothness selection criterion")
  G$family <- mgcv:::fix.family(G$family)
  Trunc_Chol <- mini.roots.si(G$S, G$off, ncol(G$X), G$rank)
  G$rS <- Trunc_Chol$B
  rank <- Trunc_Chol$rank
  Pivoting <- Trunc_Chol$pivoting
  Piv <- Trunc_Chol$piv
  rs_out <- G$rS
  if (method %in% c("REML", "P-REML", "ML", "P-ML")) {
    if (optimizer[1] == "perf") {
      warning("Reset optimizer to outer/newton")
      optimizer <- c("outer", "newton")
    }
    reml <- TRUE
  } else reml <- FALSE
  
  Ssp <- totalPenaltySpace.si(G$S, G$H, G$off, ncol(G$X)) #Create a version of totalPenaltySpace.si which will calculate derivatives
  G$Eb <- Ssp$E
  G$U1 <- cbind(Ssp$Y, Ssp$Z) #This is the T_matrix so I need D_T, D2_T
  G$Mp <- ncol(Ssp$Z)
  #outputed the scaled version of ZSZ
  Scaled_S <- Ssp$St
  #output the eigenvalues
  E_vals <- Ssp$Values
  penaltyRank <- Ssp$penaltyRank
  d_Ssp <- d_totalPenaltySpace(G$intercept, beta_star, G$X_u, D_Design, D2_Design, S_uncons=G$S_uncons, Pivoting, evals=Ssp$Values, U1=G$U1, splineRank, St=Ssp$St, S=G$S, H=G$H, off=G$off, p=ncol(G$X), ind=Ssp$ind)
  
  G$D_U1 <- d_Ssp$D_U1
  G$D2_U1 <- d_Ssp$D2_U1
  
  G$Z <- d_Ssp$Z
  G$D_Z <- d_Ssp$D_Z
  G$D2_Z <- d_Ssp$D2_Z
  G$D_ZSZ <- d_Ssp$D_S
  G$D2_ZSZ <- d_Ssp$D2_S
  
  DrS <- d_rS.si(beta_star, S=G$S, penaltyRank=rank,splineRank=splineRank,Pivoting=Pivoting,
                 D_S_Pivoted=d_Ssp$D_S_Pivoted,D2_S_Pivoted=d_Ssp$D2_S_Pivoted,Piv=Trunc_Chol$piv, b=Trunc_Chol$b, off=G$off, rS=G$rS)
  
  G$D_rS <- DrS$D_Chol
  G$D2_rS <- DrS$D2_Chol 
  G$splineRank <- splineRank
  
  G$D_X_cons <- D_Design
  G$D2_X_cons <- D2_Design
  
  k0 <- G$off[1]
  k1 <- k0 + nrow(G$S[[1]]) - 1
  k2 <- k0 + splineRank -1
  for(d1 in 1:length(beta_star))
  {
    G$D_X_cons[[d1]] <- G$D_X_cons[[d1]][,1:k2]
    G$D_X_cons[[d1]][,k0:k1] <- (G$X_u[[1]] %*% G$D_Z[[d1]]) + (D_Design[[d1]][,k0:k2] %*% G$Z)
    G$D_X_cons[[d1]] <- cbind(G$D_X_cons[[d1]][,1:k1],matrix(0,nrow=nrow(G$D_X_cons[[d1]]),ncol=(ncol(G$X)-k1) ))
    for(d2 in d1:length(beta_star))
    {
      G$D2_X_cons[[d1]][[d2]] <- G$D2_X_cons[[d1]][[d2]][,1:k2]
      
      G$D2_X_cons[[d1]][[d2]][,k0:k1] <- ((D2_Design[[d1]][[d2]][,k0:k2] %*% G$Z) + (D_Design[[d1]][,k0:k2] %*% G$D_Z[[d2]]) + 
                                            (D_Design[[d2]][,k0:k2] %*% G$D_Z[[d1]]) + (G$X_u[[1]] %*% G$D2_Z[[d1]][[d2]]))
      
      G$D2_X_cons[[d1]][[d2]] <- cbind(G$D2_X_cons[[d1]][[d2]][,1:k1],matrix(0,nrow=nrow(G$D2_X_cons[[d1]][[d2]]),ncol=(ncol(G$X)-k1) ))
    }
  }
  
  G$d_X <- D_Design
  rm(D_Design)
  G$d2_X <- D2_Design
  rm(D2_Design)
  G$beta_star <- beta_star
  rm(beta_star)
  
  G$UrS <- G$D_UrS <- G$D2_UrS <- list()
  if (length(G$S) > 0) {
    for(i in 1:length(G$S))
    {
      G$UrS[[i]] <- t(Ssp$Y) %*% G$rS[[i]]
      if(i==1)
      {
        for(d1 in 1:length(G$beta_star))
        {
          G$D2_UrS[[d1]] <- list()
          
          G$D_UrS[[d1]] <- t(d_Ssp$D_Y[[d1]]) %*% G$rS[[i]] + t(Ssp$Y) %*% G$D_rS[[d1]]
          
          for(d2 in d1:length(G$beta_star))
          {
            G$D2_UrS[[d1]][[d2]] <- (
              t(d_Ssp$D2_Y[[d1]][[d2]]) %*% G$rS[[i]] + 
                t(d_Ssp$D_Y[[d1]]) %*% G$D_rS[[d2]] +
                t(d_Ssp$D_Y[[d2]]) %*% G$D_rS[[d1]] + 
                t(Ssp$Y) %*% G$D2_rS[[d1]][[d2]]
            )
          }
        }
      }
    }
  } else i <- 0  
  if (!is.null(G$H)) {
    G$UrS[[i + 1]] <- t(Ssp$Y) %*% mroot(G$H)
  }
  outer.looping <- ((!G$am && (optimizer[1] == "outer")) || 
                      reml || method == "GACV.Cp")
  fam.name <- G$family$family[1]
  if (scale == 0) {
    if (fam.name == "binomial" || fam.name == "poisson") 
      G$sig2 <- 1
    else G$sig2 <- -1
  } else {
    G$sig2 <- scale
  }
  if (reml) {
    criterion <- method
    if (fam.name == "binomial" || fam.name == "poisson") 
      scale <- 1
    if (inherits(G$family, "extended.family") && scale <= 
        0) {
      scale <- if (is.null(G$family$scale)) 
        1
      else G$family$scale
    }
  } else {
    if (scale == 0) {
      if (fam.name == "binomial" || fam.name == "poisson") 
        scale <- 1
      else scale <- -1
    }
    if (scale > 0) 
      criterion <- "UBRE"
    else {
      if (method == "GCV.Cp") 
        criterion <- "GCV"
      else criterion <- "GACV"
    }
  }
  if (substr(fam.name, 1, 17) == "Negative Binomial") {
    scale <- 1
    if (method %in% c("GCV.Cp", "GACV.Cp")) 
      criterion <- "UBRE"
  }
  if (scale > 0) {
    if (method == "P-ML") 
      criterion <- method <- "ML"
    else if (method == "P-REML") 
      criterion <- method <- "REML"
  }
  family <- G$family
  nb.fam.reset <- FALSE
  if (outer.looping) {
    fixedSteps <- if (inherits(G$family, "extended.family")) 
      0
    else control$outerPIsteps
    if (substr(G$family$family[1], 1, 17) == "Negative Binomial") {
      scale <- G$sig2 <- 1
      G$family <- negbin(max(family$getTheta()), link = family$link)
      nb.fam.reset <- TRUE
    }
  } else fixedSteps <- control$maxit + 2
  if (!is.null(G$family$preinitialize)) 
    eval(G$family$preinitialize)
  if (length(G$sp) > 0){
    #This is where we set up initial sp
    d_lsp2 <- list()
    d2_lsp2 <- list()
    initial_lsp2 <- initial.spg.si(x=G$X, x_uncons=G$X_uncons, d_x=G$D_X_cons, d2_x=G$D2_X_cons, beta_star=G$beta_star, d_S=G$D_ZSZ, d2_S=G$D2_ZSZ,
                                   y=G$y, weights=G$w, family=G$family, S=G$S, off=G$off, L=G$L, lsp0=G$lsp0, E = G$Eb, ...)
    
    lsp2 <- log(initial_lsp2$def.sp)
    d_lsp2[[1]] <- initial_lsp2$d_def.sp
    d2_lsp2[[1]] <- initial_lsp2$d2_def.sp
    
    for(i in 2:length(lsp2))
    {
      d_lsp2[[i]] <- rep(0, length(G$beta_star))
      d2_lsp2[[i]] <- matrix(0,nrow=length(G$beta_star),ncol=length(G$beta_star))
    }
  } else lsp2 <- rep(0, 0)
  if (outer.looping && !is.null(in.out)) { #We currently do not enter here
    ok <- TRUE
    if (is.null(in.out$sp) || is.null(in.out$scale)) 
      ok <- FALSE
    if (length(in.out$sp) != length(G$sp)) 
      ok <- FALSE
    if (!ok) 
      stop("in.out incorrect: see documentation")
    lsp <- log(in.out$sp)
  } else { # We enter here for the next fitting
    if (fixedSteps > 0) {
      #We don't actually get into gam.fit.si - will have to address later 
      object <- gam.fit.si(G, family = G$family, control = control, 
                           gamma = gamma, fixedSteps = fixedSteps, ...)
      lsp <- log(object$sp)
    }
    else {
      lsp <- lsp2 #We just set lsp=lsp2
      d_lsp <- d_lsp2
      d2_lsp <- d2_lsp2
    }
  }
  if (nb.fam.reset) 
    G$family <- family
  if (outer.looping) {
    #We enter here - This is where the outer smoothing is done
    if (is.null(in.out) && length(lsp) > 0) { print("We don't get here do we???")
      ind <- lsp > lsp2 + 5
      lsp[ind] <- lsp2[ind] + 5
      d_lsp[ind] <- d_lsp2[ind] #THIS STUFF ISN'T RIGHT, NEEDS TO BE CHANGED
      d2_lsp[ind] <- d2_lsp2[ind]
      ind <- lsp < lsp2 - 5
      lsp[ind] <- lsp2[ind] - 5
      d_lsp[ind] <- d_lsp2[ind]
      d2_lsp[ind] <- d2_lsp2[ind]
    }
    null.stuff <- if (inherits(G$family, "general.family")){ 
      list()
    } else get.null.coef.si(G, ...) #This get initial spline coefficients 
    
    if (fixedSteps > 0 && is.null(in.out)){ 
      mgcv.conv <- object$mgcv.conv
    } else mgcv.conv <- NULL
    if (criterion %in% c("REML", "ML") && scale <= 0) { #We use REML
      if (fixedSteps > 0) {
        log.scale <- log(sum(object$weights * object$residuals^2)/(G$n - 
                                                                     sum(object$edf)))
      }
      else {
        if (is.null(in.out)) {
          log.scale <- log(null.stuff$null.scale/10) #This is what we do
        }
        else {
          log.scale <- log(in.out$scale)
        }
      }
      lsp <- c(lsp, log.scale)
      d_lsp[[length(lsp)]] <- rep(0,times=length(lsp)-1)
      d2_lsp[[length(lsp)]] <- matrix(0,nrow=length(lsp)-1,ncol=length(lsp)-1)
      if (!is.null(G$L)) {
        G$L <- cbind(rbind(G$L, rep(0, ncol(G$L))), c(rep(0, 
                                                          nrow(G$L)), 1))
      }
      if (!is.null(G$lsp0)) 
        G$lsp0 <- c(G$lsp0, 0)
    }
    
    if (inherits(G$family, "extended.family") && !inherits(G$family, 
                                                           "general.family") && G$family$n.theta > 0) {
      #We also don't currently enter this if statement
      th0 <- G$family$getTheta()
      nth <- length(th0)
      nlsp <- length(lsp)
      ind <- 1:nlsp + nth
      lsp <- c(th0, lsp)
      if (!is.null(G$L) && nth > 0) {
        L <- rbind(cbind(diag(nth), matrix(0, nth, ncol(G$L))), 
                   cbind(matrix(0, nrow(G$L), nth), G$L))
        G$L <- L
      }
      if (!is.null(G$lsp0)) 
        G$lsp0 <- c(th0 * 0, G$lsp0)
    } else nth <- 0
    G$null.coef <- null.stuff$null.coef
    
    #We find the derivatives of the null coefficients here
    G$d_null.coef <- d_null.coef_fn(beta_star = G$beta_star, sp = G$sp, X = G$X)
    
    #This is where we perform smoothing and fitting - we pass d_lsp and d2_lsp 
    object <- gam.outer.si(lsp, d_lsp, d2_lsp, fscale = null.stuff$null.scale, 
                           family = G$family, control = control, criterion = criterion, 
                           method = method, optimizer = optimizer, scale = scale, 
                           gamma = gamma, G = G, ...)
    if (criterion %in% c("REML", "ML") && scale <= 0) 
      object$sp <- object$sp[-length(object$sp)]
    if (inherits(G$family, "extended.family") && nth > 0) 
      object$sp <- object$sp[-(1:nth)]
    object$mgcv.conv <- mgcv.conv
  }
  if (!inherits(G$family, "extended.family") && G$intercept && 
      any(G$offset != 0)) 
    object$null.deviance <- glm(object$y ~ offset(G$offset), 
                                family = object$family, weights = object$prior.weights)$deviance
  object$method <- criterion
  object$smooth <- G$smooth
  names(object$edf) <- G$term.names
  names(object$edf1) <- G$term.names
  if (!is.null(G$family$postproc)) 
    eval(G$family$postproc)
  if (!is.null(G$P)) {
    object$coefficients <- as.numeric(G$P %*% object$coefficients)
    object$Vp <- G$P %*% object$Vp %*% t(G$P)
    object$Ve <- G$P %*% object$Ve %*% t(G$P)
    rownames(object$Vp) <- colnames(object$Vp) <- G$term.names
    rownames(object$Ve) <- colnames(object$Ve) <- G$term.names
  }
  names(object$coefficients) <- G$term.names
  object$Pivoting <- Pivoting
  object$Piv <- Piv
  object$rS <- rs_out
  object$Scaled_S <- Scaled_S 
  object$E_vals <- E_vals
  object$penaltyRank <- penaltyRank
  
  object
}

#02a
d_rS.si <- function(beta_star,S,penaltyRank,splineRank,Pivoting,D_S_Pivoted,D2_S_Pivoted,Piv,b,off,rS){
  D_Chol <- D2_Chol <- list()
  
  for(deriv_1 in 1:length(beta_star))
  {
    D_Chol[[deriv_1]] <- matrix(0,ncol=ncol(rS[[1]]), nrow=nrow(rS[[1]]))
    D2_Chol[[deriv_1]] <- list()
    for(deriv_2 in deriv_1:length(beta_star))
    {
      D2_Chol[[deriv_1]][[deriv_2]] <- matrix(0,ncol=ncol(rS[[1]]), nrow=nrow(rS[[1]]))
      Chol_out <- suppressWarnings(D2Chol(P=S[[1]][Pivoting[[1]],Pivoting[[1]]], DP_list=D_S_Pivoted, D2P=D2_S_Pivoted[[deriv_1]][[deriv_2]], np=ncol(S[[1]]),beta_star=beta_star,deriv_1=deriv_1,deriv_2=deriv_2))    
      D2_rS <- t(Chol_out$D2)[, Piv[[1]], drop = FALSE]
      D2_rS <- t(D2_rS[1:penaltyRank[1], , drop = FALSE])
      D2_Chol[[deriv_1]][[deriv_2]][off[1]:(off[1] + nrow(b[[1]]) - 1), ] <- D2_rS     
    }    
    D_rS <- t(Chol_out$D)[, Piv[[1]], drop = FALSE]
    D_rS <- t(D_rS[1:penaltyRank[1], , drop = FALSE])    
    D_Chol[[deriv_1]][off[1]:(off[1] + nrow(b[[1]]) - 1), ] <- D_rS    
  }
  list(D_Chol=D_Chol,D2_Chol=D2_Chol)
}

#02b
initial.spg.si <- function (x, x_uncons, d_x, d2_x, beta_star, d_S, d2_S,  
                            y, weights, family, S, off, L = NULL, lsp0 = NULL, 
                            type = 1, start = NULL, mustart = NULL, etastart = NULL, 
                            E = NULL, ...) {
  if (length(S) == 0) 
    return(rep(0, 0))
  nobs <- nrow(x)
  if (is.null(mustart)) 
    mukeep <- NULL
  else mukeep <- mustart
  eval(family$initialize)
  if (inherits(family, "general.family")) { 
    lbb <- family$ll(y, x, start, weights, family, deriv = 1)$lbb
    lambda <- rep(0, length(S))
    for (i in 1:length(S)) {
      ind <- off[i]:(off[i] + ncol(S[[i]]) - 1)
      lami <- 1
      dlb <- -diag(lbb[ind, ind])
      dS <- diag(S[[i]])
      ind <- rowSums(abs(S[[i]])) > max(S[[i]]) * .Machine$double.eps^0.75
      dlb <- dlb[ind]
      dS <- dS[ind]
      while (mean(dlb/(dlb + lami * dS)) > 0.4) lami <- lami * 
        5
      while (mean(dlb/(dlb + lami * dS)) < 0.4) lami <- lami/5
      lambda[i] <- lami
    }
  }
  else { #This is the part we get into  
    if (is.null(mukeep)) {
      if (!is.null(start)) 
        etastart <- drop(x %*% start)
      if (!is.null(etastart)) 
        mustart <- family$linkinv(etastart)
    }
    else mustart <- mukeep
    if (inherits(family, "extended.family")) {
      theta <- family$getTheta()
      w <- 0.5 * drop(family$Dd(y, mustart, theta, weights)$EDmu2 * 
                        family$mu.eta(family$linkfun(mustart))^2)
    }
    else w <- drop(weights * family$mu.eta(family$linkfun(mustart))^2/family$variance(mustart))    
    w <- sqrt(w)
    if (type == 1) { #Type is indeed 1, we use this
      w_d2_x <- list()
      for(i in 1:length(beta_star))
      {
        w_d2_x[[i]] <- list()
        for(j in i:length(beta_star))
        {
          w_d2_x[[i]][[j]] <- w*d2_x[[i]][[j]]
        }
      }
      d_lambda <- d_initial.sp(X=w*x, derivX=lapply(d_x, "*", w), deriv2X=w_d2_x, beta_star=beta_star, S=S, derivS=d_S, deriv2S=d2_S, off=off, expensive = FALSE, XX = FALSE) 
    }
    else {  #Will probably need to address this later
      csX <- colSums((w * x)^2)
      lambda <- rep(0, length(S))
      for (i in 1:length(S)) {
        ind <- off[i]:(off[i] + ncol(S[[i]]) - 1)
        lambda[i] <- sum(csX[ind])/sqrt(sum(S[[i]]^2))
      }
    }
  }
  if (!is.null(L)) { print("IS NULL L")
    lsp <- log(lambda)
    if (is.null(lsp0)) 
      lsp0 <- rep(0, nrow(L))
    lsp <- as.numeric(coef(lm(lsp ~ L - 1 + offset(lsp0))))
    lambda <- exp(lsp)
  }
  d_lambda
}

#02c
d_initial.sp <- function (X, derivX, deriv2X, beta_star, S, derivS, deriv2S, off, expensive = FALSE, XX = FALSE) {
  n.p <- length(S)
  if (XX) 
    expensive <- FALSE
  def.sp <- array(0, n.p)
  if (n.p) {
    ldxx <- if (XX) 
      diag(X)
    else colSums(X * X)
    ldss <- ldxx * 0
    
    d_ldxx <- list()
    d2_ldxx <- list()
    
    d_ldss <- list()
    d2_ldss <- list()
    
    d_St <- list()
    d2_St <- list()
    
    d_pen <- list()
    d2_pen <- list()
    
    d_def.sp <- vector()
    d2_def.sp <- matrix(nrow=length(beta_star),ncol=length(beta_star))
    
    for(k in 1:length(beta_star))
    {
      d2_ldxx[[k]] <- list()
      d2_ldss[[k]] <- list()
      d2_St[[k]] <- list()
      d2_pen[[k]] <- list()
      d_ldxx[[k]] <- if (XX) 
        diag(derivX[[k]])
      else colSums(2* X[,1:ncol(derivX[[k]])] * derivX[[k]])
      d_ldss[[k]] <- d_ldxx[[k]] * 0
      if (expensive){ 
        d_St[[k]] <- matrix(0, ncol(X), ncol(X))
      }
      d_pen[[k]] <- rep(FALSE, length(d_ldxx[[k]]))
      
      for(j in k:length(beta_star))
      {
        d2_ldxx[[k]][[j]] <- if (XX) 
          diag(deriv2X[[k]][[j]])
        else colSums(2* derivX[[j]] * derivX[[k]] + 2* X[,1:ncol(deriv2X[[k]][[j]])] * deriv2X[[k]][[j]])
        d2_ldss[[k]][[j]] <- d2_ldxx[[k]][[j]] * 0
        
        if (expensive){ 
          d2_St[[k]][[j]] <- matrix(0, ncol(X), ncol(X))
        }
        d2_pen[[k]][[j]] <- rep(FALSE, length(d2_ldxx[[k]][[j]]))
      }
      
    }
    
    
    if (expensive) 
      St <- matrix(0, ncol(X), ncol(X))
    pen <- rep(FALSE, length(ldxx))
    
    for (i in 1:n.p) {
      maS <- max(abs(S[[i]]))
      rsS <- rowMeans(abs(S[[i]]))
      csS <- colMeans(abs(S[[i]]))
      dS <- diag(abs(S[[i]]))
      thresh <- .Machine$double.eps^0.8 * maS
      ind <- rsS > thresh & csS > thresh & dS > thresh
      ss <- diag(S[[i]])[ind]
      start <- off[i]
      finish <- start + ncol(S[[i]]) - 1
      xx <- ldxx[start:finish]
      xx <- xx[ind]
      pen[start:finish] <- pen[start:finish] | ind
      sizeXX <- mean(xx)
      sizeS <- mean(ss)
      if (sizeS <= 0) 
        stop(gettextf("S[[%d]] matrix is not +ve definite.", 
                      i))
      def.sp[i] <- sizeXX/sizeS
      ldss[start:finish] <- ldss[start:finish] + def.sp[i] * 
        diag(S[[i]])
      if (expensive) 
        St[start:finish, start:finish] <- St[start:finish, 
                                             start:finish] + def.sp[i] * S[[i]]
      #We only need derivatives with respect to beta^\star for i=1
      if(i==1)
      {
        for(k in 1:length(beta_star))
        {
          d_xx <- d_ldxx[[k]][start:finish]
          d_xx <- d_xx[ind]
          d_ss <- diag(derivS[[k]])[ind]
          d_pen[[k]] <- d_pen[[k]][start:finish] | ind
          d_sizeXX<- mean(d_xx)
          d_sizeS <- mean(d_ss)
          d_def.sp[k] <- d_sizeXX/sizeS - (sizeXX*d_sizeS)/(sizeS^2)
          
          d_ldss[[k]] <- d_ldss[[k]][start:finish] + d_def.sp[[k]]*diag(S[[i]]) + def.sp[i]*diag(derivS[[k]])
          
          if (expensive)
            d_St[[k]][start:finish, start:finish] <- d_St[[k]][start:finish, 
                                                               start:finish] + d_def.sp[[k]] * S[[i]] + def.sp[i] *derivS[[k]]
          for(j in k:length(beta_star))
          {
            dj_xx <- d_ldxx[[j]][start:finish]
            dj_xx <- dj_xx[ind]
            dj_ss <- diag(derivS[[j]])[ind]
            dj_pen<- rep(FALSE, length(d_ldxx[[j]]))
            dj_pen <- dj_pen[[k]][start:finish] | ind
            dj_sizeXX<- mean(dj_xx)
            dj_sizeS <- mean(dj_ss)
            
            d2_xx <- d2_ldxx[[k]][[j]][start:finish]
            d2_xx <- d2_xx[ind]
            
            d2_ss <- diag(deriv2S[[k]][[j]])[ind]
            d2_pen[[k]][[j]] <- d2_pen[[k]][[j]][start:finish] | ind
            
            d2_sizeXX<- mean(d2_xx)
            d2_sizeS <- mean(d2_ss)
            
            d2_def.sp[k,j] <- d2_sizeXX/sizeS - (d_sizeS*dj_sizeS)/(sizeS^2) - (dj_sizeXX*d_sizeS + sizeXX*d2_sizeS)/(sizeS^2) + (sizeXX*d_sizeXX*2*sizeS*dj_sizeS)/(sizeS^4)       
          }
        }
      }
    }
    
    if (expensive) { #We don't use expensive but for cases I do, need to extend this
      msp <- single.sp(X, St)
      if (msp > 0) 
        def.sp <- def.sp * msp
    }
    else { 
      ind <- ldss > 0 & pen
      ldxx <- ldxx[ind]
      ldss <- ldss[ind]
      while (mean(ldxx/(ldxx + ldss)) > 0.4) {
        d2_def.sp <- d2_def.sp * 10
        d_def.sp <- d_def.sp * 10
        def.sp <- def.sp * 10
        ldss <- ldss * 10
      }
      while (mean(ldxx/(ldxx + ldss)) < 0.4) {
        d2_def.sp <- d2_def.sp/10
        d_def.sp <- d_def.sp/10
        def.sp <- def.sp/10
        ldss <- ldss/10
      }
    }
  }
  list(def.sp=as.numeric(def.sp),d_def.sp=d_def.sp,d2_def.sp=d2_def.sp)
}

#02d
get.null.coef.si <- function (G, start = NULL, etastart = NULL, mustart = NULL, ...) {
  y <- G$y
  weights <- G$w
  nobs <- G$n
  family <- G$family
  eval(family$initialize)
  y <- as.numeric(y)
  mum <- mean(y) + 0 * y
  etam <- family$linkfun(mum)
  
  null.coef <- qr.coef(qr(G$X), etam)
  null.coef[is.na(null.coef)] <- 0
  null.scale <- sum(family$dev.resids(y, mum, weights))/nrow(G$X)
  list(null.coef = null.coef, null.scale = null.scale)
}

#02e
mini.roots.si <- function (S, off, np, rank = NULL) {
  b <- list()
  m <- length(S)
  pivoting <- list()
  piv <- list()
  if (m <= 0) 
    return(list())
  B <- S
  if (is.null(rank)) 
    rank <- rep(-1, m)
  for (i in 1:m) {
    b0 <- mroot.si(S[[i]], rank = rank[i])
    b[[i]] <- b0$L
    pivoting[[i]] <- b0$pivoting 
    piv[[i]] <- b0$piv
    B[[i]] <- matrix(0, np, ncol(b[[i]]))
    B[[i]][off[i]:(off[i] + nrow(b[[i]]) - 1), ] <- b[[i]]
  }
  list(B=B,pivoting=pivoting,piv=piv,b=b,rank=rank)
}

#02f
mroot.si <- function (A, rank = NULL, method = "chol") {
  if (is.null(rank)) 
    rank <- 0
  if (!isTRUE(all.equal(A, t(A)))) 
    stop("Supplied matrix not symmetric")
  if (method == "svd") {
    um <- La.svd(A)
    if (sum(um$d != sort(um$d, decreasing = TRUE)) > 0) 
      stop("singular values not returned in order")
    if (rank < 1) {
      rank <- dim(A)[1]
      if (um$d[1] <= 0) 
        rank <- 1
      else while (rank > 0 && (um$d[rank]/um$d[1] < .Machine$double.eps || 
                               all.equal(um$u[, rank], um$vt[rank, ]) != TRUE)) rank <- rank - 
          1
      if (rank == 0) 
        stop("Something wrong - matrix probably not +ve semi definite")
    }
    d <- um$d[1:rank]^0.5
    return(t(t(um$u[, 1:rank]) * as.vector(d)))
  }
  else if (method == "chol") {
    L <- suppressWarnings(chol(A, pivot = TRUE, tol = 0))
    piv <- order(attr(L, "pivot"))
    pivoting <- attr(L, "pivot")
    r <- attr(L, "rank")
    p <- ncol(L)
    if (r < p) 
      L[(r + 1):p, (r + 1):p] <- 0
    if (rank < 1) 
      rank <- r
    L <- L[, piv, drop = FALSE]
    L <- t(L[1:rank, , drop = FALSE])
    return(list(L=L,pivoting=pivoting,piv=piv))
  }
  else stop("method not recognised.")
}

#03
totalPenaltySpace.si <- function (S, H, off, p) {
  Hscale <- sqrt(sum(H * H))
  if (Hscale == 0) 
    H <- NULL
  if (is.null(H)) 
    St <- matrix(0, p, p)
  else {
    St <- H/sqrt(sum(H * H))
    if (ncol(H) != p || nrow(H) != p) 
      stop("H has wrong dimension")
  }
  m <- length(S)
  if (m > 0) 
    for (i in 1:m) { #For now, m=1 which will need to be ammended when there is more then one smooth
      k0 <- off[i]
      k1 <- k0 + nrow(S[[i]]) - 1
      St[k0:k1, k0:k1] <- St[k0:k1, k0:k1] + S[[i]]/sqrt(sum(S[[i]] * 
                                                               S[[i]])) #St is scaled S
    }
  es <- eigen(St, symmetric = TRUE)
  ind <- es$values > max(es$values) * .Machine$double.eps^0.66
  penaltyRank <- sum(ind)
  Y <- es$vectors[, ind, drop = FALSE]
  Z <- es$vectors[, !ind, drop = FALSE]
  E <- sqrt(as.numeric(es$values[ind])) * t(Y)
  
  list(Y = Y, Z = Z, E = E, St=St, Values=es$values,penaltyRank=penaltyRank,ind=ind)
}

#04
d_totalPenaltySpace <- function(intercept, beta_star, Design, D_Design, D2_Design, S_uncons, Pivoting, evals, U1, splineRank, St ,S, H, off, p, ind) {
  D_U1 <- D2_U1 <- D_evals <- D2_St_S <- D2beta_C <- D2_Z <- D_St_S <- D2_S <- D_S <- D_S_Pivoted <- D_Z <- Dbeta_C <- D_Identity <- D2_S_Pivoted <- list()
  D_Y <- D2_Y <- D_Zevecs <- D2_Zevecs <- list()
  if(intercept)
  {
    subset <- splineRank-1
    
    for(deriv_1 in 1:length(beta_star))
    {
      D_Design[[deriv_1]] <- D_Design[[deriv_1]][,-1]
      for(deriv_2 in deriv_1:length(beta_star))
      {
        D2_Design[[deriv_1]][[deriv_2]] <- D2_Design[[deriv_1]][[deriv_2]][,-1]
      }
    }
    
    # user  system elapsed 
    # 0.016   0.000   0.016 
  } else {
    subset <- splineRank
  }
  #At this stage, S is the constrained version already
  #Already have S here as S
  #NOTE, PROBABLKY DON'T NEED TO LOOP THROUGH ALL SMOOTHS GIVEN THE FACT THAT ONHLY THE SMOOTH ASSOCIATED TO THE INDEX HAS DEPENDENCTY ON BETA
  #SO IN THE CASE OF TAKING DERIVATIVES WITH RESPECT TO BETA, I DON'T NEED THE LOOP THROUGH LENGTH(S) - We make these ammendments and have backup2 saved
  
  k0 <- off[1]
  k1 <- k0 + nrow(S[[1]]) - 1
  C <- matrix(colMeans(Design[[1]][,1:splineRank]),nrow=1,ncol=splineRank) #We pass the matrix without intercept
  
  Identity <- diag(rep(1,times=nrow(t(C))))
  
  D2_Identity <- matrix(0,ncol=nrow(t(C)),nrow=nrow(t(C)))
  I_k1 <- diag(nrow=p,ncol=p)
  
  for(deriv_1 in 1:length(beta_star))
  {
    D_Identity[[deriv_1]] <- matrix(0,ncol=nrow(t(C)),nrow=nrow(t(C)))
  }
  
  for(deriv_1 in 1:length(beta_star))
  {
    Dbeta_C[[deriv_1]] <- t(matrix(colMeans(D_Design[[deriv_1]][,1:splineRank]),nrow=1,ncol=splineRank))
  }
  
  
  for(deriv_1 in 1:length(beta_star))
  {
    D_St_S[[deriv_1]] <- matrix(0, p, p)
    
    D2beta_C[[deriv_1]] <- list()
    D2_Z[[deriv_1]] <- list()
    D2_S[[deriv_1]] <- list()
    D2_St_S[[deriv_1]] <- list()
    
    for(deriv_2 in deriv_1:length(beta_star))
    {
      D2_St_S[[deriv_1]][[deriv_2]] <- matrix(0, p, p)
      
      D2beta_C[[deriv_1]][[deriv_2]] <- matrix(colMeans(D2_Design[[deriv_1]][[deriv_2]][,1:splineRank]),nrow=1,ncol=splineRank) 
      QR_diff <- QR.Hessian(Deriv2_Design=t(D2beta_C[[deriv_1]][[deriv_2]]),Deriv_Design_list=Dbeta_C,Design=t(C),Deriv2_B=D2_Identity,Deriv_B_list=D_Identity,B=Identity,deriv_1=deriv_1,deriv_2=deriv_2)     
      Constraint_matrix <- (QR_diff$B)[,(nrow(C)+1):ncol(C)] #This will be the same both times        
      D2_Z[[deriv_1]][[deriv_2]] <- QR_diff$D2_B[,(nrow(C)+1):ncol(C)]      
    }
    D_Z[[deriv_1]] <- QR_diff$D_B[,(nrow(C)+1):ncol(C)]
    D_S[[deriv_1]] <- crossprod(D_Z[[deriv_1]],S_uncons[[1]]) %*% Constraint_matrix + (crossprod(Constraint_matrix,S_uncons[[1]])) %*% D_Z[[deriv_1]]
    D_S_Pivoted[[deriv_1]] <- D_S[[deriv_1]][Pivoting[[1]],Pivoting[[1]]]
    
    D_St_S[[deriv_1]][k0:k1, k0:k1] <- D_St_S[[deriv_1]][k0:k1, k0:k1] + 
      
      S[[1]]*(-1*sum(S[[1]]*D_S[[deriv_1]])/(( sum(S[[1]]*S[[1]]) )^(1.5))) + 
      
      D_S[[deriv_1]]/sqrt(sum(S[[1]]*S[[1]]))
    
  }
  
  for(deriv_1 in 1:length(beta_star))
  {
    D2_S_Pivoted[[deriv_1]] <- list()
    for(deriv_2 in deriv_1:length(beta_star))
    {
      D2_S[[deriv_1]][[deriv_2]] <- crossprod(D2_Z[[deriv_1]][[deriv_2]],S_uncons[[1]]) %*% Constraint_matrix + 
        crossprod(D_Z[[deriv_1]],S_uncons[[1]]) %*% D_Z[[deriv_2]] + 
        crossprod(D_Z[[deriv_2]],S_uncons[[1]]) %*% D_Z[[deriv_1]] + 
        crossprod(Constraint_matrix,S_uncons[[1]]) %*% D2_Z[[deriv_1]][[deriv_2]]
      
      D2_St_S[[deriv_1]][[deriv_2]][k0:k1, k0:k1] <- D2_St_S[[deriv_1]][[deriv_2]][k0:k1, k0:k1] + 
        (D_S[[deriv_2]]*(-sum(S[[1]]*D_S[[deriv_1]])) + S[[1]]*(-sum(D_S[[deriv_2]]*D_S[[deriv_1]] + S[[1]]*D2_S[[deriv_1]][[deriv_2]])))/(sum(S[[1]]^2)^(1.5)) -
        (S[[1]]*(-sum(S[[1]]*D_S[[deriv_1]])*(1.5)*(sum(S[[1]]^2)^(0.5))*2*sum(S[[1]]*D_S[[deriv_2]])))/(sum(S[[1]]^2)^3) +
        D2_S[[deriv_1]][[deriv_2]]/sqrt(sum(S[[1]]^2)) - 
        (D_S[[deriv_1]]*0.5*sqrt(sum(S[[1]]^2))*2*sum(S[[1]]*D_S[[deriv_2]]))/sum(S[[1]]^2)
      
      D2_S_Pivoted[[deriv_1]][[deriv_2]] <- D2_S[[deriv_1]][[deriv_2]][Pivoting[[1]],Pivoting[[1]]]  
    }
  }
  
  for(deriv_1 in 1:length(beta_star))
  {
    D_U1[[deriv_1]] <- matrix(nrow=ncol(U1),ncol=ncol(U1))
    D_evals[[deriv_1]] <- vector()
    
    #Fast version
    tmp <- D_St_S[[deriv_1]] %*% U1
    tmp2 <- crossprod(U1,D_St_S[[deriv_1]])
    for(k in 1:ncol(U1))
    {
      if(sum(abs(tmp[,k])) == 0){
        D_U1[[deriv_1]][,k] <- rep(0, times = nrow(U1))
      } else{
        D_U1[[deriv_1]][,k] <- -1* ginv( St- evals[k]*I_k1 ) %*% tmp[,k]
      }
      
      if(sum(abs(tmp2[k,])) == 0){
        D_evals[[deriv_1]][k] <- 0
      } else{
        D_evals[[deriv_1]][k] <- tmp2[k,] %*% U1[,k]
      }
    }
    
    D_Y[[deriv_1]] <- D_U1[[deriv_1]][, ind, drop = FALSE]
    D_Zevecs[[deriv_1]] <- D_U1[[deriv_1]][, !ind, drop = FALSE]
  }
  
  for(deriv_1 in 1:length(beta_star))
  {
    D2_U1[[deriv_1]] <- D2_Y[[deriv_1]] <- D2_Zevecs[[deriv_1]] <- list()
    for(deriv_2 in deriv_1:length(beta_star))
    {
      D2_U1[[deriv_1]][[deriv_2]]<- matrix(nrow=ncol(U1),ncol=ncol(U1))
      
      tmp3 <- D2_St_S[[deriv_1]][[deriv_2]] %*% U1
      tmp4 <- D_St_S[[deriv_1]] %*% D_U1[[deriv_2]]
      tmp5 <- D_St_S[[deriv_1]] %*% U1
      
      for(k in 1:ncol(U1))
      {
        if(sum(abs(tmp3[,k])) == 0 & sum(abs(tmp4[,k])) == 0 & sum(abs(tmp5[,k])) == 0 ){
          
          D2_U1[[deriv_1]][[deriv_2]][,k] <- rep(0, times = nrow(U1))
          
        } else{
          regExpr0 <- St- evals[k]*I_k1
          regExpr <- ginv(regExpr0)
          
          D2_U1[[deriv_1]][[deriv_2]][,k] <-
            -1*(
              -regExpr %*% (D_St_S[[deriv_2]] - D_evals[[deriv_2]][k]*I_k1) %*% regExpr  +
                
                tcrossprod(regExpr) %*% crossprod(D_St_S[[deriv_2]] - D_evals[[deriv_2]][k]*I_k1, (I_k1 - (regExpr0) %*% regExpr)) +
                
                tcrossprod(I_k1 - regExpr %*% (regExpr0), D_St_S[[deriv_2]] - D_evals[[deriv_2]][k]*I_k1) %*% crossprod(regExpr)
            ) %*% tmp5[,k] +
            
            -1*regExpr %*% tmp3[,k] +
            -1*regExpr %*% tmp4[,k]
        }
      }
      
      D2_Y[[deriv_1]][[deriv_2]] <- D2_U1[[deriv_1]][[deriv_2]][, ind, drop = FALSE]
      D2_Zevecs[[deriv_1]] <- D2_U1[[deriv_1]][[deriv_2]][, !ind, drop = FALSE]
    }
  }
  
  list(D_U1=D_U1, D2_U1=D2_U1, D_S_Pivoted=D_S_Pivoted, D2_S_Pivoted=D2_S_Pivoted, D_Z=D_Z, D2_Z=D2_Z, Z=Constraint_matrix, D_S=D_S, D2_S=D2_S, D_Y=D_Y, D2_Y=D2_Y, D_Zevecs=D_Zevecs, D2_Zevecs=D2_Zevecs)
  
}

#05
gam.fit.si <- function (G, start = NULL, etastart = NULL, mustart = NULL, family = gaussian(), 
                        control = gam.control(), gamma = 1, fixedSteps = (control$maxit + 
                                                                            1), ...) {
  intercept <- G$intercept
  conv <- FALSE
  n <- nobs <- NROW(G$y)
  nvars <- NCOL(G$X)
  y <- G$y
  X <- G$X
  if (nvars == 0) 
    stop("Model seems to contain no terms")
  olm <- G$am
  find.theta <- FALSE
  if (substr(family$family[1], 1, 17) == "Negative Binomial") {
    Theta <- family$getTheta()
    if (length(Theta) == 1) {
      find.theta <- FALSE
      G$sig2 <- 1
    }
    else {
      if (length(Theta) > 2) 
        warning("Discrete Theta search not available with performance iteration")
      Theta <- range(Theta)
      T.max <- Theta[2]
      T.min <- Theta[1]
      Theta <- sqrt(T.max * T.min)
      find.theta <- TRUE
    }
    nb.link <- family$link
  }
  n.S <- length(G$S)
  if (n.S > 0) {
    S.size <- 0
    for (i in 1:n.S) S.size[i] <- mean(abs(G$S[[i]]))
  }
  weights <- G$w
  n.score <- sum(weights != 0)
  offset <- G$offset
  variance <- family$variance
  dev.resids <- family$dev.resids
  aic <- family$aic
  linkinv <- family$linkinv
  linkfun <- family$linkfun
  mu.eta <- family$mu.eta
  if (!is.function(variance) || !is.function(linkinv)) 
    stop("illegal `family' argument")
  valideta <- family$valideta
  if (is.null(valideta)) 
    valideta <- function(eta) TRUE
  validmu <- family$validmu
  if (is.null(validmu)) 
    validmu <- function(mu) TRUE
  if (is.null(mustart)) {
    eval(family$initialize)
  }
  else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  if (NCOL(y) > 1) 
    stop("y must be univariate unless binomial")
  coefold <- NULL
  eta <- if (!is.null(etastart)) 
    etastart
  else if (!is.null(start)) 
    if (length(start) != nvars) 
      stop(gettextf("Length of start should equal %d and correspond to initial coefs.", 
                    nvars))
  else {
    coefold <- start
    offset + as.vector(if (NCOL(G$X) == 1) 
      G$X * start
      else G$X %*% start)
  }
  else family$linkfun(mustart)
  mu <- linkinv(eta)
  if (!(validmu(mu) && valideta(eta))) 
    stop("Can't find valid starting values: please specify some")
  devold <- sum(dev.resids(y, mu, weights))
  boundary <- FALSE
  scale <- G$sig2
  msp <- G$sp
  magic.control <- list(tol = G$conv.tol, step.half = G$max.half, 
                        rank.tol = control$rank.tol)
  for (iter in 1:(control$maxit)) {
    good <- weights > 0
    varmu <- variance(mu)[good]
    if (any(is.na(varmu))) 
      stop("NAs in V(mu)")
    if (any(varmu == 0)) 
      stop("0s in V(mu)")
    mu.eta.val <- mu.eta(eta)
    if (any(is.na(mu.eta.val[good]))) 
      stop("NAs in d(mu)/d(eta)")
    good <- (weights > 0) & (mu.eta.val != 0)
    if (all(!good)) {
      conv <- FALSE
      warning(gettextf("No observations informative at iteration %d", 
                       iter))
      break
    }
    mevg <- mu.eta.val[good]
    mug <- mu[good]
    yg <- y[good]
    weg <- weights[good]
    var.mug <- variance(mug)
    G$y <- z <- (eta - offset)[good] + (yg - mug)/mevg
    w <- sqrt((weg * mevg^2)/var.mug)
    G$w <- w
    G$X <- X[good, , drop = FALSE]
    G$sig2 <- scale
    if (sum(!is.finite(G$y)) + sum(!is.finite(G$w)) > 0) 
      stop("iterative weights or data non-finite in gam.fit - regularization may help. See ?gam.control.")
    mr <- magic(G$y, G$X, msp, G$S, G$off, L = G$L, lsp0 = G$lsp0, 
                G$rank, G$H, matrix(0, 0, ncol(G$X)), G$w, gamma = gamma, 
                G$sig2, G$sig2 < 0, ridge.parameter = control$irls.reg, 
                control = magic.control, n.score = n.score, nthreads = control$nthreads)
    G$p <- mr$b
    msp <- mr$sp
    G$sig2 <- mr$scale
    G$gcv.ubre <- mr$score
    if (find.theta) {
      mv <- magic.post.proc(G$X, mr, w = G$w^2)
      G$edf <- mv$edf
      Theta <- mgcv:::mgcv.find.theta(Theta, T.max, T.min, weights, 
                                      good, mu, mu.eta.val, G, .Machine$double.eps^0.5)
      family <- do.call("negbin", list(theta = Theta, link = nb.link))
      variance <- family$variance
      dev.resids <- family$dev.resids
      aic <- family$aic
      family$Theta <- Theta
    }
    if (any(!is.finite(G$p))) {
      conv <- FALSE
      warning(gettextf("Non-finite coefficients at iteration %d", 
                       iter))
      break
    }
    start <- G$p
    eta <- drop(X %*% start)
    mu <- linkinv(eta <- eta + offset)
    eta <- linkfun(mu)
    dev <- sum(dev.resids(y, mu, weights))
    if (control$trace) 
      message(gettextf("Deviance = %s Iterations - %d", 
                       dev, iter, domain = "R-mgcv"))
    boundary <- FALSE
    if (!is.finite(dev)) {
      if (is.null(coefold)) 
        stop("no valid set of coefficients has been found:please supply starting values", 
             call. = FALSE)
      warning("Step size truncated due to divergence", 
              call. = FALSE)
      ii <- 1
      while (!is.finite(dev)) {
        if (ii > control$maxit) 
          stop("inner loop 1; can't correct step size")
        ii <- ii + 1
        start <- (start + coefold)/2
        eta <- drop(X %*% start)
        mu <- linkinv(eta <- eta + offset)
        eta <- linkfun(mu)
        dev <- sum(dev.resids(y, mu, weights))
      }
      boundary <- TRUE
      if (control$trace) 
        cat("Step halved: new deviance =", dev, "\n")
    }
    if (!(valideta(eta) && validmu(mu))) {
      warning("Step size truncated: out of bounds.", call. = FALSE)
      ii <- 1
      while (!(valideta(eta) && validmu(mu))) {
        if (ii > control$maxit) 
          stop("inner loop 2; can't correct step size")
        ii <- ii + 1
        start <- (start + coefold)/2
        eta <- drop(X %*% start)
        mu <- linkinv(eta <- eta + offset)
        eta <- linkfun(mu)
      }
      boundary <- TRUE
      dev <- sum(dev.resids(y, mu, weights))
      if (control$trace) 
        cat("Step halved: new deviance =", dev, "\n")
    }
    if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon || 
        olm || iter >= fixedSteps) {
      conv <- TRUE
      coef <- start
      break
    }
    else {
      devold <- dev
      coefold <- coef <- start
    }
  }
  if (!conv) {
    warning("Algorithm did not converge")
  }
  if (boundary) 
    warning("Algorithm stopped at boundary value")
  eps <- 10 * .Machine$double.eps
  if (family$family[1] == "binomial") {
    if (any(mu > 1 - eps) || any(mu < eps)) 
      warning("fitted probabilities numerically 0 or 1 occurred")
  }
  if (family$family[1] == "poisson") {
    if (any(mu < eps)) 
      warning("fitted rates numerically 0 occurred")
  }
  residuals <- rep(NA, nobs)
  residuals[good] <- z - (eta - offset)[good]
  wt <- rep(0, nobs)
  wt[good] <- w^2
  wtdmu <- if (intercept) 
    sum(weights * y)/sum(weights)
  else linkinv(offset)
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)
  mv <- magic.post.proc(G$X, mr, w = G$w^2)
  G$Vp <- mv$Vb
  G$hat <- mv$hat
  G$Ve <- mv$Ve
  G$edf <- mv$edf
  G$conv <- mr$gcv.info
  G$sp <- msp
  rank <- G$conv$rank
  aic.model <- aic(y, n, mu, weights, dev) + 2 * sum(G$edf)
  if (scale < 0) {
    gcv.ubre.dev <- n.score * dev/(n.score - gamma * sum(G$edf))^2
  }
  else {
    gcv.ubre.dev <- dev/n.score + 2 * gamma * sum(G$edf)/n.score - 
      G$sig2
  }
  list(coefficients = as.vector(coef), residuals = residuals, 
       fitted.values = mu, family = family, linear.predictors = eta, 
       deviance = dev, null.deviance = nulldev, iter = iter, 
       weights = wt, prior.weights = weights, df.null = nulldf, 
       y = y, converged = conv, sig2 = G$sig2, edf = G$edf, 
       edf1 = mv$edf1, hat = G$hat, R = mr$R, boundary = boundary, 
       sp = G$sp, nsdf = G$nsdf, Ve = G$Ve, Vp = G$Vp, rV = mr$rV, 
       mgcv.conv = G$conv, gcv.ubre = G$gcv.ubre, aic = aic.model, 
       rank = rank, gcv.ubre.dev = gcv.ubre.dev, scale.estimated = (scale < 
                                                                      0))
}

#06
gam.outer.si <- function(lsp, d_lsp, d2_lsp, fscale, family, control, method, optimizer, criterion, 
                         scale, gamma, G, ...) {
  if (is.null(optimizer[2])) 
    optimizer[2] <- "newton"
  if (!optimizer[2] %in% c("newton", "bfgs", "nlm", "optim", 
                           "nlm.fd")) 
    stop("unknown outer optimization method.")
  if (length(lsp) == 0) {
    optimizer[2] <- "no.sps"
  }
  nbGetTheta <- substr(family$family[1], 1, 17) == "Negative Binomial" && 
    length(family$getTheta()) > 1
  if (optimizer[2] == "nlm.fd") {
    #Don't get in here
    if (nbGetTheta) 
      stop("nlm.fd not available with negative binomial Theta estimation")
    if (method %in% c("REML", "ML", "GACV.Cp", "P-ML", "P-REML")) 
      stop("nlm.fd only available for GCV/UBRE")
    um <- nlm(full.score, lsp, typsize = lsp, fscale = fscale, 
              stepmax = control$nlm$stepmax, ndigit = control$nlm$ndigit, 
              gradtol = control$nlm$gradtol, steptol = control$nlm$steptol, 
              iterlim = control$nlm$iterlim, G = G, family = family, 
              control = control, gamma = gamma, ...)
    lsp <- um$estimate
    object <- attr(full.score(lsp, G, family, control, gamma = gamma, 
                              ...), "full.gam.object")
    object$gcv.ubre <- um$minimum
    object$outer.info <- um
    object$sp <- exp(lsp)
    return(object)
  }
  family <- fix.family.link(family)
  family <- fix.family.var(family)
  if (method %in% c("REML", "ML", "P-REML", "P-ML")) 
    family <- fix.family.ls(family)
  if (nbGetTheta) {
    #Don't get in here
    if (!(optimizer[2] %in% c("newton", "bfgs", "no.sps"))) {
      warning("only outer methods `newton' & `bfgs' supports `negbin' family and theta selection: reset")
      optimizer[2] <- "newton"
    }
    object <- mgcv:::gam.negbin(lsp, fscale, family, control, method, 
                                optimizer, gamma, G, ...)
  }
  else if (optimizer[2] == "newton" || optimizer[2] == "bfgs") { 
    #We enter this statement and the optimizer we use is newton
    if (optimizer[2] == "bfgs"){ 
      b <- bfgs(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, 
                UrS = G$UrS, L = G$L, lsp0 = G$lsp0, offset = G$offset, 
                U1 = G$U1, Mp = G$Mp, family = family, weights = G$w, 
                control = control, gamma = gamma, scale = scale, 
                conv.tol = control$newton$conv.tol, maxNstep = control$newton$maxNstep, 
                maxSstep = control$newton$maxSstep, maxHalf = control$newton$maxHalf, 
                printWarn = FALSE, scoreType = criterion, null.coef = G$null.coef, 
                pearson.extra = G$pearson.extra, dev.extra = G$dev.extra, 
                n.true = G$n.true, Sl = G$Sl, ...)
    }
    else{
      b <- newton.si(lsp = lsp, d_lsp = d_lsp, d2_lsp = d2_lsp, beta_star = G$beta_star, splineRank = G$splineRank, rS = G$rS, d_rS = G$D_rS, d2_rS = G$D2_rS, 
                     X = G$X, d_X = G$d_X, d2_X= G$d2_X, X_uncons=G$X_uncons, d_X_cons = G$D_X_cons, d2_X_cons = G$D2_X_cons,Z = G$Z, d_Z = G$D_Z, d2_Z = G$D2_Z,
                     y = G$y, Eb = G$Eb, 
                     UrS = G$UrS, d_UrS = G$D_UrS, d2_UrS = G$D2_UrS,
                     L = G$L, lsp0 = G$lsp0, offset = G$offset, 
                     U1 = G$U1, d_U1 = G$D_U1, d2_U1 = G$D2_U1, 
                     Mp = G$Mp, family = family, weights = G$w, 
                     control = control, gamma = gamma, scale = scale, 
                     conv.tol = control$newton$conv.tol, maxNstep = control$newton$maxNstep, 
                     maxSstep = control$newton$maxSstep, maxHalf = control$newton$maxHalf, 
                     printWarn = FALSE, scoreType = criterion, null.coef = G$null.coef, d_null.coef = G$d_null.coef,
                     pearson.extra = G$pearson.extra, dev.extra = G$dev.extra, 
                     n.true = G$n.true, Sl = G$Sl, ...)
    }
    
    
    object <- b$object
    object$REML <- object$REML1 <- object$REML2 <- object$GACV <- object$D2 <- object$P2 <- object$UBRE2 <- object$trA2 <- object$GACV1 <- object$GACV2 <- object$GCV2 <- object$D1 <- object$P1 <- NULL
    object$sp <- as.numeric(exp(b$lsp))
    object$gcv.ubre <- b$score
    b <- list(conv = b$conv, iter = b$iter, grad = b$grad, 
              hess = b$hess, score.hist = b$score.hist)
    object$outer.info <- b
  }
  else {
    args <- list(X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, 
                 offset = G$offset, U1 = G$U1, Mp = G$Mp, family = family, 
                 weights = G$w, control = control, scoreType = criterion, 
                 gamma = gamma, scale = scale, L = G$L, lsp0 = G$lsp0, 
                 null.coef = G$null.coef, n.true = G$n.true)
    if (optimizer[2] == "nlm") {
      #Don't get in here
      b <- nlm(gam4objective, lsp, typsize = lsp, fscale = fscale, 
               stepmax = control$nlm$stepmax, ndigit = control$nlm$ndigit, 
               gradtol = control$nlm$gradtol, steptol = control$nlm$steptol, 
               iterlim = control$nlm$iterlim, check.analyticals = control$nlm$check.analyticals, 
               args = args, ...)
      lsp <- b$estimate
    }
    else if (optimizer[2] == "optim") {
      #Don't get in here
      b <- optim(par = lsp, fn = gam2objective, gr = gam2derivative, 
                 method = "L-BFGS-B", control = list(fnscale = fscale, 
                                                     factr = control$optim$factr, lmm = min(5, length(lsp))), 
                 args = args, ...)
      lsp <- b$par
    }
    else b <- NULL
    obj <- gam2objective(lsp, args, printWarn = TRUE, ...)
    object <- attr(obj, "full.fit")
    object$gcv.ubre <- as.numeric(obj)
    object$outer.info <- b
    object$sp <- exp(lsp)
  }
  if (scale > 0) {
    object$scale.estimated <- FALSE
    object$scale <- scale
  }
  else {
    object$scale <- object$scale.est
    object$scale.estimated <- TRUE
  }
  object$control <- control
  if (inherits(family, "general.family")) {
    #Don't get in here
    mv <- mgcv:::gam.fit5.post.proc(object, G$Sl, G$L, G$S, G$off)
    object$coefficients <- Sl.initial.repara(G$Sl, object$coefficients, 
                                             inverse = TRUE)
  }
  else mv <- mgcv:::gam.fit3.post.proc(G$X, G$L, G$lsp0, G$S, G$off,object)
  if (!is.null(mv$Vc)) 
    object$Vc <- mv$Vc
  if (!is.null(mv$edf2)) 
    object$edf2 <- mv$edf2
  object$Vp <- mv$Vb
  object$hat <- mv$hat
  object$Ve <- mv$Ve
  object$edf <- mv$edf
  object$edf1 <- mv$edf1
  object$R <- mv$R
  object$aic <- object$aic + 2 * sum(mv$edf)
  object$nsdf <- G$nsdf
  object$K <- object$D1 <- object$D2 <- object$P <- object$P1 <- object$P2 <- object$GACV <- object$GACV1 <- object$GACV2 <- object$REML <- object$REML1 <- object$REML2 <- object$GCV <- object$GCV1 <- object$GCV2 <- object$UBRE <- object$UBRE1 <- object$UBRE2 <- object$trA <- object$trA1 <- object$trA2 <- object$alpha <- object$alpha1 <- object$scale.est <- NULL
  object$sig2 <- object$scale
  object
}

#07
newton.si <- function(lsp, d_lsp, d2_lsp, beta_star, splineRank, rS, d_rS, d2_rS,
                      X, d_X, d2_X, X_uncons,d_X_cons, d2_X_cons, Z, d_Z, d2_Z, 
                      y, Eb, UrS, d_UrS, d2_UrS, L, lsp0, offset, U1, d_U1, d2_U1, Mp, family, 
                      weights, control, gamma, scale, conv.tol = 1e-06, maxNstep = 5, 
                      maxSstep = 2, maxHalf = 30, printWarn = FALSE, scoreType = "deviance", 
                      start = NULL, mustart = NULL, null.coef = rep(0, ncol(X)), d_null.coef = NULL,
                      pearson.extra, dev.extra = 0, n.true = -1, Sl = NULL, edge.correct = FALSE, 
                      ...) 
{
  d_lsp_out <- d2_lsp_out <- list()
  reml <- scoreType %in% c("REML", "P-REML", "ML", "P-ML")
  if (is.null(L)){ 
    L <- diag(length(lsp))
  } else {
    if (!inherits(L, "matrix")) 
      stop("L must be a matrix.")
    if (nrow(L) < ncol(L)) 
      stop("L must have at least as many rows as columns.")
    if (nrow(L) != length(lsp0) || ncol(L) != length(lsp)) 
      stop("L has inconsistent dimensions.")
  }
  if (is.null(lsp0)) 
    lsp0 <- rep(0, nrow(L))
  if (reml && FALSE) {
    frob.X <- sqrt(sum(X * X))
    lsp.max <- rep(NA, length(lsp0))
    for (i in 1:nrow(L)) {
      lsp.max[i] <- 16 + log(frob.X/sqrt(sum(UrS[[i]]^2))) - 
        lsp0[i]
      if (lsp.max[i] < 2) 
        lsp.max[i] <- 2
    }
  } else lsp.max <- NULL
  if (!is.null(lsp.max)) {
    lsp1.max <- coef(lm(lsp.max - lsp0 ~ L - 1))
    ind <- lsp > lsp1.max
    lsp[ind] <- lsp1.max[ind] - 1
    delta <- rti(lsp, lsp1.max)
  } else {
    delta <- lsp
  }
  check.derivs <- FALSE
  sp.trace <- FALSE
  if (check.derivs) {
    deriv <- 2
    eps <- 1e-04
    deriv.check(x = X, y = y, sp = L %*% lsp + lsp0, Eb = Eb, 
                UrS = UrS, offset = offset, U1 = U1, Mp = Mp, family = family, 
                weights = weights, deriv = deriv, control = control, 
                gamma = gamma, scale = scale, printWarn = FALSE, 
                start = start, mustart = mustart, scoreType = scoreType, 
                eps = eps, null.coef = null.coef, Sl = Sl, ...)
  }
  print("first")
  
  d2_lsp <- sapply(d2_lsp, function(x) { x[lower.tri(x,diag = FALSE)] <- t(x)[lower.tri(x,diag = FALSE)] 
  x },simplify = FALSE) #assign lower triangular to upper
  
  L0 <- diag(length(d_lsp[[1]]))
  
  d_lsp_in <- lapply(d_lsp, function(x) L0 %*% x)
  d2_lsp_in <- lapply(d2_lsp, function(x) L0 %*% x)
  
  d_lsp_out[[1]] <- d_lsp_in
  d2_lsp_out[[1]] <- d2_lsp_in
  
  b <- gam.fit3.si(init_mod_pen = NULL, coef_derivs = NULL, beta_star = beta_star, splineRank = splineRank, rS = rS, d_rS = d_rS, d2_rS = d2_rS,
                   x = X, d_X = d_X, d2_X = d2_X, X_uncons=X_uncons, d_X_cons = d_X_cons, d2_X_cons = d2_X_cons,Z = Z, d_Z = d_Z, d2_Z = d2_Z,
                   y = y, sp = L %*% lsp + lsp0, d_sp = d_lsp_in, d2_sp = d2_lsp_in, Eb = Eb,
                   UrS = UrS, d_UrS = d_UrS, d2_UrS = d2_UrS,
                   offset = offset, U1 = U1, d_U1 = d_U1, d2_U1 = d2_U1,
                   Mp = Mp, family = family,
                   weights = weights, deriv = 2, control = control, gamma = gamma,
                   scale = scale, printWarn = FALSE, start = start, mustart = mustart,
                   scoreType = scoreType, null.coef = null.coef, d_null.coef = d_null.coef, pearson.extra = pearson.extra,
                   dev.extra = dev.extra, n.true = n.true, Sl = Sl, mgcv_version = FALSE, ...)
  
  modPen_out <- b$mod_pen
  deviance_gradient <- b$deviance_gradient
  deviance_hessian <- b$deviance_hessian
  smoothing_parameters <- L %*% lsp + lsp0
  
  mustart <- b$fitted.values
  etastart <- b$linear.predictors
  start <- b$coefficients
  coefDerivs_out <- b$coef_derivs
  if (reml) {
    old.score <- score <- b$REML
    grad <- b$REML1
    hess <- b$REML2
    
    d_old.score <- d_score <- b$d_REML
    d_grad <- b$d_REML1
    d_hess <- b$d_REML2
    
    d2_old.score <- d2_score <- b$d2_REML
    d2_grad <- b$d2_REML1
    d2_hess <- b$d2_REML2 
    
    d_hess <- lapply(d_hess, function(x) {x[lower.tri(x,diag = FALSE)] <- t(x)[lower.tri(x,diag = FALSE)]
    x})
    d2_hess <- lapply(d2_hess,function(x){
      out <- lapply(x,function(y){ 
        if(!is.null(y)) { y[lower.tri(y,diag = FALSE)] <- t(y)[lower.tri(y,diag = FALSE)] } 
        y })
      out })
  } else if (scoreType == "GACV") {
    old.score <- score <- b$GACV
    grad <- b$GACV1
    hess <- b$GACV2
  } else if (scoreType == "UBRE") {
    old.score <- score <- b$UBRE
    grad <- b$UBRE1
    hess <- b$UBRE2
  } else {
    old.score <- score <- b$GCV
    grad <- b$GCV1
    hess <- b$GCV2
  }
  grad <- t(L) %*% grad
  hess <- t(L) %*% hess %*% L
  
  d_grad <- lapply(d_grad, function(x) crossprod(L,x) )
  d_hess <- lapply(d_hess, function(x) crossprod(L,x) %*% L)
  
  d2_grad <- lapply(d2_grad, function(x) lapply(x, function(y){if(!is.null(y)) crossprod(L,y)}) )
  d2_hess <- lapply(d2_hess, function(x) lapply(x, function(y){if(!is.null(y)) crossprod(L,y) %*% L }) )
  
  if (!is.null(lsp.max)) {
    rho <- rt(delta, lsp1.max)
    nr <- length(rho$rho1)
    hess <- diag(rho$rho1, nr, nr) %*% hess %*% diag(rho$rho1, 
                                                     nr, nr) + diag(rho$rho2 * grad)
    grad <- rho$rho1 * grad
  }
  if (reml){ 
    score.scale <- abs(log(b$scale.est)) + abs(score)
  } else score.scale <- b$scale.est + abs(score)
  uconv.ind <- abs(grad) > score.scale * conv.tol
  if (!sum(uconv.ind)) 
    uconv.ind <- uconv.ind | TRUE
  score.hist <- rep(NA, 200)
  qerror.thresh <- 0.8
  for (i in 1:200) {
    if (control$trace) {
      cat("\n", i, "newton max(|grad|) =", max(abs(grad)), 
          "\n")
    }
    okc <- check.derivs
    while (okc) {
      okc <- FALSE
      eps <- 1e-04
      deriv <- 2
      deriv.check(x = X, y = y, sp = L %*% lsp + lsp0, 
                  Eb = Eb, UrS = UrS, offset = offset, U1 = U1, 
                  Mp = Mp, family = family, weights = weights, 
                  deriv = deriv, control = control, gamma = gamma, 
                  scale = scale, printWarn = FALSE, etastart = etastart, 
                  start = start, scoreType = scoreType, eps = eps, 
                  null.coef = null.coef, Sl = Sl, ...)
      if (inherits(family, "general.family")) {
        eps <- 1e-06
        spe <- 0.001
        er <- deriv.check5(x = X, y = y, sp = L %*% lsp + 
                             lsp0, weights = weights, start = start, offset = offset, 
                           Mp = Mp, family = family, control = control, 
                           deriv = deriv, eps = eps, spe = spe, Sl = Sl, 
                           ...)
      }
    }
    uconv.ind <- uconv.ind & abs(grad) > max(abs(grad)) * 0.001
    hess1 <- hess[uconv.ind, uconv.ind]
    grad1 <- grad[uconv.ind]
    
    d_grad1 <- lapply(d_grad, function(x) x[uconv.ind])
    d_hess1 <- lapply(d_hess, function(x) as.matrix(x[uconv.ind,uconv.ind]))
    
    d2_grad1 <- d2_hess1 <- list()
    for(d1 in 1:length(beta_star))
    {
      d2_grad1[[d1]] <- d2_hess1[[d1]] <- list()
      for(d2 in d1:length(beta_star))
      {
        d2_grad1[[d1]][[d2]] <- d2_grad[[d1]][[d2]][uconv.ind]
        d2_hess1[[d1]][[d2]] <- d2_hess[[d1]][[d2]][uconv.ind, uconv.ind]
      }
    }
    
    eh <- eigen(hess1, symmetric = TRUE)
    d <- eh$values
    U <- eh$vectors
    
    Id <- diag(nrow=nrow(U),ncol=ncol(U))
    
    d_U <- lapply(d_hess1, function(x){
      out <- as.matrix(sapply(1:ncol(U),function(k) {
        inner <- -1* ginv( hess1- (d[k]*Id) ) %*% (x %*% U[,k])
        inner
      }))
      out
    })
    
    d_d <- lapply(d_hess1, function(x){
      out <- unlist(lapply(1:ncol(U), function(k) {
        inner <- crossprod(U,x)[k,] %*% U[,k] 
        inner
      }))
      out
    })
    
    d2_U <- d2_d <- list()
    for(d1 in 1:length(beta_star))
    {
      d2_U[[d1]] <- list()
      d2_d[[d1]] <- list()
      for(d2 in d1:length(beta_star))
      {
        d2_U[[d1]][[d2]]<- matrix(nrow=ncol(U),ncol=ncol(U))
        d2_d[[d1]][[d2]] <- vector()
        
        for(k in 1:ncol(U))
        {
          d2_U[[d1]][[d2]][,k] <- (
            -1*(-ginv(hess1- (d[k]*Id)) %*% ( d_hess1[[d2]]-(d_d[[d2]][k]*Id) ) %*% ginv(hess1- (d[k]*Id))
                + ginv(hess1- (d[k]*Id)) %*% t(ginv(hess1- (d[k]*Id))) %*% t( d_hess1[[d2]]-(d_d[[d2]][k]*Id) ) %*% ( Id-(hess1- (d[k]*Id)) %*% ginv(hess1- (d[k]*Id)) )
                + ( Id-(hess1- (d[k]*Id)) %*% ginv(hess1- (d[k]*Id)) ) %*% t( d_hess1[[d2]]-(d_d[[d2]][k]*Id) ) %*% t(ginv(hess1- (d[k]*Id))) %*% ginv(hess1- (d[k]*Id))
            ) %*% d_hess1[[d1]] %*% U[,k]
            - ginv(hess1- (d[k]*Id)) %*% d2_hess1[[d1]][[d2]] %*% U[,k]
            - ginv(hess1- (d[k]*Id)) %*% d_hess1[[d1]] %*% d_U[[d2]][,k]
          )
          
          d2_d[[d1]][[d2]][k] <-  ( (crossprod(d_U[[d2]],d_hess1[[d1]]))[k,] %*% U[,k]) + 
            ( (crossprod(U,d2_hess1[[d1]][[d2]] ))[k,] %*% U[,k] ) + 
            ( (crossprod(U,d_hess1[[d1]]))[k,] %*% d_U[[d2]][,k] )
        }
      }
    }
    
    ind <- d < 0
    pdef <- if (sum(ind) > 0){ 
      FALSE
    } else TRUE
    
    for(d1 in 1:length(beta_star))
    {
      for(d2 in d1:length(beta_star))
      {
        d2_d[[d1]][[d2]][ind] <- -d2_d[[d1]][[d2]][ind]
      }
    }
    
    d_d <- lapply(d_d,function(x){
      x[ind] <- -x[ind]
      x
    })
    
    d[ind] <- -d[ind]
    
    low.d <- max(d) * .Machine$double.eps^0.7
    ind <- d < low.d
    if (sum(ind) > 0) 
      pdef <- FALSE
    
    for(d1 in 1:length(beta_star))
    {
      for(d2 in d1:length(beta_star))
      {
        d2_d[[d1]][[d2]][ind] <- d2_d[[d1]][[d2]][which(d==max(d))] * .Machine$double.eps^0.7
      }
    }
    
    d_d <- lapply(d_d,function(x){
      x[ind] <- x[which(d==max(d))] * .Machine$double.eps^0.7
      x
    })
    
    d[ind] <- low.d
    
    ind <- d != 0
    
    for(d1 in 1:length(beta_star))
    {
      for(d2 in d1:length(beta_star))
      {
        d2_d[[d1]][[d2]][ind] <- 
          
          -d2_d[[d1]][[d2]][ind]/(d[ind]*d[ind]) + 
          
          (2 * d[ind] * d_d[[d1]][ind] * d_d[[d2]][ind] )/( d[ind]*d[ind]*d[ind]*d[ind] )
        
      }
    }
    
    d_d <- lapply(d_d,function(x){
      x[ind] <- -x[ind]/(d[ind]*d[ind])
      x
    })
    
    d[ind] <- 1/d[ind]
    
    Nstep <- 0 * grad
    Nstep[uconv.ind] <- -drop(U %*% (d * (t(U) %*% grad1)))
    
    d_Nstep <- mapply(function(d_U,d_d,d_grad1){
      out <- 0*grad
      out[uconv.ind] <- -drop( d_U %*% (d * t(U) %*% grad1) + 
                                 U %*% (d_d * t(U) %*% grad1) + 
                                 U %*% (d * t(d_U) %*% grad1) + 
                                 U %*% (d * t(U) %*% d_grad1) )
      out
      
    },d_U,d_d,d_grad1,SIMPLIFY = FALSE)
    
    d2_Nstep <- list()
    for(d1 in 1:length(beta_star)){
      
      d2_Nstep[[d1]] <- list()
      for(d2 in d1:length(beta_star))
      {
        d2_Nstep[[d1]][[d2]] <- 0 * grad
        d2_Nstep[[d1]][[d2]][uconv.ind] <- -drop( (d2_U[[d1]][[d2]] %*% (d * (t(U) %*% grad1))) + 
                                                    (d_U[[d1]] %*% (d_d[[d2]] * (t(U) %*% grad1))) + 
                                                    (d_U[[d1]] %*% (d * (t(d_U[[d2]]) %*% grad1))) + 
                                                    (d_U[[d1]] %*% (d * (t(U) %*% d_grad1[[d2]]))) +
                                                    
                                                    (d_U[[d2]] %*% (d_d[[d1]] * (t(U) %*% grad1))) +
                                                    (U %*% (d2_d[[d1]][[d2]] * (t(U) %*% grad1))) + 
                                                    (U %*% (d_d[[d1]] * (t(d_U[[d2]]) %*% grad1))) + 
                                                    (U %*% (d_d[[d1]] * (t(U) %*% d_grad1[[d2]]))) +
                                                    
                                                    (d_U[[d2]] %*% (d * (t(d_U[[d1]]) %*% grad1))) + 
                                                    (U %*% (d_d[[d2]] * (t(d_U[[d1]]) %*% grad1))) + 
                                                    (U %*% (d * (t(d2_U[[d1]][[d2]]) %*% grad1))) + 
                                                    (U %*% (d * (t(d_U[[d1]]) %*% d_grad1[[d2]]))) +
                                                    
                                                    (d_U[[d2]] %*% (d * (t(U) %*% d_grad1[[d1]]))) + 
                                                    (U %*% (d_d[[d2]] * (t(U) %*% d_grad1[[d1]]))) + 
                                                    (U %*% (d * (t(d_U[[d2]]) %*% d_grad1[[d1]]))) + 
                                                    (U %*% (d * (t(U) %*% d2_grad1[[d1]][[d2]])))
        )
      }
    }
    
    Sstep <- -grad/max(abs(grad))
    
    d_Sstep <- lapply(d_grad, function(x){
      out <- -( x/max(abs(grad)) - (grad*x[which.max(abs(grad))]*grad[which.max(abs(grad))]) /(max(abs(grad))^3) )
      
      out
    })
    
    d2_Sstep <- list()
    for(d1 in 1:length(beta_star))
    {
      d2_Sstep[[d1]] <- list()
      for(d2 in d1:length(beta_star))
      {
        d2_Sstep[[d1]][[d2]] <-  -d2_grad[[d1]][[d2]]/max(abs(grad)) +
          
          ( d_grad[[d1]]*d_grad[[d2]][which.max(abs(grad))]*grad[which.max(abs(grad))] )/(max(abs(grad))^3) +
          
          (d_grad[[d2]]*d_grad[[d1]][which.max(abs(grad))]*grad[which.max(abs(grad))] +
             grad*d2_grad[[d1]][[d2]][which.max(abs(grad))]*grad[which.max(abs(grad))] +
             grad*d_grad[[d1]][which.max(abs(grad))]*d_grad[[d2]][which.max(abs(grad))]
          )/(max(abs(grad))^3) -
          
          (grad*d_grad[[d1]][which.max(abs(grad))]*grad[which.max(abs(grad))]*d_grad[[d2]][which.max(abs(grad))]*
             grad[which.max(abs(grad))]*3)/(max(abs(grad))^5)
      }
    }
    
    ms <- max(abs(Nstep))
    
    ind_max <- which.max(abs(Nstep))
    
    d_ms <- unlist(lapply(d_Nstep, function(x){
      out <- x[ind_max] * Nstep[ind_max]/abs(Nstep[ind_max])
      out
    }))
    
    d2_ms <- list()
    for(d1 in 1:length(beta_star))
    {
      d2_ms[[d1]] <- list()
      for(d2 in d1:length(beta_star))
      {
        d2_ms[[d1]][[d2]] <- d2_Nstep[[d1]][[d2]][ind_max] * Nstep[ind_max]/abs(Nstep[ind_max]) +
          
          d_Nstep[[d1]][ind_max] * (d_Nstep[[d2]][ind_max]/ms - (Nstep[ind_max]^2 * d_Nstep[[d2]][ind_max])/(ms^3) )
      }
    }
    
    mns <- maxNstep
    if (ms > maxNstep){ print("ms > maxNstep")
      for(d1 in 1:length(beta_star))
      {
        for(d2 in d1:length(beta_star))
        {
          d2_Nstep[[d1]][[d2]] <- mns*( d2_Nstep[[d1]][[d2]]/ms - 
                                          
                                          (d_Nstep[[d1]]*d_ms[[d2]])/(ms*ms) -
                                          
                                          (d_ms[[d1]]*d_Nstep[[d2]] + Nstep*d2_ms[[d1]][[d2]])/(ms*ms) + 
                                          
                                          ( (2*ms*d_ms[[d2]]*d_ms[[d1]]*Nstep)/(ms*ms*ms*ms) ) )
        }
      }
      
      d_Nstep <- mapply(function(x,y){
        out <- mns*( x/ms - y*Nstep/(ms*ms) )
      },d_Nstep,d_ms,SIMPLIFY = FALSE)
      
      Nstep <- mns * Nstep/ms
    }
    sd.unused <- TRUE
    if (sp.trace) 
      cat(lsp, "\n")
    if (!is.null(lsp.max)) {
      delta1 <- delta + Nstep
      lsp1 <- rt(delta1, lsp1.max)$rho
      while (max(abs(lsp1 - lsp)) > maxNstep) {
        Nstep <- Nstep/2
        delta1 <- delta + Nstep
        lsp1 <- rt(delta1, lsp1.max)$rho
      }
    } else{
      d2_lsp1 <- list()
      for(i in 1:length(lsp))
      {
        d2_lsp1[[i]] <- matrix(nrow=length(beta_star),ncol=length(beta_star))
        
        for(d1 in 1:length(beta_star))
        {
          for(d2 in d1:length(beta_star))
          {
            d2_lsp1[[i]][d1,d2] <- d2_lsp[[i]][d1,d2]*exp(Nstep[i]) + 
              d_lsp[[i]][d1]*d_Nstep[[d2]][i]*exp(Nstep[i]) +
              d_lsp[[i]][d2]*d_Nstep[[d1]][i]*exp(Nstep[i]) + 
              exp(lsp[i])*d2_Nstep[[d1]][[d2]][i]*exp(Nstep[i]) + 
              exp(lsp[i])*d_Nstep[[d1]][i]*d_Nstep[[d2]][i]*exp(Nstep[i])
          }
        }
      }
      
      #This odd expression is due to the fact I worked with d_lsp = d/dbeta^\star (exp(sp)) 
      d_lsp1 <- lapply(1:length(lsp), function(i){
        out <- unlist(lapply(1:length(beta_star), function(d1){
          out2 <- d_lsp[[i]][d1]*exp(Nstep[i]) + exp(lsp[i])*d_Nstep[[d1]][i]*exp(Nstep[i])
          out2
        }))
        out
      })
      
      lsp1 <- lsp + Nstep
      
      d2_lsp1 <- sapply(d2_lsp1, function(x) { x[lower.tri(x,diag = FALSE)] <- t(x)[lower.tri(x,diag = FALSE)] 
      x },simplify = FALSE) #assign lower triangular to upper
      
      L0 <- diag(length(d_lsp[[1]]))
      
      d_lsp1_in <- lapply(d_lsp1, function(x) L0 %*% x)
      d2_lsp1_in <- lapply(d2_lsp1, function(x) L0 %*% x)
      
    }
    print("Second")
    
    d_lsp_out[[length(d_lsp_out) + 1]] <- d_lsp1_in
    d2_lsp_out[[length(d2_lsp_out) + 1]] <- d2_lsp1_in
    
    b <- gam.fit3.si(init_mod_pen = modPen_out, coef_derivs = coefDerivs_out, beta_star = beta_star, splineRank = splineRank, 
                     rS = rS, d_rS = d_rS, d2_rS = d2_rS, x = X, d_X = d_X, d2_X = d2_X, X_uncons=X_uncons, 
                     d_X_cons = d_X_cons, d2_X_cons = d2_X_cons, Z = Z, d_Z = d_Z, d2_Z = d2_Z,
                     y = y, sp = L %*% lsp1 + lsp0, d_sp = d_lsp1_in, d2_sp = d2_lsp1_in, Eb = Eb, 
                     UrS = UrS, d_UrS = d_UrS, d2_UrS = d2_UrS,
                     offset = offset, U1 = U1, d_U1 = d_U1, d2_U1 = d2_U1,
                     Mp = Mp, family = family, 
                     weights = weights, deriv = as.numeric(pdef) * 2, control = control, gamma = gamma, 
                     scale = scale, printWarn = FALSE, start = start, mustart = mustart, 
                     scoreType = scoreType, null.coef = null.coef, d_null.coef = d_null.coef, 
                     pearson.extra = pearson.extra, 
                     dev.extra = dev.extra, n.true = n.true, Sl = Sl, mgcv_version = FALSE, ...)
    
    modPen_out <- b$mod_pen
    
    deviance_gradient <- b$deviance_gradient
    deviance_hessian <- b$deviance_hessian
    smoothing_parameters <- c(smoothing_parameters,L %*% lsp1 + lsp0)
    
    pred.change <- sum(grad * Nstep) + 0.5 * t(Nstep) %*% hess %*% Nstep
    if (reml) {
      score1 <- b$REML
      d_score1 <- b$d_REML
      d2_score1 <- b$d2_REML
    } else if (scoreType == "GACV") {
      score1 <- b$GACV
    } else if (scoreType == "UBRE") {
      score1 <- b$UBRE
    } else score1 <- b$GCV
    ii <- 0
    score.change <- score1 - score
    qerror <- abs(pred.change - score.change)/(max(abs(pred.change), 
                                                   abs(score.change)) + score.scale * conv.tol)
    if (is.finite(score1) && score.change < 0 && pdef && qerror < qerror.thresh)
    { 
      
      print("is.finite(score1) && score.change < 0 && pdef && qerror < qerror.thresh")
      
      old.score <- score
      d_old.score <- d_score
      d2_old.score <- d2_score
      
      mustart <- b$fitted.values
      etastart <- b$linear.predictors
      start <- b$coefficients
      coefDerivs_out <- b$coef_derivs
      
      lsp <- lsp1
      d_lsp <- d_lsp1
      d2_lsp <- d2_lsp1
      
      if (reml) {
        score <- b$REML
        grad <- b$REML1
        hess <- b$REML2
        
        d_score <- b$d_REML
        d_grad <- b$d_REML1
        d_hess <- b$d_REML2
        
        d2_score <- b$d2_REML
        d2_grad <- b$d2_REML1
        d2_hess <- b$d2_REML2
        
        d_hess <- lapply(d_hess, function(x) {x[lower.tri(x,diag = FALSE)] <- t(x)[lower.tri(x,diag = FALSE)]
        x})
        
        d2_hess <- lapply(d2_hess,function(x){
          out <- lapply(x,function(y){
            if(!is.null(y)) { y[lower.tri(y,diag = FALSE)] <- t(y)[lower.tri(y,diag = FALSE)] }
            y })
          out })
        
      } else if (scoreType == "GACV") {
        score <- b$GACV
        grad <- b$GACV1
        hess <- b$GACV2
      } else if (scoreType == "UBRE") {
        score <- b$UBRE
        grad <- b$UBRE1
        hess <- b$UBRE2
      } else {
        score <- b$GCV
        grad <- b$GCV1
        hess <- b$GCV2
      }
      
      grad <- t(L) %*% grad
      hess <- t(L) %*% hess %*% L
      
      d_grad <- lapply(d_grad, function(x) crossprod(L,x) )
      d_hess <- lapply(d_hess, function(x) crossprod(L,x) %*% L)
      
      d2_grad <- lapply(d2_grad, function(x) lapply(x, function(y){if(!is.null(y)) crossprod(L,y)}) )
      d2_hess <- lapply(d2_hess, function(x) lapply(x, function(y){if(!is.null(y)) crossprod(L,y) %*% L }) )
      
      if (!is.null(lsp.max)) {
        delta <- delta1
        rho <- rt(delta, lsp1.max)
        nr <- length(rho$rho1)
        hess <- diag(rho$rho1, nr, nr) %*% hess %*% diag(rho$rho1, 
                                                         nr, nr) + diag(rho$rho2 * grad)
        grad <- rho$rho1 * grad
      }
    } else if (!is.finite(score1) || score1 >= score || qerror >= qerror.thresh) { 
      
      print("!is.finite(score1) || score1 >= score || qerror >=qerror.thresh")
      
      step <- Nstep
      d_step <- d_Nstep
      d2_step <- d2_Nstep
      
      while ((!is.finite(score1) || score1 >= score || 
              qerror >= qerror.thresh) && ii < maxHalf) { 
        
        print("(!is.finite(score1) || score1 >= score || 
              qerror >= qerror.thresh) && ii < maxHalf")
        if (ii == 3 && i < 10) {
          print("BAD ii fix!")
          s.length <- min(sum(step^2)^0.5, maxSstep)
          step <- Sstep * s.length/sum(Sstep^2)^0.5
          sd.unused <- FALSE
        } else{
          step <- step/2
          
          d_step <- lapply(d_step,function(x){ out <- x/2; out})
          
          for(d1 in 1:length(beta_star))
          {
            for(d2 in d1:length(beta_star))
            {
              d2_step[[d1]][[d2]] <- d2_step[[d1]][[d2]]/2
            }
          }
        }
        if (!is.null(lsp.max)) {
          delta1 <- delta + step
          lsp1 <- rt(delta1, lsp1.max)$rho
        }
        else{ 
          d2_lsp1 <- list()
          for(i in 1:length(lsp))
          {
            d2_lsp1[[i]] <- matrix(nrow=length(beta_star),ncol=length(beta_star))
            
            for(d1 in 1:length(beta_star))
            {
              for(d2 in d1:length(beta_star))
              {
                d2_lsp1[[i]][d1,d2] <- d2_lsp[[i]][d1,d2]*exp(step[i]) + 
                  d_lsp[[i]][d1]*d_step[[d2]][i]*exp(step[i]) +
                  d_lsp[[i]][d2]*d_step[[d1]][i]*exp(step[i]) + 
                  exp(lsp[i])*d2_step[[d1]][[d2]][i]*exp(step[i]) + 
                  exp(lsp[i])*d_step[[d1]][i]*d_step[[d2]][i]*exp(step[i])
              }
            }
          }
          
          #This odd expression is due to the fact I worked with d_lsp = d/dbeta^\star (exp(sp)) which is somewhat questionable now
          d_lsp1 <- lapply(1:length(lsp), function(i){
            out <- unlist(lapply(1:length(beta_star), function(d1){
              out2 <- d_lsp[[i]][d1]*exp(step[i]) + exp(lsp[i])*d_step[[d1]][i]*exp(step[i])
              out2
            }))
            out
          })
          
          lsp1 <- lsp + step
        }
        
        d2_lsp1 <- sapply(d2_lsp1, function(x) { x[lower.tri(x,diag = FALSE)] <- t(x)[lower.tri(x,diag = FALSE)] 
        x },simplify = FALSE) #assign lower triangular to upper
        
        L0 <- diag(length(d_lsp[[1]]))
        
        d_lsp1_in <- lapply(d_lsp1, function(x) L0 %*% x)
        d2_lsp1_in <- lapply(d2_lsp1, function(x) L0 %*% x)
        
        print("Third")
        
        d_lsp_out[[length(d_lsp_out) + 1]] <- d_lsp1_in
        d2_lsp_out[[length(d2_lsp_out) + 1]] <- d2_lsp1_in
        
        b1 <- gam.fit3.si(init_mod_pen = modPen_out, coef_derivs = coefDerivs_out, beta_star = beta_star, splineRank = splineRank, rS = rS, 
                          d_rS = d_rS, d2_rS = d2_rS, x = X, d_X = d_X, d2_X = d2_X, X_uncons=X_uncons, 
                          d_X_cons = d_X_cons, d2_X_cons = d2_X_cons, Z = Z, d_Z = d_Z, d2_Z = d2_Z,
                          y = y, sp = L %*% lsp1 + lsp0, d_sp = d_lsp1_in, d2_sp = d2_lsp1_in, Eb = Eb, 
                          UrS = UrS, d_UrS = d_UrS, d2_UrS = d2_UrS,
                          offset = offset, U1 = U1, d_U1 = d_U1, d2_U1 = d2_U1,
                          Mp = Mp, family = family, 
                          weights = weights, deriv = 0, control = control, gamma = gamma, 
                          scale = scale, printWarn = FALSE, start = start, mustart = mustart, 
                          scoreType = scoreType, null.coef = null.coef, d_null.coef = d_null.coef,
                          pearson.extra = pearson.extra, 
                          dev.extra = dev.extra, n.true = n.true, Sl = Sl, mgcv_version = FALSE, ...)
        
        deviance_gradient <- b1$deviance_gradient
        deviance_hessian <- b1$deviance_hessian
        smoothing_parameters <- c(smoothing_parameters,L %*% lsp1 + lsp0)
        
        pred.change <- sum(grad * step) + 0.5 * t(step) %*% 
          hess %*% step
        if (reml) {
          score1 <- b1$REML
          d_score1 <- b1$d_REML
          d2_score1 <- b1$d2_REML
        }
        else if (scoreType == "GACV") {
          score1 <- b1$GACV
        }
        else if (scoreType == "UBRE") {
          score1 <- b1$UBRE
        }
        else score1 <- b1$GCV
        score.change <- score1 - score
        qerror <- abs(pred.change - score.change)/(max(abs(pred.change), 
                                                       abs(score.change)) + score.scale * conv.tol)
        if (is.finite(score1) && score.change < 0 && 
            qerror < qerror.thresh) { 
          
          print("is.finite(score1) && score.change < 0 && qerror < qerror.thresh")
          
          if (pdef || !sd.unused) { 
            
            print("pdef || !sd.unused")
            
            d2_lsp1 <- sapply(d2_lsp1, function(x) { x[lower.tri(x,diag = FALSE)] <- t(x)[lower.tri(x,diag = FALSE)] 
            x },simplify = FALSE) #assign lower triangular to upper
            
            L0 <- diag(length(d_lsp[[1]]))
            
            d_lsp1_in <- lapply(d_lsp1, function(x) L0 %*% x)
            d2_lsp1_in <- lapply(d2_lsp1, function(x) L0 %*% x)
            
            d_lsp_out[[length(d_lsp_out) + 1]] <- d_lsp1_in
            d2_lsp_out[[length(d2_lsp_out) + 1]] <- d2_lsp1_in
            
            print("Fourth")
            
            b <- gam.fit3.si(init_mod_pen = modPen_out, coef_derivs = coefDerivs_out, beta_star = beta_star, splineRank = splineRank, rS = rS, d_rS = d_rS, d2_rS = d2_rS,
                             x = X, d_X = d_X, d2_X = d2_X, X_uncons=X_uncons, d_X_cons = d_X_cons, d2_X_cons = d2_X_cons, Z = Z, d_Z = d_Z, d2_Z = d2_Z,
                             y = y, sp = L %*% lsp1 + lsp0, d_sp = d_lsp1_in, d2_sp = d2_lsp1_in, Eb = Eb, 
                             UrS = UrS, d_UrS = d_UrS, d2_UrS = d2_UrS,
                             offset = offset, U1 = U1, d_U1 = d_U1, d2_U1 = d2_U1,
                             Mp = Mp, family = family, 
                             weights = weights, deriv = 2, control = control, gamma = gamma, 
                             scale = scale, printWarn = FALSE, start = start, mustart = mustart, 
                             scoreType = scoreType, null.coef = null.coef, d_null.coef = d_null.coef,
                             pearson.extra = pearson.extra, 
                             dev.extra = dev.extra, n.true = n.true, Sl = Sl, ...)
            
            deviance_gradient <- b$deviance_gradient
            deviance_hessian <- b$deviance_hessian
            smoothing_parameters <- c(smoothing_parameters,L %*% lsp1 + lsp0)
            
            mustart <- b$fitted.values
            etastart <- b$linear.predictors
            start <- b$coefficients
            coefDerivs_out <- b$coef_derivs
            
            old.score <- score
            d_old.score <- d_score
            d2_old.score <- d2_score
            
            lsp <- lsp1
            d_lsp <- d_lsp1
            d2_lsp <- d2_lsp1
            if (reml) {
              score <- b$REML
              grad <- b$REML1
              hess <- b$REML2
              
              d_score <- b$d_REML
              d_grad <- b$d_REML1
              d_hess <- b$d_REML2
              
              d2_score <- b$d2_REML
              d2_grad <- b$d2_REML1
              d2_hess <- b$d2_REML2
              
              d_hess <- lapply(d_hess, function(x) {x[lower.tri(x,diag = FALSE)] <- t(x)[lower.tri(x,diag = FALSE)]
              x})
              d2_hess <- lapply(d2_hess,function(x){
                out <- lapply(x,function(y){ 
                  if(!is.null(y)) { y[lower.tri(y,diag = FALSE)] <- t(y)[lower.tri(y,diag = FALSE)] } 
                  y })
                out })
            }
            else if (scoreType == "GACV") {
              score <- b$GACV
              grad <- b$GACV1
              hess <- b$GACV2
            }
            else if (scoreType == "UBRE") {
              score <- b$UBRE
              grad <- b$UBRE1
              hess <- b$UBRE2
            }
            else {
              score <- b$GCV
              grad <- b$GCV1
              hess <- b$GCV2
            }
            grad <- t(L) %*% grad
            hess <- t(L) %*% hess %*% L
            
            d_grad <- lapply(d_grad, function(x) crossprod(L,x) )
            d_hess <- lapply(d_hess, function(x) crossprod(L,x) %*% L)
            
            d2_grad <- lapply(d2_grad, function(x) lapply(x, function(y){if(!is.null(y)) crossprod(L,y)}) )
            d2_hess <- lapply(d2_hess, function(x) lapply(x, function(y){if(!is.null(y)) crossprod(L,y) %*% L }) )
            
            if (!is.null(lsp.max)) {
              delta <- delta1
              rho <- rt(delta, lsp1.max)
              nr <- length(rho$rho1)
              hess <- diag(rho$rho1, nr, nr) %*% hess %*% 
                diag(rho$rho1, nr, nr) + diag(rho$rho2 * 
                                                grad)
              grad <- rho$rho1 * grad
            }
          }
          else {
            b <- b1
            score2 <- score1
            d_score2 <- d_score1
            d2_score2 <- d2_score1
          }
          
          score1 <- score - abs(score) - 1
          d_score1 <- d_score - d_score*score/abs(score)
          
          for(d1 in 1:length(beta_star))
          {
            for(d2 in d1:length(beta_star))
            {
              d2_score1[d1,d2] <- d2_score[d1,d2] - 
                ( (d2_score[d1,d2]*score + d_score[d1]*d_score[d2])/abs(score) - 
                    (d_score[d1]*score * (d_score[d2]*score/abs(score)))/(abs(score)^2)   )
            }
          }
        }
        if (!is.finite(score1) || score1 >= score || 
            qerror >= qerror.thresh) 
          ii <- ii + 1
        
      } #The end of the while
      if (!pdef && sd.unused && ii < maxHalf){ 
        print("!pdef && sd.unused && ii < maxHalf")
        score1 <- score2
        d_score1 <- d_score2
        d2_score1 <- d2_score2
      }
  }
    if (!pdef && sd.unused) { 
      print("!pdef && sd.unused")
      
      step <- Sstep * 2
      d_step <- lapply(d_Sstep, function(x){out <- x*2;out})
      
      d2_step <- list()
      for(d1 in 1:length(beta_star))
      {
        d2_step[[d1]] <- list()
        for(d2 in d1:length(beta_star))
        {
          d2_step[[d1]][[d2]] <- d2_Sstep[[d1]][[d2]] * 2
        }
      }
      
      kk <- 0
      score2 <- NA
      d_score2 <- lapply(score2,function(x){ x <- NA ; x})
      
      d2_score2 <- list()
      for(d1 in 1:length(beta_star))
      {
        d2_score2[[d1]] <- list()
        for(d2 in d1:length(beta_star)){
          d2_score2[[d1]][[d2]] <- NA
        }
      }
      
      ok <- TRUE
      while (ok) { print("while (ok)")
        step <- step/2
        
        d_step <- lapply(d_step,function(x){ out <- x/2; out})
        
        for(d1 in 1:length(beta_star))
        {
          for(d2 in d1:length(beta_star))
          {
            d2_step[[d1]][[d2]] <- d2_step[[d1]][[d2]]/2
          }
        }
        
        kk <- kk + 1
        if (!is.null(lsp.max)) {
          delta3 <- delta + step
          lsp3 <- rt(delta3, lsp1.max)$rho
        } else{
          d2_lsp3 <- list()
          for(i in 1:length(lsp))
          {
            d2_lsp3[[i]] <- matrix(nrow=length(beta_star),ncol=length(beta_star))
            for(d1 in 1:length(beta_star))
            {
              for(d2 in d1:length(beta_star))
              {
                
                d2_lsp3[[i]][d1,d2] <-
                  d2_lsp[[i]][d1,d2]*exp(step[i]) + 
                  d_lsp[[i]][d1]*d_step[[d2]][i]*exp(step[i]) + 
                  d_lsp[[i]][d2]*d_step[[d1]][i]*exp(step[i]) + 
                  exp(lsp[i])*d2_step[[d1]][[d2]][i]*exp(step[i]) + 
                  exp(lsp[i])*d_step[[d1]][i]*d_step[[d2]][i]*exp(step[i])
              }
            }
          }
          
          d_lsp3 <- lapply(1:length(lsp), function(i){
            out <- unlist(lapply(1:length(beta_star), function(d1){
              out2 <- d_lsp[[i]][d1]*exp(step[i]) + exp(lsp[i])*d_step[[d1]][i]*exp(step[i])
              out2
            }))
            out
          })
          
          lsp3 <- lsp + step
        }
        
        d2_lsp3 <- sapply(d2_lsp3, function(x) { x[lower.tri(x,diag = FALSE)] <- t(x)[lower.tri(x,diag = FALSE)] 
        x },simplify = FALSE) #assign lower triangular to upper
        
        L0 <- diag(length(d_lsp[[1]]))
        
        d_lsp3_in <- lapply(d_lsp3, function(x) L0 %*% x)
        d2_lsp3_in <- lapply(d2_lsp3, function(x) L0 %*% x)
        
        print("Fifth")
        
        d_lsp_out[[length(d_lsp_out) + 1]] <- d_lsp3_in
        d2_lsp_out[[length(d2_lsp_out) + 1]] <- d2_lsp3_in
        
        #Updated deriv inputs
        b1 <- gam.fit3.si(init_mod_pen = modPen_out, coef_derivs = coefDerivs_out, beta_star = beta_star, splineRank = splineRank, rS = rS, d_rS = d_rS, d2_rS = d2_rS,
                          x = X, d_X = d_X, d2_X = d2_X, X_uncons=X_uncons, d_X_cons = d_X_cons, d2_X_cons = d2_X_cons, Z = Z, d_Z = d_Z, d2_Z = d2_Z,
                          y = y, sp = L %*% lsp3 + lsp0, d_sp = d_lsp3_in, d2_sp = d2_lsp3_in, Eb = Eb, 
                          UrS = UrS, d_UrS = d_UrS, d2_UrS = d2_UrS,
                          offset = offset, U1 = U1, d_U1 = d_U1, d2_U1 = d2_U1,
                          Mp = Mp, family = family, 
                          weights = weights, deriv = 0, control = control, gamma = gamma, 
                          scale = scale, printWarn = FALSE, start = start, mustart = mustart, 
                          scoreType = scoreType, null.coef = null.coef, d_null.coef = d_null.coef,
                          pearson.extra = pearson.extra, 
                          dev.extra = dev.extra, n.true = n.true, Sl = Sl, mgcv_version = FALSE, ...)
        
        deviance_gradient <- b1$deviance_gradient
        deviance_hessian <- b1$deviance_hessian
        smoothing_parameters <- c(smoothing_parameters,L %*% lsp3 + lsp0)
        
        pred.change <- sum(grad * step) + 0.5 * t(step) %*% hess %*% step
        if (reml) {
          score3 <- b1$REML
          d_score3 <- b1$d_REML
          d2_score3 <- b1$d2_REML
        } else if (scoreType == "GACV") {
          score3 <- b1$GACV
        } else if (scoreType == "UBRE") {
          score3 <- b1$UBRE
        } else score3 <- b1$GCV
        score.change <- score3 - score
        qerror <- abs(pred.change - score.change)/(max(abs(pred.change), abs(score.change)) + score.scale * conv.tol)
        if (!is.finite(score2) || (is.finite(score3) && score3 <= score2 && qerror < qerror.thresh)) {
          
          print("!is.finite(score2) || (is.finite(score3) && score3 <= score2 && qerror < qerror.thresh)")
          
          score2 <- score3
          d_score2 <- d_score3
          d2_score2 <- d2_score3
          
          lsp2 <- lsp3
          d_lsp2 <- d_lsp3
          d2_lsp2 <- d2_lsp3
          
          if (!is.null(lsp.max)) 
            delta2 <- delta3
        }
        if ((is.finite(score2) && is.finite(score3) && score2 < score && score3 > score2) || kk == 40){ 
          print("(is.finite(score2) && is.finite(score3) && score2 < score && score3 > score2) || kk == 40")
          ok <- FALSE
        }
      }
      
      if (is.finite(score2) && score2 < score1) { 
        
        print("is.finite(score2) && score2 < score1")
        
        lsp1 <- lsp2
        d_lsp1 <- d_lsp2
        d2_lsp1 <- d2_lsp2
        
        if (!is.null(lsp.max)) 
          delta1 <- delta2
        
        score1 <- score2
        d_score1 <- d_score2
        d2_score1 <- d2_score2
      }
      
      d2_lsp1 <- sapply(d2_lsp1, function(x) { x[lower.tri(x,diag = FALSE)] <- t(x)[lower.tri(x,diag = FALSE)] 
      x },simplify = FALSE) #assign lower triangular to upper
      
      L0 <- diag(length(d_lsp[[1]]))
      
      d_lsp1_in <- lapply(d_lsp1, function(x) L0 %*% x)
      d2_lsp1_in <- lapply(d2_lsp1, function(x) L0 %*% x)
      
      print("Sixth")
      
      d_lsp_out[[length(d_lsp_out) + 1]] <- d_lsp1_in
      d2_lsp_out[[length(d2_lsp_out) + 1]] <- d2_lsp1_in
      
      b <- gam.fit3.si(init_mod_pen = modPen_out, coef_derivs = coefDerivs_out, beta_star = beta_star, splineRank = splineRank, rS = rS, d_rS = d_rS, d2_rS = d2_rS,
                       x = X, d_X = d_X, d2_X = d2_X, X_uncons=X_uncons, d_X_cons = d_X_cons, d2_X_cons = d2_X_cons, Z = Z, d_Z = d_Z, d2_Z = d2_Z,
                       y = y, sp = L %*% lsp1 + lsp0, d_sp = d_lsp1_in, d2_sp = d2_lsp1_in, Eb = Eb, 
                       UrS = UrS, d_UrS = d_UrS, d2_UrS = d2_UrS,
                       offset = offset, U1 = U1, d_U1 = d_U1, d2_U1 = d2_U1,
                       Mp = Mp, family = family, 
                       weights = weights, deriv = 2, control = control, gamma = gamma, 
                       scale = scale, printWarn = FALSE, start = start, mustart = mustart, 
                       scoreType = scoreType, null.coef = null.coef, d_null.coef = d_null.coef,
                       pearson.extra = pearson.extra, 
                       dev.extra = dev.extra, n.true = n.true, Sl = Sl, mgcv_version = FALSE, ...)
      
      deviance_gradient <- b$deviance_gradient
      deviance_hessian <- b$deviance_hessian
      smoothing_parameters <- c(smoothing_parameters,L %*% lsp1 + lsp0)
      
      mustart <- b$fitted.values
      etastart <- b$linear.predictors
      start <- b$coefficients
      coefDerivs_out <- b$coef_derivs
      modPen_out <- b$mod_pen
      
      old.score <- score
      d_old.score <- d_score
      d2_old.score <- d2_score
      
      lsp <- lsp1
      d_lsp <- d_lsp1
      d2_lsp <- d2_lsp1
      
      if (reml) {
        score <- b$REML
        grad <- b$REML1
        hess <- b$REML2
        
        d_score <- b$d_REML
        d_grad <- b$d_REML1
        d_hess <- b$d_REML2
        
        d2_score <- b$d2_REML
        d2_grad <- b$d2_REML1
        d2_hess <- b$d2_REML2
        
        d_hess <- lapply(d_hess, function(x) {x[lower.tri(x,diag = FALSE)] <- t(x)[lower.tri(x,diag = FALSE)]
        x})
        
        d2_hess <- lapply(d2_hess,function(x){
          out <- lapply(x,function(y){ 
            if(!is.null(y)) { y[lower.tri(y,diag = FALSE)] <- t(y)[lower.tri(y,diag = FALSE)] } 
            y })
          out })
        
      } else if (scoreType == "GACV") {
        score <- b$GACV
        grad <- b$GACV1
        hess <- b$GACV2
      } else if (scoreType == "UBRE") {
        score <- b$UBRE
        grad <- b$UBRE1
        hess <- b$UBRE2
      } else {
        score <- b$GCV
        grad <- b$GCV1
        hess <- b$GCV2
      }
      grad <- t(L) %*% grad
      hess <- t(L) %*% hess %*% L
      
      d_grad <- lapply(d_grad, function(x) crossprod(L,x) )
      d_hess <- lapply(d_hess, function(x) crossprod(L,x) %*% L)
      
      d2_grad <- lapply(d2_grad, function(x) lapply(x, function(y){if(!is.null(y)) crossprod(L,y)}) )
      d2_hess <- lapply(d2_hess, function(x) lapply(x, function(y){if(!is.null(y)) crossprod(L,y) %*% L }) )
      
      if (!is.null(lsp.max)) {
        delta <- delta1
        rho <- rt(delta, lsp1.max)
        nr <- length(rho$rho1)
        hess <- diag(rho$rho1, nr, nr) %*% hess %*% diag(rho$rho1, 
                                                         nr, nr) + diag(rho$rho2 * grad)
        grad <- rho$rho1 * grad
      }
    }
    score.hist[i] <- score
    converged <- TRUE
    if (reml){ 
      score.scale <- abs(log(b$scale.est)) + abs(score)
    } else score.scale <- abs(b$scale.est) + abs(score)
    grad2 <- diag(hess)
    uconv.ind <- (abs(grad) > score.scale * conv.tol * 0.1) | 
      (abs(grad2) > score.scale * conv.tol * 0.1)
    if (sum(abs(grad) > score.scale * conv.tol * 5)) 
      converged <- FALSE
    if (abs(old.score - score) > score.scale * conv.tol) {
      if (converged) 
        uconv.ind <- uconv.ind | TRUE
      converged <- FALSE
    }
    if (ii == maxHalf) 
      converged <- TRUE
    if (converged) 
      break
  }
  if (ii == maxHalf) {
    ct <- "step failed"
    warning("Fitting terminated with step failure - check results carefully")
  }
  else if (i == 200) {
    ct <- "iteration limit reached"
    warning("Iteration limit reached without full convergence - check carefully")
  }
  else ct <- "full convergence"
  b$dVkk <- NULL
  if (as.logical(edge.correct) && reml) { print("as.logical(edge.correct) && reml") 
    stop("Not yet integrated")
    flat <- which(abs(grad2) < abs(grad) * 100)
    REML <- b$REML
    alpha <- if (is.logical(edge.correct)) 
      0.02
    else abs(edge.correct)
    b1 <- b
    lsp1 <- lsp
    if (length(flat)) 
      for (i in flat) {
        REML <- b1$REML + alpha
        while (b1$REML < REML) {
          lsp1[i] <- lsp1[i] - 1
          print("Seventh")
          b1 <- gam.fit3(x = X, y = y, sp = L %*% lsp1 + 
                           lsp0, Eb = Eb, UrS = UrS, offset = offset, 
                         U1 = U1, Mp = Mp, family = family, weights = weights, 
                         deriv = 0, control = control, gamma = gamma, 
                         scale = scale, printWarn = FALSE, start = start, 
                         mustart = mustart, scoreType = scoreType, 
                         null.coef = null.coef, pearson.extra = pearson.extra, 
                         dev.extra = dev.extra, n.true = n.true, Sl = Sl, 
                         ...)
        }
      }
    print("Eigth")
    b1 <- gam.fit3(x = X, y = y, sp = L %*% lsp1 + lsp0, 
                   Eb = Eb, UrS = UrS, offset = offset, U1 = U1, Mp = Mp, 
                   family = family, weights = weights, deriv = 2, control = control, 
                   gamma = gamma, scale = scale, printWarn = FALSE, 
                   start = start, mustart = mustart, scoreType = scoreType, 
                   null.coef = null.coef, pearson.extra = pearson.extra, 
                   dev.extra = dev.extra, n.true = n.true, Sl = Sl, 
                   ...)
    score1 <- b1$REML
    grad1 <- b1$REML1
    hess1 <- b1$REML2
    grad1 <- t(L) %*% grad1
    hess1 <- t(L) %*% hess1 %*% L
    if (!is.null(lsp.max)) {
      delta <- delta1
      rho <- rt(delta, lsp1.max)
      nr <- length(rho$rho1)
      hess1 <- diag(rho$rho1, nr, nr) %*% hess1 %*% diag(rho$rho1, 
                                                         nr, nr) + diag(rho$rho2 * grad1)
      grad1 <- rho$rho1 * grad1
    }
    attr(hess, "edge.correct") <- TRUE
    attr(hess, "hess1") <- hess1
    attr(hess, "db.drho1") <- b1$db.drho
    attr(hess, "dw.drho1") <- b1$dw.drho
    attr(hess, "lsp1") <- lsp1
    attr(hess, "rp") <- b1$rp
  }
  
  b$Gradient <- deviance_gradient
  b$Hessian <- deviance_hessian
  
  b$d_lsp_out <- d_lsp_out
  b$d2_lsp_out <- d2_lsp_out
  
  list(score = score, lsp = lsp, lsp.full = L %*% lsp + lsp0, 
       grad = grad, hess = hess, iter = i, conv = ct, score.hist = score.hist[!is.na(score.hist)], 
       object = b)
}

#08
gam.fit3.si <- function (init_mod_pen = NULL, coef_derivs = NULL, beta_star, splineRank, rS, d_rS, d2_rS,
                         x, d_X, d2_X, X_uncons, d_X_cons, d2_X_cons ,Z, d_Z, d2_Z, y, sp, d_sp, d2_sp, Eb, UrS = list(), d_UrS = list(), d2_UrS = list(),
                         weights = rep(1, nobs), 
                         start = NULL, etastart = NULL, mustart = NULL, offset = rep(0, 
                                                                                     nobs), U1 = diag(ncol(x)), d_U1, d2_U1,
                         Mp = -1, family = gaussian(), 
                         control = gam.control(), intercept = TRUE, deriv = 2, gamma = 1, 
                         scale = 1, printWarn = TRUE, scoreType = "REML", null.coef = rep(0, ncol(x)), d_null.coef = NULL,
                         pearson.extra = 0, dev.extra = 0, n.true = -1, 
                         Sl = NULL, mgcv_version = TRUE,...) { 
  XZ <- x
  rS0 <- rS
  if (control$trace) {
    t0 <- proc.time()
    tc <- 0
  }
  if (inherits(family, "extended.family")) {
    if (inherits(family, "general.family")) {
      return(gam.fit5(x, y, sp, Sl = Sl, weights = weights, 
                      offset = offset, deriv = deriv, family = family, 
                      control = control, Mp = Mp, start = start))
    } else return(gam.fit4(x, y, sp, Eb, UrS = UrS, weights = weights, 
                           start = start, etastart = etastart, mustart = mustart, 
                           offset = offset, U1 = U1, Mp = Mp, family = family, 
                           control = control, deriv = deriv, scale = scale, 
                           scoreType = scoreType, null.coef = null.coef, ...))
  }
  if (family$link == family$canonical){ 
    fisher <- TRUE #We use canonical link so fisher is true
  } else{ 
    fisher = FALSE
  }
  if (scale > 0){ 
    scale.known <- TRUE #scale.known is set to true
  } else{ 
    scale.known <- FALSE
  }
  if (!scale.known && scoreType %in% c("REML", "ML")) { 
    nsp <- length(sp)
    scale <- exp(sp[nsp])
    d_scale <- d_sp[[nsp]]
    d2_scale <- d2_sp[[nsp]]
    sp <- sp[-nsp]
    d_sp[[nsp]] <- d2_sp[[nsp]] <- NULL
  } else if (scale.known && scoreType %in% c("REML", "ML")) {
    d_scale <- 0*d_sp[[1]]
    d2_scale <- 0*d2_sp[[1]]  
  }
  if (!deriv %in% c(0, 1, 2)) 
    stop("unsupported order of differentiation requested of gam.fit3")
  x <- as.matrix(x)
  nSp <- length(sp)
  if (nSp == 0){ 
    deriv.sp <- 0
  } else{
    deriv.sp <- deriv 
  }
  rank.tol <- .Machine$double.eps * 100
  xnames <- dimnames(x)[[2]]
  ynames <- if (is.matrix(y)){ 
    rownames(y)
  } else{ 
    names(y)
  }
  q <- ncol(x)
  if(length(UrS)){
    grderiv <- deriv * as.numeric(scoreType %in% c("REML", "ML", "P-REML", "P-ML"))
    
    if(mgcv_version){
      rp <- mgcv:::gam.reparam(UrS, sp, grderiv) 
    } else{
      rp <- gam.reparamR.gam(UrS, sp, grderiv)
    }
    
    T <- diag(q)
    T[1:ncol(rp$Qs), 1:ncol(rp$Qs)] <- rp$Qs 
    T <- U1 %*% T
    
    if( identical( rp$Qs-diag(ncol(rp$Qs)), matrix(0,nrow=nrow(rp$Qs),ncol= ncol(rp$Qs)) ) ){
      d_T <- d_U1
      d2_T <- d2_U1
    } else{
      print("Non-identity Qs")
      
      dpen <- d_Pen(beta_star=beta_star, splineRank = splineRank, rS = UrS, d_rS = d_UrS, 
                    d2_rS =d2_UrS, T_mat = T, d_T=d_U1, d2_T=d2_U1, UrS = UrS) #These have been transformed with U1
      
      Si <- dpen$Penalty
      d_Si <- dpen$d_Penalty
      d2_Si <- dpen$d2_Penalty
      
      inputs <<- list(rS = UrS, d_rS = d_UrS, d2_rS = d2_UrS, lsp = sp, deriv = grderiv, beta_star = beta_star, 
                      d_lsp = d_sp, d2_lsp = d2_sp, Si = Si, d_Si = d_Si, d2_Si = d2_Si, dpen = dpen)
      
      d_rp <- gam.reparamR.si(rS = UrS, d_rS = d_UrS, d2_rS = d2_UrS, lsp = sp, deriv = grderiv, beta_star = beta_star, 
                              d_lsp = d_sp, d2_lsp = d2_sp, Si = Si, d_Si = d_Si, d2_Si = d2_Si, dpen = dpen)
      
      Qs <- diag(q)
      Qs[1:ncol(rp$Qs), 1:ncol(rp$Qs)] <- rp$Qs
      
      d_Qs <- d2_Qs <- list()
      for(d1 in 1:length(beta_star))
      {
        d2_Qs[[d1]] <- list()
        d_Qs[[d1]] <- matrix(0,nrow = q, ncol = q)
        d_Qs[[d1]][1:ncol(rp$Qs), 1:ncol(rp$Qs)] <- d_rp$d_Qs[[d1]]
        for(d2 in d1:length(beta_star))
        {
          d2_Qs[[d1]][[d2]] <- matrix(0,nrow = q, ncol = q)
          d2_Qs[[d1]][[d2]][1:ncol(rp$Qs), 1:ncol(rp$Qs)] <- d_rp$d2_Qs[[d1]][[d2]]
        }
      }
      
      d_T <- d2_T <- list()
      for(d1 in 1:length(beta_star))
      {
        d2_T[[d1]] <- list()
        for(d2 in d1:length(beta_star))
        {
          d2_T[[d1]][[d2]] <- d2_U1[[d1]][[d2]] %*% Qs + d_U1[[d1]] %*% d_Qs[[d2]] + 
            d_U1[[d2]] %*% d_Qs[[d1]] + U1 %*% d2_Qs[[d1]][[d2]]
        }
        d_T[[d1]] <- d_U1[[d1]] %*% Qs + U1 %*% d_Qs[[d1]]
      }
    }
    
    #Set up model matrix x, and smoothing matrix
    null.coef <- t(T) %*% null.coef 
    if(!is.null(start)){ 
      start <- t(T) %*% start 
    }
    x <- .Call(mgcv:::C_mgcv_pmmult2, x, T, 0, 0, control$nthreads) 
    rS <- list()
    for (i in 1:length(UrS)) 
    {
      rS[[i]] <- rbind(rp$rS[[i]], matrix(0, Mp, ncol(rp$rS[[i]])))  
    }
    
    nrowS <- nrow(rp$S)
    Eb <- Eb %*% T
    rows.E <- q - Mp
    Sr <- cbind(rp$E, matrix(0, nrow(rp$E), Mp))
    St <- rbind(cbind(rp$S, matrix(0, nrow(rp$S), Mp)), matrix(0,Mp, q))
    
    if(is.null(init_mod_pen)){
      #Set this up assuming identity Qs
      d_modmat <- d_mm(intercept=intercept, beta_star=beta_star, X=XZ, d_X=d_X_cons, d2_X=d2_X_cons, Z=Z, d_Z=d_Z, d2_Z=d2_Z, T_mat = U1, d_T=d_U1, d2_T=d2_U1, UrS=UrS) 
      dpen <- d_Pen(beta_star=beta_star, splineRank = splineRank, rS = UrS, d_rS = d_UrS, d2_rS =d2_UrS, T_mat = U1, d_T=d_U1, d2_T=d2_U1, UrS = UrS)
      
      Si <- dpen$Penalty
      d_Si <- dpen$d_Penalty
      d2_Si <- dpen$d2_Penalty
      mod_pen <- list(d_modmat=d_modmat,dpen=dpen)
      
      #Initialize the derivatives
      if(is.null(coef_derivs)){
        drho_coef <- d_null.coef$drho_coefstart 
        d2rho_coef <- d_null.coef$d2rho_coefstart
        d_coef <- d_null.coef$d_coefstart
        drho_d_coef <- d_null.coef$drho_d_coefstart
        d2rho_d_coef <- d_null.coef$d2rho_d_coefstart
        d2_coef <- d_null.coef$d2_coefstart
        drho_d2_coef <- d_null.coef$drho_d2_coefstart
        d2rho_d2_coef <- d_null.coef$d2rho_d2_coefstart
      }
      
      #The model matrix in this case has not been transformed with U1 so we don't need to worry about different cases of identity/non-identy Qs
    } else{
      mod_pen <- init_mod_pen
      
      #Pull out penalties
      Si <- mod_pen$dpen$Penalty
      d_Si <- mod_pen$dpen$d_Penalty
      d2_Si <- mod_pen$dpen$d2_Penalty
      
      #Pull out model matrix
      d_modmat <- mod_pen$d_modmat
      
      d_x <- d_modmat$d_X
      d2_x <- d_modmat$d2_X
      
      #Pull out the coefficients
      drho_coef <- coef_derivs$drho_coef 
      d2rho_coef <- coef_derivs$d2rho_coef
      d_coef <- coef_derivs$d_coef
      drho_d_coef <- coef_derivs$drho_d_coef
      d2rho_d_coef <- coef_derivs$d2rho_d_coef
      d2_coef <- coef_derivs$d2_coef
      drho_d2_coef <- coef_derivs$drho_d2_coef
      d2rho_d2_coef <- coef_derivs$d2rho_d2_coef
    }
    
    #Transform model matrix and smoothing matrix for non-identity Qs
    if(identical( rp$Qs-diag(ncol(rp$Qs)) , matrix(0,nrow=nrow(rp$Qs),ncol= ncol(rp$Qs)) )){
      
      #Pull out the d_X, d2_X matrices:
      if(!is.null(init_mod_pen)){
        mod_pen <- init_mod_pen 
      } 
      
      d_x <- mod_pen$d_modmat$d_X
      d2_x <- mod_pen$d_modmat$d2_X
      
      drho_S <- d2rho_S <- drho_d_S <- drho_d2_S <- d2rho_d_S <- d2rho_d2_S <- list()  
      for(d1 in 1:length(beta_star))
      {
        drho_d_S[[d1]] <- drho_d2_S[[d1]] <- d2rho_d_S[[d1]] <- d2rho_d2_S[[d1]] <- list()
        for(d2 in d1:length(beta_star))
        {
          drho_d2_S[[d1]][[d2]] <- d2rho_d2_S[[d1]][[d2]] <- list()
          for(i in 1:length(sp))
          {
            drho_S[[i]] <- exp(sp[i])*mod_pen$dpen$Penalty[[i]]
            
            drho_d_S[[d1]][[i]] <- d_sp[[i]][d1] * mod_pen$dpen$Penalty[[i]] + exp(sp[i]) * mod_pen$dpen$d_Penalty[[i]][[d1]]
            
            drho_d2_S[[d1]][[d2]][[i]] <- d2_sp[[i]][d1,d2] * mod_pen$dpen$Penalty[[i]] + d_sp[[i]][d1] * mod_pen$dpen$d_Penalty[[i]][[d2]] + 
              d_sp[[i]][d2] * mod_pen$dpen$d_Penalty[[i]][[d1]] + exp(sp[i]) * mod_pen$dpen$d2_Penalty[[i]][[d1]][[d2]]
            
            d2rho_S[[i]] <- d2rho_d_S[[d1]][[i]] <- d2rho_d2_S[[d1]][[d2]][[i]] <- list()
            for(j in i:length(sp))
            {
              d2rho_S[[i]][[j]] <- if(i==j)
              {
                exp(sp[i])*mod_pen$dpen$Penalty[[i]]
              }
              else
              {
                0*mod_pen$dpen$Penalty[[i]]
              }
              
              d2rho_d_S[[d1]][[i]][[j]] <- if(i==j)
              {
                d_sp[[i]][d1] * mod_pen$dpen$Penalty[[i]] + exp(sp[i]) * mod_pen$dpen$d_Penalty[[i]][[d1]]
              }
              else
              {
                0*mod_pen$dpen$Penalty[[i]]
              }
              
              d2rho_d2_S[[d1]][[d2]][[i]][[j]] <- if(i==j)
              {
                d2_sp[[i]][d1,d2] * mod_pen$dpen$Penalty[[i]] + d_sp[[i]][d1] * mod_pen$dpen$d_Penalty[[i]][[d2]] + 
                  d_sp[[i]][d2] * mod_pen$dpen$d_Penalty[[i]][[d1]] + exp(sp[i]) * mod_pen$dpen$d2_Penalty[[i]][[d1]][[d2]]
              }
              else
              {
                0*mod_pen$dpen$Penalty[[i]]
              }
            }
          }
        }
      }
      
      d_S <- d2_S <- list()
      for(d1 in 1:length(beta_star))
      {
        d_S[[d1]] <- matrix(0,nrow=nrow(Si[[1]]),ncol=ncol(Si[[1]]))
        d2_S[[d1]] <- list()
        for(d2 in d1:length(beta_star))
        {
          d2_S[[d1]][[d2]] <- matrix(0,nrow=nrow(Si[[1]]),ncol=ncol(Si[[1]]))
          
          for(i in 1:length(Si)) #This probably doesn't need a loop as it's zero for all apart from the first
          {
            d2_S[[d1]][[d2]] <- d2_S[[d1]][[d2]] + 
              
              d_sp[[i]][d2]*d_Si[[i]][[d1]][1:nrow(Si[[i]]),1:ncol(Si[[i]])] + 
              
              exp(sp[i])*d2_Si[[i]][[d1]][[d2]][1:nrow(Si[[i]]),1:ncol(Si[[i]])] +
              
              d2_sp[[i]][d1,d2]*Si[[i]] + 
              
              d_sp[[i]][d1]*d_Si[[i]][[d2]][1:nrow(Si[[i]]),1:ncol(Si[[i]])]
          }
        }
        for(i in 1:length(Si))
        {
          d_S[[d1]] <- d_S[[d1]] + 
            exp(sp[i])*d_Si[[i]][[d1]][1:nrow(Si[[i]]),1:ncol(Si[[i]])] + 
            d_sp[[i]][d1]*Si[[i]]  
        }
      }
      
      for(i in 1:length(sp))
      {
        for(j in i:length(sp))
        {
          d2rho_S[[i]][[j]] <- rbind(cbind(d2rho_S[[i]][[j]], matrix(0, nrowS, Mp)), matrix(0,Mp, q))
        }
        drho_S[[i]] <- rbind(cbind(drho_S[[i]], matrix(0, nrowS, Mp)), matrix(0,Mp, q))
      }
      
      for(d1 in 1:length(beta_star))
      {
        for(d2 in d1:length(beta_star))
        {
          d2_S[[d1]][[d2]] <- rbind(cbind(d2_S[[d1]][[d2]], matrix(0, nrowS, Mp)), matrix(0,Mp, q))
          
          for(i in 1:length(sp))
          {
            for(j in i:length(sp))
            {
              d2rho_d2_S[[d1]][[d2]][[i]][[j]] <- rbind(cbind(d2rho_d2_S[[d1]][[d2]][[i]][[j]], matrix(0, nrowS, Mp)), matrix(0,Mp, q))
            }
            drho_d2_S[[d1]][[d2]][[i]] <- rbind(cbind(drho_d2_S[[d1]][[d2]][[i]], matrix(0, nrowS, Mp)), matrix(0,Mp, q))
          }
        }
        
        d_S[[d1]] <- rbind(cbind(d_S[[d1]], matrix(0, nrowS, Mp)), matrix(0,Mp, q))
        
        for(i in 1:length(sp))
        {
          for(j in i:length(sp))
          {
            d2rho_d_S[[d1]][[i]][[j]] <- rbind(cbind(d2rho_d_S[[d1]][[i]][[j]], matrix(0, nrowS, Mp)), matrix(0,Mp, q))
          }
          drho_d_S[[d1]][[i]] <- rbind(cbind(drho_d_S[[d1]][[i]], matrix(0, nrowS, Mp)), matrix(0,Mp, q))
        }
      }
      
    } else{
      #Transform with Qs
      mmQs <- d_x_Qs(beta_star = beta_star, x = XZ %*% U1, d_x = d_x, d2_x = d2_x, Qs = Qs, d_Qs = d_Qs, d2_Qs = d2_Qs)
      d_x <- mmQs$d_x
      d2_x <- mmQs$d2_x
      
      #Pull out the smoothing matrices:
      d_S <- d_rp$d_S
      d2_S <- d_rp$d2_S
      drho_S <- d_rp$drho_S
      d2rho_S <- d_rp$d2rho_S
      drho_d_S <- d_rp$drho_d_S
      drho_d2_S <- d_rp$drho_d2_S
      d2rho_d_S <- d_rp$d2rho_d_S
      d2rho_d2_S <- d_rp$d2rho_d2_S
      
      for(i in 1:length(sp))
      {
        for(j in i:length(sp))
        {
          d2rho_S[[i]][[j]] <- rbind(cbind(d2rho_S[[i]][[j]], matrix(0, nrowS, Mp)), matrix(0,Mp, q))
          
        }
        drho_S[[i]] <- rbind(cbind(drho_S[[i]], matrix(0, nrowS, Mp)), matrix(0,Mp, q))
      }
      
      for(d1 in 1:length(beta_star))
      {
        for(d2 in d1:length(beta_star))
        {
          d2_S[[d1]][[d2]] <- rbind(cbind(d2_S[[d1]][[d2]], matrix(0, nrowS, Mp)), matrix(0,Mp, q))
          
          for(i in 1:length(sp))
          {
            for(j in i:length(sp))
            {
              d2rho_d2_S[[d1]][[d2]][[i]][[j]] <- rbind(cbind(d2rho_d2_S[[d1]][[d2]][[i]][[j]], matrix(0, nrowS, Mp)), matrix(0,Mp, q))
            }
            drho_d2_S[[d1]][[d2]][[i]] <- rbind(cbind(drho_d2_S[[d1]][[d2]][[i]], matrix(0, nrowS, Mp)), matrix(0,Mp, q))
          }
        }
        
        d_S[[d1]] <- rbind(cbind(d_S[[d1]], matrix(0, nrowS, Mp)), matrix(0,Mp, q))
        
        for(i in 1:length(sp))
        {
          for(j in i:length(sp))
          {
            d2rho_d_S[[d1]][[i]][[j]] <- rbind(cbind(d2rho_d_S[[d1]][[i]][[j]], matrix(0, nrowS, Mp)), matrix(0,Mp, q))
          }
          drho_d_S[[d1]][[i]] <- rbind(cbind(drho_d_S[[d1]][[i]], matrix(0, nrowS, Mp)), matrix(0,Mp, q))
        }
      }
    }
    
    #Need to transform coefficients, this inlcudes the possibility of non-identity Qs
    
    coef_trans <- if(is.null(coef_derivs)){
      null.coef
    } else{
      coef_derivs$coef
    }
    
    trans <- d_coef_trans0(beta_star = beta_star, T_mat = T, d_T=d_T, d2_T=d2_T, coef = coef_trans, drho_coef = drho_coef, 
                           d2rho_coef = d2rho_coef, d_coef = d_coef, d2_coef = d2_coef, drho_d_coef = drho_d_coef,
                           drho_d2_coef = drho_d2_coef, d2rho_d_coef = d2rho_d_coef, d2rho_d2_coef = d2rho_d2_coef, sp=sp)
    
    drho_coef <- trans$drho_coef
    d2rho_coef <- trans$d2rho_coef
    d_coef <- trans$d_coef
    drho_d_coef <- trans$drho_d_coef
    d2rho_d_coef <- trans$d2rho_d_coef
    d2_coef <- trans$d2_coef
    drho_d2_coef <- trans$drho_d2_coef
    d2rho_d2_coef <- trans$d2rho_d2_coef
    
    ###End new version
  } else {
    T <- diag(q)
    St <- matrix(0, q, q)
    rSncol <- sp <- rows.E <- Eb <- Sr <- 0
    rS <- list(0)
    rp <- list(det = 0, det1 = rep(0, 0), det2 = rep(0, 0), 
               fixed.penalty = FALSE)
  }
  iter <- 0
  coef <- rep(0, ncol(x))
  conv <- FALSE
  n <- nobs <- NROW(y)
  if (n.true <= 0) 
    n.true <- nobs
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  if (is.null(weights)) weights <- rep.int(1, nobs)
  if (is.null(offset)) offset <- rep.int(0, nobs)
  variance <- family$variance 
  dev.resids <- family$dev.resids 
  aic <- family$aic
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  if (!is.function(variance) || !is.function(linkinv)) 
    stop("illegal `family' argument")
  valideta <- family$valideta
  if (is.null(valideta)) 
    valideta <- function(eta) TRUE
  validmu <- family$validmu
  if (is.null(validmu)) 
    validmu <- function(mu) TRUE
  if (is.null(mustart)) { #This is the initialization step for the coefficients
    eval(family$initialize) 
  } else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  if (family$family == "gaussian" && family$link == "identity"){
    strictly.additive <- TRUE
  } else{ 
    strictly.additive <- FALSE
  }
  D1 <- D2 <- P <- P1 <- P2 <- trA <- trA1 <- trA2 <- GCV <- GCV1 <- GCV2 <- GACV <- GACV1 <- GACV2 <- UBRE <- UBRE1 <- UBRE2 <- REML <- REML1 <- REML2 <- NULL
  d_UBRE <- d2_UBRE <- d_UBRE1 <- d2_UBRE1 <- d_UBRE2 <- d2_UBRE2 <- d_REML <- d2_REML <- d_REML1 <- d2_REML1 <- d_REML2 <- d2_REML2 <- NULL
  if (EMPTY) { 
    eta <- rep.int(0, nobs) + offset
    if (!valideta(eta)) 
      stop("Invalid linear predictor values in empty model")
    mu <- linkinv(eta)
    if (!validmu(mu)) 
      stop("Invalid fitted means in empty model")
    dev <- sum(dev.resids(y, mu, weights))
    w <- (weights * mu.eta(eta)^2)/variance(mu)
    residuals <- (y - mu)/mu.eta(eta)
    good <- rep(TRUE, length(residuals))
    boundary <- conv <- TRUE
    coef <- numeric(0)
    iter <- 0
    V <- variance(mu)
    alpha <- dev
    trA2 <- trA1 <- trA <- 0
    if (deriv) 
      GCV2 <- GCV1 <- UBRE2 <- UBRE1 <- trA1 <- rep(0, 
                                                    nSp)
    GCV <- nobs * alpha/(nobs - gamma * trA)^2
    UBRE <- alpha/nobs - scale + 2 * gamma/n * trA
    scale.est <- alpha/(nobs - trA)
  } else { #We enter the else part of this if statement on the initialisation - We don't have an empty model
    eta <- if (!is.null(etastart)){
      etastart
    } else if (!is.null(start)){
      if (length(start) != nvars){ 
        stop(gettextf("Length of start should equal %d and correspond to initial coefs for %s", 
                      nvars, deparse(xnames)))
      } else { #We don't get into this else
        coefold <- start
        offset + as.vector(if (NCOL(x) == 1){
          x * start
        } else{
          x %*% start
        })
      }
    } else family$linkfun(mustart)
    mu <- linkinv(eta)
    boundary <- conv <- FALSE
    rV = matrix(0, ncol(x), ncol(x))
    coefold <- null.coef #During this initialisation, at this point coefold is 
    
    etaold <- null.eta <- as.numeric(x %*% null.coef + as.numeric(offset))
    old.pdev <- sum(dev.resids(y, linkinv(null.eta), weights)) + 
      t(null.coef) %*% St %*% null.coef
    
    ii <- 0
    while (!(validmu(mu) && valideta(eta))) { # We do not seem to get into this while loop
      ii <- ii + 1
      if (ii > 20) 
        stop("Can't find valid starting values: please specify some")
      if (!is.null(start)) 
        start <- start * 0.9 + coefold * 0.1
      eta <- 0.9 * eta + 0.1 * etaold
      mu <- linkinv(eta)
    }
    for (iter in 1:control$maxit) {
      good <- weights > 0
      var.val <- variance(mu)
      varmu <- var.val[good]
      if (any(is.na(varmu))) 
        stop("NAs in V(mu)")
      if (any(varmu == 0)) 
        stop("0s in V(mu)")
      mu.eta.val <- mu.eta(eta)
      if (any(is.na(mu.eta.val[good]))) 
        stop("NAs in d(mu)/d(eta)")
      good <- (weights > 0) & (mu.eta.val != 0)
      if (all(!good)) {
        conv <- FALSE
        warning(gettextf("No observations informative at iteration %d", 
                         iter))
        break
      }
      mevg <- mu.eta.val[good]
      mug <- mu[good]
      yg <- y[good]
      weg <- weights[good]
      var.mug <- var.val[good]
      if (fisher) {
        z <- (eta - offset)[good] + (yg - mug)/mevg
        w <- (weg * mevg^2)/var.mug
      } else {
        c = yg - mug
        alpha <- 1 + c * (family$dvar(mug)/var.mug + 
                            family$d2link(mug) * mevg)
        alpha[alpha == 0] <- .Machine$double.eps
        z <- (eta - offset)[good] + (yg - mug)/(mevg * 
                                                  alpha)
        w <- weg * alpha * mevg^2/var.mug
      }
      if (sum(good) < ncol(x)) 
        stop("Not enough informative observations.")
      if (control$trace) 
        t1 <- proc.time()
      
      if(iter ==1)
      {
        if(!is.null(start)){
          coef <- start
        } else{
          coef <- coefold
        }
      }
      
      d_oo <- d_pls_fit1(deriv = deriv, family = family$family, beta_star = beta_star, good = good, X = x[good, ], y=yg, 
                         d_X = d_x, d2_X = d2_x, mu = mug, variance = w, pseudodata = z, 
                         sp=as.numeric(exp(sp)), d_sp=d_sp, d2_sp=d2_sp, St=St, d_St = d_S, d2_St = d2_S,
                         drho_S = drho_S, d2rho_S = d2rho_S, drho_d_S = drho_d_S, drho_d2_S = drho_d2_S, d2rho_d_S = d2rho_d_S, d2rho_d2_S = d2rho_d2_S,
                         coef = coef, drho_coef=drho_coef, d2rho_coef=d2rho_coef, d_coef=d_coef, drho_d_coef=drho_d_coef, 
                         d2rho_d_coef=d2rho_d_coef, d2_coef=d2_coef, drho_d2_coef=drho_d2_coef, d2rho_d2_coef=d2rho_d2_coef)
      
      drho_coef <- d_oo$drho_coef 
      d2rho_coef <- d_oo$d2rho_coef
      d_coef <- d_oo$d_coef
      drho_d_coef <- d_oo$drho_d_coef
      d2rho_d_coef <- d_oo$d2rho_d_coef
      d2_coef <- d_oo$d2_coef
      drho_d2_coef <- d_oo$drho_d2_coef
      d2rho_d2_coef <- d_oo$d2rho_d2_coef 
      
      oo <- .C(mgcv:::C_pls_fit1, y = as.double(z), X = as.double(x[good, 
                                                                    ]), w = as.double(w), wy = as.double(w * z), 
               E = as.double(Sr), Es = as.double(Eb), n = as.integer(sum(good)), 
               q = as.integer(ncol(x)), rE = as.integer(rows.E), 
               eta = as.double(z), penalty = as.double(1), 
               rank.tol = as.double(rank.tol), nt = as.integer(control$nthreads), 
               use.wy = as.integer(0))
      
      if (control$trace) 
        tc <- tc + sum((proc.time() - t1)[c(1, 4)])
      if (!fisher && oo$n < 0) { #We don't get into this if statement
        z <- (eta - offset)[good] + (yg - mug)/mevg
        w <- (weg * mevg^2)/var.mug
        if (control$trace) 
          t1 <- proc.time()
        oo <- .C(mgcv:::C_pls_fit1, y = as.double(z), X = as.double(x[good, 
                                                                      ]), w = as.double(w), E = as.double(Sr), Es = as.double(Eb), 
                 n = as.integer(sum(good)), q = as.integer(ncol(x)), 
                 rE = as.integer(rows.E), eta = as.double(z), 
                 penalty = as.double(1), rank.tol = as.double(rank.tol), 
                 nt = as.integer(control$nthreads))
        if (control$trace) 
          tc <- tc + sum((proc.time() - t1)[c(1, 4)])
      }
      start <- oo$y[1:ncol(x)] 
      penalty <- oo$penalty
      eta <- drop(x %*% start)
      if (any(!is.finite(start))) {
        conv <- FALSE
        warning(gettextf("Non-finite coefficients at iteration %d", 
                         iter))
        break
      }
      mu <- linkinv(eta <- eta + offset)
      dev <- sum(dev.resids(y, mu, weights))
      if (control$trace) 
        message(gettextf("Deviance = %s Iterations - %d", 
                         dev, iter, domain = "R-mgcv"))
      boundary <- FALSE
      if (!is.finite(dev)) {
        if (is.null(coefold)) {
          if (is.null(null.coef)) 
            stop("no valid set of coefficients has been found:please supply starting values", 
                 call. = FALSE)
          coefold <- null.coef
          etaold <- null.eta
        }
        warning("Step size truncated due to divergence", 
                call. = FALSE)
        ii <- 1
        while (!is.finite(dev)) { #We don't get into this while loop - This is when deviance is not finite so this may need to be addressed
          if (ii > control$maxit) 
            stop("inner loop 1; can't correct step size")
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- (eta + etaold)/2
          mu <- linkinv(eta)
          dev <- sum(dev.resids(y, mu, weights))
        }
        boundary <- TRUE
        penalty <- t(start) %*% St %*% start 
        if (control$trace) 
          cat("Step halved: new deviance =", dev, "\n")
      }
      if (!(valideta(eta) && validmu(mu))) { #Might need to use this if statement later if we don't have valid mu, valid eta as a failsafe procedure
        warning("Step size truncated: out of bounds", 
                call. = FALSE)
        ii <- 1
        while (!(valideta(eta) && validmu(mu))) { 
          if (ii > control$maxit) 
            stop("inner loop 2; can't correct step size")
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- (eta + etaold)/2
          mu <- linkinv(eta)
        }
        boundary <- TRUE
        penalty <- t(start) %*% St %*% start
        dev <- sum(dev.resids(y, mu, weights))
        if (control$trace) 
          cat("Step halved: new deviance =", dev, "\n")
      }
      pdev <- dev + penalty
      if (control$trace) 
        message(gettextf("penalized deviance = %s", pdev, 
                         domain = "R-mgcv"))
      div.thresh <- 10 * (0.1 + abs(old.pdev)) * .Machine$double.eps^0.5
      if (pdev - old.pdev > div.thresh) { #Might need to use this if statement later as a failsafe procedure
        ii <- 1
        if (iter == 1) {
          etaold <- null.eta
          coefold <- null.coef
        }
        while (pdev - old.pdev > div.thresh) { 
          if (ii > 100) 
            stop("inner loop 3; can't correct step size")
          ii <- ii + 1
          start <- (start + coefold)/2
          eta <- (eta + etaold)/2
          mu <- linkinv(eta)
          dev <- sum(dev.resids(y, mu, weights))
          pdev <- dev + t(start) %*% St %*% start
          if (control$trace) 
            message(gettextf("Step halved: new penalized deviance = %g", 
                             pdev, "\n"))
        }
      }
      if (strictly.additive) { 
        conv <- TRUE
        coef <- start
        break
      }
      if (abs(pdev - old.pdev)/(0.1 + abs(pdev)) < control$epsilon) { 
        grad <- 2 * t(x[good, ]) %*% (w * ((x %*% start)[good] - 
                                             z)) + 2 * St %*% start
        if (max(abs(grad)) > control$epsilon * max(abs(start + 
                                                       coefold))/2) {
          old.pdev <- pdev
          coef <- coefold <- start
          etaold <- eta
        }
        else {
          conv <- TRUE
          coef <- start
          etaold <- eta
          break
        }
      }
      else { 
        old.pdev <- pdev
        coef <- coefold <- start
        etaold <- eta
      }
    }
    wdr <- dev.resids(y, mu, weights) #calculate deviance residuals
    dev <- sum(wdr) #thus calculate deviance
    wdr <- sign(y - mu) * sqrt(pmax(wdr, 0))
    good <- weights > 0
    var.val <- variance(mu) 
    varmu <- var.val[good]
    if (any(is.na(varmu))) 
      stop("NAs in V(mu)")
    if (any(varmu == 0)) 
      stop("0s in V(mu)")
    mu.eta.val <- mu.eta(eta)
    if (any(is.na(mu.eta.val[good]))) 
      stop("NAs in d(mu)/d(eta)")
    good <- (weights > 0) & (mu.eta.val != 0)
    mevg <- mu.eta.val[good]
    mug <- mu[good]
    yg <- y[good]
    weg <- weights[good]
    etag <- eta[good]
    var.mug <- var.val[good]
    if (fisher) { 
      z <- (eta - offset)[good] + (yg - mug)/mevg
      w <- (weg * mevg^2)/var.mug
      alpha <- wf <- 0
    } else {
      c <- yg - mug
      alpha <- 1 + c * (family$dvar(mug)/var.mug + family$d2link(mug) * 
                          mevg)
      alpha[alpha == 0] <- .Machine$double.eps
      z <- (eta - offset)[good] + (yg - mug)/(mevg * alpha)
      wf <- weg * mevg^2/var.mug
      w <- wf * alpha
    }
    g1 <- 1/mevg
    g2 <- family$d2link(mug)
    g3 <- family$d3link(mug)
    V <- family$variance(mug)
    V1 <- family$dvar(mug)
    V2 <- family$d2var(mug)
    if (fisher) {
      g4 <- V3 <- 0
    } else {
      g4 <- family$d4link(mug)
      V3 <- family$d3var(mug)
    }
    if (TRUE) {
      g2 <- g2/g1
      g3 <- g3/g1
      g4 <- g4/g1
      V1 <- V1/V
      V2 <- V2/V
      V3 <- V3/V
    }
    P1 <- D1 <- array(0, nSp)
    P2 <- D2 <- matrix(0, nSp, nSp)
    trA1 <- array(0, nSp)
    trA2 <- matrix(0, nSp, nSp)
    rV = matrix(0, ncol(x), ncol(x))
    dum <- 1
    if (control$trace) 
      cat("calling gdi...")
    REML <- 0
    if (scoreType %in% c("REML", "P-REML")) {
      REML <- 1
      remlInd <- 1
    } else if (scoreType %in% c("ML", "P-ML")) {
      REML <- -1
      remlInd <- 0
    }
    if (REML == 0){
      rSncol <- unlist(lapply(rS, ncol))
    } else rSncol <- unlist(lapply(UrS, ncol))
    if (control$trace) 
      t1 <- proc.time()
    
    oo <- .C(mgcv:::C_gdi1, X = as.double(x[good, ]), E = as.double(Sr), 
             Eb = as.double(Eb), rS = as.double(unlist(rS)), U1 = as.double(U1), 
             sp = as.double(exp(sp)), z = as.double(z), w = as.double(w), 
             wf = as.double(wf), alpha = as.double(alpha), mu = as.double(mug), 
             eta = as.double(etag), y = as.double(yg), p.weights = as.double(weg), 
             g1 = as.double(g1), g2 = as.double(g2), g3 = as.double(g3), 
             g4 = as.double(g4), V0 = as.double(V), V1 = as.double(V1), 
             V2 = as.double(V2), V3 = as.double(V3), beta = as.double(coef), 
             b1 = as.double(rep(0, nSp * ncol(x))), w1 = as.double(rep(0,nSp * length(z))), D1 = as.double(D1), D2 = as.double(D2), 
             P = as.double(dum), P1 = as.double(P1), P2 = as.double(P2), 
             trA = as.double(dum), trA1 = as.double(trA1), trA2 = as.double(trA2), 
             rV = as.double(rV), rank.tol = as.double(rank.tol), 
             conv.tol = as.double(control$epsilon), rank.est = as.integer(1), 
             n = as.integer(length(z)), p = as.integer(ncol(x)), 
             M = as.integer(nSp), Mp = as.integer(Mp), Enrow = as.integer(rows.E), 
             rSncol = as.integer(rSncol), deriv = as.integer(deriv.sp), 
             REML = as.integer(REML), fisher = as.integer(fisher), 
             fixed.penalty = as.integer(rp$fixed.penalty), nthreads = as.integer(control$nthreads), 
             dVkk = as.double(rep(0, nSp * nSp)))
    
    d_oo <- d_gdi1(deriv = deriv, family = family$family, scoreType = scoreType, beta_star = beta_star, 
                   good = good, y = yg, X = x[good, ], d_X = d_x, d2_X = d2_x, 
                   mu = mug, variance = w, coef = coef, drho_coef = drho_coef, 
                   d2rho_coef = d2rho_coef, d_coef = d_coef, drho_d_coef = drho_d_coef, 
                   d2rho_d_coef = d2rho_d_coef, d2_coef = d2_coef, drho_d2_coef = drho_d2_coef, 
                   d2rho_d2_coef = d2rho_d2_coef, dev = dev, sp = as.numeric(exp(sp)), d_sp = d_sp, d2_sp = d2_sp, 
                   St=St, d_St = d_S, d2_St = d2_S, drho_S = drho_S, d2rho_S = d2rho_S, drho_d_S = drho_d_S, 
                   drho_d2_S = drho_d2_S, d2rho_d_S = d2rho_d_S, d2rho_d2_S = d2rho_d2_S)
    
    if (control$trace) {
      tg <- sum((proc.time() - t1)[c(1, 4)])
      cat("done!\n")
    }
    db.drho <- if (deriv){ 
      T %*% matrix(oo$b1, ncol(x), nSp) #We conduct this step
    } else NULL
    dw.drho <- if (deriv){ 
      matrix(oo$w1, length(z), nSp) #We conduct this step
    } else NULL
    rV <- matrix(oo$rV, ncol(x), ncol(x))
    Kmat <- matrix(0, nrow(x), ncol(x))
    Kmat[good, ] <- oo$X
    coef <- oo$beta
    eta <- drop(x %*% coef + offset)
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta))) {
      coef <- start
      eta <- etaold
      mu <- linkinv(eta)
    }
    trA <- oo$trA
    
    if(control$scale.est %in% c("pearson", "fletcher", "Pearson", 
                                "Fletcher")) { #We get in here when we're looking for scale.est
      pearson <- sum(weights * (y - mu)^2/family$variance(mu))
      scale.est <- (pearson + dev.extra)/(n.true - trA)
      if (control$scale.est %in% c("fletcher", "Fletcher")) { #We then get in here
        s.bar = mean(family$dvar(mu) * (y - mu) * sqrt(weights)/family$variance(mu))
        if (is.finite(s.bar)){ #We make this adjustment to scale.est and this is our final scale.est
          scale.est <- scale.est/(1 + s.bar)
        }
      }
    } else {
      scale.est <- (dev + dev.extra)/(n.true - trA)
    }
    reml.scale <- NA
    if (scoreType %in% c("REML", "ML")) {
      ls <- family$ls(y, weights, n, scale) * n.true/nobs
      Dp <- dev + oo$conv.tol + dev.extra
      REML <- Dp/(2 * scale) - ls[1] + oo$rank.tol/2 - 
        rp$det/2 - remlInd * Mp/2 * log(2 * pi * scale)
      
      d_REML <- vector(length=length(beta_star))
      d2_REML <- matrix(nrow=length(beta_star),ncol=length(beta_star))
      
      S_inv <- solve(rp$S)
      for(d1 in 1:length(beta_star))
      {
        d_ls <- (-nobs/2)*(d_scale[d1]/scale)
        d_Dp <- d_oo$d_dev[d1] + d_oo$d_conv.tol[d1]
        
        d_REML[d1] <- d_Dp/(2 * scale) - Dp*2*d_scale[d1]/(4 * scale * scale)  - d_ls + d_oo$d_rank.tol[d1]/2 -
          sum(diag( S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] ))/2 - remlInd * Mp *  d_scale[d1] /(2 *  scale)
        
        for(d2 in d1:length(beta_star))
        {
          d2_ls <- (-nobs/2)*( d2_scale[d1,d2]/scale - (d_scale[d1]*d_scale[d2])/(scale*scale) )
          d2_Dp <- d_oo$d2_dev[d1,d2] + d_oo$d2_conv.tol[d1,d2]
          
          d2_REML[d1,d2] <- ( (d_oo$d2_dev[d1,d2] + d_oo$d2_conv.tol[d1,d2])/(2 * scale) -
                                (d_oo$d_dev[d1] + d_oo$d_conv.tol[d1])*2*d_scale[d2]/(4 * scale * scale) ) -
            (((d_oo$d_dev[d2] + d_oo$d_conv.tol[d2])*2*d_scale[d1] + Dp*2*d2_scale[d1,d2])/(4 * scale * scale) -
               Dp*2*d_scale[d1]*8*scale*d_scale[d2]/(16 * scale * scale * scale * scale)) - d2_ls +
            d_oo$d2_rank.tol[d1,d2]/2 -
            sum(diag(-S_inv %*% d_S[[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] +
                       S_inv %*% d2_S[[d1]][[d2]][1:nrow(rp$S),1:ncol(rp$S)]))/2 -
            remlInd * Mp *  ( d2_scale[d1,d2]/(2 *  scale) - d_scale[d1]*2*d_scale[d2]/(4 * scale * scale))
        }
      }
      attr(REML, "Dp") <- Dp/(2 * scale)
      if (deriv) {
        d_UBRE <- d_UBRE <- d2_UBRE <- d_UBRE1 <- d2_UBRE1 <- d_UBRE2 <- d2_UBRE2 <- NULL
        
        REML1 <- oo$D1/(2 * scale) + oo$trA1/2 - rp$det1/2 
        
        d_REML1 <- d2_REML1 <- list()
        for(d1 in 1:length(beta_star))
        {
          d_REML1[[d1]] <- vector(length=length(sp))
          for(i in 1:length(sp))
          {
            d_REML1[[d1]][i] <- ( (d_oo$drho_d_dev[[d1]][i] + d_oo$drho_d_conv.tol[[d1]][i] )/(2 * scale) ) - 
              oo$D1[i]*2*d_scale[d1]/(4 * scale * scale) +
              d_oo$drho_d_rank.tol[[d1]][i]/2 - 
              sum(diag( -S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                          S_inv %*% drho_d_S[[d1]][[i]][1:nrow(rp$S),1:ncol(rp$S)]
              ))/2
          }
          
          d2_REML1[[d1]] <- list()
          for(d2 in d1:length(beta_star))
          {
            d2_REML1[[d1]][[d2]] <- vector(length=length(sp))
            for(i in 1:length(sp))
            {
              d2_REML1[[d1]][[d2]][i] <- 
                ( (d_oo$drho_d2_dev[[d1]][[d2]][i] + d_oo$drho_d2_conv.tol[[d1]][[d2]][i] )/(2 * scale) - 
                    ((d_oo$drho_d_dev[[d1]][i] + d_oo$drho_d_conv.tol[[d1]][i])*2*d_scale[d2] )/(4 * scale * scale) ) -
                ( ((d_oo$drho_d_dev[[d2]][i] + d_oo$drho_d_conv.tol[[d2]][i])*2*d_scale[d1] + oo$D1[i]*2*d2_scale[d1,d2])/(4 * scale * scale) - 
                    (oo$D1[i]*2*d_scale[d1]*8*scale*d_scale[d2])/(16 * scale * scale * scale * scale) ) +
                d_oo$drho_d2_rank.tol[[d1]][[d2]][i]/2 - 
                sum(diag(
                  S_inv %*% d_S[[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                    -S_inv %*% d2_S[[d1]][[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                    -S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% -S_inv %*% d_S[[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                    -S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_d_S[[d2]][[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                    
                    -S_inv %*% d_S[[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_d_S[[d1]][[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                    S_inv %*% drho_d2_S[[d1]][[d2]][[i]][1:nrow(rp$S),1:ncol(rp$S)]
                ))/2 
            }
          }
        }
        
        if (deriv == 2){
          REML2 <- (matrix(oo$D2, nSp, nSp)/scale + matrix(oo$trA2, nSp, nSp) - rp$det2)/2
          d_REML2 <- list()
          d2_REML2 <- list()
          for(d1 in 1:length(beta_star))
          {
            d_REML2[[d1]] <- matrix(nrow=length(sp),ncol=length(sp))
            for(i in 1:length(sp))
            {
              for(j in i:length(sp))
              {
                d_REML2[[d1]][i,j] <- ( (d_oo$d2rho_d_dev[[d1]][i,j] + d_oo$d2rho_d_conv.tol[[d1]][i,j] )/(2 * scale) -
                                          matrix(oo$D2, nSp, nSp)[i,j]*2*d_scale[d1]/(4 * scale * scale) ) + 
                  d_oo$d2rho_d_rank.tol[[d1]][i,j]/2 - 
                  sum(diag(
                    S_inv %*% drho_S[[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                      -S_inv %*% drho_d_S[[d1]][[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                      -S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% -S_inv %*% drho_S[[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                      -S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d2rho_S[[i]][[j]][1:nrow(rp$S),1:ncol(rp$S)] +  
                      -S_inv %*% drho_S[[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_d_S[[d1]][[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                      S_inv %*% d2rho_d_S[[d1]][[i]][[j]][1:nrow(rp$S),1:ncol(rp$S)]))/2
              }
            }
            
            d2_REML2[[d1]] <- list()
            for(d2 in d1:length(beta_star))
            {
              d2_REML2[[d1]][[d2]] <- matrix(nrow=length(sp),ncol=length(sp))              
              for(i in 1:length(sp))
              {
                for(j in i:length(sp))
                {
                  d2_REML2[[d1]][[d2]][i,j] <- ( (d_oo$d2rho_d2_dev[[d1]][[d2]][i,j] + d_oo$d2rho_d2_conv.tol[[d1]][[d2]][i,j] )/(2 * scale) -
                                                   (d_oo$d2rho_d_dev[[d1]][i,j] + d_oo$d2rho_d_conv.tol[[d1]][i,j] )*2*d_scale[d2]/(4 * scale * scale) ) - 
                    ( (d_oo$d2rho_d_dev[[d2]][i,j]*2*d_scale[d1] + matrix(oo$D2, nSp, nSp)[i,j]*2*d2_scale[d1,d2])/(4 * scale * scale) - 
                        matrix(oo$D2, nSp, nSp)[i,j]*2*d_scale[d1]*8*scale*d_scale[d2]/(16 * scale * scale * scale * scale) ) +
                    d_oo$d2rho_d2_rank.tol[[d1]][[d2]][i,j]/2 - 
                    sum(diag(
                      -S_inv %*% d_S[[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        S_inv %*% drho_d_S[[d2]][[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        S_inv %*% drho_S[[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% -S_inv %*% d_S[[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        S_inv %*% drho_S[[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d2_S[[d1]][[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        S_inv %*% drho_S[[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% -S_inv %*% d_S[[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        S_inv %*% drho_S[[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_d_S[[d2]][[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        S_inv %*% d_S[[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_d_S[[d1]][[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        -S_inv %*% drho_d2_S[[d1]][[d2]][[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        -S_inv %*% drho_d_S[[d1]][[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% -S_inv %*% d_S[[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        -S_inv %*% drho_d_S[[d1]][[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_d_S[[d2]][[i]][1:nrow(rp$S),1:ncol(rp$S)] + 
                        S_inv %*% d_S[[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% -S_inv %*% drho_S[[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        -S_inv %*% d2_S[[d1]][[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% -S_inv %*% drho_S[[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        -S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d_S[[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        -S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% -S_inv %*% drho_d_S[[d2]][[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        -S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% -S_inv %*% drho_S[[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% -S_inv %*% d_S[[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        -S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% -S_inv %*% drho_S[[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_d_S[[d2]][[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        S_inv %*% d_S[[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d2rho_S[[i]][[j]][1:nrow(rp$S),1:ncol(rp$S)] +
                        -S_inv %*% d2_S[[d1]][[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d2rho_S[[i]][[j]][1:nrow(rp$S),1:ncol(rp$S)] +
                        -S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% -S_inv %*% d_S[[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d2rho_S[[i]][[j]][1:nrow(rp$S),1:ncol(rp$S)] +
                        -S_inv %*% d_S[[d1]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d2rho_d_S[[d2]][[i]][[j]][1:nrow(rp$S),1:ncol(rp$S)] +
                        S_inv %*% d_S[[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_S[[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_d_S[[d1]][[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        -S_inv %*% drho_d_S[[d2]][[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_d_S[[d1]][[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        -S_inv %*% drho_S[[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% -S_inv %*% d_S[[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_d_S[[d1]][[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        -S_inv %*% drho_S[[j]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% drho_d2_S[[d1]][[d2]][[i]][1:nrow(rp$S),1:ncol(rp$S)] +
                        -S_inv %*% d_S[[d2]][1:nrow(rp$S),1:ncol(rp$S)] %*% S_inv %*% d2rho_d_S[[d1]][[i]][[j]][1:nrow(rp$S),1:ncol(rp$S)] +
                        S_inv %*% d2rho_d2_S[[d1]][[d2]][[i]][[j]][1:nrow(rp$S),1:ncol(rp$S)]
                    ))/2
                }
              }
            }
          }
        }
        if (sum(!is.finite(REML2))) {
          stop("Non finite derivatives. Try decreasing fit tolerance! See `epsilon' in `gam.contol'")
        }
      }
      if (!scale.known && deriv) { 
        dlr.dlphi <- -Dp/(2 * scale) - ls[2] * scale - Mp/2 * remlInd
        d2lr.d2lphi <- Dp/(2 * scale) - ls[3] * scale^2 - ls[2] * scale
        d2lr.dspphi <- -oo$D1/(2 * scale)
        REML1 <- c(REML1, dlr.dlphi)
        
        for(d1 in 1:length(beta_star))
        {
          d_dlr.dlphi <- ( -(d_oo$d_dev[d1] + d_oo$d_conv.tol[d1])/(2 * scale) - 
                             -Dp*2*d_scale[d1]/(4 * scale * scale) ) - 
            ( (nobs*d_scale[d1]/(2* scale * scale) * scale + 
                 ls[2] * d_scale[d1]) )
          d_REML1[[d1]] <- c(d_REML1[[d1]],d_dlr.dlphi)
          
          for(d2 in d1:length(beta_star))
          {
            d2_dlr.dlphi <- (
              -( (d_oo$d2_dev[d1,d2] + d_oo$d2_conv.tol[d1,d2])/(2 * scale) - 
                   (d_oo$d_dev[d1] + d_oo$d_conv.tol[d1])*2*d_scale[d2]/(4 * scale * scale) ) - 
                ((-(d_oo$d_dev[d2] + d_oo$d_conv.tol[d2])*2*d_scale[d1] + -Dp*2*d2_scale[d1,d2])/(4 * scale * scale) - 
                   -Dp*2*d_scale[d1]*8*scale*d_scale[d2]/(16 * scale * scale * scale * scale) ) ) - 
              ( (nobs*d2_scale[d1,d2]/(2 * scale) - nobs*d_scale[d1]*2*d_scale[d2]/(4 * scale * scale) )  + 
                  nobs*d_scale[d2]/(2 * scale * scale) * d_scale[d1] + ls[2] * d2_scale[d1,d2] )
            
            d2_REML1[[d1]][[d2]] <- c(d2_REML1[[d1]][[d2]],d2_dlr.dlphi)
          }
        }
        if (deriv == 2) {
          REML2 <- rbind(REML2, as.numeric(d2lr.dspphi))
          REML2 <- cbind(REML2, c(as.numeric(d2lr.dspphi), d2lr.d2lphi))
          
          for(d1 in 1:length(beta_star))
          {
            d_d2lr.dspphi <- -( d_oo$drho_d_dev[[d1]] + d_oo$drho_d_conv.tol[[d1]] )/(2 * scale) - -oo$D1*2*d_scale[d1]/(4 * scale * scale)
            d_d2lr.d2lphi <- ( (d_oo$d_dev[d1] + d_oo$d_conv.tol[d1])/(2 * scale) - Dp*2*d_scale[d1]/(4 * scale * scale) ) 
            d_REML2[[d1]] <- rbind(d_REML2[[d1]], as.numeric(d_d2lr.dspphi))
            d_REML2[[d1]] <- cbind(d_REML2[[d1]], c(as.numeric(d_d2lr.dspphi), d_d2lr.d2lphi))
            
            for(d2 in d1:length(beta_star))
            {
              d2_d2lr.dspphi <- ( -( d_oo$drho_d2_dev[[d1]][[d2]] + d_oo$drho_d2_conv.tol[[d1]][[d2]] )/(2 * scale) - 
                                    -( d_oo$drho_d_dev[[d1]] + d_oo$drho_d_conv.tol[[d1]] )*2*d_scale[d2]/(4 * scale * scale) ) - 
                ( (-d_oo$drho_d_dev[[d2]]*2*d_scale[d1] + -oo$D1*2*d2_scale[d1,d2])/(4 * scale * scale) -
                    -oo$D1*2*d_scale[d1]*8*scale*d_scale[d2]/(16 * scale * scale * scale * scale) )
              
              d2_d2lr.d2lphi <- ( (d_oo$d2_dev[d1,d2] + d_oo$d2_conv.tol[d1,d2])/(2 * scale) - 
                                    (d_oo$d_dev[d1] + d_oo$d_conv.tol[d1])*2*d_scale[d2]/(4 * scale * scale) ) - 
                ( ((d_oo$d_dev[d2] + d_oo$d_conv.tol[d2])*2*d_scale[d1] + Dp*2*d2_scale[d1,d2])/(4 * scale * scale) - 
                    Dp*2*d_scale[d1]*8*scale*d_scale[d2]/(16 * scale * scale * scale * scale) ) 
              
              d2_REML2[[d1]][[d2]] <- rbind(d2_REML2[[d1]][[d2]], as.numeric(d2_d2lr.dspphi))
              d2_REML2[[d1]][[d2]] <- cbind(d2_REML2[[d1]][[d2]], c(as.numeric(d2_d2lr.dspphi), d2_d2lr.d2lphi))
              
            }
          }
        }
      }
      reml.scale <- scale
    } else if (scoreType %in% c("P-REML", "P-ML")) {
      reml.scale <- phi <- (oo$P * (nobs - Mp) + pearson.extra)/(n.true - 
                                                                   Mp)
      oo$P1 <- oo$P1 * (nobs - Mp)/(n.true - Mp)
      oo$P2 <- oo$P2 * (nobs - Mp)/(n.true - Mp)
      ls <- family$ls(y, weights, n, phi) * n.true/nobs
      Dp <- dev + oo$conv.tol + dev.extra
      K <- oo$rank.tol/2 - rp$det/2
      REML <- Dp/(2 * phi) - ls[1] + K - Mp/2 * log(2 * 
                                                      pi * phi) * remlInd
      attr(REML, "Dp") <- Dp/(2 * phi)
      if (deriv) {
        phi1 <- oo$P1
        Dp1 <- oo$D1
        K1 <- oo$trA1/2 - rp$det1/2
        REML1 <- Dp1/(2 * phi) - phi1 * (Dp/(2 * phi^2) + 
                                           Mp/(2 * phi) * remlInd + ls[2]) + K1
        if (deriv == 2) {
          phi2 <- matrix(oo$P2, nSp, nSp)
          Dp2 <- matrix(oo$D2, nSp, nSp)
          K2 <- matrix(oo$trA2, nSp, nSp)/2 - rp$det2/2
          REML2 <- Dp2/(2 * phi) - (outer(Dp1, phi1) + 
                                      outer(phi1, Dp1))/(2 * phi^2) + (Dp/phi^3 - 
                                                                         ls[3] + Mp/(2 * phi^2) * remlInd) * outer(phi1, 
                                                                                                                   phi1) - (Dp/(2 * phi^2) + ls[2] + Mp/(2 * 
                                                                                                                                                           phi) * remlInd) * phi2 + K2
        }
      }
    } else {
      P <- oo$P
      delta <- nobs - gamma * trA
      delta.2 <- delta * delta
      GCV <- nobs * dev/delta.2
      GACV <- dev/nobs + P * 2 * gamma * trA/(delta * nobs) 
      UBRE <- dev/nobs - 2 * delta * scale/nobs + scale
      
      d_UBRE <- vector(length=length(beta_star))
      d2_UBRE <- matrix(nrow=length(beta_star), ncol=length(beta_star))
      for(d1 in 1:length(beta_star))
      {
        d_delta <- -gamma * d_oo$d_hat[d1]
        d_UBRE[d1] <- d_oo$d_dev[d1]/nobs  - 2 * d_delta * scale/nobs 
        for(d2 in d1:length(beta_star))
        {
          d2_delta <- -gamma * d_oo$d2_hat[d1,d2]
          d2_UBRE[d1,d2] <- d_oo$d2_dev[d1,d2]/nobs  - 2 * d2_delta * scale/nobs 
        }
      }
      
      if (deriv) {
        trA1 <- oo$trA1
        D1 <- oo$D1
        P1 <- oo$P1
        if (sum(!is.finite(D1)) || sum(!is.finite(P1)) || 
            sum(!is.finite(trA1))) {
          stop("Non-finite derivatives. Try decreasing fit tolerance! See `epsilon' in `gam.contol'")
        }
        delta.3 <- delta * delta.2
        GCV1 <- nobs * D1/delta.2 + 2 * nobs * dev * 
          trA1 * gamma/delta.3
        GACV1 <- D1/nobs + 2 * P/delta.2 * trA1 + 2 * 
          gamma * trA * P1/(delta * nobs)
        UBRE1 <- D1/nobs + gamma * trA1 * 2 * scale/nobs
        
        d_UBRE1 <- list()
        d2_UBRE1 <- list()
        for(d1 in 1:length(beta_star))
        {
          d_UBRE1[[d1]] <- vector(length=length(sp))
          d2_UBRE1[[d1]] <- list()
          for(d2 in d1:length(beta_star))
          {
            d2_UBRE1[[d1]][[d2]] <- vector(length=length(sp))
            for(i in 1:length(sp))
            {
              d_UBRE1[[d1]][i] <- d_oo$drho_d_dev[[d1]][i]/nobs + gamma * d_oo$drho_d_hat[[d1]][i] * 2 * scale/nobs  
              d2_UBRE1[[d1]][[d2]][i] <- d_oo$drho_d2_dev[[d1]][[d2]][i]/nobs + gamma * d_oo$drho_d2_hat[[d1]][[d2]][i] * 2 * scale/nobs
            }
          }
        }
        
        if (deriv == 2) {
          trA2 <- matrix(oo$trA2, nSp, nSp)
          D2 <- matrix(oo$D2, nSp, nSp)
          P2 <- matrix(oo$P2, nSp, nSp)
          if (sum(!is.finite(D2)) || sum(!is.finite(P2)) || 
              sum(!is.finite(trA2))) {
            stop("Non-finite derivatives. Try decreasing fit tolerance! See `epsilon' in `gam.contol'")
          }
          GCV2 <- outer(trA1, D1)
          GCV2 <- (GCV2 + t(GCV2)) * gamma * 2 * nobs/delta.3 + 
            6 * nobs * dev * outer(trA1, trA1) * gamma * 
            gamma/(delta.2 * delta.2) + nobs * D2/delta.2 + 
            2 * nobs * dev * gamma * trA2/delta.3
          GACV2 <- D2/nobs + outer(trA1, trA1) * 4 * 
            P/(delta.3) + 2 * P * trA2/delta.2 + 2 * 
            outer(trA1, P1)/delta.2 + 2 * outer(P1, trA1) * 
            (1/(delta * nobs) + trA/(nobs * delta.2)) + 
            2 * trA * P2/(delta * nobs)
          GACV2 <- (GACV2 + t(GACV2)) * 0.5
          UBRE2 <- D2/nobs + 2 * gamma * trA2 * scale/nobs
          
          d_UBRE2 <- list()
          d2_UBRE2 <- list()
          for(d1 in 1:length(beta_star))
          {
            d_UBRE2[[d1]] <- matrix(nrow=length(sp),ncol=length(sp))
            d2_UBRE2[[d1]] <- list()
            for(d2 in d1:length(beta_star))
            {
              d2_UBRE2[[d1]][[d2]] <- matrix(nrow=length(sp),ncol=length(sp))
              for(i in 1:length(sp))
              {
                for(j in i:length(sp))
                {
                  d_UBRE2[[d1]][i,j] <- d_oo$d2rho_d_dev[[d1]][i,j]/nobs + gamma * d_oo$d2rho_d_hat[[d1]][i,j] * 2 * scale/nobs 
                  d2_UBRE2[[d1]][[d2]][i,j] <- d_oo$d2rho_d2_dev[[d1]][[d2]][i,j]/nobs + gamma * d_oo$d2rho_d2_hat[[d1]][[d2]][i,j] * 2 * scale/nobs
                }
              }
            }
          }
        }
      }
    }
    if (!conv && printWarn) 
      warning("Algorithm did not converge")
    if (printWarn && boundary) 
      warning("Algorithm stopped at boundary value")
    eps <- 10 * .Machine$double.eps
    if (printWarn && family$family[1] == "binomial") {
      if (any(mu > 1 - eps) || any(mu < eps)) 
        warning("fitted probabilities numerically 0 or 1 occurred")
    }
    if (printWarn && family$family[1] == "poisson") {
      if (any(mu < eps)) 
        warning("fitted rates numerically 0 occurred")
    }
    residuals <- rep.int(NA, nobs)
    residuals[good] <- z - (eta - offset)[good]
    
    trans <- d_coef_trans(beta_star = beta_star, T_mat = T, d_T=d_T, d2_T=d2_T,coef = coef,drho_coef = drho_coef, 
                          d2rho_coef = d2rho_coef, d_coef = d_coef, d2_coef = d2_coef, drho_d_coef = drho_d_coef,
                          drho_d2_coef = drho_d2_coef, d2rho_d_coef = d2rho_d_coef, d2rho_d2_coef = d2rho_d2_coef, sp=sp)
    
    
    drho_coef <- trans$drho_coef
    d2rho_coef <- trans$d2rho_coef
    d_coef <- trans$d_coef
    drho_d_coef <- trans$drho_d_coef
    d2rho_d_coef <- trans$d2rho_d_coef
    d2_coef <- trans$d2_coef
    drho_d2_coef <- trans$drho_d2_coef
    d2rho_d2_coef <- trans$d2rho_d2_coef
    
    coef <- as.numeric(T %*% coef)
    rV <- T %*% rV
    names(coef) <- xnames
  }
  names(residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames
  ww <- wt <- rep.int(0, nobs)
  if (fisher) {
    wt[good] <- w
    ww <- wt
  } else {
    wt[good] <- wf
    ww[good] <- w
  }
  names(wt) <- ynames
  names(weights) <- ynames
  names(y) <- ynames
  if (deriv && nrow(dw.drho) != nrow(x)) {
    w1 <- dw.drho
    dw.drho <- matrix(0, nrow(x), ncol(w1))
    dw.drho[good, ] <- w1
  }
  wtdmu <- if (intercept){
    sum(weights * y)/sum(weights)
  } else{ 
    linkinv(offset)
  }
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)
  aic.model <- aic(y, n, mu, weights, dev)
  if (control$trace) {
    t1 <- proc.time()
    at <- sum((t1 - t0)[c(1, 4)])
    cat("Proportion time in C: ", (tc + tg)/at, " ls:", tc/at, 
        " gdi:", tg/at, "\n")
  }
  coef_derivs <- list(coef = coef, drho_coef = drho_coef, d2rho_coef = d2rho_coef, d_coef = d_coef, drho_d_coef = drho_d_coef,
                      d2rho_d_coef = d2rho_d_coef, d2_coef = d2_coef, drho_d2_coef = drho_d2_coef, d2rho_d2_coef = d2rho_d2_coef)
  
  list(coefficients = coef, residuals = residuals, fitted.values = mu, 
       family = family, linear.predictors = eta, deviance = dev, deviance_gradient = d_oo$d_dev, deviance_hessian = d_oo$d2_dev, 
       null.deviance = nulldev, iter = iter, weights = wt, working.weights = ww, 
       prior.weights = weights, df.null = nulldf, y = y, converged = conv, 
       boundary = boundary, D1 = D1, D2 = D2, P = P, P1 = P1, 
       P2 = P2, trA = trA, trA1 = trA1, trA2 = trA2, GCV = GCV, 
       GCV1 = GCV1, GCV2 = GCV2, GACV = GACV, GACV1 = GACV1, 
       GACV2 = GACV2, UBRE = UBRE, UBRE1 = UBRE1, UBRE2 = UBRE2, 
       REML = REML, REML1 = REML1, REML2 = REML2, rV = rV, db.drho = db.drho, 
       dw.drho = dw.drho, scale.est = scale.est, reml.scale = reml.scale, 
       aic = aic.model, rank = oo$rank.est, K = Kmat,iter=iter,T_matrix = T,
       coef_derivs = coef_derivs,
       mod_pen = mod_pen, d_UBRE = d_UBRE, d2_UBRE = d2_UBRE, d_UBRE1 = d_UBRE1, d2_UBRE1 = d2_UBRE1,
       d_UBRE2 = d_UBRE2, d2_UBRE2 = d2_UBRE2,
       d_REML = d_REML, d2_REML = d2_REML, d_REML1 = d_REML1, d2_REML1 = d2_REML1,
       d_REML2 = d_REML2, d2_REML2 = d2_REML2)
}

#08a
d_mm <- function(intercept, X, d_X, d2_X, Z, d_Z, d2_Z, T_mat, d_T, d2_T, beta_star, UrS){
  for(deriv_1 in 1:length(beta_star))
  {  
    for(deriv_2 in deriv_1:length(beta_star))
    {
      d2_X[[deriv_1]][[deriv_2]] <- (
        (d2_X[[deriv_1]][[deriv_2]] %*% T_mat) + (d_X[[deriv_1]] %*% d_T[[deriv_2]]) +
          (d_X[[deriv_2]] %*% d_T[[deriv_1]]) + (X %*% d2_T[[deriv_1]][[deriv_2]])
      )
      
    }
    d_X[[deriv_1]] <- (d_X[[deriv_1]] %*% T_mat) + (X %*% d_T[[deriv_1]])
  }
  list(d_X=d_X, d2_X=d2_X)
}

#08b
d_Pen <- function(beta_star, splineRank, rS, d_rS, d2_rS, T_mat, d_T, d2_T, UrS){
  Penalty <- d_Penalty <-  d2_Penalty <- list()
  #PenaltyRank <- vector()
  #This might need addressing, UrS is already re-param under U1, but in the case where we have gam.reparam where rp$Qs is not the identity, we will need to do soemthing different because we 
  #transform U1 via rp$Qs so will need derivatives - This I think is the procedure from Wood 2011
  for(i in 1:length(rS))
  {
    d_Penalty[[i]] <- d2_Penalty[[i]] <- list()
    
    Penalty[[i]] <- rS[[i]] %*% t(rS[[i]])
    
    if(i==1)
    {  
      for(deriv_1 in 1:length(beta_star))
      {
        d2_Penalty[[i]][[deriv_1]] <- list()
        for(deriv_2 in deriv_1:length(beta_star))
        {            
          d2_Penalty[[i]][[deriv_1]][[deriv_2]] <- (
            (d2_rS[[deriv_1]][[deriv_2]] %*% t(rS[[i]])) + (d_rS[[deriv_1]] %*% t(d_rS[[deriv_2]])) +
              (d_rS[[deriv_2]] %*% t(d_rS[[deriv_1]])) + (rS[[i]] %*% t(d2_rS[[deriv_1]][[deriv_2]]))
          )
          
        }
        
        d_Penalty[[i]][[deriv_1]] <- (d_rS[[deriv_1]] %*% t(rS[[i]])) + (rS[[i]] %*% t(d_rS[[deriv_1]]))
        
      }
    }
    else if(i!=1)
    {
      deriv_rS_temp <- matrix(0,nrow=nrow(rS[[i]]),ncol=ncol(rS[[i]]))
      for(deriv_1 in 1:length(beta_star))
      {
        d2_Penalty[[i]][[deriv_1]] <- list()
        for(deriv_2 in deriv_1:length(beta_star))
        {            
          d2_Penalty[[i]][[deriv_1]][[deriv_2]] <- (
            (deriv_rS_temp %*% t(rS[[i]])) + (deriv_rS_temp %*% t(deriv_rS_temp)) +
              (deriv_rS_temp %*% t(deriv_rS_temp)) + (rS[[i]] %*% t(deriv_rS_temp))
          )
          
        }
        
        d_Penalty[[i]][[deriv_1]] <- (deriv_rS_temp %*% t(rS[[i]])) + (rS[[i]] %*% t(deriv_rS_temp))
        
      }
    }
  }
  
  list(Penalty=Penalty, d_Penalty=d_Penalty, d2_Penalty=d2_Penalty)
}

#08c no longer in use

#08d
d_pls_fit0 <- function(beta_star, good, X, d_X, d2_X, variance, pseudodata, sp, d_sp, d2_sp, St, Penalty, d_Penalty, d2_Penalty){ 
  d_term <- d_inv_term <- d_coefstart <- drho_d_coefstart <- d2rho_d_term <- d2rho_d_coefstart <- list()
  d2_inv_term <- d2_term <- d2_coefstart <- drho_d2_coefstart <- d2rho_d2_coefstart <- list()
  d2rho_coefstart <- drho_coefstart <- list()
  
  #m <- ncol(X)
  #colNames <- colnames(X) <- as.character(1:ncol(X))
  
  singularityTest <- tryCatch(
    {
      test <- solve(crossprod(X * variance,X) + St)
      out <- "non-singular"
    },
    warning = function(war) {
      out <- "singular"
    },
    error = function(err){
      out <- "singular"
    }
  )
  
  if(singularityTest == "singular"){
    m <- ncol(X)
    colNames <- colnames(X) <- as.character(1:ncol(X))
    X_Decomp <- qr(sqrt(variance)*X, complete = FALSE)
    
    fullRank <-function (x, tol = 1e-07, qrx = qr(x, tol = tol)) 
    {
      keep_cols <- drop_cols <- rank <- NULL
      d <- dim(x)
      n <- d[[1L]]
      p <- d[[2L]]
      if (n < p){
        return(t(fullRank(t(x), tol = tol))) #transpose & then perform the routine
      }
      
      rnk <- qrx$rank
      if (rnk == p){
        x_full <- x
      } else {
        keep_cols <- qrx$pivot[seq_len(rnk)]
        drop_cols <- rep(0,q)
        drop_cols <- as.numeric(1:q %in% keep_cols)
        names(drop_cols) <- 1:q
        drop_cols <- names(drop_cols)[drop_cols == 0]
        
        x <- x[, keep_cols, drop = FALSE]
      }
      return(list(x = x, drop_cols = drop_cols, keep_cols = keep_cols, rank = rnk))
    }
    
    tt <- fullRank(x = sqrt(variance) * X, tol)
    tt2 <- fullRank(x = term)
    tt3 <- solve(tt2$x)
    
    
    form <- paste0(paste0("X",colnames(X)), collapse="+")
    form <- paste("y ~",form)
    form <- as.formula(form)
    
    test_lm <- lm(form, data = data.frame(X,y))
    
    drop <- X_Decomp$pivot[-c(1:X_Decomp$rank)]
    
    if(length(drop) > 0){
      X <- X[,!colnames(X) %in% drop]
      #We need to drop unidentifiable columns from a number of objects
      d_X <- sapply(d_X, function(d_X){
        colnames(d_X) <- 1:ncol(d_X)
        d_X <- d_X[,!colnames(d_X) %in% drop]
      },simplify = FALSE)
      d2_X <- sapply(1:length(d2_X), function(i){
        out <- sapply(d2_X[[i]], function(d2_X){
          if(ncol(d2_X)>0)
          {
            colnames(d2_X) <- 1:ncol(d2_X)
            d2_X <- d2_X[,!colnames(d2_X) %in% drop]
          }
        },simplify = FALSE)
        return(out)
      },simplify = FALSE)
    }
  }
  
  term <- crossprod(X,variance*X)
  inv_term <- solve( term )
  temp <- crossprod(X,variance*pseudodata)
  
  for(i in 1:length(sp))
  {
    d2rho_coefstart[[i]] <- list()
    for(j in i:length(sp))
    {
      d2rho_coefstart[[i]][[j]] <- drop(matrix(0,nrow = length(colNames), ncol = 1))
      names(d2rho_coefstart[[i]][[j]]) <- colNames
      
      d2rho_coefstart[[i]][[j]][var_ind] <- drop(matrix(0,nrow = ncol(X), ncol = 1))
    }
    drho_coefstart[[i]] <- drop(matrix(0,nrow = length(colNames), ncol = 1))
    names(drho_coefstart[[i]]) <- colNames
    
    drho_coefstart[[i]][var_ind] <- drop(matrix(0,nrow = ncol(X), ncol = 1))
  }
  
  for(d1 in 1:length(beta_star))
  { 
    drho_d_coefstart[[d1]] <- d2rho_d_coefstart[[d1]] <- list()
    
    d_X[[d1]] <- d_X[[d1]][good, ]
    d_term[[d1]] <- crossprod(d_X[[d1]],variance*X) + crossprod(X,variance*d_X[[d1]])  
    d_inv_term[[d1]] <- -1*(inv_term %*% (d_term[[d1]] %*% inv_term))
    
    d_coefstart[[d1]] <- drop(matrix(0,nrow = length(colNames), ncol = 1))
    names(d_coefstart[[d1]]) <- colNames
    
    d_coefstart[[d1]][var_ind] <- drop((d_inv_term[[d1]] %*% temp) + (inv_term %*% crossprod(d_X[[d1]],variance*pseudodata)))
    
    for(i in 1:length(sp))
    {
      drho_d_coefstart[[d1]][[i]] <- drop(matrix(0,nrow = length(colNames), ncol = 1))
      names(drho_d_coefstart[[d1]][[i]]) <- colNames
      
      drho_d_coefstart[[d1]][[i]][var_ind] <- drop(matrix(0,nrow = ncol(X), ncol = 1))
    }
    
    for(i in 1:length(sp))
    {
      d2rho_d_coefstart[[d1]][[i]] <- list()
      for(j in i:length(sp))
      { 
        d2rho_d_coefstart[[d1]][[i]][[j]] <- drop(matrix(0,nrow = length(colNames), ncol = 1))
        names(d2rho_d_coefstart[[d1]][[i]][[j]]) <- colNames
        
        d2rho_d_coefstart[[d1]][[i]][[j]][var_ind] <- drop(matrix(0,nrow = ncol(X), ncol = 1))
      }
    }
  }  
  d2_inv_term <- d2_term <- list()
  for(d1 in 1:length(beta_star))
  {
    d2_inv_term[[d1]] <- d2_term[[d1]] <- d2_coefstart[[d1]] <- drho_d2_coefstart[[d1]] <- d2rho_d2_coefstart[[d1]] <- list()
    for(d2 in d1:length(beta_star))
    {    
      drho_d2_coefstart[[d1]][[d2]] <- d2rho_d2_coefstart[[d1]][[d2]] <- list()
      d2_X[[d1]][[d2]] <- d2_X[[d1]][[d2]][good, ]
      
      d2_term[[d1]][[d2]] <- crossprod(d2_X[[d1]][[d2]],variance*X) + 
        crossprod(d_X[[d1]],variance*d_X[[d2]]) + 
        crossprod(d_X[[d2]],variance*d_X[[d1]]) + 
        crossprod(X,variance*d2_X[[d1]][[d2]]) 
      
      d2_inv_term[[d1]][[d2]] <- -1*( d_inv_term[[d2]] %*% (d_term[[d1]] %*% inv_term) + 
                                        inv_term %*% (d2_term[[d1]][[d2]] %*% inv_term) + 
                                        inv_term %*% (d_term[[d1]] %*% d_inv_term[[d2]]) )
      
      d2_coefstart[[d1]][[d2]] <- drop(matrix(0,nrow = length(colNames), ncol = 1))
      names(d2_coefstart[[d1]][[d2]]) <- colNames
      
      d2_coefstart[[d1]][[d2]][var_ind] <- drop(d2_inv_term[[d1]][[d2]] %*% temp +
                                                  d_inv_term[[d1]] %*% (crossprod(d_X[[d2]],variance*pseudodata)) +
                                                  d_inv_term[[d2]] %*% (crossprod(d_X[[d1]],variance*pseudodata)) +
                                                  inv_term %*% (crossprod(d2_X[[d1]][[d2]],variance*pseudodata)))
      
      for(i in 1:length(sp))
      {
        drho_d2_coefstart[[d1]][[d2]][[i]] <- drop(matrix(0,nrow = length(colNames), ncol = 1))
        names(drho_d2_coefstart[[d1]][[d2]][[i]]) <- colNames
        
        drho_d2_coefstart[[d1]][[d2]][[i]][var_ind]  <- drop(matrix(0,nrow = ncol(X), ncol = 1))
      }
      
      for(i in 1:length(sp))
      {
        d2rho_d2_coefstart[[d1]][[d2]][[i]] <- list()
        for(j in i:length(sp))
        {
          d2rho_d2_coefstart[[d1]][[d2]][[i]][[j]] <- drop(matrix(0,nrow = length(colNames), ncol = 1))
          names(d2rho_d2_coefstart[[d1]][[d2]][[i]][[j]]) <- colNames
          
          d2rho_d2_coefstart[[d1]][[d2]][[i]][[j]][var_ind]  <- drop(matrix(0,nrow = ncol(X), ncol = 1))
        }
      }
    }    
  }
  
  list(drho_coefstart=drho_coefstart, d2rho_coefstart=d2rho_coefstart, d_coefstart=d_coefstart, d2_coefstart=d2_coefstart, drho_d_coefstart=drho_d_coefstart, 
       d2rho_d_coefstart=d2rho_d_coefstart, drho_d2_coefstart = drho_d2_coefstart, d2rho_d2_coefstart=d2rho_d2_coefstart)
}

#08e
d_pls_fit1 <- function(deriv = 2, family, beta_star, good, y, X, d_X, d2_X, mu, 
                       variance, pseudodata, sp, d_sp, d2_sp, St, d_St, d2_St,
                       drho_S, d2rho_S, drho_d_S, drho_d2_S, d2rho_d_S, d2rho_d2_S,
                       coef, drho_coef, d2rho_coef, d_coef, drho_d_coef, d2rho_d_coef, 
                       d2_coef, drho_d2_coef, d2rho_d2_coef){
  
  n <- nrow(X)
  t10 <- variance*pseudodata
  t1 <- crossprod(X,t10)
  
  drho_d_eta <- drho_d_mu <- drho_d_pseudodata <- drho_d_variance <- drho_d_term <- drho_d_inv_term <- d2rho_d_eta <- d2rho_d_mu <- d2rho_d_pseudodata <- d2rho_d_variance <- d2rho_d_term <- d2rho_d_inv_term <- list()
  d2_eta <- d2_mu <- d2_pseudodata <- d2_variance <- drho_d2_eta <- drho_d2_mu <- drho_d2_pseudodata <- drho_d2_variance <- drho_d2_inv_term <- d2rho_d_eta <- d2rho_d_mu <- d2rho_d2_pseudodata <- d2rho_d2_variance <- d2rho_d2_term <- d2rho_d2_inv_term <- drho_d2_term <- list()
  d2_term <- d2_inv_term <- d2rho_d2_eta <- d2rho_d2_mu <- d2rho_d2_variance <- d2rho_d2_pseudodata <- d2rho_eta <- d2rho_mu <- d2rho_variance <- d2rho_pseudodata <- d2rho_term <- d2rho_inv_term <- list()
  
  singularityTest <- tryCatch(
    {
      test <- solve(crossprod(X * variance,X) + St)
      out <- "non-singular"
    },
    warning = function(war) {
      out <- "singular"
    },
    error = function(err){
      out <- "singular"
    }
  )
  
  if(singularityTest == "singular"){
    m <- ncol(X)
    colNames <- colnames(X) <- as.character(1:ncol(X))
    X_Decomp <- qr(X, complete = FALSE)
    drop <- X_Decomp$pivot[-c(1:X_Decomp$rank)]
    
    
    if(length(drop) > 0){
      X <- X[,!colnames(X) %in% drop]
      #We need to drop unidentifiable columns from a number of objects
      d_X <- sapply(d_X, function(d_X){
        colnames(d_X) <- 1:ncol(d_X)
        d_X <- d_X[,!colnames(d_X) %in% drop]
      },simplify = FALSE)
      d2_X <- sapply(1:length(d2_X), function(i){
        out <- sapply(d2_X[[i]], function(d2_X){
          if(ncol(d2_X)>0)
          {
            colnames(d2_X) <- 1:ncol(d2_X)
            d2_X <- d2_X[,!colnames(d2_X) %in% drop]
          }
        },simplify = FALSE)
        return(out)
      },simplify = FALSE)
    }
  }
  
  inv_variance <- 1/variance
  inv_variance[inv_variance==Inf] <- 10^100
  inv_term <- solve( crossprod(X * variance,X)  + St )
  
  if(deriv){ #First derivative
    drho_eta <- lapply(drho_coef,FUN= function(input) X %*% input)
    drho_mu <- lapply(drho_eta, "*", variance)
    drho_variance <- drho_var_fn(family=family, mu, drho_mu)
    drho_pseudodata <- drho_pseudo_fn(family = family, y = y, inv_variance = inv_variance, mu = mu, drho_mu = drho_mu, drho_eta = drho_eta)
    temp <- lapply(drho_variance, FUN=function(input) crossprod(X*input ,X)  )
    drho_term <- mapply("+" ,temp ,drho_S,SIMPLIFY =F)
    drho_inv_term <- lapply(drho_term, FUN=function(input) -1*inv_term %*% input %*% inv_term )
    
    t2 <- lapply(drho_variance, FUN= function(input) crossprod(X,input*pseudodata) )
    t3 <- lapply(drho_pseudodata, FUN= function(input) crossprod(X,variance*input) )
    t2t3 <- mapply(FUN=function(input1,input2) input1+input2, t2,t3, SIMPLIFY =F)
    t11 <- mapply(FUN=function(input1,input2) input1+input2, lapply(drho_variance, FUN=function(input) input*pseudodata), lapply(drho_pseudodata, FUN=function(input) variance*input) ,SIMPLIFY =F)  
    
  }
  
  if(deriv == 2){
    t16 <- t41 <- regularExp <- list()
    for(i in 1:length(sp))
    {
      d2rho_eta[[i]] <- d2rho_mu[[i]] <- d2rho_variance[[i]] <- d2rho_pseudodata[[i]] <- d2rho_term[[i]] <- d2rho_inv_term[[i]] <- t16[[i]] <- t41[[i]] <- regularExp[[i]] <- list()
      for(j in i:length(sp))
      {
        d2rho_eta[[i]][[j]] <- tcrossprod(X,t(d2rho_coef[[i]][[j]]))
        d2rho_mu[[i]][[j]] <- d2rho_eta[[i]][[j]]*variance + drho_eta[[i]]*drho_variance[[j]] 
        d2rho_variance[[i]][[j]] <- as.vector(d2rho_var_fn(family=family, mu=mu, drho_mu_i=drho_mu[[i]], drho_mu_j=drho_mu[[j]], d2rho_mu=d2rho_mu[[i]][[j]]))
        d2rho_pseudodata[[i]][[j]] <- d2rho_pseudo_fn(family = family, y = y, inv_variance = inv_variance, mu = mu, drho_mu_i = drho_mu[[i]], drho_mu_j = drho_mu[[j]], d2rho_mu = d2rho_mu[[i]][[j]], 
                                                      d2rho_eta = d2rho_eta[[i]][[j]])
        
        t16[[i]][[j]] <- d2rho_variance[[i]][[j]]*X
        
        d2rho_term[[i]][[j]] <- crossprod(X,t16[[i]][[j]]) + d2rho_S[[i]][[j]]
        
        d2rho_inv_term[[i]][[j]] <- -1*( drho_inv_term[[j]] %*% drho_term[[i]] %*% inv_term + inv_term %*% (d2rho_term[[i]][[j]] %*% inv_term + drho_term[[i]] %*% drho_inv_term[[j]]) )
        
        regularExp[[i]][[j]] <- d2rho_variance[[i]][[j]]*pseudodata + drho_variance[[i]]*drho_pseudodata[[j]] + drho_variance[[j]]*drho_pseudodata[[i]] + variance*d2rho_pseudodata[[i]][[j]]
        t41[[i]][[j]] <- crossprod(X,regularExp[[i]][[j]])
      }
    }
  }
  
  if(sum(good) != n){
    d_X <- lapply(d_X, FUN=function(input) input[good, ])
  }
  
  d_eta <- mapply("+",lapply(d_X,FUN= function(input) tcrossprod(input,t(coef))), lapply(d_coef, FUN=function(input) tcrossprod(X,t(input))) ,SIMPLIFY =F)
  d_mu <- lapply(d_eta, "*", variance)
  d_variance <- d_var_fn(family=family, mu, d_mu)
  d_pseudodata <- d_pseudo_fn(family = family, y=y, inv_variance = inv_variance, mu = mu, d_mu = d_mu, d_eta = d_eta)
  
  temp1 <- lapply(d_X, FUN=function(input) crossprod(input*variance ,X))
  temp2 <- lapply(d_variance, FUN=function(input) crossprod(X*input ,X))
  temp3 <- lapply(d_X, FUN=function(input) crossprod(X,variance*input))
  
  d_term <- mapply(FUN=function(input1,input2,input3,input4) input1+input2+input3+input4 ,temp1 ,temp2 ,temp3 ,d_St ,SIMPLIFY =F)
  d_inv_term <- lapply(d_term, FUN=function(input) -1*inv_term %*% input %*% inv_term )
  
  t4 <- lapply(d_X, FUN=function(input) crossprod(input,t10) )
  t5 <- lapply(d_variance, FUN=function(input) crossprod(X,input*pseudodata) )
  t6 <- lapply(d_pseudodata, FUN=function(input) crossprod(X,variance*input) )
  t8 <- lapply(d_term, FUN=function(input) input %*% inv_term)
  t9 <- mapply(FUN=function(input1,input2) input1+input2, lapply(d_variance, FUN=function(input) input*pseudodata), lapply(d_pseudodata, FUN=function(input) variance*input) ,SIMPLIFY =F)  
  
  crossT2 <- t17 <- t19 <- t20 <- t21 <- t26 <- t27 <- t35 <- t36 <- t37 <- t40 <- list()
  for(d1 in 1:(length(beta_star)))
  { 
    if(deriv){
      drho_d_eta[[d1]] <- drho_d_mu[[d1]] <- drho_d_variance[[d1]] <- drho_d_pseudodata[[d1]] <- drho_d_term[[d1]] <- drho_d_inv_term[[d1]] <- d2rho_d_eta[[d1]] <- d2rho_d_mu[[d1]] <- d2rho_d_variance[[d1]] <- d2rho_d_pseudodata[[d1]] <- d2rho_d_term[[d1]] <- d2rho_d_inv_term[[d1]] <- list()
      crossT2[[d1]] <- t17[[d1]] <- t19[[d1]] <- t20[[d1]] <- t21[[d1]] <- t26[[d1]] <- t27[[d1]] <- t35[[d1]] <- t36[[d1]] <- t37[[d1]] <- t40[[d1]] <- list()
      for(i in 1:length(sp))
      {
        drho_d_eta[[d1]][[i]] <- drop( tcrossprod(d_X[[d1]] , t(drho_coef[[i]]) )  + tcrossprod(X , t(drho_d_coef[[d1]][[i]])) )
        
        drho_d_mu[[d1]][[i]] <- drho_d_eta[[d1]][[i]] * variance + d_eta[[d1]] * drho_variance[[i]]
        
        drho_d_variance[[d1]][[i]] <- as.vector(drho_d_var_fn(family=family, mu=mu, d_mu=d_mu[[d1]], drho_mu=drho_mu[[i]], drho_d_mu=drho_d_mu[[d1]][[i]]))
        drho_d_pseudodata[[d1]][[i]] <- drho_d_pseudo_fn(family = family, y = y, inv_variance = inv_variance, mu = mu, d_mu = d_mu[[d1]], drho_mu = drho_mu[[i]], 
                                                         drho_d_mu = drho_d_mu[[d1]][[i]], drho_d_eta = drho_d_eta[[d1]][[i]])
        
        drho_d_term[[d1]][[i]] <- crossprod(d_X[[d1]],drho_variance[[i]]*X) + 
          crossprod(X,drho_d_variance[[d1]][[i]]*X + drho_variance[[i]]*d_X[[d1]]) + 
          drho_d_S[[d1]][[i]]       
        
        drho_d_inv_term[[d1]][[i]] <- -1*( drho_inv_term[[i]] %*% t8[[d1]] + 
                                             inv_term %*% ((drho_d_term[[d1]][[i]] %*% inv_term) + (d_term[[d1]] %*% drho_inv_term[[i]]))  
        ) 
        
        crossT2[[d1]][[i]] <- crossprod(d_X[[d1]],t11[[i]])
        t17[[d1]][[i]] <- drho_d_variance[[d1]][[i]]*pseudodata
        t19[[d1]][[i]] <- t17[[d1]][[i]] + d_variance[[d1]]*drho_pseudodata[[i]]
        t20[[d1]][[i]] <- drho_variance[[i]]*d_pseudodata[[d1]] + variance*drho_d_pseudodata[[d1]][[i]]
        t21[[d1]][[i]] <- crossprod(X,t19[[d1]][[i]] + t20[[d1]][[i]])
        t26[[d1]][[i]] <- t17[[d1]][[i]] + drho_variance[[i]]*d_pseudodata[[d1]]
        t27[[d1]][[i]] <- d_variance[[d1]]*drho_pseudodata[[i]] + variance*drho_d_pseudodata[[d1]][[i]]
        t37[[d1]][[i]] <- crossprod(X,t26[[d1]][[i]] + t27[[d1]][[i]])
      }
    }
    
    if(deriv == 2){
      for(i in 1:length(sp))
      {
        d2rho_d_eta[[d1]][[i]] <- d2rho_d_mu[[d1]][[i]] <- d2rho_d_variance[[d1]][[i]] <- d2rho_d_pseudodata[[d1]][[i]] <- d2rho_d_term[[d1]][[i]] <- d2rho_d_inv_term[[d1]][[i]] <- list()
        t35[[d1]][[i]] <- t36[[d1]][[i]] <- t40[[d1]][[i]] <- list()
        for(j in i:length(sp))
        {
          d2rho_d_eta[[d1]][[i]][[j]] <- drop( tcrossprod(d_X[[d1]] , t(d2rho_coef[[i]][[j]])) + tcrossprod(X , t(d2rho_d_coef[[d1]][[i]][[j]])) )
          
          d2rho_d_mu[[d1]][[i]][[j]] <- d2rho_d_eta[[d1]][[i]][[j]] * variance + drho_d_eta[[d1]][[i]] * drho_variance[[j]] + drho_d_eta[[d1]][[j]] * drho_variance[[i]] + d_eta[[d1]] * d2rho_variance[[i]][[j]]
          
          d2rho_d_variance[[d1]][[i]][[j]] <- as.vector(d2rho_d_var_fn(family=family, mu=mu, d_mu=d_mu[[d1]], drho_mu_i=drho_mu[[i]], drho_mu_j=drho_mu[[j]], 
                                                                       drho_d_mu_i=drho_d_mu[[d1]][[i]],
                                                                       drho_d_mu_j=drho_d_mu[[d1]][[j]],d2rho_mu=d2rho_mu[[i]][[j]], d2rho_d_mu=d2rho_d_mu[[d1]][[i]][[j]]))
          
          d2rho_d_pseudodata[[d1]][[i]][[j]] <- d2rho_d_pseudo_fn(family = family, y = y, inv_variance = inv_variance, variance = variance, d_variance = d_variance[[d1]], 
                                                                  drho_variance_i = drho_variance[[i]], drho_variance_j = drho_variance[[j]], 
                                                                  drho_d_variance_i = drho_d_variance[[d1]][[i]], drho_d_variance_j = drho_d_variance[[d1]][[j]], 
                                                                  d2rho_variance = d2rho_variance[[i]][[j]], d2rho_d_variance = d2rho_d_variance[[d1]][[i]][[j]],
                                                                  mu = mu, d_mu = d_mu[[d1]], drho_mu_i = drho_mu[[i]], drho_mu_j = drho_mu[[j]], 
                                                                  drho_d_mu_i = drho_d_mu[[d1]][[i]], drho_d_mu_j = drho_d_mu[[d1]][[j]], d2rho_mu = d2rho_mu[[i]][[j]], 
                                                                  d2rho_d_mu = d2rho_d_mu[[d1]][[i]][[j]], d2rho_d_eta = d2rho_d_eta[[d1]][[i]][[j]])
          
          d2rho_d_term[[d1]][[i]][[j]] <- crossprod(d_X[[d1]] ,t16[[i]][[j]]) + crossprod(X ,d2rho_d_variance[[d1]][[i]][[j]]*X + d2rho_variance[[i]][[j]]*d_X[[d1]]) + d2rho_d_S[[d1]][[i]][[j]]
          
          d2rho_d_inv_term[[d1]][[i]][[j]] <- -1*( 
            d2rho_inv_term[[i]][[j]] %*% t8[[d1]] + 
              drho_inv_term[[i]] %*% (drho_d_term[[d1]][[j]] %*% inv_term + d_term[[d1]] %*% drho_inv_term[[j]]) + 
              
              drho_inv_term[[j]] %*% (drho_d_term[[d1]][[i]] %*% inv_term + d_term[[d1]] %*% drho_inv_term[[i]]) + 
              
              inv_term %*% (d2rho_d_term[[d1]][[i]][[j]] %*% inv_term + drho_d_term[[d1]][[i]] %*% drho_inv_term[[j]] +
                              drho_d_term[[d1]][[j]] %*% drho_inv_term[[i]] + d_term[[d1]] %*% d2rho_inv_term[[i]][[j]]
              ) )
          
          t36[[d1]][[i]][[j]] <- d2rho_d_variance[[d1]][[i]][[j]]*pseudodata + drho_d_variance[[d1]][[i]]*drho_pseudodata[[j]] +
            d2rho_variance[[i]][[j]]*d_pseudodata[[d1]] + drho_variance[[i]]*drho_d_pseudodata[[d1]][[j]] +
            drho_d_variance[[d1]][[j]]*drho_pseudodata[[i]] + d_variance[[d1]]*d2rho_pseudodata[[i]][[j]] +
            drho_variance[[j]]*drho_d_pseudodata[[d1]][[i]] + variance*d2rho_d_pseudodata[[d1]][[i]][[j]]
          
          t35[[d1]][[i]][[j]] <- crossprod(X,t36[[d1]][[i]][[j]])
          
          t40[[d1]][[i]][[j]] <- crossprod(d_X[[d1]],regularExp[[i]][[j]])
        }
      } 
    }
  }
  
  crossT1 <- t22 <- t23 <- t24 <- t25 <- t30 <- t38 <- t39 <- list()
  for(d1 in 1:length(beta_star))
  {
    crossT1[[d1]] <- list()
    d2_eta[[d1]] <- d2_mu[[d1]] <- d2_pseudodata[[d1]] <- d2_variance[[d1]] <- drho_d2_eta[[d1]] <- drho_d2_mu[[d1]] <- drho_d2_pseudodata[[d1]] <- drho_d2_variance[[d1]] <- drho_d2_inv_term[[d1]] <- d2rho_d2_eta[[d1]] <- d2rho_d2_mu[[d1]] <- d2rho_d2_pseudodata[[d1]] <- d2rho_d2_variance[[d1]] <- d2rho_d2_term[[d1]] <- d2rho_d2_inv_term[[d1]] <- drho_d2_term[[d1]] <- list()
    d2_term[[d1]] <- d2_inv_term[[d1]] <- d2rho_d2_variance[[d1]] <- d2rho_d2_pseudodata[[d1]] <- t22[[d1]] <- t23[[d1]] <- t24[[d1]] <- t25[[d1]] <- t30[[d1]] <- t38[[d1]] <- t39[[d1]] <- list()
    for(d2 in d1:length(beta_star))
    {  
      if(sum(good) != n){
        d2_X[[d1]][[d2]] <- d2_X[[d1]][[d2]][good, ]     
      }
      
      d2_eta[[d1]][[d2]] <- tcrossprod(d2_X[[d1]][[d2]] , t(coef)) + tcrossprod(d_X[[d1]] , t(d_coef[[d2]])) + 
        tcrossprod(d_X[[d2]] , t(d_coef[[d1]])) + tcrossprod(X , t(d2_coef[[d1]][[d2]])) 
      
      d2_mu[[d1]][[d2]] <- d2_eta[[d1]][[d2]]*variance + d_eta[[d1]]*d_variance[[d2]]
      d2_variance[[d1]][[d2]] <- as.vector(d2_var_fn(family=family, mu, d_mu_1=d_mu[[d1]], d_mu_2=d_mu[[d2]], d2_mu=d2_mu[[d1]][[d2]]))
      
      d2_pseudodata[[d1]][[d2]] <- d2_pseudo_fn(family = family, y = y, inv_variance = inv_variance, variance = variance, d_variance_1 = d_variance[[d1]], 
                                                d_variance_2 = d_variance[[d2]], d2_variance = d2_variance[[d1]][[d2]], mu = mu, 
                                                d_mu_1 = d_mu[[d1]], d_mu_2 = d_mu[[d2]], d2_mu = d2_mu[[d1]][[d2]], d2_eta = d2_eta[[d1]][[d2]])
      
      d2_term[[d1]][[d2]] <-  crossprod(d2_X[[d1]][[d2]]*variance + d_X[[d1]]*d_variance[[d2]] + d_X[[d2]]*d_variance[[d1]] + X*d2_variance[[d1]][[d2]],X) + 
        crossprod(d_X[[d1]]*variance + X*d_variance[[d1]],d_X[[d2]]) +
        crossprod(d_X[[d2]]*variance + X*d_variance[[d2]],d_X[[d1]]) + 
        crossprod(X*variance,d2_X[[d1]][[d2]]) +
        d2_St[[d1]][[d2]]
      
      d2_inv_term[[d1]][[d2]] <- -1*(d_inv_term[[d2]]  %*% t8[[d1]] + inv_term %*% (d2_term[[d1]][[d2]] %*% inv_term + d_term[[d1]] %*% d_inv_term[[d2]]) ) 
      
      crossT1[[d1]][[d2]] <- crossprod(d2_X[[d1]][[d2]],t10)
      t24[[d1]][[d2]] <- crossprod(d_X[[d1]],t9[[d2]])
      t25[[d1]][[d2]] <- crossprod(d_X[[d2]],t9[[d1]])
      t30[[d1]][[d2]] <- crossprod(X,d2_variance[[d1]][[d2]]*pseudodata + d_variance[[d1]]*d_pseudodata[[d2]] + d_variance[[d2]]*d_pseudodata[[d1]] + variance*d2_pseudodata[[d1]][[d2]])
      
      d2_coef[[d1]][[d2]] <- as.vector( d2_inv_term[[d1]][[d2]] %*% t1 + 
                                          d_inv_term[[d1]] %*% (t4[[d2]] + t5[[d2]] + t6[[d2]]) +               
                                          d_inv_term[[d2]] %*% (t4[[d1]] + t5[[d1]] + t6[[d1]]) +
                                          inv_term %*% (crossT1[[d1]][[d2]] + t24[[d1]][[d2]] + t25[[d1]][[d2]] + t30[[d1]][[d2]] )
      ) 
      
      if(deriv){
        drho_d2_eta[[d1]][[d2]] <- drho_d2_mu[[d1]][[d2]] <- drho_d2_variance[[d1]][[d2]] <- drho_d2_pseudodata[[d1]][[d2]] <- drho_d2_term[[d1]][[d2]] <- drho_d2_inv_term[[d1]][[d2]] <- list()
        d2rho_d2_eta[[d1]][[d2]] <- d2rho_d2_mu[[d1]][[d2]] <- d2rho_d2_variance[[d1]][[d2]] <- d2rho_d2_pseudodata[[d1]][[d2]] <- d2rho_d2_term[[d1]][[d2]] <- d2rho_d2_inv_term[[d1]][[d2]] <- list()
        t22[[d1]][[d2]] <- t23[[d1]][[d2]] <- t38[[d1]][[d2]] <- t39[[d1]][[d2]] <- list()
        for(i in 1:length(sp))
        {
          drho_d2_eta[[d1]][[d2]][[i]] <- tcrossprod(d2_X[[d1]][[d2]] , t(drho_coef[[i]])) + tcrossprod(d_X[[d1]] , t(drho_d_coef[[d2]][[i]])) + 
            tcrossprod(d_X[[d2]] , t(drho_d_coef[[d1]][[i]])) + tcrossprod(X , t(drho_d2_coef[[d1]][[d2]][[i]]))
          
          drho_d2_mu[[d1]][[d2]][[i]] <- drho_d2_eta[[d1]][[d2]][[i]] * variance + drho_d_eta[[d1]][[i]] * d_variance[[d2]] + d2_eta[[d1]][[d2]] * drho_variance[[i]] + d_eta[[d1]] * drho_d_variance[[d2]][[i]]        
          drho_d2_variance[[d1]][[d2]][[i]] <- as.vector(drho_d2_var_fn(family=family, mu, d_mu_1=d_mu[[d1]], d_mu_2=d_mu[[d2]], drho_mu=drho_mu[[i]], 
                                                                        drho_d_mu_1=drho_d_mu[[d1]][[i]], drho_d_mu_2=drho_d_mu[[d2]][[i]], d2_mu=d2_mu[[d1]][[d2]], 
                                                                        drho_d2_mu=drho_d2_mu[[d1]][[d2]][[i]]))
          drho_d2_pseudodata[[d1]][[d2]][[i]] <- drho_d2_pseudo_fn(family = family, y = y, inv_variance = inv_variance, 
                                                                   variance = variance, d_variance_1 = d_variance[[d1]], 
                                                                   d_variance_2 = d_variance[[d2]], drho_variance = drho_variance[[i]], 
                                                                   drho_d_variance_1 = drho_d_variance[[d1]][[i]], 
                                                                   drho_d_variance_2 = drho_d_variance[[d2]][[i]], 
                                                                   d2_variance = d2_variance[[d1]][[d2]], 
                                                                   drho_d2_variance = drho_d2_variance[[d1]][[d2]][[i]], 
                                                                   mu = mu, d_mu_1 = d_mu[[d1]], d_mu_2 = d_mu[[d2]], 
                                                                   drho_mu = drho_mu[[i]], 
                                                                   drho_d_mu_1 = drho_d_mu[[d1]][[i]], 
                                                                   drho_d_mu_2 = drho_d_mu[[d2]][[i]], 
                                                                   d2_mu = d2_mu[[d1]][[d2]], 
                                                                   drho_d2_mu = drho_d2_mu[[d1]][[d2]][[i]], 
                                                                   drho_d2_eta = drho_d2_eta[[d1]][[d2]][[i]])
          
          drho_d2_term[[d1]][[d2]][[i]] <-  crossprod(d2_X[[d1]][[d2]],drho_variance[[i]]*X) + 
            crossprod(d_X[[d1]],drho_d_variance[[d2]][[i]]*X + drho_variance[[i]]*d_X[[d2]]) +
            crossprod(d_X[[d2]],drho_d_variance[[d1]][[i]]*X + drho_variance[[i]]*d_X[[d1]]) +
            crossprod(X,drho_d2_variance[[d1]][[d2]][[i]]*X + drho_d_variance[[d1]][[i]]*d_X[[d2]] +
                        drho_d_variance[[d2]][[i]]*d_X[[d1]] + drho_variance[[i]]*d2_X[[d1]][[d2]]) +
            drho_d2_S[[d1]][[d2]][[i]]
          
          drho_d2_inv_term[[d1]][[d2]][[i]] <- -1*( drho_d_inv_term[[d2]][[i]]  %*% t8[[d1]] + 
                                                      d_inv_term[[d2]]  %*% (drho_d_term[[d1]][[i]] %*% inv_term + d_term[[d1]] %*% drho_inv_term[[i]]) +
                                                      drho_inv_term[[i]] %*% (d2_term[[d1]][[d2]] %*% inv_term + d_term[[d1]] %*% d_inv_term[[d2]]) + 
                                                      inv_term %*% (drho_d2_term[[d1]][[d2]][[i]] %*% inv_term + d2_term[[d1]][[d2]] %*% drho_inv_term[[i]] +
                                                                      drho_d_term[[d1]][[i]] %*% d_inv_term[[d2]] + d_term[[d1]] %*% drho_d_inv_term[[d2]][[i]] )
          )
          
          t22[[d1]][[d2]][[i]] <- crossprod(d_X[[d2]],t19[[d1]][[i]] + t20[[d1]][[i]])
          t23[[d1]][[d2]][[i]] <- crossprod(d2_X[[d1]][[d2]],t11[[i]])
          t38[[d1]][[d2]][[i]] <- crossprod(d_X[[d1]],t26[[d2]][[i]] + t27[[d2]][[i]])
          t39[[d1]][[d2]][[i]] <- crossprod(X,drho_d2_variance[[d1]][[d2]][[i]]*pseudodata + drho_d_variance[[d1]][[i]]*d_pseudodata[[d2]] + 
                                              d2_variance[[d1]][[d2]]*drho_pseudodata[[i]] + d_variance[[d1]]*drho_d_pseudodata[[d2]][[i]] +
                                              drho_d_variance[[d2]][[i]]*d_pseudodata[[d1]] + drho_variance[[i]]*d2_pseudodata[[d1]][[d2]] +
                                              d_variance[[d2]]*drho_d_pseudodata[[d1]][[i]] + variance*drho_d2_pseudodata[[d1]][[d2]][[i]])
          
          drho_d2_coef[[d1]][[d2]][[i]] <- drop( 
            drho_d2_inv_term[[d1]][[d2]][[i]] %*% t1 + 
              drho_d_inv_term[[d1]][[i]] %*% ( t4[[d2]] + t5[[d2]] + t6[[d2]] )  + 
              d2_inv_term[[d1]][[d2]] %*% t2t3[[i]] + 
              d_inv_term[[d1]] %*% ( crossT2[[d2]][[i]] + t37[[d2]][[i]] ) +
              drho_d_inv_term[[d2]][[i]] %*% ( t4[[d1]] + t5[[d1]] + t6[[d1]] ) +
              d_inv_term[[d2]] %*% ( crossT2[[d1]][[i]] + t21[[d1]][[i]] )+
              inv_term %*% ( t23[[d1]][[d2]][[i]] + t38[[d1]][[d2]][[i]] + t22[[d1]][[d2]][[i]] + t39[[d1]][[d2]][[i]] )  +  
              drho_inv_term[[i]] %*% ( crossT1[[d1]][[d2]] + t24[[d1]][[d2]] + t25[[d1]][[d2]] + t30[[d1]][[d2]] )  
          )
        }
      }
      
      if(deriv == 2){
        for(i in 1:length(sp))
        {
          d2rho_d2_eta[[d1]][[d2]][[i]] <- d2rho_d2_mu[[d1]][[d2]][[i]] <- d2rho_d2_variance[[d1]][[d2]][[i]] <- d2rho_d2_pseudodata[[d1]][[d2]][[i]] <- d2rho_d2_inv_term[[d1]][[d2]][[i]] <- d2rho_d2_term[[d1]][[d2]][[i]] <- list()
          for(j in i:length(sp))
          {
            d2rho_d2_eta[[d1]][[d2]][[i]][[j]] <- tcrossprod(d2_X[[d1]][[d2]] , t(d2rho_coef[[i]][[j]])) + 
              tcrossprod(d_X[[d1]] , t(d2rho_d_coef[[d2]][[i]][[j]])) + 
              tcrossprod(d_X[[d2]] , t(d2rho_d_coef[[d1]][[i]][[j]])) + 
              tcrossprod(X , t(d2rho_d2_coef[[d1]][[d2]][[i]][[j]]))
            
            d2rho_d2_mu[[d1]][[d2]][[i]][[j]] <- d2rho_d2_eta[[d1]][[d2]][[i]][[j]] * variance + d2rho_d_eta[[d1]][[i]][[j]] * d_variance[[d2]] + drho_d2_eta[[d1]][[d2]][[i]] * drho_variance[[j]] + drho_d_eta[[d1]][[i]] * drho_d_variance[[d2]][[j]] +
              drho_d2_eta[[d1]][[d2]][[j]] * drho_variance[[i]] + drho_d_eta[[d1]][[j]] * drho_d_variance[[d2]][[i]] + d2_eta[[d1]][[d2]] * d2rho_variance[[i]][[j]] + d_eta[[d1]] * d2rho_d_variance[[d2]][[i]][[j]]
            d2rho_d2_variance[[d1]][[d2]][[i]][[j]] <- as.vector(d2rho_d2_var_fn(family=family, mu=mu, d_mu_1=d_mu[[d1]], d_mu_2=d_mu[[d2]], drho_mu_i=drho_mu[[i]], drho_mu_j=drho_mu[[j]], 
                                                                                 d2_mu=d2_mu[[d1]][[d2]], d2rho_mu=d2rho_mu[[i]][[j]], drho_d_mu_i_1=drho_d_mu[[d1]][[i]], drho_d_mu_i_2=drho_d_mu[[d2]][[i]], drho_d_mu_j_1=drho_d_mu[[d1]][[j]], drho_d_mu_j_2=drho_d_mu[[d2]][[j]],
                                                                                 drho_d2_mu_i=drho_d2_mu[[d1]][[d2]][[i]], drho_d2_mu_j=drho_d2_mu[[d1]][[d2]][[j]], d2rho_d_mu_1=d2rho_d_mu[[d1]][[i]][[j]], d2rho_d_mu_2=d2rho_d_mu[[d2]][[i]][[j]], d2rho_d2_mu=d2rho_d2_mu[[d1]][[d2]][[i]][[j]]))
            d2rho_d2_pseudodata[[d1]][[d2]][[i]][[j]] <- d2rho_d2_pseudo_fn(family = family, y = y, inv_variance = inv_variance, variance = variance, drho_variance_i = drho_variance[[i]], drho_variance_j = drho_variance[[j]], d_variance_1 = d_variance[[d1]], d_variance_2 = d_variance[[d2]], 
                                                                            drho_d_variance_i_1 = drho_d_variance[[d1]][[i]], drho_d_variance_i_2 = drho_d_variance[[d2]][[i]], drho_d_variance_j_1 = drho_d_variance[[d1]][[j]], 
                                                                            drho_d_variance_j_2 = drho_d_variance[[d2]][[j]], d2rho_variance = d2rho_variance[[i]][[j]], d2_variance = d2_variance[[d1]][[d2]], d2rho_d_variance_1 = d2rho_d_variance[[d1]][[i]][[j]], 
                                                                            d2rho_d_variance_2 = d2rho_d_variance[[d2]][[i]][[j]], drho_d2_variance_i = drho_d2_variance[[d1]][[d2]][[i]], drho_d2_variance_j = drho_d2_variance[[d1]][[d2]][[j]],
                                                                            d2rho_d2_variance = d2rho_d2_variance[[d1]][[d2]][[i]][[j]], mu = mu, d_mu_1 = d_mu[[d1]], d_mu_2 = d_mu[[d2]], drho_mu_i = drho_mu[[i]], drho_mu_j = drho_mu[[j]], 
                                                                            d2_mu = d2_mu[[d1]][[d2]], d2rho_mu =d2rho_mu[[i]][[j]], drho_d_mu_i_1 = drho_d_mu[[d1]][[i]], drho_d_mu_i_2 = drho_d_mu[[d2]][[i]], drho_d_mu_j_1 = drho_d_mu[[d1]][[j]], 
                                                                            drho_d_mu_j_2 = drho_d_mu[[d2]][[j]], drho_d2_mu_i = drho_d2_mu[[d1]][[d2]][[i]], drho_d2_mu_j = drho_d2_mu[[d1]][[d2]][[j]], d2rho_d_mu_1 = d2rho_d_mu[[d1]][[i]][[j]], 
                                                                            d2rho_d_mu_2 = d2rho_d_mu[[d2]][[i]][[j]], d2rho_d2_mu = d2rho_d2_mu[[d1]][[d2]][[i]][[j]], d2rho_d2_eta = d2rho_d2_eta[[d1]][[d2]][[i]][[j]])
            
            d2rho_d2_term[[d1]][[d2]][[i]][[j]] <- crossprod(d2_X[[d1]][[d2]],t16[[i]][[j]]) + 
              crossprod(d_X[[d1]],d2rho_d_variance[[d2]][[i]][[j]]*X + d2rho_variance[[i]][[j]]*d_X[[d2]]) +  
              crossprod(d_X[[d2]],d2rho_d_variance[[d1]][[i]][[j]]*X + d2rho_variance[[i]][[j]]*d_X[[d1]]) + 
              crossprod(X,d2rho_d2_variance[[d1]][[d2]][[i]][[j]]*X + d2rho_d_variance[[d1]][[i]][[j]]*d_X[[d2]] +
                          d2rho_d_variance[[d2]][[i]][[j]]*d_X[[d1]] + d2rho_variance[[i]][[j]]*d2_X[[d1]][[d2]]) +  
              d2rho_d2_S[[d1]][[d2]][[i]][[j]]
            
            d2rho_d2_inv_term[[d1]][[d2]][[i]][[j]] <-  -1*( d2rho_d_inv_term[[d2]][[i]][[j]]  %*% t8[[d1]] + 
                                                               
                                                               drho_d_inv_term[[d2]][[i]]  %*% (drho_d_term[[d1]][[j]] %*% inv_term + d_term[[d1]] %*% drho_inv_term[[j]]) + 
                                                               drho_d_inv_term[[d2]][[j]]  %*% (drho_d_term[[d1]][[i]] %*% inv_term + d_term[[d1]] %*% drho_inv_term[[i]]) + 
                                                               d_inv_term[[d2]]  %*% (d2rho_d_term[[d1]][[i]][[j]] %*% inv_term +
                                                                                        drho_d_term[[d1]][[i]] %*% drho_inv_term[[j]] +
                                                                                        drho_d_term[[d1]][[j]] %*% drho_inv_term[[i]] +
                                                                                        d_term[[d1]] %*% d2rho_inv_term[[i]][[j]]
                                                               ) + 
                                                               d2rho_inv_term[[i]][[j]] %*% (d2_term[[d1]][[d2]] %*% inv_term +
                                                                                               d_term[[d1]] %*% d_inv_term[[d2]]
                                                               ) + 
                                                               drho_inv_term[[i]] %*% (drho_d2_term[[d1]][[d2]][[j]] %*% inv_term +
                                                                                         d2_term[[d1]][[d2]] %*% drho_inv_term[[j]] +
                                                                                         drho_d_term[[d1]][[j]] %*% d_inv_term[[d2]] +
                                                                                         d_term[[d1]] %*% drho_d_inv_term[[d2]][[j]]
                                                               ) + 
                                                               drho_inv_term[[j]] %*% (drho_d2_term[[d1]][[d2]][[i]] %*% inv_term +
                                                                                         d2_term[[d1]][[d2]] %*% drho_inv_term[[i]] +
                                                                                         drho_d_term[[d1]][[i]] %*% d_inv_term[[d2]] +
                                                                                         d_term[[d1]] %*% drho_d_inv_term[[d2]][[i]]
                                                               ) +
                                                               inv_term %*% (d2rho_d2_term[[d1]][[d2]][[i]][[j]] %*% inv_term +
                                                                               drho_d2_term[[d1]][[d2]][[i]] %*% drho_inv_term[[j]] +
                                                                               drho_d2_term[[d1]][[d2]][[j]] %*% drho_inv_term[[i]] +
                                                                               d2_term[[d1]][[d2]] %*% d2rho_inv_term[[i]][[j]] +
                                                                               d2rho_d_term[[d1]][[i]][[j]] %*% d_inv_term[[d2]] +
                                                                               drho_d_term[[d1]][[i]] %*% drho_d_inv_term[[d2]][[j]] +
                                                                               drho_d_term[[d1]][[j]] %*% drho_d_inv_term[[d2]][[i]] +
                                                                               d_term[[d1]] %*% d2rho_d_inv_term[[d2]][[i]][[j]]
                                                               )  
            )
            d2rho_d2_coef[[d1]][[d2]][[i]][[j]] <- drop( 
              d2rho_d2_inv_term[[d1]][[d2]][[i]][[j]] %*% t1 + 
                d2rho_d_inv_term[[d1]][[i]][[j]] %*% ( t4[[d2]] + t5[[d2]] + t6[[d2]] )  +     
                drho_d2_inv_term[[d1]][[d2]][[j]] %*% t2t3[[i]] +
                drho_d2_inv_term[[d1]][[d2]][[i]] %*% t2t3[[j]] +  
                drho_d_inv_term[[d1]][[i]] %*% ( crossT2[[d2]][[j]] + t21[[d2]][[j]] ) + 
                d2_inv_term[[d1]][[d2]] %*% t41[[i]][[j]] +   
                drho_d_inv_term[[d1]][[j]] %*% ( crossT2[[d2]][[i]] + t37[[d2]][[i]] ) + 
                d_inv_term[[d1]] %*% ( t40[[d2]][[i]][[j]] + t35[[d2]][[i]][[j]]) +
                d2rho_d_inv_term[[d2]][[i]][[j]] %*% ( t4[[d1]] + t5[[d1]] + t6[[d1]] ) +
                drho_d_inv_term[[d2]][[i]] %*% ( crossT2[[d1]][[j]] + t21[[d1]][[j]] ) + 
                drho_d_inv_term[[d2]][[j]] %*% ( crossT2[[d1]][[i]] + t21[[d1]][[i]] ) +
                d_inv_term[[d2]] %*% ( t40[[d1]][[i]][[j]] + t35[[d1]][[i]][[j]]) +
                drho_inv_term[[j]] %*% ( t23[[d1]][[d2]][[i]] + t38[[d1]][[d2]][[i]] + t22[[d1]][[d2]][[i]] + t39[[d1]][[d2]][[i]] ) +
                inv_term %*% ( crossprod(d2_X[[d1]][[d2]],regularExp[[i]][[j]]) + 
                                 crossprod(d_X[[d1]],t36[[d2]][[i]][[j]] ) +
                                 crossprod(d_X[[d2]],t36[[d1]][[i]][[j]] ) + 
                                 crossprod(X,d2rho_d2_variance[[d1]][[d2]][[i]][[j]]*pseudodata + drho_d2_variance[[d1]][[d2]][[i]]*drho_pseudodata[[j]] +
                                             d2rho_d_variance[[d1]][[i]][[j]]*d_pseudodata[[d2]] + drho_d_variance[[d1]][[i]]*drho_d_pseudodata[[d2]][[j]] +
                                             drho_d2_variance[[d1]][[d2]][[j]]*drho_pseudodata[[i]] + d2_variance[[d1]][[d2]]*d2rho_pseudodata[[i]][[j]] +
                                             drho_d_variance[[d1]][[j]]*drho_d_pseudodata[[d2]][[i]] + d_variance[[d1]]*d2rho_d_pseudodata[[d2]][[i]][[j]] +
                                             d2rho_d_variance[[d2]][[i]][[j]]*d_pseudodata[[d1]] + drho_d_variance[[d2]][[i]]*drho_d_pseudodata[[d1]][[j]] +
                                             d2rho_variance[[i]][[j]]*d2_pseudodata[[d1]][[d2]] + drho_variance[[i]]*drho_d2_pseudodata[[d1]][[d2]][[j]] +
                                             drho_d_variance[[d2]][[j]]*drho_d_pseudodata[[d1]][[i]] + d_variance[[d2]]*d2rho_d_pseudodata[[d1]][[i]][[j]] +
                                             drho_variance[[j]]*drho_d2_pseudodata[[d1]][[d2]][[i]] + variance*d2rho_d2_pseudodata[[d1]][[d2]][[i]][[j]] ) ) +
                d2rho_inv_term[[i]][[j]] %*% ( crossT1[[d1]][[d2]] + t24[[d1]][[d2]] + t25[[d1]][[d2]] + t30[[d1]][[d2]] ) +
                drho_inv_term[[i]] %*% ( t23[[d1]][[d2]][[j]] + crossprod(d_X[[d1]],t19[[d2]][[j]] + t20[[d2]][[j]]) + t22[[d1]][[d2]][[j]] + 
                                           crossprod(X,drho_d2_variance[[d1]][[d2]][[j]]*pseudodata + d2_variance[[d1]][[d2]]*drho_pseudodata[[j]] +
                                                       drho_d_variance[[d1]][[j]]*d_pseudodata[[d2]] + d_variance[[d1]]*drho_d_pseudodata[[d2]][[j]] +
                                                       drho_d_variance[[d2]][[j]]*d_pseudodata[[d1]] + d_variance[[d2]]*drho_d_pseudodata[[d1]][[j]] +
                                                       drho_variance[[j]]*d2_pseudodata[[d1]][[d2]] + variance*drho_d2_pseudodata[[d1]][[d2]][[j]] ) )   
            )
          }
        }
      }
    }
  }
  
  for(d1 in 1:length(beta_star))
  { 
    for(i in 1:length(sp))
    {
      if(deriv == 2){
        for(j in i:length(sp))
        {
          d2rho_d_coef[[d1]][[i]][[j]] <- drop( 
            d2rho_d_inv_term[[d1]][[i]][[j]] %*% t1 + 
              drho_d_inv_term[[d1]][[i]] %*% crossprod(X,t11[[j]]) +
              drho_d_inv_term[[d1]][[j]] %*% crossprod(X,t11[[i]]) + 
              d_inv_term[[d1]] %*% t41[[i]][[j]] +
              d2rho_inv_term[[i]][[j]] %*% (t4[[d1]] + crossprod(X,t9[[d1]])) + 
              drho_inv_term[[i]] %*% (crossT2[[d1]][[j]] + t21[[d1]][[j]] ) +
              drho_inv_term[[j]] %*% (crossT2[[d1]][[i]] + t21[[d1]][[i]] ) +
              inv_term %*% (t40[[d1]][[i]][[j]] + t35[[d1]][[i]][[j]]) 
          )
        }
      }
      
      if(deriv){
        drho_d_coef[[d1]][[i]] <- drop( 
          drho_d_inv_term[[d1]][[i]] %*% t1 + d_inv_term[[d1]] %*% t2t3[[i]] + 
            drho_inv_term[[i]] %*% (t4[[d1]] + t5[[d1]] + t6[[d1]]) + inv_term %*% (crossT2[[d1]][[i]] + t21[[d1]][[i]]) 
        )
      }
    }
    d_coef[[d1]] <- drop( d_inv_term[[d1]] %*% t1 + inv_term %*% (t4[[d1]] + t5[[d1]] + t6[[d1]]) )
  }
  
  for(i in 1:length(sp))
  {
    if(deriv == 2){
      for(j in i:length(sp))
      {
        d2rho_coef[[i]][[j]] <- drop(d2rho_inv_term[[i]][[j]] %*% t1 + drho_inv_term[[i]] %*% t2t3[[j]] + drho_inv_term[[j]] %*% t2t3[[i]] + inv_term %*% t41[[i]][[j]])
      }
    }
    if(deriv){
      drho_coef[[i]] <- drop(drho_inv_term[[i]] %*% t1 + inv_term %*% t2t3[[i]])
    }
  }
  
  return(list(drho_coef = drho_coef, d2rho_coef = d2rho_coef, d_coef = d_coef, drho_d_coef = drho_d_coef, 
              d2rho_d_coef = d2rho_d_coef, d2_coef = d2_coef, drho_d2_coef = drho_d2_coef, d2rho_d2_coef = d2rho_d2_coef))
}

#08f
d_gdi1 <- function(deriv, family, scoreType,  beta_star, good, y, X, d_X, d2_X, mu, variance, 
                   coef, drho_coef, d2rho_coef, d_coef, drho_d_coef, d2rho_d_coef, d2_coef, drho_d2_coef, d2rho_d2_coef, 
                   dev, sp, d_sp, d2_sp, St,d_St, d2_St, drho_S, d2rho_S, drho_d_S, drho_d2_S, d2rho_d_S, d2rho_d2_S){
  
  d_eta <- d_mu <- d_variance <- d_term <- d_inv_term <- drho_d_eta <- drho_d_mu <- drho_d_variance <- drho_d_term <- drho_d_inv_term <- d2rho_d_eta <- d2rho_d_mu <- d2rho_d_variance <- d2rho_d_term <- d2rho_d_inv_term <- list()
  d2_eta <- d2_mu <- d2_variance <- d2_term <- d2_inv_term <- drho_d2_eta <- drho_d2_mu <- drho_d2_variance <- drho_d2_term <- drho_d2_inv_term <- d2rho_d2_eta <- d2rho_d2_mu <- d2rho_d2_variance <- d2rho_d2_term <- d2rho_d2_inv_term <- list()
  drho_d2_dev <- d2rho_d2_dev <- d2_dev <- d2_hat <- d2_rank.tol <- d2_conv.tol <- matrix(nrow=length(beta_star),ncol=length(beta_star))
  drho_d_dev <- d2rho_d_dev <- drho_d2_dev <- d2rho_d2_dev <- d2rho_eta <- d2rho_mu <- d2rho_variance <- d2rho_term <- d2rho_inv_term <- drho_d_hat <- d2rho_d_hat <- drho_d2_hat <- d2rho_d2_hat <- drho_eta <- drho_mu <- drho_variance <- drho_term <- drho_inv_term <- list()
  d_dev <- d_hat <- vector(length=length(beta_star))
  
  
  
  n <- nrow(X)
  
  singularityTest <- tryCatch(
    {
      test <- solve(crossprod(X * variance,X) + St)
      out <- "non-singular"
    },
    warning = function(war) {
      out <- "singular"
    },
    error = function(err){
      out <- "singular"
    }
  )
  
  if(singularityTest == "singular"){
    m <- ncol(X)
    colNames <- colnames(X) <- as.character(1:ncol(X))
    X_Decomp <- qr(X, complete = FALSE)
    drop <- X_Decomp$pivot[-c(1:X_Decomp$rank)]
    
    
    if(length(drop) > 0){
      X <- X[,!colnames(X) %in% drop]
      #We need to drop unidentifiable columns from a number of objects
      d_X <- sapply(d_X, function(d_X){
        colnames(d_X) <- 1:ncol(d_X)
        d_X <- d_X[,!colnames(d_X) %in% drop]
      },simplify = FALSE)
      d2_X <- sapply(1:length(d2_X), function(i){
        out <- sapply(d2_X[[i]], function(d2_X){
          if(ncol(d2_X)>0)
          {
            colnames(d2_X) <- 1:ncol(d2_X)
            d2_X <- d2_X[,!colnames(d2_X) %in% drop]
          }
        },simplify = FALSE)
        return(out)
      },simplify = FALSE)
    }
    
  }
  
  term <- crossprod(X * variance,X)  + St
  inv_term <- solve( term )
  
  if(deriv){
    drho_eta <- lapply(drho_coef,FUN= function(input) X %*% input)
    drho_mu <- lapply(drho_eta, "*", variance)
    drho_variance <- drho_var_fn(family=family, mu, drho_mu)
    temp <- lapply(drho_variance, FUN=function(input) crossprod(X*input ,X)  )
    drho_term <- mapply("+" ,temp ,drho_S,SIMPLIFY =F)
    drho_inv_term <- lapply(drho_term, FUN=function(input) -1*inv_term %*% input %*% inv_term  )
  }
  
  if(deriv == 2){
    t16 <- list()
    for(i in 1:length(sp))
    {
      d2rho_eta[[i]] <- d2rho_mu[[i]] <- d2rho_variance[[i]] <- d2rho_term[[i]] <- d2rho_inv_term[[i]] <- t16[[i]] <- list()
      for(j in i:length(sp))
      {
        
        d2rho_eta[[i]][[j]] <- tcrossprod(X,t(d2rho_coef[[i]][[j]]))
        d2rho_mu[[i]][[j]] <- d2rho_eta[[i]][[j]]*variance + drho_eta[[i]]*drho_variance[[j]] 
        d2rho_variance[[i]][[j]] <- as.vector(d2rho_var_fn(family=family, mu=mu, drho_mu_i=drho_mu[[i]], drho_mu_j=drho_mu[[j]], d2rho_mu=d2rho_mu[[i]][[j]]))
        
        t16[[i]][[j]] <- d2rho_variance[[i]][[j]]*X
        
        d2rho_term[[i]][[j]] <- crossprod(X,t16[[i]][[j]]) + d2rho_S[[i]][[j]]
        d2rho_inv_term[[i]][[j]] <- -1*( 
          drho_inv_term[[j]] %*% drho_term[[i]] %*% inv_term + 
            inv_term %*% (d2rho_term[[i]][[j]] %*% inv_term + drho_term[[i]] %*% drho_inv_term[[j]])
        )
      }
    }
  }
  
  if(sum(good) != n){
    d_X <- lapply(d_X, FUN=function(input) input[good, ])
  }
  d_eta <- mapply("+",lapply(d_X,FUN= function(input) input %*% coef), lapply(d_coef, FUN=function(input) X %*% input ) ,SIMPLIFY =F)
  d_mu <- lapply(d_eta, "*", variance)
  d_variance <- d_var_fn(family=family, mu, d_mu)
  temp1 <- lapply(d_X, FUN=function(input) crossprod(input,variance*X)  )
  temp2 <- lapply(d_variance, FUN=function(input) crossprod(X ,input*X)  )
  temp3 <- lapply(d_X, FUN=function(input) crossprod(X ,variance*input)  )
  d_term <- mapply(FUN=function(input1,input2,input3,input4) input1+input2+input3+input4 ,temp1 ,temp2 ,temp3 ,d_St ,SIMPLIFY =F)
  d_inv_term <- lapply(d_term, FUN=function(input) -1*inv_term %*% input %*% inv_term  )
  d_dev <- d_deviance_fn(family=family, y, mu, d_mu)
  
  if (scoreType %in% c("REML", "ML")){
    d_hat <- rep(NA,times=length(d_mu))
    d_conv.tol <- mapply(FUN=function(d_St, d_coefficients) { crossprod(d_coefficients,St) %*% coef + crossprod(coef,d_St) %*% coef + crossprod(coef,St) %*% d_coefficients}, 
                         d_St=d_St, d_coefficients=d_coef)
    d_rank.tol <- unlist(lapply(d_term, FUN=function(input) sum(diag( inv_term %*% input )) )) 
  }
  else
  {
    temp1 <- lapply(d_variance, FUN=function(input) rowSums((input*X %*% inv_term)*X) )
    temp2 <- lapply(d_X, FUN=function(input) rowSums( (variance*input %*% inv_term)*X + variance*(X %*%inv_term)*input) ) 
    temp3 <- lapply(d_inv_term, FUN=function(input) rowSums( (variance*X %*% input)*X) )
    d_hat <- unlist(mapply(FUN=function(input1,input2,input3) sum(input1+input2+input3), temp1,temp2,temp3,SIMPLIFY =F)) 
    d_conv.tol <- rep(NA,times=length(d_mu))
    d_rank.tol <- rep(NA,times=length(d_mu))
  }
  
  t8 <- lapply(d_term, FUN=function(input) input %*% inv_term)
  drho_d_rank.tol <- drho_d_conv.tol <- d2rho_d_rank.tol <- d2rho_d_conv.tol <- list()
  for(d1 in 1:length(beta_star))
  {
    if(deriv){
      drho_d_eta[[d1]] <- drho_d_mu[[d1]] <- drho_d_variance[[d1]] <- drho_d_term[[d1]] <- drho_d_inv_term[[d1]] <- list()
      drho_d_rank.tol[[d1]] <- drho_d_conv.tol[[d1]] <- drho_d_hat[[d1]] <- drho_d_dev[[d1]] <- vector(length=length(sp))
      
      for(i in 1:length(sp))
      {
        drho_d_eta[[d1]][[i]] <- drop( tcrossprod(d_X[[d1]] , t(drho_coef[[i]]) )  + tcrossprod(X , t(drho_d_coef[[d1]][[i]])) )
        drho_d_mu[[d1]][[i]] <- drho_d_eta[[d1]][[i]]*variance + d_eta[[d1]]*drho_variance[[i]]
        drho_d_variance[[d1]][[i]] <- as.vector(drho_d_var_fn(family=family, mu=mu, d_mu=d_mu[[d1]], drho_mu=drho_mu[[i]], drho_d_mu=drho_d_mu[[d1]][[i]]))
        
        drho_d_term[[d1]][[i]] <- crossprod((d_X[[d1]]),drho_variance[[i]]*X) + crossprod(X,drho_d_variance[[d1]][[i]]*X + drho_variance[[i]]*d_X[[d1]]) + drho_d_S[[d1]][[i]]       
        drho_d_inv_term[[d1]][[i]] <- -1*( drho_inv_term[[i]] %*% (t8[[d1]]) + 
                                             inv_term %*% (drho_d_term[[d1]][[i]] %*% inv_term) + 
                                             inv_term %*% (d_term[[d1]] %*% drho_inv_term[[i]]) )
        drho_d_dev[[d1]][i] <- drho_d_deviance_fn(family=family, y=y, mu=mu, d_mu=d_mu[[d1]], drho_mu=drho_mu[[i]], drho_d_mu=drho_d_mu[[d1]][[i]])
        
        if (scoreType %in% c("REML", "ML")){
          drho_d_hat[[d1]][i] <- NA
          drho_d_rank.tol[[d1]][i] <- sum(diag( d_inv_term[[d1]] %*% drho_term[[i]] + inv_term %*% drho_d_term[[d1]][[i]] )) 
          
          drho_d_conv.tol[[d1]][i] <- crossprod(drho_d_coef[[d1]][[i]],St) %*% coef + 
            crossprod(d_coef[[d1]],drho_S[[i]]) %*% coef + 
            crossprod(d_coef[[d1]],St) %*% drho_coef[[i]] + 
            crossprod(drho_coef[[i]],d_St[[d1]]) %*% coef + 
            crossprod(coef,drho_d_S[[d1]][[i]]) %*% coef + 
            crossprod(coef,d_St[[d1]]) %*% drho_coef[[i]] + 
            crossprod(drho_coef[[i]],St) %*% d_coef[[d1]] + 
            crossprod(coef,drho_S[[i]]) %*% d_coef[[d1]] + 
            crossprod(coef,St) %*% drho_d_coef[[d1]][[i]]
        }
        else
        {
          drho_d_hat[[d1]][i] <- sum(rowSums( 
            (drho_d_variance[[d1]][[i]]*X %*% inv_term)*X + (drho_variance[[i]]* d_X[[d1]] %*% inv_term)*X + (drho_variance[[i]]* X %*% d_inv_term[[d1]])*X + 
              (drho_variance[[i]]* X %*% inv_term)*d_X[[d1]] + (d_variance[[d1]]*X %*% drho_inv_term[[i]])*X + (variance*d_X[[d1]] %*% drho_inv_term[[i]])*X +
              (variance*X %*% drho_d_inv_term[[d1]][[i]])*X + (variance*X %*% drho_inv_term[[i]])*d_X[[d1]] 
          ))
          drho_d_rank.tol[[d1]][i] <- NA
          drho_d_conv.tol[[d1]][i] <- NA
        }
      }
    }
    
    if(deriv == 2){
      d2rho_d_eta[[d1]] <- d2rho_d_mu[[d1]] <- d2rho_d_variance[[d1]] <- d2rho_d_term[[d1]] <- d2rho_d_inv_term[[d1]] <- list()
      d2rho_d_hat[[d1]] <- matrix(nrow=length(sp),ncol=length(sp))
      d2rho_d_dev[[d1]] <- matrix(nrow=length(sp),ncol=length(sp))
      d2rho_d_rank.tol[[d1]] <- matrix(nrow=length(sp),ncol=length(sp))
      d2rho_d_conv.tol[[d1]] <- matrix(nrow=length(sp),ncol=length(sp))
      for(i in 1:length(sp))
      {
        d2rho_d_eta[[d1]][[i]] <- d2rho_d_mu[[d1]][[i]] <- d2rho_d_variance[[d1]][[i]] <- d2rho_d_term[[d1]][[i]] <- d2rho_d_inv_term[[d1]][[i]] <- list()
        for(j in i:length(sp))
        {
          d2rho_d_eta[[d1]][[i]][[j]] <- drop( tcrossprod(d_X[[d1]] , t(d2rho_coef[[i]][[j]])) + tcrossprod(X , t(d2rho_d_coef[[d1]][[i]][[j]])) )
          d2rho_d_mu[[d1]][[i]][[j]] <- d2rho_d_eta[[d1]][[i]][[j]] * variance + drho_d_eta[[d1]][[i]] * drho_variance[[j]] + drho_d_eta[[d1]][[j]] * drho_variance[[i]] + d_eta[[d1]] * d2rho_variance[[i]][[j]]
          d2rho_d_variance[[d1]][[i]][[j]] <- as.vector(d2rho_d_var_fn(family=family, mu, d_mu=d_mu[[d1]], drho_mu_i=drho_mu[[i]], drho_mu_j=drho_mu[[j]], 
                                                                       drho_d_mu_i=drho_d_mu[[d1]][[i]], drho_d_mu_j=drho_d_mu[[d1]][[j]],d2rho_mu=d2rho_mu[[i]][[j]], d2rho_d_mu=d2rho_d_mu[[d1]][[i]][[j]]))
          
          d2rho_d_term[[d1]][[i]][[j]] <- crossprod(d_X[[d1]] ,t16[[i]][[j]]) + crossprod(X ,d2rho_d_variance[[d1]][[i]][[j]]*X + d2rho_variance[[i]][[j]]*d_X[[d1]]) + d2rho_d_S[[d1]][[i]][[j]]
          
          d2rho_d_inv_term[[d1]][[i]][[j]] <- -1*( 
            d2rho_inv_term[[i]][[j]] %*% t8[[d1]] + 
              drho_inv_term[[i]] %*% (drho_d_term[[d1]][[j]] %*% inv_term + d_term[[d1]] %*% drho_inv_term[[j]]) + 
              
              drho_inv_term[[j]] %*% (drho_d_term[[d1]][[i]] %*% inv_term + d_term[[d1]] %*% drho_inv_term[[i]]) + 
              
              inv_term %*% (d2rho_d_term[[d1]][[i]][[j]] %*% inv_term + drho_d_term[[d1]][[i]] %*% drho_inv_term[[j]] +
                              drho_d_term[[d1]][[j]] %*% drho_inv_term[[i]] + d_term[[d1]] %*% d2rho_inv_term[[i]][[j]]
              ) )
          
          d2rho_d_dev[[d1]][i,j] <- d2rho_d_deviance_fn(family=family, y=y, mu=mu, d_mu=d_mu[[d1]], drho_mu_i=drho_mu[[i]], drho_mu_j=drho_mu[[j]], d2rho_mu= d2rho_mu[[i]][[j]],
                                                        drho_d_mu_i=drho_d_mu[[d1]][[i]], drho_d_mu_j=drho_d_mu[[d1]][[j]], d2rho_d_mu=d2rho_d_mu[[d1]][[i]][[j]])
          
          if (scoreType %in% c("REML", "ML")){
            d2rho_d_hat[[d1]][i,j] <- NA
            d2rho_d_rank.tol[[d1]][i,j] <- sum(diag( 
              drho_d_inv_term[[d1]][[j]] %*% drho_term[[i]] + d_inv_term[[d1]] %*% d2rho_term[[i]][[j]] +  drho_inv_term[[j]] %*% drho_d_term[[d1]][[i]] + inv_term %*% d2rho_d_term[[d1]][[i]][[j]] 
            ))
            d2rho_d_conv.tol[[d1]][i,j] <- crossprod(d2rho_d_coef[[d1]][[i]][[j]],St) %*% coef + 
              crossprod(drho_d_coef[[d1]][[i]],drho_S[[j]]) %*% coef + 
              crossprod(drho_d_coef[[d1]][[i]],St) %*% drho_coef[[j]] + 
              crossprod(drho_d_coef[[d1]][[j]],drho_S[[i]]) %*% coef + 
              crossprod(d_coef[[d1]],d2rho_S[[i]][[j]]) %*% coef + 
              crossprod(d_coef[[d1]],drho_S[[i]]) %*% drho_coef[[j]] + 
              crossprod(drho_d_coef[[d1]][[j]],St) %*% drho_coef[[i]] + 
              crossprod(d_coef[[d1]],drho_S[[j]]) %*% drho_coef[[i]] + 
              crossprod(d_coef[[d1]],St) %*% d2rho_coef[[i]][[j]] + 
              crossprod(d2rho_coef[[i]][[j]],d_St[[d1]]) %*% coef + 
              crossprod(drho_coef[[i]],drho_d_S[[d1]][[j]]) %*% coef + 
              crossprod(drho_coef[[i]],d_St[[d1]]) %*% drho_coef[[j]] + 
              crossprod(drho_coef[[j]],drho_d_S[[d1]][[i]]) %*% coef + 
              crossprod(coef,d2rho_d_S[[d1]][[i]][[j]]) %*% coef + 
              crossprod(coef,drho_d_S[[d1]][[i]]) %*% drho_coef[[j]] + 
              crossprod(drho_coef[[j]],d_St[[d1]]) %*% drho_coef[[i]] + 
              crossprod(coef,drho_d_S[[d1]][[j]]) %*% drho_coef[[i]] + 
              crossprod(coef,d_St[[d1]]) %*% d2rho_coef[[i]][[j]] + 
              crossprod(d2rho_coef[[i]][[j]],St) %*% d_coef[[d1]] + 
              crossprod(drho_coef[[i]],drho_S[[j]]) %*% d_coef[[d1]] + 
              crossprod(drho_coef[[i]],St) %*% drho_d_coef[[d1]][[j]] + 
              crossprod(drho_coef[[j]],drho_S[[i]]) %*% d_coef[[d1]] + 
              crossprod(coef,d2rho_S[[i]][[j]]) %*% d_coef[[d1]] + 
              crossprod(coef,drho_S[[i]]) %*% drho_d_coef[[d1]][[j]] + 
              crossprod(drho_coef[[j]],St) %*% drho_d_coef[[d1]][[i]] + 
              crossprod(coef,drho_S[[j]]) %*% drho_d_coef[[d1]][[i]] + 
              crossprod(coef,St) %*% d2rho_d_coef[[d1]][[i]][[j]]
          }
          else{
            d2rho_d_hat[[d1]][i,j] <- sum(rowSums(
              (d2rho_d_variance[[d1]][[i]][[j]]*X %*% inv_term)*X + (drho_d_variance[[d1]][[i]]*X %*% drho_inv_term[[j]])*X +
                (d2rho_variance[[i]][[j]]*d_X[[d1]] %*% inv_term)*X + (d2rho_variance[[i]][[j]]*X %*% d_inv_term[[d1]])*X +
                (d2rho_variance[[i]][[j]]*X %*% inv_term)*d_X[[d1]] + (drho_variance[[i]]*d_X[[d1]] %*% drho_inv_term[[j]])*X + 
                (drho_variance[[i]]*X %*% drho_d_inv_term[[d1]][[j]])*X + (drho_variance[[i]]*X %*% drho_inv_term[[j]])*d_X[[d1]] +
                (drho_d_variance[[d1]][[j]]*X %*% drho_inv_term[[i]])*X + (d_variance[[d1]]*X %*% d2rho_inv_term[[i]][[j]])*X +
                (drho_variance[[j]]*d_X[[d1]] %*% drho_inv_term[[i]])*X + (drho_variance[[j]]*X %*% drho_d_inv_term[[d1]][[i]])*X +
                (drho_variance[[j]]*X %*% drho_inv_term[[i]])*d_X[[d1]] + (variance*d_X[[d1]] %*% d2rho_inv_term[[i]][[j]])*X + 
                (variance*X %*% d2rho_d_inv_term[[d1]][[i]][[j]])*X + (variance*X %*% d2rho_inv_term[[i]][[j]])*d_X[[d1]] 
            ))
            d2rho_d_rank.tol[[d1]][i,j] <- NA
            d2rho_d_conv.tol[[d1]][i,j] <- NA
          }
        }
      }
    }
  }
  
  drho_d2_rank.tol <- drho_d2_conv.tol <- d2rho_d2_rank.tol <- d2rho_d2_conv.tol <- list()
  for(d1 in 1:length(beta_star))
  { 
    d2_eta[[d1]] <- d2_mu[[d1]] <- d2_variance[[d1]] <- d2_term[[d1]] <- d2_inv_term[[d1]] <- drho_d2_eta[[d1]] <- drho_d2_mu[[d1]] <- drho_d2_variance[[d1]] <- drho_d2_term[[d1]] <- drho_d2_inv_term[[d1]] <- d2rho_d2_eta[[d1]] <- d2rho_d2_mu[[d1]] <- d2rho_d2_variance[[d1]] <- d2rho_d2_term[[d1]] <- d2rho_d2_inv_term[[d1]] <- list()
    drho_d2_dev[[d1]] <- d2rho_d2_dev[[d1]] <- drho_d2_rank.tol[[d1]] <- drho_d2_conv.tol[[d1]] <- d2rho_d2_rank.tol[[d1]] <- d2rho_d2_conv.tol[[d1]] <- drho_d2_hat[[d1]] <- d2rho_d2_hat[[d1]] <- list() 
    for(d2 in d1:length(beta_star))
    {  
      if(sum(good) != n){
        d2_X[[d1]][[d2]] <- d2_X[[d1]][[d2]][good, ]     
      }
      
      d2_eta[[d1]][[d2]] <- ( tcrossprod(d2_X[[d1]][[d2]] , t(coef)) + tcrossprod(d_X[[d1]] , t(d_coef[[d2]])) + 
                                tcrossprod(d_X[[d2]] , t(d_coef[[d1]])) + tcrossprod(X , t(d2_coef[[d1]][[d2]])) )
      
      d2_mu[[d1]][[d2]] <- d2_eta[[d1]][[d2]]*variance + d_eta[[d1]]*d_variance[[d2]]
      d2_variance[[d1]][[d2]] <- as.vector(d2_var_fn(family=family, mu=mu, d_mu_1=d_mu[[d1]], d_mu_2=d_mu[[d2]], d2_mu = d2_mu[[d1]][[d2]]))
      
      d2_term[[d1]][[d2]] <- crossprod(d2_X[[d1]][[d2]],variance*X) + crossprod(d_X[[d1]],d_variance[[d2]]*X + variance*d_X[[d2]]) + 
        crossprod(d_X[[d2]],d_variance[[d1]]*X + variance*d_X[[d1]]) + 
        crossprod(X,d2_variance[[d1]][[d2]]*X + d_variance[[d1]]*d_X[[d2]] + d_variance[[d2]]*d_X[[d1]] + 
                    variance*d2_X[[d1]][[d2]]) + d2_St[[d1]][[d2]]
      
      d2_inv_term[[d1]][[d2]] <- -1*(d_inv_term[[d2]]  %*% (t8[[d1]]) + inv_term %*% (d2_term[[d1]][[d2]] %*% inv_term) + inv_term %*% (d_term[[d1]] %*% d_inv_term[[d2]])) 
      
      d2_dev[d1,d2] <- d2_deviance_fn(family = family, y=y, mu=mu, d_mu_1=d_mu[[d1]], d_mu_2=d_mu[[d2]], d2_mu = d2_mu[[d1]][[d2]])
      
      if (scoreType %in% c("REML", "ML")){
        d2_hat[d1,d2] <- NA
        d2_rank.tol[d1,d2] <- sum(diag( d_inv_term[[d2]] %*% d_term[[d1]] + inv_term %*% d2_term[[d1]][[d2]] ))
        d2_conv.tol[d1,d2] <- crossprod(d2_coef[[d1]][[d2]],St) %*% coef + 
          crossprod(d_coef[[d1]],d_St[[d2]]) %*% coef + 
          crossprod(d_coef[[d1]],St) %*% d_coef[[d2]] + 
          crossprod(d_coef[[d2]],d_St[[d1]]) %*% coef + 
          crossprod(coef,d2_St[[d1]][[d2]]) %*% coef + 
          crossprod(coef,d_St[[d1]]) %*% d_coef[[d2]] + 
          crossprod(d_coef[[d2]],St) %*% d_coef[[d1]] + 
          crossprod(coef,d_St[[d2]]) %*% d_coef[[d1]] + 
          crossprod(coef,St) %*% d2_coef[[d1]][[d2]]
      }
      else{
        d2_hat[d1,d2] <- sum(rowSums(
          (d2_variance[[d1]][[d2]]*X %*% inv_term)*X  + (d_variance[[d1]]*d_X[[d2]] %*% inv_term)*X + (d_variance[[d1]]*X %*% d_inv_term[[d2]])*X + (d_variance[[d1]]*X %*% inv_term)*d_X[[d2]] +
            (d_variance[[d2]]*d_X[[d1]] %*% inv_term)*X + (variance*d2_X[[d1]][[d2]] %*% inv_term)*X + (variance*d_X[[d1]] %*% d_inv_term[[d2]])*X + (variance*d_X[[d1]] %*% inv_term)*d_X[[d2]] + 
            (d_variance[[d2]]*X %*% d_inv_term[[d1]])*X + (variance*d_X[[d2]] %*% d_inv_term[[d1]])*X + (variance*X %*% d2_inv_term[[d1]][[d2]])*X + (variance*X %*% d_inv_term[[d1]])*d_X[[d2]] +
            (d_variance[[d2]]*X %*% inv_term)*d_X[[d1]] + (variance*d_X[[d2]] %*% inv_term)*d_X[[d1]] + (variance*X %*% d_inv_term[[d2]])*d_X[[d1]] + (variance*X %*% inv_term)*d2_X[[d1]][[d2]] 
        ))
        d2_rank.tol[d1,d2] <- NA
        d2_conv.tol[d1,d2] <- NA
      }
      
      if(deriv){
        drho_d2_eta[[d1]][[d2]] <- drho_d2_mu[[d1]][[d2]] <- drho_d2_variance[[d1]][[d2]] <- drho_d2_term[[d1]][[d2]] <- drho_d2_inv_term[[d1]][[d2]] <- list()
        drho_d2_dev[[d1]][[d2]] <- drho_d2_hat[[d1]][[d2]] <- drho_d2_rank.tol[[d1]][[d2]] <- drho_d2_conv.tol[[d1]][[d2]] <- vector(length=length(sp))
        for(i in 1:length(sp))
        {
          drho_d2_eta[[d1]][[d2]][[i]] <- ( tcrossprod(d2_X[[d1]][[d2]] , t(drho_coef[[i]])) + tcrossprod(d_X[[d1]] , t(drho_d_coef[[d2]][[i]])) + 
                                              tcrossprod(d_X[[d2]] , t(drho_d_coef[[d1]][[i]])) + tcrossprod(X , t(drho_d2_coef[[d1]][[d2]][[i]])) )
          
          drho_d2_mu[[d1]][[d2]][[i]] <- drho_d2_eta[[d1]][[d2]][[i]] * variance + drho_d_eta[[d1]][[i]] * d_variance[[d2]] + d2_eta[[d1]][[d2]] * drho_variance[[i]] + d_eta[[d1]] * drho_d_variance[[d2]][[i]]
          drho_d2_variance[[d1]][[d2]][[i]] <- as.vector(drho_d2_var_fn(family=family, mu, d_mu_1=d_mu[[d1]], d_mu_2=d_mu[[d2]], drho_mu=drho_mu[[i]], 
                                                                        drho_d_mu_1=drho_d_mu[[d1]][[i]], drho_d_mu_2=drho_d_mu[[d2]][[i]], d2_mu=d2_mu[[d1]][[d2]], 
                                                                        drho_d2_mu=drho_d2_mu[[d1]][[d2]][[i]]))
          
          drho_d2_term[[d1]][[d2]][[i]] <-  crossprod(d2_X[[d1]][[d2]],drho_variance[[i]]*X) + 
            crossprod(d_X[[d1]],drho_d_variance[[d2]][[i]]*X + drho_variance[[i]]*d_X[[d2]]) +
            crossprod(d_X[[d2]],drho_d_variance[[d1]][[i]]*X + drho_variance[[i]]*d_X[[d1]]) +
            crossprod(X,drho_d2_variance[[d1]][[d2]][[i]]*X + drho_d_variance[[d1]][[i]]*d_X[[d2]] +
                        drho_d_variance[[d2]][[i]]*d_X[[d1]] + drho_variance[[i]]*d2_X[[d1]][[d2]]
            ) +
            drho_d2_S[[d1]][[d2]][[i]]
          
          drho_d2_inv_term[[d1]][[d2]][[i]] <- -1*( drho_d_inv_term[[d2]][[i]]  %*% t8[[d1]] + 
                                                      d_inv_term[[d2]]  %*% (drho_d_term[[d1]][[i]] %*% inv_term + d_term[[d1]] %*% drho_inv_term[[i]]) +
                                                      drho_inv_term[[i]] %*% (d2_term[[d1]][[d2]] %*% inv_term + d_term[[d1]] %*% d_inv_term[[d2]]) + 
                                                      inv_term %*% (drho_d2_term[[d1]][[d2]][[i]] %*% inv_term + d2_term[[d1]][[d2]] %*% drho_inv_term[[i]] +
                                                                      drho_d_term[[d1]][[i]] %*% d_inv_term[[d2]] + d_term[[d1]] %*% drho_d_inv_term[[d2]][[i]] )
          )
          drho_d2_dev[[d1]][[d2]][i] <- drho_d2_deviance_fn(family=family, y=y, mu=mu, d_mu_1=d_mu[[d1]], d_mu_2=d_mu[[d2]], drho_mu=drho_mu[[i]], 
                                                            drho_d_mu_1=drho_d_mu[[d1]][[i]], drho_d_mu_2=drho_d_mu[[d2]][[i]], d2_mu=d2_mu[[d1]][[d2]], 
                                                            drho_d2_mu=drho_d2_mu[[d1]][[d2]][[i]])
          if (scoreType %in% c("REML", "ML")){
            drho_d2_hat[[d1]][[d2]][i] <- NA
            
            drho_d2_rank.tol[[d1]][[d2]][i] <- sum(diag( drho_d_inv_term[[d2]][[i]] %*% d_term[[d1]] + 
                                                           d_inv_term[[d2]] %*% drho_d_term[[d1]][[i]] + 
                                                           drho_inv_term[[i]] %*% d2_term[[d1]][[d2]] + 
                                                           inv_term %*% drho_d2_term[[d1]][[d2]][[i]] )) 
            
            drho_d2_conv.tol[[d1]][[d2]][i] <- crossprod(drho_d2_coef[[d1]][[d2]][[i]],St) %*% coef + 
              crossprod(drho_d_coef[[d1]][[i]],d_St[[d2]]) %*% coef + 
              crossprod(drho_d_coef[[d1]][[i]],St) %*% d_coef[[d2]] + 
              crossprod(d2_coef[[d1]][[d2]],drho_S[[i]]) %*% coef + 
              crossprod(d_coef[[d1]],drho_d_S[[d2]][[i]]) %*% coef + 
              crossprod(d_coef[[d1]],drho_S[[i]]) %*% d_coef[[d2]] + 
              crossprod(d2_coef[[d1]][[d2]],St) %*% drho_coef[[i]] + 
              crossprod(d_coef[[d1]],d_St[[d2]]) %*% drho_coef[[i]] + 
              crossprod(d_coef[[d1]],St) %*% drho_d_coef[[d2]][[i]] + 
              crossprod(drho_d_coef[[d2]][[i]],d_St[[d1]]) %*% coef + 
              crossprod(drho_coef[[i]],d2_St[[d1]][[d2]]) %*% coef + 
              crossprod(drho_coef[[i]],d_St[[d1]]) %*% d_coef[[d2]] + 
              crossprod(d_coef[[d2]],drho_d_S[[d1]][[i]]) %*% coef + 
              crossprod(coef,drho_d2_S[[d1]][[d2]][[i]]) %*% coef + 
              crossprod(coef,drho_d_S[[d1]][[i]]) %*% d_coef[[d2]] + 
              crossprod(d_coef[[d2]],d_St[[d1]]) %*% drho_coef[[i]] + 
              crossprod(coef,d2_St[[d1]][[d2]]) %*% drho_coef[[i]] + 
              crossprod(coef,d_St[[d1]]) %*% drho_d_coef[[d2]][[i]] + 
              crossprod(drho_d_coef[[d2]][[i]],St) %*% d_coef[[d1]] + 
              crossprod(drho_coef[[i]],d_St[[d2]]) %*% d_coef[[d1]] + 
              crossprod(drho_coef[[i]],St) %*% d2_coef[[d1]][[d2]] + 
              crossprod(d_coef[[d2]],drho_S[[i]]) %*% d_coef[[d1]] + 
              crossprod(coef,drho_d_S[[d2]][[i]]) %*% d_coef[[d1]] + 
              crossprod(coef,drho_S[[i]]) %*% d2_coef[[d1]][[d2]] + 
              crossprod(d_coef[[d2]],St) %*% drho_d_coef[[d1]][[i]] + 
              crossprod(coef,d_St[[d2]]) %*% drho_d_coef[[d1]][[i]] + 
              crossprod(coef,St) %*% drho_d2_coef[[d1]][[d2]][[i]]
          }
          else
          {
            drho_d2_hat[[d1]][[d2]][i] <- sum(rowSums(
              (drho_d2_variance[[d1]][[d2]][[i]]*X %*% inv_term)*X + (drho_d_variance[[d1]][[i]]*d_X[[d2]] %*% inv_term)*X +     
                (drho_d_variance[[d1]][[i]]*X %*% d_inv_term[[d2]])*X + (drho_d_variance[[d1]][[i]]*X %*% inv_term)*d_X[[d2]] +
                (drho_d_variance[[d2]][[i]]*d_X[[d1]] %*% inv_term)*X  + (drho_d_variance[[d2]][[i]]*X %*% inv_term)*d_X[[d1]] + 
                (drho_d_variance[[d2]][[i]]*X %*% d_inv_term[[d1]])*X + (drho_variance[[i]]*d_X[[d2]] %*% inv_term)*d_X[[d1]] + 
                (drho_variance[[i]]*d_X[[d2]] %*% d_inv_term[[d1]])*X + (drho_variance[[i]]*X %*% d_inv_term[[d2]])*d_X[[d1]] +
                (drho_variance[[i]]*X %*% inv_term)*d2_X[[d1]][[d2]] + (drho_variance[[i]]*X %*% d2_inv_term[[d1]][[d2]])*X +
                (drho_variance[[i]]*X %*% d_inv_term[[d1]])*d_X[[d2]] + (drho_variance[[i]]*d2_X[[d1]][[d2]] %*% inv_term)*X + 
                (drho_variance[[i]]*d_X[[d1]] %*% d_inv_term[[d2]])*X + (drho_variance[[i]]*d_X[[d1]] %*% inv_term)*d_X[[d2]]  +   
                (d2_variance[[d1]][[d2]]*X %*% drho_inv_term[[i]])*X + (d_variance[[d1]]*d_X[[d2]] %*% drho_inv_term[[i]])*X + 
                (d_variance[[d1]]*X %*% drho_d_inv_term[[d2]][[i]])*X + (d_variance[[d1]]*X %*% drho_inv_term[[i]])*d_X[[d2]] +
                (variance*d2_X[[d1]][[d2]] %*% drho_inv_term[[i]])*X + (variance*d_X[[d1]] %*% drho_d_inv_term[[d2]][[i]])*X + 
                (variance*d_X[[d1]] %*% drho_inv_term[[i]])*d_X[[d2]] + (variance*d_X[[d2]] %*% drho_d_inv_term[[d1]][[i]])*X + 
                (variance*d_X[[d2]] %*% drho_inv_term[[i]])*d_X[[d1]] + (variance*X %*% drho_d2_inv_term[[d1]][[d2]][[i]])*X + 
                (variance*X %*% drho_d_inv_term[[d1]][[i]])*d_X[[d2]] + (variance*X %*% drho_d_inv_term[[d2]][[i]])*d_X[[d1]] + 
                (variance*X %*% drho_inv_term[[i]])*d2_X[[d1]][[d2]] + (d_variance[[d2]]*d_X[[d1]] %*% drho_inv_term[[i]])*X + 
                (d_variance[[d2]]*X %*% drho_d_inv_term[[d1]][[i]])*X + (d_variance[[d2]]*X %*% drho_inv_term[[i]])*d_X[[d1]]           
            ))
            drho_d2_rank.tol[[d1]][[d2]][i] <- NA
            drho_d2_conv.tol[[d1]][[d2]][i] <- NA
          }
        }
      }
      
      if(deriv == 2){
        d2rho_d2_eta[[d1]][[d2]] <- d2rho_d2_mu[[d1]][[d2]] <- d2rho_d2_variance[[d1]][[d2]] <- d2rho_d2_term[[d1]][[d2]] <- d2rho_d2_inv_term[[d1]][[d2]] <- list()
        d2rho_d2_dev[[d1]][[d2]] <- d2rho_d2_hat[[d1]][[d2]] <- d2rho_d2_rank.tol[[d1]][[d2]] <- d2rho_d2_conv.tol[[d1]][[d2]] <- matrix(nrow=length(sp),ncol=length(sp))
        for(i in 1:length(sp))
        {
          d2rho_d2_eta[[d1]][[d2]][[i]] <- d2rho_d2_mu[[d1]][[d2]][[i]] <- d2rho_d2_variance[[d1]][[d2]][[i]] <- d2rho_d2_term[[d1]][[d2]][[i]] <- d2rho_d2_inv_term[[d1]][[d2]][[i]] <- list()
          for(j in i:length(sp))
          {
            d2rho_d2_eta[[d1]][[d2]][[i]][[j]] <- ( tcrossprod(d2_X[[d1]][[d2]] , t(d2rho_coef[[i]][[j]])) + 
                                                      tcrossprod(d_X[[d1]] , t(d2rho_d_coef[[d2]][[i]][[j]])) + 
                                                      tcrossprod(d_X[[d2]] , t(d2rho_d_coef[[d1]][[i]][[j]])) + 
                                                      tcrossprod(X , t(d2rho_d2_coef[[d1]][[d2]][[i]][[j]])) )
            
            d2rho_d2_mu[[d1]][[d2]][[i]][[j]] <-  d2rho_d2_eta[[d1]][[d2]][[i]][[j]] * variance + d2rho_d_eta[[d1]][[i]][[j]] * d_variance[[d2]] + 
              drho_d2_eta[[d1]][[d2]][[i]] * drho_variance[[j]] + drho_d_eta[[d1]][[i]] * drho_d_variance[[d2]][[j]] +
              drho_d2_eta[[d1]][[d2]][[j]] * drho_variance[[i]] + drho_d_eta[[d1]][[j]] * drho_d_variance[[d2]][[i]] + 
              d2_eta[[d1]][[d2]] * d2rho_variance[[i]][[j]] + d_eta[[d1]] * d2rho_d_variance[[d2]][[i]][[j]]
            d2rho_d2_variance[[d1]][[d2]][[i]][[j]] <- as.vector(d2rho_d2_var_fn(family, mu=mu, d_mu_1=d_mu[[d1]], d_mu_2=d_mu[[d2]], drho_mu_i=drho_mu[[i]], drho_mu_j=drho_mu[[j]], 
                                                                                 d2_mu=d2_mu[[d1]][[d2]], d2rho_mu=d2rho_mu[[i]][[j]], drho_d_mu_i_1=drho_d_mu[[d1]][[i]], drho_d_mu_i_2=drho_d_mu[[d2]][[i]], 
                                                                                 drho_d_mu_j_1=drho_d_mu[[d1]][[j]], drho_d_mu_j_2=drho_d_mu[[d2]][[j]], drho_d2_mu_i=drho_d2_mu[[d1]][[d2]][[i]], drho_d2_mu_j=drho_d2_mu[[d1]][[d2]][[j]], 
                                                                                 d2rho_d_mu_1=d2rho_d_mu[[d1]][[i]][[j]], d2rho_d_mu_2=d2rho_d_mu[[d2]][[i]][[j]], d2rho_d2_mu=d2rho_d2_mu[[d1]][[d2]][[i]][[j]]))
            d2rho_d2_term[[d1]][[d2]][[i]][[j]] <- crossprod(d2_X[[d1]][[d2]],t16[[i]][[j]]) + 
              crossprod(d_X[[d1]],d2rho_d_variance[[d2]][[i]][[j]]*X + d2rho_variance[[i]][[j]]*d_X[[d2]]) +  
              crossprod(d_X[[d2]],d2rho_d_variance[[d1]][[i]][[j]]*X + d2rho_variance[[i]][[j]]*d_X[[d1]]) + 
              
              crossprod(X,d2rho_d2_variance[[d1]][[d2]][[i]][[j]]*X + d2rho_d_variance[[d1]][[i]][[j]]*d_X[[d2]] +
                          d2rho_d_variance[[d2]][[i]][[j]]*d_X[[d1]] + d2rho_variance[[i]][[j]]*d2_X[[d1]][[d2]]
              ) +  
              d2rho_d2_S[[d1]][[d2]][[i]][[j]]
            
            d2rho_d2_inv_term[[d1]][[d2]][[i]][[j]] <-  -1*( d2rho_d_inv_term[[d2]][[i]][[j]]  %*% t8[[d1]] + 
                                                               
                                                               
                                                               drho_d_inv_term[[d2]][[i]]  %*% (drho_d_term[[d1]][[j]] %*% inv_term + 
                                                                                                  d_term[[d1]] %*% drho_inv_term[[j]]) + 
                                                               drho_d_inv_term[[d2]][[j]]  %*% (drho_d_term[[d1]][[i]] %*% inv_term +
                                                                                                  d_term[[d1]] %*% drho_inv_term[[i]]) + 
                                                               d_inv_term[[d2]]  %*% (d2rho_d_term[[d1]][[i]][[j]] %*% inv_term +
                                                                                        drho_d_term[[d1]][[i]] %*% drho_inv_term[[j]] +
                                                                                        drho_d_term[[d1]][[j]] %*% drho_inv_term[[i]] +
                                                                                        d_term[[d1]] %*% d2rho_inv_term[[i]][[j]]
                                                               ) + 
                                                               d2rho_inv_term[[i]][[j]] %*% (d2_term[[d1]][[d2]] %*% inv_term +
                                                                                               d_term[[d1]] %*% d_inv_term[[d2]]
                                                               ) + 
                                                               drho_inv_term[[i]] %*% (drho_d2_term[[d1]][[d2]][[j]] %*% inv_term +
                                                                                         d2_term[[d1]][[d2]] %*% drho_inv_term[[j]] +
                                                                                         drho_d_term[[d1]][[j]] %*% d_inv_term[[d2]] +
                                                                                         d_term[[d1]] %*% drho_d_inv_term[[d2]][[j]]
                                                               ) + 
                                                               drho_inv_term[[j]] %*% (drho_d2_term[[d1]][[d2]][[i]] %*% inv_term +
                                                                                         d2_term[[d1]][[d2]] %*% drho_inv_term[[i]] +
                                                                                         drho_d_term[[d1]][[i]] %*% d_inv_term[[d2]] +
                                                                                         d_term[[d1]] %*% drho_d_inv_term[[d2]][[i]]
                                                               ) +
                                                               inv_term %*% (d2rho_d2_term[[d1]][[d2]][[i]][[j]] %*% inv_term +
                                                                               drho_d2_term[[d1]][[d2]][[i]] %*% drho_inv_term[[j]] +
                                                                               drho_d2_term[[d1]][[d2]][[j]] %*% drho_inv_term[[i]] +
                                                                               d2_term[[d1]][[d2]] %*% d2rho_inv_term[[i]][[j]] +
                                                                               d2rho_d_term[[d1]][[i]][[j]] %*% d_inv_term[[d2]] +
                                                                               drho_d_term[[d1]][[i]] %*% drho_d_inv_term[[d2]][[j]] +
                                                                               drho_d_term[[d1]][[j]] %*% drho_d_inv_term[[d2]][[i]] +
                                                                               d_term[[d1]] %*% d2rho_d_inv_term[[d2]][[i]][[j]]
                                                               )
            )
            d2rho_d2_dev[[d1]][[d2]][i,j] <- d2rho_d2_deviance_fn(family=family, y=y, mu=mu, d_mu_1 = d_mu[[d1]], d_mu_2 = d_mu[[d2]], 
                                                                  drho_mu_i=drho_mu[[i]], drho_mu_j=drho_mu[[j]], d2rho_mu=d2rho_mu[[i]][[j]], d2_mu=d2_mu[[d1]][[d2]], 
                                                                  drho_d_mu_i_1 = drho_d_mu[[d1]][[i]], drho_d_mu_i_2 = drho_d_mu[[d2]][[i]], drho_d_mu_j_1 = drho_d_mu[[d1]][[j]], drho_d_mu_j_2 = drho_d_mu[[d2]][[j]],
                                                                  drho_d2_mu_i = drho_d2_mu[[d1]][[d2]][[i]], drho_d2_mu_j = drho_d2_mu[[d1]][[d2]][[j]], 
                                                                  d2rho_d_mu_1 = d2rho_d_mu[[d1]][[i]][[j]], d2rho_d_mu_2 = d2rho_d_mu[[d2]][[i]][[j]], 
                                                                  d2rho_d2_mu = d2rho_d2_mu[[d1]][[d2]][[i]][[j]], d_eta_1 = d_eta[[d1]], d_eta_2 = d_eta[[d2]], 
                                                                  drho_eta_i=drho_eta[[i]], drho_eta_j=drho_eta[[j]], d2rho_eta=d2rho_eta[[i]][[j]], d2_eta=d2_eta[[d1]][[d2]], 
                                                                  drho_d_eta_i_1 = drho_d_eta[[d1]][[i]], drho_d_eta_i_2 = drho_d_eta[[d2]][[i]], drho_d_eta_j_1 = drho_d_eta[[d1]][[j]], drho_d_eta_j_2 = drho_d_eta[[d2]][[j]],
                                                                  drho_d2_eta_i = drho_d2_eta[[d1]][[d2]][[i]], drho_d2_eta_j = drho_d2_eta[[d1]][[d2]][[j]], 
                                                                  d2rho_d_eta_1 = d2rho_d_eta[[d1]][[i]][[j]], d2rho_d_eta_2 = d2rho_d_eta[[d2]][[i]][[j]], d2rho_d2_eta = d2rho_d2_eta[[d1]][[d2]][[i]][[j]])
            
            if (scoreType %in% c("REML", "ML")){
              d2rho_d2_hat[[d1]][[d2]][i,j] <- NA
              d2rho_d2_rank.tol[[d1]][[d2]][i,j] <- sum(diag( 
                d2rho_d_inv_term[[d1]][[i]][[j]] %*% d_term[[d2]] + drho_d_inv_term[[d1]][[i]] %*% drho_d_term[[d2]][[j]] +   
                  drho_d_inv_term[[d1]][[j]] %*% drho_d_term[[d2]][[i]] + d_inv_term[[d1]] %*% d2rho_d_term[[d2]][[i]][[j]] +  
                  d2rho_inv_term[[i]][[j]] %*% d2_term[[d1]][[d2]] + drho_inv_term[[i]] %*% drho_d2_term[[d1]][[d2]][[j]] +  
                  drho_inv_term[[j]] %*% drho_d2_term[[d1]][[d2]][[i]] + inv_term %*% d2rho_d2_term[[d1]][[d2]][[i]][[j]] 
              ))
              d2rho_d2_conv.tol[[d1]][[d2]][i,j] <- 
                crossprod(d2rho_d2_coef[[d1]][[d2]][[i]][[j]],St) %*% coef + crossprod(drho_d2_coef[[d1]][[d2]][[i]],drho_S[[j]]) %*% coef + 
                crossprod(drho_d2_coef[[d1]][[d2]][[i]],St) %*% drho_coef[[j]] + crossprod(d2rho_d_coef[[d1]][[i]][[j]],d_St[[d2]]) %*% coef + 
                crossprod(drho_d_coef[[d1]][[i]],drho_d_S[[d2]][[j]]) %*% coef + crossprod(drho_d_coef[[d1]][[i]],d_St[[d2]]) %*% drho_coef[[j]] + 
                crossprod(d2rho_d_coef[[d1]][[i]][[j]],St) %*% d_coef[[d2]] + crossprod(drho_d_coef[[d1]][[i]],drho_S[[j]]) %*% d_coef[[d2]] + 
                crossprod(drho_d_coef[[d1]][[i]],St) %*% drho_d_coef[[d2]][[j]] + crossprod(drho_d2_coef[[d1]][[d2]][[j]],drho_S[[i]]) %*% coef + 
                crossprod(d2_coef[[d1]][[d2]],d2rho_S[[i]][[j]]) %*% coef + crossprod(d2_coef[[d1]][[d2]],drho_S[[i]]) %*% drho_coef[[j]] + 
                crossprod(drho_d_coef[[d1]][[j]],drho_d_S[[d2]][[i]]) %*% coef + crossprod(d_coef[[d1]],d2rho_d_S[[d2]][[i]][[j]]) %*% coef + 
                crossprod(d_coef[[d1]],drho_d_S[[d2]][[i]]) %*% drho_coef[[j]] + crossprod(drho_d_coef[[d1]][[j]],drho_S[[i]]) %*% d_coef[[d2]] + 
                crossprod(d_coef[[d1]],d2rho_S[[i]][[j]]) %*% d_coef[[d2]] + crossprod(d_coef[[d1]],drho_S[[i]]) %*% drho_d_coef[[d2]][[j]] + 
                crossprod(drho_d2_coef[[d1]][[d2]][[j]],St) %*% drho_coef[[i]] + crossprod(d2_coef[[d1]][[d2]],drho_S[[j]]) %*% drho_coef[[i]] + 
                crossprod(d2_coef[[d1]][[d2]],St) %*% d2rho_coef[[i]][[j]] + crossprod(drho_d_coef[[d1]][[j]],d_St[[d2]]) %*% drho_coef[[i]] + 
                crossprod(d_coef[[d1]],drho_d_S[[d2]][[j]]) %*% drho_coef[[i]] + crossprod(d_coef[[d1]],d_St[[d2]]) %*% d2rho_coef[[i]][[j]] + 
                crossprod(drho_d_coef[[d1]][[j]],St) %*% drho_d_coef[[d2]][[i]] + crossprod(d_coef[[d1]],drho_S[[j]]) %*% drho_d_coef[[d2]][[i]] +
                crossprod(d_coef[[d1]],St) %*% d2rho_d_coef[[d2]][[i]][[j]] + crossprod(d2rho_d_coef[[d2]][[i]][[j]],d_St[[d1]]) %*% coef + 
                crossprod(drho_d_coef[[d2]][[i]],drho_d_S[[d1]][[j]]) %*% coef + crossprod(drho_d_coef[[d2]][[i]],d_St[[d1]]) %*% drho_coef[[j]] + 
                crossprod(d2rho_coef[[i]][[j]],d2_St[[d1]][[d2]]) %*% coef + crossprod(drho_coef[[i]],drho_d2_S[[d1]][[d2]][[j]]) %*% coef + 
                crossprod(drho_coef[[i]],d2_St[[d1]][[d2]]) %*% drho_coef[[j]] + crossprod(d2rho_coef[[i]][[j]],d_St[[d1]]) %*% d_coef[[d2]] + 
                crossprod(drho_coef[[i]],drho_d_S[[d1]][[j]]) %*% d_coef[[d2]] + crossprod(drho_coef[[i]],d_St[[d1]]) %*% drho_d_coef[[d2]][[j]] + 
                crossprod(drho_d_coef[[d2]][[j]],drho_d_S[[d1]][[i]]) %*% coef + crossprod(d_coef[[d2]],d2rho_d_S[[d1]][[i]][[j]]) %*% coef + 
                crossprod(d_coef[[d2]],drho_d_S[[d1]][[i]]) %*% drho_coef[[j]] + crossprod(drho_coef[[j]],drho_d2_S[[d1]][[d2]][[i]]) %*% coef + 
                crossprod(coef,d2rho_d2_S[[d1]][[d2]][[i]][[j]]) %*% coef + crossprod(coef,drho_d2_S[[d1]][[d2]][[i]]) %*% drho_coef[[j]] + 
                crossprod(drho_coef[[j]],drho_d_S[[d1]][[i]]) %*% d_coef[[d2]] + crossprod(coef,d2rho_d_S[[d1]][[i]][[j]]) %*% d_coef[[d2]] + 
                crossprod(coef,drho_d_S[[d1]][[i]]) %*% drho_d_coef[[d2]][[j]] + crossprod(drho_d_coef[[d2]][[j]],d_St[[d1]]) %*% drho_coef[[i]] + 
                crossprod(d_coef[[d2]],drho_d_S[[d1]][[j]]) %*% drho_coef[[i]] + crossprod(d_coef[[d2]],d_St[[d1]]) %*% d2rho_coef[[i]][[j]] + 
                crossprod(drho_coef[[j]],d2_St[[d1]][[d2]]) %*% drho_coef[[i]] + crossprod(coef,drho_d2_S[[d1]][[d2]][[j]]) %*% drho_coef[[i]] + 
                crossprod(coef,d2_St[[d1]][[d2]]) %*% d2rho_coef[[i]][[j]] + crossprod(drho_coef[[j]],d_St[[d1]]) %*% drho_d_coef[[d2]][[i]] + 
                crossprod(coef,drho_d_S[[d1]][[j]]) %*% drho_d_coef[[d2]][[i]] + crossprod(coef,d_St[[d1]]) %*% d2rho_d_coef[[d2]][[i]][[j]] + 
                crossprod(d2rho_d_coef[[d2]][[i]][[j]],St) %*% d_coef[[d1]] + crossprod(drho_d_coef[[d2]][[i]],drho_S[[j]]) %*% d_coef[[d1]] + 
                crossprod(drho_d_coef[[d2]][[i]],St) %*% drho_d_coef[[d1]][[j]] + crossprod(d2rho_coef[[i]][[j]],d_St[[d2]]) %*% d_coef[[d1]] + 
                crossprod(drho_coef[[i]],drho_d_S[[d2]][[j]]) %*% d_coef[[d1]] + crossprod(drho_coef[[i]],d_St[[d2]]) %*% drho_d_coef[[d1]][[j]] + 
                crossprod(d2rho_coef[[i]][[j]],St) %*% d2_coef[[d1]][[d2]] + crossprod(drho_coef[[i]],drho_S[[j]]) %*% d2_coef[[d1]][[d2]] + 
                crossprod(drho_coef[[i]],St) %*% drho_d2_coef[[d1]][[d2]][[j]] + crossprod(drho_d_coef[[d2]][[j]],drho_S[[i]]) %*% d_coef[[d1]] + 
                crossprod(d_coef[[d2]],d2rho_S[[i]][[j]]) %*% d_coef[[d1]] + crossprod(d_coef[[d2]],drho_S[[i]]) %*% drho_d_coef[[d1]][[j]] + 
                crossprod(drho_coef[[j]],drho_d_S[[d2]][[i]]) %*% d_coef[[d1]] + crossprod(coef,d2rho_d_S[[d2]][[i]][[j]]) %*% d_coef[[d1]] + 
                crossprod(coef,drho_d_S[[d2]][[i]]) %*% drho_d_coef[[d1]][[j]] + crossprod(drho_coef[[j]],drho_S[[i]]) %*% d2_coef[[d1]][[d2]] + 
                crossprod(coef,d2rho_S[[i]][[j]]) %*% d2_coef[[d1]][[d2]] + crossprod(coef,drho_S[[i]]) %*% drho_d2_coef[[d1]][[d2]][[j]] + 
                crossprod(drho_d_coef[[d2]][[j]],St) %*% drho_d_coef[[d1]][[i]] + crossprod(d_coef[[d2]],drho_S[[j]]) %*% drho_d_coef[[d1]][[i]] + 
                crossprod(d_coef[[d2]],St) %*% d2rho_d_coef[[d1]][[i]][[j]] + crossprod(drho_coef[[j]],d_St[[d2]]) %*% drho_d_coef[[d1]][[i]] + 
                crossprod(coef,drho_d_S[[d2]][[j]]) %*% drho_d_coef[[d1]][[i]] + crossprod(coef,d_St[[d2]]) %*% d2rho_d_coef[[d1]][[i]][[j]] + 
                crossprod(drho_coef[[j]],St) %*% drho_d2_coef[[d1]][[d2]][[i]] + crossprod(coef,drho_S[[j]]) %*% drho_d2_coef[[d1]][[d2]][[i]] + 
                crossprod(coef,St) %*% d2rho_d2_coef[[d1]][[d2]][[i]][[j]]
              
            }
            else 
            {
              d2rho_d2_hat[[d1]][[d2]][i,j] <- sum(rowSums(
                (d2rho_d2_variance[[d1]][[d2]][[i]][[j]]*X %*% inv_term)*X + (drho_d2_variance[[d1]][[d2]][[i]]*X %*% drho_inv_term[[j]])*X + 
                  (d2rho_d_variance[[d1]][[i]][[j]]*d_X[[d2]] %*% inv_term)*X + (d2rho_d_variance[[d1]][[i]][[j]]*X %*% d_inv_term[[d2]])*X +  
                  (d2rho_d_variance[[d1]][[i]][[j]]*X %*% inv_term)*d_X[[d2]] + (drho_d_variance[[d1]][[i]]*d_X[[d2]] %*% drho_inv_term[[j]])*X + 
                  (drho_d_variance[[d1]][[i]]*X %*% drho_d_inv_term[[d2]][[j]])*X + (drho_d_variance[[d1]][[i]]*X %*% drho_inv_term[[j]])*d_X[[d2]] + 
                  (d2rho_d_variance[[d2]][[i]][[j]]* d_X[[d1]] %*% inv_term)*X + (d2rho_d_variance[[d2]][[i]][[j]]*X %*% inv_term)*d_X[[d1]] + 
                  (d2rho_d_variance[[d2]][[i]][[j]]*X %*% d_inv_term[[d1]])*X + (drho_d_variance[[d2]][[i]]*d_X[[d1]] %*% drho_inv_term[[j]])*X + 
                  (drho_d_variance[[d2]][[i]]*X %*% drho_inv_term[[j]])*d_X[[d1]] + (drho_d_variance[[d2]][[i]]*X %*% drho_d_inv_term[[d1]][[j]])*X +
                  (d2rho_variance[[i]][[j]]*d_X[[d2]] %*% inv_term)*d_X[[d1]] + (d2rho_variance[[i]][[j]]*d_X[[d2]] %*% d_inv_term[[d1]])*X +
                  (d2rho_variance[[i]][[j]]*X %*% d_inv_term[[d2]])*d_X[[d1]] + (d2rho_variance[[i]][[j]]*X %*% inv_term)*d2_X[[d1]][[d2]] + 
                  (d2rho_variance[[i]][[j]]*X %*% d2_inv_term[[d1]][[d2]])*X + (d2rho_variance[[i]][[j]]*X %*% d_inv_term[[d1]])*d_X[[d2]] +
                  (d2rho_variance[[i]][[j]]*d2_X[[d1]][[d2]] %*% inv_term)*X + (d2rho_variance[[i]][[j]]*d_X[[d1]] %*% d_inv_term[[d2]])*X + 
                  (d2rho_variance[[i]][[j]]*d_X[[d1]] %*% inv_term)*d_X[[d2]] + (drho_variance[[i]]*d_X[[d2]] %*% drho_inv_term[[j]])*d_X[[d1]] + 
                  (drho_variance[[i]]*d_X[[d2]] %*% drho_d_inv_term[[d1]][[j]])*X + (drho_variance[[i]]*X %*% drho_d_inv_term[[d2]][[j]])*d_X[[d1]] + 
                  (drho_variance[[i]]*X %*% drho_inv_term[[j]])*d2_X[[d1]][[d2]] + (drho_variance[[i]]*X %*% drho_d2_inv_term[[d1]][[d2]][[j]])*X + 
                  (drho_variance[[i]]*X %*% drho_d_inv_term[[d1]][[j]])*d_X[[d2]] + (drho_variance[[i]]*d2_X[[d1]][[d2]] %*% drho_inv_term[[j]])*X + 
                  (drho_variance[[i]]*d_X[[d1]] %*% drho_d_inv_term[[d2]][[j]])*X + (drho_variance[[i]]*d_X[[d1]] %*% drho_inv_term[[j]])*d_X[[d2]] +
                  (drho_d2_variance[[d1]][[d2]][[j]]*X %*% drho_inv_term[[i]])*X + (d2_variance[[d1]][[d2]]*X %*% d2rho_inv_term[[i]][[j]])*X +
                  (drho_d_variance[[d1]][[j]]* d_X[[d2]] %*% drho_inv_term[[i]])*X + (drho_d_variance[[d1]][[j]]*X %*% drho_d_inv_term[[d2]][[i]])*X +  
                  (drho_d_variance[[d1]][[j]]*X %*% drho_inv_term[[i]])*d_X[[d2]] + (d_variance[[d1]]*d_X[[d2]] %*% d2rho_inv_term[[i]][[j]])*X  + 
                  (d_variance[[d1]]*X %*% d2rho_d_inv_term[[d2]][[i]][[j]])*X + (d_variance[[d1]]*X %*% d2rho_inv_term[[i]][[j]])*d_X[[d2]] + 
                  (drho_variance[[j]]*d2_X[[d1]][[d2]] %*% drho_inv_term[[i]])*X + (drho_variance[[j]]*d_X[[d1]] %*% drho_d_inv_term[[d2]][[i]])*X + 
                  (drho_variance[[j]]*d_X[[d1]] %*% drho_inv_term[[i]])*d_X[[d2]] + (drho_variance[[j]]*d_X[[d2]] %*% drho_d_inv_term[[d1]][[i]])*X + 
                  (drho_variance[[j]]*d_X[[d2]] %*% drho_inv_term[[i]])*d_X[[d1]] + (drho_variance[[j]]*X %*% drho_d2_inv_term[[d1]][[d2]][[i]])*X  + 
                  (drho_variance[[j]]*X %*% drho_d_inv_term[[d1]][[i]])*d_X[[d2]]  + (drho_variance[[j]]*X %*% drho_d_inv_term[[d2]][[i]])*d_X[[d1]]  + 
                  (drho_variance[[j]]*X %*% drho_inv_term[[i]])*d2_X[[d1]][[d2]]  + (variance*d2_X[[d1]][[d2]] %*% d2rho_inv_term[[i]][[j]])*X + 
                  (variance*d_X[[d1]] %*% d2rho_d_inv_term[[d2]][[i]][[j]])*X + (variance*d_X[[d1]] %*% d2rho_inv_term[[i]][[j]])*d_X[[d2]] +
                  (variance*d_X[[d2]] %*% d2rho_d_inv_term[[d1]][[i]][[j]])*X + (variance*d_X[[d2]] %*% d2rho_inv_term[[i]][[j]])*d_X[[d1]] +
                  (variance*X %*% d2rho_d2_inv_term[[d1]][[d2]][[i]][[j]])*X  + (variance*X %*% d2rho_d_inv_term[[d1]][[i]][[j]])*d_X[[d2]]  + 
                  (variance*X %*% d2rho_d_inv_term[[d2]][[i]][[j]])*d_X[[d1]]  + (variance*X %*% d2rho_inv_term[[i]][[j]])*d2_X[[d1]][[d2]]  + 
                  (drho_d_variance[[d2]][[j]]*d_X[[d1]] %*% drho_inv_term[[i]])*X + (drho_d_variance[[d2]][[j]]*X %*% drho_d_inv_term[[d1]][[i]])*X + 
                  (drho_d_variance[[d2]][[j]]*X %*% drho_inv_term[[i]])*d_X[[d1]] + (d_variance[[d2]]*d_X[[d1]] %*% d2rho_inv_term[[i]][[j]])*X +  
                  (d_variance[[d2]]*X %*% d2rho_d_inv_term[[d1]][[i]][[j]])*X + (d_variance[[d2]]*X %*% d2rho_inv_term[[i]][[j]])*d_X[[d1]]    
              ))  
              d2rho_d2_rank.tol[[d1]][[d2]][i,j] <- NA
              d2rho_d2_conv.tol[[d1]][[d2]][i,j] <- NA
            } 
          }
        } 
      }
    }
  }
  
  list(
    d_dev = d_dev, d2_dev = d2_dev, drho_d_dev = drho_d_dev, d2rho_d_dev = d2rho_d_dev, drho_d2_dev = drho_d2_dev, d2rho_d2_dev = d2rho_d2_dev,
    d_hat = d_hat, d2_hat = d2_hat, drho_d_hat = drho_d_hat, d2rho_d_hat = d2rho_d_hat, drho_d2_hat = drho_d2_hat, d2rho_d2_hat = d2rho_d2_hat, 
    d_rank.tol=d_rank.tol, d2_rank.tol=d2_rank.tol, drho_d_rank.tol = drho_d_rank.tol, d2rho_d_rank.tol = d2rho_d_rank.tol, drho_d2_rank.tol = drho_d2_rank.tol, 
    d2rho_d2_rank.tol = d2rho_d2_rank.tol, drho_d_conv.tol = drho_d_conv.tol, drho_d2_conv.tol=drho_d2_conv.tol, d2rho_d_conv.tol=d2rho_d_conv.tol,
    d_conv.tol = d_conv.tol, d2_conv.tol=d2_conv.tol, d2rho_d2_conv.tol=d2rho_d2_conv.tol
  ) 
}

d_coef_trans0 <- function(beta_star,T_mat,d_T,d2_T,coef,drho_coef, d2rho_coef, d_coef, 
                          d2_coef, drho_d_coef, drho_d2_coef, d2rho_d_coef, d2rho_d2_coef, sp){
  
  for(d1 in 1:length(beta_star)){
    for(d2 in d1:length(beta_star)){
      for(i in 1:length(sp)){
        for(j in i:length(sp)){
          # d2rho_d2_coef[[d1]][[d2]][[i]][[j]] <- drop(t(d2_T[[d1]][[d2]]) %*% d2rho_coef[[i]][[j]] + t(d_T[[d1]]) %*% d2rho_d_coef[[d2]][[i]][[j]] +    
          #                                               t(d_T[[d2]]) %*% d2rho_d_coef[[d1]][[i]][[j]] + t(T_mat) %*% d2rho_d2_coef[[d1]][[d2]][[i]][[j]])
          
          d2rho_d2_coef[[d1]][[d2]][[i]][[j]] <- drop(crossprod(d2_T[[d1]][[d2]],d2rho_coef[[i]][[j]]) + crossprod(d_T[[d1]],d2rho_d_coef[[d2]][[i]][[j]]) +
                                                        crossprod(d_T[[d2]],d2rho_d_coef[[d1]][[i]][[j]]) + crossprod(T_mat,d2rho_d2_coef[[d1]][[d2]][[i]][[j]]))
          
        }
        # drho_d2_coef[[d1]][[d2]][[i]] <- drop(t(d2_T[[d1]][[d2]]) %*% drho_coef[[i]] + t(d_T[[d1]]) %*% drho_d_coef[[d2]][[i]] +    
        #                                         t(d_T[[d2]]) %*% drho_d_coef[[d1]][[i]] + t(T_mat) %*% drho_d2_coef[[d1]][[d2]][[i]])
        
        drho_d2_coef[[d1]][[d2]][[i]] <- drop(crossprod(d2_T[[d1]][[d2]], drho_coef[[i]]) + crossprod(d_T[[d1]], drho_d_coef[[d2]][[i]]) +    
                                                crossprod(d_T[[d2]], drho_d_coef[[d1]][[i]]) + crossprod(T_mat, drho_d2_coef[[d1]][[d2]][[i]]))
        
      }
      # d2_coef[[d1]][[d2]] <- drop(t(d2_T[[d1]][[d2]]) %*% coef + t(d_T[[d1]]) %*% d_coef[[d2]] +    
      #                               t(d_T[[d2]]) %*% d_coef[[d1]] + t(T_mat) %*% d2_coef[[d1]][[d2]])
      
      d2_coef[[d1]][[d2]] <- drop(crossprod(d2_T[[d1]][[d2]], coef) + crossprod(d_T[[d1]], d_coef[[d2]]) +    
                                    crossprod(d_T[[d2]], d_coef[[d1]]) + crossprod(T_mat,d2_coef[[d1]][[d2]]))
      
    }
  }
  
  for(i in 1:length(sp))
  {
    for(j in i:length(sp))
    {
      for(d1 in 1:length(beta_star)){
        # d2rho_d_coef[[d1]][[i]][[j]] <- drop( (t(d_T[[d1]]) %*% d2rho_coef[[i]][[j]]) + (t(T_mat) %*% d2rho_d_coef[[d1]][[i]][[j]]) )
        
        d2rho_d_coef[[d1]][[i]][[j]] <- drop( crossprod(d_T[[d1]],d2rho_coef[[i]][[j]]) + crossprod(T_mat, d2rho_d_coef[[d1]][[i]][[j]]) )
      }
      # d2rho_coef[[i]][[j]] <- drop(t(T_mat) %*% d2rho_coef[[i]][[j]])
      d2rho_coef[[i]][[j]] <- drop(crossprod(T_mat, d2rho_coef[[i]][[j]]))
    }
    # drho_coef[[i]] <- drop(t(T_mat) %*% drho_coef[[i]])
    drho_coef[[i]] <- drop(crossprod(T_mat,drho_coef[[i]]))
  }
  
  for(d1 in 1:length(beta_star)){
    # d_coef[[d1]] <- drop(t(d_T[[d1]]) %*% coef + t(T_mat) %*% d_coef[[d1]])
    d_coef[[d1]] <- drop(crossprod(d_T[[d1]],coef) + crossprod(T_mat,d_coef[[d1]]))
  }
  
  list(coef=coef,drho_coef=drho_coef, d2rho_coef=d2rho_coef, d_coef=d_coef, d2_coef=d2_coef, 
       drho_d_coef=drho_d_coef, drho_d2_coef=drho_d2_coef, d2rho_d_coef=d2rho_d_coef, d2rho_d2_coef=d2rho_d2_coef)
  
}

d_coef_trans <- function(beta_star,T_mat,d_T,d2_T,coef,drho_coef, d2rho_coef, d_coef, d2_coef, 
                         drho_d_coef, drho_d2_coef, d2rho_d_coef, d2rho_d2_coef, sp){
  #We now calculate derivatives in non-reparametrized space
  for(d1 in 1:length(beta_star)){
    for(d2 in d1:length(beta_star)){
      for(i in 1:length(sp)){
        for(j in i:length(sp)){
          d2rho_d2_coef[[d1]][[d2]][[i]][[j]] <- drop(d2_T[[d1]][[d2]] %*% d2rho_coef[[i]][[j]] + d_T[[d1]] %*% d2rho_d_coef[[d2]][[i]][[j]] +    
                                                        d_T[[d2]] %*% d2rho_d_coef[[d1]][[i]][[j]] + T_mat %*% d2rho_d2_coef[[d1]][[d2]][[i]][[j]])
        }
        drho_d2_coef[[d1]][[d2]][[i]] <- drop( 
          (d2_T[[d1]][[d2]] %*% drho_coef[[i]]) + (d_T[[d1]] %*% drho_d_coef[[d2]][[i]]) +    
            (d_T[[d2]] %*% drho_d_coef[[d1]][[i]]) + T_mat %*% drho_d2_coef[[d1]][[d2]][[i]]
        )
        
      }
      d2_coef[[d1]][[d2]] <- drop( 
        (d2_T[[d1]][[d2]] %*% coef) + (d_T[[d1]] %*% d_coef[[d2]]) +    
          (d_T[[d2]] %*% d_coef[[d1]]) + T_mat %*% d2_coef[[d1]][[d2]]
      )
      
    }
  }
  
  for(i in 1:length(sp))
  {
    for(j in i:length(sp))
    {
      for(d1 in 1:length(beta_star)){
        d2rho_d_coef[[d1]][[i]][[j]] <- drop( (d_T[[d1]] %*% d2rho_coef[[i]][[j]]) + (T_mat %*% d2rho_d_coef[[d1]][[i]][[j]]) )
      }
      d2rho_coef[[i]][[j]] <- drop(T_mat %*% d2rho_coef[[i]][[j]])
    }
    drho_coef[[i]] <- drop(T_mat %*% drho_coef[[i]])
  }
  
  for(d1 in 1:length(beta_star)){
    d_coef[[d1]] <- drop( (d_T[[d1]] %*% coef) + (T_mat %*% d_coef[[d1]]) )
  }
  
  list(coef=coef,drho_coef=drho_coef, d2rho_coef=d2rho_coef, d_coef=d_coef, d2_coef=d2_coef, 
       drho_d_coef=drho_d_coef, drho_d2_coef=drho_d2_coef, d2rho_d_coef=d2rho_d_coef, d2rho_d2_coef=d2rho_d2_coef)
}

#08i
dQs <- function(beta_star,Si,penaltyRank,sp,d_Si, d_sp, d2_Si, d2_sp, dTol, rTol, drho_S, d2rho_S, drho_d_S, drho_d2_S, d2rho_d_S, d2rho_d2_S){
  Qs <- diag(penaltyRank)
  
  p <- penaltyRank
  q <- penaltyRank
  
  d_Qs <- list()
  d2_Qs <- list()
  for(deriv_1 in 1:length(beta_star))
  {
    d2_Qs[[deriv_1]] <- list()
    for(deriv_2 in deriv_1:length(beta_star))
    {
      d2_Qs[[deriv_1]][[deriv_2]] <- matrix(0,nrow=penaltyRank,ncol=penaltyRank)
    }
    d_Qs[[deriv_1]] <- matrix(0,nrow=penaltyRank,ncol=penaltyRank)
  }
  
  ##We calculate initial St, d_St, d2_St
  St <- matrix(0,nrow=p,ncol=q)
  for(k in 1:length(Si))
  {
    St <- St + exp(sp[k])*Si[[k]] 
    
  }
  
  d_St <- list()
  d2_St <- list()
  for(deriv_1 in 1:length(beta_star))
  {
    d_St[[deriv_1]] <- matrix(0,nrow=p,ncol=q)
    d2_St[[deriv_1]] <- list()
    for(deriv_2 in deriv_1:length(beta_star))
    {
      d2_St[[deriv_1]][[deriv_2]] <- matrix(0,nrow=p,ncol=q)
      for(k in 1:length(Si))
      {
        d2_St[[deriv_1]][[deriv_2]] <- d2_St[[deriv_1]][[deriv_2]] + 
          d_sp[[k]][deriv_2]*d_Si[[k]][[deriv_1]][1:p,1:q] + 
          exp(sp[k])*d2_Si[[k]][[deriv_1]][[deriv_2]][1:p,1:q] +
          d2_sp[[k]][deriv_1,deriv_2]*Si[[k]] + 
          d_sp[[k]][deriv_1]*d_Si[[k]][[deriv_2]][1:p,1:q]
      }
    }
    for(k in 1:length(Si))
    {
      d_St[[deriv_1]] <- d_St[[deriv_1]] + exp(sp[k])*d_Si[[k]][[deriv_1]][1:p,1:q] + d_sp[[k]][deriv_1]*Si[[k]]  
    }
  }
  
  ##We create global versions of St, d_St and d2_St
  m <- length(Si)
  K <- 0
  Q <- dim(St)[1]
  gamma <- rep(1,times=m)
  while(K != penaltyRank )
  {
    
    Omega <- rep(0,times = length(gamma))  
    # Step 1
    for(k in 1:length(gamma))
    {
      if(gamma[k] == 1)
      {
        Omega[k] <- sqrt(sum(Si[[k]] * Si[[k]])) * exp(sp[k])
      }
    }
    
    #Step 2
    alpha <- as.numeric( Omega >= dTol * max(Omega) ) 
    gammaPrime <- as.numeric( Omega < dTol * max(Omega) )
    
    SScaled <- matrix(0,nrow=nrow(St),ncol=ncol(St))
    #Step 3
    for(k in 1:length(alpha))
    {
      if(alpha[k] == 1)
      {
        SScaled <- SScaled + Si[[k]]/sqrt(sum(Si[[k]] * Si[[k]]))
      }
    }
    
    es <- eigen(SScaled, symmetric = TRUE)
    ind <- es$values > max(es$values) * rTol
    penaltyRank <- sum(ind)
    
    # Step 5
    SUnscaled <- matrix(0,nrow=nrow(St),ncol=ncol(St))
    for(k in 1:length(alpha))
    {
      if(alpha[k] == 1)
      {
        SUnscaled <- SUnscaled + exp(sp[k])*Si[[k]] 
      }
    }
    
    d_SUnscaled <- list()
    d2_SUnscaled <- list()
    for(deriv_1 in 1:length(beta_star))
    {
      d_SUnscaled[[deriv_1]] <- matrix(0,nrow=nrow(St),ncol=ncol(St))
      d2_SUnscaled[[deriv_1]] <- list()
      for(deriv_2 in deriv_1:length(beta_star))
      {
        d2_SUnscaled[[deriv_1]][[deriv_2]] <- matrix(0,nrow=nrow(St),ncol=ncol(St))
        for(k in 1:length(alpha))
        {
          if(alpha[k] == 1)
          {
            d2_SUnscaled[[deriv_1]][[deriv_2]] <- d2_SUnscaled[[deriv_1]][[deriv_2]] + 
              d_sp[[k]][deriv_2]*d_Si[[k]][[deriv_1]][1:p,1:q] + 
              exp(sp[k])*d2_Si[[k]][[deriv_1]][[deriv_2]][1:p,1:q] +
              d2_sp[[k]][deriv_1,deriv_2]*Si[[k]] + d_sp[[k]][deriv_1]*d_Si[[k]][[deriv_2]][1:p,1:q]
          }
        }
      }
      
      for(k in 1:length(alpha))
      {
        if(alpha[k] == 1)
        {
          d_SUnscaled[[deriv_1]] <- d_SUnscaled[[deriv_1]] + exp(sp[k])*d_Si[[k]][[deriv_1]][1:p,1:q] + d_sp[[k]][deriv_1]*Si[[k]]
        }
      }
    }
    
    es0 <- eigen(SUnscaled, symmetric = TRUE)
    
    U0 <- es0$vectors
    Ur <- es0$vectors[, 1:penaltyRank]
    Un <- es0$vectors[, -c(1:penaltyRank)]
    evals <- es0$values
    
    U <- matrix(0,nrow = (K+Q), ncol = (K+Q))
    Id <- diag(nrow=nrow(SUnscaled),ncol=ncol(SUnscaled))
    d_U0 <- list()
    d_evals <- list()
    d_Un <- list()
    for(deriv_1 in 1:length(beta_star))
    {
      d_U0[[deriv_1]] <- matrix(nrow=ncol(U0),ncol=ncol(U0))
      d_evals[[deriv_1]] <- vector()
      d_Un[[deriv_1]] <- matrix(nrow=ncol(Un),ncol=ncol(Un))
      for(k in 1:ncol(U0))
      {
        d_U0[[deriv_1]][,k] <- -1* ginv( SUnscaled - ( evals[k]*Id ) ) %*% (d_SUnscaled[[deriv_1]] %*% U0[,k])
        d_evals[[deriv_1]][k] <- (crossprod(U0,d_SUnscaled[[deriv_1]]))[k,] %*% U0[,k]
      }
      d_Un[[deriv_1]] <- d_U0[[deriv_1]][, -c(1:penaltyRank)]
    }
    
    d2_U0 <- list()
    d2_Un <- list()
    for(deriv_1 in 1:length(beta_star))
    {
      d2_U0[[deriv_1]] <- list()
      d2_Un[[deriv_1]] <- list()
      for(deriv_2 in deriv_1:length(beta_star))
      {
        d2_U0[[deriv_1]][[deriv_2]]<- matrix(nrow=ncol(U0),ncol=ncol(U0))
        for(k in 1:ncol(U0))
        {
          d2_U0[[deriv_1]][[deriv_2]][,k] <- (          
            -1*(-ginv(SUnscaled- (evals[k]*Id)) %*% ( d_SUnscaled[[deriv_2]]-(d_evals[[deriv_2]][k]*Id) ) %*% ginv(SUnscaled- (evals[k]*Id))              
                + ginv(SUnscaled- (evals[k]*Id)) %*% t(ginv(SUnscaled- (evals[k]*Id))) %*% t( d_SUnscaled[[deriv_2]]-(d_evals[[deriv_2]][k]*Id) ) %*% ( Id-(SUnscaled- (evals[k]*Id)) %*% ginv(SUnscaled- (evals[k]*Id)) )              
                + ( Id-(SUnscaled- (evals[k]*Id)) %*% ginv(SUnscaled- (evals[k]*Id)) ) %*% t( d_SUnscaled[[deriv_2]]-(d_evals[[deriv_2]][k]*Id) ) %*% t(ginv(SUnscaled- (evals[k]*Id))) %*% ginv(SUnscaled- (evals[k]*Id))             
            ) %*% d_SUnscaled[[deriv_1]] %*% U0[,k]          
            - ginv(SUnscaled- (evals[k]*Id)) %*% d2_SUnscaled[[deriv_1]][[deriv_2]] %*% U0[,k]          
            - ginv(SUnscaled- (evals[k]*Id)) %*% d_SUnscaled[[deriv_1]] %*% d_U0[[deriv_2]][,k]          
          )
        }
        d2_Un[[deriv_1]][[deriv_2]] <- d2_U0[[deriv_1]][[deriv_2]][, -c(1:penaltyRank)]
      }
    }
    d_U <- list()
    d2_U <- list()
    if(K==0)
    {
      U <- U0
      d2_U[[deriv_1]] <- list()
      for(deriv_1 in 1:length(beta_star))
      {
        for(deriv_2 in deriv_1:length(beta_star))
        {
          d2_U[[deriv_1]][[deriv_2]] <- d2_U0[[deriv_1]][[deriv_2]]
        }
        d_U[[deriv_1]] <- d_U0[[deriv_1]]
      }
      
    } else if (K > 0) {
      U[1:K,1:K] <- diag(K)
      U[(K+1):(K+Q),(K+1):(K+Q)] <- U0
      
      for(deriv_1 in 1:length(beta_star))
      {
        for(deriv_2 in deriv_1:length(beta_star))
        {
          d2_U[[deriv_1]][[deriv_2]][1:K,1:K] <- matrix(0,nrow = K,ncol = K)
          d2_U[[deriv_1]][[deriv_2]][(K+1):(K+Q),(K+1):(K+Q)] <- d2_U0[[deriv_1]][[deriv_2]]
        }
        d_U[[deriv_1]][1:K,1:K] <- matrix(0,nrow = K,ncol = K)
        d_U[[deriv_1]][(K+1):(K+Q),(K+1):(K+Q)] <- d_U0[[deriv_1]] 
      }
    }
    
    StPrime <-  crossprod(U, St) %*% U
    
    drho_SPrime <- list()
    d2rho_SPrime <- list()
    for(i in 1:length(sp))
    {
      d2rho_SPrime[[i]] <- list()
      for(j in i:length(sp))
      {
        d2rho_SPrime[[i]][[j]] <- crossprod(U, d2rho_S[[i]][[j]][1:p,1:q]) %*% U
      }
      
      drho_SPrime[[i]] <- crossprod(U, drho_S[[i]][1:p,1:q]) %*% U
    }
    
    dStPrime <- list()
    d2StPrime <- list()
    drho_d_SPrime <- list()
    d2rho_d_SPrime <- list()
    drho_d2_SPrime <- list()
    d2rho_d2_SPrime <- list()
    for(deriv_1 in 1:length(beta_star))
    {
      d2StPrime[[deriv_1]] <- list()
      drho_d_SPrime[[deriv_1]] <- list()
      d2rho_d_SPrime[[deriv_1]] <- list()
      drho_d2_SPrime[[deriv_1]] <- list()
      d2rho_d2_SPrime[[deriv_1]] <- list()
      
      for(deriv_2 in deriv_1:length(beta_star))
      {
        drho_d2_SPrime[[deriv_1]][[deriv_2]] <- list()
        d2rho_d2_SPrime[[deriv_1]][[deriv_2]] <- list()
        
        d2StPrime[[deriv_1]][[deriv_2]] <- crossprod(d2_U[[deriv_1]][[deriv_2]], St) %*% U + 
          crossprod(d_U[[deriv_1]], d_St[[deriv_2]]) %*% U + 
          crossprod(d_U[[deriv_1]], St) %*% d_U[[deriv_2]] + 
          crossprod(d_U[[deriv_2]], d_St[[deriv_1]]) %*% U + 
          crossprod(U, d2_St[[deriv_1]][[deriv_2]]) %*% U + 
          crossprod(U, d_St[[deriv_1]]) %*% d_U[[deriv_2]] + 
          crossprod(d_U[[deriv_2]], St) %*% d_U[[deriv_1]] +
          crossprod(U, d_St[[deriv_2]]) %*% d_U[[deriv_1]] +
          crossprod(U, St) %*% d2_U[[deriv_1]][[deriv_2]]
        
        for(i in 1:length(sp))
        {
          d2rho_d2_SPrime[[deriv_1]][[deriv_2]][[i]] <- list()
          for(j in i:length(sp))
          {
            d2rho_d2_SPrime[[deriv_1]][[deriv_2]][[i]][[j]] <- crossprod(d2_U[[deriv_1]][[deriv_2]], d2rho_S[[i]][[j]][1:p,1:q]) %*% U + 
              crossprod(d_U[[deriv_1]], drho_d2_S[[deriv_1]][[deriv_2]][[j]][1:p,1:q]) %*% U + 
              crossprod(d_U[[deriv_1]], d2rho_S[[i]][[j]][1:p,1:q]) %*% d_U[[deriv_2]] + 
              crossprod(d_U[[deriv_2]], d2rho_d_S[[deriv_1]][[i]][[j]][1:p,1:q]) %*% U + 
              crossprod(U, d2rho_d2_S[[deriv_1]][[deriv_2]][[i]][[j]][1:p,1:q]) %*% U + 
              crossprod(U, d2rho_d_S[[deriv_1]][[i]][[j]][1:p,1:q]) %*% d_U[[deriv_2]] + 
              crossprod(d_U[[deriv_2]], d2rho_S[[i]][[j]][1:p,1:q]) %*% d_U[[deriv_1]] +
              crossprod(U, d2rho_d_S[[deriv_2]][[i]][[j]][1:p,1:q]) %*% d_U[[deriv_1]] +
              crossprod(U, d2rho_S[[i]][[j]][1:p,1:q]) %*% d2_U[[deriv_1]][[deriv_2]]
          }
          drho_d2_SPrime[[deriv_1]][[deriv_2]][[i]] <- crossprod(d2_U[[deriv_1]][[deriv_2]], drho_S[[i]][1:p,1:q]) %*% U + 
            crossprod(d_U[[deriv_1]], drho_d_S[[deriv_2]][[i]][1:p,1:q]) %*% U + 
            crossprod(d_U[[deriv_1]], drho_S[[i]][1:p,1:q]) %*% d_U[[deriv_2]] + 
            crossprod(d_U[[deriv_2]], drho_d_S[[deriv_1]][[i]][1:p,1:q]) %*% U + 
            crossprod(U, drho_d2_S[[deriv_1]][[deriv_2]][[i]][1:p,1:q]) %*% U + 
            crossprod(U, drho_d_S[[deriv_1]][[i]][1:p,1:q]) %*% d_U[[deriv_2]] + 
            crossprod(d_U[[deriv_2]], drho_S[[i]][1:p,1:q]) %*% d_U[[deriv_1]] +
            crossprod(U, drho_d_S[[deriv_2]][[i]][1:p,1:q]) %*% d_U[[deriv_1]] +
            crossprod(U, drho_S[[i]][1:p,1:q]) %*% d2_U[[deriv_1]][[deriv_2]]
        }
        
      }
      
      dStPrime[[deriv_1]] <- crossprod(d_U[[deriv_1]], St) %*% U + 
        crossprod(U, d_St[[deriv_1]]) %*% U + 
        crossprod(U, St) %*% d_U[[deriv_1]]
      
      for(i in 1:length(sp))
      {
        d2rho_d_SPrime[[deriv_1]][[i]] <- list()
        
        for(j in i:length(sp))
        {
          d2rho_d_SPrime[[deriv_1]][[i]][[j]] <- crossprod(d_U[[deriv_1]], drho_S[[i]][1:p,1:q]) %*% U + 
            crossprod(U, d2rho_d_S[[deriv_1]][[i]][[j]][1:p,1:q]) %*% U + 
            crossprod(U, drho_S[[i]][1:p,1:q]) %*% d_U[[deriv_1]] 
        }
        drho_d_SPrime[[deriv_1]][[i]] <- crossprod(d_U[[deriv_1]], drho_S[[i]][1:p,1:q]) %*% U + 
          crossprod(U, drho_d_S[[deriv_1]][[i]][1:p,1:q]) %*% U + 
          crossprod(U, drho_S[[i]][1:p,1:q]) %*% d_U[[deriv_1]] 
      }
    }
    
    T_alpha <- matrix(0,nrow = (K+Q), ncol = (K+Q))
    T_gammaPrime <- matrix(0,nrow = (K+Q), ncol = (K+Q)) 
    for(k in 1:length(gammaPrime))
    {
      if(gammaPrime[k] == 1){
        
        for(deriv_1 in 1:length(beta_star))
        {
          for(deriv_2 in deriv_1:length(beta_star))
          {
            d2_Si[[k]][[deriv_1]][[deriv_2]] <- crossprod(d2_Un[[deriv_1]][[deriv_2]],Si[[k]]) %*% Un + 
              crossprod(d_Un[[deriv_1]],d_Si[[k]][[deriv_2]][1:p,1:q]) %*% Un + 
              crossprod(d_Un[[deriv_1]],Si[[k]]) %*% d_Un[[deriv_2]] + 
              crossprod(d_Un[[deriv_2]],d_Si[[k]][[deriv_1]][1:p,1:q]) %*% Un + 
              crossprod(Un,d2_Si[[k]][[deriv_1]][[deriv_2]][1:p,1:q]) %*% Un + 
              crossprod(Un,d_Si[[k]][[deriv_1]][1:p,1:q]) %*% d_Un[[deriv_2]] + 
              crossprod(d_Un[[deriv_2]],Si[[k]]) %*% d_Un[[deriv_1]] +
              crossprod(Un,d_Si[[k]][[deriv_2]][1:p,1:q]) %*% d_Un[[deriv_1]] +
              crossprod(Un,Si[[k]]) %*% d2_Un[[deriv_1]][[deriv_2]]
          }
          
          d_Si[[k]][[deriv_1]] <- crossprod(d_Un[[deriv_1]],Si[[k]]) %*% Un + 
            crossprod(Un,d_Si[[k]][[deriv_1]][1:p,1:q]) %*% Un + 
            crossprod(Un,Si[[k]]) %*% d_Un[[deriv_1]]  
        }
        Si[[k]] <- crossprod(Un,Si[[k]]) %*% Un  
      }
    } 
    
    for(deriv_1 in 1:length(beta_star))
    {
      for(deriv_2 in deriv_1:length(beta_star))
      {
        d2_Qs[[deriv_1]][[deriv_2]] <- d2_U[[deriv_1]][[deriv_2]] %*% Qs + d_U[[deriv_1]] %*% d_Qs[[deriv_2]] + 
          d_U[[deriv_2]] %*% d_Qs[[deriv_1]] + U %*% d2_Qs[[deriv_1]][[deriv_2]]
      }
      d_Qs[[deriv_1]] <- d_U[[deriv_1]] %*% Qs + U %*% d_Qs[[deriv_1]]
    }
    Qs <- U %*% Qs
    
    K <- K + penaltyRank
    Q <- Q - penaltyRank
    St <- StPrime
    d_St <- dStPrime
    d2_St <- d2StPrime
    drho_S <- drho_SPrime
    d2rho_S <- d2rho_SPrime
    drho_d_S <- drho_d_SPrime
    d2rho_d_S <- d2rho_d_SPrime
    drho_d2_S <- drho_d2_SPrime
    d2rho_d2_S <- d2rho_d2_SPrime
    
    gamma <- gammaPrime
  }
  
  
  list(Si = Si, Qs = Qs, d_Qs = d_Qs, d2_Qs = d2_Qs, St = St, d_St = d_St, d2_St = d2_St,
       drho_S = drho_S, d2rho_S = d2rho_S, drho_d_S = drho_d_S, d2rho_d_S = d2rho_d_S, drho_d2_S = drho_d2_S, d2rho_d2_S = d2rho_d2_S
       
  )
}

#08j not used

#Dvariance functions
drho_var_fn <- function(family, mu, drho_mu){
  if(family=="binomial")
  {
    drho_variance <- lapply( drho_mu, FUN=function(input) as.vector(input*(1-2*mu) )  )  
  }
  else if(family=="gaussian")
  {
    drho_variance <- lapply(drho_mu,FUN=function(input) as.vector(0*input) )
  }  
  drho_variance
}

d2rho_var_fn <- function(family, mu, drho_mu_i, drho_mu_j, d2rho_mu){
  if(family=="binomial")
  {
    d2rho_variance <- as.vector( d2rho_mu*(1-2*mu) - 2*drho_mu_i*drho_mu_j )
  }
  else if(family=="gaussian")
  {
    d2rho_variance <- 0*d2rho_mu
  }  
  d2rho_variance
}

d_var_fn <- function(family, mu, d_mu){
  if(family=="binomial")
  {
    d_variance <- lapply( d_mu, FUN=function(input) as.vector(input*(1-2*mu) ) )
  }
  else if(family=="gaussian")
  {
    d_variance <- lapply(d_mu,FUN=function(input) as.vector(0*input))
  }  
  d_variance
}

d2_var_fn <- function(family, mu, d_mu_1, d_mu_2, d2_mu){
  if(family=="binomial")
  {
    d2_variance <- as.vector( d2_mu*(1-2*mu) - 2*d_mu_1*d_mu_2 ) 
  }
  else if(family=="gaussian")
  {
    d2_variance <- 0*d2_mu
  }  
  d2_variance
}

drho_d_var_fn <- function(family, mu, d_mu, drho_mu, drho_d_mu){
  if(family=="binomial")
  {
    drho_d_variance <- as.vector( drho_d_mu*(1-2*mu) -2*d_mu*drho_mu )
  }
  else if(family=="gaussian")
  {
    drho_d_variance <- 0*drho_d_mu
  }  
  drho_d_variance
}

#Corrected version
d2rho_d_var_fn <- function(family, mu, d_mu, drho_mu_i, drho_mu_j, drho_d_mu_i, drho_d_mu_j, d2rho_mu, d2rho_d_mu){
  if(family=="binomial")
  {
    d2rho_d_variance <- as.vector( d2rho_d_mu*(1-2*mu) -2*d2rho_mu*d_mu - 2*drho_d_mu_i*drho_mu_j - 2*drho_mu_i*drho_d_mu_j ) 
  }
  else if(family=="gaussian")
  {
    d2rho_d_variance <- 0*d2rho_d_mu
  }  
  d2rho_d_variance
}

drho_d2_var_fn <- function(family, mu, d_mu_1, d_mu_2, drho_mu, drho_d_mu_1, drho_d_mu_2, d2_mu, drho_d2_mu){
  if(family=="binomial")
  {
    drho_d2_variance <- as.vector( drho_d2_mu*(1-2*mu) -2*d2_mu*drho_mu
                                   -2*drho_d_mu_1*d_mu_2 -2*d_mu_1*drho_d_mu_2 )
  }
  else if(family=="gaussian")
  {
    drho_d2_variance <- 0*drho_d2_mu
  }  
  drho_d2_variance
}

#Corrected version
d2rho_d2_var_fn <- function(family, mu, d_mu_1, d_mu_2, drho_mu_i, drho_mu_j, d2_mu, d2rho_mu, drho_d_mu_i_1, drho_d_mu_i_2, drho_d_mu_j_1, drho_d_mu_j_2,
                            drho_d2_mu_i, drho_d2_mu_j, d2rho_d_mu_1, d2rho_d_mu_2, d2rho_d2_mu){
  if(family=="binomial")
  {
    d2rho_d2_variance <- as.vector( (d2rho_d2_mu*(1-2*mu) + d2rho_d_mu_1*(-2*d_mu_2)) - 
                                      (2*d2rho_d_mu_2*d_mu_1 + 2*d2rho_mu*d2_mu) - 
                                      (2*drho_d2_mu_i*drho_mu_j + 2*drho_d_mu_i_1*drho_d_mu_j_2) - 
                                      (2*drho_d_mu_i_2*drho_d_mu_j_1 + 2*drho_mu_i*drho_d2_mu_j) ) 
  }
  else if(family=="gaussian")
  {
    d2rho_d2_variance <- 0*d2rho_d2_mu
  }  
  d2rho_d2_variance
}

#Dpseudodata functions
drho_pseudo_fn <- function(family, y, inv_variance, mu, drho_mu, drho_eta){
  if(family=="binomial")
  {
    temp <- lapply(drho_mu, FUN=function(input)  -input*inv_variance - ( ( y-mu )*(input-2*mu*input) )*(inv_variance*inv_variance))
    drho_pseudodata <- mapply("+" ,temp ,drho_eta,SIMPLIFY =F)
  }
  else if(family=="gaussian")
  {
    drho_pseudodata <- lapply(drho_mu,FUN=function(input) as.vector(0*input) )
  }  
  drho_pseudodata
}

d2rho_pseudo_fn <- function(family, y, inv_variance, mu, drho_mu_i, drho_mu_j, d2rho_mu, d2rho_eta){
  if(family=="binomial")
  {
    d2rho_pseudodata <-  ( -d2rho_mu*inv_variance ) - 
      
      (-drho_mu_i)*( drho_mu_j-2*mu*drho_mu_j )*(inv_variance*inv_variance)  -   
      
      ( (y-mu)*( d2rho_mu-2*drho_mu_i*drho_mu_j-2*mu*d2rho_mu ) + (-drho_mu_j)*( drho_mu_i-2*mu*drho_mu_i) )*(inv_variance*inv_variance)  + 
      
      (y-mu)*(drho_mu_i-2*mu*drho_mu_i)*2*(drho_mu_j-2*mu*drho_mu_j)*(inv_variance*inv_variance*inv_variance) + 
      
      d2rho_eta
  }
  else if(family=="gaussian")
  {
    d2rho_pseudodata <- 0*d2rho_mu
  }  
  d2rho_pseudodata
}

d_pseudo_fn <- function(family, y, inv_variance, mu, d_mu, d_eta){
  if(family=="binomial")
  {
    temp <- lapply(d_mu, FUN=function(input)  -input*inv_variance - ( ( y-mu )*(input-2*mu*input) )*(inv_variance*inv_variance))
    d_pseudodata <- mapply("+" ,temp ,d_eta,SIMPLIFY =F)
  }
  else if(family=="gaussian")
  {
    d_pseudodata <- lapply(d_mu,FUN=function(input) as.vector(0*input))
  }  
  d_pseudodata
}

d2_pseudo_fn <- function(family, y, inv_variance, variance, d_variance_1, d_variance_2, d2_variance, mu, d_mu_1, d_mu_2, d2_mu, d2_eta){
  if(family=="binomial")
  {
    d2_pseudodata <- -d2_mu*inv_variance + 
      
      d_mu_1*d_variance_2*inv_variance*inv_variance - 
      
      ( (y-mu)*d2_variance - d_mu_2*d_variance_1 )*inv_variance*inv_variance + 
      
      (y-mu)*d_variance_1*2*variance*d_variance_2*inv_variance*inv_variance*inv_variance*inv_variance + 
      
      d2_eta
  }
  else if(family=="gaussian")
  {
    d2_pseudodata <- 0*d2_mu
  }  
  d2_pseudodata
}

drho_d_pseudo_fn <- function(family, y, inv_variance, mu, d_mu, drho_mu, drho_d_mu, drho_d_eta){
  if(family=="binomial")
  {
    drho_d_pseudodata <- -drho_d_mu*inv_variance + 
      
      d_mu*( drho_mu-2*mu*drho_mu )*( inv_variance*inv_variance ) -  
      
      ( (y-mu)*(drho_d_mu -2*drho_mu*d_mu-2*mu*drho_d_mu) - drho_mu*(d_mu-2*mu*d_mu) )*( inv_variance*inv_variance ) +
      
      (y-mu)*(d_mu-2*mu*d_mu)*2*(mu-mu*mu)*(drho_mu-2*mu*drho_mu)*(inv_variance * inv_variance * inv_variance * inv_variance)  + 
      
      drho_d_eta 
  }
  else if(family=="gaussian")
  {
    drho_d_pseudodata <- 0*drho_d_mu
  }  
  drho_d_pseudodata
}

d2rho_d_pseudo_fn <- function(family, y, inv_variance, variance, d_variance, drho_variance_i, drho_variance_j, drho_d_variance_i, drho_d_variance_j, d2rho_variance, d2rho_d_variance,
                              mu, d_mu, drho_mu_i, drho_mu_j, drho_d_mu_i, drho_d_mu_j, d2rho_mu, d2rho_d_mu , d2rho_d_eta){
  if(family=="binomial")
  {
    f1 <- -1*drho_d_mu_i
    d_f1 <- -1*d2rho_d_mu
    
    d_g1 <- drho_variance_j
    f2 <- d_mu*drho_variance_i  
    d_f2 <- drho_d_mu_j*drho_variance_i + d_mu*d2rho_variance 
    
    d_g2 <- 2*variance*drho_variance_j
    f3 <- (y-mu)*drho_d_variance_i - drho_mu_i*d_variance    
    d_f3 <- -drho_mu_j*drho_d_variance_i + (y-mu)*d2rho_d_variance - d2rho_mu*d_variance - drho_mu_i*drho_d_variance_j
    
    d_g3 <- d_g2
    f4 <- (y-mu)*d_variance*2*variance*drho_variance_i
    d_f4 <- -drho_mu_j*d_variance*2*variance*drho_variance_i + (y-mu)*drho_d_variance_j*2*variance*drho_variance_i +
      (y-mu)*d_variance*2*drho_variance_j*drho_variance_i + (y-mu)*d_variance*2*variance*d2rho_variance
    
    d_g4 <- 4*variance*variance*variance*drho_variance_j 
    
    d2rho_d_pseudodata <- ( d_f1*inv_variance - f1*d_g1*inv_variance*inv_variance ) + 
      
      ( d_f2*inv_variance*inv_variance - f2*d_g2*inv_variance*inv_variance*inv_variance*inv_variance ) -
      
      ( d_f3*inv_variance*inv_variance - f3*d_g3*inv_variance*inv_variance*inv_variance*inv_variance ) + 
      
      ( d_f4*inv_variance*inv_variance*inv_variance*inv_variance - f4*d_g4*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance ) + 
      
      d2rho_d_eta
  }
  else if(family=="gaussian")
  {
    d2rho_d_pseudodata <- 0*d2rho_d_mu
  }  
  d2rho_d_pseudodata
}

drho_d2_pseudo_fn <- function(family, y, inv_variance, variance, d_variance_1, d_variance_2, drho_variance, drho_d_variance_1, drho_d_variance_2, d2_variance, drho_d2_variance, 
                              mu, d_mu_1, d_mu_2, drho_mu, drho_d_mu_1, drho_d_mu_2, d2_mu, drho_d2_mu, drho_d2_eta){
  if(family=="binomial")
  {
    f1 <- -drho_d_mu_1
    d_f1 <- -drho_d2_mu
    
    d_g1 <- d_variance_2
    f2 <- d_mu_1*drho_variance
    d_f2 <- d2_mu*drho_variance + d_mu_1*drho_d_variance_2 
    
    d_g2 <- 2*variance*d_variance_2
    f3 <- (y-mu)*drho_d_variance_1 - drho_mu*d_variance_1
    d_f3 <- -d_mu_2*drho_d_variance_1 + (y-mu)*drho_d2_variance - drho_d_mu_2*d_variance_1- drho_mu*d2_variance
    
    d_g3 <- d_g2
    f4 <- (y-mu)*d_variance_1*2*variance*drho_variance  
    d_f4 <- -d_mu_2*d_variance_1*2*variance*drho_variance + (y-mu)*d2_variance*2*variance*drho_variance + (y-mu)*d_variance_1*2*d_variance_2*drho_variance + (y-mu)*d_variance_1*2*variance*drho_d_variance_2
    
    d_g4 <- 4*variance*variance*variance*d_variance_2
    
    drho_d2_pseudodata <-  ( d_f1*inv_variance - f1*d_g1*inv_variance*inv_variance ) +
      
      ( d_f2*inv_variance*inv_variance - f2*d_g2*inv_variance*inv_variance*inv_variance*inv_variance ) - 
      
      ( d_f3*inv_variance*inv_variance - f3*d_g3*inv_variance*inv_variance*inv_variance*inv_variance ) +
      
      ( d_f4*inv_variance*inv_variance*inv_variance*inv_variance - f4*d_g4*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance ) +
      
      drho_d2_eta
  }
  else if(family=="gaussian")
  {
    drho_d2_pseudodata <- 0*drho_d2_mu
  }  
  drho_d2_pseudodata
}

d2rho_d2_pseudo_fn <- function(family, y, inv_variance, variance, drho_variance_i, drho_variance_j, d_variance_1, d_variance_2, drho_d_variance_i_1,
                               drho_d_variance_i_2, drho_d_variance_j_1, drho_d_variance_j_2, d2rho_variance, d2_variance, d2rho_d_variance_1, d2rho_d_variance_2,
                               drho_d2_variance_i, drho_d2_variance_j, d2rho_d2_variance,
                               mu, d_mu_1, d_mu_2, drho_mu_i, drho_mu_j, d2_mu, d2rho_mu, 
                               drho_d_mu_i_1, 
                               drho_d_mu_i_2, drho_d_mu_j_1, drho_d_mu_j_2,
                               drho_d2_mu_i, drho_d2_mu_j, d2rho_d_mu_1, d2rho_d_mu_2, d2rho_d2_mu, d2rho_d2_eta){
  if(family=="binomial")
  {
    
    f11 <- -drho_d2_mu_i
    drho_f11 <- -d2rho_d2_mu
    
    drho_g11 <- drho_variance_j
    f12 <- -drho_d_mu_i_1*d_variance_2
    drho_f12 <- -d2rho_d_mu_1*d_variance_2 -drho_d_mu_i_1*drho_d_variance_j_2
    
    drho_g12 <- 2*variance*drho_variance_j
    f21 <-  d2_mu*drho_variance_i + d_mu_1*drho_d_variance_i_2
    drho_f21 <- drho_d2_mu_j*drho_variance_i + d2_mu*d2rho_variance + drho_d_mu_j_1*drho_d_variance_i_2 + d_mu_1*d2rho_d_variance_2
    
    drho_g21 <- 2*variance*drho_variance_j
    f22 <- d_mu_1*drho_variance_i * 2*variance*d_variance_2
    
    drho_f22 <- drho_d_mu_j_1*drho_variance_i*2*variance*d_variance_2 + d_mu_1*d2rho_variance*2*variance*d_variance_2 + d_mu_1*drho_variance_i*2*drho_variance_j*d_variance_2 + d_mu_1*drho_variance_i*2*variance*drho_d_variance_j_2
    
    drho_g22 <- 4*variance*variance*variance*drho_variance_j
    
    f31 <- -d_mu_2*drho_d_variance_i_1 + (y-mu)*drho_d2_variance_i - drho_d_mu_i_2*d_variance_1 - drho_mu_i*d2_variance
    drho_f31 <- -drho_d_mu_j_2*drho_d_variance_i_1 - d_mu_2*d2rho_d_variance_1 + -drho_mu_j*drho_d2_variance_i + (y-mu)*d2rho_d2_variance - 
      d2rho_d_mu_2*d_variance_1 - drho_d_mu_i_2*drho_d_variance_j_1 - d2rho_mu*d2_variance- drho_mu_i*drho_d2_variance_j 
    
    drho_g31 <- 2*variance*drho_variance_j
    
    d_g3_temp <- 2*variance*d_variance_2
    drho_d_g3_temp <- 2*drho_variance_j*d_variance_2+ 2*variance*drho_d_variance_j_2
    
    f3_temp <- (y-mu)*drho_d_variance_i_1 - drho_mu_i*d_variance_1
    drho_f3_temp <- -drho_mu_j*drho_d_variance_i_1 + (y-mu)*d2rho_d_variance_1 - d2rho_mu*d_variance_1 - drho_mu_i*drho_d_variance_j_1
    
    f32 <- ( (y-mu)*drho_d_variance_i_1 - drho_mu_i*d_variance_1 ) * 2*variance*d_variance_2
    drho_f32 <- drho_f3_temp*d_g3_temp + f3_temp*drho_d_g3_temp
    
    drho_g32 <- 4*variance*variance*variance*drho_variance_j 
    f41 <-  -d_mu_2*d_variance_1*2*variance*drho_variance_i + (y-mu)*d2_variance*2*variance*drho_variance_i + (y-mu)*d_variance_1*2*d_variance_2*drho_variance_i + (y-mu)*d_variance_1*2*variance*drho_d_variance_i_2 
    
    drho_f41 <- -drho_d_mu_j_2*d_variance_1*2*variance*drho_variance_i - d_mu_2*drho_d_variance_j_1*2*variance*drho_variance_i + -d_mu_2*d_variance_1*2*drho_variance_j*drho_variance_i - 
      d_mu_2*d_variance_1*2*variance*d2rho_variance + -drho_mu_j*d2_variance*2*variance*drho_variance_i + (y-mu)*drho_d2_variance_j*2*variance*drho_variance_i +  
      (y-mu)*d2_variance*2*drho_variance_j*drho_variance_i + (y-mu)*d2_variance*2*variance*d2rho_variance + -drho_mu_j*d_variance_1*2*d_variance_2*drho_variance_i + 
      (y-mu)*drho_d_variance_j_1*2*d_variance_2*drho_variance_i + (y-mu)*d_variance_1*2*drho_d_variance_j_2*drho_variance_i + (y-mu)*d_variance_1*2*d_variance_2*d2rho_variance +
      -drho_mu_j*d_variance_1*2*variance*drho_d_variance_i_2 + (y-mu)*drho_d_variance_j_1*2*variance*drho_d_variance_i_2 + (y-mu)*d_variance_1*2*drho_variance_j*drho_d_variance_i_2 + (y-mu)*d_variance_1*2*variance*d2rho_d_variance_2 
    
    drho_g41 <- 4*variance*variance*variance*drho_variance_j
    
    f42 <- (y-mu)*d_variance_1*2*variance*drho_variance_i*4*variance*variance*variance*d_variance_2
    
    drho_f42 <- -drho_mu_j*d_variance_1*2*variance*drho_variance_i*4*variance*variance*variance*d_variance_2 +  (y-mu)*drho_d_variance_j_1*2*variance*drho_variance_i*4*variance*variance*variance*d_variance_2+
      (y-mu)*d_variance_1*2*drho_variance_j*drho_variance_i*4*variance*variance*variance*d_variance_2 + (y-mu)*d_variance_1*2*variance*d2rho_variance*4*variance*variance*variance*d_variance_2+
      (y-mu)*d_variance_1*2*variance*drho_variance_i*4*3*variance*variance*drho_variance_j*d_variance_2 + (y-mu)*d_variance_1*2*variance*drho_variance_i*4*variance*variance*variance*drho_d_variance_j_2 
    
    drho_g42 <- 8*variance*variance*variance*variance*variance*variance*variance*drho_variance_j 
    
    d2rho_d2_pseudodata <- ( ( drho_f11*inv_variance - f11*drho_g11*inv_variance*inv_variance ) - 
                               
                               ( drho_f12*inv_variance*inv_variance - f12*drho_g12*inv_variance*inv_variance*inv_variance*inv_variance ) ) + 
      
      ( ( drho_f21*inv_variance*inv_variance - f21*drho_g21*inv_variance*inv_variance*inv_variance*inv_variance ) - 
          
          ( drho_f22*inv_variance*inv_variance*inv_variance*inv_variance - 
              
              f22*drho_g22*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance ) ) -              
      
      ( ( drho_f31*inv_variance*inv_variance - f31*drho_g31*inv_variance*inv_variance*inv_variance*inv_variance ) - 
          
          ( drho_f32*inv_variance*inv_variance*inv_variance*inv_variance - 
              
              f32*drho_g32*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance ) ) +              
      
      ( ( drho_f41*inv_variance*inv_variance*inv_variance*inv_variance - 
            
            f41*drho_g41*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance ) - 
          
          ( drho_f42*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance - 
              
              f42*drho_g42*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance*inv_variance ) ) + 
      
      d2rho_d2_eta 
    
  }
  else if(family=="gaussian")
  {
    d2rho_d2_pseudodata <- 0*d2rho_d2_mu
  }  
  d2rho_d2_pseudodata
}

#Ddeviance functions
d_deviance_fn <- function(family, y, mu, d_mu){
  if(family=="binomial")
  {
    d_dev <- unlist(lapply(d_mu, FUN=function(input) -2*sum( ( y*input )/mu  - ( (1-y)*input )/(1-mu) ) ))
  }
  else if(family=="gaussian")
  {
    d_dev <- unlist(lapply(d_mu, FUN=function(input) sum( 2*(y-mu)*(-input) ) ) )
  }  
  d_dev
}

d2_deviance_fn <- function(family, y, mu, d_mu_1, d_mu_2, d2_mu){
  if(family=="binomial")
  {
    d2_dev <- -2*sum(
      y*(d2_mu/mu - (d_mu_1*d_mu_2)/(mu*mu) )
      -
        (1-y)*(d2_mu/(1-mu) - (d_mu_1*(-d_mu_2))/((1-mu)*(1-mu)) )
    )
  }
  else if(family=="gaussian")
  {
    d2_dev <- sum( 2*(-d_mu_2)*(-d_mu_1) + 2*(y-mu)*(-d2_mu) )
  }  
  d2_dev
}

drho_d_deviance_fn <- function(family, y, mu, d_mu, drho_mu, drho_d_mu){
  if(family=="binomial")
  {
    drho_d_dev <- -2*sum( 
      y*( drho_d_mu/mu - d_mu*drho_mu/(mu*mu) ) -
        
        (1-y)*( drho_d_mu/(1-mu) - d_mu*(-drho_mu)/((1-mu)*(1-mu)) ) 
    )
  }
  else if(family=="gaussian")
  {
    drho_d_dev <- sum(  2*(-drho_mu)*(-d_mu) + 2*(y-mu)*(-drho_d_mu) )
  }  
  drho_d_dev
}

d2rho_d_deviance_fn <- function(family, y, mu, d_mu, drho_mu_i, drho_mu_j, d2rho_mu, drho_d_mu_i, drho_d_mu_j, d2rho_d_mu){
  if(family=="binomial")
  {
    d2rho_d_dev <- -2*sum( 
      y*( 
        d2rho_d_mu/mu - 
          drho_d_mu_i*drho_mu_j/(mu*mu)  -
          
          (drho_d_mu_j*drho_mu_i + d_mu*d2rho_mu)/(mu*mu)  + 
          d_mu*drho_mu_i*2*mu*drho_mu_j/(mu*mu*mu*mu)        
      ) -
        
        (1-y)*( 
          d2rho_d_mu/(1-mu) - 
            drho_d_mu_i*(-drho_mu_j)/((1-mu)*(1-mu)) -
            
            (drho_d_mu_j*(-drho_mu_i) + d_mu*(-d2rho_mu))/((1-mu)*(1-mu)) + 
            d_mu*(-drho_mu_i)*2*(1-mu)*(-drho_mu_j)/((1-mu)*(1-mu)*(1-mu)*(1-mu))  
        ) 
    )
    
  }
  else if(family=="gaussian")
  {
    d2rho_d_dev <- sum(  
      2*(-d2rho_mu)*(-d_mu) + 2*(-drho_mu_i)*(-drho_d_mu_j) + 
        2*(-drho_mu_j)*(-drho_d_mu_i) + 2*(y-mu)*(-d2rho_d_mu) 
    )
  }  
  d2rho_d_dev
}

drho_d2_deviance_fn <- function(family, y, mu, d_mu_1, d_mu_2, drho_mu, drho_d_mu_1, drho_d_mu_2, d2_mu, drho_d2_mu){
  if(family=="binomial")
  {
    drho_d2_dev <- -2*sum(
      y*(
        drho_d2_mu/mu - 
          d2_mu*drho_mu/(mu*mu) -
          
          (drho_d_mu_1*d_mu_2 + d_mu_1*drho_d_mu_2)/(mu*mu) +
          d_mu_1*d_mu_2*2*mu*drho_mu/(mu*mu*mu*mu)
      )
      - 
        (1-y)*(
          drho_d2_mu/(1-mu) -
            d2_mu*(-drho_mu)/((1-mu)*(1-mu)) -
            
            (drho_d_mu_1*(-d_mu_2) + d_mu_1*(-drho_d_mu_2))/((1-mu)*(1-mu)) +
            d_mu_1*(-d_mu_2)*2*(1-mu)*(-drho_mu)/((1-mu)*(1-mu)*(1-mu)*(1-mu))
        )
    )
    
  }
  else if(family=="gaussian")
  {
    drho_d2_dev <- sum( 2*(-drho_d_mu_2)*(-d_mu_1) + 2*(-d_mu_2)*(-drho_d_mu_1) + 
                          2*(-drho_mu)*(-d2_mu) + 2*(y-mu)*(-drho_d2_mu) )
  }  
  drho_d2_dev
}

d2rho_d2_deviance_fn <- function(family, y, mu, d_mu_1, d_mu_2, drho_mu_i, drho_mu_j, d2rho_mu, d2_mu, drho_d_mu_i_1, drho_d_mu_i_2, drho_d_mu_j_1, drho_d_mu_j_2,
                                 drho_d2_mu_i, drho_d2_mu_j, d2rho_d_mu_1, d2rho_d_mu_2, d2rho_d2_mu,
                                 d_eta_1, d_eta_2, drho_eta_i, drho_eta_j, d2rho_eta, d2_eta, drho_d_eta_i_1, drho_d_eta_i_2, drho_d_eta_j_1, drho_d_eta_j_2,
                                 drho_d2_eta_i, drho_d2_eta_j, d2rho_d_eta_1, d2rho_d_eta_2, d2rho_d2_eta ){
  
  if(family=="binomial")
  {
    d2rho_d2_dev <- -2*sum( 
      y*(
        d2rho_d2_eta*(1-mu) + d2rho_d_eta_1*(-d_mu_2) + 
          drho_d2_eta_i*(-drho_mu_j) + drho_d_eta_i_1*(-drho_d_mu_j_2) +
          drho_d2_eta_j*( -drho_mu_i ) + drho_d_eta_j_1*( -drho_d_mu_i_2 ) +                                      
          d2_eta*( -d2rho_mu ) + d_eta_1*( -d2rho_d_mu_2 ))
      +
        (1-y)*( 
          (-d2rho_d2_eta)*mu + (-d2rho_d_eta_1)*d_mu_2 +
            (-drho_d2_eta_i)*drho_mu_j + (-drho_d_eta_i_1)*drho_d_mu_j_2 +
            (-drho_d2_eta_j)*( drho_mu_i ) + (-drho_d_eta_j_1)*( drho_d_mu_i_2 ) +                                              
            (-d2_eta)*( d2rho_mu )  + (-d_eta_1)*( d2rho_d_mu_2 ) )                                                
      
    )
    
  }
  else if(family=="gaussian")
  {
    d2rho_d2_dev <- sum( 2*(-d2rho_d_mu_2)*(-d_mu_1) + 2*(-drho_d_mu_i_2)*(-drho_d_mu_j_1) + 
                           2*(-drho_d_mu_j_2)*(-drho_d_mu_i_1) + 2*(-d_mu_2)*(-d2rho_d_mu_1) + 
                           2*(-d2rho_mu)*(-d2_mu) + 2*(-drho_mu_i)*(-drho_d2_mu_j) + 
                           2*(-drho_mu_j)*(-drho_d2_mu_i) + 2*(y-mu)*(-d2rho_d2_mu) 
    )
  }  
  d2rho_d2_dev
}

#Main functions
smooth.construct.pb.smooth.spec <- function (object, data, knots) {
  if (length(object$p.order) == 1) 
    m <- rep(object$p.order, 2)
  else m <- object$p.order
  m[is.na(m)] <- 2
  object$p.order <- m
  if (object$bs.dim < 0) 
    object$bs.dim <- max(10, m[1] + 1)
  nk <- object$bs.dim - m[1]
  if (nk <= 0) 
    stop("basis dimension too small for b-spline order")
  if (length(object$term) != 1) 
    stop("Basis only handles 1D smooths")
  x <- data[[object$term]]
  k <- knots[[object$term]]
  
  #Set global variables for gradient testing
  
  if (is.null(k)) {
    xl <- min(x) 
    xu <- max(x)
  }
  else if (length(k) == 2) {
    xl <- min(x) 
    xu <- max(x)
    if (xl > min(x) || xu < max(x)) 
      stop("knot range does not include data")
  }
  if (is.null(k) || length(k) == 2) {
    xr <- xu - xl
    xl <- xl - xr * 0.001
    xu <- xu + xr * 0.001
    dx <- (xu - xl)/(nk - 1)
    k <- seq(xl - dx * (m[1] + 1), xu + dx * (m[1] + 1), 
             length = nk + 2 * m[1] + 2)
  }
  else {
    if (length(k) != nk + 2 * m[1] + 2) 
      stop(paste("there should be ", nk + 2 * m[1] + 2, 
                 " supplied knots"))
  }
  object$X <- splines::spline.des(k, x, m[1] + 2, x * 0, outer.ok=T)$design
  if (!is.null(k)) {
    if (sum(colSums(object$X) == 0) > 0) 
      warning("knot range is so wide that there is *no* information about some basis coefficients")
  }
  if (length(unique(x)) < object$bs.dim) 
    warning("basis dimension is larger than number of unique covariates")
  S <- diag(object$bs.dim)
  if (m[2]) 
    for (i in 1:m[2]) S <- diff(S)
  object$S <- list(t(S) %*% S)
  object$S[[1]] <- (object$S[[1]] + t(object$S[[1]]))/2
  object$rank <- object$bs.dim - m[2]
  object$null.space.dim <- m[2]
  object$knots <- k
  object$nk <- nk
  object$m <- m
  class(object) <- "pb.smooth"
  object
}

Predict.matrix.pb.smooth <- function(object,data){ 
  x <- data[[object$term]]
  m <- object$m;     # spline order (3=cubic)
  #k<-object$knots    # knot locations 
  nk <- object$nk
  xl <- min(x)
  xu <- max(x)
  xr <- xu - xl
  xl <- xl - xr * 0.001
  xu <- xu + xr * 0.001
  dx <- (xu - xl)/(nk - 1)
  k <- seq(xl - dx * (m[1] + 1), xu + dx * (m[1] + 1), 
           length = nk + 2 * m[1] + 2)
  X <- splines::spline.des(k, x, m[1] + 2, x * 0,outer.ok=T)$design
  X
}

Dbeta <- function(beta_star){ #The argument it takes is beta rather than beta_star. We need dbeta_star_i/d_beta_i, where beta_star = f(beta)
  d<-(length(beta_star)+1)
  dbeta<-matrix(nrow=(d),ncol=(d-1))
  #we have all_derivs[[1]] as zero vector 
  beta <- c(sqrt( 1-sum((beta_star/sqrt((1+sum(beta_star^2))))^2) ), beta_star/sqrt((1+sum(beta_star^2))))
  
  for(t in 2:d) #Differentiate with respect to beta_star[t]
  {
    for(i in 2:d) #p indexes through the elements of beta_1, beta_2, ... , beta_d
    {
      if(i!=1 & i==t)
      {
        dbeta[i,(t-1)] <- 1/sqrt(1 + sum(beta_star^2)) - beta_star[i-1]*beta_star[i-1]/((1 + sum(beta_star^2))^(3/2))
      }
      else if(i!=1 & i!=t)
      {
        dbeta[i,(t-1)] <- - beta_star[i-1]*beta_star[t-1]/((1 + sum(beta_star^2))^(3/2))
      }
    }
    dbeta[1,(t-1)] <- -( sum(2*beta[2:d]*dbeta[2:d,(t-1)]) )/( 2*sqrt( 1-sum( beta[2:d]^2 ) ) )
  }
  dbeta
}

D2beta <- function(beta_star,dbeta){ #The argument it takes is beta rather than beta_star. We need dbeta_star_i/d_beta_i, where beta_star = f(beta)
  d<-(length(beta_star)+1)
  d2beta<-array(dim=c((d),(d-1),(d-1)) )
  #we have all_derivs[[1]] as zero vector 
  beta <- c(sqrt( 1-sum((beta_star/sqrt((1+sum(beta_star^2))))^2) ), beta_star/sqrt((1+sum(beta_star^2))))
  
  for(alpha in 2:d)
  {
    for(sigma in 2:d)
    {
      for(t in 2:d)
      {
        if(t!=1 & t==alpha)
        {
          if(alpha==sigma)
          {
            d2beta[t,(alpha-1),(sigma-1)] <- -beta_star[(sigma-1)]*(1 + sum(beta_star^2))^(-3/2) - 
              
              2*beta_star[t-1]/((1 + sum(beta_star^2))^(3/2)) +
              
              beta_star[t-1]*beta_star[t-1]*3*beta_star[sigma-1]/((1 + sum(beta_star^2))^(5/2))
          }
          else if (alpha != sigma)
          {
            d2beta[t,(alpha-1),(sigma-1)] <- -beta_star[(sigma-1)]*(1 + sum(beta_star^2))^(-3/2) + 
              
              beta_star[t-1]*beta_star[t-1]*3*beta_star[sigma-1]/((1 + sum(beta_star^2))^(5/2))
          }
          
        }
        else if(t!=1 & t!=alpha)
        {
          if(alpha == sigma & sigma == t)
          {
            d2beta[t,(alpha-1),(sigma-1)] <- - ( 
              (beta_star[(alpha-1)] + beta_star[(t-1)])/((1 + sum(beta_star^2))^(3/2))  -
                
                beta_star[(t-1)]*beta_star[(alpha-1)]*beta_star[(sigma-1)]*3/((1 + sum(beta_star^2))^(5/2))
            )
          }
          else if(alpha == sigma & sigma != t)
          {
            d2beta[t,(alpha-1),(sigma-1)] <- - ( 
              beta_star[t-1]/((1 + sum(beta_star^2))^(3/2))  -
                
                beta_star[t-1]*beta_star[alpha-1]*beta_star[sigma-1]*3/((1 + sum(beta_star^2))^(5/2))
            )
          }
          else if (alpha != sigma & sigma == t)
          {
            d2beta[t,(alpha-1),(sigma-1)] <- - ( 
              beta_star[alpha-1]/((1 + sum(beta_star^2))^(3/2))  -
                
                beta_star[t-1]*beta_star[alpha-1]*beta_star[sigma-1]*3/((1 + sum(beta_star^2))^(5/2))
            )
          }
          else if (alpha != sigma & sigma != t)
          {
            d2beta[t,(alpha-1),(sigma-1)] <- - ( 
              - beta_star[t-1]*beta_star[alpha-1]*beta_star[sigma-1]*3/((1 + sum(beta_star^2))^(5/2))
            )
          }
          
        }
      }
      
      d2beta[1,(alpha-1),(sigma-1)] <- -( sum(dbeta[2:d,(sigma-1)]*dbeta[2:d,(alpha-1)] + beta[2:d]* d2beta[2:d,(alpha-1),(sigma-1)]) )/( sqrt( 1-sum( beta[2:d]^2 ) ) ) -
        
        ( -sum(beta[2:d]*dbeta[2:d,(alpha-1)]) * -sum(beta[2:d]*dbeta[2:d,(sigma-1)])  )/( ( 1-sum( beta[2:d]^2 ) )^(3/2) )
      
    }
  }  
  d2beta
}

#Updated version 
QR.Differentiation <- function(Deriv_Design, Design, Deriv_B, B){
  nr <- nrow(Design)
  nc <- ncol(Design)
  q <- ncol(B)
  p <- nrow(B)
  #The first iteration of the householder matrices, H, is H^0 = Deriv_Design, namely the matrix of element wise derivatives of the spline 
  #transformed design matrix, dependent on the index parameter vector beta
  
  #At each iteration of the recursive differentation, we need the householder matrix and the matrix of derivatives from the iteration before.
  #Thus, to initialize, we need to pass this function the design matrix and the matrix of derivatives. In my case, I need the spline transformation 
  #matrix as well as the matrix of the derivatives of the spline transformation
  
  #H^0 = Design
  #D(H^0) = Deriv_Design
  
  #Deriv_B will be the matrix of zeros when B is the identity - B being the identity gives us Q as the result
  
  H <- Design
  D_H <- Deriv_Design
  
  
  # Householder transformations
  for (j in seq_len(nc))
  {
    H_copy <- H
    t <- seq.int(j, nr)
    t2 <- seq.int(j,p)
    sigma <- sum(H[t,j]^2)
    D_sigma <- 2*sum(D_H[t,j]*H[t,j])
    s <- sqrt(sigma)
    D_s <- D_sigma * 0.5 * sigma^(-0.5)
    diag_ej <- H[j,j]
    #note that D_diag_ej <- D_H[j,j]
    D_diag_ej<- D_H[j,j]
    gamma <- 1.0 / (sigma + abs(s * diag_ej))
    ##We differentiate the denmominator of gamma which we call D_gamma_denom
    D_gamma_denom <- D_sigma + (( (D_s * diag_ej) + (s * D_diag_ej) )*(s * diag_ej))/abs(s * diag_ej)
    #D_gamma <- -( D_sigma + ((s*diag_ej)/abs(s * diag_ej))*(D_s*diag_ej + s*D_diag_ej  )  ) / ( (sigma + abs(s * diag_ej))^2 )    
    kappa <- if (diag_ej < 0) s else -s
    D_kappa <- if (diag_ej < 0) D_s else -D_s
    
    g <- 1/gamma   
    D_g <- D_sigma + ((D_diag_ej *s + diag_ej * D_s) * diag_ej * s)/abs(s * diag_ej)
    q <- seq.int(j+1, nr)
    
    H[j,j] <- H[j,j] - kappa
    if (j < nc)
      for (k in seq.int(j+1, nc))
      {
        yPrime <- sum(H[t,j] * H[t,k]) * gamma
        H[t,k] <- H[t,k] - H[t,j] * yPrime             
      }
    
    #######We add this additional loop which shall return Q      
    #Note that at this stage, H has been updated to the extent that the H[j,j]=H[j,j]-kappa (so not their final forMulation) however the rest
    #of H has the final forMulation
    
    #Note that at this stage, D_H will need to have been updated
    
    t3 <- seq.int(j+1,p)   
    for (k2 in 1:ncol(B))
    {
      yPrime <- sum(H[t2,j] * B[t2,k2]) * gamma      
      f_yPrime <- sum(H[t2,j] * B[t2,k2])       
      D_f_yPrime <- sum( (D_H[t3,j] * B[t3,k2]) + (H[t3,j] * Deriv_B[t3,k2]) ) + ( (D_diag_ej - D_kappa) * B[j,k2] ) + ( (diag_ej-kappa) * Deriv_B[j,k2] )
      
      D_yPrime <- D_f_yPrime/(1/gamma) - (f_yPrime * D_gamma_denom)/((1/gamma)^2)
      
      #t3=j
      Deriv_B[j,k2] <- Deriv_B[j,k2] - (D_diag_ej-D_kappa)*yPrime - (diag_ej-kappa)*D_yPrime
      #t3 =/= j
      Deriv_B[t3,k2] <- Deriv_B[t3,k2]- D_H[t3,j]*yPrime - H[t3,j]*D_yPrime
      
      B[t2,k2] <- B[t2,k2] - H[t2,j] * yPrime
      
    }
    
    
    #This Must be done after computing Q (B)
    H[j,j] <- kappa
    
    #This gives me the householder matrix for the current iteration, I just need to calculate the householder matrix derivatives for this iteration
    #We retain the same value of j
    #We shall calculate the components seperately, particularly for the two instances of the quotient rule we need to calculate
    #Need to remember that in my written derivatives I differentiated d/dbeta^star{H^(j+1)}, where as here I'm using H^(j), so need to replace 
    #j+1's with j's
    
    #j: 1 to nc
    #t: j to nc
    #k: j+1 to nc
    #p: j to nc <=> p: t to nc
    #q: j to nc <=> q: t to nc
    
    if (j < nc)
      for (k in seq.int(j+1, nc))
      { 
        #use f_1 when t=j
        f1 <- (diag_ej-kappa) * ( sum(H_copy[q,j]*H_copy[q,k]) + (diag_ej-kappa)*H_copy[j,k] )
        D_f1 <- (D_diag_ej-D_kappa) * ( sum(H_copy[q,j]*H_copy[q,k]) + (diag_ej-kappa)*H_copy[j,k] )   + (diag_ej-kappa) * (  sum( D_H[q,j]*H_copy[q,k] + H_copy[q,j]*D_H[q,k] ) + (D_diag_ej-D_kappa) * H_copy[j,k]  + (diag_ej-kappa)*D_H[j,k] )       
        #use f_2 when t>j
        f2 <- H_copy[q,j]* ( sum(H_copy[q,j]*H_copy[q,k]) + (diag_ej-kappa)*H_copy[j,k] )
        D_f2 <- D_H[q,j] * ( sum(H_copy[q,j]*H_copy[q,k]) + (diag_ej-kappa)*H_copy[j,k] ) + H_copy[q,j] * (  sum( D_H[q,j]*H_copy[q,k] + H_copy[q,j]*D_H[q,k] ) + (D_diag_ej-D_kappa) * H_copy[j,k]  + (diag_ej-kappa)*D_H[j,k] ) 
        
        #when t=j
        D_H[j,k] <- D_H[j,k] - ( D_f1/g - (f1*D_g)/(g^2) )
        #when t>j
        D_H[q,k] <- D_H[q,k] - ( D_f2/g - (f2*D_g)/(g^2) )
        
      }
    D_H[j,j] <- D_kappa      
  } # end Householder
  
  # set zeros in the lower triangular side of X (which stores) 
  # not really necessary, this is just to return R for illustration
  for (i in seq_len(ncol(H)))
    H[seq.int(i+1, nr),i] <- 0
  
  for (i in seq_len(ncol(D_H)))
    D_H[seq.int(i+1, nr),i] <- 0 
  
  list(R=H[1:nc,1:nc], D_R=D_H[1:nc,1:nc],B=B, D_B= Deriv_B)  
}

bspline <- function(x,k,i,m=2){ if (m==-1) # base of recursion
{ 
  res <- as.numeric(x<k[i+1]&x>=k[i])
} else{ # construct from call to lower order basis
  z0 <- (x-k[i])/(k[i+m+1]-k[i])
  z1 <- (k[i+m+2]-x)/(k[i+m+2]-k[i+1])
  res <- z0*bspline(x,k,i,m-1)+ z1*bspline(x,k,i+1,m-1)
}
  res
}

Dbspline <- function(x,dx,k,dk,i,m=2){
  if (m==-1) # base of recursion
  { 
    dres <- rep(0,times=length(x))
  }
  else
  {
    z0 <- (x-k[i])/(k[i+m+1]-k[i])
    z1 <- (k[i+m+2]-x)/(k[i+m+2]-k[i+1])
    
    dz0 <- (dx-dk[i])/(k[i+m+1]-k[i]) - ((x-k[i])*(dk[i+m+1]-dk[i]))/((k[i+m+1]-k[i])^2)
    dz1 <- (dk[i+m+2]-dx)/(k[i+m+2]-k[i+1]) - ((k[i+m+2]-x)*(dk[i+m+2]-dk[i+1]))/((k[i+m+2]-k[i+1])^2)
    
    dres <- dz0*bspline(x=x,k=k,i=i,m=(m-1)) + z0*Dbspline(x=x,dx=dx,k=k,dk=dk,i=i,m=(m-1)) + dz1*bspline(x=x,k,i=(i+1),m=(m-1)) + z1*Dbspline(x=x,dx=dx,k=k,dk=dk,i=(i+1),m=(m-1))
  }  
  dres
}

D2bspline <- function(x,dx,d2x,k,dk,d2k,i,m=2,deriv_1,deriv_2){  
  if (m==-1) # base of recursion
  { 
    d2res <- rep(0,times=length(x))
  }
  else
  {
    z0 <- (x-k[i])/(k[i+m+1]-k[i])
    z1 <- (k[i+m+2]-x)/(k[i+m+2]-k[i+1])
    
    dz0_deriv_1 <- (dx[,deriv_1]-dk[,deriv_1][i])/(k[i+m+1]-k[i]) - ((x-k[i])*(dk[,deriv_1][i+m+1]-dk[,deriv_1][i]))/((k[i+m+1]-k[i])^2)
    dz1_deriv_1 <- (dk[,deriv_1][i+m+2]-dx[,deriv_1])/(k[i+m+2]-k[i+1]) - ((k[i+m+2]-x)*(dk[,deriv_1][i+m+2]-dk[,deriv_1][i+1]))/((k[i+m+2]-k[i+1])^2)
    
    
    dz0_deriv_2 <- (dx[,deriv_2]-dk[,deriv_2][i])/(k[i+m+1]-k[i]) - ((x-k[i])*(dk[,deriv_2][i+m+1]-dk[,deriv_2][i]))/((k[i+m+1]-k[i])^2)
    dz1_deriv_2 <- (dk[,deriv_2][i+m+2]-dx[,deriv_2])/(k[i+m+2]-k[i+1]) - ((k[i+m+2]-x)*(dk[,deriv_2][i+m+2]-dk[,deriv_2][i+1]))/((k[i+m+2]-k[i+1])^2)
    
    f1z0 <- dx[,deriv_1]-dk[,deriv_1][i]
    df1z0 <- d2x[[deriv_1]][,deriv_2] - d2k[[deriv_1]][,deriv_2][i]
    
    g1z0 <- k[i+m+1]-k[i]
    dg1z0 <- dk[,deriv_2][i+m+1] - dk[,deriv_2][i]
    
    f2z0 <- ((x-k[i])*(dk[,deriv_1][i+m+1]-dk[,deriv_1][i]))
    df2z0 <- ((dx[,deriv_2]-dk[,deriv_2][i])*(dk[,deriv_1][i+m+1]-dk[,deriv_1][i])) + ((x-k[i])*(d2k[[deriv_1]][,deriv_2][i+m+1]-d2k[[deriv_1]][,deriv_2][i]))
    
    g2z0 <- (k[i+m+1]-k[i])^2
    dg2z0 <- 2*(k[i+m+1]-k[i])*(dk[,deriv_2][i+m+1]-dk[,deriv_2][i])
    ######
    f1z1 <- dk[,deriv_1][i+m+2]-dx[,deriv_1]
    df1z1 <- d2k[[deriv_1]][,deriv_2][i+m+2] - d2x[[deriv_1]][,deriv_2]
    
    g1z1 <- k[i+m+2]-k[i+1]
    dg1z1 <- dk[,deriv_2][i+m+2] - dk[,deriv_2][i+1]
    
    f2z1 <- ((k[i+m+2]-x)*(dk[,deriv_1][i+m+2]-dk[,deriv_1][i+1]))
    df2z1 <- ((dk[,deriv_2][i+m+2]-dx[,deriv_2])*(dk[,deriv_1][i+m+2]-dk[,deriv_1][i+1])) + ((k[i+m+2]-x)*(d2k[[deriv_1]][,deriv_2][i+m+2]-d2k[[deriv_1]][,deriv_2][i+1]))
    
    g2z1 <- (k[i+m+2]-k[i+1])^2
    dg2z1 <- 2*(k[i+m+2]-k[i+1])*(dk[,deriv_2][i+m+2]-dk[,deriv_2][i+1])
    
    d2z0 <- (df1z0/g1z0) - ((f1z0*dg1z0)/(g1z0^2)) - ( (df2z0/g2z0) - ((f2z0*dg2z0)/(g2z0^2)) )
    d2z1 <- (df1z1/g1z1) - ((f1z1*dg1z1)/(g1z1^2)) - ( (df2z1/g2z1) - ((f2z1*dg2z1)/(g2z1^2)) )
    
    
    d2res <- (d2z0*bspline(x=x,k=k,i=i,m=(m-1)) + dz0_deriv_1*Dbspline(x=x,dx=dx[,deriv_2],k=k,dk=dk[,deriv_2],i=i,m=(m-1)) +
                dz0_deriv_2*Dbspline(x=x,dx=dx[,deriv_1],k=k,dk=dk[,deriv_1],i=i,m=(m-1)) + z0*D2bspline(x=x,dx=dx,d2x=d2x,k=k,dk=dk,d2k=d2k,i=i,m=(m-1),deriv_1=deriv_1,deriv_2=deriv_2) + 
                d2z1*bspline(x=x,k,i=(i+1),m=(m-1)) + dz1_deriv_1*Dbspline(x=x,dx=dx[,deriv_2],k,dk=dk[,deriv_2],i=(i+1),m=(m-1)) +
                dz1_deriv_2*Dbspline(x=x,dx=dx[,deriv_1],k=k,dk=dk[,deriv_1],i=(i+1),m=(m-1)) + z1*D2bspline(x=x,dx=dx,d2x=d2x,k=k,dk=dk,d2k=d2k,i=(i+1),m=(m-1),deriv_1=deriv_1,deriv_2=deriv_2) 
    )
  }  
  d2res
}

createKnots <- function(x, spline_rank, spline_degree){
  nk <- spline_rank - spline_degree
  xu <- max(x)
  xl <- min(x)
  xr <- xu - xl
  xl <- xl - xr * 0.001
  xu <- xu + xr * 0.001
  dx <- (xu - xl)/(nk - 1)
  knots <- seq(xl - dx * (2 + 1), xu + dx * (2 + 1), length = nk + 2 * 2 + 2)
  return(knots)
}

createDknots <- function(beta_star, x, dx, spline_rank, spline_degree, knots){
  
  dxl <- dx[which.min(x)]
  dxu <- dx[which.max(x)]
  ndx <- spline_rank-spline_degree
  bdeg <- spline_degree
  
  dxr <- dxu - dxl 
  dxl <- dxl - dxr * 0.001
  dxu <- dxu + dxr * 0.001
  ddx <- (dxu - dxl)/(ndx - 1)
  Dknots <- seq(dxl - ddx * (bdeg + 1), dxu + ddx * (bdeg + 1), length = ndx + 2 * bdeg + 2)
  return(Dknots)
}

createD2knots <- function(beta_star, x, d2x, spline_rank, spline_degree, knots){
  
  D2knots <- matrix(nrow=length(knots),ncol=length(beta_star))
  for(d2 in 1:length(beta_star))
  {
    d2xl <- d2x[which.min(x),d2] 
    d2xu <- d2x[which.max(x),d2] 
    ndx <- spline_rank-spline_degree
    bdeg <- spline_degree
    
    d2xr <- d2xu - d2xl #ndx = splinerank - bdeg
    d2xl <- d2xl - d2xr * 0.001
    d2xu <- d2xu + d2xr * 0.001
    d2dx <- (d2xu - d2xl)/(ndx - 1)
    
    D2knots[,d2] <- seq(d2xl - d2dx * (bdeg + 1), d2xu + d2dx * (bdeg + 1), length = ndx + 2 * bdeg + 2)
  }
  
  return(D2knots)
}

random_row_perMute <- function(Design,X,y){
  rownames(Design) <- as.numeric(1:nrow(Design))
  rownames(X) <- as.numeric(1:nrow(X))
  names(y) <- as.numeric(1:length(y))
  
  nc <- ncol(Design)
  
  Design_Perm <- Design
  rownames(Design_Perm) <- rownames(Design)
  
  X_Perm <- X
  rownames(X_Perm) <- rownames(X)
  
  y_Perm <- y
  names(y_Perm) <- names(y)
  
  row_numbers <- vector()
  row_storage <- list()
  for(i in 1:nc)
  {
    if(Design[(nc-i+1),(nc-i+1)]==0) #We work backwards from entry Design[nc,nc]
    {
      non_zero <- Design[Design[,(nc-i+1)]>0, ,drop=FALSE]
      row_numbers[i] <- sample(rownames(non_zero), 1)
      row_storage[[i]] <- non_zero[rownames(non_zero)==row_numbers[i], ,drop=FALSE]
      
      #Design
      Design_Perm[(nc-i+1),] <- non_zero[rownames(non_zero)==row_numbers[i], ,drop=FALSE]
      Design_Perm[row_numbers[i],] <- Design[(nc-i+1),]
      
      
      #X
      X_Perm[(nc-i+1),] <- X[rownames(X)==row_numbers[i], ,drop=FALSE]
      X_Perm[row_numbers[i],] <- X[(nc-i+1),]
      
      #y
      y_Perm[(nc-i+1)] <- y[names(y)==row_numbers[i]]
      y_Perm[row_numbers[i]] <- y[(nc-i+1)]
      
    }
    else 
    {
      row_numbers[i] <- 0
      row_storage[[i]] <- NULL
    }
    
  }
  
  list(Design_Perm=Design_Perm,X_Perm=X_Perm,y_Perm=y_Perm,rownumbers=rev(row_numbers))
}

Apply_perMutations <- function(Object,Input,is.matrix=TRUE) {
  rownumbers <- rev(Object$rownumbers) #rownumbers are passed in reverse for readability so I flip them back
  nc <- length(rownumbers)
  if(is.matrix)
  {
    rownames(Input) <- 1:nrow(Input)
    
    Input_PerMuted <- Input
    rownames(Input_PerMuted) <- rownames(Input)
    
    for(i in 1:nc)
    {
      if(as.numeric(rownumbers[i])>0) #We work backwards from entry Design[nc,nc]
      {
        Input_PerMuted[(nc-i+1),] <- Input[rownames(Input)==as.numeric(rownumbers[i]),]
        Input_PerMuted[as.numeric(rownumbers[i]),] <- Input[(nc-i+1),]      
      }    
    }
  }
  if(!is.matrix)
  {
    names(Input) <- 1:length(Input)
    Input_PerMuted <- Input
    names(Input_PerMuted) <- names(Input)
    for(i in 1:nc)
    {
      if(as.numeric(rownumbers[i])>0) #We work backwards from entry Design[nc,nc]
      {
        Input_PerMuted[(nc-i+1)] <- Input[names(Input)==as.numeric(rownumbers[i])]
        Input_PerMuted[as.numeric(rownumbers[i])] <- Input[(nc-i+1)]
        
      }
    }
  }
  
  Input_PerMuted
}

DChol <- function(P,DP,np) {
  A <- matrix(0,nrow=np,ncol=np)
  D <- matrix(0,nrow=np,ncol=np)
  for(i in 1:np)
  {
    if(i==1)
    {
      A[i,i] <- sqrt( P[i,i] )
      #D[i,i] <- ( DP[i,i] )/(A[i,i]/2) This doesn't seem to work
      D[i,i] <- 0.5 * DP[i,i] * P[i,i]^(-0.5)
      
      for(j in i:np)
      {
        A[j,i] <- P[j,i]/A[i,i]
        temp <- DP[j,i]
        D[j,i] <- temp/A[i,i] - ( D[i,i]/A[i,i] )*A[j,i]
      }
      
    }
    if(i>1)
    {
      k <- seq(1:(i-1))
      A[i,i] <- sqrt( P[i,i]- sum(A[i,k]^2) )
      #D[i,i] <- ( DP[i,i] - 2*sum( D[i,k]*A[i,k] ) )/(A[i,i]/2)
      D[i,i] <- 0.5*( P[i,i] - sum(A[i,k]^2) )^(-0.5) * ( DP[i,i] - 2*sum( A[i,k]*D[i,k] ) )
      
      for(j in i:np)
      {
        A[j,i] <- (P[j,i] - sum( A[j,k]*A[i,k] ))/A[i,i]
        temp <- DP[j,i] - sum( D[j,k]*A[i,k] + A[j,k]*D[i,k] )
        D[j,i] <- temp/A[i,i] - ( D[i,i]/A[i,i] )*A[j,i]
      }
    }          
  }
  list(A=A,D=D)
}

#Fixed some bugs
#P is the matrix we're decomposing, DP is the partial derivative, and P lives in the space of np x np matrices
D2Chol <- function(P,DP_list,D2P,np,beta_star,deriv_1,deriv_2) {
  #DP_list needs to be the derivatives of P with respect to beta_star in corresponding order
  A <- matrix(0,nrow=np,ncol=np)
  D <- list()
  for(i in 1:length(beta_star))
  {
    D[[i]] <- matrix(0,nrow=np,ncol=np)
  }
  
  DP <- DP_list
  D2 <- matrix(0,nrow=np,ncol=np)  
  
  for(i in 1:np)
  {
    if(i==1)
    {
      D2[i,i] <- -0.25*P[i,i]^(-1.5)*DP[[deriv_1]][i,i]*DP[[deriv_2]][i,i] + 0.5*D2P[i,i]*P[i,i]^(-0.5)
      
      D[[deriv_1]][i,i] <- 0.5 * DP[[deriv_1]][i,i] * P[i,i]^(-0.5)  #These are fine
      D[[deriv_2]][i,i] <- 0.5 * DP[[deriv_2]][i,i] * P[i,i]^(-0.5)
      
      A[i,i] <- sqrt( P[i,i] )
      
      for(j in i:np)
      {
        A[j,i] <- P[j,i]/A[i,i]
        
        #Update first derivatives
        # D[[deriv_1]][j,i] <- DP[[deriv_1]][j,i]/A[i,i] - ( D[[deriv_1]][i,i]/A[i,i] )*A[j,i]
        
        #Slight bug in this previously:
        D[[deriv_1]][j,i] <- DP[[deriv_1]][j,i]/A[i,i] - D[[deriv_1]][i,i]*P[j,i]/(A[i,i]^2)
        
        if(deriv_1 != deriv_2)
        {
          # D[[deriv_2]][j,i] <- DP[[deriv_2]][j,i]/A[i,i] - ( D[[deriv_2]][i,i]/A[i,i] )*A[j,i]
          
          #Slight bug in this previously
          D[[deriv_2]][j,i] <- DP[[deriv_2]][j,i]/A[i,i] - ( D[[deriv_2]][i,i]*P[j,i]/(A[i,i]^2) )
        }
        
        #Update second derivatives
        # D2[j,i] <- D2P[j,i]/A[i,i] - ( DP[[deriv_1]][j,i]*D[[deriv_2]][i,i] )/(A[i,i]^2) - 
        #   
        #   ( D[[deriv_2]][j,i]*D[[deriv_1]][i,i] + A[j,i]*D2[i,i] )/(A[i,i]) + ( A[j,i]*D[[deriv_1]][i,i]*D[[deriv_2]][i,i] )/(A[i,i]^2)
        
        #Slight bug in this previously:
        D2[j,i] <- D2P[j,i]/A[i,i] - ( DP[[deriv_1]][j,i]*D[[deriv_2]][i,i] )/(A[i,i]^2) - 
          ( DP[[deriv_2]][j,i]*D[[deriv_1]][i,i] + P[j,i]*D2[i,i] )/(A[i,i]^2) + ( P[j,i]*D[[deriv_1]][i,i]*2*A[i,i]*D[[deriv_2]][i,i] )/(A[i,i]^4)
      }
      
    }
    if(i>1)
    {
      k <- seq(1:(i-1))
      
      D2[i,i] <- ( 
        (D2P[i,i] - 2*sum(D[[deriv_1]][i,k]*D[[deriv_2]][i,k] + A[i,k]*D2[i,k])) * 0.5 *(( P[i,i]-sum(A[i,k]^2) )^(-0.5))
        -0.25*(DP[[deriv_1]][i,i]-2*sum(A[i,k]*D[[deriv_1]][i,k]))*(DP[[deriv_2]][i,i]-2*sum(A[i,k]*D[[deriv_2]][i,k]))*(( P[i,i]-sum(A[i,k]^2) )^(-1.5))
      )
      
      D[[deriv_1]][i,i] <- 0.5*((P[i,i]-sum(A[i,k]^2))^(-0.5))*(DP[[deriv_1]][i,i]-2*sum(A[i,k]*D[[deriv_1]][i,k]))
      D[[deriv_2]][i,i] <- 0.5*((P[i,i]-sum(A[i,k]^2))^(-0.5))*(DP[[deriv_2]][i,i]-2*sum(A[i,k]*D[[deriv_2]][i,k]))
      
      A[i,i] <- sqrt( P[i,i]- sum(A[i,k]^2) )
      for(j in i:np)
      {
        A[j,i] <- (P[j,i] - sum( A[j,k]*A[i,k] ))/A[i,i]
        
        #Update first derivatives
        temp <- DP[[deriv_1]][j,i] - sum( D[[deriv_1]][j,k]*A[i,k] + A[j,k]*D[[deriv_1]][i,k] )
        # D[[deriv_1]][j,i] <- temp/A[i,i] - ( D[[deriv_1]][i,i]/A[i,i] )*A[j,i]
        
        #Minor bug
        D[[deriv_1]][j,i] <- temp/A[i,i] - ( D[[deriv_1]][i,i]*(P[j,i] - sum(A[j,k]*A[i,k]))/(A[i,i]^2) )
        
        if(deriv_1 != deriv_2)
        {
          temp2 <- DP[[deriv_2]][j,i] - sum( D[[deriv_2]][j,k]*A[i,k] + A[j,k]*D[[deriv_2]][i,k] )
          # D[[deriv_2]][j,i] <- temp2/A[i,i] - ( D[[deriv_2]][i,i]/A[i,i] )*A[j,i]
          
          #Minor bug
          D[[deriv_2]][j,i] <- temp2/A[i,i] - ( D[[deriv_2]][i,i]*(P[j,i] - sum(A[j,k]*A[i,k]))/(A[i,i]^2) )
        }
        
        ##Update second derivatives
        D2[j,i] <- (
          ( D2P[j,i] - sum( D2[j,k]*A[i,k] + D[[deriv_1]][j,k]*D[[deriv_2]][i,k] + D[[deriv_2]][j,k]*D[[deriv_1]][i,k] + A[j,k]*D2[i,k] ) )/A[i,i] -
            ( (DP[[deriv_1]][j,i] - sum( D[[deriv_1]][j,k]*A[i,k] + A[j,k]*D[[deriv_1]][i,k] ) )*D[[deriv_2]][i,i] )/(A[i,i]^2) -
            ( (DP[[deriv_2]][j,i] - sum( D[[deriv_2]][j,k]*A[i,k] + A[j,k]*D[[deriv_2]][i,k] ) )*D[[deriv_1]][i,i] + ( P[j,i] - sum(A[j,k]*A[i,k]) )*D2[i,i] )/(A[i,i]^2) +
            ( ( (P[j,i] - sum(A[j,k]*A[i,k])) * D[[deriv_1]][i,i] ) * 2 * A[i,i] * D[[deriv_2]][i,i] )/(A[i,i]^4)  
        )
        
      }
    }          
  }
  list(A=A,D=D[[deriv_1]],D2=D2)
}

#Updated version 
QR.Hessian <- function(Deriv2_Design,Deriv_Design_list,Design,Deriv2_B,Deriv_B_list, B, deriv_1,deriv_2) { #Deriv_Design_list should be the list containing all of the derivatives of Design with respect to beta_star
  #deriv_1, deriv_2 corresponds to the hessian we're calculating, where we first take derivatives wrt beta^star_deriv_1 and then beta^star_deriv_2
  #Similar for Deriv_B_list
  nr <- nrow(Design)
  nc <- ncol(Design)
  q <- ncol(B)
  p <- nrow(B)
  #The first iteration of the householder matrices, H, is H^0 = Deriv_Design, namely the matrix of element wise derivatives of the spline 
  #transformed design matrix, dependent on the index parameter vector beta
  
  #At each iteration of the recursive differentation, we need the householder matrix and the matrix of derivatives from the iteration before.
  #Thus, to initialize, we need to pass this function the design matrix and the matrix of derivatives. In my case, I need the spline transformation 
  #matrix as well as the matrix of the derivatives of the spline transformation
  
  #H^0 = Design
  #D(H^0) = Deriv_Design
  
  #Deriv_B will be the matrix of zeros when B is the identity - B being the identity gives us Q as the result
  
  H <- Design
  D_H <- Deriv_Design_list
  D2_H <- Deriv2_Design
  
  Deriv_B <- Deriv_B_list
  
  # Householder transformations
  for (j in seq_len(nc))
  {
    H_copy <- H
    t <- seq.int(j, nr)
    t2 <- seq.int(j,p)
    sigma <- sum(H[t,j]^2)
    
    D_sigma <- list()
    D_sigma[[deriv_1]] <- 2*sum(D_H[[deriv_1]][t,j]*H[t,j])
    D_sigma[[deriv_2]] <- 2*sum(D_H[[deriv_2]][t,j]*H[t,j])
    
    D2_sigma <- 2*sum(D2_H[t,j]*H[t,j] + D_H[[deriv_1]][t,j]*D_H[[deriv_2]][t,j] )    
    
    s <- sqrt(sigma)
    D_s <- list()
    D_s[[deriv_1]] <- D_sigma[[deriv_1]] * 0.5 * sigma^(-0.5)
    D_s[[deriv_2]] <- D_sigma[[deriv_2]] * 0.5 * sigma^(-0.5)
    
    D2_s <- 0.5 * ( D2_sigma*sigma^(-0.5) - ( 0.5*(D_sigma[[deriv_1]]*D_sigma[[deriv_2]] )*sigma^(-1.5) ) )
    
    
    diag_ej <- H[j,j]
    
    D_diag_ej <- list()
    D_diag_ej[[deriv_1]] <- D_H[[deriv_1]][j,j]
    D_diag_ej[[deriv_2]] <- D_H[[deriv_2]][j,j]
    
    D2_diag_ej <- D2_H[j,j]
    
    
    gamma <- 1.0 / (sigma + abs(s * diag_ej))
    ##We differentiate the denmominator of gamma which we call D_gamma_denom
    D_gamma_denom <- list()
    D_gamma_denom[[deriv_1]] <- D_sigma[[deriv_1]] + (( (D_s[[deriv_1]] * diag_ej) + (s * D_diag_ej[[deriv_1]]) )*(s * diag_ej))/abs(s * diag_ej)
    D_gamma_denom[[deriv_2]] <- D_sigma[[deriv_2]] + (( (D_s[[deriv_2]] * diag_ej) + (s * D_diag_ej[[deriv_2]]) )*(s * diag_ej))/abs(s * diag_ej)
    
    
    D2_gamma_denom <-
      (
        D2_sigma + 
          
          ( ( (D2_s * diag_ej) + (D_s[[deriv_1]]*D_diag_ej[[deriv_2]]) + (D_s[[deriv_2]]*D_diag_ej[[deriv_1]]) + (s * D2_diag_ej) )*(s * diag_ej) + ((D_s[[deriv_1]] * diag_ej) + (s * D_diag_ej[[deriv_1]]))*((D_s[[deriv_2]] * diag_ej) + (s * D_diag_ej[[deriv_2]])))/abs(s * diag_ej) 
        
        - ( ( ( (D_s[[deriv_1]] * diag_ej) + (s * D_diag_ej[[deriv_1]]) )*( (D_s[[deriv_2]] * diag_ej) + (s * D_diag_ej[[deriv_2]]) )*((s * diag_ej)^2) )/((abs(s * diag_ej))^3) )
      )
    
    kappa <- if (diag_ej < 0) s else -s
    
    D_kappa <- list()
    D_kappa[[deriv_1]] <- if (diag_ej < 0) D_s[[deriv_1]] else -D_s[[deriv_1]]
    D_kappa[[deriv_2]] <- if (diag_ej < 0) D_s[[deriv_2]] else -D_s[[deriv_2]]
    
    
    D2_kappa <- if (diag_ej < 0) D2_s else -D2_s
    
    g <- 1/gamma  
    
    D_g <- list()
    D_g[[deriv_1]] <- D_sigma[[deriv_1]] + ((D_diag_ej[[deriv_1]] *s + diag_ej * D_s[[deriv_1]]) * diag_ej * s)/abs(s * diag_ej)
    D_g[[deriv_2]] <- D_sigma[[deriv_2]] + ((D_diag_ej[[deriv_2]] *s + diag_ej * D_s[[deriv_2]]) * diag_ej * s)/abs(s * diag_ej)
    
    q <- seq.int(j+1, nr)
    
    H[j,j] <- H[j,j] - kappa
    if (j < nc)
      for (k in seq.int(j+1, nc))
      {
        yPrime <- sum(H[t,j] * H[t,k]) * gamma
        H[t,k] <- H[t,k] - H[t,j] * yPrime             
      }
    
    #######We add this additional loop which shall return Q      
    #Note that at this stage, H has been updated to the extent that the H[j,j]=H[j,j]-kappa (so not their final forMulation) however the rest
    #of H has the final forMulation
    
    #Note that at this stage, D_H will need to have been updated
    
    t3 <- seq.int(j+1,p)   
    for (k2 in 1:ncol(B))
    {
      yPrime <- sum(H[t2,j] * B[t2,k2]) * gamma
      
      f_yPrime <- sum(H[t2,j] * B[t2,k2])
      
      D_f_yPrime <- list()
      D_f_yPrime[[deriv_1]] <- sum( (D_H[[deriv_1]][t3,j] * B[t3,k2]) + (H[t3,j] * Deriv_B[[deriv_1]][t3,k2]) ) + ( (D_diag_ej[[deriv_1]] - D_kappa[[deriv_1]]) * B[j,k2] ) + ( (diag_ej-kappa) * Deriv_B[[deriv_1]][j,k2] )
      if(deriv_1!=deriv_2)
      {
        D_f_yPrime[[deriv_2]] <- sum( (D_H[[deriv_2]][t3,j] * B[t3,k2]) + (H[t3,j] * Deriv_B[[deriv_2]][t3,k2]) ) + ( (D_diag_ej[[deriv_2]] - D_kappa[[deriv_2]]) * B[j,k2] ) + ( (diag_ej-kappa) * Deriv_B[[deriv_2]][j,k2] )
      }
      D2_f_yPrime <- sum( (D2_H[t3,j] * B[t3,k2]) + (D_H[[deriv_1]][t3,j]*Deriv_B[[deriv_2]][t3,k2]) + (D_H[[deriv_2]][t3,j]*Deriv_B[[deriv_1]][t3,k2]) + (H[t3,j]*Deriv2_B[t3,k2]) ) + (( D2_diag_ej-D2_kappa )*B[j,k2] ) + ( (D_diag_ej[[deriv_1]]-D_kappa[[deriv_1]])*Deriv_B[[deriv_2]][j,k2] ) + ( (D_diag_ej[[deriv_2]]-D_kappa[[deriv_2]])*Deriv_B[[deriv_1]][j,k2] ) + ( (diag_ej-kappa)*Deriv2_B[j,k2] )
      
      D_yPrime <- list()
      D_yPrime[[deriv_1]] <- D_f_yPrime[[deriv_1]]/(1/gamma) - (f_yPrime * D_gamma_denom[[deriv_1]])/((1/gamma)^2)
      if(deriv_1!=deriv_2)
      {
        D_yPrime[[deriv_2]] <- D_f_yPrime[[deriv_2]]/(1/gamma) - (f_yPrime * D_gamma_denom[[deriv_2]])/((1/gamma)^2)
      }
      D2_yPrime <- D2_f_yPrime/(1/gamma) - (D_f_yPrime[[deriv_1]] * D_gamma_denom[[deriv_2]])/((1/gamma)^2) - (D_f_yPrime[[deriv_2]]*D_gamma_denom[[deriv_1]] + f_yPrime*D2_gamma_denom)/((1/gamma)^2) + ( 2*f_yPrime*(D_gamma_denom[[deriv_1]]*D_gamma_denom[[deriv_2]])*(1/gamma) )/((1/gamma)^4)
      
      ##Update the second derivatives first
      #t3=j
      Deriv2_B[j,k2] <- Deriv2_B[j,k2] - (D2_diag_ej-D2_kappa)*yPrime - (D_diag_ej[[deriv_1]]-D_kappa[[deriv_1]])*D_yPrime[[deriv_2]] - (D_diag_ej[[deriv_2]]-D_kappa[[deriv_2]])*D_yPrime[[deriv_1]]- (diag_ej-kappa)*D2_yPrime
      #t3 =/= j
      Deriv2_B[t3,k2] <- Deriv2_B[t3,k2] - D2_H[t3,j]*yPrime - D_H[[deriv_1]][t3,j]*D_yPrime[[deriv_2]] - D_H[[deriv_2]][t3,j]*D_yPrime[[deriv_1]] - H[t3,j]*D2_yPrime
      
      #Update the first derivatives
      #t3=j
      Deriv_B[[deriv_1]][j,k2] <- Deriv_B[[deriv_1]][j,k2] - (D_diag_ej[[deriv_1]]-D_kappa[[deriv_1]])*yPrime - (diag_ej-kappa)*D_yPrime[[deriv_1]]
      #t3 =/= j
      Deriv_B[[deriv_1]][t3,k2] <- Deriv_B[[deriv_1]][t3,k2]- D_H[[deriv_1]][t3,j]*yPrime - H[t3,j]*D_yPrime[[deriv_1]]
      
      if(deriv_1!=deriv_2)
      {
        #t3=j
        Deriv_B[[deriv_2]][j,k2] <- Deriv_B[[deriv_2]][j,k2] - (D_diag_ej[[deriv_2]]-D_kappa[[deriv_2]])*yPrime - (diag_ej-kappa)*D_yPrime[[deriv_2]]
        #t3 =/= j
        Deriv_B[[deriv_2]][t3,k2] <- Deriv_B[[deriv_2]][t3,k2]- D_H[[deriv_2]][t3,j]*yPrime - H[t3,j]*D_yPrime[[deriv_2]]
      }
      
      
      #Update B
      B[t2,k2] <- B[t2,k2] - H[t2,j] * yPrime
      
    }
    
    
    #This Must be done after computing Q (B)
    H[j,j] <- kappa
    
    #This gives me the householder matrix for the current iteration, I just need to calculate the householder matrix derivatives for this iteration
    #We retain the same value of j
    #We shall calculate the components seperately, particularly for the two instances of the quotient rule we need to calculate
    #Need to remember that in my written derivatives I differentiated d/dbeta^star{H^(j+1)}, where as here I'm using H^(j), so need to replace 
    #j+1's with j's
    
    #j: 1 to nc
    #t: j to nc
    #k: j+1 to nc
    #p: j to nc <=> p: t to nc
    #q: j to nc <=> q: t to nc
    
    
    
    if (j < nc)
      for (k in seq.int(j+1, nc))
      { 
        #use f_1 when t=j
        f1 <- (diag_ej-kappa) * ( sum(H_copy[q,j]*H_copy[q,k]) + (diag_ej-kappa)*H_copy[j,k] )
        D_f1 <- (D_diag_ej[[deriv_1]]-D_kappa[[deriv_1]]) * ( sum(H_copy[q,j]*H_copy[q,k]) + (diag_ej-kappa)*H_copy[j,k] )   + (diag_ej-kappa) * (  sum( D_H[[deriv_1]][q,j]*H_copy[q,k] + H_copy[q,j]*D_H[[deriv_1]][q,k] ) + (D_diag_ej[[deriv_1]]-D_kappa[[deriv_1]]) * H_copy[j,k]  + (diag_ej-kappa)*D_H[[deriv_1]][j,k] )       
        #use f_2 when t>j
        f2 <- H_copy[q,j]* ( sum(H_copy[q,j]*H_copy[q,k]) + (diag_ej-kappa)*H_copy[j,k] )
        D_f2 <- D_H[[deriv_1]][q,j] * ( sum(H_copy[q,j]*H_copy[q,k]) + (diag_ej-kappa)*H_copy[j,k] ) + H_copy[q,j] * (  sum( D_H[[deriv_1]][q,j]*H_copy[q,k] + H_copy[q,j]*D_H[[deriv_1]][q,k] ) + (D_diag_ej[[deriv_1]]-D_kappa[[deriv_1]]) * H_copy[j,k]  + (diag_ej-kappa)*D_H[[deriv_1]][j,k] ) 
        
        #when t=j
        D_H[[deriv_1]][j,k] <- D_H[[deriv_1]][j,k] - ( D_f1/g - (f1*D_g[[deriv_1]])/(g^2) )
        #when t>j
        D_H[[deriv_1]][q,k] <- D_H[[deriv_1]][q,k] - ( D_f2/g - (f2*D_g[[deriv_1]])/(g^2) )
        
      }
    if(deriv_1!=deriv_2)
    {
      if (j < nc)
        for (k in seq.int(j+1, nc))
        { 
          #use f_1 when t=j
          f1 <- (diag_ej-kappa) * ( sum(H_copy[q,j]*H_copy[q,k]) + (diag_ej-kappa)*H_copy[j,k] )
          D_f1 <- (D_diag_ej[[deriv_2]]-D_kappa[[deriv_2]]) * ( sum(H_copy[q,j]*H_copy[q,k]) + (diag_ej-kappa)*H_copy[j,k] )   + (diag_ej-kappa) * (  sum( D_H[[deriv_2]][q,j]*H_copy[q,k] + H_copy[q,j]*D_H[[deriv_2]][q,k] ) + (D_diag_ej[[deriv_2]]-D_kappa[[deriv_2]]) * H_copy[j,k]  + (diag_ej-kappa)*D_H[[deriv_2]][j,k] )       
          #use f_2 when t>j
          f2 <- H_copy[q,j]* ( sum(H_copy[q,j]*H_copy[q,k]) + (diag_ej-kappa)*H_copy[j,k] )
          D_f2 <- D_H[[deriv_2]][q,j] * ( sum(H_copy[q,j]*H_copy[q,k]) + (diag_ej-kappa)*H_copy[j,k] ) + H_copy[q,j] * (  sum( D_H[[deriv_2]][q,j]*H_copy[q,k] + H_copy[q,j]*D_H[[deriv_2]][q,k] ) + (D_diag_ej[[deriv_2]]-D_kappa[[deriv_2]]) * H_copy[j,k]  + (diag_ej-kappa)*D_H[[deriv_2]][j,k] ) 
          
          #when t=j
          D_H[[deriv_2]][j,k] <- D_H[[deriv_2]][j,k] - ( D_f1/g - (f1*D_g[[deriv_2]])/(g^2) )
          #when t>j
          D_H[[deriv_2]][q,k] <- D_H[[deriv_2]][q,k] - ( D_f2/g - (f2*D_g[[deriv_2]])/(g^2) )
          
        }
      D_H[[deriv_2]][j,j] <- D_kappa[[deriv_2]]
    }
    
    D_H[[deriv_1]][j,j] <- D_kappa[[deriv_1]]
    
  } # end Householder
  
  # set zeros in the lower triangular side of X (which stores) 
  # not really necessary, this is just to return R for illustration
  for (i in seq_len(ncol(H)))
    H[seq.int(i+1, nr),i] <- 0
  
  for (i in seq_len(ncol(D_H[[deriv_1]])))
    D_H[[deriv_1]][seq.int(i+1, nr),i] <- 0 
  
  if(deriv_1!=deriv_2)
  {
    for (i in seq_len(ncol(D_H[[deriv_1]])))
      D_H[[deriv_2]][seq.int(i+1, nr),i] <- 0 
  }
  
  
  list(R=H[1:nc,1:nc], D_R=D_H[[deriv_1]][1:nc,1:nc],B=B, D_B= Deriv_B[[deriv_1]], D2_B=Deriv2_B)  
}

######GAM functions:
gam.reparamR.gam <- function (rS, lsp, deriv) 
{
  q <- nrow(rS[[1]])
  rSncol <- unlist(lapply(rS, ncol))
  M <- length(lsp)
  if (length(rS) > M){
    fixed.penalty <- TRUE
  } else fixed.penalty <- FALSE
  
  d.tol <- .Machine$double.eps^0.3
  r.tol <- .Machine$double.eps^0.75
  
  #Form Si:
  Si <- list()
  for(i in 1:length(lsp)){
    Si[[i]] <- rS[[i]] %*% t(rS[[i]])
  }
  
  oo <- get_stableS.gam(rS = rS, deriv = deriv, dTol = d.tol, rTol = r.tol, lsp = lsp, Si = Si)
  
  S <- oo$S
  S <- (S + t(S)) * 0.5
  p <- abs(diag(S))^0.5
  p[p == 0] <- 1
  St <- t(t(S/p)/p)
  St <- (St + t(St)) * 0.5
  E <- t(mroot(St, rank = q)$L * p)
  Qs <- oo$Qs
  #Just re-assign rS to the outputs:
  rS <- oo$rS
  
  if (deriv > 0){
    det1 <- oo$det1
  } else{ 
    det1 <- NULL
  }
  if (deriv > 1){
    det2 <- oo$det2
  } else{ 
    det2 <- NULL
  }
  list(S = S, E = E, Qs = Qs, rS = rS, det = oo$det, det1 = det1, 
       det2 = det2, fixed.penalty = fixed.penalty)
}

get_stableS.gam <- function(beta_star, deriv, dTol = .Machine$double.eps^0.3, 
                            rTol = .Machine$double.eps^0.75, lsp, Si, rS){
  
  Si_Q <- Si
  rank <-p <- q <- ncol(Si[[1]])
  sp <- exp(lsp)
  iter <- 0
  
  drho_det <- NULL
  d2rho_det <- NULL
  
  #If we don't want derivatives we do the simple version of the routine which just returns Qs
  #Initialize parameters:
  Qs <- diag(rank)
  m <- length(Si)
  K <- 0
  Q <- rank
  gamma <- rep(1,times=m)
  
  #Form the total penalty matrix
  S <- matrix(0,nrow=p,ncol=q)
  for(k in 1:length(Si))
  {
    S <- S + sp[k]*Si[[k]]
  }
  
  while(1)
  {
    Omega <- rep(0,times = length(gamma))  
    # Step 1
    for(k in 1:length(gamma))
    {
      if(gamma[k] == 1)
      {
        Omega[k] <- sqrt(sum(Si[[k]] * Si[[k]])) * sp[k]
      }
    }
    
    #Step 2
    alpha <- as.numeric( Omega >= dTol * max(Omega) & gamma == 1) 
    gamma_prime <- as.numeric( Omega < dTol * max(Omega) & gamma == 1)
    
    S_Scaled <- matrix(0,nrow=Q,ncol=Q)
    #Step 3
    for(k in 1:length(alpha))
    {
      if(alpha[k] == 1)
      {
        S_Scaled <- S_Scaled + Si[[k]]/sqrt(sum(Si[[k]] * Si[[k]]))
      }
    }
    
    es <- eigen(S_Scaled, symmetric = TRUE)
    ind <- es$values > max(es$values) * rTol
    penaltyRank <- sum(ind)
    
    #Step 4
    if(Q == penaltyRank){
      break
    }
    
    # Step 5
    S_Unscaled <- matrix(0,nrow=Q,ncol=Q)
    for(k in 1:m)
    {
      if(alpha[k] == 1)
      {
        S_Unscaled <- S_Unscaled + sp[k]*Si[[k]]
      }
    }
    
    es0 <- eigen(S_Unscaled, symmetric = TRUE)
    U0 <- es0$vectors
    evals <- es0$values
    Ur <- es0$vectors[, 1:penaltyRank]
    Un <- es0$vectors[, -c(1:penaltyRank)]
    
    U <- matrix(0,nrow = (K+Q), ncol = (K+Q))
    Id <- diag(nrow=p,ncol=q)
    
    
    if(K==0)
    {
      U <- U0
    } else if (K > 0) {
      U[1:K,1:K] <- diag(K)
      U[(K+1):(K+Q),(K+1):(K+Q)] <- U0
    }
    
    StPrime <-  crossprod(U, S) %*% U
    
    for(k in 1:length(gamma_prime))
    {
      if(gamma_prime[k] == 1){
        Si[[k]] <- crossprod(Un,Si[[k]]) %*% Un  
      }
    } 
    
    Qs <- U %*% Qs
    K <- K + penaltyRank
    Q <- Q - penaltyRank
    S <- StPrime
    gamma <- gamma_prime
  }
  
  #Transform the rS
  for(i in 1:m){
    rS[[i]] <- t(Qs) %*% rS[[i]]
  }
  
  #Transform the Si
  Si_Q <- list()
  for(i in 1:m){
    Si_Q[[i]] <- rS[[i]] %*% t(rS[[i]])
  }
  
  #Calculate log of determinant
  det <- log(det(S))
  
  #Deriv positive then we calculate derivatives of log determinant
  if(deriv > 0){
    #We calculate derivatives of log|S|
    drho_det <- vector()
    d2rho_det <- matrix(0,length(sp),length(sp))
    
    inv_S <- solve(S)
    
    for(i in 1:length(sp)){
      for(j in 1:length(sp)){
        
        k_delta <- if(i==j) 1 else 0
        
        d2rho_det[i,j] <- k_delta * sp[i] * sum(diag(inv_S %*% Si_Q[[i]])) - 
          sp[i] * sp[j] * sum(diag(inv_S %*% Si_Q[[i]] %*% inv_S %*% Si_Q[[j]]))
      }
      drho_det[i] <- sp[i] * sum(diag(inv_S %*% Si_Q[[i]]))
    }
  }
  
  return(list(Si = Si_Q, Qs = Qs, S = S, det = det, det1 = drho_det, det2 = d2rho_det, rS = rS))
}


####Single index versions
gam.reparamR.si <- function(rS, d_rS, d2_rS, lsp, deriv, beta_star = NULL, d_lsp = NULL, d2_lsp = NULL, 
                            Si = NULL, d_Si = NULL, d2_Si = NULL, dpen){
  q <- nrow(rS[[1]])
  rSncol <- unlist(lapply(rS, ncol))
  M <- length(lsp)
  if (length(rS) > M){
    fixed.penalty <- TRUE
  } else fixed.penalty <- FALSE
  
  d.tol <- .Machine$double.eps^0.3
  r.tol <- .Machine$double.eps^0.75
  
  drho_S <- d2rho_S <- drho_d_S <- drho_d2_S <- d2rho_d_S <- d2rho_d2_S <- list()
  for(d1 in 1:length(beta_star))
  {
    drho_d_S[[d1]] <- drho_d2_S[[d1]] <- d2rho_d_S[[d1]] <- d2rho_d2_S[[d1]] <- list()
    for(d2 in d1:length(beta_star))
    {
      drho_d2_S[[d1]][[d2]] <- list()
      d2rho_d2_S[[d1]][[d2]] <- list()
      
      for(i in 1:length(lsp))
      {
        drho_S[[i]] <- exp(lsp[i])*dpen$Penalty[[i]]
        
        drho_d_S[[d1]][[i]] <- d_lsp[[i]][d1] * dpen$Penalty[[i]] + exp(lsp[i]) * dpen$d_Penalty[[i]][[d1]]
        
        drho_d2_S[[d1]][[d2]][[i]] <- (
          d2_lsp[[i]][d1,d2] * dpen$Penalty[[i]] + d_lsp[[i]][d1] * dpen$d_Penalty[[i]][[d2]] + 
            d_lsp[[i]][d2] * dpen$d_Penalty[[i]][[d1]] + exp(lsp[i]) * dpen$d2_Penalty[[i]][[d1]][[d2]]
        )
        
        d2rho_S[[i]] <- d2rho_d_S[[d1]][[i]] <- d2rho_d2_S[[d1]][[d2]][[i]] <- list()
        for(j in i:length(lsp))
        {
          d2rho_S[[i]][[j]] <- if(i==j)
          {
            exp(lsp[i])*dpen$Penalty[[i]]
          }
          else
          {
            0*dpen$Penalty[[i]]
          }
          
          
          d2rho_d_S[[d1]][[i]][[j]] <- if(i==j)
          {
            d_lsp[[i]][d1] * dpen$Penalty[[i]] + exp(lsp[i]) * dpen$d_Penalty[[i]][[d1]]
          }
          else
          {
            0*dpen$Penalty[[i]]
          }
          
          
          d2rho_d2_S[[d1]][[d2]][[i]][[j]] <- if(i==j)
          {
            (
              d2_lsp[[i]][d1,d2] * dpen$Penalty[[i]] + d_lsp[[i]][d1] * dpen$d_Penalty[[i]][[d2]] + 
                d_lsp[[i]][d2] * dpen$d_Penalty[[i]][[d1]] + exp(lsp[i]) * dpen$d2_Penalty[[i]][[d1]][[d2]]
            )
          }
          else
          {
            0*dpen$Penalty[[i]]
          }
          
        } 
      }
    }
  }
  
  #Inputs:
  # beta_star = beta_star
  # deriv = deriv
  # dTol = d.tol
  # rTol = r.tol
  # rS = rS
  # d_rS = d_rS
  # d2_rS = d2_rS
  # lsp = lsp
  # d_lsp = d_lsp
  # d2_lsp = d2_lsp
  # Si = Si
  # d_Si = d_Si
  # d2_Si = d2_Si
  # drho_S = drho_S
  # d2rho_S = d2rho_S
  # drho_d_S = drho_d_S
  # drho_d2_S = drho_d2_S
  # d2rho_d_S = d2rho_d_S
  # d2rho_d2_S = d2rho_d2_S
  
  oo <- get_stableS.si(beta_star = beta_star, deriv = deriv, dTol = d.tol, rTol = r.tol, rS = rS, d_rS = d_rS, d2_rS = d2_rS,
                       lsp = lsp, d_lsp = d_lsp, d2_lsp = d2_lsp, Si = Si, d_Si = d_Si,
                       d2_Si = d2_Si, drho_S = drho_S, d2rho_S = d2rho_S, drho_d_S = drho_d_S,
                       drho_d2_S = drho_d2_S, d2rho_d_S = d2rho_d_S, d2rho_d2_S = d2rho_d2_S)
  S <- oo$S
  S <- (S + t(S)) * 0.5
  
  d_Si <- oo$d_Si
  d2_Si <- oo$d2_Si
  d_S <- oo$d_S
  d2_S <- oo$d2_S
  drho_S <- oo$drho_S
  d2rho_S <- oo$d2rho_S
  drho_d_S <- oo$drho_d_S
  drho_d2_S <- oo$drho_d2_S
  d2rho_d_S <- oo$d2rho_d_S
  d2rho_d2_S <- oo$d2rho_d2_S
  for(i in 1:length(lsp)){
    for(j in i:length(lsp)){
      for(d1 in 1:length(beta_star)){
        for(d2 in d1:length(beta_star)){
          d2rho_d2_S[[d1]][[d2]][[i]][[j]] <- (d2rho_d2_S[[d1]][[d2]][[i]][[j]] + t(d2rho_d2_S[[d1]][[d2]][[i]][[j]])) * 0.5
        }
        d2rho_d_S[[d1]][[i]][[j]] <- (d2rho_d_S[[d1]][[i]][[j]] + t(d2rho_d_S[[d1]][[i]][[j]])) * 0.5
      }
      d2rho_S[[i]][[j]] <- (d2rho_S[[i]][[j]] + t(d2rho_S[[i]][[j]])) * 0.5
    }
    drho_S[[i]] <- (drho_S[[i]] + t(drho_S[[i]])) * 0.5
    
    for(d1 in 1:length(beta_star)){
      for(d2 in d1:length(beta_star)){
        drho_d2_S[[d1]][[d2]][[i]] <- (drho_d2_S[[d1]][[d2]][[i]] + t(drho_d2_S[[d1]][[d2]][[i]])) * 0.5
      }
      drho_d_S[[d1]][[i]] <- (drho_d_S[[d1]][[i]] + t(drho_d_S[[d1]][[i]])) * 0.5
    }
  }
  
  for(d1 in 1:length(beta_star)){
    for(d2 in d1:length(beta_star)){
      d2_S[[d1]][[d2]] <- (d2_S[[d1]][[d2]] + t(d2_S[[d1]][[d2]])) * 0.5
    }
    d_S[[d1]] <- (d_S[[d1]] + t(d_S[[d1]])) * 0.5
  }
  
  p <- abs(diag(S))^0.5
  p[p == 0] <- 1
  St <- t(t(S/p)/p)
  St <- (St + t(St)) * 0.5
  E <- t(mroot(St, rank = q)$L * p)
  Qs <- oo$Qs
  d_Qs <- oo$d_Qs
  d2_Qs <- oo$d2_Qs
  rS <- oo$rS
  d_rS <- oo$d_rS
  d2_rS <- oo$d2_rS
  
  if(deriv > 0){
    det1 <- oo$det1
    d_det1 <- oo$d_det1
    d2_det1 <- oo$d2_det1
  } else{ 
    det1 <- NULL
    d_det1 <- 0
    d2_det1 <- 0
  }
  
  if(deriv > 1){
    det2 <- oo$det2
    d_det2 <- oo$d_det2
    d2_det2 <- oo$d2_det2
  } else{ 
    det2 <- NULL
    d_det2 <- 0
    d2_det2 <- 0
  }
  
  list(S = S, E = E, Qs = Qs, d_Qs = d_Qs, d2_Qs = d2_Qs, rS = rS, d_rS = d_rS, d2_rS = d2_rS, d_Si = d_Si, d2_Si = d2_Si,
       d_S = d_S, d2_S = d2_S, drho_S = drho_S, d2rho_S = d2rho_S, drho_d_S = drho_d_S,
       drho_d2_S = drho_d2_S, d2rho_d_S = d2rho_d_S, d2rho_d2_S = d2rho_d2_S,
       det = oo$det, d_det = oo$d_det, d2_det = oo$d2_det,
       det1 = det1, d_det1 = d_det1, d2_det1 = d2_det1,
       det2 = det2, d_det2 = d_det2, d2_det2 = d2_det2,
       fixed.penalty = fixed.penalty)
}

get_stableS.si <- function(beta_star, deriv, dTol = .Machine$double.eps^0.3, 
                           rTol = .Machine$double.eps^0.75, rS, d_rS, d2_rS,
                           lsp, d_lsp, d2_lsp, Si, d_Si, d2_Si, drho_S, 
                           d2rho_S, drho_d_S, drho_d2_S, d2rho_d_S, d2rho_d2_S){
  
  Si_Q <- Si
  d_Si_Q <- d_Si
  d2_Si_Q <- d2_Si
  rank <-p <- q <- ncol(Si[[1]])
  sp <- exp(lsp)
  iter <- 0
  
  drho_det <- NULL
  d2rho_det <- NULL
  det <- d_det <- d2_det <- drho_det <- drho_d_det <- drho_d2_det <- d2rho_det <- d2rho_d_det <- d2rho_d2_det <- NULL
  
  #Initialize parameters:
  Qs <- diag(rank)
  m <- length(Si)
  K <- 0
  Q <- rank
  gamma <- rep(1,times=m)
  
  d_Qs <- d2_Qs <- list()
  for(d1 in 1:length(beta_star))
  {
    d2_Qs[[d1]] <- list()
    for(d2 in d1:length(beta_star))
    {
      d2_Qs[[d1]][[d2]] <- matrix(0,nrow=rank,ncol=rank)
    }
    d_Qs[[d1]] <- matrix(0,nrow=rank,ncol=rank)
  }
  
  #Form the total penalty matrix
  S <- matrix(0,nrow=p,ncol=q)
  for(k in 1:length(Si))
  {
    S <- S + sp[k]*Si[[k]]
  }
  
  d_S <- d2_S <- list()
  for(d1 in 1:length(beta_star))
  {
    d_S[[d1]] <- matrix(0,nrow=p,ncol=q)
    d2_S[[d1]] <- list()
    for(d2 in d1:length(beta_star))
    {
      d2_S[[d1]][[d2]] <- matrix(0,nrow=p,ncol=q)
      for(k in 1:length(Si))
      {
        d2_S[[d1]][[d2]] <- d2_S[[d1]][[d2]] + 
          d_lsp[[k]][d2]*d_Si[[k]][[d1]] + sp[k]*d2_Si[[k]][[d1]][[d2]] +
          d2_lsp[[k]][d1,d2]*Si[[k]] + d_lsp[[k]][d1]*d_Si[[k]][[d2]]
      }
    }
    for(k in 1:m)
    {
      d_S[[d1]] <- d_S[[d1]] + sp[k]*d_Si[[k]][[d1]] + d_lsp[[k]][d1]*Si[[k]]
    }
  }
  
  while(1){
    Omega <- rep(0,times = length(gamma))  
    # Step 1
    for(k in 1:length(gamma))
    {
      if(gamma[k] == 1)
      {
        Omega[k] <- sqrt(sum(Si[[k]] * Si[[k]])) * sp[k]
      }
    }
    
    #Step 2
    alpha <- as.numeric( Omega >= dTol * max(Omega) & gamma == 1) 
    gamma_prime <- as.numeric( Omega < dTol * max(Omega) & gamma == 1)
    
    S_Scaled <- matrix(0,nrow=Q,ncol=Q)
    #Step 3
    for(k in 1:length(alpha))
    {
      if(alpha[k] == 1)
      {
        S_Scaled <- S_Scaled + Si[[k]]/sqrt(sum(Si[[k]] * Si[[k]]))
      }
    }
    
    es <- eigen(S_Scaled, symmetric = TRUE)
    ind <- es$values > max(es$values) * rTol
    penaltyRank <- sum(ind)
    
    #Step 4
    if(Q == penaltyRank){
      break
    }
    
    # Step 5
    S_Unscaled <- matrix(0,nrow=Q,ncol=Q)
    for(k in 1:m)
    {
      if(alpha[k] == 1)
      {
        S_Unscaled <- S_Unscaled + sp[k]*Si[[k]]
      }
    }
    
    d_S_Unscaled <- d2_S_Unscaled <- list()
    for(d1 in 1:length(beta_star))
    {
      d_S_Unscaled[[d1]] <- matrix(0,nrow=Q,ncol=Q)
      d2_S_Unscaled[[d1]] <- list()
      for(d2 in d1:length(beta_star))
      {
        d2_S_Unscaled[[d1]][[d2]] <- matrix(0,nrow=Q,ncol=Q)
        for(k in 1:m)
        {
          if(alpha[k] == 1)
          {
            d2_S_Unscaled[[d1]][[d2]] <- d2_S_Unscaled[[d1]][[d2]] + 
              d_lsp[[k]][d2]*d_Si[[k]][[d1]] + sp[k]*d2_Si[[k]][[d1]][[d2]] +
              d2_lsp[[k]][d1,d2]*Si[[k]] + d_lsp[[k]][d1]*d_Si[[k]][[d2]]
          }
        }
      }
      
      for(k in 1:m)
      {
        if(alpha[k] == 1)
        {
          d_S_Unscaled[[d1]] <- d_S_Unscaled[[d1]] + sp[k]*d_Si[[k]][[d1]] + d_lsp[[k]][d1]*Si[[k]]
        }
      }
    }
    
    es0 <- eigen(S_Unscaled, symmetric = TRUE)
    U0 <- es0$vectors
    evals <- es0$values
    Ur <- es0$vectors[, 1:penaltyRank]
    Un <- es0$vectors[, -c(1:penaltyRank)]
    
    U <- matrix(0,nrow = (K+Q), ncol = (K+Q))
    Id <- diag(nrow=Q,ncol=Q)
    d_U0 <- d_evals <- d_Un <- list()
    for(d1 in 1:length(beta_star))
    {
      d_U0[[d1]] <- matrix(nrow=nrow(U0),ncol=ncol(U0))
      d_evals[[d1]] <- vector()
      d_Un[[d1]] <- matrix(nrow=nrow(Un),ncol=ncol(Un))
      for(k in 1:ncol(U0))
      {
        d_U0[[d1]][,k] <- -1* ginv(S_Unscaled - evals[k]*Id) %*% (d_S_Unscaled[[d1]] %*% U0[,k])
        d_evals[[d1]][k] <- crossprod(U0,d_S_Unscaled[[d1]])[k,] %*% U0[,k]
      }
      d_Un[[d1]] <- d_U0[[d1]][, -c(1:penaltyRank)]
    }
    
    d2_U0 <- d2_Un <- list()
    for(d1 in 1:length(beta_star))
    {
      d2_U0[[d1]] <- d2_Un[[d1]] <- list()
      for(d2 in d1:length(beta_star))
      {
        d2_U0[[d1]][[d2]]<- matrix(nrow=ncol(U0),ncol=ncol(U0))
        for(k in 1:ncol(U0))
        {
          d2_U0[[d1]][[d2]][,k] <- (
            -1*(
              -ginv(S_Unscaled- evals[k]*Id) %*% ( d_S_Unscaled[[d2]]-d_evals[[d2]][k]*Id ) %*% ginv(S_Unscaled- evals[k]*Id) + 
                tcrossprod(ginv(S_Unscaled- evals[k]*Id),ginv(S_Unscaled- evals[k]*Id)) %*% crossprod(( Id-(S_Unscaled- evals[k]*Id)),ginv(S_Unscaled- evals[k]*Id)) +
                tcrossprod( ( Id-(S_Unscaled- evals[k]*Id) %*% ginv(S_Unscaled- evals[k]*Id) ), d_S_Unscaled[[d2]]-d_evals[[d2]][k]*Id ) %*% crossprod(ginv(S_Unscaled- evals[k]*Id), ginv(S_Unscaled- evals[k]*Id))
            ) %*% d_S_Unscaled[[d1]] %*% U0[,k]          
            - ginv(S_Unscaled- evals[k]*Id) %*% d2_S_Unscaled[[d1]][[d2]] %*% U0[,k]          
            - ginv(S_Unscaled- evals[k]*Id) %*% d_S_Unscaled[[d1]] %*% d_U0[[d2]][,k]          
          )
        }
        d2_Un[[d1]][[d2]] <- d2_U0[[d1]][[d2]][, -c(1:penaltyRank)]
      }
    }
    
    d_U <- d2_U <- list()
    if(K==0) {
      U <- U0
      d_U <- d_U0
      d2_U <- d2_U0
    } else if (K > 0) {
      U[1:K,1:K] <- diag(K)
      U[(K+1):(K+Q),(K+1):(K+Q)] <- U0
      
      for(d1 in 1:length(beta_star))
      {
        d2_U[[d1]] <- list()
        for(d2 in d1:length(beta_star))
        {
          d2_U[[d1]][[d2]] <- matrix(0,nrow = (K+Q), ncol = (K+Q))
          d2_U[[d1]][[d2]][(K+1):(K+Q),(K+1):(K+Q)] <- d2_U0[[d1]][[d2]]
        }
        d_U[[d1]] <- matrix(0,nrow = (K+Q), ncol = (K+Q))
        d_U[[d1]][(K+1):(K+Q),(K+1):(K+Q)] <- d_U0[[d1]] 
      }
    }
    
    StPrime <-  crossprod(U, S) %*% U
    
    drho_SPrime <- d2rho_SPrime <- list()
    for(i in 1:length(sp))
    {
      d2rho_SPrime[[i]] <- list()
      for(j in i:length(sp))
      {
        d2rho_SPrime[[i]][[j]] <- crossprod(U, d2rho_S[[i]][[j]]) %*% U
      }
      drho_SPrime[[i]] <- crossprod(U, drho_S[[i]]) %*% U
    }
    
    dStPrime <- d2StPrime <- drho_d_SPrime <- d2rho_d_SPrime <- drho_d2_SPrime <- d2rho_d2_SPrime <- list()
    for(d1 in 1:length(beta_star))
    {
      d2StPrime[[d1]] <- drho_d_SPrime[[d1]] <- d2rho_d_SPrime[[d1]] <- drho_d2_SPrime[[d1]] <- d2rho_d2_SPrime[[d1]] <- list()
      for(d2 in d1:length(beta_star))
      {
        drho_d2_SPrime[[d1]][[d2]] <- d2rho_d2_SPrime[[d1]][[d2]] <- list()
        
        d2StPrime[[d1]][[d2]] <- crossprod(d2_U[[d1]][[d2]], S) %*% U + 
          crossprod(d_U[[d1]], d_S[[d2]]) %*% U + 
          crossprod(d_U[[d1]], S) %*% d_U[[d2]] +
          
          crossprod(d_U[[d2]], d_S[[d1]]) %*% U + 
          crossprod(U, d2_S[[d1]][[d2]]) %*% U + 
          crossprod(U, d_S[[d1]]) %*% d_U[[d2]] + 
          
          crossprod(d_U[[d2]], S) %*% d_U[[d1]] +
          crossprod(U, d_S[[d2]]) %*% d_U[[d1]] + 
          crossprod(U, S) %*% d2_U[[d1]][[d2]]
        
        for(i in 1:length(sp))
        {
          d2rho_d2_SPrime[[d1]][[d2]][[i]] <- list()
          for(j in i:length(sp))
          {
            d2rho_d2_SPrime[[d1]][[d2]][[i]][[j]] <- crossprod(d2_U[[d1]][[d2]], d2rho_S[[i]][[j]]) %*% U + 
              crossprod(d_U[[d1]], d2rho_d_S[[d2]][[i]][[j]]) %*% U + 
              crossprod(d_U[[d1]], d2rho_S[[i]][[j]]) %*% d_U[[d2]] + 
              crossprod(d_U[[d2]], d2rho_d_S[[d1]][[i]][[j]]) %*% U + 
              crossprod(U, d2rho_d2_S[[d1]][[d2]][[i]][[j]]) %*% U + 
              crossprod(U, d2rho_d_S[[d1]][[i]][[j]]) %*% d_U[[d2]] + 
              crossprod(d_U[[d2]], d2rho_S[[i]][[j]]) %*% d_U[[d1]] +
              crossprod(U, d2rho_d_S[[d2]][[i]][[j]]) %*% d_U[[d1]] + 
              crossprod(U, d2rho_S[[i]][[j]]) %*% d2_U[[d1]][[d2]]
          }
          drho_d2_SPrime[[d1]][[d2]][[i]] <- crossprod(d2_U[[d1]][[d2]], drho_S[[i]]) %*% U + 
            crossprod(d_U[[d1]], drho_d_S[[d2]][[i]]) %*% U + 
            crossprod(d_U[[d1]], drho_S[[i]]) %*% d_U[[d2]] + 
            crossprod(d_U[[d2]], drho_d_S[[d1]][[i]]) %*% U + 
            crossprod(U, drho_d2_S[[d1]][[d2]][[i]]) %*% U + 
            crossprod(U, drho_d_S[[d1]][[i]]) %*% d_U[[d2]] + 
            crossprod(d_U[[d2]], drho_S[[i]]) %*% d_U[[d1]] +
            crossprod(U, drho_d_S[[d2]][[i]]) %*% d_U[[d1]] + 
            crossprod(U, drho_S[[i]]) %*% d2_U[[d1]][[d2]]
        }
      }
      
      dStPrime[[d1]] <- crossprod(d_U[[d1]], S) %*% U + 
        crossprod(U, d_S[[d1]]) %*% U + 
        crossprod(U, S) %*% d_U[[d1]]
      
      for(i in 1:length(sp))
      {
        d2rho_d_SPrime[[d1]][[i]] <- list()
        for(j in i:length(sp))
        {
          d2rho_d_SPrime[[d1]][[i]][[j]] <- crossprod(d_U[[d1]], d2rho_S[[i]][[j]]) %*% U + 
            crossprod(U, d2rho_d_S[[d1]][[i]][[j]]) %*% U + 
            crossprod(U, d2rho_S[[i]][[j]]) %*% d_U[[d1]]
        }
        drho_d_SPrime[[d1]][[i]] <- crossprod(d_U[[d1]], drho_S[[i]]) %*% U + 
          crossprod(U, drho_d_S[[d1]][[i]]) %*% U + 
          crossprod(U, drho_S[[i]]) %*% d_U[[d1]] 
      }
    }
    
    for(k in 1:length(gamma_prime))
    {
      if(gamma_prime[k] == 1){
        
        for(d1 in 1:length(beta_star))
        {
          for(d2 in d1:length(beta_star))
          {
            d2_Si[[k]][[d1]][[d2]] <- crossprod(d2_Un[[d1]][[d2]],Si[[k]]) %*% Un + 
              crossprod(d_Un[[d1]],d_Si[[k]][[d2]]) %*% Un + 
              crossprod(d_Un[[d1]],Si[[k]]) %*% d_Un[[d2]] + 
              
              crossprod(d_Un[[d2]],d_Si[[k]][[d1]]) %*% Un + 
              crossprod(Un,d2_Si[[k]][[d1]][[d2]]) %*% Un + 
              crossprod(Un,d_Si[[k]][[d1]]) %*% d_Un[[d2]] + 
              
              crossprod(d_Un[[d2]],Si[[k]]) %*% d_Un[[d1]] +
              crossprod(Un,d_Si[[k]][[d2]]) %*% d_Un[[d1]] + 
              crossprod(Un,Si[[k]]) %*% d2_Un[[d1]][[d2]]
          }
          
          d_Si[[k]][[d1]] <- crossprod(d_Un[[d1]],Si[[k]]) %*% Un + 
            crossprod(Un,d_Si[[k]][[d1]]) %*% Un + 
            crossprod(Un,Si[[k]]) %*% d_Un[[d1]]  
        }
        Si[[k]] <- crossprod(Un,Si[[k]]) %*% Un  
      }
    } 
    
    for(d1 in 1:length(beta_star))
    {
      for(d2 in d1:length(beta_star))
      {
        d2_Qs[[d1]][[d2]] <- d2_U[[d1]][[d2]] %*% Qs + d_U[[d1]] %*% d_Qs[[d2]] + d_U[[d2]] %*% d_Qs[[d1]] + U %*% d2_Qs[[d1]][[d2]]
      }
      d_Qs[[d1]] <- d_U[[d1]] %*% Qs + U %*% d_Qs[[d1]]
    }
    Qs <- U %*% Qs
    
    K <- K + penaltyRank
    Q <- Q - penaltyRank
    S <- StPrime
    d_S <- dStPrime
    d2_S <- d2StPrime
    drho_S <- drho_SPrime
    d2rho_S <- d2rho_SPrime
    drho_d_S <- drho_d_SPrime
    d2rho_d_S <- d2rho_d_SPrime
    drho_d2_S <- drho_d2_SPrime
    d2rho_d2_S <- d2rho_d2_SPrime
    gamma <- gamma_prime
  }
  
  #Calculate log of determinant
  det <- log(det(S))
  
  for(d1 in 1:length(beta_star)){
    for(d2 in d1:length(beta_star)){
      d2_rS[[d1]][[d2]] <- crossprod(d2_Qs[[d1]][[d2]],rS[[1]]) +
        crossprod(d_Qs[[d1]],d_rS[[d2]]) +
        crossprod(d_Qs[[d2]],d_rS[[d1]]) +
        crossprod(Qs,d2_rS[[d1]][[d2]])
    }
  }
  
  for(d1 in 1:length(beta_star)){
    d_rS[[d1]] <- crossprod(d_Qs[[d1]], rS[[1]]) + crossprod(Qs,d_rS[[d1]])
  }
  
  #Transform the rS and derivatives
  for(i in 1:m){
    rS[[i]] <- crossprod(Qs,rS[[i]])
  }
  
  #Transform the Si & beta_star derivatives
  Si_Q <- list()
  d_Si_Q <- d_Si
  d2_Si_Q <- d2_Si
  for(i in 1:m){
    Si_Q[[i]] <- tcrossprod(rS[[i]],rS[[i]])
    if(i == 1){
      for(d1 in 1:length(beta_star)){
        for(d2 in d1:length(beta_star)){
          d2_Si_Q[[i]][[d1]][[d2]] <- tcrossprod(d2_rS[[d1]][[d2]],rS[[i]]) +
            tcrossprod(d_rS[[d1]],d_rS[[d2]]) +
            tcrossprod(d_rS[[d2]],d_rS[[d1]]) +
            tcrossprod(rS[[i]],d2_rS[[d1]][[d2]])
        }
        d_Si_Q[[i]][[d1]] <- tcrossprod(d_rS[[d1]],rS[[i]]) + tcrossprod(rS[[i]],d_rS[[d1]])
      }
    }
  }
  
  if(deriv > 0){
    #We calculate derivatives of log|S|
    drho_det <- vector()
    d2rho_det <- matrix(0,length(sp),length(sp))
    d_det <- vector()
    d2_det <- matrix(0,length(beta_star),length(beta_star))
    drho_d_det <- d2rho_d_det <- drho_d2_det <- d2rho_d2_det <- list()
    inv_S <- solve(S)
    
    for(i in 1:length(sp)){
      drho_d_det[[i]] <- vector()
      d2rho_d_det[[i]] <- d2rho_d2_det[[i]] <- list()
      drho_d2_det[[i]] <- matrix(0,length(beta_star),length(beta_star))
      for(j in i:length(sp)){
        d2rho_d_det[[j]] <- d2rho_d2_det[[j]] <- list()
        
        d2rho_d_det[[i]][[j]] <- d2rho_d_det[[j]][[i]] <- vector()
        d2rho_d2_det[[i]][[j]] <- d2rho_d2_det[[j]][[i]] <- matrix(0,length(beta_star),length(beta_star))
        
        k_delta <- if(i==j) 1 else 0
        d2rho_det[i,j] <- k_delta * sp[i] * sum(diag(inv_S %*% Si_Q[[i]])) - 
          sp[i] * sp[j] * sum(diag(inv_S %*% Si_Q[[i]] %*% inv_S %*% Si_Q[[j]]))
        
        d2rho_det[j,i] <- d2rho_det[i,j]
        
        for(d1 in 1:length(beta_star)){
          d2rho_d_det[[i]][[j]][d1] <- sum(diag(
            inv_S %*% drho_S[[j]][1:p,1:q] %*% inv_S %*% d_S[[d1]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
              -inv_S %*% drho_d_S[[d1]][[j]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
              -inv_S %*% d_S[[d1]][1:p,1:q] %*% -inv_S %*% drho_S[[j]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
              -inv_S %*% d_S[[d1]][1:p,1:q] %*% inv_S %*% d2rho_S[[i]][[j]][1:p,1:q] +
              -inv_S %*% drho_S[[j]][1:p,1:q] %*% inv_S %*% drho_d_S[[d1]][[i]][1:p,1:q] +
              inv_S %*% d2rho_d_S[[d1]][[i]][[j]][1:p,1:q]))
          
          #d2rho_d_det[[j]][[i]][d1] <- d2rho_d_det[[i]][[j]][d1]
          
          for(d2 in d1:length(beta_star)){
            d2rho_d2_det[[i]][[j]][d1,d2] <- sum(diag(
              -inv_S %*% d_S[[d2]][1:p,1:q] %*% inv_S %*% drho_S[[j]][1:p,1:q] %*% inv_S %*% d_S[[d1]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
                inv_S %*% drho_d_S[[d2]][[j]][1:p,1:q] %*% inv_S %*% d_S[[d1]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
                inv_S %*% drho_S[[j]][1:p,1:q] %*% -inv_S %*% d_S[[d2]][1:p,1:q] %*% inv_S %*% d_S[[d1]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
                inv_S %*% drho_S[[j]][1:p,1:q] %*% inv_S %*% d2_S[[d1]][[d2]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
                inv_S %*% drho_S[[j]][1:p,1:q] %*% inv_S %*% d_S[[d1]][1:p,1:q] %*% -inv_S %*% d_S[[d2]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
                inv_S %*% drho_S[[j]][1:p,1:q] %*% inv_S %*% d_S[[d1]][1:p,1:q] %*% inv_S %*% drho_d_S[[d2]][[i]][1:p,1:q] +
                inv_S %*% d_S[[d2]][1:p,1:q] %*% inv_S %*% drho_d_S[[d1]][[j]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
                -inv_S %*% drho_d2_S[[d1]][[d2]][[j]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
                -inv_S %*% drho_d_S[[d1]][[j]][1:p,1:q] %*% -inv_S %*% d_S[[d2]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
                -inv_S %*% drho_d_S[[d1]][[j]][1:p,1:q] %*% inv_S %*% drho_d_S[[d2]][[i]][1:p,1:q] + 
                inv_S %*% d_S[[d2]][1:p,1:q] %*% inv_S %*% d_S[[d1]][1:p,1:q] %*% -inv_S %*% drho_S[[j]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
                -inv_S %*% d2_S[[d1]][[d2]][1:p,1:q] %*% -inv_S %*% drho_S[[j]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
                -inv_S %*% d_S[[d1]][1:p,1:q] %*% inv_S %*% d_S[[d2]][1:p,1:q] %*% inv_S %*% drho_S[[j]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
                -inv_S %*% d_S[[d1]][1:p,1:q] %*% -inv_S %*% drho_d_S[[d2]][[j]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
                -inv_S %*% d_S[[d1]][1:p,1:q] %*% -inv_S %*% drho_S[[j]][1:p,1:q] %*% -inv_S %*% d_S[[d2]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
                -inv_S %*% d_S[[d1]][1:p,1:q] %*% -inv_S %*% drho_S[[j]][1:p,1:q] %*% inv_S %*% drho_d_S[[d2]][[i]][1:p,1:q] +
                inv_S %*% d_S[[d2]][1:p,1:q] %*% inv_S %*% d_S[[d1]][1:p,1:q] %*% inv_S %*% d2rho_S[[i]][[j]][1:p,1:q] +
                -inv_S %*% d2_S[[d1]][[d2]][1:p,1:q] %*% inv_S %*% d2rho_S[[i]][[j]][1:p,1:q] +
                -inv_S %*% d_S[[d1]][1:p,1:q] %*% -inv_S %*% d_S[[d2]][1:p,1:q] %*% inv_S %*% d2rho_S[[i]][[j]][1:p,1:q] +
                -inv_S %*% d_S[[d1]][1:p,1:q] %*% inv_S %*% d2rho_d_S[[d2]][[i]][[j]][1:p,1:q] +
                inv_S %*% d_S[[d2]][1:p,1:q] %*% inv_S %*% drho_S[[j]][1:p,1:q] %*% inv_S %*% drho_d_S[[d1]][[i]][1:p,1:q] +
                -inv_S %*% drho_d_S[[d2]][[j]][1:p,1:q] %*% inv_S %*% drho_d_S[[d1]][[i]][1:p,1:q] +
                -inv_S %*% drho_S[[j]][1:p,1:q] %*% -inv_S %*% d_S[[d2]][1:p,1:q] %*% inv_S %*% drho_d_S[[d1]][[i]][1:p,1:q] +
                -inv_S %*% drho_S[[j]][1:p,1:q] %*% inv_S %*% drho_d2_S[[d1]][[d2]][[i]][1:p,1:q] +
                -inv_S %*% d_S[[d2]][1:p,1:q] %*% inv_S %*% d2rho_d_S[[d1]][[i]][[j]][1:p,1:q] +
                inv_S %*% d2rho_d2_S[[d1]][[d2]][[i]][[j]][1:p,1:q]))
            
            d2rho_d2_det[[i]][[j]][d2,d1] <- d2rho_d2_det[[i]][[j]][d1,d2]
          }
        }
      }
      
      drho_det[i] <- sp[i] * sum(diag(inv_S %*% Si_Q[[i]]))
      
      for(d1 in 1:length(beta_star)){
        for(d2 in d1:length(beta_star)){
          drho_d2_det[[i]][d1,d2] <- sum(diag(
            inv_S %*% d_S[[d2]][1:p,1:q] %*% inv_S %*% d_S[[d1]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
              -inv_S %*% d2_S[[d1]][[d2]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
              -inv_S %*% d_S[[d1]][1:p,1:q] %*% -inv_S %*% d_S[[d2]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
              -inv_S %*% d_S[[d1]][1:p,1:q] %*% inv_S %*% drho_d_S[[d2]][[i]][1:p,1:q] +
              -inv_S %*% d_S[[d2]][1:p,1:q] %*% inv_S %*% drho_d_S[[d1]][[i]][1:p,1:q] +
              inv_S %*% drho_d2_S[[d1]][[d2]][[i]][1:p,1:q]))
          
          drho_d2_det[[i]][d2,d1] <- drho_d2_det[[i]][d1,d2]
        }
        drho_d_det[[i]][d1] <- sum(diag(-inv_S %*% d_S[[d1]][1:p,1:q] %*% inv_S %*% drho_S[[i]][1:p,1:q] +
                                          inv_S %*% drho_d_S[[d1]][[i]][1:p,1:q]))
      }
    }
    
    for(d1 in 1:length(beta_star)){
      for(d2 in d1:length(beta_star)){
        d2_det[d1,d2] <- sum(diag(-inv_S %*% d_S[[d2]][1:p,1:q] %*% inv_S %*% d_S[[d1]][1:p,1:q] + 
                                    inv_S %*% d2_S[[d1]][[d2]][1:p,1:q]))
        
        d2_det[d2,d1] <- d2_det[d1,d2]
      }
      d_det[d1] <- sum(diag(inv_S %*% d_S[[d1]][1:p,1:q]))
    }
  }
  
  return(list(Si = Si_Q, d_Si = d_Si_Q, d2_Si = d2_Si_Q, Qs = Qs, d_Qs = d_Qs, d2_Qs = d2_Qs, S = S, 
              d_S = d_S, d2_S = d2_S, rS = rS, d_rS = d_rS, d2_rS = d2_rS, 
              drho_S = drho_S, d2rho_S = d2rho_S, drho_d_S = drho_d_S, d2rho_d_S = d2rho_d_S, 
              drho_d2_S = drho_d2_S, d2rho_d2_S = d2rho_d2_S, 
              det = det, d_det = d_det, d2_det = d2_det,
              det1 = drho_det, d_det1 = drho_d_det, d2_det1 = drho_d2_det,
              det2 = d2rho_det, d_det2 = d2rho_d_det , d2_det2 = d2rho_d2_det))
}


mroot <- function (A, rank = NULL, method = "chol"){
  if (is.null(rank)) 
    rank <- 0
  if (!isTRUE(all.equal(A, t(A)))) 
    stop("Supplied matrix not symmetric")
  if (method == "svd") {
    um <- La.svd(A)
    if (sum(um$d != sort(um$d, decreasing = TRUE)) > 0) 
      stop("singular values not returned in order")
    if (rank < 1) {
      rank <- dim(A)[1]
      if (um$d[1] <= 0) 
        rank <- 1
      else while (rank > 0 && (um$d[rank]/um$d[1] < .Machine$double.eps || 
                               all.equal(um$u[, rank], um$vt[rank, ]) != TRUE)) rank <- rank - 
          1
      if (rank == 0) 
        stop("Something wrong - matrix probably not +ve semi definite")
    }
    d <- um$d[1:rank]^0.5
    return(t(t(um$u[, 1:rank]) * as.vector(d)))
  }
  else if (method == "chol") {
    L <- suppressWarnings(chol(A, pivot = TRUE, tol = 0))
    piv <- order(attr(L, "pivot"))
    pivoting <- attr(L, "pivot")
    r <- attr(L, "rank")
    p <- ncol(L)
    if (r < p) 
      L[(r + 1):p, (r + 1):p] <- 0
    if (rank < 1) 
      rank <- r
    L <- L[, piv, drop = FALSE]
    L <- t(L[1:rank, , drop = FALSE])
    return(list(L=L,pivoting=pivoting,piv=piv))
  }
  else stop("method not recognised.")
}

d_x_Qs <- function(beta_star,x, d_x, d2_x, Qs, d_Qs, d2_Qs){
  
  for(d1 in 1:length(beta_star))
  {  
    for(d2 in d1:length(beta_star))
    {
      d2_x[[d1]][[d2]] <- d2_x[[d1]][[d2]] %*% Qs + d_x[[d1]] %*% d_Qs[[d2]] +
        (d_x[[d2]] %*% d_Qs[[d1]]) + (x %*% d2_Qs[[d1]][[d2]])
    }
    d_x[[d1]] <- d_x[[d1]] %*% Qs + x %*% d_Qs[[d1]]
  }
  list(d_x = d_x, d2_x = d2_x)
}



