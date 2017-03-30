# CauseAndCorrelation
Package of functions for the path analysis [summer school course](http://www.billshipley.recherche.usherbrooke.ca/summer%20school%20path%20analysis.htm) by [Bill Shipley](http://www.billshipley.recherche.usherbrooke.ca/)

## Installation

First, make sure that you have the [latest version of R](https://cran.r-project.org/) on your computer, or at least a very recent version. To install the package from GitHub, type the following commands in your [R](https://cran.r-project.org/)  or [RStudio](https://www.rstudio.com/products/RStudio/#Desktop) console:

```r
# install package dependencies first
install.packages("lavaan")
install.packages("ggm")

# install and load devtools to be able to install packages from GitHub with install_github
install.packages("devtools")
library(devtools)

# install CauseAndCorrelation from Bill's GitHub
install_github("BillShipley/CauseAndCorrelation")
library(CauseAndCorrelation)
?Causal.Inference

```
If everything worked, the last command should have opened the help file for the `Causal.Inference` function. 


## If things do not work


### Step 1

Read the error messages and make sure all packages dependencies are installed and loaded, especially package `lavaan` and `ggm`. If a message says that a package could not be loaded, try installing it manually by typing:

```r
# manually installing a dependencies
install.packages("packagename")
```
until all packages are installed. A current bug in `install_github` on Windows prevents the installation of package dependencies of dependencies (`lavaan` and `ggm`).


### Step 2

Although this should not be required, for certain packages with compiled code, [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and [MiKTeX](https://miktex.org/) need to be installed on Windows to be able to build source packages. For Mac users, Xcode is required and can be installed through the apple store. Here is a more detailed list of [prerequisites](https://support.rstudio.com/hc/en-us/articles/200486498-Package-Development-Prerequisites) for building source packages for Windows, Mac and Linux.


### Step 3

In case something goes wrong with the package installation and the previous instructions do not work, here is the code for every function in the package. The code can be copied and pasted in the R console to get the definition of all functions. However, in this case, help files won't be available.


```r
Causal.Inference<-function (dat, alpha.reject = 0.05, write.result = T) 
{
    x <- apply(dat, 1, function(x) {
        sum(is.na(x)) == 0
    })
    Pcor.prob <- function(dat, x, y, Q) {
        if (is.na(x) | is.na(y)) 
            stop("ERROR in Pcor.prob of Causal.Inference; x or y missing")
        if (sum(is.na(Q)) == 0 & any(x == Q)) {
            cat("x= ", x, " Q= ", Q, "\n")
            stop("ERROR in Pcor.prob of Causal.Inference; x=Q")
        }
        if (sum(is.na(Q)) == 0 & any(y == Q)) 
            stop("ERROR in Pcor.prob of Causal.Inference; y=Q")
        if (x == y) 
            stop("ERROR in Causal.Inference; x=y")
        if (sum(is.na(Q)) == 0) {
            new.mat <- cbind(dat[, x], dat[, y], dat[, Q])
            n.cond <- length(Q)
        }
        else {
            new.mat <- cbind(dat[, x], dat[, y])
            n.cond <- 0
        }
        n <- dim(new.mat)[1]
        df <- n - 2 - length(n.cond)
        inv <- try(solve(var(new.mat, na.rm = T)))
        r.value <- -1 * inv[1, 2]/sqrt(inv[1, 1] * inv[2, 2])
        t.value <- r.value * sqrt((n - 2) * (1 - r.value^2))
        2 * (1 - pt(abs(t.value), df))
    }
    pairs.with.edge <- function(cgraph) {
        nvars <- dim(cgraph)[2]
        com <- combn(1:nvars, 2)
        ncombs <- dim(com)[2]
        keep <- rep(1, ncombs)
        for (i in 1:ncombs) {
            if (cgraph[com[1, i], com[2, i]] == 0) {
                com[1, i] <- com[2, i] <- 0
            }
        }
        com[, com[1, ] > 0]
    }
    find.possible.Q <- function(nvars, x, y) {
        z <- 1:nvars
        z[x] <- z[y] <- 0
        z[z > 0]
    }
    nvars <- dim(dat)[2]
    cgraph <- matrix(1, nrow = nvars, ncol = nvars)
    diag(cgraph) <- rep(0, nvars)
    do.pairs <- pairs.with.edge(cgraph)
    n.pairs <- dim(do.pairs)[2]
    if (n.pairs > 0) {
        for (j in 1:n.pairs) {
            p <- Pcor.prob(dat, x = do.pairs[1, j], y = do.pairs[2, 
                j], Q = NA)
            if (p > alpha.reject) {
                cgraph[do.pairs[1, j], do.pairs[2, j]] <- 0
                cgraph[do.pairs[2, j], do.pairs[1, j]] <- 0
            }
        }
    }
    max.order <- nvars - 2
    for (i in 1:max.order) {
        do.pairs <- pairs.with.edge(cgraph)
        if (is.vector(do.pairs)) 
            do.pairs <- matrix(do.pairs, ncol = 1)
        n.pairs <- dim(do.pairs)[2]
        if (n.pairs > 0) {
            for (j in 1:n.pairs) {
                Q <- find.possible.Q(nvars = nvars, x = do.pairs[1, 
                  j], y = do.pairs[2, j])
                if (length(Q) > 1) 
                  x <- combn(Q, i)
                else x <- as.matrix(Q, 1, 1)
                for (k in 1:length(x[1, ])) {
                  x1 <- do.pairs[1, j]
                  y1 <- do.pairs[2, j]
                  Qcond <- x[, k]
                  p <- Pcor.prob(dat, x = x1, y = y1, Q = Qcond)
                  if (p > alpha.reject) {
                    cgraph[do.pairs[1, j], do.pairs[2, j]] <- 0
                    cgraph[do.pairs[2, j], do.pairs[1, j]] <- 0
                  }
                }
            }
        }
    }
    triplets <- combn(1:nvars, 3)
    n.triplets <- dim(triplets)[2]
    for (i in 1:n.triplets) {
        X <- Y <- Z <- 0
        if (cgraph[triplets[1, i], triplets[2, i]] > 0 & cgraph[triplets[2, 
            i], triplets[3, i]] > 0 & cgraph[triplets[1, i], 
            triplets[3, i]] == 0) {
            X <- triplets[1, i]
            Y <- triplets[2, i]
            Z <- triplets[3, i]
        }
        if (cgraph[triplets[1, i], triplets[3, i]] > 0 & cgraph[triplets[1, 
            i], triplets[2, i]] == 0 & cgraph[triplets[3, i], 
            triplets[2, i]] > 0) {
            X <- triplets[1, i]
            Z <- triplets[2, i]
            Y <- triplets[3, i]
        }
        if (cgraph[triplets[1, i], triplets[3, i]] > 0 & cgraph[triplets[1, 
            i], triplets[2, i]] > 0 & cgraph[triplets[2, i], 
            triplets[3, i]] == 0) {
            Y <- triplets[1, i]
            X <- triplets[2, i]
            Z <- triplets[3, i]
        }
        if (X > 0 & Y > 0 & Z > 0) {
            flag <- 0
            p <- Pcor.prob(dat, x = X, y = Z, Q = Y)
            if (p > alpha.reject & nvars > 3) 
                flag <- 0
            if (p > alpha.reject & nvars == 3) 
                flag <- 1
            if (nvars > 3) {
                var.set <- (1:nvars)[-c(X, Y, Z)]
                corder <- 1
                ncond <- length(var.set)
                while (flag == 0 & corder <= ncond) {
                  if (ncond == 1) 
                    cset <- matrix(var.set, 1, 1)
                  if (ncond > 1) 
                    cset <- combn(var.set, corder)
                  ncset <- dim(cset)[2]
                  for (i2 in 1:ncset) {
                    p <- Pcor.prob(dat, x = X, y = Z, Q = c(Y, 
                      cset[, i2]))
                    if (p > alpha.reject) 
                      flag <- 1
                  }
                  corder <- corder + 1
                }
            }
            if (flag == 0) 
                cgraph[X, Y] <- cgraph[Z, Y] <- 2
        }
    }
    EPA.write <- function(cgraph, dat) {
        nvars <- dim(cgraph)[1]
        if (!is.null(names(dat))) 
            var.names <- names(dat)
        if (is.null(names(dat))) 
            var.names <- 1:nvars
        npossible <- factorial(nvars)/(factorial(nvars - 2) * 
            2)
        count <- 0
        for (i in 1:(nvars - 1)) {
            for (j in (i + 1):nvars) {
                if (cgraph[i, j] > 0 | cgraph[j, i] > 0) 
                  count <- count + 1
                if (count > npossible) 
                  return("ERROR")
                if (cgraph[i, j] == 1 & cgraph[j, i] == 1) {
                  cat(var.names[i], "--", var.names[j], "\n")
                }
                if (cgraph[i, j] == 2 & cgraph[j, i] == 1) {
                  cat(var.names[i], "->", var.names[j], "\n")
                }
                if (cgraph[j, i] == 2 & cgraph[i, j] == 1) {
                  cat(var.names[i], "<-", var.names[j], "\n")
                }
                if (cgraph[j, i] == 2 & cgraph[i, j] == 2) {
                  cat(var.names[i], "<->", var.names[j], "\n")
                }
            }
        }
        out <- apply(cgraph, 2, sum)
        for (i in 1:nvars) if (out[i] == 0) 
            cat(var.names[i], "-- none\n")
    }
    if (write.result) 
        EPA.write(cgraph, dat)
    if (!write.result) 
        cgraph
}

 
dsep.test<-function (amat, S, n, only.null = F) 
{
    pval <- function(r, q, n) {
        df = n - 2 - q
        tval <- r * sqrt(df)/sqrt(1 - r * r)
        2 * (1 - pt(abs(tval), df))
    }
    l <- basiSet(amat)
    k <- length(l)
    p <- rep(0, k)
    if (!only.null) 
        cat("Individual d-sep claims in basis set", "\n")
    for (i in 1:k) {
        r <- pcor(l[[i]], S)
        q <- length(l[[i]]) - 2
        p[i] <- pval(r, q, n)
        if (is.nan(p[i])) 
            return(list(r = r, q = q, n = n, p = p[i]))
        if (!only.null) 
            cat(l[[i]][1], "_||_", l[[i]][2], "|{", l[[i]][-c(1, 
                2)], "} r=", round(r, 3), " p=", round(p[i], 
                3), "\n")
    }
    ctest <- -2 * sum(log(p))
    df <- 2 * k
    pv <- 1 - pchisq(ctest, df)
    if (only.null) 
        pv
    else list(ctest = ctest, df = df, pvalue = pv)
}

 
EPA<-function (dat, alpha.reject = 0.05, write.result = T) 
{
    x <- apply(dat, 1, function(x) {
        sum(is.na(x)) == 0
    })
    Pcor.prob <- function(dat, x, y, Q) {
        if (sum(is.na(Q)) == 0) {
            new.mat <- cbind(dat[, x], dat[, y], dat[, Q])
            n.cond <- length(Q)
        }
        else {
            new.mat <- cbind(dat[, x], dat[, y])
            n.cond <- 0
        }
        n <- dim(new.mat)[1]
        df <- n - 2 - length(n.cond)
        inv <- try(solve(var(new.mat, na.rm = T)))
        r.value <- -1 * inv[1, 2]/sqrt(inv[1, 1] * inv[2, 2])
        t.value <- r.value * sqrt((n - 2) * (1 - r.value^2))
        2 * (1 - pt(abs(t.value), df))
    }
    pairs.with.edge <- function(cgraph) {
        com <- combn(1:nvars, 2)
        ncombs <- dim(com)[2]
        keep <- rep(1, ncombs)
        for (i in 1:ncombs) {
            if (cgraph[com[1, i], com[2, i]] == 0) {
                com[1, i] <- com[2, i] <- 0
            }
        }
        com[, com[1, ] > 0]
    }
    find.possible.Q <- function(nvars, x, y) {
        z <- 1:nvars
        z[x] <- z[y] <- 0
        z[z > 0]
    }
    nvars <- dim(dat)[2]
    cgraph <- matrix(1, nrow = nvars, ncol = nvars)
    diag(cgraph) <- rep(0, nvars)
    do.pairs <- pairs.with.edge(cgraph)
    n.pairs <- dim(do.pairs)[2]
    if (n.pairs > 0) {
        for (j in 1:n.pairs) {
            p <- Pcor.prob(dat, x = do.pairs[1, j], y = do.pairs[2, 
                j], Q = NA)
            if (p > alpha.reject) {
                cgraph[do.pairs[1, j], do.pairs[2, j]] <- 0
                cgraph[do.pairs[2, j], do.pairs[1, j]] <- 0
            }
        }
    }
    max.order <- nvars - 2
    for (i in 1:max.order) {
        do.pairs <- pairs.with.edge(cgraph)
        if (is.vector(do.pairs)) 
            do.pairs <- matrix(do.pairs, ncol = 1)
        n.pairs <- dim(do.pairs)[2]
        if (n.pairs > 0) {
            for (j in 1:n.pairs) {
                Q <- find.possible.Q(nvars, x = do.pairs[1, j], 
                  y = do.pairs[2, j])
                x <- combn(Q, i)
                for (k in 1:length(x[1, ])) {
                  x1 <- do.pairs[1, j]
                  y1 <- do.pairs[2, j]
                  Qcond <- x[, k]
                  p <- Pcor.prob(dat, x = x1, y = y1, Q = Qcond)
                  if (p > alpha.reject) {
                    cgraph[do.pairs[1, j], do.pairs[2, j]] <- 0
                    cgraph[do.pairs[2, j], do.pairs[1, j]] <- 0
                  }
                }
            }
        }
    }
    triplets <- combn(1:nvars, 3)
    n.triplets <- dim(triplets)[2]
    for (i in 1:n.triplets) {
        X <- Y <- Z <- 0
        if (cgraph[triplets[1, i], triplets[2, i]] > 0 & cgraph[triplets[2, 
            i], triplets[3, i]] > 0 & cgraph[triplets[1, i], 
            triplets[3, i]] == 0) {
            X <- triplets[1, i]
            Y <- triplets[2, i]
            Z <- triplets[3, i]
        }
        if (cgraph[triplets[1, i], triplets[3, i]] > 0 & cgraph[triplets[1, 
            i], triplets[2, i]] == 0 & cgraph[triplets[3, i], 
            triplets[2, i]] > 0) {
            X <- triplets[1, i]
            Z <- triplets[2, i]
            Y <- triplets[3, i]
        }
        if (cgraph[triplets[1, i], triplets[3, i]] > 0 & cgraph[triplets[1, 
            i], triplets[2, i]] > 0 & cgraph[triplets[2, i], 
            triplets[3, i]] == 0) {
            Y <- triplets[1, i]
            X <- triplets[2, i]
            Z <- triplets[3, i]
        }
        if (X > 0 & Y > 0 & Z > 0) {
            var.set <- (1:nvars)[-c(X, Y, Z)]
            flag <- 0
            p <- Pcor.prob(dat, x = X, y = Z, Q = Y)
            if (p > alpha.reject) 
                flag <- 0
            corder <- 1
            ncond <- length(var.set)
            while (flag == 0 & corder <= ncond) {
                if (ncond == 1) 
                  cset <- matrix(var.set, 1, 1)
                if (ncond > 1) 
                  cset <- combn(var.set, corder)
                ncset <- dim(cset)[2]
                for (i2 in 1:ncset) {
                  p <- Pcor.prob(dat, x = X, y = Z, Q = c(Y, 
                    cset[, i2]))
                  if (p > alpha.reject) 
                    flag <- 1
                }
                corder <- corder + 1
            }
            if (flag == 0) 
                cgraph[X, Y] <- cgraph[Z, Y] <- 2
        }
    }
    EPA.write <- function(cgraph, dat) {
        nvars <- dim(cgraph)[1]
        if (!is.null(names(dat))) 
            var.names <- names(dat)
        if (is.null(names(dat))) 
            var.names <- 1:nvars
        npossible <- factorial(nvars)/(factorial(nvars - 2) * 
            2)
        count <- 0
        for (i in 1:(nvars - 1)) {
            for (j in (i + 1):nvars) {
                if (cgraph[i, j] > 0 | cgraph[j, i] > 0) 
                  count <- count + 1
                if (count > npossible) 
                  return("ERROR")
                if (cgraph[i, j] == 1 & cgraph[j, i] == 1) {
                  cat(var.names[i], "--", var.names[j], "\n")
                }
                if (cgraph[i, j] == 2 & cgraph[j, i] == 1) {
                  cat(var.names[i], "->", var.names[j], "\n")
                }
                if (cgraph[j, i] == 2 & cgraph[i, j] == 1) {
                  cat(var.names[i], "<-", var.names[j], "\n")
                }
                if (cgraph[j, i] == 2 & cgraph[i, j] == 2) {
                  cat(var.names[i], "<->", var.names[j], "\n")
                }
            }
        }
        out <- apply(cgraph, 2, sum)
        for (i in 1:nvars) if (out[i] == 0) 
            cat(var.names[i], "-- none\n")
    }
    if (write.result) 
        EPA.write(cgraph, dat)
    if (!write.result) 
        cgraph
}

 
Exploratory.path.analysis<-function (dat, upper.bound = 0.5, significance.level = 0.05) 
{
    library(ggm)
    EPA.write <- function(cgraph, dat) {
        nvars <- dim(cgraph)[1]
        if (!is.null(names(dat))) 
            var.names <- names(dat)
        if (is.null(names(dat))) 
            var.names <- 1:nvars
        npossible <- factorial(nvars)/(factorial(nvars - 2) * 
            2)
        count <- 0
        for (i in 1:(nvars - 1)) {
            for (j in (i + 1):nvars) {
                if (cgraph[i, j] > 0 | cgraph[j, i] > 0) 
                  count <- count + 1
                if (count > npossible) 
                  return("ERROR")
                if (cgraph[i, j] == 1 & cgraph[j, i] == 1) {
                  cat(var.names[i], "--", var.names[j], "\n")
                }
                if (cgraph[i, j] == 2 & cgraph[j, i] == 1) {
                  cat(var.names[i], "->", var.names[j], "\n")
                }
                if (cgraph[j, i] == 2 & cgraph[i, j] == 1) {
                  cat(var.names[i], "<-", var.names[j], "\n")
                }
                if (cgraph[j, i] == 2 & cgraph[i, j] == 2) {
                  cat(var.names[i], "<->", var.names[j], "\n")
                }
            }
        }
        out <- apply(cgraph, 2, sum)
        for (i in 1:nvars) if (out[i] == 0) 
            cat(var.names[i], "-- none\n")
    }
    VCV <- var(dat)
    n.obs <- dim(dat)[1]
    nvars <- dim(VCV)[2]
    explore.range <- c(0.01, seq(0.05, upper.bound, 0.05))
    n.ranges <- length(explore.range)
    old <- matrix(-1, nvars, nvars)
    for (i in 1:n.ranges) {
        cgraph <- Causal.Inference(dat, alpha.reject = explore.range[i], 
            write.result = F)
        same.cgraph <- sum(cgraph != old)
        if (same.cgraph > 0) {
            old <- cgraph
            dag <- orient.graph(cgraph, nvars)
            amat <- dag
            amat[dag == 1] <- 0
            amat[dag == 2] <- 1
            if (!isAcyclic(amat)) 
                next
            dimnames(amat) <- list(names(dat), names(dat))
            p.val <- round(dsep.test(amat = amat, S = VCV, n = n.obs, 
                only.null = T), 4)
            if (p.val > significance.level) {
                cat("Partially oriented graph:", "\n")
                Causal.Inference(dat, alpha.reject = explore.range[i], 
                  write.result = T)
                cat("Null probability for partially oriented graph ", 
                  as.character(p.val), "\n")
                cat("obtained at boundary level ", as.character(explore.range[i]), 
                  "\n")
                cat("An equivalent DAG:", "\n")
                EPA.write(dag, dat)
                cat("________", "\n", "\n")
            }
        }
    }
}
```