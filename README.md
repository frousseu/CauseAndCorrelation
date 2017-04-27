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
  
  
 gen.data<-function (n = 100)  
 { 
     A <- rnorm(n) 
     B <- 0.5 * A + rnorm(n, 0, sqrt(1 - 0.5^2)) 
     C <- 0.5 * B + rnorm(n, 0, sqrt(1 - 0.5^2)) 
     D <- 0.5 * B + rnorm(n, 0, sqrt(1 - 0.5^2)) 
     E <- 0.5 * C + 0.5 * D + rnorm(n, 0, sqrt(1 - 2 * 0.5^2)) 
     data.frame(A = A, B = B, C = C, D = D, E = E) 
 } 
  
  
 gen.perms<-function (n = 5)  
 { 
     all <- expand.grid(p1 = 1:n, p2 = 1:n, p3 = 1:n, stringsAsFactors = FALSE) 
     perms <- all[apply(all, 1, function(x) { 
         length(unique(x)) == 3 
     }), ] 
     perms 
 } 
  
  
 MCX2<-function (model.df, n.obs, model.chi.square, n.sim = 10000)  
 { 
     x <- (-1 + sqrt(1 + 8 * model.df))/2 
     if ((x - as.integer(x)) == 0)  
         v <- x 
     if ((x - as.integer(x)) > 0 & (x - as.integer(x)) < 1)  
         v <- as.integer(x) + 1 
     if ((x - as.integer(x)) > 1)  
         return("error") 
     c.value <- v * (v + 1)/2 - model.df 
     MCX2 <- rep(NA, n.sim) 
     for (i in 1:n.sim) { 
         dat <- matrix(rnorm(n.obs * v), ncol = v) 
         obs.VCV <- var(dat) 
         model.VCV <- diag(v) 
         diag(model.VCV)[1:c.value] <- diag(obs.VCV)[1:c.value] 
         MCX2[i] <- (n.obs - 1) * (log(det(model.VCV)) + sum(diag(obs.VCV) *  
             (1/diag(model.VCV))) - log(det(obs.VCV)) - v) 
     } 
     prob <- sum(MCX2 >= model.chi.square)/n.sim 
     x <- seq(0, max(MCX2)) 
     theoretical.prob <- dchisq(x, model.df) 
     hist(MCX2, freq = F, ylab = "proportion of simulations",  
         xlab = "Maximum likelihood chi-square statistic", main = "Monte Carlo simulations",  
         ylim = c(0, max(theoretical.prob)), sub = paste(as.character(model.df),  
             " df")) 
     lines(x, theoretical.prob, lty = 2) 
     lines(x = c(model.chi.square, model.chi.square), y = c(0,  
         1), lwd = 2) 
     legend(x = "topright", legend = "theoretical X2 distribution",  
         lty = 2) 
     list(MCprobability = prob, MLprobability = 1 - pchisq(model.chi.square,  
         model.df)) 
 } 
  
  
 Observed.Equivalent.DAG<-function (full.DAG, latents = NA)  
 { 
     library(ggm) 
     pairs.without.edge <- function(my.graph) { 
         nvars <- dim(my.graph)[2] 
         com <- combn(1:nvars, 2) 
         ncombs <- dim(com)[2] 
         keep <- rep(T, ncombs) 
         for (i in 1:ncombs) { 
             if (my.graph[com[1, i], com[2, i]] != 0 | my.graph[com[2,  
                 i], com[1, i]] != 0) { 
                 com[1, i] <- com[2, i] <- 0 
                 keep[i] <- F 
             } 
         } 
         matrix(com[, keep], ncol = sum(keep)) 
     } 
     find.possible.Q <- function(nvars, x, y) { 
         z <- 1:nvars 
         z[x] <- z[y] <- 0 
         z[z > 0] 
     } 
     dag.name <- function(amat, n) { 
         rownames(amat)[n] 
     } 
     full.vars <- row.names(full.DAG) 
     full.vars.index <- 1:length(full.vars) 
     n.observed <- length(full.vars) - length(latents) 
     observed.DAG <- full.DAG 
     observed.vars <- full.vars 
     observed.vars.index <- full.vars.index 
     for (i in 1:length(latents)) { 
         observed.vars[latents[i] == full.vars] <- NA 
         observed.vars.index[latents[i] == full.vars] <- NA 
         observed.DAG[latents[i] == full.vars, ] <- NA 
         observed.DAG[, latents[i] == full.vars] <- NA 
     } 
     cat("the original DAG is:", "\n") 
     total.n.vars <- dim(full.DAG)[2] 
     for (i in 1:(total.n.vars - 1)) { 
         for (j in (i + 1):total.n.vars) { 
             if (full.DAG[i, j] == 1 & full.DAG[j, i] == 0)  
                 cat(full.vars[i], "->", full.vars[j], "\n") 
             if (full.DAG[i, j] == 0 & full.DAG[j, i] == 1)  
                 cat(full.vars[j], "->", full.vars[i], "\n") 
         } 
     } 
     if (sum(is.na(latents)) > 0) { 
         return(cat("There are no latents; the DAG doesn't change ",  
             "\n")) 
     } 
     if (sum(is.na(latents)) == 0) { 
         cat("latent variable(s): ", latents, "\n") 
         n.latents <- length(latents) 
         for (i in 1:n.latents) { 
             ok <- F 
             for (j in 1:length(full.vars)) if (latents[i] ==  
                 full.vars[j])  
                 ok <- T 
             if (!ok)  
                 return("ERROR: latent variable name not in the DAG") 
         } 
     } 
     cat("_____________________", "\n") 
     observed.vars <- observed.vars[!is.na(observed.vars)] 
     observed.vars.index <- observed.vars.index[!is.na(observed.vars.index)] 
     if (n.observed <= 0)  
         return(cat("No observed variables", "\n")) 
     if (n.observed == 1)  
         return(cat("Only one observed variable", "\n")) 
     if (n.observed == 2)  
         return(cat("Only two observed variables", "\n")) 
     observed.DAG <- observed.DAG[observed.vars.index, observed.vars.index] 
     if (n.observed <= 0) { 
         return(cat("All variables are latent; there is no equivalent observed DAG",  
             "\n")) 
     } 
     pairs.to.test <- pairs.without.edge(observed.DAG) 
     n.pairs.to.test <- dim(pairs.to.test)[2] 
     n.remaining <- length(observed.vars) - 2 
     if (n.pairs.to.test <= 0) { 
         return(cat("Since there are only two observed variables, nothing further will be done",  
             "\n")) 
     } 
     add.edge <- matrix(NA, nrow = 2, ncol = n.pairs.to.test) 
     kount <- 0 
     for (i in 1:n.pairs.to.test) { 
         is.pair.dsep <- F 
         possible.Q <- find.possible.Q(n.observed, pairs.to.test[1,  
             i], pairs.to.test[2, i]) 
         first.var <- observed.vars.index[pairs.to.test[1, i]] 
         second.var <- observed.vars.index[pairs.to.test[2, i]] 
         test <- dSep(amat = full.DAG, first = dag.name(full.DAG,  
             first.var), second = dag.name(full.DAG, second.var),  
             cond = NULL) 
         if (test) { 
             is.pair.dsep <- T 
             next 
         } 
         if (sum(is.na(possible.Q) == 0)) { 
             n.possible.Q <- length(possible.Q) 
             for (j in 1:n.possible.Q) { 
                 Q <- combn(possible.Q, j) 
                 if (j == n.possible.Q)  
                   Q <- matrix(possible.Q, nrow = j, ncol = 1) 
                 n.Q <- dim(Q)[2] 
                 first.var <- observed.vars.index[pairs.to.test[1,  
                   i]] 
                 second.var <- observed.vars.index[pairs.to.test[2,  
                   i]] 
                 for (k in 1:n.Q) { 
                   cond.vars <- as.vector(observed.vars.index[Q[,  
                     k]]) 
                   test <- dSep(amat = full.DAG, first = dag.name(full.DAG,  
                     first.var), second = dag.name(full.DAG, second.var),  
                     cond = dag.name(full.DAG, cond.vars)) 
                   if (test) { 
                     is.pair.dsep <- T 
                     break 
                   } 
                 } 
             } 
         } 
         if (!is.pair.dsep) { 
             kount <- kount + 1 
             add.edge[1, kount] <- pairs.to.test[1, i] 
             add.edge[2, kount] <- pairs.to.test[2, i] 
         } 
     } 
     cgraph <- matrix(0, n.observed, n.observed, dimnames = list(observed.vars,  
         observed.vars)) 
     for (i in 1:(n.observed - 1)) { 
         for (j in (i + 1):n.observed) { 
             if (observed.DAG[i, j] == 1 & observed.DAG[j, i] ==  
                 0) { 
                 cgraph[i, j] <- 2 
                 cgraph[j, i] <- 1 
             } 
             if (observed.DAG[j, i] == 1 & observed.DAG[i, j] ==  
                 0) { 
                 cgraph[j, i] <- 2 
                 cgraph[i, j] <- 1 
             } 
         } 
     } 
     for (i in 1:kount) { 
         cgraph[add.edge[1, i], add.edge[2, i]] <- cgraph[add.edge[2,  
             i], add.edge[1, i]] <- 1 
     } 
     cat("Equivalent partially oriented graph involving only the observed variables:",  
         "\n") 
     ind.vars <- rep(T, n.observed) 
     for (i in 1:(n.observed - 1)) { 
         for (j in (i + 1):n.observed) { 
             if (cgraph[i, j] == 2 & cgraph[j, i] == 1)  
                 cat(observed.vars[i], "->", observed.vars[j],  
                   "\n") 
             if (cgraph[i, j] == 1 & cgraph[j, i] == 2)  
                 cat(observed.vars[j], "->", observed.vars[i],  
                   "\n") 
             if (cgraph[i, j] == 1 & cgraph[j, i] == 1)  
                 cat(observed.vars[i], "--", observed.vars[j],  
                   "\n") 
             if (cgraph[i, j] > 0)  
                 ind.vars[i] <- F 
             if (cgraph[j, i] > 0)  
                 ind.vars[j] <- F 
         } 
     } 
     for (i in 1:n.observed) if (ind.vars[i])  
         cat(observed.vars[i], "--", observed.vars[i], "\n") 
     if (n.observed < 3)  
         return() 
     final.dag <- orient.graph(cgraph, n.observed) 
     new.dag <- matrix(0, n.observed, n.observed) 
     for (i in 1:(n.observed - 1)) { 
         for (j in (i + 1):n.observed) { 
             if (final.dag[i, j] == 2 & final.dag[j, i] == 1)  
                 new.dag[i, j] <- 1 
             if (final.dag[i, j] == 1 & final.dag[j, i] == 2)  
                 new.dag[j, 1] <- 1 
         } 
     } 
     test.acyclic <- isAcyclic(new.dag) 
     if (!test.acyclic) { 
         return(cat("No possible DAG; orientation requires at least one latent",  
             "\n")) 
     } 
     if (test.acyclic)  
         cat("One *possible* equivalent DAG involving only the observed variables:",  
             "\n") 
     for (i in 1:(n.observed - 1)) { 
         for (j in (i + 1):n.observed) { 
             if (final.dag[i, j] == 2 & final.dag[j, i] == 1)  
                 cat(observed.vars[i], "->", observed.vars[j],  
                   "\n") 
             if (final.dag[i, j] == 1 & final.dag[j, i] == 2)  
                 cat(observed.vars[j], "->", observed.vars[i],  
                   "\n") 
             if (final.dag[i, j] == 1 & final.dag[j, i] == 1)  
                 cat(observed.vars[i], "--", observed.vars[j],  
                   "\n") 
         } 
     } 
     for (i in 1:n.observed) if (ind.vars[i])  
         cat(observed.vars[i], "--", observed.vars[i], "\n") 
 } 
  
  
 orient.graph<-function (cgraph, nvars)  
 { 
     find.undirected.edge <- function(cgraph, nvars) { 
         if (nvars < 3)  
             stop("error in find.undirected") 
         x1 <- y1 <- 0 
         for (i in 1:(nvars - 1)) { 
             for (j in (i + 1):nvars) { 
                 if (cgraph[i, j] == 1 & cgraph[j, i] == 1) { 
                   x1 <- i 
                   y1 <- j 
                   return(data.frame(x1 = x1, y1 = y1)) 
                 } 
             } 
         } 
         data.frame(x1 = x1, y1 = y1) 
     } 
     old.graph <- cgraph 
     is.change1 <- is.change2 <- T 
     while (is.change1 | is.change2) { 
         new.graph <- orient.phaseII.1(old.graph, nvars) 
         if (sum(old.graph != new.graph) == 0)  
             is.change1 <- F 
         if (sum(old.graph != new.graph) > 0)  
             is.change1 <- T 
         new.edge <- find.undirected.edge(new.graph, nvars) 
         if (new.edge$x1 > 0 & new.edge$y1 > 0) { 
             new.graph[new.edge$x1, new.edge$y1] <- 2 
             old.graph <- new.graph 
             is.change2 <- T 
         } 
         if (new.edge$x1 == 0 & new.edge$y1 == 0) { 
             old.graph <- new.graph 
             is.change2 <- F 
         } 
     } 
     amat <- new.graph 
     amat[new.graph == 1] <- 0 
     amat[new.graph == 2] <- 1 
     new.graph 
 } 
  
  
 orient.phaseII.1<-function (cgraph, nvars)  
 { 
     if (nvars < 3)  
         return(cgraph) 
     var.trips <- combn(1:nvars, 3) 
     ntrips <- dim(var.trips)[2] 
     all.triplets <- expand.grid(p1 = 1:3, p2 = 1:3, p3 = 1:3,  
         stringsAsFactors = FALSE) 
     triplet.perms <- all.triplets[apply(all.triplets, 1, function(x) { 
         length(unique(x)) == 3 
     }), ] 
     n.triplet.perms <- dim(triplet.perms)[1] 
     if (nvars > 3) { 
         var.quads <- combn(1:nvars, 4) 
         nquads <- dim(var.quads)[2] 
         all.quads <- expand.grid(p1 = 1:4, p2 = 1:4, p3 = 1:4,  
             p4 = 1:4, stringsAsFactors = FALSE) 
         quad.perms <- all.quads[apply(all.quads, 1, function(x) { 
             length(unique(x)) == 4 
         }), ] 
         n.quad.perms <- dim(quad.perms)[1] 
     } 
     old <- cgraph 
     test <- 1 
     while (test != 0) { 
         for (i in 1:ntrips) { 
             for (j in 1:n.triplet.perms) { 
                 x <- triplet.perms[j, 1] 
                 y <- triplet.perms[j, 2] 
                 z <- triplet.perms[j, 3] 
                 if (cgraph[var.trips[x, i], var.trips[y, i]] ==  
                   2 & cgraph[var.trips[y, i], var.trips[z, i]] ==  
                   1 & cgraph[var.trips[z, i], var.trips[y, i]] ==  
                   1 & cgraph[var.trips[x, i], var.trips[z, i]] ==  
                   0)  
                   cgraph[var.trips[y, i], var.trips[z, i]] <- 2 
                 if (cgraph[var.trips[x, i], var.trips[y, i]] ==  
                   2 & cgraph[var.trips[y, i], var.trips[z, i]] ==  
                   2 & cgraph[var.trips[x, i], var.trips[z, i]] ==  
                   1 & cgraph[var.trips[z, i], var.trips[x, i]] ==  
                   1)  
                   cgraph[var.trips[x, i], var.trips[z, i]] <- 2 
             } 
         } 
         if (nvars < 4)  
             return(cgraph) 
         for (i in 1:nquads) { 
             for (j in 1:n.quad.perms) { 
                 w <- quad.perms[j, 1] 
                 x <- quad.perms[j, 2] 
                 y <- quad.perms[j, 3] 
                 z <- quad.perms[j, 4] 
                 if (cgraph[var.quads[w, i], var.quads[x, i]] ==  
                   1 & cgraph[var.quads[x, i], var.quads[w, i]] ==  
                   1 & cgraph[var.quads[w, i], var.quads[y, i]] ==  
                   1 & cgraph[var.quads[y, i], var.quads[w, i]] ==  
                   1 & cgraph[var.quads[w, i], var.quads[z, i]] ==  
                   1 & cgraph[var.quads[z, i], var.quads[w, i]] ==  
                   1 & cgraph[var.quads[x, i], var.quads[z, i]] ==  
                   2 & cgraph[var.quads[z, i], var.quads[x, i]] ==  
                   1 & cgraph[var.quads[y, i], var.quads[z, i]] ==  
                   2 & cgraph[var.quads[z, i], var.quads[y, i]] ==  
                   1)  
                   cgraph[var.quads[w, i], var.quads[z, i]] <- 2 
                 if (cgraph[var.quads[w, i], var.quads[x, i]] ==  
                   1 & cgraph[var.quads[x, i], var.quads[w, i]] ==  
                   1 & cgraph[var.quads[w, i], var.quads[y, i]] ==  
                   1 & cgraph[var.quads[y, i], var.quads[w, i]] ==  
                   1 & cgraph[var.quads[w, i], var.quads[z, i]] >  
                   0 & cgraph[var.quads[z, i], var.quads[w, i]] >  
                   0 & cgraph[var.quads[x, i], var.quads[z, i]] ==  
                   1 & cgraph[var.quads[z, i], var.quads[x, i]] ==  
                   2 & cgraph[var.quads[y, i], var.quads[z, i]] ==  
                   2 & cgraph[var.quads[z, i], var.quads[y, i]] ==  
                   1)  
                   cgraph[var.quads[w, i], var.quads[x, i]] <- 2 
             } 
             test <- sum(old != cgraph) 
             if (test != 0)  
                 old <- cgraph 
         } 
     } 
     cgraph 
 } 
  
  
 Pcor.prob<-function (dat, x, y, Q)  
 { 
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
  
  
 shipley.test2<-function (amat, S, n)  
 { 
     pval <- function(r, q, n) { 
         df = n - 2 - q 
         tval <- r * sqrt(df)/sqrt(1 - r * r) 
         2 * (1 - pt(abs(tval), df)) 
     } 
     l <- basiSet(amat) 
     k <- length(l) 
     p <- rep(0, k) 
     cat("testing individual d-sep claims in basis set", "\n") 
     for (i in 1:k) { 
         r <- pcor(l[[i]], S) 
         q <- length(l[[i]]) - 2 
         p[i] <- pval(r, q, n) 
         if (is.nan(p[i]))  
             return(list(r = r, q = q, n = n, p = p[i])) 
         cat(l[[i]][1], "_||_", l[[i]][2], "|{", l[[i]][-c(1,  
             2)], "} r=", round(r, 3), " p=", round(p[i], 3),  
             "\n") 
     } 
     ctest <- -2 * sum(log(p)) 
     df <- 2 * k 
     pv <- 1 - pchisq(ctest, df) 
     list(ctest = ctest, df = df, pvalue = pv) 
 } 
  
  
 vanishing.tetrads<-function (dat, sig = 0.05)  
 { 
     get.3.equations <- function(tet.vector) { 
         mat <- matrix(NA, ncol = 8, nrow = 3) 
         mat[1, ] <- cbind(tet.vector[1], tet.vector[2], tet.vector[3],  
             tet.vector[4], tet.vector[1], tet.vector[4], tet.vector[2],  
             tet.vector[3]) 
         mat[2, ] <- cbind(tet.vector[1], tet.vector[3], tet.vector[2],  
             tet.vector[4], tet.vector[1], tet.vector[4], tet.vector[2],  
             tet.vector[3]) 
         mat[3, ] <- cbind(tet.vector[1], tet.vector[3], tet.vector[2],  
             tet.vector[4], tet.vector[1], tet.vector[2], tet.vector[3],  
             tet.vector[4]) 
         mat 
     } 
     test.stat <- function(dat, triplet) { 
         t.vars <- sort(triplet[1:4]) 
         r <- var(dat, na.rm = T) 
         tao <- r[triplet[1], triplet[2]] * r[triplet[3], triplet[4]] -  
             r[triplet[5], triplet[6]] * r[triplet[7], triplet[8]] 
         D13 <- det(r[c(triplet[1], triplet[3]), c(triplet[1],  
             triplet[3])]) 
         D24 <- det(r[c(triplet[2], triplet[4]), c(triplet[2],  
             triplet[4])]) 
         D <- det(r[triplet[1:4], triplet[1:4]]) 
         N <- dim(dat)[1] 
         tao.var <- (D13 * D24 * (N + 1)/(N - 1) - D) * (1/(N -  
             2)) 
         if (tao.var <= 0) { 
             cat("triplet: ", triplet, "\n") 
             cat("variance of tao is ", tao.var, "\n") 
             cat("tao.var<=0. D=", D, "D13=", D13, "D24=", D24,  
                 "\n") 
             stop() 
         } 
         z <- tao/sqrt(tao.var) 
         list(triplet = triplet, VCV = r, tao = tao, tao.var = tao.var,  
             z = z, prob = 2 * (1 - pnorm(abs(z)))) 
     } 
     get.choke.points <- function(vec) { 
         tetrad <- matrix(vec, ncol = 2, byrow = T) 
         all.comb <- cbind(c(vec[1], vec[1], vec[1], vec[2], vec[2],  
             vec[3]), c(vec[2], vec[3], vec[4], vec[3], vec[4],  
             vec[4])) 
         chokes <- rep(T, 6) 
         for (j in 1:4) { 
             for (i in 1:6) { 
                 if (sum(tetrad[j, ] == all.comb[i, c(1, 2)]) ==  
                   2)  
                   chokes[i] <- F 
                 if (sum(tetrad[j, ] == all.comb[i, c(2, 1)]) ==  
                   2)  
                   chokes[i] <- F 
             } 
         } 
         list(tetrad = tetrad, all.comb = all.comb, choke.points = all.comb[chokes,  
             ]) 
     } 
     nvars <- dim(dat)[2] 
     tetrad.quadriplets <- combn(1:nvars, 4) 
     ntetrads <- dim(tetrad.quadriplets)[2] 
     z <- prob <- rep(NA, ntetrads * 3) 
     count <- 0 
     for (i in 1:ntetrads) { 
         triplets <- get.3.equations(tetrad.quadriplets[, i]) 
         for (j in 1:3) { 
             count <- count + 1 
             temp <- test.stat(dat, triplets[j, ]) 
             z[count] <- temp$z 
             prob[count] <- temp$prob 
             if (prob[count] <= sig)  
                 cat("triplet:", triplets[j, ], " does not vanish (p=",  
                   prob[count], ") \n\n") 
             if (prob[count] > sig) { 
                 chokes <- get.choke.points(triplets[j, ]) 
                 cat("triplet:", triplets[j, ], "  vanishes (p=",  
                   prob[count], ") \n") 
                 cat("If there is a saturated dependency graph for the four variables (via EPA):",  
                   triplets[j, 1], triplets[j, 2], triplets[j,  
                     3], triplets[j, 4], "\n") 
                 cat("then there is at least one latent common cause of either (",  
                   chokes$choke.points[1, 1], ",", chokes$choke.points[1,  
                     2], ") and/or of (", chokes$choke.points[2,  
                     1], ",", chokes$choke.points[2, 2], ")\n\n") 
             } 
         } 
     } 
 } 
 
```
