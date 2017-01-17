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
