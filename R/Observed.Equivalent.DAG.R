Observed.Equivalent.DAG<-function (full.DAG, latents = NA) 
{
    library(ggm)
    pairs.without.edge <- function(my.graph) {
        nvars <- dim(my.graph)[2]
        com <- combn(1:nvars, 2)
        ncombs <- dim(com)[2]
        keep <- rep(1, ncombs)
        for (i in 1:ncombs) {
            if (my.graph[com[1, i], com[2, i]] != 0 | my.graph[com[2, 
                i], com[1, i]] != 0) {
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
    }
    cat("_____________________", "\n")
    observed.vars <- observed.vars[!is.na(observed.vars)]
    observed.vars.index <- observed.vars.index[!is.na(observed.vars.index)]
    observed.DAG <- observed.DAG[observed.vars.index, observed.vars.index]
    if (dim(observed.DAG)[2] == 0) {
        return(cat("All variables are latent; there is no equivalent observed DAG", 
            "\n"))
    }
    pairs.to.test <- pairs.without.edge(observed.DAG)
    n.pairs.to.test <- dim(pairs.to.test)[2]
    n.remaining <- length(observed.vars) - 2
    if (n.remaining <= 0) {
        cat("The equivalent DAG is ", "\n")
        return(observed.DAG)
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
