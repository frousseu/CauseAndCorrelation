Observed.Equivalent.DAG<-function (full.DAG, latents = NA) 
{
# needs the ggm library
# full.DAG is a binary (0/1) matrix produced from DAG() function in ggm
#
# combn gives all unique combinations of a vector, taken
# 2 at a time, and outputs a matrix with each unique combination in a
# column and the total number of columns equal to the total number of
# unique combinations.
#
#pairs.without.edge outputs, as a matrix, the number of pairs of variables in my.graph
# that don't share an edge, with one column
# per pair and the two rows giving the variables in the pair.
    pairs.without.edge <- function(my.graph) {
        nvars<-dim(my.graph)[2]
        com <- combn(1:nvars, 2)
        ncombs <- dim(com)[2]
        keep <- rep(T, ncombs)
        for (i in 1:ncombs) {
# if(there is an edge between this pair) then remove from com
            if (my.graph[com[1, i], com[2, i]] != 0 |
                my.graph[com[2, i], com[1, i]]!=0) {
                com[1, i] <- com[2, i] <- 0
                keep[i]<-F
            }
        }
        matrix(com[, keep],ncol=sum(keep))
    }
#
# find.possible.Q outputs a vector listing all other variables from 1:nvars
# except x and y
    find.possible.Q <- function(nvars, x, y) {
        z <- 1:nvars
        z[x] <- z[y] <- 0
        z[z > 0]
    }
#
# converts indices of variables to variable names; for dSep
dag.name<-function (amat,n) 
{
rownames(amat)[n]
}
#
full.vars<-row.names(full.DAG)
full.vars.index<-1:length(full.vars)
n.observed<-length(full.vars)-length(latents)
observed.DAG<-full.DAG
observed.vars<-full.vars
observed.vars.index<-full.vars.index
for(i in 1:length(latents)){
 observed.vars[latents[i]==full.vars]<-NA
 observed.vars.index[latents[i]==full.vars]<-NA
 observed.DAG[latents[i]==full.vars,]<-NA
 observed.DAG[,latents[i]==full.vars]<-NA
}
cat("the original DAG is:","\n")
total.n.vars<-dim(full.DAG)[2]
for(i in 1:(total.n.vars-1)){
 for(j in (i+1):total.n.vars){
  if(full.DAG[i,j]==1 & full.DAG[j,i]==0)cat(full.vars[i],"->",full.vars[j],"\n")
  if(full.DAG[i,j]==0 & full.DAG[j,i]==1)cat(full.vars[j],"->",full.vars[i],"\n")
 }
}
if(sum(is.na(latents))>0){
 return(cat("There are no latents; the DAG doesn't change ","\n"))
}
if(sum(is.na(latents))==0){
 cat("latent variable(s): ",latents,"\n")
 n.latents<-length(latents)
 for(i in 1:n.latents){
  ok<-F
  for(j in 1:length(full.vars))if(latents[i]==full.vars[j])ok<-T
  if(!ok)return("ERROR: latent variable name not in the DAG")
 }
}
cat("_____________________","\n")
observed.vars<-observed.vars[!is.na(observed.vars)]
observed.vars.index<-observed.vars.index[!is.na(observed.vars.index)]
#cat("full.vars:",full.vars,"\n")
#cat("full.vars.index:",full.vars.index,"\n")
#cat("observed.vars:",observed.vars,"\n")
#cat("observed.vars.index:",observed.vars.index,"\n","\n")
#
# construct initial observed DAG by removing latents and conserving directed
# edges between pairs of observed variables
#
# HERE
if(n.observed<=0)return(cat("No observed variables","\n"))
if(n.observed==1)return(cat("Only one observed variable","\n"))
if(n.observed==2)return(cat("Only two observed variables","\n"))

#cat("observed.DAG","\n")
observed.DAG<-observed.DAG[observed.vars.index,observed.vars.index]
#HERE
#for(i in 1:n.observed)cat(observed.DAG[i,],"\n")
#
# now find those observed variable pairs that don't share an edge in
# the full DAG.
#
if(n.observed<=0){
 return(cat("All variables are latent; there is no equivalent observed DAG","\n"))
}
#
pairs.to.test<-pairs.without.edge(observed.DAG)
#print(pairs.to.test)
#HERE
#cat("pairs.to.test; i.e. without edge obs index not full index:","\n")
n.pairs.to.test<-dim(pairs.to.test)[2]
#print(n.pairs.to.test)
#HERE
#for(i in 1:n.pairs.to.test)cat(pairs.to.test[,i],"\n")
n.remaining<-length(observed.vars)-2
# if all observed variables share an edge then return...
#HERE IS A CHANGE
#if(n.remaining<=0){
if(n.pairs.to.test<=0){
 return(cat("Since there are only two observed variables, nothing further will be done","\n"))
}
add.edge<-matrix(NA,nrow=2,ncol=n.pairs.to.test)
# for each pair (i) to test, determine dsep in full graph given only the observed
#
kount<-0
# i cycles over each pair that are not adjacent...
for(i in 1:n.pairs.to.test){
 is.pair.dsep<-F
# get those other observed variables in graph except this pair...
 possible.Q<-find.possible.Q(n.observed,pairs.to.test[1,i],pairs.to.test[2,i])
#cat("i=",i,"\n")
# Do do unconditional dseparation...
# i.e. conditional order=0
 first.var<-observed.vars.index[pairs.to.test[1,i]]
 second.var<-observed.vars.index[pairs.to.test[2,i]]
 test<-dSep(amat=full.DAG,first=dag.name(full.DAG,first.var),second=dag.name(full.DAG,second.var),cond=NULL)
# if first.var is dsep from second.var then there is no edge between them;
 if(test){
  is.pair.dsep<-T
#cat("pair ",dag.name(full.DAG,first.var),dag.name(full.DAG,second.var),"are unconditionally dsep","\n")
#  break
   next
 }
# if here then there are potential conditional variables to consider
# so cycle through all possible conditional orders...
 if(sum(is.na(possible.Q)==0)){
  n.possible.Q<-length(possible.Q)
#cat("for pair ",i,"possible.Q is:",possible.Q,"\n")
#cat("index for these vars is",observed.vars.index[possible.Q],"\n")
#
#now, determine, using observed.vars.index[possible.Q], if the pair are dsep
# in the full graph
# j gives the conditional order for a given pair
  for(j in 1:n.possible.Q){
#cat("conditional order is: ",j,"\n")
# Q has column = different combinations and rows=elements in each combination
   Q<-combn(possible.Q,j)
#cat("before..., Q= ",Q,"\n")
   if(j==n.possible.Q)Q<-matrix(possible.Q,nrow=j,ncol=1)
#cat("i=",i,"j=",j," Q= ",Q,"\n")
   n.Q<-dim(Q)[2]
#cat("There are ",n.Q," combinations at this conditional order","\n")
   first.var<-observed.vars.index[pairs.to.test[1,i]]
#cat("first.var=",first.var,"i.e. pairs.to.test[1,i]=",
#   pairs.to.test[1,i],"dag name=",dag.name(full.DAG,first.var),"\n")
   second.var<-observed.vars.index[pairs.to.test[2,i]]
#cat("second.var=",second.var,"i.e. pairs.to.test[2,i]=",
#    pairs.to.test[2,i],"dag.name=",dag.name(full.DAG,second.var),"\n")
# k cycles through these different combinations
   for(k in 1:n.Q){
    cond.vars<-as.vector(observed.vars.index[Q[,k]])
#cat("cond.vars=",cond.vars,"dag name(s)=",dag.name(full.DAG,cond.vars),"\n")
    test<-dSep(amat=full.DAG,first=dag.name(full.DAG,first.var),second=dag.name(full.DAG,second.var),
     cond=dag.name(full.DAG,cond.vars))
#cat("dSep?",test,"\n")
# if first.var dsep from second.var then there is no edge...
    if(test){
     is.pair.dsep<-T
     break
    }
   }
  }
 }
#cat(" for pair (",pairs.to.test[1,i],pairs.to.test[2,i],"), is.pair.dsep=",is.pair.dsep,"\n")
 if(!is.pair.dsep){
  kount<-kount+1
  add.edge[1,kount]<-pairs.to.test[1,i]
  add.edge[2,kount]<-pairs.to.test[2,i]
 }
}
# convert observed DAG to a partially oriented graph
cgraph<-matrix(0,n.observed,n.observed,dimnames=list(observed.vars,observed.vars))
for(i in 1:(n.observed-1)){
 for(j in (i+1):n.observed){
  if(observed.DAG[i,j]==1 & observed.DAG[j,i]==0){
   cgraph[i,j]<-2
   cgraph[j,i]<-1
  }
  if(observed.DAG[j,i]==1 & observed.DAG[i,j]==0){
   cgraph[j,i]<-2
   cgraph[i,j]<-1
  } 
 }
}
for(i in 1:kount){
 cgraph[add.edge[1,i],add.edge[2,i]]<-cgraph[add.edge[2,i],add.edge[1,i]]<-1
} 
cat("Equivalent partially oriented graph involving only the observed variables:","\n")
ind.vars<-rep(T,n.observed)
for(i in 1:(n.observed-1)){
 for(j in (i+1):n.observed){
  if(cgraph[i,j]==2 & cgraph[j,i]==1)cat(observed.vars[i],"->",observed.vars[j],"\n")
  if(cgraph[i,j]==1 & cgraph[j,i]==2)cat(observed.vars[j],"->",observed.vars[i],"\n")
  if(cgraph[i,j]==1 & cgraph[j,i]==1)cat(observed.vars[i],"--",observed.vars[j],"\n")
  if(cgraph[i,j]>0)ind.vars[i]<-F
  if(cgraph[j,i]>0)ind.vars[j]<-F
 }
}
for(i in 1:n.observed)if(ind.vars[i])cat(observed.vars[i],"--",observed.vars[i],"\n")

#HERE
#cat("before final.dag","\n")
if(n.observed<3)return()
final.dag<-orient.graph(cgraph,n.observed)
#HERE
#cat("after final.dag","\n")
# convert back to binary matrix and test for cycles...
new.dag<-matrix(0,n.observed,n.observed)
for(i in 1:(n.observed-1)){
 for(j in (i+1):n.observed){
  if(final.dag[i,j]==2 & final.dag[j,i]==1)new.dag[i,j]<-1
  if(final.dag[i,j]==1 & final.dag[j,i]==2)new.dag[j,1]<-1
 }
}
test.acyclic<-isAcyclic(new.dag)
if(!test.acyclic){
 return(cat("No possible DAG; orientation requires at least one latent","\n"))
}

if(test.acyclic)cat("One *possible* equivalent DAG involving only the observed variables:","\n")
for(i in 1:(n.observed-1)){
 for(j in (i+1):n.observed){
  if(final.dag[i,j]==2 & final.dag[j,i]==1)cat(observed.vars[i],"->",observed.vars[j],"\n")
  if(final.dag[i,j]==1 & final.dag[j,i]==2)cat(observed.vars[j],"->",observed.vars[i],"\n")
  if(final.dag[i,j]==1 & final.dag[j,i]==1)cat(observed.vars[i],"--",observed.vars[j],"\n")
 }
}
for(i in 1:n.observed)if(ind.vars[i])cat(observed.vars[i],"--",observed.vars[i],"\n")
}
