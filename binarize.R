# R Code for Calculating Network Measures
#Calculation of Brainerd-Robinson similarity coefficients:
#  Input: frequency table sites (or other units of analysis) as rows and type categories as columns.
#Example format:
#  Type 1	Type 2	Type 3
#Site 1	551	331	25
#Site 2	10	300	62
#Site 3	60	112	0
#Site 4	12	10	14

sim.calc <- function(x) { # define function to calculate similarity scores
  x <- prop.table(as.matrix(x),1)*100 # convert frequency table into row % matrix
  rd <- dim(x)[1] # determine number of rows (sites) in table
  results <- matrix(0,rd,rd) # create matrix with rows and columns for each site
  for (s1 in 1:rd) { # open loop for s1 in 1 to the number of rows
    for (s2 in 1:rd) { # open loop for s2 in 1 to the number of rows
      x1Temp <- as.numeric(x[s1, ]) # lookup values in matrix x at row s1
      x2Temp <- as.numeric(x[s2, ]) # lookup values in matrix x at row s2
      results[s1,s2] <- 200 - (sum(abs(x1Temp - x2Temp)))}} # similarity of row s1 to s2
  row.names(results) <- row.names(x) # add site names to row labels
  colnames(results) <- row.names(x) # add site names to column labels
  results <- results/200 # rescale similarity to range from 0 to 1
  return(results)} # return the result of this function

sim.mat <- sim.calc(mydata) # run sim.calc function on frequency table (mydata)

#Example sim.mat from sim.calc function using sample table above:
#  Site 1	Site 2	Site 3	Site 4
#Site 1	1	0.42	0.71	0.64
#Site 2	0.42	1	0.68	0.47
#Site 3	0.71	0.68	1	0.61
#Site 4	0.64	0.47	0.61	1

#Creation of Binary Network:
#  Requires R packages statnet and sna
#To install package, type the following at the R console. This also installs the sna package.

#install.packages('statnet')

#Input: similarity matrix produced in the previous step - sim.mat

library(statnet)# initialize statnet library
cutoff <- 0.5 # set threshold for creating a tie to 0.5 similarity
bin <- event2dichot(sim.mat,method='absolute',thresh=cutoff) # binarize sim.mat at cutoff
rownames(bin) <- rownames(sim.mat) # apply row names from sim.mat to bin
colnames(bin) <- colnames(sim.mat) # apply column names from sim.mat to bin
net <- network(bin,directed=F) # create undirected network from binarized matrix bin
plot(net,displaylabels=T) # show network graph with sites labeled


#The procedure above creates and displays a binary network from a similarity matrix (See above). To select a different binarization threshold, simply change the value for the cutoff variable. 

#Calculation of Centrality Scores for Binary Networks:
  
#  Input: R network produced in the previous step - net

net.stats <- function(y){ # define function to calculate centrality scores
  dg <- as.matrix(degree(y,gmode='graph')) # calculate degree centrality
  eg <- as.matrix(evcent(y)) # calculate eigenvector centrality
  eg <- sqrt((eg^2)*length(eg)) # scale so sum of squared scores = number of nodes
  bw <- betweenness(y,gmode='graph') # calculate betweenness centrality
  output <- cbind(dg,eg,bw) # combine centrality scores into matrix
  rownames(output) <- rownames(as.matrix(y)) # add row names to matrix
  colnames(output) <- c('dg','eg','bw') # add column names to matrix
  return (output)} # return results of this function

cent.scores <- net.stats(net) # run centrality scores for network (net)

#Example cent.scores from net.stats function using sample table above:
#  dg	eg	bw
#Site 1	2	1.05	0
#Site 2	1	0.56	0
#Site 3	3	1.22	2
#Site 4	2	1.05	0


#Calculation of Centrality Scores for Weighted Networks:
#  Requires R packages tnet and igraph
#To install packages, type the following at the R console (will install both packages).

#install.packages('tnet')

#Input: similarity matrix produced above - sim.mat

net.stats.wt <- function(y){ # define function to calculate weighted centrality scores
  dg.wt <- as.matrix(rowSums(y)-1) # calculate weighted degree centrality
  eg.wt <- as.matrix(evcent(y)) # calculate weighted eigenvector centrality
  eg.wt <- sqrt((eg.wt^2)*length(eg.wt)) # scale so sum of squared scores = number of nodes
  output <- cbind(dg.wt,eg.wt) # combine dg.wt and eg.wt centrality scores into matrix
  rownames(output) <- rownames(as.matrix(y)) # add row names to matrix
  colnames(output) <- c('dg','eg') # add column names to matrix
  return (output)} # return results of this function

temp.scores <- net.stats.wt(sim.mat) # run centrality scores for similarity matrix 

library(tnet) # initialize tnet library

sim.mat.t <- as.tnet(sim.mat) # convert similarity matrix to tnet object
bw.wt <- betweenness_w(sim.mat,directed=F,alpha=1) # calculate weighted betweenness
bw.wt <- as.matrix(bw.wt[,2]) # remove node labels
colnames(bw.wt) <- c('bw.wt') # add column name
cent.scores.wt <- cbind(temp.scores,bw.wt) # create matrix of weighted centrality scores

detach(package:tnet, unload=TRUE) # unload tnet package (required for other analyses)

#Example cent.scores.wt from net.stats.wt function using sample table above:
  
#  dg.wt	eg.wt	bw.wt
#Site 1	1.77	1.01	0
#Site 2	1.57	0.91	0
#Site 3	2.00	1.09	0
#Site 4	1.72	0.98	0


#Calculation of network centralization:
#  Requires R packages statnet

library(statnet) # initialize statnet library

## calculate centralization measures for binary network (net) created above
centralization(net,degree,normalize=T) #calculate binary degree centralization
centralization(net,betweenness,normalize=T) #calculate binary betweenness centralization
centralization(net,evcent,normalize=T) #calculate binary eigenvector centralization

## calculate centralization measures for weighted network (sim.mat) created above
centralization(sim.mat,degree,normalize=T) #calculate weighted degree centralization
centralization(sim.mat,evcent,normalize=T) #calculate weighted eigenvector centralization

bw.cent <- function(x){ # define function for calculating
  Cstar <- max(x) # determine maximum weighted betweenness centrality score in x
  Csum <- (Cstar-x) # create vector or maximum centrality score minus each individual score
  num <- 2*(sum(Csum)) # calculate numerator of centralization equation
  den <- ((length(x)-1)^2)*(length(x)-2) # calculate denominator of centralization equation
  out <- num/den # calculate centralization score
  return(out)} # output result of this function

# calculate betweenness centralization for bw.wt variable in cent.scores.wt (column 3)
bw.cent(cent.scores.wt[,3]) 
