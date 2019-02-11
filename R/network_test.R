library(reshape2)
library(network)
library(igraph)
library(visNetwork)

df.orig <- data.frame(list(healthy=c(.75,.15,.1), cancer=c(0,.7,.3), death=c(0,0,1)), row.names=c('x', 'y', 'z'))

nodeIds <- seq(length(df.orig))
nodeNames <- names(df.orig)
nodes <- data.frame(list(id=nodeIds, label=nodeNames), stringsAsFactors = F)
nodes$test <- ifelse(nodes$id==1, 'red', 'blue')

dfm <- as.matrix(df.orig, rownames.force = F)
dimnames(dfm)[[2]] <- NULL
edges <- melt(dfm)
names(edges) <- c('from', 'to', 'weight')
  edges <- edges[edges$weight != 0,]

net <- network(edges, vertex.attr = nodes, matrix.type = 'edgelist', ignore.eval = F)
plot(net, vertex.cex = 3)

# net_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = T)
# plot(net_igraph, edge.arrow.size = 0.2, edge.label=edges$weight)

edges$label <- edges$weight
p <- visNetwork(nodes, edges) %>% 
  visEdges(arrows='from', length=200)
print(p)