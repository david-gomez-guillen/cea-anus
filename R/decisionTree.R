library(jsonlite)
library(ggplot2)
library(profvis)
library(dirmult)
library(yaml)
library(profvis)
library(visNetwork)
library(xlsx)

Node <- setRefClass('Node',
                    fields=c('id', 'name', 'info', 'probs', 'out'),
                    methods=list(
                      parseNodeValue = function(value, context=list()) {
                        if (length(value) > 0 && 
                            !is.null(value) && !is.na(value) &&
                            tryCatch(is.na(as.numeric(value)), warning=function(c){TRUE})) {
                          if (! value %in% names(context) || is.null(context[[value]])) {
                            stop(paste0('"', value, '" is not defined in context'))
                          }
                          retValue <- context[[value]]
                        } else {retValue <- value}
                        return(as.numeric(retValue))
                      },
                      parseProbs = function(context=list()) {
                        actualOutProbs <- probs
                        if (length(actualOutProbs) == 1) {
                          # Reference to distribution sample
                          actualOutProbs <- parseNodeValue(actualOutProbs, context)
                        } else {
                          actualOutProbs <- sapply(actualOutProbs, function(p) {
                            if (p == '#') return(p)
                            else parseNodeValue(p, context)
                          }, USE.NAMES = FALSE)
                        }
                        
                        if (any(actualOutProbs == '#')) {
                          actualOutProbs[actualOutProbs == '#'] = 1 - sum(as.numeric(actualOutProbs[actualOutProbs != '#']))
                          actualOutProbs <- as.numeric(actualOutProbs)
                        } 
                        
                        EPSILON <- .01
                        if (abs(1 - sum(actualOutProbs)) > EPSILON) {
                          stop(paste0('Probabilities from node "', name, '" don\'t add up to one'))
                        } 
                        
                        return(actualOutProbs)
                      },
                      getExpectedMeasure = function(measure, context=list(), use.only.leaves=FALSE) {
                        # If accumulative = TRUE, the expected value of the measure will be the weighted mean
                        # of all nodes. If false, the weighted mean will only consider the leaves.
                        measureValue <- info[[measure]]
                        m <- parseNodeValue(measureValue, context)
                        if ((use.only.leaves && length(out) > 0) || length(m) == 0 || is.na(m)) {
                          # warning(paste0('"', measureValue, '" does not exist for node "', name, '"; assuming 0'))
                          currentMeasure <- 0
                        }
                        else {
                          currentMeasure <- m
                        }
                        
                        if (length(out) > 0) {
                          actualOutProbs <- parseProbs(context)
                          outMeasure <- sum(sapply(seq_along(out), 
                                                        function(i) {
                                                          child <- out[[i]]
                                                          p <- as.numeric(actualOutProbs[[i]])
                                                          return(p * child$getExpectedMeasure(measure, context))
                                                        }
                          )
                          )
                          currentMeasure <- currentMeasure + outMeasure
                        }
                        return(currentMeasure)
                      },
                      getExpectedMeasureBatch = function(measure, context=list(), use.only.leaves=FALSE) {
                        if (length(context) == 0) {
                          num_simulations <- 1
                        } else if (is.matrix(context[[1]])) {
                          num_simulations <- nrow(context[[1]])
                        } else {
                          num_simulations <- length(context[[1]])
                        }
                        num_samples <- min(sapply(context, function(v) length(v)))
                        results <- sapply(seq(num_simulations), function(i) {
                          iterationContext <- sapply(context, function(v) {
                            if (length(v) != num_samples) {
                              v
                            } else if (is.matrix(v)) {
                              v[i,]
                            } else {
                              v[[i]]
                            }
                          })
                          iterResults <- getExpectedMeasure(measure, iterationContext, use.only.leaves)
                          return(iterResults)
                        })
                        return(results)
                      },
                      getExpectedCost = function(context=list()) {
                        return(getExpectedMeasureBatch('cost', context, use.only.leaves = FALSE))
                      },
                      getExpectedEffectiveness = function(context=list()) {
                        return(getExpectedMeasureBatch('utility', context, use.only.leaves = TRUE))
                      },
                      getCostEffectiveness = function(context=list()) {
                        return(getExpectedCost(context) / getExpectedEffectiveness(context))
                      },
                      getNetMonetaryBenefit = function(wtp, context=list()) {
                        return(getExpectedEffectiveness(context) * wtp - getExpectedCost(context))
                      },
                      getNetHealthBenefit = function(wtp, context=list()) {
                        return(getExpectedEffectiveness(context) - getExpectedCost(context) / wtp)
                      },
                      getUsedVariables = function() {
                        vars <- c(probs)
                        vars <- append(vars, c(info[['cost']], info[['utility']]))
                        vars <- vars[vars != '#' & !is.na(vars) & suppressWarnings(is.na(as.numeric(vars)))]
                        for(node in out) {
                          vars <- c(vars, node$getUsedVariables())
                        }
                        return(vars[!duplicated(vars)])
                      },
                      runMarkovSimulation = function(context=list(), max.steps=100) {
                        if (length(context) == 0) {
                          num_simulations <- 1
                        } else if (is.matrix(context[[1]])) {
                          num_simulations <- nrow(context[[1]])
                        } else {
                          num_simulations <- length(context[[1]])
                        }
                        results <- lapply(seq(num_simulations), function(i) {
                          iterationContext <- sapply(context, function(v) {
                            if (is.matrix(v)) {
                              v[i,]
                            } else {
                              v[[i]]
                            }
                          })
                          iterResults <- as.environment(.runSingleMarkovSimulation(context=iterationContext, max.steps=max.steps))
                          return(iterResults)
                        })
                        if (!is.null(attribute)) {
                          results <- sapply(results, function(r) get(attribute, envir=r))
                        }
                        return(results)
                      },
                      .runSingleMarkovSimulation = function(context=list(), max.steps=100, step=0) {
                        if (length(out) == 0 || 
                            step == max.steps || 
                            any(probs==1)) {
                          return(.self)
                        }
                        else {
                          pbs <- parseProbs(context)
                          randomChildIndex <- sample(length(pbs), size=1, prob=pbs)
                          return(out[[randomChildIndex]]$.runSingleMarkovSimulation(context, step=step+1, max.steps=max.steps))
                        }
                      },
                      runSimulation = function(context=list(), max.steps=100, attribute=NULL) {
                        if (length(context) == 0) {
                          num_simulations <- 1
                        } else if (is.matrix(context[[1]])) {
                          num_simulations <- nrow(context[[1]])
                        } else {
                          num_simulations <- length(context[[1]])
                        }
                        results <- lapply(seq(num_simulations), function(i) {
                          iterationContext <- sapply(context, function(v) {
                            if (is.matrix(v)) {
                              v[i,]
                            } else {
                              v[[i]]
                            }
                          })
                          iterResults <- as.environment(.runSingleSimulation(context=iterationContext, max.steps=max.steps))
                          return(iterResults)
                        })
                        if (!is.null(attribute)) {
                          results <- sapply(results, function(r) get(attribute, envir=r))
                        }
                        return(results)
                      },
                      .runSingleSimulation = function(context=list(), max.steps=100, step=0) {
                        if (length(out) == 0 || 
                            step == max.steps || 
                            any(probs==1)) {
                          return(.self)
                        }
                        else {
                          pbs <- parseProbs(context)
                          randomChildIndex <- sample(length(pbs), size=1, prob=pbs)
                          return(out[[randomChildIndex]]$.runSingleSimulation(context, step=step+1, max.steps=max.steps))
                        }
                      },
                      toString = function(level = 0, prob = NULL) {
                        displayedTree <- paste0(paste0(rep('--',level), collapse=''), name)
                        text <- c()
                        if (!is.null(prob)) {
                          text <- c(text, paste0('p = ', prob))
                        }
                        if (!is.null(info[['cost']])) {
                          text <- c(text, paste0('cost = ', info[['cost']]))
                        }
                        if (!is.null(info[['utility']])) {
                          text <- c(text, paste0('u = ', info[['utility']]))
                        }
                        if (length(text) > 0) {
                          displayedTree <- paste0(displayedTree, ' (', paste0(text, collapse = ', '), ')')
                        }
                        for(i in seq_along(out)) {
                          child <- out[[i]]
                          if (length(probs) > 1) {
                            childProb <- probs[[i]]
                          } else {
                            childProb <- paste0(probs, '[', i, ']')
                          }
                          displayedTree <- paste(displayedTree, child$toString(level=level+1, prob=childProb), sep='\n')
                        }
                        return(displayedTree)
                      },
                      getNodes = function() {
                        nodes <- .getNodes()
                        for(node in nodes) {
                          node$info['visited'] <- NULL
                        }
                        return(nodes)
                      },
                      .getNodes = function() {
                        if (!is.null(.self$info[['visited']])) {
                          return(list())
                        }
                        .self$info['visited'] <- TRUE
                        nodes <- .self
                        for(child in out) {
                          nodes <- append(nodes, child$.getNodes())
                        }
                        return(nodes)
                      },
                      getEdgeDisplayInfo = function(context=NULL) {
                        nodes <- getNodes()
                        for(node in nodes)
                          node$info['visited'] <- NULL
                        
                        edges <- .getEdgeDisplayInfo(context)
                        for(node in nodes)
                          node$info['visited'] <- NULL
                        return(edges)
                      },
                      .getEdgeDisplayInfo = function(context=NULL) {
                        if (!is.null(.self$info[['visited']]))
                          return(list())
                        .self$info['visited'] <- TRUE
                        edges <- data.frame()
                        if (!is.null(context)) {
                          actualProbs <- parseProbs(context)
                        } else {
                          actualProbs <- probs
                        }
                        for(i in seq_along(out)) {
                          child <- out[[i]]
                          label <- ''
                          if (!is.null(child$info[['in_transition']])) {
                            label <- paste(strwrap(child$info['in_transition'], 20), '\n', collapse='\n')
                          }
                          if (length(actualProbs) > 1) {
                            if (actualProbs[[i]] == 0) next  # If prob=0 we don't draw the edge
                            label <- paste0(label, '[', actualProbs[[i]], ']')
                          } else {
                            label <- paste0(label, '[', actualProbs, '(', i, ')', ']')
                          }
                          edges <- rbind(edges, data.frame(from=id, to=child$id, label=label, stringsAsFactors = F))
                          edges <- rbind(edges, child$.getEdgeDisplayInfo(context))
                        }
                        return(edges)
                      },
                      show = function(context=NULL, spacing=400) {
                        nodes <- getNodes()
                        names <- sapply(nodes, function(n) n$name)
                        ids <- sapply(nodes, function(n) n$id)
                        nodes <- data.frame(id=ids, label=names, stringsAsFactors = F)
                        edges <- getEdgeDisplayInfo(context)
                        if (info[['model']] == 'tree') {
                          nodes$shape <- 'box'
                          p <- visNetwork(nodes, edges) %>% 
                            visEdges(arrows='to', length=200) %>% 
                            visIgraphLayout() %>%
                            visHierarchicalLayout(sortMethod = 'directed', nodeSpacing = spacing) 
                        } else if (info[['model']] == 'markov') {
                          nodes$shape <- 'circle'
                          p <- visNetwork(nodes, edges) %>% 
                            visEdges(arrows='to', length=300)
                        }
                        print(p)
                      },
                      compareStrategies = function(..., context=list()) {
                        if (length(context) > 0) {
                          samples <- length(context[[1]])
                          ls <- lapply(context, function(v) {
                            if (!is.null(nrow(v))) {
                              return(nrow(v))
                            } else {
                              length(v)
                            }
                          })
                          if (any(ls != samples)) {
                            stop('Context variables are not all equal in length')
                          }
                        } else {
                          samples <- 1
                        }
                        
                        fullResults <- data.frame()
                        num_samples <- sapply(context, function(v) length(v))
                        for(i in seq(samples)) {
                          currentContext <- sapply(context, function(v) {
                            if (is.matrix(v) || (length(v) != num_samples)) {
                              return(v[i,])
                            } else {
                              return(v[[i]])
                            }
                          })
                          
                          baseCost <- getExpectedCost(currentContext)
                          baseEff <- getExpectedEffectiveness(currentContext)
                          baseCE <- getCostEffectiveness(currentContext)
                          
                          displayedContext <- lapply(currentContext, function(e) {if (length(e) == 1) e else {paste(e, collapse = ', ')}})
                          
                          results <- data.frame(
                            append(
                              displayedContext,
                              list(
                                strategy='base',
                                C=baseCost,
                                IC=NA,
                                E=baseEff,
                                IE=NA,
                                CE=baseCE,
                                ICER=NA
                              )))
                          alternatives <- list(...)
                          for (alt in names(alternatives[[1]])) {
                            altTree <- alternatives[[1]][[alt]]
                            c <- altTree$getExpectedCost(currentContext)
                            ic <- c - baseCost
                            e <- altTree$getExpectedEffectiveness(currentContext)
                            ie <- e - baseEff
                            ce <- altTree$getCostEffectiveness(currentContext)
                            if (ie != 0) {
                              icer <- ic / ie
                            }
                            else {
                              icer <- NA
                            }
                            alternativeResults <- data.frame(
                              append(
                                displayedContext,
                                list(
                                  strategy=alt,
                                  C=c,
                                  IC=ic,
                                  E=e,
                                  IE=ie,
                                  CE=ce,
                                  ICER=icer
                                )))
                            results <- rbind(results, alternativeResults)
                          }
                          fullResults <- rbind(fullResults, results)
                        }
                        
                        fullSummary <- fullResults[,c('strategy', 'C', 'IC', 'E', 'IE', 'CE', 'ICER')]
                        fullSummary.mean <- aggregate(fullSummary, by=list(group=fullSummary$strategy), function(g) {
                          if (!any(is.numeric(g))) NA else mean(g)
                          })
                        fullSummary.mean$strategy <- fullSummary.mean$group
                        fullSummary.sd <- aggregate(fullSummary, by=list(group=fullSummary$strategy), function(g) {
                          if (!any(is.numeric(g))) NA else sd(g)
                        })
                        fullSummary.sd$strategy <- fullSummary.sd$group
                        
                        fullSummary.mean <- fullSummary.mean[order(fullSummary.mean$E),]
                        fullSummary.mean <- rbind(
                          fullSummary.mean[fullSummary.mean$strategy=='base',],
                          fullSummary.mean[fullSummary.mean$strategy!='base',]
                        )
                        fullSummary.mean$domination <- 'undominated'
                        for(index in seq(2,nrow(fullSummary.mean))) {
                          if (fullSummary.mean[index, 'IE'] < 0) {
                            fullSummary.mean[index, 'domination'] <- 'absolute'
                          } else if (any(fullSummary.mean[index, 'ICER'] > fullSummary.mean[min(index+1, nrow(fullSummary.mean)):nrow(fullSummary.mean), 'ICER'])) {
                            fullSummary.mean[index, 'domination'] <- 'extended'
                          }
                        }
                        
                        plot.df <- fullSummary.mean
                        plt <- ggplot(plot.df, aes(x=E, y=C)) +
                          geom_line(data=plot.df[plot.df$domination=='undominated',]) +
                          geom_point(data=plot.df[plot.df$domination=='absolute',], color='red', size=5) +
                          geom_point(data=plot.df[plot.df$domination=='extended',], color='blue', size=5) +
                          geom_point(size=3, aes(shape=strategy)) +
                          xlab('Effectiveness') +
                          ylab('Cost')
                        print(plt)
                        
                        return(fullSummary.mean)
                      }
                    ))

parseYAML <- function(filePath, nextId=0) {
  treeDir <- dirname(filePath)
  treeSpec <- yaml.load_file(filePath)
  nodes <- list()
  # childrenNodes <- c()
  
  parseProbs <- function(probStr) {
    pbs <- strsplit(probStr, ',')[[1]]
    pbs <- trimws(pbs)
    pbs <- ifelse(pbs=='_', '#', pbs)
    return(pbs)
  }
  
  parseNode <- function(node) {
    name <- names(node)
    node <- node[[1]]
    attributes <- names(node)
    attributes <- attributes[!attributes %in% c('name', 'children', 'probs', 'include')]
    
    newNode <- Node(id=nextId,
                    name=name,
                    info=list(),
                    probs=node$probs)
    nextId <<- nextId + 1
    
    for(val in attributes) {
      newNode$info[val] <- node[val]
    }
    newNode$info['model'] <- 'tree'
    
    if (length(node$children) > 0) {
      newNode$probs <- parseProbs(node$probs)
      newNode$out <- lapply(node$children, function(child) {
        isIncluded <- 'include' %in% names(child[[1]]) && 
          !is.null(child[[1]][['include']]) && 
          !is.na(child[[1]][['include']])
        childNode <- parseNode(child)
        info <- childNode$info
        if (isIncluded) {
          childNode <- parseYAML(paste0(treeDir, '/', child[[1]]['include'][[1]]), nextId=nextId)
        }
        childNode$info <- modifyList(childNode$info, info)
        childNode
      })
      names(newNode$out) <- lapply(newNode$out, function(child) {child$name})
    } else {
      newNode$probs <- numeric(0)
      newNode$out <- list()
    }
    newNode
  }
  
  tree <- parseNode(treeSpec)
  
  tree
}

parseExcel <- function(filePath) {
  wb <- loadWorkbook(filePath)
  sheetnames <- names(getSheets(wb))
  
  parseProbs <- function(probStr) {
    pbs <- strsplit(probStr, ',')[[1]]
    pbs <- trimws(pbs)
    pbs <- ifelse(pbs=='_', '#', pbs)
    return(pbs)
  }
  
  nodeInfo <- read.xlsx(filePath, sheetName = 'nodes', stringsAsFactors=F)
  
  edge_sheets <- sheetnames[startsWith(sheetnames, 'edges')]
  for (net in edge_sheets) {
    edges <- read.xlsx(filePath, sheetName = net, stringsAsFactors=F)
    if (any(names(edges)[2:length(names(edges))] != edges[[1]]))
      stop('Row and column states are not equal')
    row.names(edges) <- unlist(edges[1])
    edges[1] <- NULL
    
    nodes <- list()
    for(i in row.names(nodeInfo)) {
      node <- nodeInfo[i,]
      newNode <- Node(id=i,
                      name=node$name,
                      info=list(),
                      out=list(),
                      probs=node$probs)
      newNode$info['model'] <- 'markov'
      nodes <- append(nodes, newNode)
    }
    
    for(i in seq_along(nodes)) {
      for(j in seq_along(nodes)) {
        nodes[[i]]$probs[j] <- edges[i,j]
        nodes[[i]]$out <- append(nodes[[i]]$out, nodes[[j]])
      }
    }
    return(nodes[[1]])
  }
}