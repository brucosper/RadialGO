#' TODO
#' Change color scale to better show the differences between small p-values
#' Fix edge colors not showing (seems like a bug in DiagrammeR or Graphviz)
#' Add documentation
#' Fix naming conventions
#' Cluster nodes by cellular component localization
#' Add edge attributes to new edges after reducing
#'
require(GO.db, quietly = TRUE)
require(DataCombine, quietly = TRUE)
require(DiagrammeR, quietly = TRUE)
require(searchable, quietly = TRUE)


#' Generate the subgraph implied by the list of nodes passed in.
#' @param enrichment_results A named numerical vector containing GO IDs as names
#'                           and the corresponding score
#' @return A data frame containing one column, named "id", with the nodes needed
#'         for building the graph
#' @examples
#' scores <- c(0.5, 0.2, 0.001)
#' names(scores) <- c("GO:0008150", "GO:1901360", "GO:0006139")
#' getSubgraphNodes(scores)
getSubgraphNodes <- function(enrichment_results){
  nodes <- data.frame(id=character(), stringsAsFactors = FALSE)
  for (i in names(enrichment_results)){
    nodes <- rbind(nodes, data.frame(id=i, stringsAsFactors = FALSE))
    ancestors <- as.list(GOBPANCESTOR[[i]])
    ancestors <- ancestors[!is.na(ancestors)]
    ancestors <- ancestors[which(ancestors != "all")]
    if (length(ancestors) > 0){
      for (ancestor in ancestors){
        newRow <- data.frame(id=ancestor, stringsAsFactors = FALSE)
        nodes <- rbind(nodes, newRow)
      }
    }
  }
  return(unique(nodes))
}

#' Create a text label describing the GO term and its score
#' @param goTerm A string containing the GO ID to use
#' @param scores A named numerical vector containing GO IDs as names
#'               and the corresponding score
#' @return A string containing the label text
#' @examples <- c(0.4)
#' scores <- c(0.5, 0.2, 0.001)
#' names(scores) <- c("GO:0008150", "GO:1901360", "GO:0006139")
#' buildLabel(names(scores[1]), scores)
buildLabel <- function(goTerm, scores){
  label <- paste(c(as.character(gsub(" ", "\n", Term(GOTERM[[goTerm]]))),
                                "\n\n", goTerm,
                                "\n\n", as.character(round(scores[[goTerm]], 4))),
                 collapse='')
  return(label)
}

#' Builds a node data frame from the top GO terms using the enrichment results
#' @param enrichment_results A named numerical vector containing GO IDs as names
#'                           and the corresponding score
#' @param top The number of top nodes to use to build the graph
#' @examples
#' scores <- c(0.5, 0.2, 0.001)
#' names(scores) <- c("GO:0008150", "GO:1901360", "GO:0006139")
#' buildNodeDF(scores, 2)
buildNodeDF <- function(enrichment_results, top){
  n <- data.frame(nodes=character(), type=character(),
                  label=character(), style=character(),
                  color=character(), shape=character(),
                  data=numeric(), stringsAsFactors = FALSE)
  nodeList <- getSubgraphNodes(enrichment_results[1:top])
  for (i in seq(along=1:dim(nodeList)[1])){
    nodeLabel <- buildLabel(nodeList[i, ], enrichment_results)
    if (as.character(nodeList[i, ]) %in% names(enrichment_results)){
      colorIndex <- (enrichment_results[nodeList[i, ]] *
                     length(enrichment_results))
      nodeColor <- colorRampPalette(c("red", "grey", "white"),
                                    space="rgb")(length(enrichment_results))[colorIndex]
      newNode <- data.frame(id=as.numeric(substr(nodeList[i, ], 4,
                                                 nchar(nodeList[i, ]))),
                            type="normal",
                            label= nodeLabel, style="filled",
                            fillcolor=nodeColor,
                            shape="circle", data=0,
                            fontsize = 80, fontcolor = "black",
                            stringsAsFactors = FALSE)
    } else {
      # as.character can be removed?
      newNode <- data.frame(id=as.numeric(substr(as.character(nodeList[i, ]),
                                                 4, nchar(nodeList[i, ]))),
                            type="normal",
                            label= nodeLabel, style="filled", fillcolor="white",
                            shape="circle", data=0,
                            fontsize = 80, fontcolor = "black",
                            stringsAsFactors = FALSE)
    }
    n <- rbind(n, newNode)
  }
  return(n)
}

#' Rebuild a GO ID string from its number
#' @param goID A string or number containing the GO ID number
#' @return The full GO ID
#' @examples
#'  rebuildGOID(8150)
#'
rebuildGOID <- function(goID){
  return(paste(c("GO:", rep("0",(7-nchar(goID))), goID), collapse=''))
}

#' Determine what color the edge should take in the graph.
#' @param edgeRelationships A string containing a relationship between GO terms
#' @return A string containing the color the edge should take
#' @examples
#' getEdgeColor("is_a")
#'
getEdgeColor <- function(edgeRelationships){
  edgeColor <- ""
  switch(edgeRelationships,
         "is_a" = edgeColor <- "blue",
         "part_of" = edgeColor <- "yellow",
         "has_part" = edgeColor <- "black",
         "regulates" = edgeColor <- "purple",
         "negatively_regulates" = edgeColor <- "green",
         "positively_regulates" = edgeColor <- "red")
  return(edgeColor)
}


#' Build a data frame of edges from the nodes. NOTE: Due to a bug
#' in the way DiagrammeR renders Graphviz objects, edges currently do not render
#' with the correct stylings. To see that they are specified correctly, one can
#' render the graph with the output="visNetwork" option, which shows correct
#' edges but loses the node layout
#' @param nodeDF A data frame containing at least columns id, label,
#' @return A data frame of edges (from, to, rel)
#' @examples
#' scores <- c(0.5, 0.2, 0.001)
#' names(scores) <- c("GO:0008150", "GO:1901360", "GO:0006139")
#' buildEdgeDF(buildNodeDF(scores, 2))
#'
buildEdgeDF <- function(nodeDF){
  edges <- data.frame(from=character(), to=character(), rel=character(),
                      style=character(), stringsAsFactors = FALSE)
  terms <- searchable::invert(Term(GOTERM))
  listOfNodes <- character()
  for (i in nodeDF[["id"]]){
    # build list of node GO IDs in the nodeDF here, to check with children
    go_id <- rebuildGOID(i)
    listOfNodes <- c(listOfNodes, go_id)
  }
  for (i in nodeDF[["id"]]){
    go_id <- rebuildGOID(i)
    children <- as.list(GOBPCHILDREN[[go_id]])
    children <- children[!is.na(children)]
    if (length(children) > 0){
      in_common <- unique(intersect(children, listOfNodes))
      for (j in in_common){
        edgeColor <- ""
        lookUp <- searchable::invert(GOBPCHILDREN[[go_id]])
        edgeColor <- getEdgeColor(lookUp[[j]])
        edgeLabel <- lookUp[[j]]
        edges <- rbind(edges, data.frame(from=as.numeric(substr(j, 4, nchar(j))),
                                         to=as.numeric(substr(go_id, 4,
                                                              nchar(go_id))),
                                         rel=lookUp[[j]],
                                         arrowsize=5,
                                         color=edgeColor,
                                         headlabel=edgeLabel))
      }
    }
  }
  return(edges)
}

#' Removes nodes from the graph which have a p-value above the specified cutoff,
#' collapsing the edges
#' @param graph The graph to reduce
#' @param cutoff The cutoff to be used (nodes above this cutoff will be removed)
#' @return A graph with the nodes removed
#'
reduceGraph <- function(graph, cutoff){
  edgeDF <- get_edge_df(graph)
  nodeDF <- get_node_df(graph)
  for(i in nodeDF[["id"]]){
    if(i == 8150){
      # GO:0008150 is the biological_process node (root)
      next
    } else {
      splitString <- strsplit(nodeDF[nodeDF$id == i, "label"], "\n")
      score <- as.numeric(splitString[[1]][length(splitString[[1]])])
      if(score > cutoff){
        # we're deleting this node
        edgesIn <- edgeDF[edgeDF$to == as.numeric(i), ]
        edgesOut <- edgeDF[edgeDF$from == as.numeric(i), ]
        if(nrow(edgesIn) != 0 && nrow(edgesOut) != 0){
          for(j in seq(along=1:nrow(edgesIn))){
            for(k in seq(along=1:nrow(edgesOut))){
              edgeAttrs <- edgesIn[edgesIn$from == edgesIn[j, ]$from, ]
              edgeAttrs <- edgeAttrs[edgeAttrs$to == edgesIn[j, ]$to, ]
              edgeAttrs <- edgeAttrs[c(5:ncol(edgeAttrs))]
              newEdge <- data.frame(from = edgesIn[j, ]$from,
                                    to = edgesOut[k, ]$to,
                                    rel=edgesOut[k, ]$rel,
                                    stringsAsFactors = FALSE)
              newEdge <- merge(newEdge, edgeAttrs)
              edgeDF <- combine_edfs(edgeDF, newEdge)
            }
          }
        }
        #remove old edges
        edgeDF <- edgeDF[edgeDF$to != as.numeric(i), ]
        edgeDF <- edgeDF[edgeDF$from != as.numeric(i), ]
        nodeDF <- nodeDF[nodeDF$id != i, ]
      }
    }
  }
  return(create_graph(nodes_df = nodeDF, edges_df = edgeDF))
}

#' Generate the graph visualization of GO enrichment scores
#' @param scores A named numeric vector with GO IDs and corresponding scores
#' @param top The number of nodes to use when generating the subgraph
#' @param cutoff The cutoff to be used for hiding nodes
#' @return A graph object to be used with render_graph()
#' @examples
#' scores <- c(0.5, 0.2, 0.001)
#' names(scores) <- c("GO:0008150", "GO:1901360", "GO:0006139")
#' generateGraph(scores, 2, 0.3)
generateGraph <- function(scores, top, cutoff){
  nodes <- buildNodeDF(scores, top)
  edges <- buildEdgeDF(nodes)


  gr <- create_graph(nodes_df = nodes, edges_df = edges, directed=TRUE)
  gr <- reduceGraph(gr, cutoff)
  gr <- add_global_graph_attrs(gr, attr="layout", value="twopi",
                               attr_type = "graph")
  gr <- add_global_graph_attrs(gr, attr="root", value="8150",
                               attr_type = "graph")
  gr <- add_global_graph_attrs(gr, attr="ranksep", value="15",
                               attr_type = "graph")
  gr <- add_global_graph_attrs(gr, attr="nodesep", value="20",
                               attr_type = "graph")
  gr <- add_global_graph_attrs(gr, attr="rank", value="same",
                               attr_type = "graph")
  gr <- add_global_graph_attrs(gr, attr="splines", value="curved",
                               attr_type = "graph")
  gr <- add_global_graph_attrs(gr, attr="forcelabels", value="true",
                               attr_type = "graph")
  gr <- add_global_graph_attrs(gr, attr="orientation", value="[lL]*",
                               attr_type = "graph")
  gr <- add_global_graph_attrs(gr, attr="overlap", value="true",
                               attr_type = "graph")
  gr <- add_global_graph_attrs(gr, attr="fixedsize", value="true",
                               attr_type = "node")
  gr <- add_global_graph_attrs(gr, attr="width", value="12",
                               attr_type = "node")
  gr <- add_global_graph_attrs(gr, attr="height", value="12",
                               attr_type = "node")
  gr <- add_global_graph_attrs(gr, attr="penwidth", value="8",
                               attr_type = "edge")

  return(gr)
}

go_ids<-readRDS("example_ids.rds")
graph <- generateGraph(go_ids, 5, 0.2)
render_graph(graph, width=1280, height=1024)
# [END]
