#' TODO
#' Change color scale to better show the differences between small p-values
#' Find some way to remove node overlap while keeping the radial shape intact
#' Fix edge colors not showing (seems like a bug in DiagrammeR or Graphviz)
#' Add documentation
#' Fix naming conventions
#' Cluster nodes by cellular component localization

require(GO.db, quietly = TRUE)
require(DataCombine, quietly = TRUE)
require(DiagrammeR, quietly = TRUE)
require(searchable, quietly = TRUE)

go_ids<-readRDS("example_ids.rds")

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

buildLabel <- function(goTerm, scores){
  label <- paste(c(as.character(gsub(" ", "\n", Term(GOTERM[[goTerm]]))),
                                "\n\n", as.character(round(scores[[goTerm]], 4)),
                                "\n\n", goTerm),
                 collapse='')
  return(label)
}

buildNodeDF <- function(enrichment_results, top){
  n <- data.frame(nodes=character(), type=character(),
                  label=character(), style=character(),
                  color=character(), shape=character(),
                  data=numeric(), stringsAsFactors = FALSE)
  nodeList <- getSubgraphNodes(enrichment_results[1:top])
  for (i in seq(along=dim(nodeList)[1])){
    nodeLabel <- buildLabel(nodeList[i, ], go_ids)
    if (as.character(nodeList[i, ]) %in% names(enrichment_results)){
      colorIndex <- (enrichment_results[nodeList[i, ]] *
                     length(enrichment_results))
      nodeColor <- heat.colors(length(enrichment_results))[colorIndex]
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

rebuildGOID <- function(goID){
  return(paste(c("GO:", rep("0",(7-nchar(goID))), goID), collapse=''))
}

getEdgeColor <- function(edgeRelationships){
  edgeColor <- ""
  switch(edgeRelationships,
         "is_a" = edgeColor <- "blue",
         "part_of" = edgeColor <- "yellow",
         "has_part" = edgeColor <- "black",
         "regulates" = edgeColor <- "purple",
         "negatively_regulates" = color <- "green",
         "positively_regulates" = color <- "red")
  return(edgeColor)
}

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
        edgeStyle <- ""
        edgeColor <- ""
        lookUp <- searchable::invert(GOBPCHILDREN[[go_id]])
        edgeColor <- getEdgeColor(lookUp[[j]])
        edges <- rbind(edges, data.frame(from=as.numeric(substr(j, 4, nchar(j))),
                                         to=as.numeric(substr(go_id, 4,
                                                              nchar(go_id))),
                                         rel=lookUp[[j]],
                                         arrowsize=5,
                                         color=edgeColor))
      }
    }
  }
  return(edges)
}

generateGraph <- function(scores){
  nodes <- buildNodeDF(go_ids, 5)
  edges <- buildEdgeDF(nodes)
  gr <- create_graph(nodes_df = nodes,
                     edges_df = edges,
                     directed=TRUE)

  gr <- add_global_graph_attrs(gr, attr="layout", value="twopi",
                               attr_type = "graph")
  gr <- add_global_graph_attrs(gr, attr="root", value="8150",
                               attr_type = "graph")
  gr <- add_global_graph_attrs(gr, attr="ranksep", value="20",
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
  gr <- add_global_graph_attrs(gr, attr="overlap", value="prism0",
                               attr_type = "graph")
  gr <- add_global_graph_attrs(gr, attr="fixedsize", value="true",
                               attr_type = "node")
  gr <- add_global_graph_attrs(gr, attr="width", value="12",
                               attr_type = "node")
  gr <- add_global_graph_attrs(gr, attr="height", value="12",
                               attr_type = "node")
  gr <- add_global_graph_attrs(gr, attr="penwidth", value="8",
                               attr_type = "edge")

  render_graph(gr, width=1280, height=1024)
}
# [END]
