# TODO
# change code to take in list of all nodes with scores, rather than just top x
# color nodes according to scores
# make nodes equidistant
# add relationships edges (part_of, etc.)
# cluster by cellular component
# add documentation

require(GO.db, quietly = TRUE)
require(DataCombine, quietly = TRUE)
require(igraph, quietly = TRUE)
require(DiagrammeR, quietly = TRUE)
require(searchable, quietly = TRUE)

go_ids<-readRDS("example_ids.rds")

getSubgraphNodes <- function(enrichment_results){
  nodes <- data.frame(id=character(), stringsAsFactors = FALSE)
  for(i in enrichment_results[["id"]]){
    nodes <- rbind(nodes, data.frame(id=i))
    ancestors <- as.list(GOBPANCESTOR[[i]])
    ancestors <- ancestors[!is.na(ancestors)]
    ancestors <- ancestors[which(ancestors != "all")]
    if(length(ancestors) > 0){
      for(ancestor in ancestors){
        newRow <- data.frame(id=ancestor, stringsAsFactors = FALSE)
        nodes <- rbind(nodes, newRow)
      }
    }
  }
  return(unique(nodes))
}

buildNodeDF <- function(enrichment_results, top){
  n <- data.frame(nodes=character(), type=character(),
                  label=character(), style=character(),
                  color=character(), shape=character(),
                  data=numeric(), stringsAsFactors = FALSE)
  nodeList <- getSubgraphNodes(enrichment_results[1:top, ])
  names(nodeList) <- c("nodes")
  for(i in 1:dim(nodeList)){
    if(as.character(nodeList[i, ]) %in% enrichment_results[["id"]]){
      colorIndex <- (as.numeric(enrichment_results[enrichment_results$id == as.character(nodeList[i, ]), "score"])) * 100 + 1
      nodeColor <- heat.colors(100)[colorIndex]
      newNode <- data.frame(id=as.numeric(substr(as.character(nodeList[i, ]), 4, nchar(as.character(nodeList[i, ])))),
                 type="normal",
                 label= gsub(" ", "\n", Term(GOTERM[[as.character(nodeList[i, ])]])), style="filled",
                 fillcolor=nodeColor,
                 shape="circle", data=0,
                 fontsize = 40, fontcolor = "black",
                 stringsAsFactors = FALSE)
    }
    else{
      newNode <- data.frame(id=as.numeric(substr(as.character(nodeList[i, ]), 4, nchar(as.character(nodeList[i, ])))),
                            type="normal",
                            label= gsub(" ", "\n", Term(GOTERM[[as.character(nodeList[i, ])]])), style="filled", fillcolor="white",
                            shape="circle", data=0,
                            fontsize = 40, fontcolor = "black",
                            stringsAsFactors = FALSE)
    }
    n <- rbind(n, newNode)
  }
  #nodes <- rbind(nodes, nodeList)
  # TODO the above code is losing the empty columns, need to add default values for them
  return(n)
  # add node types, colour, etc here
}



buildEdgeDF <- function(nodeDF){
  edges <- data.frame(from=character(), to=character(), rel=character(), stringsAsFactors = FALSE)
  terms <- invert(Term(GOTERM))
  listOfNodes <- character()
  for(i in nodeDF[["label"]]){
    # build list of node GO IDs in the nodeDF here, to check with children
    listOfNodes <- c(listOfNodes, terms[[gsub("\n", " ", i)]])
  }
  for(i in nodeDF[["label"]]){
    go_id <- terms[[gsub("\n", " ", i)]]
    children <- as.list(GOBPCHILDREN[[go_id]])
    children <- children[!is.na(children)]
    if(length(children) > 0){
      in_common <- unique(intersect(children, listOfNodes))
      for(j in in_common){
        edges <- rbind(edges, data.frame(from=as.numeric(substr(j, 4, nchar(j))),
                                         to=as.numeric(substr(go_id, 4, nchar(go_id))),
                                         rel="leading_to"))
      }
    }
  }
  return(edges)
}

# need to pass in full go_ids with scores, and an integer (top x to use for subgraph generation)
nodes <- buildNodeDF(go_ids, 10)
edges <- buildEdgeDF(nodes)
gr <- create_graph(nodes_df = nodes,
                   edges_df = edges,
                   directed=TRUE)
gr <- add_global_graph_attrs(gr, attr="layout", value="twopi", attr_type = "graph")
gr <- add_global_graph_attrs(gr, attr="root", value="8150", attr_type = "graph")
gr <- add_global_graph_attrs(gr, attr="ranksep", value="10", attr_type = "graph")
gr <- add_global_graph_attrs(gr, attr="nodesep", value="1.2", attr_type = "graph")
gr <- add_global_graph_attrs(gr, attr="rank", value="max", attr_type = "graph")
gr <- add_global_graph_attrs(gr, attr="smoothing", value="true", attr_type = "graph")
gr <- add_global_graph_attrs(gr, attr="splines", value="curved", attr_type = "graph")
gr <- add_global_graph_attrs(gr, attr="forcelabels", value="true", attr_type = "graph")
gr <- add_global_graph_attrs(gr, attr="fixedsize", value="true", attr_type = "node")
gr <- add_global_graph_attrs(gr, attr="width", value="5", attr_type = "node")
gr <- add_global_graph_attrs(gr, attr="height", value="5", attr_type = "node")


render_graph(gr)
