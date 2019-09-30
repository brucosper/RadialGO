# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

# want to pass to top function (enrichment analysis results with scores, number of top nodes to use when determining graph)


require(GO.db, quietly = TRUE)
require(DataCombine, quietly = TRUE)
require(igraph, quietly = TRUE)
require(DiagrammeR, quietly = TRUE)

go_ids<-readRDS("example_ids.rds")

# topGO enrichment analysis returns a list of GO IDs with a score (KS or Fischer)
# need to grab nodes, organize with edges (relationships) and plot with visNetwork
# TODO maybe use topGO's functions which already get the relationships
# see: https://rdrr.io/bioc/topGO/src/R/topGOviz.R

#get list of children for each node, so we can build edges


get_edges <- function(ids){
  edge_list<-data.frame(from=character(), to=character(), stringsAsFactors=FALSE)
  for(id in ids){
    children <- as.list(GOBPCHILDREN[[id]])
    children <- children[!is.na(children)]
    if(length(children) > 0){
      in_common <- unique(intersect(children, ids))
      for(i in in_common){
        newrow <- data.frame(from=i, to=id)
        edge_list <- InsertRow(edge_list, newrow)
      }
    }
  }
  return(edge_list)
}

get_nodes <- function(ids){
  nodes <- data.frame(id=character(), stringsAsFactors = FALSE)
  add_bp <- data.frame(id="GO:0008150") # add root BP term so it's the first one (root)
  nodes <- InsertRow(nodes, add_bp)
  for (i in ids){
      newrow <- data.frame(id=i, stringsAsFactors = FALSE)
      nodes <- InsertRow(nodes, newrow)
      ancestors <- as.list(GOBPANCESTOR[[i]])
      ancestors <- ancestors[!is.na(ancestors)]
      ancestors <- ancestors[which(ancestors != "all")]
      if(length(ancestors) > 0){
        for(ancestor in ancestors){
            newrow2 <- data.frame(id=ancestor, stringsAsFactors = FALSE)
            nodes <- InsertRow(nodes, newrow2)
        }
      }
  }
  return(unique(nodes))
}

buildGraphvizGraph <- function(nodes, edges){
  out <- "digraph Enrichment {
          graph [layout=twopi, overlap = true, splines=curved, fontsize = 1, nodesep=10, ranksep=4]
          node [shape=circle, fillcolor = orange, style = filled, width=2]
          "
  for(i in nodes){
    title <- c("'", gsub(" ", "\n", Term(GOTERM[[as.character(i)]])), "'")
    title <- paste(title, collapse="")
    n <- paste(title, collapse='')
    out <- paste(c(out, n, collapse="; "))
  }
  out <- paste(c(out, "edge [color=grey] "), collapse='')
  for(i in 1:dim(edges)[1]){
    fr <- c("'", gsub(" ", "\n", Term(GOTERM[[edges[i, ]["from"]]])), "'")
    fr <- paste(fr, collapse='')
    to <- c("'", gsub(" ", "\n", Term(GOTERM[[edges[i, ]["to"]]])), "'")
    to <- paste(to, collapse='')
    edg <- c(fr, to)
    edg <- paste(edg, collapse=" -> ")
    out <- paste(c(out, edg), collapse=" ")
  }
  out <- paste(out, "}", collapse="")
  return(out)
}

getSubgraphNodes <- function(enrichment_results){
  nodes<-data.frame(id=character(), stringsAsFactors = FALSE)
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

buildNodeDF <- function(enrichment_results){
  nodes <- data.frame(nodes=character(), type=character(),
                      label=character(), style=character(),
                      color=character(), shape=character(),
                      data=character())
  nodeList <- getSubgraphNodes(enrichment_results)
  names(nodeList) <- c("nodes")
  nodes <- rbind(nodes, nodeList)
  # TODO the above code is losing the empty columns, need to add default values for them
  return(nodes)
  # add node types, colour, etc here
}



buildEdgeDF <- function(nodeDF){
  edges <- data.frame(from=character(), to=character(), rel=character())
  for(i in nodeDF[["nodes"]]){
    children <- as.list(GOBPCHILDREN[[i]])
    children <- children[!is.na(children)]
    if(length(children) > 0){
      in_common <- unique(intersect(children, nodeDF[["nodes"]]))
      for(j in in_common){
        edges <- rbind(edges, data.frame(from=j, to=i, rel="leading_to"))
      }
    }
  }
  return(edges)
}

nodes <- buildNodeDF(go_ids)
edges <- buildEdgeDF(nodes)
print(nodes)
print(edges)
gr <- create_graph(nodes_df = nodes, edges_df = edges, directed=TRUE)
render_graph(gr)

#DiagrammeR::grViz(gr)

# TODO
# change colour according to p-value or score
# make nodes equidistant
# add relationships edges (part_of, etc.)
# cluster by cellular component
# change buildgraph to use diagrammer
