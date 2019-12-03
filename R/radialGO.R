# TODO Fix graph size
# TODO Add warning() and stop() for incorrect user input
# TODO Use message() to log
# TODO Use drop=F when subsetting DFs
# TODO Change color scale to better show the differences between small p-values
# TODO Fix edge colors not showing (seems like a bug in DiagrammeR or Graphviz)
# TODO Cluster nodes by cellular component localization


#' This function converts a vector of values("z") to a vector of color
#' levels. One must define the number of colors. The limits of the color
#' scale("zlim") or the break points for the color changes("breaks") can
#' also be defined. when breaks and zlim are defined, breaks overrides zlim.
#'
#' Function copied from https://menugget.blogspot.com/2011/09/converting-values-to-color-levels.html
val2col<-function(z, zlim, col = heat.colors(12), breaks){
   if(!missing(breaks)){
         if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
     }
   if(missing(breaks) & !missing(zlim)){
         breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
     }
   if(missing(breaks) & missing(zlim)){
         zlim <- range(z, na.rm=TRUE)
         zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
         zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
         breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
     }
   colorlevels <- col[((as.vector(z)-breaks[1])/(range(breaks)[2]-range(breaks)[1]))*(length(breaks)-1)+1] # assign colors to heights for each point
   colorlevels
}


#' Generate the subgraph implied by the list of nodes passed in.
#'
#' @param enrichment_results A named numerical vector containing GO IDs as names
#'                           and the corresponding score
#'
#' @return A data frame containing one column, named "id", with the nodes needed
#'         for building the graph
#'
#' @examples
#' \dontrun{
#' scores <- c(0.5, 0.2, 0.001)
#' names(scores) <- c("GO:0008150", "GO:1901360", "GO:0006139")
#' getSubgraphNodes(scores)
#' }
#'
#' @importFrom GO.db GOBPANCESTOR
#'
getSubgraphNodes <- function(enrichment_results){
  nodes <- data.frame(id=character(), stringsAsFactors = FALSE)
  for (i in enrichment_results[["GO.ID"]]){
    nodes <- rbind(nodes, data.frame(id=i, stringsAsFactors = FALSE))
    ancestors <- as.list(GO.db::GOBPANCESTOR[[i]]) # get ancestors of node
    ancestors <- ancestors[!is.na(ancestors)] #  remove NA fromlist
    ancestors <- ancestors[which(ancestors != "all")] # remove "all" ancestors
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
#'
#' @param goTerm A string containing the GO ID to use
#'
#' @param scores A named numerical vector containing GO IDs as names
#'               and the corresponding score
#'
#' @return A string containing the label text
#'
#' @examples
#' \dontrun{
#' scores <- c(0.5, 0.2, 0.001)
#' names(scores) <- c("GO:0008150", "GO:1901360", "GO:0006139")
#' buildLabel(names(scores[1]), scores)
#' }
#'
#' @importFrom AnnotationDbi Term
#'
buildLabel <- function(goTerm, scores){
  label <- paste(c(as.character(gsub(" ", "\n", AnnotationDbi::Term(GO.db::GOTERM[[goTerm]]))),
                                "\n\n", goTerm,
                                "\n\n", as.character(round(scores[scores$GO.ID == goTerm, "P.value"], 4))),
                 collapse='')
  return(label)
}

#' Builds a node data frame from the top GO terms using the enrichment results
#'
#' @param enrichment_results A named numerical vector containing GO IDs as names
#'                           and the corresponding score
#'
#' @param top The number of top nodes to use to build the graph
#'
#' @importFrom grDevices colorRampPalette
#'
buildNodeDF <- function(enrichment_results, top){
  n <- data.frame(nodes=character(), type=character(),
                  label=character(), style=character(),
                  color=character(), shape=character(),
                  data=numeric(), stringsAsFactors = FALSE)
  nodeList <- getSubgraphNodes(enrichment_results[with(enrichment_results, order(P.value)), ][1:top, ])
  # get node ids and p-values for nodeList
  checkList <- enrichment_results[enrichment_results$GO.ID %in% nodeList[ ,1], ]
  # order checklist
  checkList <- checkList[with(checkList, order(P.value)), ]
  # fix rownames to correspond to order
  rownames(checkList) <- seq(1:nrow(checkList))

  for (i in seq(along=1:dim(nodeList)[1])){
    nodeLabel <- buildLabel(nodeList[i, ], enrichment_results)
    if (as.character(nodeList[i, ]) %in% enrichment_results[[1]]){

      # TODO modify this --- get proper colours for p-value
      # use colours for nodes from nodeList rather than 1 to minimum
      # select colour by (P.value == x's p-value) rather than transforming it
      # colorIndex <- ceiling(enrichment_results[enrichment_results$GO.ID == nodeList[i, ], "P.value"] *
      #              nrow(enrichment_results))
      #checkList <- nodeList[with(nodeList, order("P.value")), ]
      #colorIndex <- as.numeric(rownames(checkList[checkList$GO.ID == nodeList[i, ], ]))
      #nodeColor <- grDevices::colorRampPalette(c("red", "yellow", "white"),
      #                                         space="rgb")(nrow(nodeList))[colorIndex]

      # TODO get p-values for nodeList, so we have both
      #nodeColor <- val2col(checkList[ ,2], c(1, (max(checkList[ ,2])*100)+1))[as.numeric(rownames(nodeList[i, ]))]
      #print(c(1, (max(checkList[ ,2])*100)+1))
      #print(val2col(checkList[ ,2], c(1, (max(checkList[ ,2])*100)+1)))
      #print(nodeColor)
      index <- as.numeric(rownames(checkList[checkList$GO.ID == nodeList[i, ], ]))
      print(index)
      nodeColor <- hcl(checkList[,2]/checkList[,2][1])[index]
      #print(hcl(checkList[,2]/checkList[,2][1]))
      #print(nodeColor)
      if(as.numeric(substr(nodeList[i, ], 4, nchar(nodeList[i, ]))) == 8150){
        nodeColor <- "gray"
      }

      newNode <- data.frame(id=as.numeric(substr(nodeList[i, ], 4, nchar(nodeList[i, ]))),
                            type="normal",
                            label= nodeLabel, style="filled",
                            fillcolor=nodeColor,
                            shape="circle", data=0,
                            fontsize = 80, fontcolor = "black",
                            stringsAsFactors = FALSE)
    } else {
      newNode <- data.frame(id=as.numeric(substr(as.character(nodeList[i, ]),
                                                 4, nchar(nodeList[i, ]))),
                            type="normal",
                            label= nodeLabel, style="filled",
                            fillcolor="white",
                            shape="circle", data=0,
                            fontsize = 80, fontcolor = "black",
                            stringsAsFactors = FALSE)
    }
    n <- rbind(n, newNode)
  }
  return(n)
}

#' Rebuild a GO ID string from its number
#'
#' @param goID A string or number containing the GO ID number
#'
#' @return The full GO ID
#'
#' @examples
#' rebuildGOID(8150)
#'
#' @export
rebuildGOID <- function(goID){
  return(paste(c("GO:", rep("0",(7-nchar(goID))), goID), collapse=''))
}

#' Determine what color the edge should take in the graph.
#'
#' @param edgeRelationships A string containing a relationship between GO terms
#'
#' @return A string containing the color the edge should take
#'
#' @examples
#' \dontrun{
#' getEdgeColor("is_a")
#' }
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
#'
#' @param nodeDF A data frame containing at least columns id, label,
#'
#' @return A data frame of edges (from, to, rel)
#'
#' @examples
#' \dontrun{
#' scores <- c(0.5, 0.2, 0.001)
#' names(scores) <- c("GO:0008150", "GO:1901360", "GO:0006139")
#' buildEdgeDF(buildNodeDF(scores, 2))
#' }
#' @importFrom AnnotationDbi Term
#' @import GO.db
buildEdgeDF <- function(nodeDF){
  edges <- data.frame(from=character(), to=character(), rel=character(),
                      style=character(), stringsAsFactors = FALSE)
  terms <- searchable::invert(AnnotationDbi::Term(GO.db::GOTERM))
  listOfNodes <- character()
  for (i in nodeDF[["id"]]){
    # build list of node GO IDs in the nodeDF here, to check with children
    go_id <- rebuildGOID(i)
    listOfNodes <- c(listOfNodes, go_id)
  }
  for (i in nodeDF[["id"]]){
    go_id <- rebuildGOID(i)
    children <- as.list(GO.db::GOBPCHILDREN[[go_id]])
    children <- children[!is.na(children)]
    if (length(children) > 0){
      in_common <- unique(intersect(children, listOfNodes))
      for (j in in_common){
        edgeColor <- ""
        lookUp <- searchable::invert(GO.db::GOBPCHILDREN[[go_id]])
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
#'
#' @param graph The graph to reduce
#'
#' @param cutoff The cutoff to be used (nodes above this cutoff will be removed)
#'
#' @return A graph with the nodes removed
#'
#' @importFrom DiagrammeR combine_edfs
reduceGraph <- function(graph, cutoff){
  edgeDF <- DiagrammeR::get_edge_df(graph)
  nodeDF <- DiagrammeR::get_node_df(graph)
  for(i in nodeDF[["id"]]){
    if(i == 8150){
      # GO:0008150 is the biological_process node (root)
      # we want to remove duplicate edges
      next
    } else {
      # get the p-value for this node from the label
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
              # check if edge is duplicate
              newEdge <- merge(newEdge, edgeAttrs)
              edgeDF <- DiagrammeR::combine_edfs(edgeDF, newEdge)
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
  # remove duplicated edges
  edgeDF <- edgeDF[!duplicated(edgeDF[,c("from", "to", "rel")]), ]
  return(DiagrammeR::create_graph(nodes_df = unique(nodeDF), edges_df = unique(edgeDF)))
}

#' Generate the graph visualization of GO enrichment scores
#'
#' @param scores A named numeric vector with GO IDs and corresponding scores
#'
#' @param top The number of nodes to use when generating the subgraph
#'
#' @param cutoff The cutoff to be used for hiding nodes
#'
#' @return A graph object to be used with render_graph()
#'
#' @examples
#' \dontrun{
#' scores <- c(0.5, 0.2, 0.001)
#' names(scores) <- c("GO:0008150", "GO:1901360", "GO:0006139")
#' generateGraph(scores, 2, 0.3)
#' }
#'
#' @export
#'
#' @import DiagrammeR
#'
generateGraph <- function(scores, top, cutoff){
  nodes <- buildNodeDF(scores, top)
  edges <- buildEdgeDF(nodes)


  gr <- DiagrammeR::create_graph(nodes_df = nodes, edges_df = edges, directed=TRUE)
  gr <- reduceGraph(gr, cutoff)
  gr <- DiagrammeR::add_global_graph_attrs(gr, attr="layout", value="twopi",
                               attr_type = "graph")
  gr <- DiagrammeR::add_global_graph_attrs(gr, attr="root", value="8150",
                               attr_type = "graph")
  gr <- DiagrammeR::add_global_graph_attrs(gr, attr="ranksep", value="20",
                               attr_type = "graph")
  gr <- DiagrammeR::add_global_graph_attrs(gr, attr="nodesep", value="20",
                               attr_type = "graph")
  gr <- DiagrammeR::add_global_graph_attrs(gr, attr="rank", value="same",
                               attr_type = "graph")
  gr <- DiagrammeR::add_global_graph_attrs(gr, attr="splines", value="curved",
                               attr_type = "graph")
  gr <- DiagrammeR::add_global_graph_attrs(gr, attr="forcelabels", value="true",
                               attr_type = "graph")
  gr <- DiagrammeR::add_global_graph_attrs(gr, attr="orientation", value="[lL]*",
                               attr_type = "graph")
  gr <- DiagrammeR::add_global_graph_attrs(gr, attr="overlap", value="true",
                               attr_type = "graph")
  gr <- DiagrammeR::add_global_graph_attrs(gr, attr="fixedsize", value="true",
                               attr_type = "node")
  gr <- DiagrammeR::add_global_graph_attrs(gr, attr="width", value="14",
                               attr_type = "node")
  gr <- DiagrammeR::add_global_graph_attrs(gr, attr="height", value="14",
                               attr_type = "node")
  gr <- DiagrammeR::add_global_graph_attrs(gr, attr="penwidth", value="1",
                               attr_type = "edge")
  gr <- DiagrammeR::add_global_graph_attrs(gr, attr="fontsize", value="15",
                                           attr_type = "node")

  return(gr)
}


# Shiny app --------------------------------------------------------------------
# Define UI ----
ui <- fluidPage(
  titlePanel("RadialGO - Gene Ontology enrichment visualization"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", h3("Import GO enrichment results"),
                placeholder = "No file selected", accept= c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")),
      numericInput("top",
                   h3("Choose number of top nodes to use when generating graph"),
                   value = 10),
      numericInput("cutoff",
                   h3("Choose cutoff p-value for nodes to be displayed"),
                   value = 0.5, step=0.01)
    ),
    mainPanel(
      h3("Output"),
      grVizOutput("image"),
    )
  )
)

#' Shiny app server function
#'
#' @param ui The UI object
#'
#' @param server The server object
#'
#' @import shiny
#'
server <- function(input, output) {
  output$image <- renderGrViz({
    inFile <- input$file

    if (is.null(inFile))
      return(NULL)

    address <- gsub("/", "\\\\", inFile$datapath)
    graph <- generateGraph(read.csv(address), input$top, input$cutoff)
    grViz(generate_dot(graph), width="100%", height="1000px")
  })
}


shinyApp(ui = ui, server = server)

# End of Shiny app -------------------------------------------------------------

# [END]
