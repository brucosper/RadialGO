#' ALL GO enrichment data
#'
#' Data from a GO enrichment analysis of dataset of Acute Lymphoblastic Leukemia
#' using topGO. This dataset was retrieved from the "ALL" library.
#'
#' @docType data
#'
#' @usage data(go_ids)
#'
#' @format A named numerical vector containing GO IDs as names and their score
#' in the enrichment analysis as values.
#'
#' @keywords datasets
#'
#' @references Chiaretti, S. et al. Gene expression profile of adult T-cell
#' acute lymphocytic leukemia identifies distinct subsets of patients with
#' different response to therapy and survival. Blood 103, 2771–2778 (2004).
#'
#' @examples
#' data(go_ids)
#' \dontrun{
#' render_graph(generateGraph(go_ids, 5, 0.2))
#' }
"go_ids"
data(go_ids, envir=environment())
