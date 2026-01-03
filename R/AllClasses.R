#' @rdname createGretool
#' @export
#'
setClass(
  Class = "Gretool",
  contains = "SummarizedExperiment",
  slots = list(
    counts = "list",
    norm_factors = "list",
    norm_metrics = "data.frame",
    norm_score = "data.frame",
    enrichment = "list",
    enrichment_filtered = "list",
    parameter = "list"
  )
)

#' Gretool object and constructor
#' 
#' @description \code{Gretool} object extends the \code{SummarizedExperiment} class. 
#' The \code{createGretool} is a easy constructor of \code{Gretool} object
#' 
#' @param data A un-normalized count data matrix of shape n x p, where n is the 
#'   number of samples and p is the number of features. 
#' @param col.data \code{data.frame} with at least two columns (indicate condition and enrich groups). 
#'   Rows of \code{col.data} correspond to columns of \code{data}. 
#' @param spike.in.prefix A character specify the prefix of spike-in id, e.g., "^FB" stands for fly spike-in id, default: NULL. 
#' @param input.id Input library id, must be consistent with the enrich column of \code{col.data}, e.g., "Input". 
#' @param enrich.id Enrich library id, must be consistent with the enrich column of \code{col.data}, e.g., "Enrich".  
#' @param synthetic.id Vector of synthetic RNA id, e.g. c("Syn1","Syn2"), default: NULL. 
#' 
#' @return Gretool object
#' 
#' @details Description of each slot: \cr
#'   \code{assay} \code{SummarizedExperiment::Assays} object, contains all counts. \cr
#'   \code{counts} list for holding raw/normalized counts of sample and spike_in. \cr
#'   \code{norm_factors} list for holding normalization factors of sample and spike_in. \cr
#'   \code{norm_metrics} data frame with normalization methods in row and metrics in columns. \cr
#'   \code{norm_score} data frame with normalization methods in row and scores in columns. \cr
#'   \code{enrichment} list for holding all differential analysis results. \cr
#'   \code{enrichment_filtered} lists for holding filtered differential analysis results. \cr
#'   \code{parameter} list of parameters.\cr
#'   
#' @export
#' 
#' @importFrom methods validObject
#'
createGretool <- function(data, 
                        col.data, 
                        spike.in.prefix = NULL, 
                        input.id = "Input", 
                        enrich.id = "Enrich",
                        synthetic.id = NULL
                        ) {
  # check if the order of rows in `col.data` matches to the order of columns in `data`
  if (!identical(rownames(col.data), colnames(data))) {
    stop("The rownames of `col.data` and the colnames of `data` must be identical.")
  }
  
  # parameters 
  params <- list(
    spike.in.prefix = spike.in.prefix,
    synthetic.id = synthetic.id,
    input.id = input.id,
    enrich.id = enrich.id
  )
  
  # create Assay object
  ## count_assay is a list for holding raw/normalized counts of sample and spike_in
  count_assay <- list(sample = list(), spike_in = list())
  ## enrichment_assay and enrichment_filtered_assay are lists for holding all/filtered differential analysis results
  enrichment_assay <- enrichment_filtered_assay <- list(sample=list(), spike_in=list())
  
  # factor slot
  norm_factors <- list(sample=list(), spike_in=list())
  
  # SummarizedExperiment object
  ## rowData for mapping gene id 
  rowDf <- S4Vectors::DataFrame(GeneID=rownames(data)) 
  
  # if spike-in prefix provided
  if (!is.null(spike.in.prefix)) {
    rowDf$SpikeIn <- rowDf$GeneID %in% grep(spike.in.prefix,rownames(data),value=TRUE)
  } else {
    rowDf$SpikeIn <- rep(FALSE, nrow(rowDf))
  }
  # if synthetic RNA id provided
  if (!is.null(synthetic.id)) {
    rowDf$Synthetic <- rowDf$GeneID %in% grep(paste(synthetic.id,collapse = "|"),rowDf$GeneID,value=TRUE)
  } else {
    rowDf$Synthetic <- rep(FALSE, nrow(rowDf))
  }
  ## colData for mapping samples
  colDf <- S4Vectors::DataFrame(
    id = colnames(data),
    col.data
  )
  # check condition column
  if (is.null(col.data$condition)) {
    stop("col.data must have a condition column to specify the grouping of samples.")
  } else {
    colDf$replicate <- countReplicate(colDf$condition)
  }
  
  # check enrich column
  if (is.null(col.data$enrich)) {
    stop("col.data must have a enrich column to specify the enrichment grouping of samples.")
  }
  
  if (is.null(col.data$batch)) colDf$batch <- NA_character_
  
  # create SummarizedExperiment object
  se <- SummarizedExperiment::SummarizedExperiment(S4Vectors::SimpleList(as.matrix(data)),
                                                   colData = colDf,
                                                   rowData = rowDf)
  # create object
  object <- methods::new("Gretool", 
                       se,
                       counts = count_assay,
                       norm_factors = norm_factors,
                       norm_metrics = data.frame(),
                       norm_score = data.frame(),
                       enrichment = enrichment_assay,
                       enrichment_filtered = enrichment_filtered_assay,
                       parameter = params
                       )
  
  # separate the raw counts
  object@counts$sample[["Raw"]] <- data[grep(paste(c(spike.in.prefix, synthetic.id), collapse = "|"), rownames(data), invert = TRUE),]
  object@counts$spike_in[["Raw"]] <- data[.safe_grep(spike.in.prefix, rownames(data)),]
  
  validObject(object)
  return(object)
}

setValidity("Gretool", function(object) {
  
  # check input.id and enrich.id match the enrich grouping
  if (!all(c(object@parameter$input.id, object@parameter$enrich.id) %in% object$enrich)) {
    return("The `input.id` and/or `enrich.id` must be the same as the enrich column of `col.data`.")
  }
  
  # check spike.in.prefix match the rownames
  if (!is.null(object@parameter$spike.in.prefix) & length(.safe_grep(object@parameter$spike.in.prefix, rownames(object))) < 1) {
    return("The `spike.in.prefix` does not match the rownames in the count matrix.")
  }

  # check synthetic.id match the rownames
  if (!is.null(object@parameter$synthetic.id) & !all(object@parameter$synthetic.id %in% rownames(object))) {
    return("The `synthetic.id` are not presented in the rownames of count matrix.")
  }
  
  # check length of condition match the sample size
  if (length(object$condition) != ncol(object)) {
    return("The number of elements in the condition column of `col.data` does not match the sample size.")
  }
  
  # check length of enrich match the sample size
  if (length(object$enrich) != ncol(object)) {
    return("The number of elements in the enrich column of `col.data` does not match the sample size.")
  }
  
  # check length of batch match the sample size
  if (!all(is.na(object$batch)) & length(object$batch) != ncol(object)) {
    return("The number of elements in the batch column of `col.data` does not match the sample size.")
  }
  
})
