#' Normalization and assessment
#'  
#' @param object Gretool object.
#' @param auto Whether to automatically select negative control, positive 
#'   evaluation and negative evaluation genes, default: TRUE. 
#' @param return.norm Whether to return normalized counts in object. By default,
#'   not return normalized counts for reducing memory costs. 
#' @param n.neg.control Number of negative control genes for RUV normalization, default: 1000. 
#' @param n.pos.eval Number of positive evaluation genes for wanted variation assessment, default: 500.
#' @param n.neg.eval Number of negative evaluation genes for unwanted variation assessment, default: 500.
#' @param neg.control Vector of negative control genes' id for RUV normalization, default: NULL. 
#' @param pos.eval Vector of positive evaluation genes' id for wanted variation assessment, default: NULL.
#' @param neg.eval Vector of negative evaluation genes' id for unwanted variation assessment, default: NULL.
#' @param scaling.method Vector of scaling methods that are applied to the data.
#'   Available methods are: \code{c("TC", "UQ", "TMM", "DESeq", "PossionSeq")}. 
#'   Select one or multiple methods. By default all scaling methods will be applied.
#' @param ruv.norm Whether to perform RUV normalization. 
#' @param ruv.k The number of factors of unwanted variation to be estimated from the data, default: 1.
#' @param ruv.drop The number of singular values to drop in the estimation of 
#'   unwanted variation, default: 0.  
#' @param eval.pam.k Integer or vector of integers indicates the number of 
#'   clusters for PAM clustering in performance evaluation, default: 2:6. 
#' @param eval.pc.n Integer indicates the evaluation metrics will be calculated 
#'   in the first nth PCs, default: 3.
#'
#' @return Gretool object.
#' @export
#' 
#' @importFrom SummarizedExperiment rowData assay
#' @importFrom stats setNames
Gretool <- function(object,
                  auto = TRUE, return.norm = FALSE,
                  n.neg.control = 1000, n.pos.eval = 500, n.neg.eval = 500,
                  neg.control = NULL, pos.eval = NULL, neg.eval = NULL,
                  scaling.method = c("TC", "UQ", "TMM", "DESeq", "PossionSeq"),
                  ruv.norm = TRUE, ruv.k = 1, ruv.drop = 0,
                  eval.pam.k = 2:6, eval.pc.n = 3) {
  
  # check ruv.k, eval.pam.k and eval.pc.n are least than the number of samples
  if (!all(c(ruv.k, eval.pam.k, eval.pc.n) < ncol(object))) {
    stop("Number of `ruv.k`, `eval.pam.k` and `eval.pc.n` should not exceed the number of samples.")
  }
  
  # retrieve parameters from object
  bio.group <- object$condition
  enrich.group <- object$enrich
  
  if (!any(is.na(object$batch))) {
    batch.group <- object$batch
  } else {
    batch.group <- NULL
  }
  
  spike.in.prefix <- object@parameter$spike.in.prefix
  input.id <- object@parameter$input.id
  enrich.id <- object@parameter$enrich.id
  synthetic.id <- object@parameter$synthetic.id
  
  # create group matrix
  sc_mat <-  CreateGroupMatrix(bio.group)
  enrich_mat <- CreateGroupMatrix(enrich.group)
  
  # get counts
  data <- SummarizedExperiment::assay(object)
  # sample counts
  counts_nsp <- data[grep(paste(c(spike.in.prefix, synthetic.id), collapse = "|"), rownames(data), invert = TRUE),]
  # spike-in counts
  counts_sp <- data[.safe_grep(spike.in.prefix, rownames(data)),]
  
  ## gene selection 
  if (auto) {
    genes.ls <- GeneSelection(object, 
                              n.neg.control = n.neg.control,
                              n.pos.eval = n.pos.eval,
                              n.neg.eval = n.neg.eval)
    
    neg.control.set <- genes.ls[["NegControl"]]
    pos.eval.set <- genes.ls[["PosEvaluation"]]
    neg.eval.set <- genes.ls[["NegEvaluation"]]
    
  } else {
    
    if (!all(neg.control %in% rownames(data))) {
      stop("`neg.control` are not presented in the rownames of count matrix.")
    } else {
      neg.control.set <- neg.control  
    }
    
    if (!is.null(pos.eval) & !all(pos.eval %in% rownames(data))) {
      stop("`pos.eval` are not presented in the rownames of count matrix.")
    } else {
      pos.eval.set <- pos.eval  
    }
    
    if (!is.null(neg.eval) & !all(neg.eval %in% rownames(data))) {
      stop("`neg.eval` are not presented in the rownames of count matrix.")
    } else {
      neg.eval.set <- neg.eval  
    }
    
  }
  
  # save gene set to object
  SummarizedExperiment::rowData(object)$NegControl <- SummarizedExperiment::rowData(object)$GeneID %in% neg.control.set
  SummarizedExperiment::rowData(object)$NegEvaluation <- SummarizedExperiment::rowData(object)$GeneID %in% neg.eval.set
  SummarizedExperiment::rowData(object)$PosEvaluation <- SummarizedExperiment::rowData(object)$GeneID %in% pos.eval.set
  
  ## apply normalization 
  cat("Apply normalization...\n")
  norm.ls <- ApplyNormalization(data,
                                scaling.method = scaling.method, 
                                ruv.norm = ruv.norm, ruv.k = ruv.k, ruv.drop = ruv.drop,
                                spike.in.prefix = spike.in.prefix,
                                synthetic.id = synthetic.id,
                                # below parameters are generated inside Gretool function
                                control.idx = neg.control.set, 
                                sc.idx = sc_mat, 
                                enrich.idx = enrich_mat)
  
  ## assessment 
  bio_group_index <- as.numeric(factor(bio.group, levels=unique(bio.group)))
  enrich_group_index <- as.numeric(factor(enrich.group, levels=unique(enrich.group)))
  if (!is.null(batch.group)) {
    batch_group_index <- as.numeric(factor(batch.group, levels=unique(batch.group)))
  } else {
    batch_group_index <- NULL
  }
  cat("Perform assessment...\n")
  norm.eval <- AssessNormalization(norm.ls,
                                   eval.pam.k = eval.pam.k,
                                   eval.pc.n = eval.pc.n,
                                   batch.group = batch_group_index,
                                   # below parameters are created inside function
                                   bio.group = bio_group_index, 
                                   enrich.group = enrich_group_index, 
                                   pos.eval.set = pos.eval.set,
                                   neg.eval.set = neg.eval.set)
  
  ## save metrics to object
  object@norm_metrics <- norm.eval$metrics
  ## save score to object
  object@norm_score <- norm.eval$score
  ## add run parameter to object
  parameter.run <- list(
    n.neg.control = n.neg.control,
    n.pos.eval = n.pos.eval,
    n.neg.eval = n.neg.eval,
    scaling.method = scaling.method,
    ruv.norm = ruv.norm,
    ruv.k = ruv.k,
    ruv.drop = ruv.drop,
    eval.pam.k = eval.pam.k,
    eval.pc.n = eval.pc.n,
    auto = auto,
    return.norm = return.norm
  )
  object@parameter <- c(object@parameter, parameter.run)
  
  # store normalization method names in object
  norm.methods <- names(norm.ls)
  object@counts$sample <- stats::setNames(vector("list", length(norm.methods)), nm = norm.methods)
  object@norm_factors$sample <- stats::setNames(vector("list", length(norm.methods)), nm = norm.methods)
  
  # store 'Raw' count matrix
  Counts(object, slot = "sample", method = "Raw") <- counts_nsp
  Counts(object, slot = "spike_in", method = "Raw") <- counts_sp
  
  if (return.norm == TRUE) {
    # store normalized counts
    object@counts$sample <- lapply(norm.ls, function(i) { i$dataNorm })
    # store normalization factors
    object@norm_factors$sample <- lapply(norm.ls, function(i) { i[c("normFactor", "adjustFactor", "alpha")] })
    # rename factor slot
    for (i in 1:length(object@norm_factors$sample)) { 
      names(object@norm_factors$sample[[i]]) <- c("normFactor", "adjustFactor", "alpha")
      }
  } 
  
  validObject(object)
  return(object)
}

#' Applies normalization on sequencing data
#'
#' @param data A un-normalized count data matrix of shape n x p, where n is the number of samples and p is the number of features. 
#' @param scaling.method Vector of normalization methods that are applied to the data.
#'   Available methods are: \code{c("TC", "UQ", "TMM", "DESeq", "PossionSeq")}. 
#'   Select one or multiple methods. By default all normalization methods will be applied.
#' @param ruv.norm Whether to perform RUV normalization. 
#' @param ruv.k The number of factors of unwanted variation to be estimated from the data.
#' @param ruv.drop The number of singular values to drop in the estimation of 
#'   unwanted variation, default drop the first singular value that represent the 
#'   difference between enrichment and input. 
#' @param control.idx Vector of the negative control genes for RUV normalization. 
#' @param sc.idx A numeric matrix specifying the replicate samples for which to 
#'   compute the count differences used to estimate the factors of unwanted variation.
#' @param enrich.idx Matrix with two rows indicating the column index of 
#'   enrichment and input samples in the raw/normalized count data matrix. 
#'   The first row is the column index of input and the second row is the 
#'   column index of enrichment samples.
#' @param spike.in.prefix A character specify the prefix of spike-in id. 
#' @param synthetic.id Character or vector of string specifying the name of synthetic RNAs.  
#' 
#' @return List of objects containing normalized data and associated normalization factors. 
#' @export
#' 
#' @importFrom pbapply pblapply
ApplyNormalization <- function(data, 
                               scaling.method = c("TC", "UQ", "TMM", "DESeq", "PossionSeq"),
                               ruv.norm = TRUE, 
                               ruv.k = 1, 
                               ruv.drop = 0, 
                               control.idx = NULL,
                               sc.idx = NULL,
                               enrich.idx = NULL,
                               spike.in.prefix = NULL,
                               synthetic.id = NULL) {
  
  scaling.method <- match.arg(scaling.method,
                              choices = c("TC", "UQ", "TMM", "DESeq", "PossionSeq"),
                              several.ok = TRUE)
  
  # scaling
  if (is.null(spike.in.prefix)) {
    cat("- Scaling... \n")
    data.scale <- pbapply::pblapply(1:length(scaling.method), function(i) {
      norm.f <- get(paste0("norm", scaling.method[i])) 
      data.norm <- norm.f(data)
    })
    names(data.scale) <- scaling.method
    data.norm <- data.scale
    data.raw <- list(dataNorm = data, normFactor = rep(1, ncol(data)))
    data.norm[["Raw"]] <- data.raw
  } else {
    data.spike.in <- data[grep(spike.in.prefix, rownames(data)),]
    data.non.spike.in <- data[grep(paste(c(spike.in.prefix, synthetic.id), collapse = "|"), rownames(data), invert = TRUE),]
    
    cat("- Scaling... \n")
    data.scale <- pbapply::pblapply(1:length(scaling.method), function(i) {
      norm.f <- get(paste0("norm", scaling.method[i])) 
      data.norm.spike.in <- norm.f(data.spike.in)
      data.norm.non.spike.in <- norm.f(data.non.spike.in)
      
      # return non-spike-in and negative control counts
      dataNorm <- rbind(data.norm.non.spike.in$dataNorm, data.norm.spike.in$dataNorm[control.idx,])
      
      return(list(
        dataNorm = dataNorm,
        normFactor = data.norm.non.spike.in$normFactor
      ))
    })
    names(data.scale) <- scaling.method
    data.norm <- data.scale
    data.raw <- list(dataNorm = data[c(grep(spike.in.prefix, rownames(data), invert=TRUE,value=TRUE),control.idx),],
                     normFactor = rep(1, ncol(data)))
    data.norm[["Raw"]] <- data.raw
  }
  
  # RUV normalization
  if (ruv.norm & !is.null(control.idx)) {
    # stop when provided ruv.k larger than sample size
    if (ruv.k > ncol(data)) {
      stop("Number of `ruv.k` must not exceed the number of samples.")
    }
    
    cat("- Regression-based normalization... \n")
    # generate all integrated methods
    # only perform RUVs, when enrichment is the only covariate of interest
    if (identical(sc.idx, enrich.idx)) {
      integrated.methods <- paste(rep(paste(rep(c("Raw",scaling.method),each=2), c("RUVg","RUVs"), sep="_"), each=ruv.k), paste0("k", 1:ruv.k), sep="_")
    } else {
      integrated.methods <- paste(rep(paste(rep(c("Raw",scaling.method),each=3), c("RUVg","RUVs","RUVse"), sep="_"), each=ruv.k), paste0("k", 1:ruv.k), sep="_")
    }
    
    ruv.ls <- pbapply::pblapply(1:length(integrated.methods), function(i) {
      # method.curr[1]: scaling; method.curr[2]: RUV; method.curr[3]: number of k
      method.curr <- unlist(strsplit(integrated.methods[i], split = "_"))
      
      # get current scaled data and scaling factors
      data.curr <- data.norm[[method.curr[1]]]$dataNorm
      normFactor.curr <- data.norm[[method.curr[1]]]$normFactor
      
      # apply all RUV
      if (method.curr[2] %in% c("RUVg", "RUVs", "RUVse")) {
        # switch sc.idx
        sc.idx <- switch(method.curr[2],
                         "RUVg" = NULL,
                         "RUVs" = sc.idx,
                         "RUVse" = enrich.idx)
        ruv.curr <- normRUV(data.curr,
                            control.idx = control.idx,
                            sc.idx = sc.idx,
                            method = method.curr[2],
                            k = as.numeric(gsub("k", "", method.curr[3])))
        ruv.curr$normFactor <- normFactor.curr
      }
      return(ruv.curr)
    })
    names(ruv.ls) <- integrated.methods
    
    data.norm <- c(data.norm, ruv.ls)
    
    # if use spike-in 
    # return only non-spike-in counts
    if (!is.null(spike.in.prefix)) {
      data.norm <- lapply(data.norm, function(x) {
        x$dataNorm <- x$dataNorm[!rownames(x$dataNorm) %in% control.idx,]
        return(x)
      })
    } 
  }
  return(data.norm)
}

#' Assessments of normalization performance
#'
#' @param data.ls List containing normalized counts and adjust factors for 
#'   adjusting unwanted variation. Output of \code{ApplyNormalization}. 
#' @param bio.group Vector of index indicating the column index of samples of 
#'   each biological groups in the raw/normalized count data matrix. 
#' @param enrich.group Vector of index indicating the column index of 
#'   enrichment and input samples in the raw/normalized count data matrix. 
#' @param batch.group Vector of index indicating the column index of 
#'   each batch groups in the raw/normalized count data matrix. 
#' @param eval.pam.k Integer or vector of integers indicates the number of 
#'   clusters for PAM clustering in performance evaluation, default: 2:6. 
#' @param eval.pc.n Integer indicates the evaluation metrics will be calculated 
#'   in the first nth PCs, default: 3.
#' @param log Whether to perform log2-transformation with 1 offset on data matrix, 
#'   default: TRUE. 
#' @param pos.eval.set Vector of genes id.
#' @param neg.eval.set Vector of genes id.
#'
#' @return List containing the metrics matrix and the ranking matrix, 
#'   both sorted by the score of methods from top to bottom. 
#' @export
#'
#' @importFrom stats dist cor lm na.omit var
#' @importFrom cluster silhouette 
#' @importFrom MatrixGenerics rowMedians colMedians colIQRs
#' @importFrom fpc pamk
#' @importFrom pbapply pblapply
AssessNormalization <- function(data.ls, 
                                bio.group = NULL, 
                                enrich.group = NULL, 
                                batch.group = NULL,
                                eval.pam.k = 2:6, 
                                eval.pc.n = 3, 
                                log = TRUE, 
                                pos.eval.set = NULL, 
                                neg.eval.set = NULL) {
  
  metrics.ls <- pbapply::pblapply(1:length(data.ls), function(i) {
    data <- as.matrix(data.ls[[i]]$dataNorm)
    # Clustering properties
    if (log) {
      data.log <- log2(data + 1)
    } else {
      data.log <- data
    }
    # remove constant genes
    data.log <- data.log[apply(data.log, 1, var, na.rm=TRUE) !=0, ]
    # PCA on expression matrix
    pca.expr <- prcomp(scale(t(data.log)))
    # compute right singular value by svd
    expr_sv <- svd(scale(t(data.log), center = TRUE, scale = TRUE),
                   nu = eval.pc.n, nv = 0)$u
    # calculate euclidean distance in the space of first k PCs (default: 3)
    dist.pca.expr <- dist(scale(pca.expr$x[, 1:eval.pc.n]), method = "euclidean")
    # dist.pca.expr <- dist(expr_sv, method = "euclidean")
    # silhouette width
    if (length(bio.group) == ncol(data)) {
      bio_sil <- mean(cluster::silhouette(bio.group, dist.pca.expr)[,"sil_width"])
    } else {
      bio_sil <- 0
    }
    if (length(enrich.group) == ncol(data)) {
      en_sil <- mean(cluster::silhouette(enrich.group, dist.pca.expr)[,"sil_width"])
    } else {
      en_sil <- 0
    }
    if (length(batch.group) == ncol(data)) {
      batch_sil <- mean(cluster::silhouette(batch.group, dist.pca.expr)[,"sil_width"])
    } else {
      batch_sil <- 0
    }
    
    prk <- fpc::pamk(pca.expr$x[,1:eval.pc.n], krange=eval.pam.k) # PAM clustering with user specified k
    pam_sil <- prk$pamobject$silinfo$avg.width
    
    # Global distribution properties
    data.log.rle <- data.log - rowMedians(data.log)
    # Mean squared Median RLE
    rle_med <- mean(colMedians(data.log.rle)^2)
    # Variance of IQR of RLE
    rle_iqr <- var(colIQRs(data.log.rle))
    
    # Association with control genes
    # wanted factors from positive set 
    if (!is.null(pos.eval.set)) {
      wv_factors <- svd(scale(t(data.log[rownames(data.log) %in% pos.eval.set,]), center = TRUE, scale = TRUE),
                        nu = eval.pc.n, nv = 0)$u
      # weighted coefficient of determination
      wv_cor <- 1 - sum(unlist(apply(expr_sv, 2, function(y) {
        lm(y ~ wv_factors)$residual
      })) ^ 2) / sum(scale(expr_sv, scale = FALSE) ^ 2)
    } else {
      wv_cor <- 0
    }
    
    # unwanted factors from negative set
    if (!is.null(neg.eval.set)) {
      uv_factors <- svd(scale(t(data.log[rownames(data.log) %in% neg.eval.set,]), center = TRUE, scale = TRUE),
                        nu = eval.pc.n, nv = 0)$u
      # weighted coefficient of determination
      uv_cor <- 1 - sum(unlist(apply(expr_sv, 2, function(y) {
        lm(y ~ uv_factors)$residual
      })) ^ 2) / sum(scale(expr_sv, scale = FALSE) ^ 2)
    } else {
      uv_cor <- 0
    }
    
    metrics <- c(
      BIO_SIM = bio_sil,
      EN_SIM = en_sil,
      BATCH_SIM = batch_sil,
      PAM_SIM = pam_sil,
      RLE_MED = rle_med,
      RLE_IQR = rle_iqr,
      WV_COR = wv_cor,
      UV_COR = uv_cor
    )
  })
  
  # reduce list of metrics into table
  # with methods in row and measures in column
  metrics <- data.frame(do.call(rbind, metrics.ls))
  rownames(metrics) <- names(data.ls)
  
  # multiplying by +/- 1 so that large values correspond to good performance
  score <- t(t(metrics) * c(1,1,-1,1,-1,-1,1,-1))  # BIO_SIM,EN_SIM,BATCH_SIM,PAM_SIM,RLE_MED,RLE_IQR,WV_COR,UV_COR
  # rank score
  ranked_score <- apply(na.omit(score), 2, rank, ties.method = "min")
  # mean score rank
  if (is.null(dim(ranked_score))) {
    mean_score_rank <- ranked_score
  } else {
    # if score all 1, remove it before scoring
    metrics.keep <- colSums(ranked_score==1) != nrow(ranked_score)
    mean_score_rank <- rowMeans(ranked_score[,metrics.keep])
  }
  
  ranked_score <- as.data.frame(ranked_score)
  ranked_score$SCORE <- mean_score_rank
  ranked_score <- ranked_score[order(mean_score_rank, decreasing = TRUE),]
  
  metrics <- metrics[order(mean_score_rank, decreasing = TRUE), ]
  
  return(list(
    metrics = metrics,
    score = ranked_score
  ))
}

#' Gene set selection 
#'
#' @description Select negative control genes for RUV normalization, positive 
#' and negative evaluation genes for assessment. 
#' 
#' @param object Gretool object
#' @param n.neg.control Number of negative control genes for RUV normalization, default: 1000. 
#' @param n.pos.eval Number of positive evaluation genes for wanted variation assessment, default: 500.
#' @param n.neg.eval Number of negative evaluation genes for unwanted variation assessment, default: 500.
#'
#' @return list of genes
#' @export
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom utils head
#' @importFrom stats model.matrix as.formula
GeneSelection <- function(object,
                          n.neg.control = 1000,
                          n.pos.eval = 500,
                          n.neg.eval = 500) {
  
  # retrieve parameters from object
  bio.group <- object$condition
  enrich.group <- object$enrich
  spike.in.prefix <- object@parameter$spike.in.prefix
  input.id <- object@parameter$input.id
  enrich.id <- object@parameter$enrich.id
  synthetic.id <- object@parameter$synthetic.id
  
  # get counts
  data <- SummarizedExperiment::assay(object)
  if (!is.null(spike.in.prefix)) {
    counts_nsp <- data[grep(spike.in.prefix, rownames(data), invert = TRUE),]
    counts_sp <- data[grep(spike.in.prefix, rownames(data)),]
  } else {
    counts_nsp <- data
  }
  
  cat("Gene set selection for normalization and assessment...\n")
  ### 1. negative control genes for RUV
  cat(paste("- The number of negative control genes for normalization:",n.neg.control,"\n"))
  #designMat <- stats::model.matrix(~0+enrich.group)
  if (!is.null(spike.in.prefix)) {
    deg.en <- edgeRDE(counts_sp,
                      group = enrich.group,
                      design.formula = stats::as.formula("~0+condition"),
                      contrast.df = data.frame(Group1=enrich.id, Group2=input.id)
    )
  } else {
    # find negative control genes from sample counts
    deg.en <- edgeRDE(counts_nsp,
                      group = enrich.group,
                      design.formula = stats::as.formula("~0+condition"),
                      contrast.df = data.frame(Group1=enrich.id, Group2=input.id)
    )
  }
  
  # top 1000 (default) non-sig de 
  res_tab <- deg.en$res.ls[[paste(enrich.id, input.id, sep="_")]]
  # res_tab <- subset(res_tab, FDR > 0.05)
  neg.control.set <- head(res_tab[order(res_tab$FDR, decreasing = TRUE),]$GeneID, n=n.neg.control)
  
  ### 2. positive evaluation genes (default 500)
  # if provided, preclude synthetic RNA from evaluation set 
  cat(paste("- The number of positive evaluation genes:",n.pos.eval,"\n"))
  deg.en <- edgeRDE(counts_nsp[!rownames(counts_nsp) %in% synthetic.id,],
                    group = enrich.group,
                    design.formula = stats::as.formula("~0+condition"),
                    contrast.df = data.frame(Group1=enrich.id, Group2=input.id)
  )
  res_tab <- deg.en$res.ls[[paste(enrich.id, input.id, sep="_")]]
  pos.eval.set <- head(res_tab[order(res_tab$FDR),]$GeneID, n=n.pos.eval)
  
  ### 3. negative evaluation genes (default 500)
  # if provided, preclude synthetic RNA from evaluation set 
  cat(paste("- The number of negative evaluation genes:",n.neg.eval,"\n"))
  de.all <- edgeRDE(counts_nsp[!rownames(counts_nsp) %in% synthetic.id,],
                    group = bio.group,
                    design.formula = stats::as.formula("~condition"),
                    coef = 2:length(unique(bio.group))
  )
  res_tab <- de.all$res.ls[[1]]
  neg.eval.set <- head(res_tab[order(res_tab$FDR, decreasing = TRUE),]$GeneID, n = n.neg.eval)
  
  return(list("NegControl" = neg.control.set, 
              "PosEvaluation" = pos.eval.set, 
              "NegEvaluation" = neg.eval.set))
}
