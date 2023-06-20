ExpandAxes <- function(x, y){
  x.range <- range(x)
  y.range <- range(y)
  x.len <- abs(x.range[2] - x.range[1])
  y.len <- abs(y.range[2] - y.range[1])
  if (x.len>y.len){
    diff <- x.len-y.len
    y.range[1] <- y.range[1]-round(diff/2,digits = 0)
    y.range[2] <- y.range[2]+round(diff/2,digits = 0)
    return(list(x.range=x.range, y.range=y.range))
  } else if (y.len>x.len){
    diff <- y.len-x.len
    x.range[1] <- x.range[1]-round(diff/2,digits = 0)
    x.range[2] <- x.range[2]+round(diff/2,digits = 0)
    return(list(x.range=x.range, y.range=y.range))
  } else {
    return(list(x.range=x.range, y.range=y.range))
  }
  
}


PlotRNAFeature <- function(s.data,feature){
  require(tidyverse)
  require(Signac)
  require(Seurat)
  new.axes <-  ExpandAxes(x=s.data@reductions$umap@cell.embeddings[,1], y=s.data@reductions$umap@cell.embeddings[,2])
  DefaultAssay(s.data) <- "RNA"
  f.p <- FeaturePlot(s.data, features = convert_feature_identity(s.data, assay = "RNA", feature.format = "symbol", features = c(feature))) + ylim(new.axes$y.range) + xlim(new.axes$x.range) + labs(title="")
  return(f.p)
}


convert_feature_identity <- function(object, assay, features, feature.format = "symbol") {
  require(tidyverse)
  require(Signac)
  require(Seurat)
  #'
  #' Converts ENS ID -> gene symbol and vice versa
  #' Returns a vector of length(features) of either matches or NAs, in corresponding indices.
  #'
  #' Assumes libraries dplyr and seuratObject.
  #' Moreover, requires seuratObject[["assay"]] to contain df/tbl of ENS-symbol correspondences.
  #'
  # Protective tests
  if (!any(feature.format %in% c("ens", "symbol"))) {
    stop("Feature format should be a sting: either 'symbol' or 'ens'")
  }
  if (length(features) == 0) {
    stop("No features found in argument 'features'")
  }
  if (feature.format == "ens" & !all(grepl("*ENS", features))) {
    message("Warning: Found non-ENS ID for argument feature format 'ens'")
  }
  if (feature.format == "symbol" & any(grepl("*ENS", features))) {
    message("Warning: Found ENS ID for argument feature format 'symbol'")
  }
  # Diverging execution: case if provided features are ENSEBL IDs => conversion to symbols
  if (feature.format == "ens") {
    
    object.features <- object[[assay]][[]] %>%
      rownames_to_column(var = "gene_id") %>%
      as_tibble() %>%
      dplyr::select("gene_id", "feature_symbol") %>%
      dplyr::filter(gene_id %in% features)
    match.index <- match(features, object.features$gene_id, nomatch = NA)
    v.out <- sapply(match.index, function (i) { ifelse(is.na(i), NA, object.features$feature_symbol[i])})
    
    succ.matches <- sum(!is.na(v.out))
    n.tot <- length(features)
    
    msg <- str_interp("Instance: Found matching for ${succ.matches} features out of total ${n.tot} provided features")
    
    message(msg)
    
    return (v.out)
  }
  
  # Case: otherwise provided symbols => conversion to ENS IDs
  object.features <- object[[assay]][[]] %>%
    rownames_to_column(var = "gene_id") %>%
    as_tibble() %>%
    dplyr::select("gene_id", "feature_symbol") %>%
    dplyr::filter(feature_symbol %in% features)
  
  match.index <- match(features, object.features$feature_symbol, nomatch = NA)
  v.out <- sapply(match.index, function (i) { ifelse(is.na(i), NA, object.features$gene_id[i])})
  
  sprintf("Instance: Found matching for %d features out of total %d provided features", sum(!is.na(v.out)), length(features)) %>%
    print()
  
  return (v.out)
}


ATAC.QC.UMAPs <- function(s.data){
  require(patchwork)
  require(tidyverse)
  require(Signac)
  require(Seurat)
  Idents(s.data)<-"replicate"
  
  new.axes <-  ExpandAxes(x=s.data@reductions$umap@cell.embeddings[,1], y=s.data@reductions$umap@cell.embeddings[,2])
  p.1 <- DimPlot(s.data, reduction = "umap", label = FALSE, shuffle = TRUE) + ylim(new.axes$y.range) + xlim(new.axes$x.range) + guides(color=guide_legend(title="Replicate")) + labs(title="Biological replicates")
  p.2 <- FeaturePlot(s.data, features = "pct_reads_in_peaks") + ylim(new.axes$y.range) + xlim(new.axes$x.range)
  p.3 <- FeaturePlot(s.data, features = "peak_region_fragments") + ylim(new.axes$y.range) + xlim(new.axes$x.range)
  return(p.1 + p.2/p.3)
}


`%notin%` = Negate(`%in%`)
