# Analysis: extra ROC experiments
library(conflicted)
conflicted::conflicts_prefer(
  dplyr::select,
  dplyr::filter,
  dplyr::desc,
  dplyr::rename,
  dplyr::mutate,
  dplyr::where,
  .quiet = TRUE)
library(tidyverse)
library(pROC)
library(xgboost)
library(FactoMineR)

ps <- function(x, l=2){print(str(x, max.level=l))}

`%nin%` <- negate(`%in%`)

#' Assume that non-numeric columns are id_vars
filter_profile <- function(mat, motif_ids, id_cols=c('obs_id', 'subject_id', 'group')){
  mat[, c(id_cols, motif_ids)]
}

#' I want to compare motif thresholds, pcn over AUC.
#' Put dot in front b/c this is a macro
.comp_mtd <- function(df, dr_fx, sample_ids, xgb_fx){
  dr <- dr_fx(df)
  data_labels <- split_patients(sample_ids) # What's a better name for this?
  xgb <- run_xgb(feat, xgb_fx)
  xgb$roc
}


#' To allow prefiltering of motif matches before passing to binary feature
#' creation (raw data is over 5 million matches!)
motifs_at_thresh <- function(diffs, thresh, thresh_col){
  x <- diffs %>% filter(.data[[thresh_col]] >= thresh)
  x$motif_id
}

#' Here is a basic concept I should have clued into a long time ago. We
#' identify the motifs using the training set only, therefore to produce a
#' motif feature for any downstream analysis we should filter the dataset
#' (whether train, test, pooled) by the list of those ids 
filt_thresh <- function(df, diffs, thresh){
  filter(df, motif_id %in% motifs_at_thresh(thresh, diffs))
}

### Heirarchical Clustering Suite

#' Get hclustering of binary features, uses jaccard distance (ade4 method 1)
#' @returns d.hc: distance.hclustered
run_hclust <- function(mat, obs_id='obs_id', hc_method='ward.D2'){
  mat %>%
    column_to_rownames(obs_id) %>%
    select(where(is.numeric)) %>%
    ade4::dist.binary(method=1) %>%
    hclust(method=hc_method)
}

#' @returns df with cols label, color
#' Specifically for day0, day14 labels
get_label_colors <- function(labels, colors){
  labels %>%
    data.frame(label = .) %>%
    mutate(color = case_when(
      endsWith(label, '0') ~ colors[[1]],
      endsWith(label, '14') ~ colors[[2]]
    ))
}

#' Color dendrogram nodes/edges based on dataset factors
#' @param var_col: df with cols label, color
#' @returns fx to be used in dendrapply
construct_dendapply <- function(var_col, lwd=2){
  return ( function(n) {
      if(is.leaf(n)) {
        a <- attributes(n)
        col <- var_col[var_col$label == a$label, 'color']
        attr(n, 'label') <- a$label
        attr(n, "nodePar") <- list(
          col = col,
          cex = 0, # hide symbols at tip of leaf
          font = 4,
          lab.col = col)
        attr(n, 'edgePar') <- list( lwd = lwd, col = col)
      }
      return(n)
  })
}

#' Even tho this is a macro, I might want to keep it in src since it
#' demonstrates the very nifty dendrapply(construct_dendapply) usage
#' @param mot: df, m*n where feat vars are numeric
#' @examples
#' To save this plot use dev.copy/off.
#' ppi <- 120
#' dev.copy(png, height=11*ppi, width=11*ppi, res=ppi, 'plot.png') 
#' dev.off()
runphylo <- function(mat, hclust_fx=run_hclust){
  lwd <- 2
  d.hc <- hclust_fx(mat)
  var_col <- get_label_colors(d.hc$labels, colors)
  dend  <- d.hc %>%
    as.dendrogram() %>%
    dendrapply(construct_dendapply(var_col, lwd=lwd))
  par(lwd = lwd, family = 'monospaced', ps = 12, font = 2 )
  # plots the dendrogram as circular, obj returned is still 'regular' dendrogram
  # can plot it like a fan with plot(dend), etc
  dendextend::circlize_dendrogram(dend)
}

#' dev copy/off is the only way I've found to save the circlize_dendrogram plots
#' Settings used for the poster (big figs): 11 inches x 120 ppi
#' Have since learned that cowplot::save_fig is much more pleasant to use, so
#' should use that to set aspect ratio if need to tweak dims for manuscript
.saveplot <- function(filename){
  outpath <- "export/dendro"
  dir.create(outpath, showWarnings = FALSE)
  ppi <- 120
  dev.copy(png, height=11*ppi, width=11*ppi, res=ppi, file.path(outpath, filename))
  dev.off()
}

### End Hierarchical Clustering Suite

#' Get the roc of a single dim red'd feature for bundling into a ggroc plot
#' Supply pcn within the dr_fx argument
#' Globals: xgb_function 
get_dr_feature <- function(df, dr_fx, sample.ids){
  feat <- dr_fx(df) %>%
    split_patients(sample.ids)
  xgb <- run_xgb(feat, xgb_function)
  roc <- xgb$roc
  return(roc)
}


# I do not love this fx
.build_roc_feats <- function(n_dims, feat_lens, mm.file, npx_all){
  features <- list()
  npx_features <- list()

  for (len in feat_lens) {
    # These stupid names are going to become important down the line, so keep motif-len
    # for non red and npx-red-len-dim for red
    mot <- mm.file %>% bin_mots_for_tree(len)
    features[[paste0('motif-', len)]] <-
      mot %>%
      split_patients(
        sample.ids,
        paste0('len = ', len)
      )

    npx_features[[paste0('npx', len)]] <-
      npx_all %>%
      select(any_of(c('label', 'subject_id', npx_sig))) %>%
      split_patients(
        sample.ids,
        paste0('len = ', len)
      )

    for (dim in n_dims) {
      res.pca <- mot %>%
        select(where(is.numeric)) %>%
        PCA(ncp = dim, scale.unit=FALSE, graph=FALSE)

      U <- res.pca$svd$U %>%
        data.frame() %>%
        cbind(mot %>% select(!where(is.double)), .)

      features[[paste('motif-red', len, dim, sep='-')]] <-
        mot %>%
        reduce_dims(dim) %>%
        split_patients(
          sample.ids,
          paste0('len = ', len, ', pcn = ', n_dims)
        )

      npx_features[[paste('npx-red', len, dim, sep='-')]] <-
        npx_all %>%
        select(any_of(c('label', 'subject_id', npx_sig))) %>%
        reduce_dims(dim) %>%
        split_patients(
          sample.ids,
          paste0('len = ', len, ' PCN = ', dim)
        )
    }
  }
  list(features, npx_features)
}

# Moved here to be able to edit and update w/ so()
.run_rocs <- function(features, group, palette, pal_offset=3){
  rocs <- lapply(features, function(x) run_xgb(x, xgb_function)) %>%
    map(~.[['roc']])
  aucs <- rocs$rocs %>% map(., ~ . [['auc']][1])

  list(
    rocs =  rocs,
    group = group,
    palette = get_palette(features, palette, pal_offset)
  )
}

#' Used for the gene ontology dotplots
groupGO_result <-  function(gene_list, go.level, key = 'ENSEMBLPROT', filt=TRUE){
  library(org.Hs.eg.db)
  ggo = clusterProfiler::groupGO(
    gene = gene_list,
    OrgDb = 'org.Hs.eg.db',
    keyType = key,
    ont = 'BP',
    level = go.level,
    readable = T)
  df <-  ggo@result %>%
    arrange(desc(Count)) %>%
    relocate(Description, Count) %>%
    # GeneRatio is a chr, which excel will maul into a date
    separate_wider_delim(
      cols = GeneRatio,
      delim = '/',
      names = c('x','y')) %>%
    # GeneRatio is numeric, safe for excel
    mutate(GeneRatio = as.numeric(x) / as.numeric(y)) %>%
    select(-c(x,y)) %>%
    relocate(GeneRatio, .after=Count)
  if (filt) {df = df %>% filter(Count > 0) }
  return(df)
}

# @returns res.pca
pca_fx <- function(mat){
  mat %>%
    select(where(is.numeric)) %>%
    FactoMineR::PCA(graph=FALSE)
}

#' @returns int, first PCN where cum var > thresh
find_pcn <- function(eig, var_thresh = 0.95){
  eig <- eig %>%
    data.frame() %>%
    mutate(pcn = row_number()) %>%
    filter_at(
      vars(starts_with('cumulative')),
      all_vars(. > 95)
    )
  eig$pcn[1]
}

# Quick barplot of from res.pca$eig
plot_eig <- function(eig){
  barplot(eig[,1], main="Eigenvalues", names.arg=1:nrow(eig))
}

#' Reduce cols of dtype double to vectors of n_dims len
#' Prepend back non-double cols
#' @returns data.frame
reduce_dims <- function(df, n_dims, scale=FALSE, plot_eigen=FALSE) {
  if (n_dims == 0) {
    return(df)
  }
  pc <- df %>%
    select(where(is.numeric)) %>%
    FactoMineR::PCA(ncp = n_dims, scale.unit=scale, graph=FALSE)

  if (plot_eigen) {
    barplot(pc$eig[,1],main="Eigenvalues",names.arg=1:nrow(pc$eig))
  }

  # Combine pcs with the other variables
  pc$svd$U %>%
    data.frame() %>%
    cbind( df %>% select(!where(is.double)), .)
}

#' Get rid of the few color in the palette, as it is usually extremely faint
#' @param palette: chr, color_brewer palette name
get_palette <- function(feature_list, palette, offset=1){
  n <- length(feature_list) + offset
  scales::brewer_pal(palette = palette)(n)[1+offset:n]
}

#' Get rid of the few color in the palette, as it is usually extremely faint
#' @param palette: chr, color_brewer palette name
get_palette2 <- function(len, palette, offset=1){
  n <- len + offset
  scales::brewer_pal(palette = palette)(n)[1+offset:n]
}

#' add item name with value to each child list in a nested list
add_subitem <- function(nested_list, name, value){
  for (i in 1:length(nested_list)){
    nested_list[[i]][[name]] <- value
  }
  nested_list
}

#' Filter df by subject_id to return train/test split
#'
#' @param df data to filter
#' @param sample.ids int vec of sample ids
#' @param negative If TRUE, return df w/ subject_ids NOT in set
#' @returns data.frame
#' @export
filter_ids <- function(df, sample.ids, negative=FALSE) {
  if (negative) {
    df = dplyr::filter(df, !subject_id %in% sample.ids)
  } else {
    df = dplyr::filter(df, subject_id %in% sample.ids)
  }
  return(df)
}

#' get motif profile specifically for decision trees
#' @param n_motifs: int, returns all motifs if null
#' @returns data.frame
bin_mots_for_tree <- function(motif_profile, n_motifs=NULL){
  motif_profile %>% 
    mutate(
      label = case_when(
          group == 'day0' ~ 0,
          group == 'day14' ~ 1),
      subject_id = as.character(subject_id)) %>%
    mutate(label = as.character(label)) %>%
    relocate(label, subject_id) %>%
    select(label, subject_id, matches('^[0-9]'))
}


#' Get table of most different motifs
#'
#' Wrap that allows returning a subset the diffs table
#' For train/test split, calculate a single diffs table using the train set.
#' Then use that table to calculate the motif profile feature for both sets
#' @param df matched/detected motifs
#' @param n_most_diff return this number of motifs, with ties. If NULL, return
#'   all motifs
#' @returns data.frame diffs_table
#' @export
diffs_table <- function(df, n_most_diff=NULL) {
  diffs <- get_diffs_df(df)
  if (is.numeric(n_most_diff)) {
    diffs <- diffs %>%
      slice_max(order_by = abs_diff, n = n_most_diff, with_ties = TRUE)}
  return(diffs)
}

#' Calculate the difference in motif presence between groups.
#' Specifically for zhong2021 dataset
#' DAY 0 - DAY 14 (negative diffs are enriched in disease condition)
#'
#' @param df matched/detected motifs
#' @returns df containing motif_id, counts per group, and diff of counts,
#'   arranged by abs_diff desc
get_diffs_df <- function(df, group='group', motif_id='motif_id') {
  # TODO: Ensure This is only for two groups. Dynamically find the values for
  # the groups and make sure that the output is clear which was subtracted from
  # which (ex. instead of phenotype can have group_a_dominant, etc)
  diffs = df %>%
    group_by(.data[[group]]) %>%
    count(.data[[motif_id]]) %>%
    pivot_wider(names_from = group, values_from = 'n')

  # Missing value means zero count
  diffs[is.na(diffs)] = 0
  diffs %>%
    mutate(
      diff = day0 - day14,
      abs_diff = abs(diff),
      phenotype = case_when(
        diff > 0 ~ 'disease',
        diff < 0 ~ 'recovery'),
      freq = day0 + day14) %>%
    arrange(desc(abs_diff))
}

#' Make motif profile into a vector of binary vals
#'
#' @param m motif_diffs table
#' @param columns column names to pull from the motif file
#' @returns a wide dataframe where each row is a sample's motif profile
#' @export
motif_profile_bin <- function(m, columns = c('obs_id', 'subject_id', 'group', 'motif_id')) {
  # long df only contains instances where motifs are present. after pivoting
  # wide, the motifs absent in a sample will be NA. assign those to 0
  m$value = 1
  m <- m %>%
    select(all_of(columns))
    pivot_wider(id_cols=c('obs_id', 'subject_id', 'group'), names_from='motif_id', values_from='value') %>%
    mutate_all(~replace_na(.,0))
}

#' code zhong2021 for classifier
group_to_label <- function(df) {
  df %>%
    mutate(
    label = case_when(
      group == 'day0' ~ '0',
      group == 'day14' ~ '1',)) %>%
    relocate(label) %>%
    select(!group)
}

#' get differentially expressed protein profile (zhong2021 s6) specifically for decision trees
#' @returns data.frame
npx_for_tree <- function(table_s6_file, npx_assembled_file){
  read.csv(npx_assembled_file) %>%
    filter(protein %in% read.csv(table_s6_file)$Protein) %>%
    select(subject_id, group, protein, npx) %>%
    # deal with the 4 or so compound IDs
    distinct(protein, subject_id, group, .keep_all=TRUE) %>%
    pivot_wider(
      id_cols = c(subject_id, group),
      names_from = protein, values_from = npx) %>%
    group_to_label() %>%
    mutate(subject_id = as.character(subject_id))
}

#' use these names to filter table s4 vector
get_sig_proteins <- function(table_s6) {
  if (is.character(table_s6)){
    table_s6 <- read.csv(table_s6)
  }
  table_s6 %>%
    rename_with(tolower)  %>%
    distinct(protein) %>%
    unlist() %>%
    as.vector() # otherwise will have elements named 'protein<n>'
}

npx_assembled_to_matrix <- function(npx_assembled) {
  npx_assembled %>%
    rename_with(tolower) %>%
    select(symbol, subject_id, group, npx) %>%
    # deal with the 4 or so compound IDs
    distinct(symbol, subject_id, group, .keep_all=TRUE) %>%
    pivot_wider(id_cols = c(subject_id, group),
                names_from = symbol,
                values_from = npx) %>%
    group_to_label() %>%
    mutate(subject_id = as.character(subject_id))
}

#' split patients specifically for xgboost classifier/trees
#' @returns list( name, train, test)
split_patients = function(dat, sample.ids, name=NULL) {
  # could do this as a train_mask column in the future
	datalabel <- function(D){
		out = list()
		out$data = D %>%
			select(where(is.numeric)) %>%
			as.matrix()
		out$label = D %>%
			mutate(label = as.numeric(label)) %>%
			select(label)
		return(out)
	}
  D = list()
  D$name = name
  D$train = dat %>%
		filter_ids(sample.ids) %>%
		datalabel()
  D$test = dat %>%
		filter_ids(sample.ids, negative=T) %>%
		datalabel()
  return(D)
}


# model and rocs ----
#' feed in ex M$train
xgb_function <- function(data, label) {
  # in python can use binary:logistic and pred with predict_proba
  model = xgboost(data = data, label = label,
                max.depth = 2, eta = 0.1, nthread = 2, nrounds = 30,
                objective = "binary:logistic",
                verbose = 0)
  return(model)
}

getpred <- function(model, test) {
  pred = predict(model, test$data)
  err = mean(as.numeric(pred > 0.5) != test$label)
  return(pred)
}

#' run xgboost, get roc, importance
#'
#' @param D list(train/test$data/label)
run_xgb <- function(D, xgb_fx) {
	# ghetto out
	res = list()
	model = xgb_fx(D$train$data, D$train$label[[1]])
	nb = xgb.importance(model = model)
	pred = getpred(model,D$test)
	r = roc(D$test$label[[1]], pred,
	        direction = '<', levels = c(0, 1))
	res$importance = nb
	res$model = model
	res$roc = r
	res$pred = pred
	return(res)
}

#' @param X: list( name, train, test )
get_roc <- function(X) {
	X$xgb <- run_xgb(X, xgb_function)
  roc <- list()
	rocname <- paste0( X$name, '\nAUC = ', X$xgb$roc$auc, '\n')
  roc[[rocname]] <- X$xgb$roc
  return(roc)
}

get_roclist <- function(feature_list){
  roclist = list()
  for (i in 1:length(feature_list)){
    X <- feature_list[[i]]
    roc <- get_roc(X)
    roclist <- append(roclist, roc)
  }
  # sort by auc for plot legend/palette to look nice
  # auc is the 9th element in each sublist of roclist
  # https://stackoverflow.com/questions/24203361/r-sorting-list-by-values-in-nested-list
  roclist <- roclist[
    order(sapply(roclist, `[[`, i=9))
  ]
  roclist
}

#' @params mm.file, sample.ids pulled from global
#' @params name list(name = name)
#' @params n_dims, reduce to ndims w/ pca, disabled if 0
#' @returns list(name, train, test)
get_motif_split <- function(n, name=NULL, n_dims=0) {
  df <- mm.file %>%
    bin_mots_for_tree(n)
  if (n_dims > 0) {
    df <- reduce_dims(df, n_dims)
  }
  df %>%
    split_patients(sample.ids, name)
}

