library(tidyverse)
library(gplots)
library(enrichplot)
library(ggpubr)

filter_same_cols <- function(meta_df) {
  same <- list()
  for (col in colnames(meta_df)) {
    if (length(unique(meta_df[[col]])) == 1) {
      same[[col]] <- unique(meta_df[[col]])
      meta_df[[col]] <- NULL
    }
    if ("SRA" %in% colnames(meta_df)) {
      meta_df$SRA <- NULL
    }
    if ("BioSample" %in% colnames(meta_df)) {
      meta_df$BioSample <- NULL
    }
    if ("Reanalyzed.by" %in% colnames(meta_df)) {
      meta_df$Reanalyzed.by <- NULL
    }
  }
  return(list(meta_df, same))
}

meta_table_summary <- function(meta_df) {
  column <- colnames(meta_df)
  column_type <- as.vector(sapply(meta_df, typeof))
  df <- data.frame(column = column, column_type = column_type)
  num_sum <- data.frame(mean = sapply(meta_df, mean, na.rm = TRUE),
                        sds = sapply(meta_df, sd, na.rm = TRUE)) %>%
    mutate(summary = str_glue('{round(mean, 3)}(+/-{round(sds, 3)})'))
  df$numeric_summary <- num_sum$summary
  df$numeric_summary <- gsub('NA(+/-NA)', 'NA', df$numeric_summary, fixed = TRUE)
  df_nonnum <- data.frame(sum = sapply(meta_df, function(x) {length(unique(x))}))
  df$number_of_unique_value <- df_nonnum$sum
  return(df)
}

meta_common_summary <- function(lst) {
  for (name in names(lst)) {
    value <- lst[[name]]
    cat(sprintf('%s:\t%s\n', name, value))
  }
}

plot_sample <- function(meta_df, x_name, y_name, plot_type, group_by=NULL){

  if (plot_type == "scatter plot") {
    if (is.null(group_by)) {
      p <- ggplot(meta_df, aes(x = !!sym(x_name), y = !!sym(y_name))) +
        geom_point()
      return(p)
    }
    else {
      print("reach scatter plot with group")
      return(ggplot(meta_df, aes(x = !!sym(x_name), y = !!sym(y_name), color = !!sym(group_by))) +
               geom_point())
    }
  }
  else if (plot_type == "violin plot") {
    print("reach violin plot")
    if (is.null(group_by)) {
      return(ggplot(meta_df, aes(x = !!sym(x_name), y = !!sym(y_name))) +
               geom_violin()+geom_jitter(width=0.2, alpha=0.5))
    }
    else {
      return(ggplot(meta_df, aes(x = !!sym(x_name), y = !!sym(y_name), color = !!sym(group_by)))+
               geom_violin()+geom_jitter(width=0.2, alpha=0.5))
        }
  }
  else if (plot_type == "line plot") {
    print("reach line plot")
    if (is.null(group_by)) {
      print("reach line plot without group")
      return(ggplot(meta_df, aes(x = !!sym(x_name), y = !!sym(y_name))) +
        geom_line())
    }
    else {
      print("reach line plot with group")
      return(ggplot(meta_df, aes(x = !!sym(x_name), y = !!sym(y_name), color = !!sym(group_by))) +
        geom_line())
    }
  }
}

make_count_tibble <- function(count_df){
  n_gene <- nrow(count_df)
  n_sample <- ncol(count_df)
  n_zero <- apply(count_df, 1, function(x) sum(x != 0)*100/n_sample)
  df_var <- apply(count_df, 1, var)
  df_med <- apply(count_df, 1, median)
  processed_df <- count_df
  processed_df$var <- df_var
  processed_df$med <- df_med
  processed_df$non_zero <- n_zero
  processed_df$var_rank <- rank(-processed_df$var)*100/n_gene
  processed_df <- tibble::rownames_to_column(processed_df, "gene")%>%
    dplyr::as_tibble()
  return(processed_df)
}

anno_filter <- function (count_tibble, var_percentage, non_zero_percentage) {
  df <- count_tibble %>% mutate(keep = ifelse(var_rank <= var_percentage & non_zero >= non_zero_percentage, TRUE, FALSE))
  return(df)
}

plot_filtered_counts <- function (anno_tibble){
  p1 <- ggplot(anno_tibble, aes(x = var, y = med)) +
    geom_point(aes(color = keep)) +
    scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#FFC2C2")) +
    labs(x = "Variance(log10)", y = "Median count (log10)",color = "Pass filters") +
    scale_x_log10() +
    scale_y_log10()
  p2 <- ggplot(anno_tibble, aes(x = non_zero, y = med)) +
    geom_point(aes(color = keep)) +
    scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#FFC2C2")) +
    labs(x = "Percentage of non-zero counts", y = "Median count (log10)",color = "Pass filters") +
    scale_y_log10()
  p <- ggarrange(p1, p2, ncol = 2, nrow = 1)
  return(p)
}

plot_count_heatmap <- function (anno_tibble){
    anno_tibble <- dplyr::filter(anno_tibble, keep == TRUE) %>%
      dplyr::select(-var_rank, -var, -med, -non_zero,-keep) %>%
      tibble::column_to_rownames("gene") %>% as.matrix()

    anno_tibble <- log10(anno_tibble + 1)
    plt <- heatmap.2(anno_tibble, scale = "row",
                     col = bluered(50),trace = "none",
                     density.info = "none",keysize = 1)
    return(plt)
}

calc_pca <- function (anno_tibble){
  anno_tibble <- dplyr::filter(anno_tibble, keep == TRUE) %>%
      dplyr::select(-var_rank, -var, -med, -non_zero,-keep) %>%
      tibble::column_to_rownames("gene") %>% as.matrix()

  pca_results <- prcomp(scale(t(anno_tibble)), center=FALSE, scale=FALSE)
  return(summary(pca_results))
}

plot_pca <- function (pca_df, pc_x, pc_y, meta_df=NULL, join_col=NULL, group_by=NULL){
    pca_val <- pca_df$x %>% as.data.frame()
    pca_summary <- pca_df$importance %>% as.data.frame()
    pca_val <- pca_val[, c(pc_x, pc_y)]
    x_label <- paste0(pc_x, " (", round(pca_summary[2, pc_x]*100, 2), "%)")
    y_label <- paste0(pc_y, " (", round(pca_summary[2, pc_y]*100, 2), "%)")
    pca_val$type <- ifelse(substr(row.names(pca_val), 1, 1) == "H", "H", "C")
    if(!is.null(meta_df) & !is.null(join_col) & !is.null(group_by)){
        merge_df <- pca_val %>% rownames_to_column(join_col)%>%
            left_join(meta_df, by = join_col)
        p <- ggplot(merge_df, aes_string(x = !!sym(pc_x), y = !!sym(pc_y), color = !!sym(group_by))) +
          geom_point() +
          labs(x = x_label, y = y_label)
    }
    else{
       p <- ggplot(pca_val, aes_string(x = pc_x, y = pc_y, color = "type")) +
          geom_point() +
          labs(x = x_label, y = y_label)
    }
    return(p)
}

make_de_tibble <- function (de_df){
  de_tibble <- tibble::rownames_to_column(de_df, "gene_id") %>%
    dplyr::as_tibble()
  return(de_tibble)
}

format_number <- function(df) {
  df <- df %>% mutate_if(is.numeric, ~format(., scientific = TRUE, digits = 3))
  return(df)
}

volcano_plot <- function(dataf, x_name, y_name, p_col,slider, color1, color2) {
  ret <-  ggplot(dataf, aes(x = !!sym(x_name), y = -log10(!!sym(y_name)))) +
          geom_point(aes(color = -log10(!!sym(p_col))> -slider)) +
          theme_light()+
          labs(x = x_name,
               y = paste0("-Log10(",y_name , ")"),
               color = paste0(y_name , "<1x10^",slider))+
          scale_color_manual(values = c(color1, color2))+
          theme(legend.position = "bottom")
  return(ret)
}

make_gene_list <- function (de_df, gene_list,gene_rank){
  de_df <- de_df %>% dplyr::select(gene_list, gene_rank)
  names(de_df) <- c("g_id", "rank_value")
  de_df$g_id <- sub("\\.\\d+", "", de_df$g_id)
  return(de_df)
}

kegg_bar_plot <- function (gene_tibble){
  glist <- bitr(gene_tibble$g_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  kk <- enrichKEGG(gene = glist$ENTREZID,
                  organism = 'hsa',
                  pvalueCutoff = 0.05)
  p <-barplot(kk, showCategory = 20)
  return(p)
}

kegg_map_plot <- function (gene_tibble){
  glist <- bitr(gene_tibble$g_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  kk <- enrichKEGG(gene = glist$ENTREZID,
                  organism = 'hsa',
                  pvalueCutoff = 0.05)
  edo <- pairwise_termsim(kk)
  p <- emapplot(edo,cex_label_category=0.6)
  return(p)
}