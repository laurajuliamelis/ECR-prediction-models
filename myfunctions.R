nimpute <- function(DAT, method = "any") {
  if (method == "avg") {
    p <- mean(sapply(DAT, function(x) mean(is.na(x))))
  } else if (method == "avg_gt_0") {
    p <- sapply(DAT, function(x) mean(is.na(x)))
    p <- mean(p[p > 0])
  } else if (method == "any") {
    p <- mean(!complete.cases(DAT))
  }
  return(max(5, round(100*p)))
}

CalculateAucFromDxy <- function(validate) {
  stopifnot(class(validate) == "validate") # Test if the object is correct
  
  ## Calculate AUCs from Dxy's
  aucs <- (validate["Dxy", c("index.orig","training","test","optimism","index.corrected")])/2 + 0.5
  n <- validate["Dxy", c("n")]                ## Get n
  res <- rbind(validate, AUC = c(aucs, n))    ## Combine as result
  res["AUC","optimism"] <- res["AUC","optimism"] - 0.5  ## Fix optimism
  res  ## Return results
}

naive.auc <- function(DF, i) {
  ## DF       # Original full sample
  ## DF[i,]   # Bootstrap sample (use index i as row numbers)
  
  ## Develop model on bootstrap sample
  boot_model = train(formula, data=DF[i,], method="naive_bayes", tuneLength=2, 
                     trControl = fitControl, metric="ROC")
  
  ## Measure AUC on bootstrap sample using bootstrap model (training set AUC)
  boot_predicted     <- predict(boot_model, type = "prob")[,2]
  boot_roc           <- roc(DF[i,"Resposta"] ~ boot_predicted)
  boot_auc           <- as.numeric(boot_roc$auc)
  
  ## Measure AUC on original sample using bootstrap model (testing set AUC)
  original_predicted <- predict(boot_model, type = "prob", newdata = DF)[,2]
  original_roc       <- roc(DF[ ,"Resposta"] ~ original_predicted) # No index i
  original_auc       <- as.numeric(original_roc$auc)
  
  ## training set AUC - testing set AUC = optimism in AUC
  auc_optimism       <- (boot_auc - original_auc)

  c(training = boot_auc, testing = original_auc, optimism = auc_optimism)
}


gbm.auc <- function(DF, i) {
  ## DF       # Original full sample
  ## DF[i,]   # Bootstrap sample (use index i as row numbers)
  
  ## Develop model on bootstrap sample
  boot_model = train(formula, data=DF[i,], method="gbm", tuneLength=2, 
                     trControl = fitControl, metric="ROC", verbose=F)
  
  ## Measure AUC on bootstrap sample using bootstrap model (training set AUC)
  boot_predicted     <- predict(boot_model, type = "prob")[,2]
  boot_roc           <- roc(DF[i,"Resposta"] ~ boot_predicted)
  boot_auc           <- as.numeric(boot_roc$auc)
  
  ## Measure AUC on original sample using bootstrap model (testing set AUC)
  original_predicted <- predict(boot_model, type = "prob", newdata = DF)[,2]
  original_roc       <- roc(DF[ ,"Resposta"] ~ original_predicted) # No index i
  original_auc       <- as.numeric(original_roc$auc)
  
  ## training set AUC - testing set AUC = optimism in AUC
  auc_optimism       <- (boot_auc - original_auc)
  
  c(training = boot_auc, testing = original_auc, optimism = auc_optimism)
}


rf.auc <- function(DF, i) {
  ## DF       # Original full sample
  ## DF[i,]   # Bootstrap sample (use index i as row numbers)
  
  ## Develop model on bootstrap sample
  boot_model = train(formula, data=DF[i,], method="rf", tuneLength=2, 
                     trControl = fitControl, metric="ROC")
  
  ## Measure AUC on bootstrap sample using bootstrap model (training set AUC)
  boot_predicted     <- predict(boot_model, type = "prob")[,2]
  boot_roc           <- roc(DF[i,"Resposta"] ~ boot_predicted)
  boot_auc           <- as.numeric(boot_roc$auc)
  
  ## Measure AUC on original sample using bootstrap model (testing set AUC)
  original_predicted <- predict(boot_model, type = "prob", newdata = DF)[,2]
  original_roc       <- roc(DF[ ,"Resposta"] ~ original_predicted) # No index i
  original_auc       <- as.numeric(original_roc$auc)
  
  ## training set AUC - testing set AUC = optimism in AUC
  auc_optimism       <- (boot_auc - original_auc)
  
  c(training = boot_auc, testing = original_auc, optimism = auc_optimism)
}


svm.auc <- function(DF, i) {
  ## DF       # Original full sample
  ## DF[i,]   # Bootstrap sample (use index i as row numbers)
  
  ## Develop model on bootstrap sample
  boot_model = train(formula, data=DF[i,], method="svmLinear", tuneLength=2, 
                     trControl = fitControl, metric="ROC")
  
  ## Measure AUC on bootstrap sample using bootstrap model (training set AUC)
  boot_predicted     <- predict(boot_model, type = "prob")[,2]
  boot_roc           <- roc(DF[i,"Resposta"] ~ boot_predicted)
  boot_auc           <- as.numeric(boot_roc$auc)
  
  ## Measure AUC on original sample using bootstrap model (testing set AUC)
  original_predicted <- predict(boot_model, type = "prob", newdata = DF)[,2]
  original_roc       <- roc(DF[ ,"Resposta"] ~ original_predicted) # No index i
  original_auc       <- as.numeric(original_roc$auc)
  
  ## training set AUC - testing set AUC = optimism in AUC
  auc_optimism       <- (boot_auc - original_auc)
  
  c(training = boot_auc, testing = original_auc, optimism = auc_optimism)
}


get_model_stats = function(x, precision=60) {
  
  # remember old number formatting function
  # (which would round and transforms p-values to formats like "<0.01")
  old_format_np = rms::formatNP
  # substitute it with a function which will print out as many digits as we want
  assignInNamespace("formatNP", function(x, ...) formatC(x, format="f", digits=precision), "rms")
  
  # remember old width setting
  old_width = options('width')$width
  # substitute it with a setting making sure the table will not wrap
  options(width=old_width + 4 * precision)
  
  # actually print the data and capture it
  cap = capture.output(print(x))
  
  # restore original settings
  options(width=old_width)
  assignInNamespace("formatNP", old_format_np, "rms")
  
  #model stats
  stats = c()
  stats$R2.adj = str_match(cap, "R2 adj\\s+ (\\d\\.\\d+)") %>% na.omit() %>% .[, 2] %>% as.numeric()
  
  #coef stats lines
  coef_lines = cap[which(str_detect(cap, "Coef\\s+S\\.E\\.")):(length(cap) - 1)]
  
  #parse
  coef_lines_table = suppressWarnings(readr::read_table(coef_lines %>% stringr::str_c(collapse = "\n")))
  colnames(coef_lines_table)[1] = "Predictor"
  
  list(
    stats = stats,
    coefs = coef_lines_table
  )
}

ody_plot_violindotbox2 <- function (data, x, y, no_violin = FALSE, compare = FALSE, p_adj = "fdr", 
                                    brackets_pos = 1.08, ...) {
  rlang::check_installed(c("ggpubr", "gghalves"))
  if (!is.factor(data[[x]])) {
    data[[x]] <- factor(data[[x]])
  }
  if (no_violin) {
    p <- ggplot2::ggplot(data, ggplot2::aes(.data[[x]], .data[[y]])) + 
      ggplot2::geom_boxplot(outliers = FALSE) + 
      ggplot2::geom_jitter(shape = 16, position = position_jitter(0.05), color = "grey30", size=0.8)
  }
  else {
    p <- ggplot2::ggplot(data, ggplot2::aes(.data[[x]], .data[[y]], color=.data[[x]])) + 
      gghalves::geom_half_violin(side = "r") +
      scale_color_manual(values=c("#A6CEE3", "#1F78B4"))+ 
      ggplot2::geom_boxplot(width = 0.15, outliers = FALSE, fill=c("#A6CEE3", "#1F78B4"), color="black") + 
      gghalves::geom_half_point(side = "l", alpha = 0.8,position=position_nudge(0.03),transformation = position_jitter(0.07), color = "grey40", size=0.8) +
      scale_fill_manual(c("#A6CEE3", "#1F78B4"))
  }
  if (compare) {
    stats <- ggpubr::compare_means(formula(glue::glue("{y} ~ {x}")), 
                                   data = data, p.adjust.method = p_adj)
    y_pos <- max(dplyr::pull(data, .data[[y]]), na.rm = TRUE) * 
      brackets_pos
    p <- p + ggpubr::geom_bracket(data = stats, ggplot2::aes(xmin = .data[["group1"]], 
                                                             xmax = .data[["group2"]], label = gtsummary::style_pvalue(.data[["p.adj"]])), 
                                  y.position = y_pos, ...)
  }
  data <- dplyr::select(data, .data[[x]], .data[[y]])
  attr(p, "data") <- data
  p
}