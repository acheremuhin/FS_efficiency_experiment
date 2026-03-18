library(readr)
dependence_metrics <- function(x, y, which = 1:11, nbins_mi = 10) {
  # --- checks ---
  if (!is.numeric(x) || !is.numeric(y)) stop("x and y must be numeric vectors.")
  if (length(x) != length(y)) stop("x and y must have the same length.")
  if (any(!is.finite(which)) || any(which < 1) || any(which > 11)) {
    stop("which must contain integers from 1 to 11.")
  }
  which <- unique(as.integer(which))
  
  # --- clean ---
  x2 <- as.numeric(unlist(x))
  y2 <- as.numeric(unlist(y))
  ok <- is.finite(x2) & is.finite(y2)
  x2 <- x2[ok]; y2 <- y2[ok]
  if (length(x2) < 3) stop("Too few finite observations after cleaning.")
  if (sd(x2) == 0 || sd(y2) == 0) stop("One of the variables is constant (sd=0).")
  
  # --- package requirements per metric ---
  pkg_by_metric <- list(
    `4`  = "Hmisc",
    `5`  = "TauStar",
    `6`  = "energy",
    `7`  = "infotheo",
    `8`  = "minerva",
    `9`  = "XICOR",
    `10` = "acepack",
    `11` = "dHSIC"
  )
  
  need_pkgs <- unique(unlist(pkg_by_metric[as.character(which)]))
  need_pkgs <- need_pkgs[!is.na(need_pkgs)]
  
  for (p in need_pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p)
    }
  }
  
  # --- compute only requested metrics ---
  out_metric <- character(0)
  out_value  <- numeric(0)
  
  add <- function(name, val) {
    out_metric <<- c(out_metric, name)
    out_value  <<- c(out_value, as.numeric(val))
  }
  
  if (1 %in% which) add("Pearson correlation",  cor(x2, y2, method = "pearson"))
  if (2 %in% which) add("Spearman correlation", cor(x2, y2, method = "spearman"))
  if (3 %in% which) add("Kendall tau",          cor(x2, y2, method = "kendall"))
  
  if (4 %in% which) {
    xy <- cbind(x = x2, y = y2)
    hd <- Hmisc::hoeffd(xy)
    add("Hoeffding's D", hd$D[1, 2])
  }
  
  if (5 %in% which) {
    add("Bergsma–Dassios tau", TauStar::tStar(x2, y2))
  }
  
  if (6 %in% which) {
    add("Distance correlation", energy::dcor(x2, y2))
  }
  
  if (7 %in% which) {
    xb <- infotheo::discretize(x2, disc = "equalfreq", nbins = nbins_mi)
    yb <- infotheo::discretize(y2, disc = "equalfreq", nbins = nbins_mi)
    add("Mutual information", infotheo::mutinformation(xb, yb))
  }
  
  if (8 %in% which) {
    add("MIС", minerva::mine(x2, y2)$MIC)
  }
  
  if (9 %in% which) {
    xi_obj <- XICOR::xicor(x2, y2)
    xi_val <- if (is.list(xi_obj) && !is.null(xi_obj$xi)) as.numeric(xi_obj$xi) else as.numeric(xi_obj)
    add("Chatterjee xi", xi_val)
  }
  
  if (10 %in% which) {
    ace_fit <- acepack::ace(x2, y2)
    add("Renyi maximal correlation", cor(ace_fit$tx, ace_fit$ty, use = "complete.obs"))
  }
  
  if (11 %in% which) {
    hsic <- dHSIC::dhsic.test(as.matrix(x2), as.matrix(y2))
    add("HSIC statistic", hsic$statistic)
  }
  
  data.frame(metric = out_metric, value = out_value, stringsAsFactors = FALSE)
}
# ============================================================
# dependence_matrix_rank_2:
# returns list:
# [[1]] values matrix (numeric)   : rows = predictors, cols = metrics
# [[2]] ranks  matrix (numeric)   : rows = predictors, cols = metrics
# [[3]] y vector (original df[[y_name]])
# [[4]] X matrix (numeric)        : rows = observations, cols = predictors
# ============================================================

dependence_matrix_rank_2 <- function(df, y_name, which = 1:11, nbins_mi = 10) {
  if (!is.data.frame(df)) stop("df must be a data.frame.")
  if (!is.character(y_name) || length(y_name) != 1) stop("y_name must be a single column name.")
  if (!(y_name %in% names(df))) stop("y_name not found in df.")
  if (any(which < 1) || any(which > 11)) stop("which must contain integers from 1 to 11.")
  which <- unique(as.integer(which))
  
  y <- df[[y_name]]
  if (!is.numeric(y)) stop("Dependent variable must be numeric.")
  
  # numeric predictors (exclude y)
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  x_cols <- setdiff(num_cols, y_name)
  if (length(x_cols) == 0) stop("No numeric predictors found.")
  
  # metric names (consistent with dependence_metrics)
  test_x <- df[[x_cols[1]]]
  metric_cols <- dependence_metrics(test_x, y, which = which, nbins_mi = nbins_mi)$metric
  
  # -------- values ----------
  values_mat <- matrix(
    NA_real_,
    nrow = length(x_cols),
    ncol = length(metric_cols),
    dimnames = list(x_cols, metric_cols)
  )
  
  for (i in seq_along(x_cols)) {
    x <- df[[x_cols[i]]]
    tmp <- dependence_metrics(x, y, which = which, nbins_mi = nbins_mi)
    for (k in seq_len(nrow(tmp))) {
      values_mat[i, tmp$metric[k]] <- tmp$value[k]
    }
  }
  
  # -------- ranks ----------
  rank_strongest_first <- function(v, use_abs = FALSE) {
    vv <- if (use_abs) abs(v) else v
    r <- rep(NA_real_, length(vv))
    ok <- is.finite(vv)
    r[ok] <- rank(-vv[ok], ties.method = "average")
    r
  }
  
  is_abs_metric <- function(name) {
    name %in% c("Pearson correlation", "Spearman correlation", "Kendall tau")
  }
  
  ranks_mat <- matrix(
    NA_real_,
    nrow = nrow(values_mat),
    ncol = ncol(values_mat),
    dimnames = dimnames(values_mat)
  )
  
  for (mc in metric_cols) {
    ranks_mat[, mc] <- rank_strongest_first(values_mat[, mc], use_abs = is_abs_metric(mc))
  }
  
  # -------- X matrix (original predictors) ----------
  # Rows = observations, Cols = predictors in the same order as rownames(values_mat)
  X_mat <- as.matrix(df[, x_cols, drop = FALSE])
  # Ensure numeric storage (in case of odd numeric-like classes)
  X_mat <- apply(X_mat, 2, function(z) as.numeric(unlist(z)))
  X_mat <- as.matrix(X_mat)
  colnames(X_mat) <- x_cols
  
  # Return list with X as 4th element
  list(values_mat, ranks_mat, y, X_mat)
}

# ============================================================
# select_features_from_dm
# Input:
#   dm      : list from dependence_matrix_rank_2()
#             [[1]] values_mat (predictors x metrics)
#             [[2]] ranks_mat  (predictors x metrics)
#             [[3]] y vector
#             [[4]] X_mat      (observations x predictors)
#   variant : 1..5 (selection rule)
#   metric  : metric name (string) OR column index in dm[[1]]
#   param   : variant-specific parameter:
#             1) k (integer) OR fraction in (0,1)
#             2) threshold (numeric)
#             3) quantile in (0,1) OR top-percent in (1,100]
#             4) B (integer) number of shadow generations (optional; overrides B)
#             5) NULL (not used)
#   B       : default number of shadow generations for variant 4
#   q_shadow: quantile over max-shadow scores across generations (e.g. 0.95)
#   elbow_method: "maxdist" or "maxdrop"
# Output:
#   character vector of selected predictor names
# ============================================================

select_features_from_dm <- function(dm,
                                    variant,
                                    metric = 1,
                                    param = NULL,
                                    B = 50,
                                    q_shadow = 0.95,
                                    elbow_method = c("maxdist", "maxdrop")) {
  elbow_method <- match.arg(elbow_method)
  
  if (!is.list(dm) || length(dm) < 3) stop("dm must be a list from dependence_matrix_rank_2().")
  
  values_mat <- dm[[1]]
  ranks_mat  <- dm[[2]]
  y          <- dm[[3]]
  X_mat      <- if (length(dm) >= 4) dm[[4]] else NULL
  
  if (is.null(dim(values_mat)) || is.null(dim(ranks_mat)))
    stop("dm[[1]] and dm[[2]] must be matrices.")
  
  vars <- rownames(values_mat)
  if (is.null(vars)) vars <- as.character(seq_len(nrow(values_mat)))
  
  # --- metric selection ---
  if (is.character(metric)) {
    if (!(metric %in% colnames(values_mat))) stop("metric name not found in dm[[1]] columns.")
    mc <- metric
  } else {
    metric <- as.integer(metric)
    if (metric < 1 || metric > ncol(values_mat)) stop("metric index out of range.")
    mc <- colnames(values_mat)[metric]
  }
  
  v <- as.numeric(values_mat[, mc])
  use_abs <- mc %in% c("Pearson correlation", "Spearman correlation", "Kendall tau")
  v_eff <- if (use_abs) abs(v) else v
  
  # sorted (strongest first)
  ord <- order(v_eff, decreasing = TRUE, na.last = TRUE)
  vars_ord <- vars[ord]
  v_ord <- v_eff[ord]
  m <- length(v_ord)
  
  # ------------------------------------------------------------
  # 1) Cardinal selection: top-k or top-fraction
  # ------------------------------------------------------------
  if (variant == 1) {
    if (is.null(param)) stop("variant=1 requires param (k or fraction).")
    if (!is.numeric(param) || length(param) != 1) stop("param must be a single number.")
    if (param > 0 && param < 1) {
      k <- ceiling(param * m)
    } else {
      k <- as.integer(param)
    }
    k <- max(1, min(k, m))
    return(vars_ord[seq_len(k)])
  }
  
  # ------------------------------------------------------------
  # 2) Threshold selection: fixed threshold in metric scale
  # ------------------------------------------------------------
  if (variant == 2) {
    if (is.null(param)) stop("variant=2 requires param (threshold).")
    if (!is.numeric(param) || length(param) != 1) stop("param must be a single number.")
    t <- as.numeric(param)
    return(vars[which(v_eff >= t)])
  }
  
  # ------------------------------------------------------------
  # 3) Relative threshold:
  #    - param in (0,1): quantile q (keep >= Q_q)
  #    - param in (1,100]: top param% features
  # ------------------------------------------------------------
  if (variant == 3) {
    if (is.null(param)) stop("variant=3 requires param (quantile in (0,1) OR percent in (1,100]).")
    if (!is.numeric(param) || length(param) != 1) stop("param must be a single number.")
    p <- as.numeric(param)
    
    if (p > 0 && p < 1) {
      thr <- as.numeric(stats::quantile(v_eff, probs = p, na.rm = TRUE, type = 7))
      return(vars[which(v_eff >= thr)])
    } else if (p > 1 && p <= 100) {
      k <- ceiling((p / 100) * m)
      k <- max(1, min(k, m))
      return(vars_ord[seq_len(k)])
    } else {
      stop("For variant=3: param must be in (0,1) (quantile) or (1,100] (percent).")
    }
  }
  
  # ------------------------------------------------------------
  # 4) Adaptive shadow threshold (4b) — FAST VERSION:
  #    Permute y once per generation; compute ONLY selected metric mc for each X_j.
  #    Threshold = quantile(q_shadow) of max shadow scores across B generations.
  # ------------------------------------------------------------
  if (variant == 4) {
    if (!is.null(param)) {
      if (!is.numeric(param) || length(param) != 1) stop("param must be a single number for variant=4 (B).")
      B <- as.integer(param)
    }
    if (is.null(X_mat)) stop("variant=4 requires dm[[4]] = X_mat (observations x predictors).")
    if (is.null(colnames(X_mat))) stop("dm[[4]] must have colnames.")
    if (!all(vars %in% colnames(X_mat))) stop("X_mat colnames must include all predictors (rownames(dm[[1]])).")
    
    # align X columns order to vars
    X <- X_mat[, vars, drop = FALSE]
    
    # y as numeric vector
    yv0 <- as.numeric(unlist(y))
    
    # --- compute ONLY selected metric mc (same names as dependence_metrics) ---
    calc_metric_single <- function(xvec, yvec, mc, use_abs, nbins_mi) {
      x2 <- as.numeric(unlist(xvec))
      y2 <- as.numeric(unlist(yvec))
      ok <- is.finite(x2) & is.finite(y2)
      x2 <- x2[ok]; y2 <- y2[ok]
      if (length(x2) < 3) return(NA_real_)
      if (sd(x2) == 0 || sd(y2) == 0) return(NA_real_)
      
      val <- switch(
        mc,
        "Pearson correlation" = cor(x2, y2, method = "pearson"),
        "Spearman correlation" = cor(x2, y2, method = "spearman"),
        "Kendall tau" = cor(x2, y2, method = "kendall"),
        
        "Hoeffding's D" = {
          if (!requireNamespace("Hmisc", quietly = TRUE)) install.packages("Hmisc")
          hd <- Hmisc::hoeffd(cbind(x = x2, y = y2))
          hd$D[1, 2]
        },
        
        "Bergsma–Dassios tau" = {
          if (!requireNamespace("TauStar", quietly = TRUE)) install.packages("TauStar")
          TauStar::tStar(x2, y2)
        },
        
        "Distance correlation" = {
          if (!requireNamespace("energy", quietly = TRUE)) install.packages("energy")
          energy::dcor(x2, y2)
        },
        
        "Mutual information" = {
          if (!requireNamespace("infotheo", quietly = TRUE)) install.packages("infotheo")
          xb <- infotheo::discretize(x2, disc = "equalfreq", nbins = nbins_mi)
          yb <- infotheo::discretize(y2, disc = "equalfreq", nbins = nbins_mi)
          infotheo::mutinformation(xb, yb)
        },
        
        "MIС" = {
          if (!requireNamespace("minerva", quietly = TRUE)) install.packages("minerva")
          minerva::mine(x2, y2)$MIC
        },
        
        "Chatterjee xi" = {
          if (!requireNamespace("XICOR", quietly = TRUE)) install.packages("XICOR")
          xi_obj <- XICOR::xicor(x2, y2)
          if (is.list(xi_obj) && !is.null(xi_obj$xi)) xi_obj$xi else xi_obj
        },
        
        "Renyi maximal correlation" = {
          if (!requireNamespace("acepack", quietly = TRUE)) install.packages("acepack")
          ace_fit <- acepack::ace(x2, y2)
          cor(ace_fit$tx, ace_fit$ty, use = "complete.obs")
        },
        
        "HSIC statistic" = {
          if (!requireNamespace("dHSIC", quietly = TRUE)) install.packages("dHSIC")
          hs <- dHSIC::dhsic.test(as.matrix(x2), as.matrix(y2))
          hs$statistic
        },
        
        stop("Unknown metric name mc: ", mc)
      )
      
      val <- as.numeric(val)
      if (use_abs) val <- abs(val)
      val
    }
    
    max_shadow <- numeric(B)
    
    for (b in seq_len(B)) {
      # permute y ONCE
      yv <- sample(yv0, size = length(yv0), replace = FALSE)
      
      sh_scores <- rep(NA_real_, ncol(X))
      for (j in seq_len(ncol(X))) {
        sh_scores[j] <- calc_metric_single(X[, j], yv, mc, use_abs, nbins_mi)
      }
      max_shadow[b] <- max(sh_scores, na.rm = TRUE)
    }
    
    thr <- as.numeric(stats::quantile(max_shadow, probs = q_shadow, na.rm = TRUE, type = 7))
    return(vars[which(v_eff >= thr)])
  }
  # ------------------------------------------------------------
  # 5) Elbow method:
  #    - maxdrop: pick k at max adjacent drop s_k - s_{k+1}
  #    - maxdist: max distance to straight line from (1,s1) to (m,sm)
  # ------------------------------------------------------------
  if (variant == 5) {
    s <- v_ord
    ok <- is.finite(s)
    s <- s[ok]
    vo <- vars_ord[ok]
    n <- length(s)
    if (n < 3) return(vo)
    
    if (elbow_method == "maxdrop") {
      drops <- s[-n] - s[-1]
      k <- which.max(drops)
      return(vo[seq_len(k)])
    } else {
      xk <- seq_len(n)
      x1 <- 1; y1 <- s[1]
      x2 <- n; y2 <- s[n]
      num <- abs((y2 - y1) * xk - (x2 - x1) * s + x2 * y1 - y2 * x1)
      den <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
      d <- num / den
      k <- which.max(d)
      return(vo[seq_len(k)])
    }
  }
  
  stop("variant must be one of 1,2,3,4,5.")
}

# ============================================================
# Experiment runner: all metrics (1:11) x rules (1:5)
# Outputs LONG table: metric_id / rule_id / metric_name / indicator / value
# No "operations" metric.
# Truth: predictors with prefix "NS" are influencing, "NNS" are non-influencing.
# ============================================================

# ---- truth from prefixes NS / NNS ----
get_truth_sets <- function(vars) {
  true_pos <- vars[grepl("^NS",  vars)]
  true_neg <- vars[grepl("^NNS", vars)]
  list(pos = true_pos, neg = true_neg)
}

# ---- quality metrics for selected set ----
fs_quality <- function(selected, vars_all) {
  tr <- get_truth_sets(vars_all)
  P <- tr$pos; N <- tr$neg
  
  if (length(P) + length(N) == 0) {
    return(list(
      tanimoto_distance = NA_real_,
      pct_true_influencing_found = NA_real_,
      pct_true_noninfluencing_found = NA_real_,
      pct_missing_influencing = NA_real_,
      pct_included_noninfluencing = NA_real_
    ))
  }
  
  # Confusion wrt truth
  TP <- length(intersect(selected, P))
  FN <- length(setdiff(P, selected))
  TN <- length(setdiff(N, selected))
  
  # included non-influencing among known NNS
  FP_NNS <- length(intersect(selected, N))
  
  # Tanimoto/Jaccard distance to "true influencing set"
  A <- unique(selected)
  B <- unique(P)
  denom <- length(A) + length(B) - length(intersect(A, B))
  jac_sim <- if (denom == 0) NA_real_ else length(intersect(A, B)) / denom
  tanimoto_dist <- if (is.na(jac_sim)) NA_real_ else 1 - jac_sim
  
  pct_tp <- if (length(P) == 0) NA_real_ else 100 * TP / length(P)
  pct_fn <- if (length(P) == 0) NA_real_ else 100 * FN / length(P)
  
  pct_tn <- if (length(N) == 0) NA_real_ else 100 * TN / length(N)
  pct_fp <- if (length(N) == 0) NA_real_ else 100 * FP_NNS / length(N)
  
  list(
    tanimoto_distance = tanimoto_dist,
    pct_true_influencing_found = pct_tp,
    pct_true_noninfluencing_found = pct_tn,
    pct_missing_influencing = pct_fn,
    pct_included_noninfluencing = pct_fp
  )
}

# ---- run one (metric_id, rule_id) case ----
make_formula <- function(y_name, ss) {
  if (length(ss) == 0) {
    as.formula(paste0(y_name, " ~ 1"))
  } else {
    reformulate(termlabels = ss, response = y_name)
  }
}

run_one_case <- function(df, y_name, dm,
                         metric_id, rule_id,
                         param_list,
                         B_default = 50, q_shadow = 0.95,
                         seed = 1L) {
  
  values_mat <- dm[[1]]
  vars_all <- rownames(values_mat)
  metric_name <- colnames(values_mat)[metric_id]
  
  # choose param for this rule (can be NULL)
  param <- param_list[[as.character(rule_id)]]
  
  # seed for reproducibility (important for variant=4 permutations)
  set.seed(seed + 1000L * metric_id + 10L * rule_id)
  
  # memory snapshot before (proxy)
  gc_before <- gc()
  mem_before <- sum(gc_before[, "used"])
  
  # time + selection + model
  t <- system.time({
    ss <- tryCatch(
      select_features_from_dm(dm,
                              variant = rule_id,
                              metric = metric_name,
                              param = param,
                              B = B_default,
                              q_shadow = q_shadow),
      error = function(e) character(0)
    )
    
    fml <- fml <- make_formula(y_name, ss)
    fit <- tryCatch(lm(fml, data = df), error = function(e) NULL)
  })
  
  # memory snapshot after (proxy)
  gc_after <- gc()
  mem_after <- sum(gc_after[, "used"])
  mem_delta <- mem_after - mem_before
  
  # model size (proxy, bytes)
  fit_size <- if (is.null(fit)) NA_real_ else as.numeric(object.size(fit))
  
  # selection quality
  q <- fs_quality(ss, vars_all)
  
  # (optional but useful) how many selected
  indicators <- c(
    time_sec = as.numeric(t["elapsed"]),
    memory_delta_gc_units = as.numeric(mem_delta),      # gc "used" units (not bytes)
    model_object_size_bytes = as.numeric(fit_size),
    tanimoto_distance = as.numeric(q$tanimoto_distance),
    pct_true_influencing_found = as.numeric(q$pct_true_influencing_found),
    pct_true_noninfluencing_found = as.numeric(q$pct_true_noninfluencing_found),
    pct_missing_influencing = as.numeric(q$pct_missing_influencing),
    pct_included_noninfluencing = as.numeric(q$pct_included_noninfluencing),
    n_selected = as.numeric(length(ss))
  )
  
  data.frame(
    metric_id = metric_id,
    rule_id = rule_id,
    metric_name = metric_name,
    indicator = names(indicators),
    value = as.numeric(indicators),
    stringsAsFactors = FALSE
  )
}

# ---- full grid runner ----
run_full_grid <- function(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          nbins_mi = 10,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=50, `5`=NULL),
                          B_default = 50, q_shadow = 0.95,
                          seed = 1L) {
  
  # compute dm ONCE for all metrics requested
  dm <- dependence_matrix_rank_2(df, y_name = y_name, which = which_metrics, nbins_mi = nbins_mi)
  
  out <- list()
  idx <- 1L
  for (m in which_metrics) {
    for (v in rules) {
      out[[idx]] <- run_one_case(df, y_name, dm,
                                 metric_id = m,
                                 rule_id = v,
                                 param_list = param_list,
                                 B_default = B_default,
                                 q_shadow = q_shadow,
                                 seed = seed)
      idx <- idx + 1L
    }
  }
  do.call(rbind, out)
}

# ---------------- Example ----------------
#### Первая часть ####
I_Beta_1 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/I_Beta_1.csv")
I_Beta_2 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/I_Beta_2.csv")
I_Beta_3 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/I_Beta_3.csv")
I_Beta_4 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/I_Beta_4.csv")
I_Beta_5 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/I_Beta_5.csv")
I_Beta_6 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/I_Beta_6.csv")
I_Beta_7 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/I_Beta_7.csv")
I_Beta_8 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/I_Beta_8.csv")
I_Beta_9 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/I_Beta_9.csv")
I_Beta_10 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/I_Beta_10.csv")
df <- I_Beta_1
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                           which_metrics = 1:11,
                           rules = 1:5,
                           param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                           q_shadow = 0.95,
                           seed = 1)
file_name <- "I_Beta_1_results.csv"
write.csv(res_long, file_name, row.names = FALSE)

df <- I_Beta_2
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "I_Beta_2_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- I_Beta_3
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "I_Beta_3_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- I_Beta_4
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "I_Beta_4_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- I_Beta_5
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "I_Beta_5_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- I_Beta_6
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "I_Beta_6_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- I_Beta_7
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "I_Beta_7_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- I_Beta_8
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "I_Beta_8_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- I_Beta_9
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "I_Beta_9_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- I_Beta_10
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "I_Beta_10_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
#### Вторая часть ####
II_SymBeta_1 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/II_SymBeta_1.csv")
II_SymBeta_2 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/II_SymBeta_2.csv")
II_SymBeta_3 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/II_SymBeta_3.csv")
II_SymBeta_4 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/II_SymBeta_4.csv")
II_SymBeta_5 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/II_SymBeta_5.csv")
II_SymBeta_6 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/II_SymBeta_6.csv")
II_SymBeta_7 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/II_SymBeta_7.csv")
II_SymBeta_8 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/II_SymBeta_8.csv")
II_SymBeta_9 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/II_SymBeta_9.csv")
II_SymBeta_10 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/II_SymBeta_10.csv")
df <- II_SymBeta_1
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "II_SymBeta_1_results.csv"
write.csv(res_long, file_name, row.names = FALSE)

df <- II_SymBeta_2
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "II_SymBeta_2_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- II_SymBeta_3
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "II_SymBeta_3_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- II_SymBeta_4
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "II_SymBeta_4_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- II_SymBeta_5
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "II_SymBeta_5_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- II_SymBeta_6
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "II_SymBeta_6_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- II_SymBeta_7
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "II_SymBeta_7_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- II_SymBeta_8
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "II_SymBeta_8_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- II_SymBeta_9
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "II_SymBeta_9_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- II_SymBeta_10
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "II_SymBeta_10_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
#### Третья часть ####
III_Gamma_1 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/III_Gamma_1.csv")
III_Gamma_2 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/III_Gamma_2.csv")
III_Gamma_3 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/III_Gamma_3.csv")
III_Gamma_4 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/III_Gamma_4.csv")
III_Gamma_5 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/III_Gamma_5.csv")
III_Gamma_6 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/III_Gamma_6.csv")
III_Gamma_7 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/III_Gamma_7.csv")
III_Gamma_8 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/III_Gamma_8.csv")
III_Gamma_9 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/III_Gamma_9.csv")
III_Gamma_10 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/III_Gamma_10.csv")
df <- III_Gamma_1
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "III_Gamma_1_results.csv"
write.csv(res_long, file_name, row.names = FALSE)

df <- III_Gamma_2
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "III_Gamma_2_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- III_Gamma_3
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "III_Gamma_3_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- III_Gamma_4
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "III_Gamma_4_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- III_Gamma_5
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "III_Gamma_5_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- III_Gamma_6
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "III_Gamma_6_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- III_Gamma_7
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "III_Gamma_7_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- III_Gamma_8
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "III_Gamma_8_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- III_Gamma_9
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "III_Gamma_9_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- III_Gamma_10
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "III_Gamma_10_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
#### Четвертая часть ####
IV_LogNormal_1 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/IV_LogNormal_1.csv")
IV_LogNormal_2 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/IV_LogNormal_2.csv")
IV_LogNormal_3 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/IV_LogNormal_3.csv")
IV_LogNormal_4 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/IV_LogNormal_4.csv")
IV_LogNormal_5 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/IV_LogNormal_5.csv")
IV_LogNormal_6 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/IV_LogNormal_6.csv")
IV_LogNormal_7 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/IV_LogNormal_7.csv")
IV_LogNormal_8 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/IV_LogNormal_8.csv")
IV_LogNormal_9 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/IV_LogNormal_9.csv")
IV_LogNormal_10 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/IV_LogNormal_10.csv")
df <- IV_LogNormal_1
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "IV_LogNormal_1_results.csv"
write.csv(res_long, file_name, row.names = FALSE)

df <- IV_LogNormal_2
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "IV_LogNormal_2_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- IV_LogNormal_3
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "IV_LogNormal_3_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- IV_LogNormal_4
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "IV_LogNormal_4_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- IV_LogNormal_5
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "IV_LogNormal_5_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- IV_LogNormal_6
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "IV_LogNormal_6_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- IV_LogNormal_7
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "IV_LogNormal_7_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- IV_LogNormal_8
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "IV_LogNormal_8_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- IV_LogNormal_9
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "IV_LogNormal_9_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- IV_LogNormal_10
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "IV_LogNormal_10_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
#### Пятая часть ####
V_InvGamma_1 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/V_InvGamma_1.csv")
V_InvGamma_2 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/V_InvGamma_2.csv")
V_InvGamma_3 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/V_InvGamma_3.csv")
V_InvGamma_4 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/V_InvGamma_4.csv")
V_InvGamma_5 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/V_InvGamma_5.csv")
V_InvGamma_6 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/V_InvGamma_6.csv")
V_InvGamma_7 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/V_InvGamma_7.csv")
V_InvGamma_8 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/V_InvGamma_8.csv")
V_InvGamma_9 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/V_InvGamma_9.csv")
V_InvGamma_10 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/V_InvGamma_10.csv")
df <- V_InvGamma_1
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "V_InvGamma_1_results.csv"
write.csv(res_long, file_name, row.names = FALSE)

df <- V_InvGamma_2
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "V_InvGamma_2_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- V_InvGamma_3
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "V_InvGamma_3_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- V_InvGamma_4
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "V_InvGamma_4_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- V_InvGamma_5
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "V_InvGamma_5_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- V_InvGamma_6
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "V_InvGamma_6_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- V_InvGamma_7
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "V_InvGamma_7_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- V_InvGamma_8
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "V_InvGamma_8_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- V_InvGamma_9
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "V_InvGamma_9_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- V_InvGamma_10
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "V_InvGamma_10_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
#### Шестая часть ####
VI_BetaPrime_1 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VI_BetaPrime_1.csv")
VI_BetaPrime_2 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VI_BetaPrime_2.csv")
VI_BetaPrime_3 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VI_BetaPrime_3.csv")
VI_BetaPrime_4 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VI_BetaPrime_4.csv")
VI_BetaPrime_5 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VI_BetaPrime_5.csv")
VI_BetaPrime_6 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VI_BetaPrime_6.csv")
VI_BetaPrime_7 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VI_BetaPrime_7.csv")
VI_BetaPrime_8 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VI_BetaPrime_8.csv")
VI_BetaPrime_9 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VI_BetaPrime_9.csv")
VI_BetaPrime_10 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VI_BetaPrime_10.csv")
df <- VI_BetaPrime_1
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VI_BetaPrime_1_results.csv"
write.csv(res_long, file_name, row.names = FALSE)

df <- VI_BetaPrime_2
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VI_BetaPrime_2_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- VI_BetaPrime_3
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VI_BetaPrime_3_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- VI_BetaPrime_4
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VI_BetaPrime_4_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- VI_BetaPrime_5
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VI_BetaPrime_5_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- VI_BetaPrime_6
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VI_BetaPrime_6_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- VI_BetaPrime_7
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VI_BetaPrime_7_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- VI_BetaPrime_8
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VI_BetaPrime_8_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- VI_BetaPrime_9
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VI_BetaPrime_9_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- VI_BetaPrime_10
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VI_BetaPrime_10_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
#### Седьмая часть ####
VII_t_1 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VII_t_1.csv")
VII_t_2 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VII_t_2.csv")
VII_t_3 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VII_t_3.csv")
VII_t_4 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VII_t_4.csv")
VII_t_5 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VII_t_5.csv")
VII_t_6 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VII_t_6.csv")
VII_t_7 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VII_t_7.csv")
VII_t_8 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VII_t_8.csv")
VII_t_9 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VII_t_9.csv")
VII_t_10 <- read_csv("D:/Наука/Статьи и работы начатые/Отбор признаков/К007 - СПбГУ/Сгенерированные БД/VII_t_10.csv")
df <- VII_t_1
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VII_t_1_results.csv"
write.csv(res_long, file_name, row.names = FALSE)

df <- VII_t_2
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VII_t_2_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- VII_t_3
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VII_t_3_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- VII_t_4
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VII_t_4_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- VII_t_5
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VII_t_5_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- VII_t_6
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VII_t_6_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- VII_t_7
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VII_t_7_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- VII_t_8
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VII_t_8_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- VII_t_9
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VII_t_9_results.csv"
write.csv(res_long, file_name, row.names = FALSE)
df <- VII_t_10
y_name <- "target"
res_long <- run_full_grid(df, y_name,
                          which_metrics = 1:11,
                          rules = 1:5,
                          param_list = list(`1`=3, `2`=0.2, `3`=10, `4`=5, `5`=NULL),
                          q_shadow = 0.95,
                          seed = 1)
file_name <- "VII_t_10_results.csv"
write.csv(res_long, file_name, row.names = FALSE)