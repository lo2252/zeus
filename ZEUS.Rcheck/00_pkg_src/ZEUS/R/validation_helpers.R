# Validation helpers ------------------------------------------------------

#' Summarize protocol-level agreement between ZEUS and a reference export
#'
#' @description
#' Computes match counts and percentages for logical protocol-comparison columns
#' such as `label_match`, `wavelength_match`, and `nd_match`.
#'
#' @param compare_df Data frame containing one or more logical match columns.
#' @param match_cols Character vector of logical columns to summarize.
#'
#' @return A data frame with one row per match column and columns
#'   `comparison`, `n_match`, `n_total`, and `pct_match`.
#' @export
zeus_summarize_protocol_validation <- function(compare_df,
                                               match_cols = c(
                                                 "label_match",
                                                 "wavelength_match",
                                                 "nd_match",
                                                 "match"
                                               )) {
  if (!is.data.frame(compare_df)) {
    stop("`compare_df` must be a data frame.", call. = FALSE)
  }

  match_cols <- intersect(match_cols, names(compare_df))

  if (length(match_cols) == 0L) {
    stop("No requested `match_cols` were found in `compare_df`.", call. = FALSE)
  }

  out <- lapply(match_cols, function(col) {
    vals <- compare_df[[col]]

    if (!is.logical(vals)) {
      stop("Column `", col, "` must be logical.", call. = FALSE)
    }

    n_total <- sum(!is.na(vals))
    n_match <- sum(vals, na.rm = TRUE)

    data.frame(
      comparison = col,
      n_match = n_match,
      n_total = n_total,
      pct_match = if (n_total > 0L) 100 * n_match / n_total else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

#' Validate protocol-level agreement
#'
#' @description
#' Checks that all non-missing protocol comparison values pass. This is a small
#' assertion helper for validation scripts and package tests that need a clear
#' pass/fail result instead of only printed diagnostics.
#'
#' @param compare_df Data frame containing logical match columns.
#' @param match_cols Character vector of logical columns to validate.
#' @param error Logical; if `TRUE`, fail with an error when any comparison does
#'   not pass.
#'
#' @return The summary data frame produced by
#'   `zeus_summarize_protocol_validation()`, with an added logical `passed`
#'   column.
#' @export
zeus_validate_protocol_agreement <- function(compare_df,
                                             match_cols = c(
                                               "label_match",
                                               "wavelength_match",
                                               "nd_match",
                                               "match"
                                             ),
                                             error = TRUE) {
  out <- zeus_summarize_protocol_validation(
    compare_df = compare_df,
    match_cols = match_cols
  )

  out$passed <- out$n_total > 0L & out$n_match == out$n_total

  if (isTRUE(error) && any(!out$passed)) {
    failed <- out$comparison[!out$passed]
    stop(
      "Protocol validation failed for: ",
      paste(failed, collapse = ", "),
      call. = FALSE
    )
  }

  out
}

.zeus_require_columns <- function(df, cols, arg_name) {
  missing_cols <- setdiff(cols, names(df))

  if (length(missing_cols) > 0L) {
    stop(
      "`", arg_name, "` is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
}

.zeus_assert_numeric_column <- function(df, col, arg_name) {
  if (!is.numeric(df[[col]])) {
    stop(
      "`", arg_name, "$", col, "` must be numeric.",
      call. = FALSE
    )
  }
}

.zeus_has_duplicate_keys <- function(df, by) {
  any(duplicated(df[by]))
}

.zeus_duplicate_key_message <- function(df, by) {
  dup <- df[duplicated(df[by]) | duplicated(df[by], fromLast = TRUE), by, drop = FALSE]
  dup <- utils::head(unique(dup), 5L)

  rows <- apply(dup, 1L, function(x) {
    paste(paste(by, x, sep = "="), collapse = ", ")
  })

  paste(rows, collapse = "; ")
}

.zeus_cor_complete <- function(x, y) {
  ok <- stats::complete.cases(x, y)

  if (sum(ok) < 2L) {
    return(NA_real_)
  }

  suppressWarnings(stats::cor(x[ok], y[ok]))
}

.zeus_agreement_summary <- function(dat) {
  finite_pairs <- is.finite(dat$zeus_value) & is.finite(dat$reference_value)

  data.frame(
    n_points = nrow(dat),
    n_complete = sum(finite_pairs),
    mean_diff = mean(dat$diff, na.rm = TRUE),
    mean_abs_diff = mean(dat$abs_diff, na.rm = TRUE),
    rmse = sqrt(mean(dat$sq_diff, na.rm = TRUE)),
    max_abs_diff = if (all(is.na(dat$abs_diff))) NA_real_ else max(dat$abs_diff, na.rm = TRUE),
    cor_value = .zeus_cor_complete(dat$zeus_value, dat$reference_value),
    stringsAsFactors = FALSE
  )
}

#' Compare ZEUS and reference waveforms point-by-point
#'
#' @description
#' Joins ZEUS and reference waveform tables by shared key columns and computes
#' point-level differences together with overall and group-level agreement
#' summaries.
#'
#' @param zeus_df Data frame containing ZEUS waveform values.
#' @param reference_df Data frame containing reference waveform values, such as
#'   an Origin export reshaped to long format.
#' @param by Character vector of join columns. Defaults to `stim_label` and
#'   `time_ms`.
#' @param zeus_value_col Name of the ZEUS value column.
#' @param reference_value_col Name of the reference value column.
#' @param group_col Optional grouping column for group-level summaries. Default
#'   is `stim_label`.
#'
#' @return A named list with:
#' \describe{
#'   \item{point_compare}{Joined point-level comparison table.}
#'   \item{overall_summary}{One-row data frame of overall agreement metrics.}
#'   \item{group_summary}{Group-level agreement metrics.}
#'   \item{coverage_summary}{One-row data frame describing join coverage.}
#' }
#' @export
zeus_compare_waveforms <- function(zeus_df,
                                   reference_df,
                                   by = c("stim_label", "time_ms"),
                                   zeus_value_col = "zeus_value",
                                   reference_value_col = "reference_value",
                                   group_col = "stim_label") {
  if (!is.data.frame(zeus_df) || !is.data.frame(reference_df)) {
    stop("`zeus_df` and `reference_df` must both be data frames.", call. = FALSE)
  }

  if (!is.character(by) || length(by) == 0L || any(!nzchar(by))) {
    stop("`by` must be a non-empty character vector.", call. = FALSE)
  }

  required_zeus <- unique(c(by, zeus_value_col))
  required_ref <- unique(c(by, reference_value_col))

  .zeus_require_columns(zeus_df, required_zeus, "zeus_df")
  .zeus_require_columns(reference_df, required_ref, "reference_df")
  .zeus_assert_numeric_column(zeus_df, zeus_value_col, "zeus_df")
  .zeus_assert_numeric_column(reference_df, reference_value_col, "reference_df")

  if (!is.null(group_col) &&
      !group_col %in% names(zeus_df) &&
      !group_col %in% names(reference_df)) {
    stop("`group_col` was not found in either input data frame.", call. = FALSE)
  }

  if (.zeus_has_duplicate_keys(zeus_df, by)) {
    stop(
      "`zeus_df` contains duplicate join keys. Examples: ",
      .zeus_duplicate_key_message(zeus_df, by),
      call. = FALSE
    )
  }

  if (.zeus_has_duplicate_keys(reference_df, by)) {
    stop(
      "`reference_df` contains duplicate join keys. Examples: ",
      .zeus_duplicate_key_message(reference_df, by),
      call. = FALSE
    )
  }

  zeus_work <- zeus_df
  ref_work <- reference_df

  names(zeus_work)[names(zeus_work) == zeus_value_col] <- "zeus_value"
  names(ref_work)[names(ref_work) == reference_value_col] <- "reference_value"

  point_compare <- merge(
    zeus_work,
    ref_work,
    by = by,
    all = FALSE,
    sort = FALSE
  )

  if (nrow(point_compare) == 0L) {
    stop("No matched rows were found after joining the waveform tables.", call. = FALSE)
  }

  coverage_summary <- data.frame(
    n_zeus_rows = nrow(zeus_work),
    n_reference_rows = nrow(ref_work),
    n_matched_rows = nrow(point_compare),
    pct_zeus_matched = 100 * nrow(point_compare) / nrow(zeus_work),
    pct_reference_matched = 100 * nrow(point_compare) / nrow(ref_work),
    stringsAsFactors = FALSE
  )

  point_compare$diff <- point_compare$zeus_value - point_compare$reference_value
  point_compare$abs_diff <- abs(point_compare$diff)
  point_compare$sq_diff <- point_compare$diff^2

  overall_summary <- .zeus_agreement_summary(point_compare)

  if (is.null(group_col) || !group_col %in% names(point_compare)) {
    group_summary <- overall_summary
  } else {
    split_groups <- split(point_compare, point_compare[[group_col]])

    group_summary <- do.call(
      rbind,
      lapply(names(split_groups), function(g) {
        out <- .zeus_agreement_summary(split_groups[[g]])
        out[[group_col]] <- g
        out
      })
    )

    group_summary <- group_summary[, c(
      group_col,
      "n_points",
      "n_complete",
      "mean_diff",
      "mean_abs_diff",
      "rmse",
      "max_abs_diff",
      "cor_value"
    )]
    rownames(group_summary) <- NULL
  }

  list(
    point_compare = point_compare,
    overall_summary = overall_summary,
    group_summary = group_summary,
    coverage_summary = coverage_summary
  )
}

#' Validate response agreement metrics
#'
#' @description
#' Applies package-level response validation thresholds to one or more agreement
#' summary rows, such as the `overall_summary` or `group_summary` returned by
#' `zeus_compare_waveforms()`.
#'
#' @param summary_df Data frame with agreement columns `n_points`,
#'   `n_complete`, `mean_abs_diff`, `rmse`, `max_abs_diff`, and `cor_value`.
#' @param min_n_complete Optional minimum complete point count.
#' @param min_correlation Optional minimum Pearson correlation.
#' @param max_mean_abs_diff Optional maximum mean absolute difference.
#' @param max_rmse Optional maximum root-mean-square error.
#' @param max_abs_diff Optional maximum absolute point difference.
#' @param label Human-readable label used in error messages.
#' @param error Logical; if `TRUE`, fail with an error when any check does not
#'   pass.
#'
#' @return A data frame containing one row per threshold check and a logical
#'   `passed` column.
#' @export
zeus_validate_response_agreement <- function(summary_df,
                                             min_n_complete = NULL,
                                             min_correlation = NULL,
                                             max_mean_abs_diff = NULL,
                                             max_rmse = NULL,
                                             max_abs_diff = NULL,
                                             label = "response",
                                             error = TRUE) {
  if (!is.data.frame(summary_df)) {
    stop("`summary_df` must be a data frame.", call. = FALSE)
  }

  required <- c(
    "n_points",
    "mean_abs_diff",
    "rmse",
    "max_abs_diff",
    "cor_value"
  )
  .zeus_require_columns(summary_df, required, "summary_df")

  if (!"n_complete" %in% names(summary_df)) {
    summary_df$n_complete <- summary_df$n_points
  }

  checks <- list()

  finite_min <- function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }

    min(x, na.rm = TRUE)
  }

  finite_max <- function(x) {
    if (all(is.na(x))) {
      return(NA_real_)
    }

    max(x, na.rm = TRUE)
  }

  all_at_least <- function(x, threshold) {
    all(!is.na(x)) && all(x >= threshold)
  }

  all_at_most <- function(x, threshold) {
    all(!is.na(x)) && all(x <= threshold)
  }

  add_check <- function(metric, observed, threshold, passed) {
    checks[[length(checks) + 1L]] <<- data.frame(
      label = label,
      metric = metric,
      observed = observed,
      threshold = threshold,
      passed = passed,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(min_n_complete)) {
    add_check(
      "n_complete",
      finite_min(summary_df$n_complete),
      min_n_complete,
      all_at_least(summary_df$n_complete, min_n_complete)
    )
  }

  if (!is.null(min_correlation)) {
    add_check(
      "cor_value",
      finite_min(summary_df$cor_value),
      min_correlation,
      all_at_least(summary_df$cor_value, min_correlation)
    )
  }

  if (!is.null(max_mean_abs_diff)) {
    add_check(
      "mean_abs_diff",
      finite_max(summary_df$mean_abs_diff),
      max_mean_abs_diff,
      all_at_most(summary_df$mean_abs_diff, max_mean_abs_diff)
    )
  }

  if (!is.null(max_rmse)) {
    add_check(
      "rmse",
      finite_max(summary_df$rmse),
      max_rmse,
      all_at_most(summary_df$rmse, max_rmse)
    )
  }

  if (!is.null(max_abs_diff)) {
    add_check(
      "max_abs_diff",
      finite_max(summary_df$max_abs_diff),
      max_abs_diff,
      all_at_most(summary_df$max_abs_diff, max_abs_diff)
    )
  }

  if (length(checks) == 0L) {
    stop("At least one validation threshold must be supplied.", call. = FALSE)
  }

  out <- do.call(rbind, checks)

  if (isTRUE(error) && any(!out$passed)) {
    failed <- out$metric[!out$passed]
    stop(
      "Response validation failed for `",
      label,
      "`: ",
      paste(failed, collapse = ", "),
      call. = FALSE
    )
  }

  out
}
