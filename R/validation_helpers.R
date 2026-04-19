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

  required_zeus <- unique(c(by, zeus_value_col))
  required_ref <- unique(c(by, reference_value_col))

  missing_zeus <- setdiff(required_zeus, names(zeus_df))
  missing_ref <- setdiff(required_ref, names(reference_df))

  if (length(missing_zeus) > 0L) {
    stop(
      "`zeus_df` is missing required columns: ",
      paste(missing_zeus, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(missing_ref) > 0L) {
    stop(
      "`reference_df` is missing required columns: ",
      paste(missing_ref, collapse = ", "),
      call. = FALSE
    )
  }

  if (!group_col %in% names(zeus_df) && !group_col %in% names(reference_df)) {
    stop("`group_col` was not found in either input data frame.", call. = FALSE)
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

  point_compare$diff <- point_compare$zeus_value - point_compare$reference_value
  point_compare$abs_diff <- abs(point_compare$diff)
  point_compare$sq_diff <- point_compare$diff^2

  summarize_one <- function(dat) {
    data.frame(
      n_points = nrow(dat),
      mean_diff = mean(dat$diff, na.rm = TRUE),
      mean_abs_diff = mean(dat$abs_diff, na.rm = TRUE),
      rmse = sqrt(mean(dat$sq_diff, na.rm = TRUE)),
      max_abs_diff = max(dat$abs_diff, na.rm = TRUE),
      cor_value = stats::cor(
        dat$zeus_value,
        dat$reference_value,
        use = "complete.obs"
      ),
      stringsAsFactors = FALSE
    )
  }

  overall_summary <- summarize_one(point_compare)

  if (!group_col %in% names(point_compare)) {
    group_summary <- overall_summary
  } else {
    split_groups <- split(point_compare, point_compare[[group_col]])

    group_summary <- do.call(
      rbind,
      lapply(names(split_groups), function(g) {
        out <- summarize_one(split_groups[[g]])
        out[[group_col]] <- g
        out
      })
    )

    group_summary <- group_summary[, c(
      group_col,
      "n_points",
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
    group_summary = group_summary
  )
}
