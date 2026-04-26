# CSV export --------------------------------------------------------------

#' Export ZEUS data sources as CSV files
#'
#' @description
#' Writes one CSV per ZEUS data source using a shared base path. If `csv_path`
#' ends in `.csv`, that extension is removed before suffixes are added.
#'
#' Supported components are:
#' `raw`, `traces_280`, `traces_70`, `photocell`, `stimresp_qc`,
#' and `stimresp_settings`.
#'
#' For `raw`, objects of class `ABF` or `zeus_abf_raw` are converted to ZEUS
#' long format with [abf_as_df_long()] before export.
#'
#' @param x A named list or ZEUS object containing one or more supported
#'   components.
#' @param csv_path Base path used for the output file names.
#'
#' @return Invisibly returns a named character vector of written file paths.
#' @export
zeus_export_csv_bundle <- function(x, csv_path) {
  if (!is.list(x)) {
    stop("`x` must be a named list or ZEUS object.", call. = FALSE)
  }

  if (!is.character(csv_path) || length(csv_path) != 1L || !nzchar(csv_path)) {
    stop("`csv_path` must be a single non-empty character string.", call. = FALSE)
  }

  supported_names <- c(
    "raw",
    "traces_280",
    "traces_70",
    "photocell",
    "stimresp_qc",
    "stimresp_settings"
  )

  present_names <- intersect(names(x), supported_names)

  if (length(present_names) == 0L) {
    stop(
      "`x` does not contain any supported export components: ",
      paste(supported_names, collapse = ", "),
      call. = FALSE
    )
  }

  base_path <- sub("\\.csv$", "", csv_path, ignore.case = TRUE)
  out_dir <- dirname(base_path)

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  written_paths <- character(0)

  for (name in present_names) {
    export_df <- zeus_prepare_csv_export(x[[name]], name)
    out_file <- paste0(base_path, "_", name, ".csv")
    utils::write.csv(export_df, file = out_file, row.names = FALSE)
    written_paths[name] <- out_file
  }

  invisible(written_paths)
}

#' Prepare one ZEUS component for CSV export
#'
#' @param x Object to export.
#' @param name Component name.
#'
#' @return A data frame.
#' @keywords internal
zeus_prepare_csv_export <- function(x, name) {
  if (identical(name, "raw")) {
    if (inherits(x, "zeus_abf_raw")) {
      return(abf_as_df_long(x$raw))
    }

    if (inherits(x, "ABF")) {
      return(abf_as_df_long(x))
    }
  }

  if (identical(name, "stimresp_settings")) {
    return(flatten_zeus_settings(x))
  }

  if (is.null(x)) {
    return(data.frame())
  }

  if (!is.data.frame(x)) {
    stop("Component `", name, "` is not exportable as a data frame.", call. = FALSE)
  }

  x
}

#' Flatten a named settings list for CSV export
#'
#' @param x Named list of settings.
#'
#' @return A one-row data frame.
#' @keywords internal
flatten_zeus_settings <- function(x) {
  if (!is.list(x) || is.null(names(x))) {
    stop("`x` must be a named list.", call. = FALSE)
  }

  cols <- list()

  for (name in names(x)) {
    value <- x[[name]]

    if (length(value) <= 1L) {
      cols[[name]] <- value
    } else {
      value_names <- paste0(name, "_", seq_along(value))
      cols[value_names] <- as.list(value)
    }
  }

  as.data.frame(cols, check.names = FALSE)
}
