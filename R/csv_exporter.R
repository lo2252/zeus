# Export helpers ----------------------------------------------------------

.zeus_supported_export_names <- function() {
  c(
    "raw",
    "traces_280",
    "traces_70",
    "photocell",
    "stimresp_qc",
    "stimresp_settings",
    "peak_statistics"
  )
}

.zeus_collect_export_components <- function(x, components = NULL) {
  if (!is.list(x)) {
    stop("`x` must be a named list or ZEUS object.", call. = FALSE)
  }

  supported_names <- .zeus_supported_export_names()
  if (!is.null(components)) {
    components <- intersect(as.character(components), supported_names)
    if (length(components) == 0L) {
      stop(
        "`components` must include at least one supported export component: ",
        paste(supported_names, collapse = ", "),
        call. = FALSE
      )
    }
    supported_names <- components
  }

  present_names <- intersect(names(x), supported_names)

  if (!("peak_statistics" %in% present_names)) {
    auto_peak_statistics <- zeus_prepare_peak_statistics_export(x)

    if (!is.null(auto_peak_statistics)) {
      x[["peak_statistics"]] <- auto_peak_statistics
      present_names <- c(present_names, "peak_statistics")
    }
  }

  if (length(present_names) == 0L) {
    stop(
      "`x` does not contain any supported export components: ",
      paste(supported_names, collapse = ", "),
      call. = FALSE
    )
  }

  out <- list()
  for (name in present_names) {
    out[[name]] <- zeus_prepare_csv_export(x[[name]], name)
  }

  out
}

# CSV export --------------------------------------------------------------

#' Export ZEUS data sources as CSV files
#'
#' @description
#' Writes one CSV per ZEUS data source using a shared base path. If `csv_path`
#' ends in `.csv`, that extension is removed before suffixes are added.
#'
#' Supported components are:
#' `raw`, `traces_280`, `traces_70`, `photocell`, `stimresp_qc`,
#' `stimresp_settings`, and `peak_statistics`.
#'
#' For `raw`, objects of class `ABF` or `zeus_abf_raw` are converted to ZEUS
#' long format with [abf_as_df_long()] before export.
#'
#' @param x A named list or ZEUS object containing one or more supported
#'   components.
#' @param csv_path Base path used for the output file names.
#' @param components Optional character vector of supported ZEUS components to
#'   include. Defaults to all supported components present in `x`.
#'
#' @return Invisibly returns a named character vector of written file paths.
#' @export
zeus_export_csv_bundle <- function(x, csv_path, components = NULL) {
  if (!is.character(csv_path) || length(csv_path) != 1L || !nzchar(csv_path)) {
    stop("`csv_path` must be a single non-empty character string.", call. = FALSE)
  }

  base_path <- sub("\\.csv$", "", csv_path, ignore.case = TRUE)
  out_dir <- dirname(base_path)

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  written_paths <- character(0)
  export_components <- .zeus_collect_export_components(x, components = components)

  for (name in names(export_components)) {
    export_df <- export_components[[name]]
    out_file <- paste0(base_path, "_", name, ".csv")
    utils::write.csv(export_df, file = out_file, row.names = FALSE)
    written_paths[name] <- out_file
  }

  invisible(written_paths)
}

# Excel export ------------------------------------------------------------

.zeus_excel_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&apos;", x, fixed = TRUE)
  x
}

.zeus_excel_col_ref <- function(idx) {
  letters_out <- character(0)

  while (idx > 0) {
    remainder <- (idx - 1L) %% 26L
    letters_out <- c(intToUtf8(65L + remainder), letters_out)
    idx <- (idx - 1L) %/% 26L
  }

  paste0(letters_out, collapse = "")
}

.zeus_excel_sheet_names <- function(names_in) {
  cleaned <- gsub("[\\\\/*?:\\[\\]]", "_", names_in)
  cleaned <- substr(cleaned, 1L, 31L)
  cleaned[!nzchar(cleaned)] <- "Sheet"

  out <- character(length(cleaned))
  seen <- character(0)

  for (i in seq_along(cleaned)) {
    candidate <- cleaned[[i]]
    counter <- 1L

    while (candidate %in% seen) {
      suffix <- paste0("_", counter)
      candidate <- paste0(substr(cleaned[[i]], 1L, max(1L, 31L - nchar(suffix))), suffix)
      counter <- counter + 1L
    }

    out[[i]] <- candidate
    seen <- c(seen, candidate)
  }

  out
}

.zeus_excel_cell_xml <- function(value, ref) {
  if (length(value) == 0L || is.na(value)) {
    return(paste0('<c r="', ref, '"/>'))
  }

  if (inherits(value, "POSIXt")) {
    value <- format(value, tz = "UTC", usetz = FALSE)
  }

  if (inherits(value, "Date")) {
    value <- as.character(value)
  }

  if (is.logical(value)) {
    value <- ifelse(value, "TRUE", "FALSE")
  }

  if (is.factor(value)) {
    value <- as.character(value)
  }

  if (is.numeric(value) && is.finite(value)) {
    return(paste0('<c r="', ref, '"><v>', format(value, scientific = FALSE, trim = TRUE), '</v></c>'))
  }

  paste0(
    '<c r="', ref, '" t="inlineStr"><is><t>',
    .zeus_excel_escape(as.character(value)),
    "</t></is></c>"
  )
}

.zeus_excel_sheet_xml <- function(df, sheet_name) {
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  headers <- names(df)

  row_xml <- character(0)

  if (length(headers) > 0L) {
    header_cells <- vapply(
      seq_along(headers),
      function(i) .zeus_excel_cell_xml(headers[[i]], paste0(.zeus_excel_col_ref(i), "1")),
      character(1)
    )
    row_xml <- c(row_xml, paste0('<row r="1">', paste0(header_cells, collapse = ""), "</row>"))
  }

  if (nrow(df) > 0L) {
    for (r in seq_len(nrow(df))) {
      row_idx <- r + 1L
      cells <- vapply(
        seq_along(headers),
        function(i) {
          .zeus_excel_cell_xml(df[[i]][[r]], paste0(.zeus_excel_col_ref(i), row_idx))
        },
        character(1)
      )
      row_xml <- c(row_xml, paste0('<row r="', row_idx, '">', paste0(cells, collapse = ""), "</row>"))
    }
  }

  paste0(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">',
    "<sheetData>",
    paste0(row_xml, collapse = ""),
    "</sheetData>",
    "</worksheet>"
  )
}

.zeus_zip_files <- function(zipfile, files, root_dir = getwd()) {
  old_wd <- setwd(root_dir)
  on.exit(setwd(old_wd), add = TRUE)

  status <- utils::zip(zipfile = zipfile, files = files, flags = "-r9X")

  if (!file.exists(zipfile) || !identical(status, 0L)) {
    stop("Failed to create zip archive: ", zipfile, call. = FALSE)
  }

  invisible(zipfile)
}

.zeus_write_simple_xlsx <- function(sheets, xlsx_path) {
  tmp_dir <- tempfile("zeus_xlsx_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  dir.create(file.path(tmp_dir, "_rels"))
  dir.create(file.path(tmp_dir, "docProps"))
  dir.create(file.path(tmp_dir, "xl"))
  dir.create(file.path(tmp_dir, "xl", "_rels"))
  dir.create(file.path(tmp_dir, "xl", "worksheets"))

  sheet_names <- .zeus_excel_sheet_names(names(sheets))

  worksheets <- vapply(
    seq_along(sheets),
    function(i) .zeus_excel_sheet_xml(sheets[[i]], sheet_names[[i]]),
    character(1)
  )

  writeLines(
    c(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
      '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
      '<Default Extension="xml" ContentType="application/xml"/>',
      '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
      '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
      '<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>',
      '<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>',
      paste0(
        vapply(
          seq_along(worksheets),
          function(i) paste0(
            '<Override PartName="/xl/worksheets/sheet', i,
            '.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
          ),
          character(1)
        ),
        collapse = ""
      ),
      "</Types>"
    ),
    file.path(tmp_dir, "[Content_Types].xml")
  )

  writeLines(
    c(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
      '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>',
      '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" Target="docProps/core.xml"/>',
      '<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" Target="docProps/app.xml"/>',
      "</Relationships>"
    ),
    file.path(tmp_dir, "_rels", ".rels")
  )

  writeLines(
    c(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:dcmitype="http://purl.org/dc/dcmitype/" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">',
      "<dc:creator>ZEUS</dc:creator>",
      "<cp:lastModifiedBy>ZEUS</cp:lastModifiedBy>",
      "</cp:coreProperties>"
    ),
    file.path(tmp_dir, "docProps", "core.xml")
  )

  writeLines(
    c(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties" xmlns:vt="http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes">',
      "<Application>ZEUS</Application>",
      "</Properties>"
    ),
    file.path(tmp_dir, "docProps", "app.xml")
  )

  workbook_sheets <- vapply(
    seq_along(sheet_names),
    function(i) paste0(
      '<sheet name="', .zeus_excel_escape(sheet_names[[i]]),
      '" sheetId="', i, '" r:id="rId', i, '"/>'
    ),
    character(1)
  )

  writeLines(
    c(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
      "<sheets>",
      paste0(workbook_sheets, collapse = ""),
      "</sheets>",
      "</workbook>"
    ),
    file.path(tmp_dir, "xl", "workbook.xml")
  )

  workbook_rels <- vapply(
    seq_along(sheet_names),
    function(i) paste0(
      '<Relationship Id="rId', i,
      '" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet',
      i, '.xml"/>'
    ),
    character(1)
  )

  writeLines(
    c(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
      paste0(workbook_rels, collapse = ""),
      paste0(
        '<Relationship Id="rId', length(sheet_names) + 1L,
        '" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
      ),
      "</Relationships>"
    ),
    file.path(tmp_dir, "xl", "_rels", "workbook.xml.rels")
  )

  writeLines(
    c(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">',
      '<fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>',
      '<fills count="1"><fill><patternFill patternType="none"/></fill></fills>',
      '<borders count="1"><border/></borders>',
      '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>',
      '<cellXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/></cellXfs>',
      "</styleSheet>"
    ),
    file.path(tmp_dir, "xl", "styles.xml")
  )

  for (i in seq_along(worksheets)) {
    writeLines(worksheets[[i]], file.path(tmp_dir, "xl", "worksheets", paste0("sheet", i, ".xml")))
  }

  rel_files <- list.files(tmp_dir, recursive = TRUE, all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  .zeus_zip_files(
    zipfile = normalizePath(xlsx_path, winslash = "/", mustWork = FALSE),
    files = rel_files,
    root_dir = tmp_dir
  )
}

#' Export ZEUS data sources as an Excel workbook
#'
#' @description
#' Writes one worksheet per ZEUS data source into a single `.xlsx` workbook.
#'
#' Supported components are:
#' `raw`, `traces_280`, `traces_70`, `photocell`, `stimresp_qc`,
#' `stimresp_settings`, and `peak_statistics`.
#'
#' @param x A named list or ZEUS object containing one or more supported
#'   components.
#' @param xlsx_path Output path for the `.xlsx` workbook.
#' @param components Optional character vector of supported ZEUS components to
#'   include. Defaults to all supported components present in `x`.
#'
#' @return Invisibly returns `xlsx_path`.
#' @export
zeus_export_excel_workbook <- function(x, xlsx_path, components = NULL) {
  if (!is.character(xlsx_path) || length(xlsx_path) != 1L || !nzchar(xlsx_path)) {
    stop("`xlsx_path` must be a single non-empty character string.", call. = FALSE)
  }

  out_dir <- dirname(xlsx_path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  sheets <- .zeus_collect_export_components(x, components = components)
  .zeus_write_simple_xlsx(sheets = sheets, xlsx_path = xlsx_path)
  invisible(xlsx_path)
}

#' Prepare combined peak statistics for CSV export
#'
#' @param x Object to export.
#'
#' @return A data frame or `NULL` when no statistics export can be derived.
#' @keywords internal
zeus_prepare_peak_statistics_export <- function(x) {
  if (is.list(x) && "peak_statistics" %in% names(x) && is.data.frame(x$peak_statistics)) {
    return(x$peak_statistics)
  }

  if (inherits(x, "zeus_stimresp") || inherits(x, "zeus_abf")) {
    return(zeus_summarize_peak_statistics(x)$combined_export)
  }

  NULL
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
