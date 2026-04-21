# ZEUS Shiny App ---------------------------------------------------------

.zeus_app_www_dir <- function() {
  candidates <- c()

  if (requireNamespace("pkgload", quietly = TRUE)) {
    candidates <- c(
      candidates,
      tryCatch(pkgload::pkg_path("app", "www"), error = function(e) ""),
      tryCatch(pkgload::pkg_path("inst", "app", "www"), error = function(e) "")
    )
  }

  candidates <- c(
    candidates,
    system.file("app/www", package = "ZEUS")
  )

  candidates <- unique(stats::na.omit(candidates[nzchar(candidates)]))
  existing <- candidates[dir.exists(candidates)]

  if (length(existing) == 0L) {
    ""
  } else {
    existing[[1]]
  }
}

.zeus_app_logo_src <- function() {
  www_dir <- .zeus_app_www_dir()

  if (!nzchar(www_dir)) {
    return(NULL)
  }

  logo_path <- file.path(www_dir, "ZEUS_Logo6.png")
  if (!file.exists(logo_path)) {
    return(NULL)
  }

  "zeus-www/ZEUS_Logo6.png"
}

.zeus_app_css_path <- function() {
  www_dir <- .zeus_app_www_dir()

  if (!nzchar(www_dir)) {
    return(NULL)
  }

  css_path <- file.path(www_dir, "zeus-app.css")
  if (!file.exists(css_path)) {
    return(NULL)
  }

  css_path
}

.zeus_app_available_wavelengths <- function(x) {
  if (is.null(x) || is.null(x$traces_70) || !is.data.frame(x$traces_70)) {
    return(character(0))
  }

  stim_label <- x$traces_70$stim_label
  if (is.null(stim_label)) {
    return(character(0))
  }

  wavelength <- stringr::str_extract(as.character(stim_label), "^\\S+")
  wavelength <- unique(stats::na.omit(wavelength[nzchar(wavelength)]))
  wavelength <- wavelength[!identical(wavelength, "White") & wavelength != "White"]

  if (length(wavelength) == 0L) {
    character(0)
  } else {
    wavelength
  }
}

.zeus_app_peak_summary_table <- function(stats_list) {
  if (is.null(stats_list) || is.null(stats_list$key_statistics)) {
    return(data.frame())
  }

  stats_list$key_statistics |>
    dplyr::mutate(
      mean_peak_mv = round(.data$mean_peak_mv, 3),
      sd_peak_mv = round(.data$sd_peak_mv, 3),
      sem_peak_mv = round(.data$sem_peak_mv, 3),
      mean_latency_ms = round(.data$mean_latency_ms, 3),
      sd_latency_ms = round(.data$sd_latency_ms, 3)
    ) |>
    dplyr::rename(
      Protocol = protocol_id,
      `Wave Type` = peak_type,
      `Sample Count` = n,
      `Mean Peak (mV)` = mean_peak_mv,
      `Peak SD (mV)` = sd_peak_mv,
      `Peak SEM (mV)` = sem_peak_mv,
      `Median Peak (mV)` = median_peak_mv,
      `Minimum Peak (mV)` = min_peak_mv,
      `Maximum Peak (mV)` = max_peak_mv,
      `Mean Latency (ms)` = mean_latency_ms,
      `Latency SD (ms)` = sd_latency_ms
    )
}

.zeus_app_export_paths <- function(directory, file_stem) {
  if (!is.character(directory) || length(directory) != 1L || !nzchar(directory)) {
    stop("`directory` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.character(file_stem) || length(file_stem) != 1L || !nzchar(file_stem)) {
    stop("`file_stem` must be a single non-empty character string.", call. = FALSE)
  }

  list(
    mean_plot = file.path(directory, paste0(file_stem, "_mean_waveform.png")),
    spectral_plot = file.path(directory, paste0(file_stem, "_spectral_waveform.png")),
    intensity_plot = file.path(directory, paste0(file_stem, "_intensity_response.png")),
    csv_base = file.path(directory, file_stem)
  )
}

.zeus_app_export_stem <- function(input, item = NULL) {
  file_stem <- trimws(.zeus_app_if_null(input$export_name, ""))

  if (!nzchar(file_stem) && !is.null(item) && !is.null(item$file_name)) {
    file_stem <- tools::file_path_sans_ext(item$file_name)
  }

  if (!nzchar(file_stem)) {
    file_stem <- "zeus_export"
  }

  .zeus_app_sanitize_stem(file_stem)
}

.zeus_app_write_export_bundle <- function(bundle_dir,
                                          file_stem,
                                          mean_plot,
                                          spectral_plot,
                                          intensity_plot,
                                          x) {
  dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)

  export_paths <- .zeus_app_export_paths(bundle_dir, file_stem)

  ggplot2::ggsave(export_paths$mean_plot, plot = mean_plot, width = 11, height = 6.5, dpi = 320, bg = "white")
  ggplot2::ggsave(export_paths$spectral_plot, plot = spectral_plot, width = 14, height = 9, dpi = 320, bg = "white")
  ggplot2::ggsave(export_paths$intensity_plot, plot = intensity_plot, width = 11, height = 6.5, dpi = 320, bg = "white")
  zeus_export_csv_bundle(x, export_paths$csv_base)

  invisible(export_paths)
}

.zeus_app_placeholder_plot <- function(message) {
  ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_text(
      label = message,
      colour = "#1F2A33",
      family = "Georgia",
      size = 5
    ) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_void()
}

.zeus_app_import_settings <- function(input) {
  list(
    erg_channel = input$erg_channel,
    pc_channel = input$pc_channel,
    zero_baseline = isTRUE(input$zero_baseline),
    smooth_n = as.integer(input$smooth_n),
    align_to_stimulus = input$align_to_stimulus,
    stimresp_zero_baseline = isTRUE(input$stimresp_zero_baseline),
    stimresp_baseline_window_ms = c(input$stimresp_baseline_start, input$stimresp_baseline_end),
    stimresp_exclude_noisy = isTRUE(input$stimresp_exclude_noisy),
    stimresp_noise_window_ms = c(input$stimresp_noise_start, input$stimresp_noise_end),
    stimresp_noise_ratio_cutoff = input$stimresp_noise_ratio_cutoff,
    stimresp_runmean_k = as.integer(input$stimresp_runmean_k),
    stimresp_sg_smooth = isTRUE(input$stimresp_sg_smooth),
    stimresp_sg_n = as.integer(input$stimresp_sg_n),
    stimresp_sg_p = as.integer(input$stimresp_sg_p)
  )
}

.zeus_app_attach_source_info <- function(x, file_name = NULL) {
  if (is.null(x) || !is.list(x)) {
    return(x)
  }

  x$source_file_name <- file_name
  x
}

.zeus_app_can_compare_raw <- function(x, data_slot = "traces_70") {
  if (is.null(x) || !is.list(x) || !data_slot %in% names(x)) {
    return(FALSE)
  }

  df <- x[[data_slot]]
  is.data.frame(df) && "value_raw" %in% names(df)
}

.zeus_app_available_slots <- function(x) {
  if (is.null(x) || !is.list(x)) {
    return("traces_70")
  }

  out <- c()
  if (is.data.frame(x$traces_70)) {
    out <- c(out, "traces_70")
  }
  if (is.data.frame(x$traces_280)) {
    out <- c(out, "traces_280")
  }

  if (length(out) == 0L) {
    "traces_70"
  } else {
    out
  }
}

.zeus_app_sanitize_stem <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)

  if (!nzchar(x)) {
    "file"
  } else {
    x
  }
}

.zeus_app_if_null <- function(x, y) {
  if (is.null(x)) y else x
}

.zeus_app_default_protocols <- function(n_files, analysis_mode = "single") {
  if (identical(analysis_mode, "full_fish")) {
    defaults <- c("C1", rep("C0", max(0L, n_files - 1L)))
    defaults[seq_len(n_files)]
  } else {
    rep("C1", n_files)
  }
}

.zeus_app_collect_protocols <- function(input, files, analysis_mode = "single") {
  if (is.null(files) || nrow(files) == 0L) {
    return(character(0))
  }

  if (!identical(analysis_mode, "full_fish")) {
    return(rep(.zeus_app_if_null(input$protocol_single, "C1"), nrow(files)))
  }

  defaults <- .zeus_app_default_protocols(nrow(files), analysis_mode = "full_fish")
  vapply(
    seq_len(nrow(files)),
    function(i) {
      .zeus_app_if_null(input[[paste0("protocol_file_", i)]], defaults[[i]])
    },
    character(1)
  )
}

.zeus_app_validate_protocols <- function(protocols, analysis_mode = "single") {
  if (!identical(analysis_mode, "full_fish")) {
    return(NULL)
  }

  if (length(protocols) != 3L) {
    return("Full Fish Analysis requires exactly 3 ABF files.")
  }

  protocol_counts <- table(protocols)
  c1_n <- if ("C1" %in% names(protocol_counts)) unname(protocol_counts[["C1"]]) else 0L
  c0_n <- if ("C0" %in% names(protocol_counts)) unname(protocol_counts[["C0"]]) else 0L

  if (!identical(c1_n, 1L) || !identical(c0_n, 2L)) {
    return("Full Fish Analysis requires exactly 1 C1 file and 2 C0 files.")
  }

  NULL
}

.zeus_app_make_item_labels <- function(protocols) {
  protocol_index <- ave(seq_along(protocols), protocols, FUN = seq_along)
  protocol_n <- ave(seq_along(protocols), protocols, FUN = length)

  vapply(
    seq_along(protocols),
    function(i) {
      if (protocol_n[[i]] <= 1L) {
        protocols[[i]]
      } else {
        paste(protocols[[i]], LETTERS[[protocol_index[[i]]]])
      }
    },
    character(1)
  )
}

.zeus_app_find_item <- function(bundle, id) {
  if (is.null(bundle) || is.null(bundle$items) || length(bundle$items) == 0L) {
    return(NULL)
  }

  idx <- match(id, vapply(bundle$items, `[[`, character(1), "id"))
  if (is.na(idx)) {
    return(NULL)
  }

  bundle$items[[idx]]
}

.zeus_app_choose_directory <- function(current = "") {
  current <- if (is.character(current) && length(current) == 1L && nzchar(current)) {
    current
  } else {
    path.expand("~")
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && isTRUE(rstudioapi::isAvailable())) {
    out <- tryCatch(
      rstudioapi::selectDirectory(caption = "Select export directory", path = current),
      error = function(e) ""
    )

    if (is.character(out) && length(out) == 1L && nzchar(out)) {
      return(out)
    }
  }

  sysname <- Sys.info()[["sysname"]]

  if (identical(sysname, "Darwin")) {
    out <- tryCatch(
      system2(
        "osascript",
        c(
          "-l", "JavaScript",
          "-e",
          paste0(
            "var app = Application.currentApplication();",
            "app.includeStandardAdditions = true;",
            "app.chooseFolder({withPrompt: 'Select export directory', defaultLocation: Path('",
            normalizePath(current, winslash = "/", mustWork = FALSE),
            "')}).toString();"
          )
        ),
        stdout = TRUE,
        stderr = FALSE
      ),
      error = function(e) character(0)
    )

    out <- trimws(out)
    if (length(out) > 0L && nzchar(out[[1]])) {
      return(out[[1]])
    }

    return("")
  }

  if (identical(sysname, "Windows")) {
    out <- tryCatch(utils::choose.dir(default = current, caption = "Select export directory"),
      error = function(e) character(0)
    )
    if (length(out) > 0L && nzchar(out[[1]])) {
      return(out[[1]])
    }

    return("")
  }

  if (nzchar(Sys.which("zenity"))) {
    out <- tryCatch(
      system2("zenity", c("--file-selection", "--directory", "--filename", current), stdout = TRUE, stderr = FALSE),
      error = function(e) character(0)
    )
    out <- trimws(out)
    if (length(out) > 0L && nzchar(out[[1]])) {
      return(out[[1]])
    }
  }

  ""
}

.zeus_app_peak_settings <- function(input) {
  list(
    baseline_window_ms = c(input$peak_baseline_start, input$peak_baseline_end),
    response_window_ms = c(input$peak_response_start, input$peak_response_end),
    awave_window_ms = c(input$awave_window_start, input$awave_window_end),
    stimulus_onset_ms = input$stimulus_onset_ms,
    same_sign = isTRUE(input$same_sign),
    time_reference = input$time_reference,
    dwave_window_ms = c(input$dwave_window_start, input$dwave_window_end)
  )
}

.zeus_app_theme <- function() {
  bslib::bs_theme(
    version = 5,
    bg = "#EAF7FB",
    fg = "#12344C",
    primary = "#157FA8",
    secondary = "#2F5D76",
    success = "#2F8C8C",
    info = "#5BB7D4",
    base_font = "Georgia",
    heading_font = "Palatino Linotype",
    code_font = "Courier New"
  )
}

.zeus_app_table <- function(df) {
  DT::datatable(
    df,
    rownames = FALSE,
    extensions = "Buttons",
    options = list(
      dom = "Bfrtip",
      buttons = c("copy", "csv"),
      pageLength = 8,
      autoWidth = TRUE,
      scrollX = TRUE
    ),
    class = "compact stripe hover"
  )
}

.zeus_app_sidebar <- function() {
  shiny::tagList(
    shiny::div(
      class = "zeus-sidebar-block",
      shiny::radioButtons(
        "analysis_mode",
        "Analysis mode",
        choices = c("Single File" = "single", "Full Fish Analysis" = "full_fish"),
        selected = "single",
        inline = TRUE
      ),
      shiny::fileInput(
        "abf_file",
        "ABF File",
        accept = c(".abf"),
        multiple = TRUE
      ),
      shiny::uiOutput("protocol_assignment_ui"),
      shiny::actionButton("import_file", "Import File", class = "btn-zeus-primary")
    ),
    bslib::accordion(
      id = "zeus_controls",
      open = c("Peak Settings", "Import Settings", "Plot Settings", "Export"),
      bslib::accordion_panel(
        "Peak Settings",
        shiny::selectInput(
          "time_reference",
          "Time reference",
          choices = c("absolute", "stimulus"),
          selected = "absolute"
        ),
        shiny::numericInput("stimulus_onset_ms", "Stimulus onset (ms)", value = 400, step = 1),
        shiny::checkboxInput("same_sign", "Force consistent B-wave sign", value = TRUE),
        shiny::fluidRow(
          shiny::column(6, shiny::numericInput("peak_baseline_start", "Baseline start (ms)", value = 300)),
          shiny::column(6, shiny::numericInput("peak_baseline_end", "Baseline end (ms)", value = 400))
        ),
        shiny::fluidRow(
          shiny::column(6, shiny::numericInput("peak_response_start", "Response start (ms)", value = 400)),
          shiny::column(6, shiny::numericInput("peak_response_end", "Response end (ms)", value = 700))
        ),
        shiny::fluidRow(
          shiny::column(6, shiny::numericInput("awave_window_start", "A-wave start (ms)", value = 400)),
          shiny::column(6, shiny::numericInput("awave_window_end", "A-wave end (ms)", value = 700))
        ),
        shiny::fluidRow(
          shiny::column(6, shiny::numericInput("dwave_window_start", "D-wave start (ms)", value = 700)),
          shiny::column(6, shiny::numericInput("dwave_window_end", "D-wave end (ms)", value = 1000))
        )
      ),
      bslib::accordion_panel(
        "Import Settings",
        shiny::textInput("erg_channel", "ERG Channel", value = "ERG DAM80"),
        shiny::textInput("pc_channel", "Photocell Channel", value = "Photocell"),
        shiny::checkboxInput("zero_baseline", "Zero raw baseline", value = FALSE),
        shiny::numericInput("smooth_n", "Raw smoothing window", value = 1, min = 1, step = 1),
        shiny::selectInput(
          "align_to_stimulus",
          "Align to stimulus",
          choices = c("protocol", "photocell"),
          selected = "protocol"
        ),
        shiny::checkboxInput("stimresp_zero_baseline", "Zero StimResp baseline", value = TRUE),
        shiny::fluidRow(
          shiny::column(6, shiny::numericInput("stimresp_baseline_start", "Baseline start (ms)", value = 300)),
          shiny::column(6, shiny::numericInput("stimresp_baseline_end", "Baseline end (ms)", value = 400))
        ),
        shiny::checkboxInput("stimresp_exclude_noisy", "Exclude noisy sweeps", value = TRUE),
        shiny::fluidRow(
          shiny::column(6, shiny::numericInput("stimresp_noise_start", "Noise start (ms)", value = 300)),
          shiny::column(6, shiny::numericInput("stimresp_noise_end", "Noise end (ms)", value = 1000))
        ),
        shiny::numericInput("stimresp_noise_ratio_cutoff", "Noise ratio cutoff", value = 1.5, min = 0, step = 0.1),
        shiny::numericInput("stimresp_runmean_k", "StimResp running mean", value = 16, min = 1, step = 1),
        shiny::checkboxInput("stimresp_sg_smooth", "Savitzky-Golay smoothing", value = TRUE),
        shiny::fluidRow(
          shiny::column(6, shiny::numericInput("stimresp_sg_n", "SG window", value = 101, min = 3, step = 2)),
          shiny::column(6, shiny::numericInput("stimresp_sg_p", "SG polynomial", value = 2, min = 1, step = 1))
        )
      ),
      bslib::accordion_panel(
        "Plot Settings",
        shiny::uiOutput("plot_settings_ui")
      ),
      bslib::accordion_panel(
        "Export",
        shiny::textInput("export_name", "Export file name", value = ""),
        shiny::downloadButton("download_csv_bundle", "Download CSV Bundle", class = "btn-zeus-secondary"),
        shiny::tags$div(
          style = "font-size: 0.9rem; color: #6b7280; margin-top: 0.5rem;",
          "Your browser will ask where to save the export zip file."
        ),
        shiny::tags$hr(),
        shiny::actionButton("close_app", "Close App", class = "btn-zeus-secondary")
      )
    )
  )
}

#' Build the ZEUS analysis application
#'
#' Creates the ZEUS Shiny application for importing ABF files, reviewing
#' waveform plots, inspecting peak statistics, and exporting plots plus CSV
#' outputs.
#'
#' @return A `shiny.appobj`.
#' @export
zeus_app <- function() {
  needed <- c("shiny", "bslib", "DT")
  missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_pkgs) > 0L) {
    stop(
      "The ZEUS app requires these packages: ",
      paste(missing_pkgs, collapse = ", "),
      call. = FALSE
    )
  }

  shiny::shinyApp(
    ui = {
      www_dir <- .zeus_app_www_dir()
      if (nzchar(www_dir)) {
        shiny::addResourcePath("zeus-www", www_dir)
      }

      css_path <- .zeus_app_css_path()
      logo_src <- .zeus_app_logo_src()

      shiny::fluidPage(
        theme = .zeus_app_theme(),
        shiny::tags$head(
          if (!is.null(css_path)) shiny::includeCSS(css_path)
        ),
        shiny::div(
          class = "zeus-app-shell",
          shiny::div(
              class = "zeus-hero",
              shiny::div(
                class = "zeus-hero-brand",
                if (!is.null(logo_src)) shiny::tags$img(src = logo_src, alt = "ZEUS logo", class = "zeus-logo"),
                shiny::div(
                  class = "zeus-hero-copy",
                  shiny::h1("ZEUS"),
                  shiny::p("Import ERG recordings, review waveform panels, inspect peak statistics, and export report-ready outputs.")
                )
            ),
            shiny::uiOutput("import_summary")
          ),
          bslib::layout_sidebar(
            sidebar = bslib::sidebar(
              width = 360,
              open = "desktop",
              .zeus_app_sidebar()
            ),
            shiny::tagList(
              shiny::uiOutput("selected_file_ui"),
              shiny::tabsetPanel(
                id = "zeus_main_tabs",
                type = "tabs",
                shiny::tabPanel(
                  "Home",
                  shiny::uiOutput("home_tab_ui")
                ),
                shiny::tabPanel(
                  "Spectral Waveform",
                  shiny::uiOutput("spectral_tab_ui")
                ),
                shiny::tabPanel(
                  "Intensity Response",
                  shiny::uiOutput("intensity_tab_ui")
                ),
                shiny::tabPanel(
                  "Peak Statistics",
                  shiny::uiOutput("statistics_tab_ui")
                )
              )
            )
          )
        )
      )
    },
    server = function(input, output, session) {
      imported_data <- shiny::reactiveVal(NULL)
      import_status <- shiny::reactiveVal("Waiting for file import.")
      
      output$protocol_assignment_ui <- shiny::renderUI({
        files <- input$abf_file
        mode <- .zeus_app_if_null(input$analysis_mode, "single")

        if (identical(mode, "single")) {
          return(
            shiny::selectInput(
              "protocol_single",
              "Protocol",
              choices = c("C0", "C1"),
              selected = .zeus_app_if_null(input$protocol_single, "C1")
            )
          )
        }

        if (is.null(files) || nrow(files) == 0L) {
          return(
            shiny::tags$div(
              style = "font-size: 0.92rem; color: #5b6b75;",
              "Full Fish Analysis requires 3 files assigned as 1 C1 and 2 C0."
            )
          )
        }

        defaults <- .zeus_app_default_protocols(nrow(files), analysis_mode = "full_fish")
        shiny::tagList(
          shiny::tags$div(
            style = "font-size: 0.92rem; color: #5b6b75; margin-bottom: 0.5rem;",
            "Assign 1 C1 and 2 C0 files for each fish."
          ),
          lapply(seq_len(nrow(files)), function(i) {
            shiny::selectInput(
              inputId = paste0("protocol_file_", i),
              label = paste("Protocol for", files$name[[i]]),
              choices = c("C0", "C1"),
              selected = .zeus_app_if_null(input[[paste0("protocol_file_", i)]], defaults[[i]])
            )
          })
        )
      })

      output$plot_settings_ui <- shiny::renderUI({
        bundle <- imported_data()

        if (is.null(bundle) || length(bundle$items) == 0L) {
          return(
            shiny::tags$div(
              style = "font-size: 0.95rem; color: #6b7280;",
              "Import a file to load plot options."
            )
          )
        }

        available_slots <- unique(unlist(lapply(bundle$items, function(item) {
          .zeus_app_available_slots(item$data)
        })))

        selected_slot <- if (!is.null(input$mean_data_slot) && input$mean_data_slot %in% available_slots) {
          input$mean_data_slot
        } else {
          available_slots[[1]]
        }
        can_compare <- any(vapply(bundle$items, function(item) {
          .zeus_app_can_compare_raw(item$data, selected_slot)
        }, logical(1)))

        shiny::tagList(
          shiny::selectInput(
            "mean_data_slot",
            "Mean waveform source",
            choices = stats::setNames(available_slots, available_slots),
            selected = selected_slot
          ),
          shiny::checkboxInput(
            "compare_raw",
            "Compare raw and smoothed mean traces",
            value = isTRUE(input$compare_raw) && can_compare
          ),
          if (!can_compare) {
            shiny::tags$div(
              style = "font-size: 0.9rem; color: #6b7280; margin-top: -0.25rem;",
              "Raw-vs-smoothed comparison is only available when the selected data includes raw traces."
            )
          },
          shiny::checkboxInput(
            "include_overall",
            "Include overall mean trace",
            value = if (!is.null(input$include_overall)) isTRUE(input$include_overall) else TRUE
          ),
          shiny::checkboxInput(
            "overlay_markers",
            "Overlay waveform markers",
            value = if (!is.null(input$overlay_markers)) isTRUE(input$overlay_markers) else FALSE
          ),
          shiny::checkboxInput(
            "include_photocell",
            "Include photocell overlay",
            value = if (!is.null(input$include_photocell)) isTRUE(input$include_photocell) else TRUE
          ),
          shiny::checkboxInput(
            "use_se",
            "Show SE on intensity response",
            value = if (!is.null(input$use_se)) isTRUE(input$use_se) else TRUE
          )
        )
      })

      output$import_summary <- shiny::renderUI({
        bundle <- imported_data()

        if (is.null(bundle) || length(bundle$items) == 0L) {
          return(
            shiny::div(
              class = "zeus-import-summary",
              shiny::div(shiny::strong("Status:"), shiny::span(import_status())),
              shiny::div(shiny::span("Import an ABF file to populate the waveform and statistics tabs."))
            )
          )
        }

        summary_rows <- lapply(bundle$items, function(item) {
          stim_n <- if (!is.null(item$data$traces_70)) dplyr::n_distinct(item$data$traces_70$stim_index) else NA_integer_
          sweep_n <- if (!is.null(item$data$traces_280)) dplyr::n_distinct(item$data$traces_280$sweep) else NA_integer_

          shiny::tags$div(
            shiny::strong(paste0(item$label, ":")),
            shiny::span(paste(item$file_name, "| stimuli", stim_n, "| sweeps", sweep_n))
          )
        })

        shiny::div(
          class = "zeus-import-summary",
          shiny::div(shiny::strong("Status:"), shiny::span(import_status())),
          shiny::div(shiny::strong("Analysis:"), shiny::span(bundle$analysis_label)),
          shiny::div(shiny::strong("Files loaded:"), shiny::span(length(bundle$items))),
          summary_rows
        )
      })

      output$selected_file_ui <- shiny::renderUI({
        bundle <- imported_data()

        if (is.null(bundle) || length(bundle$items) <= 1L) {
          return(NULL)
        }

        choices <- stats::setNames(
          vapply(bundle$items, `[[`, character(1), "id"),
          vapply(bundle$items, `[[`, character(1), "title")
        )

        shiny::div(
          class = "zeus-card",
          style = "padding: 0.85rem 1rem; margin-bottom: 1rem;",
          shiny::selectInput(
            "selected_file_id",
            "Selected file",
            choices = choices,
            selected = .zeus_app_if_null(input$selected_file_id, bundle$items[[1]]$id)
          )
        )
      })

      selected_item <- shiny::reactive({
        bundle <- imported_data()
        if (is.null(bundle) || length(bundle$items) == 0L) {
          return(NULL)
        }

        if (length(bundle$items) == 1L) {
          return(bundle$items[[1]])
        }

        selected_id <- .zeus_app_if_null(input$selected_file_id, bundle$items[[1]]$id)
        item <- .zeus_app_find_item(bundle, selected_id)

        if (is.null(item)) bundle$items[[1]] else item
      })

      output$home_tab_ui <- shiny::renderUI({
        item <- selected_item()

        if (is.null(item)) {
          return(
            bslib::card(
              class = "zeus-card",
              full_screen = TRUE,
              bslib::card_header("Mean Waveform"),
              shiny::plotOutput("home_placeholder_plot", height = "520px")
            )
          )
        }

        shiny::tagList(
          bslib::card(
            class = "zeus-card",
            full_screen = TRUE,
            bslib::card_header(paste(item$title, "Mean Waveform")),
            shiny::uiOutput("mean_wavelength_ui"),
            shiny::div(
              class = "zeus-card-actions",
              shiny::downloadButton("download_mean_plot", "Download Mean ERG", class = "btn-zeus-utility")
            ),
            shiny::plotOutput("mean_plot", height = "520px")
          ),
          bslib::card(
            class = "zeus-card",
            bslib::card_header(paste(item$title, "Average Peaks")),
            DT::DTOutput("average_peaks_table")
          )
        )
      })

      output$spectral_tab_ui <- shiny::renderUI({
        item <- selected_item()

        if (is.null(item)) {
          return(
            bslib::card(
              class = "zeus-card",
              full_screen = TRUE,
              bslib::card_header("Spectral Waveform"),
              shiny::plotOutput("spectral_placeholder_plot", height = "760px")
            )
          )
        }

        bslib::card(
          class = "zeus-card",
          full_screen = TRUE,
          bslib::card_header(paste(item$title, "Spectral Waveform")),
          shiny::div(
            class = "zeus-card-actions",
            shiny::downloadButton("download_spectral_plot", "Download Spectral Waveform", class = "btn-zeus-utility")
          ),
          shiny::plotOutput("spectral_plot", height = "760px")
        )
      })

      output$intensity_tab_ui <- shiny::renderUI({
        item <- selected_item()

        if (is.null(item)) {
          return(
            bslib::card(
              class = "zeus-card",
              full_screen = TRUE,
              bslib::card_header("Intensity Response"),
              shiny::plotOutput("intensity_placeholder_plot", height = "620px")
            )
          )
        }

        bslib::card(
          class = "zeus-card",
          full_screen = TRUE,
          bslib::card_header(paste(item$title, "Intensity Response")),
          shiny::div(
            class = "zeus-card-actions",
            shiny::downloadButton("download_intensity_plot", "Download Intensity Response", class = "btn-zeus-utility")
          ),
          shiny::plotOutput("intensity_plot", height = "620px")
        )
      })

      output$statistics_tab_ui <- shiny::renderUI({
        item <- selected_item()

        if (is.null(item)) {
          return(
            bslib::card(
              class = "zeus-card",
              bslib::card_header("Summary Tables"),
              shiny::tags$p("Import a file to review peak statistics.")
            )
          )
        }

        bslib::navset_card_pill(
          title = paste(item$title, "Summary Tables"),
          bslib::nav_panel("Key Statistics", DT::DTOutput("key_statistics_table")),
          bslib::nav_panel("By ND", DT::DTOutput("by_nd_table")),
          bslib::nav_panel("By Wavelength", DT::DTOutput("by_wavelength_table")),
          bslib::nav_panel("Combined Export", DT::DTOutput("combined_statistics_table"))
        )
      })

      output$home_placeholder_plot <- shiny::renderPlot({
        .zeus_app_placeholder_plot("Import an ABF file to view the mean waveform.")
      }, res = 144)

      output$spectral_placeholder_plot <- shiny::renderPlot({
        .zeus_app_placeholder_plot("Import an ABF file to view the spectral waveform panel.")
      }, res = 144)

      output$intensity_placeholder_plot <- shiny::renderPlot({
        .zeus_app_placeholder_plot("Import an ABF file to view the intensity-response plot.")
      }, res = 144)

      shiny::observeEvent(input$import_file, {
        shiny::req(input$abf_file)

        files <- input$abf_file
        mode <- .zeus_app_if_null(input$analysis_mode, "single")
        protocols <- .zeus_app_collect_protocols(input, files, analysis_mode = mode)
        validation_error <- .zeus_app_validate_protocols(protocols, analysis_mode = mode)

        if (!is.null(validation_error)) {
          imported_data(NULL)
          import_status(validation_error)
          shiny::showNotification(validation_error, type = "error", duration = NULL)
          return()
        }

        if (identical(mode, "single") && nrow(files) != 1L) {
          message_text <- "Single File mode requires exactly 1 ABF file."
          imported_data(NULL)
          import_status(message_text)
          shiny::showNotification(message_text, type = "error", duration = NULL)
          return()
        }

        settings <- .zeus_app_import_settings(input)
        import_status("Importing file...")
        item_labels <- .zeus_app_make_item_labels(protocols)

        result <- tryCatch(
          shiny::withProgress(message = "Importing ABF file", value = 0.05, {
            imported_items <- vector("list", length = nrow(files))

            for (i in seq_len(nrow(files))) {
              shiny::incProgress(
                amount = 0.8 / nrow(files),
                detail = paste("Reading", files$name[[i]])
              )

              x <- do.call(
                zeus_read_abf,
                c(
                  list(
                    path = files$datapath[[i]],
                    protocol = protocols[[i]]
                  ),
                  settings
                )
              )

              x <- .zeus_app_attach_source_info(x, file_name = files$name[[i]])
              imported_items[[i]] <- list(
                id = paste0("file", i),
                label = item_labels[[i]],
                title = paste(item_labels[[i]], files$name[[i]], sep = " | "),
                file_name = files$name[[i]],
                protocol = protocols[[i]],
                data = x
              )
            }

            shiny::incProgress(0.15, detail = "Preparing plots and summaries")
            list(
              analysis_mode = mode,
              analysis_label = if (identical(mode, "full_fish")) "Full Fish Analysis" else "Single File",
              items = imported_items
            )
          }),
          error = function(e) e
        )

        if (inherits(result, "error")) {
          imported_data(NULL)
          import_status(paste("Import failed:", conditionMessage(result)))
          shiny::showNotification(
            paste("Import failed:", conditionMessage(result)),
            type = "error",
            duration = NULL
          )
          return()
        }

        imported_data(result)
        import_status("Import complete.")
        shiny::updateCheckboxInput(session, "compare_raw", value = FALSE)

        shiny::showNotification("ABF file imported successfully.", type = "message")
      }, ignoreInit = TRUE)

      shiny::observeEvent(input$close_app, {
        shiny::stopApp(invisible(NULL))
      }, ignoreInit = TRUE)

      current_stats <- shiny::reactive({
        item <- selected_item()
        shiny::req(item)

        do.call(
          zeus_summarize_peak_statistics,
          c(list(x = item$data), .zeus_app_peak_settings(input))
        )
      })

      output$mean_wavelength_ui <- shiny::renderUI({
        item <- selected_item()
        if (is.null(item) || !identical(item$protocol, "C0")) {
          return(NULL)
        }

        choices <- .zeus_app_available_wavelengths(item$data)
        if (length(choices) == 0L) {
          choices <- "450"
        }
        selected_choice <- if ("450" %in% choices) "450" else choices[[1]]

        shiny::selectInput(
          "mean_wavelength",
          "C0 waveform block",
          choices = choices,
          selected = .zeus_app_if_null(input$mean_wavelength, selected_choice)
        )
      })

      output$mean_plot <- shiny::renderPlot({
        item <- selected_item()
        if (is.null(item)) {
          return(.zeus_app_placeholder_plot("Import an ABF file to view the mean waveform."))
        }

        slots <- .zeus_app_available_slots(item$data)
        selected_slot <- if (!is.null(input$mean_data_slot) && input$mean_data_slot %in% slots) {
          input$mean_data_slot
        } else {
          slots[[1]]
        }
        wavelength_select <- if (identical(item$protocol, "C0")) input$mean_wavelength else NULL
        compare_raw <- isTRUE(input$compare_raw) &&
          .zeus_app_can_compare_raw(item$data, selected_slot)

        zeus_plot_mean_waveform(
          x = item$data,
          data_slot = selected_slot,
          compare_raw = compare_raw,
          include_overall = isTRUE(input$include_overall),
          overlay_markers = isTRUE(input$overlay_markers),
          include_photocell = isTRUE(input$include_photocell),
          wavelength_select = wavelength_select,
          a_window = c(input$awave_window_start, input$awave_window_end),
          b_window = c(input$peak_response_start, input$peak_response_end),
          d_window = c(input$dwave_window_start, input$dwave_window_end)
        )
      }, res = 144)

      output$spectral_plot <- shiny::renderPlot({
        item <- selected_item()
        if (is.null(item)) {
          return(.zeus_app_placeholder_plot("Import an ABF file to view the spectral waveform panel."))
        }

        zeus_plot_spectral_waveform(
          x = item$data,
          include_photocell = isTRUE(input$include_photocell)
        )
      }, res = 144)

      output$intensity_plot <- shiny::renderPlot({
        item <- selected_item()
        if (is.null(item)) {
          return(.zeus_app_placeholder_plot("Import an ABF file to view the intensity-response plot."))
        }

        zeus_plot_intensity_response(
          x = item$data,
          use_se = isTRUE(input$use_se)
        )
      }, res = 144)

      output$average_peaks_table <- DT::renderDT({
        .zeus_app_table(.zeus_app_peak_summary_table(current_stats()))
      })

      output$key_statistics_table <- DT::renderDT({
        .zeus_app_table(current_stats()$key_statistics)
      })

      output$by_nd_table <- DT::renderDT({
        .zeus_app_table(current_stats()$by_nd)
      })

      output$by_wavelength_table <- DT::renderDT({
        .zeus_app_table(current_stats()$by_wavelength)
      })

      output$combined_statistics_table <- DT::renderDT({
        .zeus_app_table(current_stats()$combined_export)
      })

      output$download_csv_bundle <- shiny::downloadHandler(
        filename = function() {
          file_stem <- trimws(input$export_name)
          if (!nzchar(file_stem)) {
            file_stem <- "zeus_export"
          }
          paste0(file_stem, "_bundle.zip")
        },
        content = function(file) {
          bundle <- imported_data()
          shiny::req(bundle)
          file_stem <- trimws(input$export_name)
          if (!nzchar(file_stem)) {
            file_stem <- "zeus_export"
          }

          shiny::withProgress(message = "Preparing CSV export", value = 0, {
            bundle_dir <- file.path(tempdir(), paste0(file_stem, "_bundle"))
            unlink(bundle_dir, recursive = TRUE, force = TRUE)
            shiny::incProgress(0.1, detail = "Creating plot and CSV files")

            for (item in bundle$items) {
              slots <- .zeus_app_available_slots(item$data)
              selected_slot <- if (!is.null(input$mean_data_slot) && input$mean_data_slot %in% slots) {
                input$mean_data_slot
              } else {
                slots[[1]]
              }
              wavelength_select <- if (identical(item$protocol, "C0")) {
                input[[paste0("mean_wavelength_", item$id)]]
              } else {
                NULL
              }
              compare_raw <- isTRUE(input$compare_raw) &&
                .zeus_app_can_compare_raw(item$data, selected_slot)

              item_dir <- if (length(bundle$items) > 1L) {
                file.path(bundle_dir, .zeus_app_sanitize_stem(item$label))
              } else {
                bundle_dir
              }

              .zeus_app_write_export_bundle(
                bundle_dir = item_dir,
                file_stem = .zeus_app_sanitize_stem(tools::file_path_sans_ext(item$file_name)),
                mean_plot = zeus_plot_mean_waveform(
                  x = item$data,
                  data_slot = selected_slot,
                  compare_raw = compare_raw,
                  include_overall = isTRUE(input$include_overall),
                  overlay_markers = isTRUE(input$overlay_markers),
                  include_photocell = isTRUE(input$include_photocell),
                  wavelength_select = wavelength_select
                ),
                spectral_plot = zeus_plot_spectral_waveform(
                  x = item$data,
                  include_photocell = isTRUE(input$include_photocell)
                ),
                intensity_plot = zeus_plot_intensity_response(
                  x = item$data,
                  use_se = isTRUE(input$use_se)
                ),
                x = item$data
              )
            }

            shiny::incProgress(0.65, detail = "Compressing CSV bundle")
            zip_files <- list.files(bundle_dir, recursive = TRUE)
            tmp_zip <- tempfile(pattern = paste0(file_stem, "_"), fileext = ".zip")
            .zeus_zip_files(
              zipfile = tmp_zip,
              files = zip_files,
              root_dir = bundle_dir
            )

            shiny::incProgress(0.2, detail = "Finalizing download")
            ok <- file.copy(tmp_zip, file, overwrite = TRUE)
            if (!isTRUE(ok) || !file.exists(file)) {
              stop("Failed to prepare ZIP download.", call. = FALSE)
            }
            shiny::incProgress(0.05)
          })
        },
        contentType = "application/zip"
      )

      output$download_mean_plot <- shiny::downloadHandler(
        filename = function() {
          item <- selected_item()
          shiny::req(item)
          paste0(.zeus_app_export_stem(input, item), "_mean_waveform.png")
        },
        content = function(file) {
          item <- selected_item()
          shiny::req(item)

          slots <- .zeus_app_available_slots(item$data)
          selected_slot <- if (!is.null(input$mean_data_slot) && input$mean_data_slot %in% slots) {
            input$mean_data_slot
          } else {
            slots[[1]]
          }
          wavelength_select <- if (identical(item$protocol, "C0")) input$mean_wavelength else NULL
          compare_raw <- isTRUE(input$compare_raw) &&
            .zeus_app_can_compare_raw(item$data, selected_slot)

          plot_obj <- zeus_plot_mean_waveform(
            x = item$data,
            data_slot = selected_slot,
            compare_raw = compare_raw,
            include_overall = isTRUE(input$include_overall),
            overlay_markers = isTRUE(input$overlay_markers),
            include_photocell = isTRUE(input$include_photocell),
            wavelength_select = wavelength_select,
            a_window = c(input$awave_window_start, input$awave_window_end),
            b_window = c(input$peak_response_start, input$peak_response_end),
            d_window = c(input$dwave_window_start, input$dwave_window_end)
          )

          ggplot2::ggsave(file, plot = plot_obj, width = 11, height = 6.5, dpi = 320, bg = "white")
        },
        contentType = "image/png"
      )

      output$download_spectral_plot <- shiny::downloadHandler(
        filename = function() {
          item <- selected_item()
          shiny::req(item)
          paste0(.zeus_app_export_stem(input, item), "_spectral_waveform.png")
        },
        content = function(file) {
          item <- selected_item()
          shiny::req(item)

          plot_obj <- zeus_plot_spectral_waveform(
            x = item$data,
            include_photocell = isTRUE(input$include_photocell)
          )

          ggplot2::ggsave(file, plot = plot_obj, width = 14, height = 9, dpi = 320, bg = "white")
        },
        contentType = "image/png"
      )

      output$download_intensity_plot <- shiny::downloadHandler(
        filename = function() {
          item <- selected_item()
          shiny::req(item)
          paste0(.zeus_app_export_stem(input, item), "_intensity_response.png")
        },
        content = function(file) {
          item <- selected_item()
          shiny::req(item)

          plot_obj <- zeus_plot_intensity_response(
            x = item$data,
            use_se = isTRUE(input$use_se)
          )

          ggplot2::ggsave(file, plot = plot_obj, width = 11, height = 6.5, dpi = 320, bg = "white")
        },
        contentType = "image/png"
      )

      shiny::outputOptions(output, "download_csv_bundle", suspendWhenHidden = FALSE)
      shiny::outputOptions(output, "download_mean_plot", suspendWhenHidden = FALSE)
      shiny::outputOptions(output, "download_spectral_plot", suspendWhenHidden = FALSE)
      shiny::outputOptions(output, "download_intensity_plot", suspendWhenHidden = FALSE)
    }
  )
}

#' Run the ZEUS Shiny application
#'
#' Launches the ZEUS Shiny app for interactive file import, visualization,
#' statistics review, and export.
#'
#' @param launch.browser Logical passed to [shiny::runApp()]. Default is `TRUE`.
#' @param display.mode Character passed to [shiny::runApp()]. Default is
#'   `"normal"`.
#' @param max_upload_mb Maximum upload size allowed by Shiny, in megabytes.
#'   Default is `100`.
#' @param ... Additional arguments passed to [shiny::runApp()].
#'
#' @return Invisibly returns the running Shiny app object.
#' @export
run_zeus_app <- function(launch.browser = TRUE,
                         display.mode = "normal",
                         max_upload_mb = 100,
                         ...) {
  old_max <- getOption("shiny.maxRequestSize")
  options(shiny.maxRequestSize = as.numeric(max_upload_mb) * 1024^2)
  on.exit(options(shiny.maxRequestSize = old_max), add = TRUE)

  shiny::runApp(
    appObj = zeus_app(),
    launch.browser = launch.browser,
    display.mode = display.mode,
    ...
  )
}
