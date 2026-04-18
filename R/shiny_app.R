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

.zeus_app_write_export_bundle <- function(bundle_dir,
                                          file_stem,
                                          export_format,
                                          mean_plot,
                                          spectral_plot,
                                          intensity_plot,
                                          x) {
  dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)

  export_paths <- .zeus_app_export_paths(bundle_dir, file_stem)

  ggplot2::ggsave(export_paths$mean_plot, plot = mean_plot, width = 11, height = 6.5, dpi = 320, bg = "white")
  ggplot2::ggsave(export_paths$spectral_plot, plot = spectral_plot, width = 14, height = 9, dpi = 320, bg = "white")
  ggplot2::ggsave(export_paths$intensity_plot, plot = intensity_plot, width = 11, height = 6.5, dpi = 320, bg = "white")

  if (identical(export_format, "excel")) {
    zeus_export_excel_workbook(x, file.path(bundle_dir, paste0(file_stem, "_data.xlsx")))
  } else {
    zeus_export_csv_bundle(x, export_paths$csv_base)
  }

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
    bg = "#F6F1E8",
    fg = "#14253D",
    primary = "#C42E4E",
    secondary = "#253B6F",
    success = "#6A8C69",
    info = "#385D8A",
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
      shiny::fileInput(
        "abf_file",
        "ABF File",
        accept = c(".abf")
      ),
      shiny::selectInput(
        "protocol",
        "Protocol",
        choices = c("C0", "C1"),
        selected = "C1"
      ),
      shiny::actionButton("import_file", "Import File", class = "btn-zeus-primary")
    ),
    bslib::accordion(
      id = "zeus_controls",
      open = c("Import Settings", "Plot Settings", "Peak Settings", "Export"),
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
        shiny::selectInput(
          "mean_data_slot",
          "Mean waveform source",
          choices = c("traces_70", "traces_280"),
          selected = "traces_70"
        ),
        shiny::uiOutput("mean_wavelength_ui"),
        shiny::uiOutput("compare_raw_ui"),
        shiny::checkboxInput("include_overall", "Include overall mean trace", value = TRUE),
        shiny::checkboxInput("overlay_markers", "Overlay waveform markers", value = FALSE),
        shiny::checkboxInput("include_photocell", "Include photocell overlay", value = TRUE),
        shiny::checkboxInput("use_se", "Show SE on intensity response", value = TRUE)
      ),
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
        "Export",
        shiny::radioButtons(
          "export_format",
          "Data export format",
          choices = c("CSV bundle" = "csv", "Excel workbook" = "excel"),
          selected = "csv",
          inline = TRUE
        ),
        shiny::textInput("export_name", "Export file name", value = ""),
        shiny::downloadButton("download_bundle", "Download Data and Plots", class = "btn-zeus-secondary"),
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
            shiny::tabsetPanel(
              id = "zeus_main_tabs",
              type = "tabs",
              shiny::tabPanel(
                "Home",
                bslib::card(
                  class = "zeus-card",
                  full_screen = TRUE,
                  bslib::card_header("Mean Waveform"),
                  shiny::plotOutput("mean_plot", height = "520px")
                ),
                bslib::card(
                  class = "zeus-card",
                  bslib::card_header("Average Peaks"),
                  DT::DTOutput("average_peaks_table")
                )
              ),
              shiny::tabPanel(
                "Spectral Waveform",
                bslib::card(
                  class = "zeus-card",
                  full_screen = TRUE,
                  bslib::card_header("Spectral Waveform"),
                  shiny::plotOutput("spectral_plot", height = "760px")
                )
              ),
              shiny::tabPanel(
                "Intensity Response",
                bslib::card(
                  class = "zeus-card",
                  full_screen = TRUE,
                  bslib::card_header("Intensity Response"),
                  shiny::plotOutput("intensity_plot", height = "620px")
                )
              ),
              shiny::tabPanel(
                "Peak Statistics",
                bslib::navset_card_pill(
                  title = "Summary Tables",
                  bslib::nav_panel("Key Statistics", DT::DTOutput("key_statistics_table")),
                  bslib::nav_panel("By ND", DT::DTOutput("by_nd_table")),
                  bslib::nav_panel("By Wavelength", DT::DTOutput("by_wavelength_table")),
                  bslib::nav_panel("Combined Export", DT::DTOutput("combined_statistics_table"))
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

      output$mean_wavelength_ui <- shiny::renderUI({
        if (is.null(imported_data()) && !identical(input$protocol, "C0")) {
          return(NULL)
        }

        if (!identical(input$protocol, "C0")) {
          return(NULL)
        }

        choices <- .zeus_app_available_wavelengths(imported_data())
        if (length(choices) == 0L) {
          choices <- "450"
        }

        shiny::selectInput(
          "mean_wavelength",
          "C0 waveform block",
          choices = choices,
          selected = choices[[1]]
        )
      })

      output$compare_raw_ui <- shiny::renderUI({
        data_slot <- if (!is.null(input$mean_data_slot)) input$mean_data_slot else "traces_70"
        can_compare <- .zeus_app_can_compare_raw(imported_data(), data_slot)

        shiny::tagList(
          shiny::checkboxInput(
            "compare_raw",
            "Compare raw and smoothed mean traces",
            value = FALSE
          ),
          if (!can_compare) {
            shiny::tags$div(
              style = "font-size: 0.9rem; color: #6b7280; margin-top: -0.25rem;",
              "Raw-vs-smoothed comparison is only available when the selected data includes raw traces."
            )
          }
        )
      })

      output$import_summary <- shiny::renderUI({
        x <- imported_data()

        if (is.null(x)) {
          return(
            shiny::div(
              class = "zeus-import-summary",
              shiny::div(shiny::strong("Status:"), shiny::span(import_status())),
              shiny::div(shiny::span("Import an ABF file to populate the waveform and statistics tabs."))
            )
          )
        }

        file_label <- if (!is.null(x$source_file_name) && nzchar(x$source_file_name)) {
          x$source_file_name
        } else {
          basename(if (!is.null(x$path)) x$path else "")
        }
        protocol_label <- if (!is.null(x$protocol)) x$protocol else input$protocol
        stim_n <- if (!is.null(x$traces_70)) dplyr::n_distinct(x$traces_70$stim_index) else NA_integer_
        sweep_n <- if (!is.null(x$traces_280)) dplyr::n_distinct(x$traces_280$sweep) else NA_integer_

        shiny::div(
          class = "zeus-import-summary",
          shiny::div(shiny::strong("Status:"), shiny::span(import_status())),
          shiny::div(shiny::strong("Loaded file:"), shiny::span(file_label)),
          shiny::div(shiny::strong("Protocol:"), shiny::span(protocol_label)),
          shiny::div(shiny::strong("Stimuli:"), shiny::span(stim_n)),
          shiny::div(shiny::strong("Sweeps:"), shiny::span(sweep_n))
        )
      })

      shiny::observeEvent(input$import_file, {
        shiny::req(input$abf_file)

        settings <- .zeus_app_import_settings(input)
        import_status("Importing file...")

        result <- tryCatch(
          shiny::withProgress(message = "Importing ABF file", value = 0.15, {
            shiny::incProgress(0.2, detail = "Reading ABF data")

            x <- do.call(
              zeus_read_abf,
              c(
                list(
                  path = input$abf_file$datapath,
                  protocol = input$protocol
                ),
                settings
              )
            )

            shiny::incProgress(0.8, detail = "Preparing plots and summaries")
            .zeus_app_attach_source_info(x, file_name = input$abf_file$name)
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

        choices <- .zeus_app_available_wavelengths(imported_data())
        if (length(choices) > 0L) {
          shiny::updateSelectInput(
            session,
            "mean_wavelength",
            choices = choices,
            selected = choices[[1]]
          )
        }

        shiny::showNotification("ABF file imported successfully.", type = "message")
      }, ignoreInit = TRUE)

      shiny::observeEvent(input$close_app, {
        shiny::stopApp(invisible(NULL))
      }, ignoreInit = TRUE)

      peak_statistics <- shiny::reactive({
        x <- imported_data()
        shiny::req(x)

        do.call(
          zeus_summarize_peak_statistics,
          c(list(x = x), .zeus_app_peak_settings(input))
        )
      })

      mean_plot_obj <- shiny::reactive({
        x <- imported_data()
        if (is.null(x)) {
          return(.zeus_app_placeholder_plot("Import an ABF file to view the mean waveform."))
        }

        wavelength_select <- if (identical(x$protocol, "C0")) input$mean_wavelength else NULL
        compare_raw <- isTRUE(input$compare_raw) &&
          .zeus_app_can_compare_raw(x, input$mean_data_slot)

        zeus_plot_mean_waveform(
          x = x,
          data_slot = input$mean_data_slot,
          compare_raw = compare_raw,
          include_overall = isTRUE(input$include_overall),
          overlay_markers = isTRUE(input$overlay_markers),
          include_photocell = isTRUE(input$include_photocell),
          wavelength_select = wavelength_select
        )
      })

      spectral_plot_obj <- shiny::reactive({
        x <- imported_data()
        if (is.null(x)) {
          return(.zeus_app_placeholder_plot("Import an ABF file to view the spectral waveform panel."))
        }

        zeus_plot_spectral_waveform(
          x = x,
          include_photocell = isTRUE(input$include_photocell)
        )
      })

      intensity_plot_obj <- shiny::reactive({
        x <- imported_data()
        if (is.null(x)) {
          return(.zeus_app_placeholder_plot("Import an ABF file to view the intensity-response plot."))
        }

        zeus_plot_intensity_response(
          x = x,
          use_se = isTRUE(input$use_se)
        )
      })

      output$mean_plot <- shiny::renderPlot({
        mean_plot_obj()
      }, res = 144)

      output$spectral_plot <- shiny::renderPlot({
        spectral_plot_obj()
      }, res = 144)

      output$intensity_plot <- shiny::renderPlot({
        intensity_plot_obj()
      }, res = 144)

      output$average_peaks_table <- DT::renderDT({
        .zeus_app_table(.zeus_app_peak_summary_table(peak_statistics()))
      })

      output$key_statistics_table <- DT::renderDT({
        .zeus_app_table(peak_statistics()$key_statistics)
      })

      output$by_nd_table <- DT::renderDT({
        .zeus_app_table(peak_statistics()$by_nd)
      })

      output$by_wavelength_table <- DT::renderDT({
        .zeus_app_table(peak_statistics()$by_wavelength)
      })

      output$combined_statistics_table <- DT::renderDT({
        .zeus_app_table(peak_statistics()$combined_export)
      })

      output$download_bundle <- shiny::downloadHandler(
        filename = function() {
          file_stem <- trimws(input$export_name)
          if (!nzchar(file_stem)) {
            file_stem <- "zeus_export"
          }
          paste0(file_stem, "_bundle.zip")
        },
        content = function(file) {
        x <- imported_data()
        shiny::req(x)
          file_stem <- trimws(input$export_name)
          if (!nzchar(file_stem)) {
            file_stem <- "zeus_export"
          }

          bundle_dir <- file.path(tempdir(), paste0(file_stem, "_bundle"))
          unlink(bundle_dir, recursive = TRUE, force = TRUE)

          shiny::withProgress(message = "Preparing export zip", value = 0.15, {
            .zeus_app_write_export_bundle(
              bundle_dir = bundle_dir,
              file_stem = file_stem,
              export_format = input$export_format,
              mean_plot = mean_plot_obj(),
              spectral_plot = spectral_plot_obj(),
              intensity_plot = intensity_plot_obj(),
              x = x
            )
            shiny::incProgress(0.7)

            zip_files <- list.files(bundle_dir)
            .zeus_zip_files(
              zipfile = file,
              files = zip_files,
              root_dir = bundle_dir
            )
            shiny::incProgress(0.15)
          })
        }
      )
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
