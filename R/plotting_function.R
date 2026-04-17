# Helper, Normalize Zeus Plotting Data ------------------------------------
#' Normalize ZEUS plotting input to a standard waveform data frame
#'
#' @description
#' Internal helper used by ZEUS plotting functions to standardize input objects
#' before plotting.
#'
#' This function accepts:
#' \itemize{
#'   \item a `zeus_stimresp` object,
#'   \item a waveform data frame, or
#'   \item a list-like ZEUS object containing an embedded waveform data frame.
#' }
#'
#' The helper resolves the appropriate waveform table, standardizes common
#' alternate column names used in newer ZEUS import pipelines, reconstructs
#' missing time columns when possible, and optionally filters by channel.
#'
#' Standardized output column names include:
#' \itemize{
#'   \item `time`
#'   \item `time_ms`
#'   \item `value`
#'   \item `stim_nd`
#'   \item `stim_label`
#'   \item `stim_index`
#'   \item `wavelength`
#' }
#'
#' @param x Input object. May be a `zeus_stimresp` object, a data frame,
#'   or a list-like ZEUS object containing waveform data.
#' @param data_slot Character string indicating which slot to use when `x`
#'   is a `zeus_stimresp` object. One of `"traces_70"` or `"traces_280"`.
#' @param channel_filter Optional character string specifying a channel to
#'   retain. If `NULL`, no channel filtering is applied.
#' @param require_cols Character vector of required column names after
#'   standardization.
#' @param allow_stimresp Logical; if `TRUE`, `zeus_stimresp` objects are
#'   handled explicitly.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{df_plot}{A standardized waveform data frame ready for plotting.}
#'   \item{df_photocell_source}{A data frame that can be used to derive a
#'   photocell trace when available.}
#' }
#'
#' @keywords internal
.zeus_prepare_plot_df <- function(
    x,
    data_slot = c("traces_70", "traces_280"),
    channel_filter = NULL,
    require_cols = c("time", "value"),
    allow_stimresp = TRUE
) {
  data_slot <- match.arg(data_slot)

  df_plot <- NULL
  df_photocell_source <- NULL

  # Case 1: zeus_stimresp object
  if (isTRUE(allow_stimresp) && inherits(x, "zeus_stimresp")) {
    if (!data_slot %in% names(x)) {
      stop("`data_slot` not found in `x`.", call. = FALSE)
    }

    df_plot <- x[[data_slot]]

    if (!is.null(x$photocell)) {
      df_photocell_source <- x$photocell
    } else if ("traces_280" %in% names(x) && is.data.frame(x$traces_280)) {
      df_photocell_source <- x$traces_280
    } else {
      df_photocell_source <- df_plot
    }
  }

  # Case 2: already a data frame
  else if (is.data.frame(x)) {
    df_plot <- x
    df_photocell_source <- x
  }

  # Case 3: generic list-like ZEUS object with an embedded waveform table
  else if (is.list(x)) {
    candidate_slots <- c(
      "df_long",
      "data",
      "plot_data",
      "traces",
      "traces_280",
      "traces_70",
      "waveforms"
    )

    for (nm in candidate_slots) {
      if (nm %in% names(x) && is.data.frame(x[[nm]])) {
        df_plot <- x[[nm]]
        df_photocell_source <- x[[nm]]
        break
      }
    }

    # Fallback: first data frame found in the list
    if (is.null(df_plot)) {
      df_candidates <- x[vapply(x, is.data.frame, logical(1))]
      if (length(df_candidates) > 0) {
        df_plot <- df_candidates[[1]]
        df_photocell_source <- df_plot
      }
    }
  }

  if (is.null(df_plot) || !is.data.frame(df_plot)) {
    stop(
      paste0(
        "`x` must be a data frame, a `zeus_stimresp` object, ",
        "or a ZEUS list object containing a waveform data frame."
      ),
      call. = FALSE
    )
  }

  # Standardize alternate names used by ZEUS import pipelines

  # Time
  if (!("time" %in% names(df_plot))) {
    if ("time_ms" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::mutate(time = .data$time_ms / 1000)
    } else if ("time_s" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::rename(time = .data$time_s)
    } else if ("t" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::rename(time = .data$t)
    }
  }

  # Signal/value
  if (!("value" %in% names(df_plot))) {
    if ("signal" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::rename(value = .data$signal)
    } else if ("mean_value" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::rename(value = .data$mean_value)
    } else if ("response" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::rename(value = .data$response)
    }
  }

  # Raw signal/value
  if (!("value_raw" %in% names(df_plot))) {
    if ("signal_raw" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::rename(value_raw = .data$signal_raw)
    } else if ("raw_value" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::rename(value_raw = .data$raw_value)
    } else if ("response_raw" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::rename(value_raw = .data$response_raw)
    }
  }

  # Sweep
  if (!("sweep" %in% names(df_plot))) {
    if ("sweep_id" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::rename(sweep = .data$sweep_id)
    }
  }

  # Stimulus ND
  if (!("stim_nd" %in% names(df_plot))) {
    if ("expected_stim_nd" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::rename(stim_nd = .data$expected_stim_nd)
    }
  }

  # Wavelength
  if (!("wavelength" %in% names(df_plot))) {
    if ("stim_wavelength" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::rename(wavelength = .data$stim_wavelength)
    }
  }

  # Stimulus label
  if (!("stim_label" %in% names(df_plot))) {
    if ("expected_stim_label" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::rename(stim_label = .data$expected_stim_label)
    } else if ("stimulus_label" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::rename(stim_label = .data$stimulus_label)
    }
  }

  # Stimulus index
  if (!("stim_index" %in% names(df_plot))) {
    if ("stimulus_index" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::rename(stim_index = .data$stimulus_index)
    }
  }

  # Reconstruct time if still missing
  if (!("time" %in% names(df_plot))) {
    if (all(c("sample_index", "dt") %in% names(df_plot))) {
      df_plot <- df_plot |>
        dplyr::mutate(time = (.data$sample_index - 1) * .data$dt)
    } else if (all(c("sample", "dt") %in% names(df_plot))) {
      df_plot <- df_plot |>
        dplyr::mutate(time = (.data$sample - 1) * .data$dt)
    }
  }

  # Final time_ms creation
  if (!("time_ms" %in% names(df_plot)) && "time" %in% names(df_plot)) {
    df_plot <- df_plot |>
      dplyr::mutate(time_ms = zeus_time_to_ms(.data$time))
  }

  # Standard fallback fields
  if (!("stim_nd" %in% names(df_plot))) {
    df_plot <- df_plot |>
      dplyr::mutate(stim_nd = NA_real_)
  }

  if (!("wavelength" %in% names(df_plot))) {
    df_plot <- df_plot |>
      dplyr::mutate(wavelength = NA_character_)
  }

  if (!("stim_label" %in% names(df_plot))) {
    if ("stim_nd" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::mutate(stim_label = as.character(.data$stim_nd))
    }
  }

  if (!("stim_index" %in% names(df_plot))) {
    if ("stim_label" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::group_by(.data$stim_label) |>
        dplyr::mutate(stim_index = dplyr::cur_group_id()) |>
        dplyr::ungroup()
    } else {
      df_plot <- df_plot |>
        dplyr::mutate(stim_index = 1L)
    }
  }

  # Required column check
  missing_cols <- setdiff(require_cols, names(df_plot))
  if (length(missing_cols) > 0) {
    stop(
      "Plot data is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Optional channel filtering
  if (!is.null(channel_filter) && "channel" %in% names(df_plot)) {
    df_plot <- df_plot |>
      dplyr::filter(.data$channel == channel_filter)

    if (nrow(df_plot) == 0) {
      stop("No data left after applying `channel_filter`.", call. = FALSE)
    }
  }

  list(
    df_plot = df_plot,
    df_photocell_source = df_photocell_source
  )
}


# Photocell scaling helper ------------------------------------------------

#' Scale and position a photocell trace for Origin-like display
#'
#' Internal helper used by ZEUS plotting functions to normalize a raw photocell
#' trace so that it appears below the ERG waveform — matching the display
#' convention used in Origin StimResp outputs.
#'
#' The photocell is normalized to \code{[0, 1]}, scaled to occupy
#' `photocell_relative_height * erg_range` in Y units, and shifted so its
#' baseline sits just below `erg_ymin`.
#'
#' @param df_photocell_source A data frame containing photocell data, expected
#'   to have `time_ms` and `signal` columns.
#' @param photocell_filter Pattern passed to [grepl()] to identify the
#'   photocell channel when a `channel` column is present.
#' @param erg_ymin Numeric. Minimum Y value of the ERG signal (used for
#'   positioning).
#' @param erg_range Numeric. Total Y range of the ERG signal (used for
#'   scaling).
#' @param photocell_relative_height Fraction of `erg_range` that the photocell
#'   pulse should occupy. Default is `1.1`.
#'
#' @return A data frame with `time_ms` and `signal` columns ready to overlay
#'   on a plot, or `NULL` if no valid photocell data could be prepared.
#' @keywords internal
.zeus_scale_photocell <- function(
    df_photocell_source,
    photocell_filter = "Photocell",
    erg_ymin,
    erg_range,
    photocell_relative_height = 1.1
) {
  if (is.null(df_photocell_source) || !is.data.frame(df_photocell_source)) {
    return(NULL)
  }

  pc <- df_photocell_source

  # Ensure time column
  if (!("time" %in% names(pc))) {
    if ("time_ms" %in% names(pc)) {
      pc <- pc |> dplyr::mutate(time = .data$time_ms / 1000)
    } else {
      return(NULL)
    }
  }

  # Ensure value column
  if (!("value" %in% names(pc))) {
    if ("signal" %in% names(pc)) {
      pc <- pc |> dplyr::rename(value = .data$signal)
    } else if ("mean_value" %in% names(pc)) {
      pc <- pc |> dplyr::rename(value = .data$mean_value)
    } else {
      return(NULL)
    }
  }

  # Ensure time_ms column
  if (!("time_ms" %in% names(pc))) {
    pc <- pc |> dplyr::mutate(time_ms = zeus_time_to_ms(.data$time))
  }

  # Filter to photocell channel when channel column is present
  if ("channel" %in% names(pc)) {
    pc <- pc |> dplyr::filter(grepl(photocell_filter, .data$channel, fixed = TRUE))
  }

  if (nrow(pc) == 0) return(NULL)

  # Average all sweeps to a single representative photocell trace
  pc <- pc |>
    dplyr::group_by(.data$time_ms) |>
    dplyr::summarise(signal = mean(.data$value, na.rm = TRUE), .groups = "drop")

  if (nrow(pc) == 0) return(NULL)

  pc_min <- min(pc$signal, na.rm = TRUE)
  pc_max <- max(pc$signal, na.rm = TRUE)
  pc_range <- pc_max - pc_min

  if (!is.finite(erg_range) || erg_range <= 0 ||
      !is.finite(pc_range) || pc_range <= 0) {
    return(NULL)
  }

  # Target: photocell occupies photocell_relative_height * erg_range in Y.
  # Baseline is positioned slightly below erg_ymin (0.1 % gap).
  target_height <- photocell_relative_height * erg_range
  pc_baseline   <- erg_ymin - 0.001 * erg_range

  pc |>
    dplyr::mutate(
      signal = ((signal - pc_min) / pc_range) * target_height + pc_baseline
    )
}

.zeus_waveform_palette <- function(n) {
  n <- max(1L, as.integer(n))

  grDevices::hcl.colors(n, palette = "Dark 3")
}

.zeus_waveform_theme <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        size = ggplot2::rel(1.28),
        colour = "#0E3254",
        margin = ggplot2::margin(b = 12)
      ),
      axis.title = ggplot2::element_text(
        face = "bold",
        colour = "#1F2A33"
      ),
      axis.text = ggplot2::element_text(
        colour = "#33424D"
      ),
      legend.title = ggplot2::element_text(
        face = "bold",
        colour = "#1F2A33"
      ),
      legend.text = ggplot2::element_text(
        size = ggplot2::rel(0.88),
        colour = "#33424D"
      ),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0.18, "cm"),
      legend.key.width = grid::unit(1.2, "cm"),
      legend.key.height = grid::unit(0.45, "cm"),
      legend.background = ggplot2::element_rect(fill = "white", color = "#D5DEE7", linewidth = 0.3),
      legend.margin = ggplot2::margin(6, 6, 6, 6),
      panel.border = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(linewidth = 0.3, colour = "#DDE5EC"),
      panel.grid.minor.y = ggplot2::element_line(linewidth = 0.2, colour = "#EEF3F7"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 0.5, color = "#2A3741"),
      axis.ticks = ggplot2::element_line(linewidth = 0.4, color = "#2A3741"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}


# Mean ERG waveform plot across sweeps ------------------------------------

#' Plot mean ERG waveforms from ZEUS StimResp output
#'
#' @description
#' Visualizes averaged ERG waveforms from ZEUS StimResp output.
#'
#' This function works with either:
#' \itemize{
#'   \item a `zeus_stimresp` object returned by `build_stimresp()` or
#'     `zeus_read_abf()`, or
#'   \item a long-format data frame containing waveform data.
#' }
#'
#' Stimulus-specific mean waveforms are drawn as thinner color-coded lines,
#' while the overall mean waveform is drawn as a thicker line.
#' Time is displayed in milliseconds.
#'
#' If marker overlay is requested, the function computes reference lines from
#' the overall mean waveform itself so the lines align with the visible mean
#' trace. The reference lines correspond to:
#' \itemize{
#'   \item A-wave trough
#'   \item B-wave peak
#'   \item D-wave peak
#' }
#'
#' @param x A `zeus_stimresp` object or a long-format waveform data frame.
#' @param data_slot For `zeus_stimresp` objects, which component to plot.
#'   Default is `"traces_70"`. Can also be `"traces_280"`.
#' @param channel_filter Optional channel to filter. Ignored if the selected
#'   data slot does not contain a `channel` column.
#' @param compare_raw Logical; if `TRUE`, plots both raw and smoothed means.
#'   Requires `value_raw` in the plotted data. Default is `FALSE`.
#' @param include_overall Logical; if `TRUE`, includes an overall mean waveform
#'   across all plotted traces.
#' @param overlay_markers Logical; if `TRUE`, overlays vertical reference lines
#'   at the A-wave trough, B-wave peak, and D-wave peak locations derived from
#'   the overall mean waveform.
#' @param include_photocell Logical; if `TRUE`, adds the mean photocell trace
#'   when available.
#' @param photocell_filter Character string used to identify the photocell
#'   channel when plotting from a data frame or from `traces_280`.
#' @param photocell_color Color for the photocell trace. Default `"black"`.
#' @param photocell_auto_scale Logical; if `TRUE` (default), rescales and
#'   positions the photocell trace below the ERG signal following the Origin
#'   StimResp convention: the pulse height equals
#'   `photocell_relative_height * ERG_range` and the baseline sits slightly
#'   below the minimum ERG value.
#' @param photocell_relative_height Fraction of the ERG amplitude range that
#'   the photocell pulse should occupy when `photocell_auto_scale = TRUE`.
#'   Default is `1.1` (110% of the ERG range).
#' @param a_window Numeric length-2 vector giving the A-wave search interval in
#'   milliseconds.
#' @param b_window Numeric length-2 vector giving the B-wave search interval in
#'   milliseconds.
#' @param d_window Numeric length-2 vector giving the D-wave search interval in
#'   milliseconds.
#' @param stim_levels Optional vector giving the desired order of stimulus labels
#'   in the legend. If `NULL`, uses the order present in the data.
#' @param color_by One of `"stim_nd"`, `"stim_label"`, or `"wavelength"`.
#'   Defaults to `"stim_nd"` so that traces with the same neutral-density level
#'   are averaged together regardless of wavelength — matching the Origin
#'   StimResp grouping for spectral (C0) protocols.  Use `"stim_label"` to keep
#'   each wavelength × ND combination as a separate trace.
#' @param overall_color Color for the overall mean line.
#' @param marker_color Color for marker reference lines.
#' @param base_size Base font size for the plot theme.
#'
#' @return A ggplot object.
#' @export
zeus_plot_mean_waveform <- function(
    x,
    data_slot = c("traces_70", "traces_280"),
    channel_filter = NULL,
    compare_raw = FALSE,
    include_overall = TRUE,
    overlay_markers = FALSE,
    include_photocell = TRUE,
    photocell_filter = "Photocell",
    photocell_color = "black",
    photocell_auto_scale = TRUE,
    photocell_relative_height = 1.1,
    a_window = c(400, 700),
    b_window = c(400, 700),
    d_window = c(700, 1000),
    stim_levels = NULL,
    color_by = c("stim_nd", "stim_label", "wavelength"),
    overall_color = "red",
    marker_color = "black",
    base_size = 12
) {
  data_slot <- match.arg(data_slot)
  color_by <- match.arg(color_by)

  prepared <- .zeus_prepare_plot_df(
    x = x,
    data_slot = data_slot,
    channel_filter = channel_filter,
    require_cols = c("time", "value"),
    allow_stimresp = TRUE
  )

  df_plot <- prepared$df_plot
  df_photocell_source <- prepared$df_photocell_source

  if (!("stim_label" %in% names(df_plot))) {
    if ("stim_nd" %in% names(df_plot)) {
      df_plot <- df_plot |>
        dplyr::mutate(stim_label = as.character(.data$stim_nd))
    } else {
      stop(
        "Plot data must contain either `stim_label` or `stim_nd`.",
        call. = FALSE
      )
    }
  }

  if (isTRUE(compare_raw) && !"value_raw" %in% names(df_plot)) {
    stop(
      "`compare_raw = TRUE` requires a `value_raw` column in the plotted data.",
      call. = FALSE
    )
  }

  # Legend grouping
  df_plot <- df_plot |>
    dplyr::mutate(
      color_group = dplyr::case_when(
        color_by == "stim_label" ~ as.character(.data$stim_label),
        color_by == "stim_nd" ~ as.character(.data$stim_nd),
        color_by == "wavelength" ~ as.character(.data$wavelength),
        TRUE ~ as.character(.data$stim_label)
      )
    )

  if (is.null(stim_levels)) {
    stim_levels_chr <- unique(df_plot$color_group)
  } else {
    stim_levels_chr <- as.character(stim_levels)
  }

  df_plot <- df_plot |>
    dplyr::mutate(
      color_group = factor(.data$color_group, levels = stim_levels_chr)
    )

  # Sequential palette:

  plot_colors <- .zeus_waveform_palette(length(stim_levels_chr))
  names(plot_colors) <- stim_levels_chr

  # Aggregate plotted data
  # Only include metadata cols that are at the SAME level as color_by to
  # avoid creating spurious fine-grained sub-groups (e.g. when color_by =
  # "stim_nd", including stim_label in the grouping would create 70 separate
  # lines for C0 instead of one averaged line per ND level).
  extra_meta_cols <- switch(
    color_by,
    stim_label = intersect(c("stim_nd", "wavelength"), names(df_plot)),
    stim_nd    = intersect(c("stim_nd"),                names(df_plot)),
    wavelength = intersect(c("wavelength"),              names(df_plot)),
    character(0)
  )

  df_by_stim <- df_plot |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c("time_ms", "color_group", extra_meta_cols)))) |>
    dplyr::summarise(
      Smoothed = mean(.data$value, na.rm = TRUE),
      Raw = if (isTRUE(compare_raw)) mean(.data$value_raw, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )

  if (isTRUE(compare_raw)) {
    df_by_stim <- df_by_stim |>
      tidyr::pivot_longer(
        cols = c("Raw", "Smoothed"),
        names_to = "signal_type",
        values_to = "signal"
      )
  } else {
    df_by_stim <- df_by_stim |>
      dplyr::rename(signal = .data$Smoothed) |>
      dplyr::mutate(signal_type = "Smoothed") |>
      dplyr::select(-.data$Raw)
  }

  # Overall mean
  df_overall <- NULL

  if (isTRUE(include_overall)) {
    df_overall <- df_plot |>
      dplyr::group_by(.data$time_ms) |>
      dplyr::summarise(
        Smoothed = mean(.data$value, na.rm = TRUE),
        Raw = if (isTRUE(compare_raw)) mean(.data$value_raw, na.rm = TRUE) else NA_real_,
        .groups = "drop"
      )

    if (isTRUE(compare_raw)) {
      df_overall <- df_overall |>
        tidyr::pivot_longer(
          cols = c("Raw", "Smoothed"),
          names_to = "signal_type",
          values_to = "signal"
        )
    } else {
      df_overall <- df_overall |>
        dplyr::rename(signal = .data$Smoothed) |>
        dplyr::mutate(signal_type = "Smoothed") |>
        dplyr::select(-.data$Raw)
    }

    df_overall <- df_overall |>
      dplyr::mutate(
        color_group = factor("Overall Mean", levels = c(stim_levels_chr, "Overall Mean"))
      )
  }

  # Photocell — Origin-like: scaled to fraction of ERG range, positioned below
  df_photocell <- NULL

  if (isTRUE(include_photocell) && !is.null(df_photocell_source)) {
    erg_ymin  <- min(df_by_stim$signal, na.rm = TRUE)
    erg_ymax  <- max(df_by_stim$signal, na.rm = TRUE)
    erg_range <- erg_ymax - erg_ymin

    if (isTRUE(photocell_auto_scale)) {
      df_photocell <- .zeus_scale_photocell(
        df_photocell_source      = df_photocell_source,
        photocell_filter         = photocell_filter,
        erg_ymin                 = erg_ymin,
        erg_range                = erg_range,
        photocell_relative_height = photocell_relative_height
      )
    } else {
      # Manual path: still extract and average, but do not rescale.
      pc <- df_photocell_source
      if (!("time_ms" %in% names(pc)) && "time" %in% names(pc)) {
        pc <- pc |> dplyr::mutate(time_ms = zeus_time_to_ms(.data$time))
      }
      if (!("value" %in% names(pc)) && "signal" %in% names(pc)) {
        pc <- pc |> dplyr::rename(value = .data$signal)
      }
      if (all(c("time_ms", "value") %in% names(pc))) {
        if ("channel" %in% names(pc)) {
          pc <- pc |> dplyr::filter(grepl(photocell_filter, .data$channel, fixed = TRUE))
        }
        if (nrow(pc) > 0) {
          df_photocell <- pc |>
            dplyr::group_by(.data$time_ms) |>
            dplyr::summarise(signal = mean(.data$value, na.rm = TRUE), .groups = "drop")
        }
      }
    }
  }

  # Build plot
  p <- ggplot2::ggplot()

  if (isTRUE(compare_raw)) {
    p <- p +
      ggplot2::geom_line(
        data = df_by_stim,
        ggplot2::aes(
          x = .data$time_ms,
          y = .data$signal,
          color = .data$color_group,
          alpha = .data$signal_type,
          group = interaction(.data$color_group, .data$signal_type)
        ),
        linewidth = 0.45,
        lineend = "round"
      ) +
      ggplot2::scale_alpha_manual(
        name = "Signal",
        values = c("Raw" = 0.30, "Smoothed" = 0.95),
        breaks = c("Raw", "Smoothed")
      )
  } else {
    p <- p +
      ggplot2::geom_line(
        data = df_by_stim,
        ggplot2::aes(
          x = .data$time_ms,
          y = .data$signal,
          color = .data$color_group,
          group = .data$color_group
        ),
        linewidth = 0.50,
        alpha = 0.95,
        lineend = "round"
      )
  }

  if (isTRUE(include_overall) && !is.null(df_overall)) {
    if (isTRUE(compare_raw)) {
      p <- p +
        ggplot2::geom_line(
          data = dplyr::filter(df_overall, .data$signal_type == "Raw"),
          ggplot2::aes(x = .data$time_ms, y = .data$signal),
          color = overall_color,
          linewidth = 1.00,
          alpha = 0.30,
          inherit.aes = FALSE,
          show.legend = FALSE,
          lineend = "round"
        )
    }

    p <- p +
      ggplot2::geom_line(
        data = dplyr::filter(df_overall, .data$signal_type == "Smoothed"),
        ggplot2::aes(
          x = .data$time_ms,
          y = .data$signal,
          color = .data$color_group
        ),
        linewidth = 1.35,
        alpha = 1,
        inherit.aes = FALSE,
        lineend = "round"
      )
  }

  if (!is.null(df_photocell)) {
    p <- p +
      ggplot2::geom_line(
        data = df_photocell,
        ggplot2::aes(
          x = .data$time_ms,
          y = .data$signal,
          color = "Photocell"
        ),
        linewidth = 0.30,
        alpha = 0.95,
        inherit.aes = FALSE,
        lineend = "round"
      )
  }

  # Marker overlay
  if (isTRUE(overlay_markers) && isTRUE(include_overall) && !is.null(df_overall)) {
    overall_for_markers <- df_overall |>
      dplyr::filter(.data$signal_type == "Smoothed") |>
      dplyr::select(.data$time_ms, .data$signal) |>
      dplyr::arrange(.data$time_ms)

    if (nrow(overall_for_markers) > 0) {
      a_marker <- overall_for_markers |>
        dplyr::filter(.data$time_ms >= a_window[1], .data$time_ms <= a_window[2])

      if (nrow(a_marker) > 0) {
        a_marker <- a_marker |>
          dplyr::slice_min(order_by = .data$signal, n = 1, with_ties = FALSE) |>
          dplyr::transmute(
            marker_type = "A-wave trough",
            time_ms = .data$time_ms
          )
      } else {
        a_marker <- tibble::tibble(marker_type = "A-wave trough", time_ms = NA_real_)
      }

      b_marker <- tibble::tibble(marker_type = "B-wave peak", time_ms = NA_real_)

      if (!is.na(a_marker$time_ms[1])) {
        b_candidates <- overall_for_markers |>
          dplyr::filter(
            .data$time_ms >= max(b_window[1], a_marker$time_ms[1]),
            .data$time_ms <= b_window[2]
          )

        if (nrow(b_candidates) > 0) {
          b_marker <- b_candidates |>
            dplyr::slice_max(order_by = .data$signal, n = 1, with_ties = FALSE) |>
            dplyr::transmute(
              marker_type = "B-wave peak",
              time_ms = .data$time_ms
            )
        }
      }

      if (is.na(b_marker$time_ms[1])) {
        b_candidates <- overall_for_markers |>
          dplyr::filter(.data$time_ms >= b_window[1], .data$time_ms <= b_window[2])

        if (nrow(b_candidates) > 0) {
          b_marker <- b_candidates |>
            dplyr::slice_max(order_by = .data$signal, n = 1, with_ties = FALSE) |>
            dplyr::transmute(
              marker_type = "B-wave peak",
              time_ms = .data$time_ms
            )
        }
      }

      d_marker <- tibble::tibble(marker_type = "D-wave peak", time_ms = NA_real_)

      d_candidates <- overall_for_markers |>
        dplyr::filter(.data$time_ms >= d_window[1], .data$time_ms <= d_window[2])

      if (nrow(d_candidates) > 0) {
        d_window_width <- diff(d_window)
        d_trough_end <- d_window[1] + (d_window_width * 0.4)

        d_trough_candidates <- d_candidates |>
          dplyr::filter(.data$time_ms <= d_trough_end)

        if (nrow(d_trough_candidates) > 0) {
          d_trough <- d_trough_candidates |>
            dplyr::slice_min(order_by = .data$signal, n = 1, with_ties = FALSE)

          d_peak_candidates <- d_candidates |>
            dplyr::filter(.data$time_ms > d_trough$time_ms)

          if (nrow(d_peak_candidates) > 0) {
            d_marker <- d_peak_candidates |>
              dplyr::slice_max(order_by = .data$signal, n = 1, with_ties = FALSE) |>
              dplyr::transmute(
                marker_type = "D-wave peak",
                time_ms = .data$time_ms
              )
          }
        }

        if (is.na(d_marker$time_ms[1])) {
          d_marker <- d_candidates |>
            dplyr::slice_max(order_by = .data$signal, n = 1, with_ties = FALSE) |>
            dplyr::transmute(
              marker_type = "D-wave peak",
              time_ms = .data$time_ms
            )
        }
      }

      marker_lines <- dplyr::bind_rows(a_marker, b_marker, d_marker) |>
        dplyr::filter(!is.na(.data$time_ms)) |>
        dplyr::mutate(
          marker_type = factor(
            .data$marker_type,
            levels = c("A-wave trough", "B-wave peak", "D-wave peak")
          )
        )

      if (nrow(marker_lines) > 0) {
        p <- p +
          ggplot2::geom_vline(
            data = marker_lines,
            ggplot2::aes(
              xintercept = .data$time_ms,
              linetype = .data$marker_type
            ),
            color = marker_color,
            linewidth = 0.8,
            alpha = 0.95,
            inherit.aes = FALSE
          ) +
          ggplot2::scale_linetype_manual(
            name = "Reference lines",
            values = c(
              "A-wave trough" = "dotted",
              "B-wave peak" = "longdash",
              "D-wave peak" = "dotdash"
            )
          )
      }
    }
  }

  # Final formatting
  color_values <- c(plot_colors)
  color_breaks <- stim_levels_chr

  if (isTRUE(include_overall)) {
    color_values <- c(color_values, "Overall Mean" = overall_color)
    color_breaks <- c(color_breaks, "Overall Mean")
  }

  if (!is.null(df_photocell)) {
    color_values <- c(color_values, "Photocell" = photocell_color)
    color_breaks <- c(color_breaks, "Photocell")
  }

  p +
    ggplot2::scale_color_manual(
      values = color_values,
      breaks = color_breaks,
      drop = FALSE
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        title = "Trace Key",
        order = 1,
        override.aes = list(
          linewidth = rep(1.35, length(color_breaks)),
          alpha = 1,
          linetype = "solid"
        )
      ),
      linetype = ggplot2::guide_legend(
        title = "Reference lines",
        order = 2,
        override.aes = list(
          color = marker_color,
          linewidth = 0.9,
          alpha = 1
        )
      ),
      alpha = ggplot2::guide_legend(
        title = "Signal",
        order = 3
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 6),
      expand = ggplot2::expansion(mult = c(0.01, 0.03))
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = 8),
      minor_breaks = function(lims, ...) scales::pretty_breaks(n = 14)(lims),
      expand = ggplot2::expansion(mult = c(0.03, 0.08))
    ) +
    ggplot2::labs(
      title = "Mean ERG Waveform",
      x = "Time (ms)",
      y = "Response (\u00B5V)"
    ) +
    .zeus_waveform_theme(base_size = base_size)
}


# Intensity Response by ND -----------------------------------------------------

#' Plot ERG intensity-response curves for A-, B-, and D-wave features
#'
#' @description
#' This function generates a multi-panel intensity-response plot for ERG data,
#' displaying A-wave trough amplitude, B-wave peak amplitude, and D-wave peak
#' amplitude as a function of stimulus ND.
#'
#' The function operates on a `zeus_stimresp` object returned by
#' [zeus_read_abf()], extracting per-trace waveform measurements via
#' [extract_irrad_wl_amp()] before computing summary statistics.
#'
#' The resulting plot shows mean response values for each wave type across
#' stimulus ND levels, optionally grouped by a categorical variable such as
#' treatment group. Standard error bars can also be included.
#'
#' @param x A `zeus_stimresp` object (output of [zeus_read_abf()]).
#' @param group_col Optional character string specifying a grouping variable
#'   present in `x$traces_70` (e.g., \code{"treatment_group"}). If provided,
#'   separate curves are plotted for each group.
#' @param use_se Logical; if \code{TRUE}, standard error bars (mean ± SE) are
#'   displayed. Default is \code{TRUE}.
#' @param base_size Numeric base font size for the plot theme. Default is 11.
#'
#' @details
#' The three panels correspond to:
#' \itemize{
#'   \item A-wave trough amplitude (`awave_mv`)
#'   \item B-wave peak amplitude (`amp_mv`)
#'   \item D-wave peak amplitude (`dwave_mv`)
#' }
#'
#' Stimulus ND is displayed on the x-axis in descending order so that larger ND
#' values (dimmer stimuli) appear first.
#'
#' @return A \code{ggplot2} object showing intensity-response curves for each
#'   wave type.
#'
#' @export
zeus_plot_intensity_response <- function(
    x,
    group_col = NULL,
    use_se = TRUE,
    base_size = 11
) {
  if (!inherits(x, "zeus_stimresp")) {
    stop("`x` must be a `zeus_stimresp` object returned by `zeus_read_abf()`.",
         call. = FALSE)
  }

  # Extract per-trace waveform measurements using the current API
  features_df <- extract_irrad_wl_amp(x)

  if (!is.null(group_col) && !group_col %in% names(features_df)) {
    stop("`group_col` '", group_col, "' not found in extracted features.", call. = FALSE)
  }

  # Reshape: A-wave = awave_mv, B-wave = amp_mv, D-wave = dwave_mv
  select_cols <- c("stim_nd", "awave_mv", "amp_mv", "dwave_mv")

  if (!is.null(group_col)) {
    select_cols <- c(select_cols, group_col)
  }

  df_plot <- features_df |>
    dplyr::select(dplyr::any_of(select_cols)) |>
    tidyr::pivot_longer(
      cols = dplyr::any_of(c("awave_mv", "amp_mv", "dwave_mv")),
      names_to = "wave",
      values_to = "response"
    ) |>
    dplyr::mutate(
      wave = dplyr::case_match(
        .data$wave,
        "awave_mv" ~ "A-wave",
        "amp_mv"   ~ "B-wave",
        "dwave_mv" ~ "D-wave"
      )
    ) |>
    dplyr::filter(!is.na(.data$response), !is.na(.data$stim_nd))

  # Summary stats
  group_vars <- c("stim_nd", "wave")
  if (!is.null(group_col)) group_vars <- c(group_vars, group_col)

  df_summary <- df_plot |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(.data$response, na.rm = TRUE),
      se = stats::sd(.data$response, na.rm = TRUE) / sqrt(n),
      .groups = "drop"
    )

  # Plot
  p <- ggplot2::ggplot(
    df_summary,
    ggplot2::aes(x = .data$stim_nd, y = .data$mean)
  )

  if (!is.null(group_col)) {
    p <- p +
      ggplot2::aes(
        color = .data[[group_col]],
        group = interaction(.data[[group_col]], .data$wave)
      )
  } else {
    p <- p +
      ggplot2::aes(group = .data$wave)
  }

  p <- p +
    ggplot2::geom_line(linewidth = 1.0, lineend = "round") +
    ggplot2::geom_point(size = 2.2)

  if (isTRUE(use_se)) {
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = .data$mean - .data$se,
          ymax = .data$mean + .data$se
        ),
        width = 0.1,
        linewidth = 0.4
      )
  }

  p +
    ggplot2::facet_wrap(~wave, nrow = 1, scales = "free_y") +
    ggplot2::scale_x_reverse() +
    ggplot2::labs(
      x = "Stimulus ND",
      y = "Response (\u00B5V)",
      title = "Wave Intensity Response",
      color = if (!is.null(group_col)) group_col else NULL
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = if (!is.null(group_col)) "right" else "none",
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.line = ggplot2::element_line(linewidth = 0.5),
      axis.ticks = ggplot2::element_line(linewidth = 0.4),
      panel.spacing = grid::unit(1.2, "lines")
    )
}


# Spectral waveform panel (10 wavelength blocks) --------------------------------

#' Plot ERG waveforms split into per-wavelength-block facets
#'
#' @description
#' Generates a publication-ready multi-panel waveform plot in which each panel
#' corresponds to one wavelength block of the C0 (spectral) protocol — or, for
#' single-wavelength (C1/white-light) data, one run block.  Within each panel,
#' individual traces for each neutral-density (ND) level are drawn and
#' color-coded using a perceptually uniform sequential palette (lightest = most
#' attenuated / highest ND; darkest = least attenuated / lowest ND).
#'
#' This layout mirrors the Origin "spectral waveform" panel: the 10 wavelength
#' blocks form a 2-row × 5-column grid. When photocell data are available a
#' representative photocell pulse is overlaid at the bottom of every panel,
#' positioned and scaled to match the Origin StimResp display convention.
#'
#' The time scale is shown underneath every individual panel (requires
#' ggplot2 >= 3.4.0). Each panel's strip label is composed of two lines:
#' \enumerate{
#'   \item The wavelength prefix (C0) or block number (C1).
#'   \item A metadata tag of the form `"Block N | C0/C1 | filename"`.
#' }
#' When `x` is a `zeus_stimresp` object the protocol and filename are
#' auto-detected; they can also be supplied explicitly via `protocol_label`
#' and `file_label`.
#'
#' @param x A `zeus_stimresp` object (output of [zeus_read_abf()]) or a
#'   long-format waveform data frame containing at least `time`, `value`,
#'   `stim_label`, `stim_nd`, and `block_index`.
#' @param data_slot For `zeus_stimresp` objects, which component to use.
#'   Default is `"traces_70"` (averaged traces). Can also be `"traces_280"`.
#' @param channel_filter Optional channel to filter. Ignored when the selected
#'   data slot has no `channel` column.
#' @param stim_levels Optional numeric vector giving the desired ND levels to
#'   display and their legend order (highest = dimmest first).  If `NULL`,
#'   all ND levels present in the data are used in descending order.
#' @param include_photocell Logical; if `TRUE` (default), overlays the mean
#'   photocell trace on every panel when photocell data are available.
#' @param photocell_filter Character pattern used to identify the photocell
#'   channel when a `channel` column is present. Default is `"Photocell"`.
#' @param photocell_color Color for the photocell overlay line. Default
#'   `"black"`.
#' @param photocell_relative_height Fraction of the global ERG amplitude range
#'   that the photocell pulse occupies. Default is `0.20`.
#' @param protocol_label Character string shown in every panel label to
#'   identify the recording protocol (e.g. `"C0"` or `"C1"`). When `x` is a
#'   `zeus_stimresp` object this is auto-detected from `x$protocol`; supply
#'   this argument to override or when passing a plain data frame.
#' @param file_label Character string shown in every panel label to identify
#'   the source file (the `.abf` extension is excluded). When `x` is a
#'   `zeus_stimresp` object this is auto-detected from `x$path`; supply this
#'   argument to override or when passing a plain data frame.
#' @param facet_ncol Number of columns in the facet grid. Default is `5`.
#' @param base_size Base font size for the plot theme. Default is `10`.
#'
#' @return A `ggplot2` object.
#' @export
zeus_plot_spectral_waveform <- function(
    x,
    data_slot = c("traces_70", "traces_280"),
    channel_filter = NULL,
    stim_levels = NULL,
    include_photocell = TRUE,
    photocell_filter = "Photocell",
    photocell_color = "black",
    photocell_relative_height = 1.1,
    protocol_label = NULL,
    file_label = NULL,
    facet_ncol = 5L,
    base_size = 10
) {
  data_slot  <- match.arg(data_slot)
  facet_ncol <- as.integer(facet_ncol)

  prepared <- .zeus_prepare_plot_df(
    x = x,
    data_slot = data_slot,
    channel_filter = channel_filter,
    require_cols = c("time", "value"),
    allow_stimresp = TRUE
  )

  df_plot              <- prepared$df_plot
  df_photocell_source  <- prepared$df_photocell_source

  if (!("stim_label" %in% names(df_plot))) {
    stop("`x` must contain `stim_label`.", call. = FALSE)
  }

  if (!("stim_nd" %in% names(df_plot))) {
    stop("`x` must contain `stim_nd`.", call. = FALSE)
  }

  # Extract the wavelength prefix: first non-whitespace token of stim_label.
  # Expected formats: "650A 4.0" -> "650A", "570 5.0" -> "570",
  #                   "White 6.0" -> "White".
  df_plot <- df_plot |>
    dplyr::mutate(
      .wl_prefix = stringr::str_extract(.data$stim_label, "^\\S+")
    )

  unique_prefixes <- unique(df_plot$.wl_prefix)

  if (length(unique_prefixes) > 1L) {
    # C0-like: multiple wavelength prefixes — one per block already.
    if ("block_index" %in% names(df_plot)) {
      prefix_block_map <- df_plot |>
        dplyr::distinct(.data$.wl_prefix, .data$block_index) |>
        dplyr::arrange(.data$block_index)
      compound_labels <- prefix_block_map$.wl_prefix
      names(compound_labels) <- prefix_block_map$.wl_prefix
      label_order <- prefix_block_map$.wl_prefix
    } else {
      compound_labels <- unique_prefixes
      names(compound_labels) <- unique_prefixes
      label_order <- unique_prefixes
    }

    df_plot <- df_plot |>
      dplyr::mutate(
        block_label = factor(
          compound_labels[.data$.wl_prefix],
          levels = unique(compound_labels[label_order])
        )
      )
  } else {
    # C1-like: single wavelength prefix, keep labels concise while still
    # distinguishing repeated runs.
    if (!("block_index" %in% names(df_plot))) {
      stop(
        "Cannot determine facet panels: `block_index` column is required ",
        "when all stim_labels share the same wavelength prefix.",
        call. = FALSE
      )
    }
    block_order <- sort(unique(df_plot$block_index))
    wl_label <- unique_prefixes[[1]]

    block_label_vec <- vapply(
      block_order,
      function(b) paste0(wl_label, " Run ", b),
      character(1L)
    )
    names(block_label_vec) <- as.character(block_order)

    df_plot <- df_plot |>
      dplyr::mutate(
        block_label = factor(
          block_label_vec[as.character(.data$block_index)],
          levels = block_label_vec
        )
      )
  }

  # ND levels for color ordering (highest ND = dimmest stimulus shown first).
  if (is.null(stim_levels)) {
    nd_levels <- sort(unique(df_plot$stim_nd), decreasing = TRUE)
  } else {
    nd_levels <- as.numeric(stim_levels)
  }
  nd_levels_chr <- as.character(nd_levels)

  # Aggregate: mean across any repeated rows for the same (block, nd, time).
  df_agg <- df_plot |>
    dplyr::mutate(
      nd_group = factor(as.character(.data$stim_nd), levels = nd_levels_chr)
    ) |>
    dplyr::group_by(.data$block_label, .data$nd_group, .data$time_ms) |>
    dplyr::summarise(
      signal = mean(.data$value, na.rm = TRUE),
      .groups = "drop"
    )

  # Publication-quality sequential palette:
  # highest ND (dimmest stimulus) → lightest blue;
  # lowest  ND (brightest stimulus) → darkest navy.
  n_nd <- length(nd_levels_chr)
  plot_colors <- stats::setNames(.zeus_waveform_palette(n_nd), nd_levels_chr)

  # Photocell — scale once globally, then replicate across all panels.
  df_photocell_panels <- NULL

  if (isTRUE(include_photocell) && !is.null(df_photocell_source)) {
    erg_ymin  <- min(df_agg$signal, na.rm = TRUE)
    erg_ymax  <- max(df_agg$signal, na.rm = TRUE)
    erg_range <- erg_ymax - erg_ymin

    df_pc_scaled <- .zeus_scale_photocell(
      df_photocell_source       = df_photocell_source,
      photocell_filter          = photocell_filter,
      erg_ymin                  = erg_ymin,
      erg_range                 = erg_range,
      photocell_relative_height = photocell_relative_height
    )

    if (!is.null(df_pc_scaled)) {
      # tidyr::crossing() efficiently duplicates the scaled photocell trace for
      # every facet panel so geom_line draws it in each subplot.
      panel_levels <- levels(df_agg$block_label)
      df_photocell_panels <- tidyr::crossing(
        df_pc_scaled,
        block_label = factor(panel_levels, levels = panel_levels)
      )
    }
  }

  # Facet control: show x-axis on ALL panels (ggplot2 >= 3.4.0).
  # Construct the facet call defensively; if the installed ggplot2 is older
  # the function falls back silently and x-axes are shown only on the bottom
  # row of the grid.
  facet_layer <- tryCatch(
    ggplot2::facet_wrap(~block_label, ncol = facet_ncol, axes = "all_x"),
    error = function(e) {
      message(
        "zeus_plot_spectral_waveform: `axes = 'all_x'` requires ggplot2 >= 3.4.0; ",
        "falling back to x-axis on bottom row only."
      )
      ggplot2::facet_wrap(~block_label, ncol = facet_ncol)
    }
  )

  # Build the plot -----------------------------------------------------------
  p <- ggplot2::ggplot(
    df_agg,
    ggplot2::aes(
      x     = .data$time_ms,
      y     = .data$signal,
      color = .data$nd_group,
      group = .data$nd_group
    )
  ) +
    ggplot2::geom_line(linewidth = 0.45, lineend = "round", alpha = 0.95) +
    facet_layer

  if (!is.null(df_photocell_panels)) {
    p <- p +
      ggplot2::geom_line(
        data = df_photocell_panels,
        ggplot2::aes(
          x = .data$time_ms,
          y = .data$signal,
          color = "Photocell"
        ),
        linewidth   = 0.28,
        alpha       = 0.85,
        inherit.aes = FALSE,
        lineend     = "round"
      )
  }

  color_values <- plot_colors
  if (!is.null(df_photocell_panels)) {
    color_values <- c(color_values, "Photocell" = photocell_color)
  }

  p +
    ggplot2::scale_color_manual(
      values = color_values,
      name   = "Trace Key",
      drop   = FALSE
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        title = "Trace Key",
        override.aes = list(linewidth = 1.2, alpha = 1),
        order = 1
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 4),
      expand = ggplot2::expansion(mult = c(0.01, 0.03))
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = 8),
      minor_breaks = function(lims, ...) scales::pretty_breaks(n = 14)(lims),
      expand = ggplot2::expansion(mult = c(0.05, 0.08))
    ) +
    ggplot2::labs(
      title = "Spectral ERG Waveforms by Wavelength Block",
      x     = "Time (ms)",
      y     = "Response (\u00B5V)"
    ) +
    .zeus_waveform_theme(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "#F3F6F8", color = "#C7D3DD", linewidth = 0.5),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(0.92),
        colour = "#1F2A33",
        lineheight = 0.98
      ),
      panel.spacing.x = grid::unit(0.6, "lines"),
      panel.spacing.y = grid::unit(0.8, "lines")
    )
}
