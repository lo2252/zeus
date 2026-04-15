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
#' @param photocell_color Color for the photocell trace.
#' @param photocell_auto_scale Logical; if `TRUE`, rescales the photocell trace
#'   relative to the ERG signal.
#' @param photocell_relative_height Relative height multiplier used when
#'   `photocell_auto_scale = TRUE`.
#' @param a_window Numeric length-2 vector giving the A-wave search interval in
#'   milliseconds.
#' @param b_window Numeric length-2 vector giving the B-wave search interval in
#'   milliseconds.
#' @param d_window Numeric length-2 vector giving the D-wave search interval in
#'   milliseconds.
#' @param stim_levels Optional vector giving the desired order of stimulus labels
#'   in the legend. If `NULL`, uses the order present in the data.
#' @param color_by One of `"stim_label"`, `"stim_nd"`, or `"wavelength"`.
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
    photocell_relative_height = 1.0,
    a_window = c(400, 700),
    b_window = c(400, 700),
    d_window = c(700, 1000),
    stim_levels = NULL,
    color_by = c("stim_label", "stim_nd", "wavelength"),
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

  # Color palette
  base_palette <- c(
    "#332288", "#88CCEE", "#44AA99", "#117733",
    "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"
  )
  plot_colors <- rep(base_palette, length.out = length(stim_levels_chr))
  names(plot_colors) <- stim_levels_chr

  # Aggregate plotted data
  group_cols <- c("time_ms", "color_group")
  optional_cols <- intersect(c("stim_label", "stim_nd", "wavelength"), names(df_plot))

  df_by_stim <- df_plot |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, optional_cols)))) |>
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

  # Photocell
  df_photocell <- NULL

  if (isTRUE(include_photocell) && !is.null(df_photocell_source)) {
    if (is.data.frame(df_photocell_source)) {
      if (!("time" %in% names(df_photocell_source))) {
        if ("time_ms" %in% names(df_photocell_source)) {
          df_photocell_source <- df_photocell_source |>
            dplyr::mutate(time = .data$time_ms / 1000)
        }
      }

      if (!("value" %in% names(df_photocell_source))) {
        if ("signal" %in% names(df_photocell_source)) {
          df_photocell_source <- df_photocell_source |>
            dplyr::rename(value = .data$signal)
        } else if ("mean_value" %in% names(df_photocell_source)) {
          df_photocell_source <- df_photocell_source |>
            dplyr::rename(value = .data$mean_value)
        }
      }

      if (!("time_ms" %in% names(df_photocell_source)) && "time" %in% names(df_photocell_source)) {
        df_photocell_source <- df_photocell_source |>
          dplyr::mutate(time_ms = zeus_time_to_ms(.data$time))
      }

      if (all(c("time", "value") %in% names(df_photocell_source))) {
        if ("channel" %in% names(df_photocell_source)) {
          df_photocell <- df_photocell_source |>
            dplyr::filter(grepl(photocell_filter, .data$channel))
        } else {
          df_photocell <- df_photocell_source
        }

        if (nrow(df_photocell) > 0) {
          df_photocell <- df_photocell |>
            dplyr::group_by(.data$time_ms) |>
            dplyr::summarise(
              signal = mean(.data$value, na.rm = TRUE),
              .groups = "drop"
            )

          if (isTRUE(photocell_auto_scale)) {
            erg_peak <- max(abs(df_by_stim$signal), na.rm = TRUE)
            photocell_peak <- max(abs(df_photocell$signal), na.rm = TRUE)

            if (is.finite(erg_peak) && is.finite(photocell_peak) && photocell_peak > 0) {
              scale_factor <- (erg_peak * photocell_relative_height) / photocell_peak
            } else {
              scale_factor <- 1
            }

            df_photocell <- df_photocell |>
              dplyr::mutate(signal = .data$signal * scale_factor)
          }
        } else {
          df_photocell <- NULL
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
          group = if ("stim_label" %in% names(df_by_stim))
            interaction(.data$stim_label, .data$color_group, .data$signal_type)
          else
            interaction(.data$color_group, .data$signal_type)
        ),
        linewidth = 0.55,
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
          group = if ("stim_label" %in% names(df_by_stim))
            interaction(.data$stim_label, .data$color_group)
          else
            .data$color_group
        ),
        linewidth = 0.60,
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
        ggplot2::aes(x = .data$time_ms, y = .data$signal),
        color = photocell_color,
        linewidth = 0.5,
        alpha = 0.95,
        inherit.aes = FALSE,
        show.legend = FALSE,
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
  color_values <- c(plot_colors, "Overall Mean" = overall_color)
  color_breaks <- c(stim_levels_chr, if (isTRUE(include_overall)) "Overall Mean")

  p +
    ggplot2::scale_color_manual(
      values = color_values,
      breaks = color_breaks,
      drop = FALSE
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        title = switch(
          color_by,
          stim_label = "Stimulus",
          stim_nd = "Stimulus ND",
          wavelength = "Wavelength"
        ),
        order = 1,
        override.aes = list(
          linewidth = c(rep(1.1, length(color_breaks) - 1), 1.6),
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
      breaks = scales::pretty_breaks(n = 6),
      expand = ggplot2::expansion(mult = c(0.03, 0.08))
    ) +
    ggplot2::labs(
      title = "Mean ERG Waveform",
      x = "Time (ms)",
      y = expression("Response (" * mu * "V)")
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 10)
      ),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = grid::unit(0.25, "cm"),
      legend.key.width = grid::unit(0.9, "cm"),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(linewidth = 0.6),
      axis.ticks = ggplot2::element_line(linewidth = 0.5)
    )
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
      se = stats::sd(.data$response, na.rm = TRUE) / sqrt(.data$n),
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
      y = expression("Response (" * mu * "V)"),
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
