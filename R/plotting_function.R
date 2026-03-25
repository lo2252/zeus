#' Plot mean ERG waveform across sweeps with optional reference lines
#'
#' Computes and visualizes the mean waveform across sweeps for each stimulus
#' level, along with an overall mean waveform across all sweeps.
#'
#' Stimulus-specific mean waveforms are drawn as thinner color-coded lines,
#' while the overall mean waveform is drawn as a thicker red line by default.
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
#' @param df_long ZEUS long-format data frame.
#' @param channel_filter Optional channel to filter (for example, `"ERG DAM80"`).
#' @param compare_raw Logical; if `TRUE`, plots both raw and smoothed means.
#'   Requires `value_raw` in `df_long`. Default is `FALSE`.
#' @param include_overall Logical; if `TRUE`, includes an overall mean waveform
#'   across all sweeps.
#' @param overlay_markers Logical; if `TRUE`, overlays vertical reference lines
#'   at the A-wave trough, B-wave peak, and D-wave peak locations derived from
#'   the overall mean waveform.
#' @param a_window Numeric length-2 vector giving the A-wave search interval in
#'   the same units as `df_long$time`.
#' @param b_window Numeric length-2 vector giving the B-wave search interval in
#'   the same units as `df_long$time`.
#' @param d_window Numeric length-2 vector giving the D-wave search interval in
#'   the same units as `df_long$time`.
#' @param stim_levels Optional vector giving the desired order of `stim_nd`
#'   values in the legend. If `NULL`, numeric values are ordered ascending.
#' @param overall_color Color for the overall mean line. Default is `"red"`.
#' @param marker_color Color for the marker reference lines. Default is
#'   `"black"`.
#' @param base_size Base font size for the plot theme. Default is `12`.
#'
#' @return A ggplot object.
#' @export
zeus_plot_mean_waveform <- function(
    df_long,
    channel_filter = NULL,
    compare_raw = FALSE,
    include_overall = TRUE,
    overlay_markers = TRUE,
    a_window = c(0.4, 0.7),
    b_window = c(0.4, 0.7),
    d_window = c(0.7, 1.0),
    stim_levels = NULL,
    overall_color = "red",
    marker_color = "black",
    base_size = 12
) {

  # 1. Input validation
  required_cols <- c("time", "value", "stim_nd")
  missing_cols <- setdiff(required_cols, names(df_long))

  if (length(missing_cols) > 0) {
    stop(
      "df_long is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  df_plot <- df_long

  if (!is.null(channel_filter)) {
    if (!"channel" %in% names(df_plot)) {
      stop(
        "`channel_filter` provided but no `channel` column found.",
        call. = FALSE
      )
    }

    df_plot <- df_plot |>
      dplyr::filter(.data$channel == channel_filter)

    if (nrow(df_plot) == 0) {
      stop(
        "No data left after applying `channel_filter`.",
        call. = FALSE
      )
    }
  }

  if (isTRUE(compare_raw) && !"value_raw" %in% names(df_plot)) {
    stop(
      "`compare_raw = TRUE` requires a `value_raw` column in `df_long`.",
      call. = FALSE
    )
  }

  # 2. Stimulus ordering
  if (is.null(stim_levels)) {
    stim_unique <- unique(df_plot$stim_nd)

    if (is.numeric(stim_unique)) {
      stim_levels <- sort(stim_unique, decreasing = FALSE)
    } else {
      stim_levels <- sort(as.character(stim_unique), decreasing = FALSE)
    }
  }

  stim_levels_chr <- as.character(stim_levels)

  # 3. Color palette
  # Publication-friendly and colorblind-aware
  nd_palette <- c(
    "#332288", "#88CCEE", "#44AA99", "#117733",
    "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"
  )
  nd_colors <- rep(nd_palette, length.out = length(stim_levels_chr))
  names(nd_colors) <- stim_levels_chr

  # 4. Aggregate data
  df_by_stim <- df_plot |>
    dplyr::group_by(.data$time, .data$stim_nd) |>
    dplyr::summarise(
      Smoothed = mean(.data$value, na.rm = TRUE),
      Raw = if (isTRUE(compare_raw)) mean(.data$value_raw, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      stim_nd = factor(as.character(.data$stim_nd), levels = stim_levels_chr)
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

  if (isTRUE(include_overall)) {
    df_overall <- df_plot |>
      dplyr::group_by(.data$time) |>
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
        stim_nd = factor("Overall Mean", levels = c(stim_levels_chr, "Overall Mean"))
      )
  } else {
    df_overall <- NULL
  }

  # 5. Build plot
  p <- ggplot2::ggplot()

  if (isTRUE(compare_raw)) {
    p <- p +
      ggplot2::geom_line(
        data = df_by_stim,
        ggplot2::aes(
          x = .data$time * 1000,
          y = .data$signal,
          color = .data$stim_nd,
          alpha = .data$signal_type,
          group = interaction(.data$stim_nd, .data$signal_type)
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
          x = .data$time * 1000,
          y = .data$signal,
          color = .data$stim_nd,
          group = .data$stim_nd
        ),
        linewidth = 0.60,
        alpha = 0.95,
        lineend = "round"
      )
  }

  if (isTRUE(include_overall)) {
    if (isTRUE(compare_raw)) {
      p <- p +
        ggplot2::geom_line(
          data = dplyr::filter(df_overall, .data$signal_type == "Raw"),
          ggplot2::aes(
            x = .data$time * 1000,
            y = .data$signal
          ),
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
          x = .data$time * 1000,
          y = .data$signal,
          color = .data$stim_nd
        ),
        linewidth = 1.35,
        alpha = 1,
        inherit.aes = FALSE,
        lineend = "round"
      )
  }

  # 6. Marker overlay
  if (isTRUE(overlay_markers) && isTRUE(include_overall)) {
    overall_for_markers <- df_overall |>
      dplyr::filter(.data$signal_type == "Smoothed") |>
      dplyr::select(.data$time, .data$signal)

    get_marker <- function(data, window, type, extreme_fn) {
      data |>
        dplyr::filter(.data$time >= window[1], .data$time <= window[2]) |>
        extreme_fn(order_by = .data$signal, n = 1, with_ties = FALSE) |>
        dplyr::transmute(
          marker_type = type,
          time_ms = .data$time * 1000
        )
    }

    marker_lines <- dplyr::bind_rows(
      get_marker(overall_for_markers, a_window, "A-wave trough", dplyr::slice_min),
      get_marker(overall_for_markers, b_window, "B-wave peak", dplyr::slice_max),
      get_marker(overall_for_markers, d_window, "D-wave peak", dplyr::slice_max)
    ) |>
      dplyr::mutate(
        marker_type = factor(
          .data$marker_type,
          levels = c("A-wave trough", "B-wave peak", "D-wave peak")
        )
      )

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

  # 7. Final styling
  color_values <- c(nd_colors, "Overall Mean" = overall_color)
  color_breaks <- c(stim_levels_chr, if (isTRUE(include_overall)) "Overall Mean")

  p +
    ggplot2::scale_color_manual(
      values = color_values,
      breaks = color_breaks,
      drop = FALSE
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        title = "Stimulus ND",
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
      title = "Mean ERG Waveform by Stimulus Level",
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


# Mean Plotting by ND -----------------------------------------------------
#' Plot mean ERG waveform for each ND level with wave identification
#'
#' @param df_long ZEUS long-format data frame.
#' @param channel_filter Optional channel to filter.
#' @param compare_raw Logical; if TRUE, plots raw and smoothed means.
#' @param stim_levels Optional desired order of stim_nd values.
#' @param free_y Logical; if TRUE, each facet gets its own y scale.
#' @param overlay_markers Logical; if TRUE, add A/B/D reference lines.
#' @param a_window Numeric length-2 vector for A-wave search interval.
#' @param b_window Numeric length-2 vector for B-wave search interval.
#' @param d_window Numeric length-2 vector for D-wave search interval.
#' @param marker_color Color for reference lines.
#' @param base_size Base font size.
#'
#' @return A ggplot object.
#' @export
zeus_plot_mean_by_nd <- function(
    df_long,
    channel_filter = NULL,
    compare_raw = FALSE,
    stim_levels = NULL,
    free_y = FALSE,
    overlay_markers = TRUE,
    a_window = c(0.4, 0.7),
    b_window = c(0.4, 0.7),
    d_window = c(0.7, 1.0),
    marker_color = "black",
    base_size = 12
) {

  required_cols <- c("time", "value", "stim_nd")
  missing_cols <- setdiff(required_cols, names(df_long))

  if (length(missing_cols) > 0) {
    stop(
      "df_long is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  df_plot <- df_long

  if (!is.null(channel_filter)) {
    if (!"channel" %in% names(df_plot)) {
      stop("`channel_filter` provided but no `channel` column found.", call. = FALSE)
    }

    df_plot <- df_plot |>
      dplyr::filter(.data$channel == channel_filter)

    if (nrow(df_plot) == 0) {
      stop("No data left after applying `channel_filter`.", call. = FALSE)
    }
  }

  if (isTRUE(compare_raw) && !"value_raw" %in% names(df_plot)) {
    stop("`compare_raw = TRUE` requires a `value_raw` column.", call. = FALSE)
  }

  if (is.null(stim_levels)) {
    stim_unique <- unique(df_plot$stim_nd)

    if (is.numeric(stim_unique)) {
      stim_levels <- sort(stim_unique, decreasing = FALSE)
    } else {
      stim_levels <- sort(as.character(stim_unique), decreasing = FALSE)
    }
  }

  stim_levels_chr <- as.character(stim_levels)

  if (isTRUE(compare_raw)) {
    df_summary <- df_plot |>
      dplyr::group_by(.data$time, .data$stim_nd) |>
      dplyr::summarise(
        Raw = mean(.data$value_raw, na.rm = TRUE),
        Smoothed = mean(.data$value, na.rm = TRUE),
        .groups = "drop"
      ) |>
      tidyr::pivot_longer(
        cols = c("Raw", "Smoothed"),
        names_to = "signal_type",
        values_to = "signal"
      ) |>
      dplyr::mutate(
        stim_nd = factor(as.character(.data$stim_nd), levels = stim_levels_chr)
      )
  } else {
    df_summary <- df_plot |>
      dplyr::group_by(.data$time, .data$stim_nd) |>
      dplyr::summarise(
        signal = mean(.data$value, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        signal_type = "Smoothed",
        stim_nd = factor(as.character(.data$stim_nd), levels = stim_levels_chr)
      )
  }

  marker_lines <- NULL

  if (isTRUE(overlay_markers)) {
    marker_source <- if (isTRUE(compare_raw)) {
      df_summary |>
        dplyr::filter(.data$signal_type == "Smoothed")
    } else {
      df_summary
    }

    get_marker <- function(data, window, label, extreme_fn) {
      data |>
        dplyr::filter(.data$time >= window[1], .data$time <= window[2]) |>
        dplyr::group_by(.data$stim_nd) |>
        extreme_fn(order_by = .data$signal, n = 1, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::transmute(
          stim_nd = .data$stim_nd,
          marker_type = label,
          time_ms = .data$time * 1000
        )
    }

    marker_lines <- dplyr::bind_rows(
      get_marker(marker_source, a_window, "A-wave trough", dplyr::slice_min),
      get_marker(marker_source, b_window, "B-wave peak", dplyr::slice_max),
      get_marker(marker_source, d_window, "D-wave peak", dplyr::slice_max)
    ) |>
      dplyr::mutate(
        marker_type = factor(
          .data$marker_type,
          levels = c("A-wave trough", "B-wave peak", "D-wave peak")
        )
      )
  }

  p <- ggplot2::ggplot(
    df_summary,
    ggplot2::aes(x = .data$time * 1000, y = .data$signal)
  )

  if (isTRUE(compare_raw)) {
    p <- p +
      ggplot2::geom_line(
        ggplot2::aes(linetype = .data$signal_type),
        linewidth = 0.7,
        color = "black"
      ) +
      ggplot2::scale_linetype_manual(
        name = "Signal",
        values = c("Raw" = "solid", "Smoothed" = "longdash")
      )
  } else {
    p <- p +
      ggplot2::geom_line(
        linewidth = 0.8,
        color = "red"
      )
  }

  if (isTRUE(overlay_markers) && !is.null(marker_lines)) {
    p <- p +
      ggplot2::geom_vline(
        data = marker_lines,
        ggplot2::aes(
          xintercept = .data$time_ms,
          linetype = .data$marker_type
        ),
        color = marker_color,
        linewidth = 0.7,
        alpha = 0.95,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_linetype_manual(
        name = if (isTRUE(compare_raw)) "Line type" else "Reference lines",
        values = c(
          "A-wave trough" = "dotted",
          "B-wave peak" = "longdash",
          "D-wave peak" = "dotdash",
          "Raw" = "solid",
          "Smoothed" = "longdash"
        ),
        breaks = if (isTRUE(compare_raw)) {
          c("Raw", "Smoothed", "A-wave trough", "B-wave peak", "D-wave peak")
        } else {
          c("A-wave trough", "B-wave peak", "D-wave peak")
        }
      )
  }

  p +
    ggplot2::facet_wrap(
      ~stim_nd,
      scales = if (isTRUE(free_y)) "free_y" else "fixed"
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 6),
      expand = ggplot2::expansion(mult = c(0.01, 0.03))
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = 5),
      expand = ggplot2::expansion(mult = c(0.03, 0.08))
    ) +
    ggplot2::labs(
      title = "Mean ERG Waveform by ND Level",
      x = "Time (ms)",
      y = expression("Response (" * mu * "V)")
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}
