# ZEUS Protocol Object ----------------------------------------------------

#' Create a ZEUS protocol object
#'
#' @param id Character scalar protocol identifier.
#' @param description Character scalar protocol description.
#' @param builder Function returning a protocol table.
#'
#' @return An object of class `"zeus_protocol"`.
#' @keywords internal
zeus_protocol <- function(id, description, builder) {
  if (!is.character(id) || length(id) != 1L || !nzchar(id)) {
    stop("'id' must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.character(description) || length(description) != 1L) {
    stop("'description' must be a single character string.", call. = FALSE)
  }

  if (!is.function(builder)) {
    stop("'builder' must be a function.", call. = FALSE)
  }

  structure(
    list(
      id = id,
      description = description,
      builder = builder
    ),
    class = "zeus_protocol"
  )
}


# Protocol Registry -------------------------------------------------------

.zeus_protocol_registry <- new.env(parent = emptyenv())


# Register Protocol -------------------------------------------------------

#' Register a ZEUS protocol
#'
#' @param protocol A `"zeus_protocol"` object.
#'
#' @return Invisibly returns the protocol.
#' @keywords internal
register_zeus_protocol <- function(protocol) {
  if (!inherits(protocol, "zeus_protocol")) {
    stop("'protocol' must inherit from class 'zeus_protocol'.", call. = FALSE)
  }

  assign(protocol$id, protocol, envir = .zeus_protocol_registry)
  invisible(protocol)
}


# Get Registered Protocol -------------------------------------------------

#' Get a registered ZEUS protocol
#'
#' @param id Character scalar protocol identifier.
#'
#' @return A `"zeus_protocol"` object.
#' @export
get_zeus_protocol <- function(id) {
  if (!is.character(id) || length(id) != 1L || !nzchar(id)) {
    stop("'id' must be a single non-empty character string.", call. = FALSE)
  }

  if (!exists(id, envir = .zeus_protocol_registry, inherits = FALSE)) {
    stop("Protocol '", id, "' is not registered.", call. = FALSE)
  }

  get(id, envir = .zeus_protocol_registry, inherits = FALSE)
}


# List Protocols ----------------------------------------------------------

#' List registered ZEUS protocols
#'
#' @return Character vector of registered protocol identifiers.
#' @export
list_zeus_protocols <- function() {
  ls(envir = .zeus_protocol_registry, all.names = FALSE)
}


# C0 Protocol -------------------------------------------------------------

#' Build the ZEUS C0 spectral protocol table
#'
#' Creates the 70-condition spectral protocol table used by the C0 protocol.
#' The design consists of 10 wavelength blocks with 7 neutral-density levels
#' per block.
#'
#' @return Tibble with 70 rows and protocol metadata.
#' @export
protocol_table_C0 <- function() {
  make_block <- function(wavelength, stim_nd) {
    tibble::tibble(
      wavelength = as.character(wavelength),
      stim_nd = as.numeric(stim_nd)
    )
  }

  dplyr::bind_rows(
    make_block(650, c(4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0)),
    make_block(570, c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0)),
    make_block(490, c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0)),
    make_block(410, c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0)),
    make_block(330, c(4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0)),
    make_block(650, c(4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5)),
    make_block(610, c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0)),
    make_block(530, c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0)),
    make_block(450, c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0)),
    make_block(370, c(5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0))
  ) |>
    dplyr::mutate(
      protocol_id = "C0",
      protocol_variant = "p4",
      stim_type = "flash",
      stim_index = dplyr::row_number(),
      block_index = rep(seq_len(10L), each = 7L),
      within_block_index = rep(seq_len(7L), times = 10L),
      .before = 1
    )
}


# C1 Protocol -------------------------------------------------------------

#' Build the ZEUS C1 white-light protocol table
#'
#' Creates the 70-condition white-light protocol table used by the C1 protocol.
#' The design consists of 10 runs with 7 neutral-density levels per run.
#'
#' @return Tibble with 70 rows and protocol metadata.
#' @export
protocol_table_C1 <- function() {
  nd_values <- c(6.0, 5.5, 5.0, 4.5, 4.0, 3.5, 3.0)

  tibble::tibble(
    protocol_id = "C1",
    protocol_variant = "white",
    stim_type = "flash",
    stim_index = seq_len(70L),
    block_index = rep(seq_len(10L), each = 7L),
    within_block_index = rep(seq_len(7L), times = 10L),
    wavelength = rep("White", 70L),
    stim_nd = rep(nd_values, times = 10L)
  )
}


# Protocol Table Dispatcher -----------------------------------------------

#' Build a ZEUS protocol table by protocol identifier
#'
#' @param protocol Protocol identifier. One of `"C0"` or `"C1"`.
#'
#' @return A protocol table tibble.
#' @keywords internal
make_protocol_table <- function(protocol = c("C0", "C1")) {
  protocol <- match.arg(protocol)

  switch(
    protocol,
    C0 = protocol_table_C0(),
    C1 = protocol_table_C1()
  )
}


# Sweep Order Expansion ---------------------------------------------------

#' Expand a protocol table to raw sweep order
#'
#' Replicates each stimulus row according to the number of technical replicates
#' per stimulus and adds sweep-order metadata.
#'
#' @param protocol_tbl A protocol table tibble.
#' @param repeats_per_stim Integer number of technical replicates per stimulus.
#'
#' @return Tibble with one row per raw sweep.
#' @export
expand_protocol_repeats <- function(protocol_tbl, repeats_per_stim = 4L) {
  if (!is.data.frame(protocol_tbl)) {
    stop("'protocol_tbl' must be a data frame.", call. = FALSE)
  }

  required <- c(
    "protocol_id",
    "protocol_variant",
    "stim_type",
    "stim_index",
    "block_index",
    "within_block_index",
    "wavelength",
    "stim_nd"
  )
  missing_cols <- setdiff(required, names(protocol_tbl))

  if (length(missing_cols) > 0L) {
    stop(
      "'protocol_tbl' must contain columns: ",
      paste(required, collapse = ", "),
      call. = FALSE
    )
  }

  repeats_per_stim <- as.integer(repeats_per_stim)

  if (!is.finite(repeats_per_stim) || repeats_per_stim < 1L) {
    stop("'repeats_per_stim' must be a positive integer.", call. = FALSE)
  }

  protocol_tbl |>
    tibble::as_tibble() |>
    tidyr::uncount(weights = repeats_per_stim, .id = "tech_rep") |>
    dplyr::mutate(
      sweep_in_protocol = dplyr::row_number(),
      .before = 1
    )
}


# Label Traces with Protocol Metadata -------------------------------------

#' Label ERG traces with stimulus protocol metadata
#'
#' Maps each raw sweep to a protocol-defined stimulus condition and joins the
#' corresponding protocol metadata onto a long ERG trace table.
#'
#' @param erg_df Long ERG data containing at least `sweep`, `time`, and `value`.
#' @param protocol Protocol identifier. One of `"C0"` or `"C1"`.
#' @param repeats_per_stim Integer number of technical replicates per stimulus.
#'
#' @return Long ERG data with protocol-derived stimulus metadata added,
#'   including `stim_index`, `stim_label`, `protocol_id`, `protocol_variant`,
#'   `stim_type`, `block_index`, `within_block_index`, `wavelength`, `stim_nd`,
#'   `tech_rep`, and `sweep_in_protocol`.
#' @keywords internal
label_stimulus_protocol <- function(erg_df,
                                    protocol = c("C0", "C1"),
                                    repeats_per_stim = 4L) {
  protocol <- match.arg(protocol)
  repeats_per_stim <- as.integer(repeats_per_stim)

  needed <- c("sweep", "time", "value")
  missing_cols <- setdiff(needed, names(erg_df))

  if (length(missing_cols) > 0L) {
    stop(
      "`erg_df` must contain columns: ",
      paste(needed, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.finite(repeats_per_stim) || repeats_per_stim < 1L) {
    stop("'repeats_per_stim' must be a positive integer.", call. = FALSE)
  }

  protocol_tbl <- make_protocol_table(protocol = protocol)

  sweep_tbl <- expand_protocol_repeats(
    protocol_tbl = protocol_tbl,
    repeats_per_stim = repeats_per_stim
  ) |>
    dplyr::mutate(
      sweep = .data$sweep_in_protocol
    )

  erg_df |>
    dplyr::left_join(
      sweep_tbl,
      by = "sweep"
    ) |>
    dplyr::mutate(
      stim_label = purrr::pmap_chr(
        list(.data$wavelength, .data$stim_nd, .data$stim_index, .data$protocol_id),
        make_stimresp_label
      )
    )
}