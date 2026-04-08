
# ZEUS Protocol Object Function -------------------------------------------


#'  ZEUS Protocol Object
#'
#'  @param id Character scalar protocol id.
#'  @param desciption Scalar description.
#'  @param builder Function returning a protocol table.
#'
#'  @return An object of class "zeus_protocol"
#'  @keywords internal
zeus_protocol <- function(id, description, builder) {
  if (!is.character(id) || length(id) != 1L || !nzchar(id)) {
    stop("'id' must be a single non-empty character string.", call. = FALSE )
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


# Environment Developer
.zeus_protocol_registry <- new.env(parent = emptyenv())

# Register Protocol -------------------------------------------------------


#' Register a ZEUS protocol
#'
#' @param protocol A "zeus_protocol" object.
#'
#' @return Invisily returns the protocol.
#' @keywords internal
register_zeus_protocol <- function(protocol) {
  if (!inherits(protocol, "zeus_protocol")) {
    stop("'protocol' must be from class 'zeus_protocol'.", call. = FALSE)
  }
  assign(protocol$id, protocol, envir = .zeus_protocol_registry)
  invisible(protocol)
}


# Calls Registered Protocol -----------------------------------------------


#' Get a resisted ZEUS protocol
#'
#' @param id Character scalar protocol id.
#'
#' @return A "zeus_protocol" object.
#' @export
get_zeus_protocol <- function(id){
  if (!exists(id, envir = .zeus_protocol_registry, inherits = FALSE)) {
    stop(
      "Protocol '", id, "' is not registered.",
      call. = FALSE
    )
  }

  get(id, envir = .zeus_protocol_registry, inherits = FALSE)
}


# List Protocols ----------------------------------------------------------

#' List registered ZEUS protocols
#'
#' @return Character vector of registered protocol ids.
#' @export
list_zeus_protocols <- function() {
  ls(envir = .zeus_protocol_registry, all.names = FALSE)
}


# C0 - Protocol -----------------------------------------------------------

#' Build the ZEUS C0 spectral protocol table
#'
#' 10 wavelength blocks x 7 ND levels = 70 stimulus conditions.
#'
#' @return Tibble with 70 rows.
#' @export
protocol_table_C0 <- function() {
  make_block <- function(wavelength, stim_nd) {
    tibble::tibble(
      wavelength = as.character(wavelength),
      stim_nd = stim_nd
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
      protocol_variant = "p4", # Based on origin documentation
      stim_type = "flash",
      stim_index = dplyr::row_number(),
      block_index = rep(seq_len(10L), each = 7L),
      within_block_index = rep(seq_len(7L), times = 10L),
      .before = 1
    )
}


# C1 - Protocol -----------------------------------------------------------


#' Build the ZEUS C1 white-light protocol table
#'
#' 10 runs x 7 ND levels = 70 stimulus conditions.
#'
#' @return Tibble with 70 rows.
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


# Sweep Order Protocol ----------------------------------------------------

#' Expand a 70-row protocol table to raw sweep order
#'
#' @param protocol_tbl A 70-row protocol table.
#' @param repeats_per_stim Integer number of technical replicates.
#'
#' @return Tibble with one row per raw sweep.
#' @export
expand_protocol_repeats <- function(protocol_tbl, repeats_per_stim = 4L) {
  if (!is.data.frame(protocol_tbl)) {
    stop("'protocol_tbl' must be a data frame.", call. = FALSE)
  }

  repeats_per_stim <- as.integer(repeats_per_stim)

  if (repeats_per_stim < 1L) {
    stop("'repeats_per_stim' must be >= 1.", call. = FALSE)
  }

  protocol_tbl |>
    tidyr::uncount(weights = repeats_per_stim, .id = "tech_rep") |>
    dplyr::mutate(
      sweep_in_protocol = dplyr::row_number(),
      .before = 1
    )
}






