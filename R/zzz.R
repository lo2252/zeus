#' Register built-in ZEUS protocols
#'
#' Registers the built-in spectral (`C0`) and white-light (`C1`) protocol
#' definitions when the package is loaded.
#'
#' @keywords internal
#' @noRd
.onLoad <- function(libname, pkgname) {
  register_zeus_protocol(
    zeus_protocol(
      id = "C0",
      description = "Spectral light protocol",
      builder = function(...) protocol_table_C0()
    )
  )

  register_zeus_protocol(
    zeus_protocol(
      id = "C1",
      description = "White-light protocol",
      builder = function(...) protocol_table_C1()
    )
  )
}
