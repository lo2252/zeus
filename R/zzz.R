#' Default imported protocol
#' 
#' Imports the pre-defines spectral (C0) and white-light (C1) protocols. 
.onLoad <- function(libname, pkgname) {
  register_zeus_protocol(
    new_zeus_protocol(
      id = "C0",
      description = "Spectral Light Protocol",
      builder = function(...) protocol_table_C0()
    )
  )
  
  register_zeus_protocol(
    new_zeus_protocol(
      id = "C1",
      description = "White-light Protocol",
      builder = function(...) protocol_table_C1()
    )
  )
}