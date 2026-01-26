#' Importing a file in inst/extdata for testing and demonstration
#' @keywords internal


zeus_extdata <- function(...){
  p <- system.file("extdata", ..., package = "ZEUS")
  if (nzchar(p)) return(p)

  root <- tryCatch(pkgload::pkg_path(), error = function(e) NULL)
  if (!is.null(root)) {
    p2 <- file.path(root, "inst", "extdata", ...)
    if (file.exists(p2)) return(p2)
  }
  "" # returns an empty string
}
