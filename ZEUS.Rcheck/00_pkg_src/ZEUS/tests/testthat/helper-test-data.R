.zeus_test_cache <- new.env(parent = emptyenv())

zeus_test_abf_path <- function(filename) {
  path <- zeus_extdata(filename)

  if (!nzchar(path) || !file.exists(path)) {
    testthat::skip(paste(filename, "not found in inst/extdata"))
  }

  path
}

zeus_test_obj <- function(protocol) {
  protocol <- match.arg(protocol, c("C0", "C1"))
  key <- paste0("obj_", protocol)

  if (!exists(key, envir = .zeus_test_cache, inherits = FALSE)) {
    filename <- if (identical(protocol, "C0")) "26225005.abf" else "26225004.abf"
    path <- zeus_test_abf_path(filename)

    assign(
      key,
      zeus_read_abf(path, protocol = protocol),
      envir = .zeus_test_cache
    )
  }

  get(key, envir = .zeus_test_cache, inherits = FALSE)
}
