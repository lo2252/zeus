if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Package 'pkgload' is required to run the development app.", call. = FALSE)
}

options(shiny.maxRequestSize = 100 * 1024^2)

pkgload::load_all(path = "..", export_all = FALSE)
ZEUS::run_zeus_app()
