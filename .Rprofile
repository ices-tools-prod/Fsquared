
# LOAD user .Rprofile if it exists
if(file.exists(file.path(Sys.getenv("HOME"), ".Rprofile")))
  source(file.path(Sys.getenv("HOME"), ".Rprofile"))

# USE progressr
if (interactive() && requireNamespace("progressr", quietly = TRUE)) {
  ## Enable global progression updates
  if (getRversion() >= "4.0.0") progressr::handlers(global = TRUE)
 
  ## In RStudio Console, or not?
  if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM"))) {
    options(progressr.handlers = progressr::handler_rstudio)
  } else {
    options(progressr.handlers = progressr::handler_progress)
  }
}
