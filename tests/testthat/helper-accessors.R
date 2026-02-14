get_fun <- function(name) {
  # Prefer exported function if present
  if (exists(name,
             where = asNamespace("mpmaggregate"),
             inherits = FALSE)) {
    return(get(name, envir = asNamespace("mpmaggregate")))
  }
  # Else fall back to internal
  getFromNamespace(name, "mpmaggregate")
}
