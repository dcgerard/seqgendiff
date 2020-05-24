# Environment
message_env <- new.env(parent = emptyenv())
message_env$poisthinwarning <- FALSE
message_env$optmatchwarning <- FALSE

message_fun <- function(x) {
  if (x == "optmatch" & !message_env$optmatchwarning) {
    message(paste0("Note that optmatch uses a strange non-standard license:\n",
                   "https://cran.r-project.org/package=optmatch/LICENSE\n\n",
                   "If this doesn't work for you, try permute_method = \"hungarian\"\n\n",
                   "This message is displayed once per R session."))
    message_env$optmatchwarning <- TRUE
  } else if (x == "poisthin" & !message_env$poisthinwarning) {
    message(paste0(
      "poisthin() is now defunct.\n",
      "Please use select_counts() and thin_2group() instead.\n",
      "This message is displayed once per R session"))
    message_env$poisthinwarning <- TRUE
  }
}
