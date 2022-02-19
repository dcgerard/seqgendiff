# Environment
message_env <- new.env(parent = emptyenv())
message_env$poisthinwarning <- FALSE

message_fun <- function(x) {
  if (x == "poisthin" & !message_env$poisthinwarning) {
    message(paste0(
      "poisthin() is now defunct.\n",
      "Please use select_counts() and thin_2group() instead.\n",
      "This message is displayed once per R session"))
    message_env$poisthinwarning <- TRUE
  }
}
