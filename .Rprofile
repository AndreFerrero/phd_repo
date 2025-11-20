# Load core packages automatically
core_packages <- c("here", "tidyverse", "evd", "copula", "rmarkdown")

for (pkg in core_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

message("Core packages loaded: ", paste(core_packages, collapse = ", "))
