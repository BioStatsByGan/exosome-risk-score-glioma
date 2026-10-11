options(repos = c(CRAN = "https://cloud.r-project.org"), warn = -1)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("MCPcounter", quietly = TRUE)) {
  remotes::install_github("ebecht/MCPcounter", ref = "master", subdir = "Source",
                          upgrade = "never", dependencies = TRUE)
}
cat("MCPcounter installed:", requireNamespace("MCPcounter", quietly = TRUE), "\n")
