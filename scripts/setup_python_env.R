# One-time setup: creates the project's reticulate virtualenv and installs
# pinned Python dependencies (python/requirements.txt). Re-run safely; it
# reuses the existing venv if present.

if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate", repos = "https://cloud.r-project.org")
}

venv_name <- "proteoformtracker-py"

if (!reticulate::virtualenv_exists(venv_name)) {
  reticulate::virtualenv_create(venv_name)
}

reticulate::virtualenv_install(
  venv_name,
  requirements = "python/requirements.txt"
)

reticulate::use_virtualenv(venv_name, required = TRUE)
reticulate::py_config()
message("Python environment ready: ", venv_name)
