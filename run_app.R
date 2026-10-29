# run_app.R

# 1. Point explicitly to the local R-Portable library if it exists
app_dir <- getwd()
portable_lib <- file.path(app_dir, "R-Portable", "library")

if (dir.exists(portable_lib)) {
  .libPaths(c(portable_lib, .libPaths()))
}

# 2. Load required packages (rstudioapi is omitted)
library(shiny)
library(GLM3r)
library(glmtools)
library(ggplot2)
library(scales)
library(shinyjs)

# 3. Launch Shiny app
options(shiny.launch.browser = TRUE)
shiny::runApp(appDir = app_dir, launch.browser = TRUE)