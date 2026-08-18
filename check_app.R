required <- c("shiny", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  stop(
    "Install missing packages first: ",
    paste(missing, collapse = ", ")
  )
}

project_directory <- normalizePath(getwd())
app_directory <- file.path(project_directory, "app")

parse(file.path(app_directory, "app.R"))

dat <- readRDS(
  file.path(
    app_directory,
    "data",
    "hmp_v35_five_habitats_40_each.rds"
  )
)
stopifnot(
  is.matrix(dat$counts),
  is.matrix(dat$phylum_counts),
  is.matrix(dat$family_counts),
  ncol(dat$counts) == nrow(dat$metadata),
  identical(colnames(dat$phylum_counts), rownames(dat$metadata)),
  identical(colnames(dat$family_counts), rownames(dat$metadata)),
  identical(colnames(dat$counts), rownames(dat$metadata)),
  all(c("HMP_BODY_SUBSITE", "HABITAT") %in% names(dat$metadata)),
  length(unique(dat$metadata$HABITAT)) == 5
)

# Evaluate app.R with its own directory as the temporary working
# directory, matching the behavior of shiny::runApp("app").
app <- source(
  file.path(app_directory, "app.R"),
  local = new.env(parent = globalenv()),
  chdir = TRUE
)$value
stopifnot(inherits(app, "shiny.appobj"))

message("App syntax, data alignment, and Shiny object checks passed.")
