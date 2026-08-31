## =============================================================================
## utils/functions.R
##
## All reusable functions for this project.
##
## RULE: this file DEFINES things. It must never DO things — no reading data,
## no fitting models, no writing files, no library() side effects you don't
## need. You should be able to source() it at any time, safely, twice.
##
## Sourced by every script in scripts/ and by utils/workflow.R.
## =============================================================================


## -----------------------------------------------------------------------------
## Reading: portability helpers
## -----------------------------------------------------------------------------

#' Read a Reassembly CSV export robustly
#'
#' The exports from the CM Google Sheet and from collaborators are not uniform:
#' some are comma-separated, some semicolon; several have a UTF-8 BOM at the
#' start of the first column name (which silently turns `Plot_ID` into
#' `\ufeffPlot_ID`); several have Windows line endings. On a machine with a
#' German locale, `read.csv2()` vs `read.csv()` also disagree about the decimal
#' mark. All of that is invisible until a join mysteriously returns zero rows.
#'
#' readr handles the BOM and line endings, and we set the delimiter and locale
#' explicitly so the result does not depend on whose computer it runs on.
#'
#' @param path Path to the file, built with here::here().
#' @param delim Field delimiter. Reassembly exports are usually ";".
#' @param decimal_mark "." or ",". State it; never rely on the system locale.
#' @return A tibble.
read_reassembly_csv <- function(path, delim = ";", decimal_mark = ".") {
  stopifnot(file.exists(path))
  readr::read_delim(
    file = path,
    delim = delim, # if it is , or ; or..
    locale = readr::locale(decimal_mark = decimal_mark, encoding = "UTF-8"),
    na = c("", "NA", "#N/A", "NULL"), # many different NAs found in the original files
    trim_ws = TRUE,
    show_col_types = FALSE,
    progress = FALSE
  )
}


#' Stop with a helpful message if required input files are missing
#'
#' Data are not tracked in git (see README). A new user gets a clear list of
#' what to copy out of Nextcloud, instead of a cryptic error 40 lines later.
#'
#' @param paths Named character vector of file paths.
require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(
      "Missing input file(s). Copy them from the Reassembly Nextcloud into data/raw/:\n",
      paste0("  - ", basename(missing), collapse = "\n"),
      call. = FALSE
    )
  }
  invisible(TRUE)
}


## -----------------------------------------------------------------------------
## Cleaning: known, documented corrections to the Reassembly data
##
## Each of these encodes a real problem found in the files. Keeping them as
## named functions (rather than inline edits) means the correction is
## documented, testable, and applied identically everywhere.
## -----------------------------------------------------------------------------

#' Fix the known species-name typo
#'
#' "Bahuinia pichinchensis" (5 rows) vs "Bauhinia pichinchensis" (335 rows).
fix_species_typo <- function(x) {
  dplyr::recode(x, "Bahuinia pichinchensis" = "Bauhinia pichinchensis")
}


#' Fix the column rotation in plot XP217-PLR
#'
#' Every plot except XP217-PLR stores Block as a number and Subblock as a
#' letter. XP217-PLR uses a cyclic rotation of all four position columns:
#'   raw Block -> true Subblock
#'   raw Subblock -> true Location
#'   raw Location -> true Tube
#'   raw Tube -> true Block
#' Verified consistent across all 209 rows of that plot.
#'
#' NOTE: still to be confirmed against Julie's field notes for XP217-PLR.
#'
#' @param df Data frame with columns No.Plot/Block/Subblock/Location/Tube.
fix_xp217_rotation <- function(df,
                               plot_col = "No.Plot",
                               target_plot = "XP217-PLR") {
  needed <- c(plot_col, "Block", "Subblock", "Location", "Tube")
  stopifnot(all(needed %in% names(df)))
  
  rows <- df[[plot_col]] == target_plot
  if (!any(rows, na.rm = TRUE)) {
    warning("fix_xp217_rotation(): no rows matched ", target_plot, " — nothing changed.")
    return(df)
  }
  
  out <- df
  out$Block[rows]    <- df$Tube[rows]
  out$Subblock[rows] <- df$Block[rows]
  out$Location[rows] <- df$Subblock[rows]
  out$Tube[rows]     <- df$Location[rows]
  out
}


#' Split the T-REX unique plant code into its parts
#'
#' Codes have the form SPECIES-MTn-Pnn, e.g. "CH-MT4-P29". Parsing the code is
#' safer than trusting the free-text species column (which has the typo above),
#' and it recovers the mother-tree ID needed for the genotype models.
add_code_parts <- function(df, code_col = "unique_code") {
  stopifnot(code_col %in% names(df))
  df |>
    dplyr::mutate(
      species_acronym = stringr::str_extract(.data[[code_col]], "^[A-Z]+"),
      mother_tree     = stringr::str_extract(.data[[code_col]], "(?<=-MT)\\d+")
    )
}


## -----------------------------------------------------------------------------
## Response variables
## -----------------------------------------------------------------------------

#' Smithson & Verkuilen (2006) transform for beta regression
#'
#' Beta regression requires y strictly inside (0, 1). Herbivory severity has
#' genuine zeros. This nudges 0s and 1s just inside the bounds without
#' distorting the rest of the distribution.
#'
#' @param p Proportion in [0, 1].
#' @param n Sample size used for the transform.
sv_transform <- function(p, n) {
  (p * (n - 1) + 0.5) / n
}


## -----------------------------------------------------------------------------
## Saving: keep output paths consistent and reproducible
## -----------------------------------------------------------------------------

#' Save a figure to output/figures/ with consistent dimensions
save_figure <- function(plot, filename, width = 7, height = 5, dpi = 300) {
  path <- here::here("output", "figures", filename)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = dpi)
  message("  saved figure: ", filename)
  invisible(path)
}


#' Save a model summary or table to output/results/
save_result <- function(object, filename) {
  path <- here::here("output", "results", filename)
  if (grepl("\\.csv$", filename)) {
    readr::write_csv(object, path)
  } else {
    utils::capture.output(print(object), file = path)
  }
  message("  saved result: ", filename)
  invisible(path)
}


#' Save an intermediate dataset to data/clean/ or data/processed/
save_data <- function(df, filename, stage = c("clean", "processed")) {
  stage <- match.arg(stage)
  path <- here::here("data", stage, filename)
  readr::write_csv(df, path)
  message("  saved ", stage, " data: ", filename, " (", nrow(df), " rows)")
  invisible(path)
}