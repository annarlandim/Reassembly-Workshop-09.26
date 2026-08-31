## =============================================================================
## 01_clean_data.R
##
## IN:  data/raw/     (never modified)
## OUT: data/clean/   (corrected, renamed, typed)
##      data/processed/ (joined, analysis-ready)
##
## Run via: source(here::here("utils", "workflow.R"))
## =============================================================================

library(dplyr)
library(stringr)

source(here::here("utils", "functions.R"))


## -----------------------------------------------------------------------------
## Inputs
## -----------------------------------------------------------------------------
## Note every path is built with here::here(). Nothing here depends on your
## working directory, your username, or your operating system.

paths <- list(
  herbivory = here::here("data", "raw", "Censo_Herbivory.csv"),
  censo     = here::here("data", "raw", "Censo_clean_27042026.csv"),
  traits    = here::here("data", "raw", "FOTODISK_LSA_master.xlsx"),
  landcover = here::here("data", "raw", "Reassembly_Landcover_FNewell_mod.csv")
)

require_files(unlist(paths))


## -----------------------------------------------------------------------------
## Read
## -----------------------------------------------------------------------------
## The traits workbook has a merged section label in row 1 above the real
## header ("leaf_penetrometer" spanning several columns), so skip it explicitly
## rather than letting readxl guess.

herb_raw   <- read_reassembly_csv(paths$herbivory, delim = ",")
censo_raw  <- read_reassembly_csv(paths$censo,     delim = ",")
land_raw   <- read_reassembly_csv(paths$landcover, delim = ";")
traits_raw <- readxl::read_excel(
  paths$traits, sheet = "LSA_master", skip = 1,
  na = c("", "N/A", "NA")   # cells literally containing "N/A" text, not blank —
  # readxl's default na="" misses these and silently
  # imports the whole column as character
)


## -----------------------------------------------------------------------------
## Clean: herbivory
## -----------------------------------------------------------------------------

herb_clean <- herb_raw |>
  rename(
    unique_code        = Unique_code,
    species_raw        = Species,
    plot_id            = `No.Plot`,
    n_leaves           = `#Leaves (H)`,
    n_leaves_damaged   = `# Leaves with damage`,
    total_damage       = `Total damage`
  ) |>
  mutate(species_raw = fix_species_typo(species_raw)) |>
  fix_xp217_rotation(plot_col = "plot_id") |>
  rename(block = Block, subblock = Subblock) |>   # must come AFTER fix_xp217_rotation(),
  # which looks for "Block"/"Subblock" by name
  add_code_parts() |>
  distinct(unique_code, .keep_all = TRUE) |>   # 25 exact duplicate rows in raw
  mutate(leaf_count_flag = !is.na(n_leaves) & !is.na(n_leaves_damaged) &
           n_leaves_damaged > n_leaves)

if (any(herb_clean$leaf_count_flag, na.rm = TRUE)) {
  bad <- herb_clean |> filter(leaf_count_flag)
  warning(
    "More damaged leaves than total leaves recorded for: ",
    paste(bad$unique_code, collapse = ", "),
    " — n_leaves_damaged > n_leaves, which is impossible. Kept in herb_clean ",
    "and flagged (not silently corrected — unclear which number is wrong). ",
    "Excluded downstream wherever incidence is computed (see dat_inc in ",
    "03_glmm_herbivory.R)."
  )
}
save_data(herb_clean, "herbivory_clean.csv", stage = "clean")

## -----------------------------------------------------------------------------
## Clean: censo
## -----------------------------------------------------------------------------
## CAUTION: the column named `Survival` is a MORTALITY flag. Survival == 1
## corresponds 1:1 to Date_mortality being filled in. Renaming it here means
## nobody downstream has to remember the flip.

censo_clean <- censo_raw |>
  rename(
    unique_code    = Unique_code,
    plot_id        = Plot_ID,
    block          = Block,
    subblock       = Subblock,
    date_mortality = Date_mortality,
    died           = Survival           # <- renamed on purpose, see above
  ) |>
  mutate(
    date_planted = as.Date(Date_treelet_planted, format = "%B %d %Y"),
    date_control = as.Date(Date_tree_Control, format = "%B %d %Y"),
    days_at_risk = as.numeric(date_control - date_planted)
  ) |>
  add_code_parts() |>
  select(unique_code, species_acronym, mother_tree,
         plot_id, block, subblock, died, date_mortality, days_at_risk)

save_data(censo_clean, "censo_clean.csv", stage = "clean")


## -----------------------------------------------------------------------------
## Clean: traits (leaf level -> plant level)
## -----------------------------------------------------------------------------

traits_clean <- traits_raw |>
  rename(
    toughness  = mean_pen,
    chl_ug_cm2 = `µg/cm2`
  ) |>
  group_by(unique_code, species_acronym) |>
  summarise(
    toughness      = mean(toughness, na.rm = TRUE),
    chl_ug_cm2     = mean(chl_ug_cm2, na.rm = TRUE),
    fresh_weight_g = mean(fresh_weight_g, na.rm = TRUE),
    dry_weight_g   = mean(dry_weight_g, na.rm = TRUE),
    .groups = "drop"
  )

save_data(traits_clean, "traits_clean.csv", stage = "clean")


## -----------------------------------------------------------------------------
## Clean: plot-level landscape covariates
## -----------------------------------------------------------------------------
## FLAG, unresolved: XP201-PR and XP202-PR have Age_2026 == 0 but
## Regeneration_year == 2024 (implying age 2), and both are marked
## Need_check == "Correct". Every other plot is internally consistent.
## Confirm with Felicity / CM before using Age_2026 as a predictor.

land_clean <- land_raw |>
  rename(
    plot_id         = Plot_ID,
    forest_1km      = Forest_1km,
    distance_forest = Distance_forest,
    distance_edge   = Distance_edge,
    age_2026        = Age_2026,
    elevation       = Elevation,
    treatment       = Treatment
  ) |>
  mutate(
    regen_year_known = !is.na(Regeneration_year) & Regeneration_year != 0,
    implied_age = ifelse(regen_year_known,
                         2026 - as.numeric(Regeneration_year),
                         NA_real_),
    ## Flag A: a real, known regeneration year that doesn't match the
    ## recorded age (e.g. XP201-PR, XP202-PR).
    age_flag_mismatch = regen_year_known & age_2026 != implied_age,
    ## Flag B: no regeneration year recorded (0/unknown) AND the plot still
    ## carries the old-growth ceiling age (50), despite not being labelled
    ## old-growth (e.g. XP220-PLR). NOT assumed benign — surfaced for
    ## Felicity to confirm whether this plot's age/treatment are correct.
    age_flag_unverified_og_age = !regen_year_known &
      age_2026 == 50 &
      treatment != "Old-Growth",
    age_flag = age_flag_mismatch | age_flag_unverified_og_age
  ) |>
  select(plot_id, treatment, age_2026, implied_age,
         age_flag_mismatch, age_flag_unverified_og_age, age_flag,
         forest_1km, distance_forest, distance_edge, elevation)

if (any(land_clean$age_flag_mismatch, na.rm = TRUE)) {
  warning(
    "Age/regeneration-year mismatch in plot(s): ",
    paste(land_clean$plot_id[land_clean$age_flag_mismatch], collapse = ", "),
    " — recorded age doesn't match what the regeneration year implies. ",
    "Confirm with Felicity/CM."
  )
}
if (any(land_clean$age_flag_unverified_og_age, na.rm = TRUE)) {
  warning(
    "Unverified old-growth-ceiling age in non-old-growth plot(s): ",
    paste(land_clean$plot_id[land_clean$age_flag_unverified_og_age], collapse = ", "),
    " — age = 50 with no regeneration year recorded, but not labelled ",
    "Old-Growth. Do not assume this is a benign placeholder — confirm ",
    "with Felicity whether the age or the treatment label is correct."
  )
}

save_data(land_clean, "landcover_clean.csv", stage = "clean")


## -----------------------------------------------------------------------------
## Process: join into the analysis table
## -----------------------------------------------------------------------------

dat <- herb_clean |>
  select(unique_code, species_acronym, mother_tree, plot_id,
         block, subblock, n_leaves, n_leaves_damaged, total_damage,
         leaf_count_flag) |>
  left_join(censo_clean |> select(unique_code, died, days_at_risk), by = "unique_code") |>
  left_join(traits_clean, by = c("unique_code", "species_acronym")) |>
  left_join(land_clean,   by = "plot_id") |>
  mutate(
    n_leaves_undamaged = pmax(n_leaves - n_leaves_damaged, 0),
    severity_prop      = total_damage / 100
  ) |>
  mutate(
    severity_prop_sv = sv_transform(severity_prop, n = dplyr::n())
  )

save_data(dat, "trex_analysis.csv", stage = "processed")

message("  trait data available for ",
        length(unique(dat$species_acronym[!is.na(dat$toughness)])),
        " of ", length(unique(dat$species_acronym)), " species")
