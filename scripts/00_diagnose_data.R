## =============================================================================
## 02_diagnose_data.R  (save locally as 00_diagnose_data.R if you want it to
##                      sort before 01 — see note in the project README)
##
## IN:  data/raw/   (deliberately — see note below)
## OUT: console diagnostics only. Writes nothing, feeds nothing downstream.
##
## WHY THIS READS RAW, NOT CLEAN:
## Every problem below is already FIXED by 01_clean_data.R. This script exists
## to show HOW each one was found — using ordinary, boring tools (str, table,
## head) rather than anything clever. It's a teaching/audit script, not part
## of the pipeline: nothing downstream depends on it, and deleting it breaks
## nothing.
##
## The lesson: none of these problems announce themselves. Every one produces
## code that runs without error and gives you a wrong answer. str() and
## table() are usually enough to catch them.
##
## Run via: source(here::here("utils", "workflow.R"))
## =============================================================================

library(dplyr)   # only used for require_files() below; everything else is base R

source(here::here("utils", "functions.R"))

paths <- list(
  herbivory = here::here("data", "raw", "Censo_Herbivory.csv"),
  censo     = here::here("data", "raw", "Censo_clean_27042026.csv"),
  traits    = here::here("data", "raw", "FOTODISK_LSA_master.xlsx")
)
require_files(unlist(paths))

herb_raw  <- read_reassembly_csv(paths$herbivory, delim = ",")
censo_raw <- read_reassembly_csv(paths$censo,     delim = ",")


## -----------------------------------------------------------------------------
## CHECK 0: column types — the cheapest, most-skipped check there is
## -----------------------------------------------------------------------------
## A column that should be numeric but isn't means stray text somewhere.
## str() shows you this in one line, before you've done anything else.

cat("\n--- CHECK 0: column types in the traits workbook ---\n")

## Naive read: readxl's default na = "" only treats truly blank cells as
## missing, so a cell containing the literal text "N/A" is kept as text.
traits_naive <- readxl::read_excel(paths$traits, sheet = "LSA_master", skip = 1)

str(traits_naive)
## traits_naive$mean_pen chr [1:642] "..." ...   <- should say "num", not "chr". Something's wrong.

## Try converting it. If R can't, it tells you so directly:
as.numeric(traits_naive$mean_pen)
## Warning message: NAs introduced by coercion

## table() on just the values that failed shows you exactly what's there:
table(traits_naive$mean_pen[is.na(as.numeric(traits_naive$mean_pen))])
##  N/A
##    9

## The fix: tell readxl "N/A" text also means missing.
traits_fixed <- readxl::read_excel(paths$traits, sheet = "LSA_master", skip = 1,
                                   na = c("", "N/A", "NA"))
str(traits_fixed$mean_pen)
## num [1:642] ...   <- fixed. This is what 01_clean_data.R does.


## -----------------------------------------------------------------------------
## CHECK 1: duplicated identifiers
## -----------------------------------------------------------------------------
## If your unique ID isn't unique, every merge downstream silently duplicates
## rows and your sample size quietly inflates.

cat("\n--- CHECK 1: duplicated Unique_code in herbivory ---\n")

nrow(herb_raw)
length(unique(herb_raw$Unique_code))
sum(duplicated(herb_raw$Unique_code))

dup_codes <- herb_raw$Unique_code[duplicated(herb_raw$Unique_code)]
herb_raw[herb_raw$Unique_code %in% dup_codes,
         c("Unique_code", "No.Plot", "Block", "Subblock", "Location")]

## STILL TO CHECK: are these exact duplicate rows, or the same plant measured
## twice with different values? 01_clean_data.R keeps the first via distinct(),
## which is only defensible if they're exact duplicates:
sum(duplicated(herb_raw))   # fully identical rows in the whole file


## -----------------------------------------------------------------------------
## CHECK 2: free-text categories — table() them, always
## -----------------------------------------------------------------------------
## A typo creates a phantom extra species. table() shows every distinct value
## with its count, side by side — a typo usually stands out immediately.

cat("\n--- CHECK 2: species names ---\n")

table(herb_raw$Species)
## look for near-duplicate names, e.g. "Bahuinia pichinchensis" (5) sitting
## right next to "Bauhinia pichinchensis" (335) — that's a typo, not two
## species.

length(unique(herb_raw$Species))   # 17, but only 16 species were planted


## -----------------------------------------------------------------------------
## CHECK 3: does a column mean the same thing in every row?
## -----------------------------------------------------------------------------
## table() of a single column, on its own, can reveal an inconsistency you'd
## never think to test for directly.

cat("\n--- CHECK 3: Block/Subblock encoding ---\n")

table(herb_raw$Block)
## 1, 2, 3, 4, 5 are expected (the T-REX block number) -- but "A", "B", "C"
## show up too. Those are subblock letters. Something's swapped, somewhere.

## head() on one normal plot and the odd one, side by side:
head(herb_raw[herb_raw$`No.Plot` == "XP201-PR",
              c("Block", "Subblock", "Location", "Tube")])
head(herb_raw[herb_raw$`No.Plot` == "XP217-PLR",
              c("Block", "Subblock", "Location", "Tube")])
## XP217-PLR has Block and Subblock (and Location and Tube) reversed relative
## to every other plot. Confirmed with table() below: it's the ONLY plot where
## Block holds letters.

table(herb_raw$`No.Plot`[herb_raw$Block %in% c("A", "B", "C")])


## -----------------------------------------------------------------------------
## CHECK 4: does the data match what people TOLD you about it?
## -----------------------------------------------------------------------------
## Julie's message said fieldwork covered blocks 1-4 and skipped block 5.
## table() answers this directly.

cat("\n--- CHECK 4: does block 5 exist? ---\n")

table(herb_raw$`No.Plot`[herb_raw$Block == "5"])
## block "5" shows up in MOST plots, not none. Neither Julie's account nor
## the raw data is obviously wrong on its own -- this is a question to ask
## (of Julie / CM), not a discrepancy to silently resolve.


## -----------------------------------------------------------------------------
## CHECK 5: does a column mean what its NAME says?
## -----------------------------------------------------------------------------
## The most dangerous item in this script. `Survival` is a MORTALITY flag.
## glmer(Survival ~ ..., family = binomial()) fits perfectly, warns about
## nothing, and returns every effect backwards.

cat("\n--- CHECK 5: is `Survival` really survival? ---\n")

table(Survival = censo_raw$Survival,
      has_mortality_date = !is.na(censo_raw$Date_mortality))
## Survival == 1 lines up almost exactly with "has a mortality date". It
## means DIED. Renamed to `died` in 01_clean_data.R.


cat("\n--- diagnostics complete ---\n")