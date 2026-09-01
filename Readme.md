# Reassembly Statistics Workshop — T-REX analyses

Two things to take from this repo:

1.  **Analyses:** GLMMs on real T-REX herbivory and mortality data.
2.  **Structure**: how the project is organised so it runs on anyone's computer.

------------------------------------------------------------------------

## Getting started

1.  **Open the project by double-clicking `reassembly-workshop.Rproj`.** This sets the working directory correctly and activates the package library.

2.  **Get the data.** The raw files are *not* in this repository (see [Data policy](#data-policy)). Copy these four from the Reassembly Nextcloud into `data/raw/`:

    | File | From |
    |------------------------------------|------------------------------------|
    | `Censo_Herbivory.csv` | Julie — herbivory census |
    | `Censo_clean_27042026.csv` | Edith / CM — censo, mortality, size |
    | `FOTODISK_LSA_master.xlsx` | Hanna — leaf traits (LSA-2050) |
    | `Reassembly_Landcover_FNewell_mod.csv` | Felicity — plot landscape covariates |

3.  **Install packages.** Opening the `.Rproj` activates `renv`, which will tell you if anything's missing. Then:

    ``` r
    renv::restore()
    ```

4.  **You can run everything:**

    ``` r
    source(here::here("utils", "workflow.R"))
    ```

    That single line reproduces every figure and table in `output/`. If it fails, something is broken.

------------------------------------------------------------------------

## The one rule

**Here you shouldn't use `setwd()`.**

``` r
setwd("~/Desktop/MyAnalysis")                  # breaks on every other machine
read.csv("/Users/anna/Documents/data/x.csv")   # same problem
```

Every path here is built from the project root instead:

``` r
here::here("data", "raw", "Censo_clean_27042026.csv")
```

`here::here()` finds the root by looking for the `.Rproj` file. The same line works on your laptop, a colleague's Windows machine, and a cluster.

------------------------------------------------------------------------

## Structure

```         
.
├── data/
│   ├── raw/         # exactly as received. READ-ONLY. Never edit these.
│   ├── clean/       # same information, but corrected (types, typos, known errors)
│   └── processed/   # analysis-ready: joined tables, derived variables
├── utils/
│   ├── functions.R  # reusable functions. Defines things; never DOES things.
│   └── workflow.R   # runs the whole project, in order, from scratch
├── scripts/         # numbered, each does ONE thing
├── output/
│   ├── figures/     # generated — safe to delete, workflow.R rebuilds them
│   └── results/     # model summaries, tables
├── renv.lock        # exact package versions
└── README.md
```

**Why `raw/` is separate and read-only.** It's the only copy nobody has touched. If a cleaning decision turns out to be wrong six months from now, you re-run the cleaning script and get a corrected `clean/`. If you edit the file in place, that history is gone forever.

**Why `clean/` and `processed/` are separate.** `clean/` is the *same information*, corrected. `processed/` is *new information* — tables joined together, variables derived for a specific model. Keeping them apart means changing an analysis decision doesn't force you to re-verify your data corrections.

------------------------------------------------------------------------

## The scripts

Run in order. Each reads its input from disk and writes its output to disk — none of them depends on objects left in your environment by an earlier one.

### `00_diagnose_data.R` — how the problems were found

Reads `data/raw/` and writes nothing. This is an audit script, not part of the pipeline. Every problem it shows is already fixed by `01_clean_data.R`.

It exists because **none of these problems announce themselves.** Every one produces code that runs without error and gives you a wrong answer.

### `01_clean_data.R` — raw → clean → processed

Applies every correction, and **flags every problem it can't resolve** as a warning that fires on each run (see [Open questions](#open-questions-still-unresolved)).

### `02_glmm_herbivory.R` — two responses, three exercises

-   **Incidence** (proportion of leaves damaged) → binomial: successes out of a known number of trials
-   **Severity** (mean % area lost) → beta: a continuous proportion with no trial count

Collapsing both into one Gaussian `lm()` on "% damage" is the most common mistake in the herbivory literature. A plant can be attacked often but lightly, or rarely but severely and one averaged number can't tell those apart.

### `03_glmm_mortality.R` — same methods, but different question

Binary survival data, with random effect: `(1 | species_acronym:mother_tree)` is the genotype question T-REX's mixed-maternal-line design was built to test.

------------------------------------------------------------------------

## The exercises

Each is built the same way: **the wrong model first, then the right one**, so you can see what changes. In your copy the answers are included, but try covering them and writing the corrected model yourself first.

### Exercise 1 — does your random-effects structure match your design?

T-REX treatments are applied to whole **subblocks**. Sixteen trees in a subblock are not sixteen independent replicates. Fit it both ways and watch the standard errors.

The takeaway is not "always add random effects." It's: *what is the unit you actually manipulated or sampled?* That's what determines the structure.

### Exercise 2 — effective n

Leaf traits exist for 10 of 16 species. However many rows your data frame has, a species-level trait effect is replicated at the level of **species**, so the effective n is 10, not 3,400.

`(1 | species_acronym)` belongs in the model whether or not it looks "significant." Its job is to stop the trait term claiming power it doesn't have.

### Exercise 3

The dispersion story, in order:

| Step | Result |
|------------------------------------|------------------------------------|
| `testDispersion()` | ratio 1.07 — looks mild |
| `plot(simulateResiduals(...))` | dense outliers stacked at 0 and 1 |
| `testOutliers(type = "bootstrap")` | 6.1% observed vs 0.3% expected — **20× excess** |
| `AIC()` binomial vs beta-binomial | **ΔAIC ≈ 4551** |

The mild-looking global ratio was a bad guide, because trial sizes range from 1 to 3,587 leaves. Run the check.

------------------------------------------------------------------------

## Open questions (still unresolved) {#open-questions-still-unresolved}

These fire as warnings every time you run the pipeline. Ideally, in the future, they would be resolved with the data owners.

| What | Who to ask |
|------------------------------------|------------------------------------|
| 4 plots where `Date_tree_Control` predates planting for nearly every plant | Edith / CM |
| XP201-PR, XP202-PR: `Age_2026 = 0` but `Regeneration_year = 2024` | Felicity |
| XP220-PLR: age 50 but labelled late-regeneration, no regeneration year | Felicity |
| `CH-MT5-P53`: the raw data itself says the species ID is wrong | censo maintainer |
| Does block 5 exist? Field notes say blocks 1–4; the data disagree | Julie / CM |
| 3 plants with "MUERTA EN EL TRAYECTO" notes, counted as ordinary deaths | censo maintainer |

Known limitation, not resolved: follow-up time varies from 1 to 381 days between plants. It was tested as an offset in the mortality model but destabilised fitting, so it's excluded. Species with longer observation windows may appear more prone to mortality than they are.

------------------------------------------------------------------------

## Conventions worth copying

-   Scripts are **numbered and ordered**. Each writes to disk; the next reads from disk. No script depends on your environment.
-   Functions live in `utils/functions.R`, never in analysis scripts. Same three lines twice? Make it a function.
-   **Adding a package:** `install.packages("x")` → `renv::snapshot()` → commit `renv.lock`. Skipping the last two steps is the most common way this breaks for everyone else.
-   `output/` is disposable. If you can't regenerate it by re-running the workflow, it doesn't belong there.
-   Filenames: lowercase, underscores, no spaces, no accents. No dates in filenames — git tracks history for you.

------------------------------------------------------------------------

## If something breaks

**Most errors in a session like this come from running chunks out of order**, so an old object is sitting in your environment. Before debugging anything else:

``` r
rm(list = ls())
source(here::here("utils", "workflow.R"))
```

If it still fails from a clean state, that's a real bug worth chasing.
