## =============================================================================
## utils/workflow.R
##
## Runs the entire project, in order, from a fresh R session:
##
##     source(here::here("utils", "workflow.R"))
##
## If this fails, the project is broken. That is the whole point — it is a
## reproducibility test you can run any time, not just a convenience.
##
## It is deliberately boring: no analysis lives here. This file only decides
## WHAT runs and in WHICH ORDER. The analysis lives in scripts/.
## =============================================================================


## -----------------------------------------------------------------------------
## 1. Find the project root
## -----------------------------------------------------------------------------
## here::here() locates the folder containing the .Rproj file. Every path in
## this project is built from it, so nothing depends on your working directory
## or on where your home folder is.
##
## Bootstrap problem: we need `here` before renv::restore() (step 2) gets a
## chance to install it. So install it directly here, just this once, outside
## renv's control — it's a tiny, stable package and this is the one dependency
## the project needs before it can manage its own dependencies.

if (!requireNamespace("here", quietly = TRUE)) {
  message("Installing 'here' (needed to locate the project root)...")
  install.packages("here", repos = "https://cloud.r-project.org")
}

message("Project root: ", here::here())


## -----------------------------------------------------------------------------
## 2. Check the environment
## -----------------------------------------------------------------------------
## renv does the real work here: opening the .Rproj already sourced
## renv/activate.R via .Rprofile, which activates this project's private
## library. All that's left is to sync it against renv.lock. This is a no-op
## (fast) if you're already in sync — safe to leave in every run.
##
## ONE-TIME SETUP (do this yourself, once, in this project):
##   renv::init()      # creates renv/activate.R, .Rprofile, first renv.lock
##   renv::snapshot()  # after installing the packages this project needs
## Commit renv.lock (and .Rprofile) to git. Never commit renv/library/.

if (file.exists(here::here("renv.lock"))) {
  if (!requireNamespace("renv", quietly = TRUE)) {
    message("Installing 'renv'...")
    install.packages("renv", repos = "https://cloud.r-project.org")
  }
  renv::restore(prompt = FALSE)
  
  ## Beyond restoring what's already locked: check whether the CODE now uses
  ## something renv.lock doesn't know about yet — the "someone added
  ## library(newpkg), installed it, and forgot to snapshot" case. This is a
  ## warning, not a stop(): the current run already succeeded (everything
  ## locked did restore), this only means new additions aren't committed yet.
  message("\nChecking renv.lock against what the project code actually uses...")
  renv::status(project = here::here())
  message(
    "(If the report above lists anything as not-yet-in-the-lockfile: ",
    "run renv::snapshot() and commit the updated renv.lock.)\n"
  )
} else {
  stop(
    "No renv.lock found. This project expects one — see the ONE-TIME SETUP ",
    "note above utils/workflow.R, or ask whoever set up the repo to run ",
    "renv::snapshot() and commit renv.lock.",
    call. = FALSE
  )
}

## Include below whichever packages are missing and need to be installed
install.packages(c("DHARMa", "glmmTMB", "lme4", "performance"))
renv::snapshot()


## -----------------------------------------------------------------------------
## 3. Make sure the output folders exist
## -----------------------------------------------------------------------------
## output/ is git-ignored, so a fresh clone has no such folders. Create them
## rather than letting ggsave() fail on a missing directory.

for (d in list(
  c("data", "clean"),
  c("data", "processed"),
  c("output", "figures"),
  c("output", "results")
)) {
  dir.create(do.call(here::here, as.list(d)), recursive = TRUE, showWarnings = FALSE)
}


## -----------------------------------------------------------------------------
## 4. Load the project's own functions
## -----------------------------------------------------------------------------

source(here::here("utils", "functions.R"))


## -----------------------------------------------------------------------------
## 5. Run the analysis scripts, in order
## -----------------------------------------------------------------------------
## Each script is self-contained: it reads its input from disk and writes its
## output to disk. None of them depends on objects left in the environment by
## an earlier script. That is what makes it possible to re-run just one.

scripts <- c(
  "01_clean_data.R",      # data/raw/  -> data/clean/
  "02_glmm_herbivory.R",  # data/processed/ -> models, output/results/
  "03_glmm_mortality.R"   # data/processed/ -> models, output/results/
)

run_script <- function(filename) {
  path <- here::here("scripts", filename)
  if (!file.exists(path)) {
    stop("Script not found: ", filename, call. = FALSE)
  }
  message("\n", strrep("=", 70))
  message("RUNNING: ", filename)
  message(strrep("=", 70))
  
  started <- Sys.time()
  ## local = new.env() runs each script in its own environment, so a stray
  ## object in script 01 cannot silently satisfy a missing object in script 03.
  ## If a script only works when run after another, that is a bug you want to
  ## find here rather than three months from now.
  source(path, local = new.env(), echo = FALSE)
  elapsed <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1)
  
  message("DONE: ", filename, " (", elapsed, "s)")
  invisible(TRUE)
}

invisible(lapply(scripts, run_script))


## -----------------------------------------------------------------------------
## 6. Record the environment this run happened in
## -----------------------------------------------------------------------------
## Six months from now, "it gave a different number" is answerable only if you
## know which R and which package versions produced the original.

session_path <- here::here("output", "results", "session_info.txt")
utils::capture.output(sessionInfo(), file = session_path)

message("\n", strrep("=", 70))
message("WORKFLOW COMPLETE")
message("Outputs in: ", here::here("output"))
message("Environment recorded in: output/results/session_info.txt")
message(strrep("=", 70))