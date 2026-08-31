## =============================================================================
## 03_glmm_herbivory.R
##
## IN:  data/processed/trex_analysis.csv
## OUT: output/figures/  (exploratory + results plots)
##      output/results/  (model summaries, comparisons)
##
##   EXERCISE 1  — matching random effects to the experimental design
##   EXERCISE 1B — adding age: interaction, convergence, and
##                 checking dispersion 
##   EXERCISE 2  — effective n when the predictor varies at the species level
##
## Run via: source(here::here("utils", "workflow.R"))
## =============================================================================

library(dplyr)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(performance)

source(here::here("utils", "functions.R"))

dat <- readr::read_csv(
  here::here("data", "processed", "trex_analysis.csv"),
  show_col_types = FALSE
) |>
  mutate(
    plot_id         = factor(plot_id),
    block           = factor(block),
    subblock        = factor(subblock),
    species_acronym = factor(species_acronym),
    mother_tree     = factor(mother_tree)
  )

theme_set(theme_bw(base_size = 12))


## -----------------------------------------------------------------------------
## Two responses (not one, as it may be assumed!)
## -----------------------------------------------------------------------------
##   INCIDENCE — what proportion of a plant's leaves were attacked at all?
##               successes out of a known number of trials -> binomial distribution
##   SEVERITY  — given damage, how much leaf area was lost?
##               continuous proportion, no trial count -> beta distribution (bounded)

dat_inc <- dat |>
  filter(!is.na(n_leaves), !is.na(n_leaves_damaged), n_leaves > 0) |>
  filter(!leaf_count_flag) |>   # excludes the 2 rows where damaged > total —
  # see the warning in 01_clean_data.R
  mutate(incidence = n_leaves_damaged / n_leaves)

dat_sev <- dat |> filter(!is.na(severity_prop_sv))

message("incidence: ", nrow(dat_inc), " plants | severity: ", nrow(dat_sev), " plants")
message("NOTE: ~28% of the herbivory census is incomplete (Julie's message) — ",
        "these models use complete cases only.")


## =============================================================================
## PART A — EXPLORATORY PLOTS
##
## Always look at the data before modelling. Each plot below previews a
## modelling decision made later.
## =============================================================================

## A1. The two responses are NOT the same variable ----------------------------
p_two_responses <- dat_inc |>
  filter(!is.na(severity_prop)) |>
  ggplot(aes(x = incidence, y = severity_prop)) +
  geom_point(alpha = 0.25, colour = "grey30") +
  geom_smooth(method = "loess", formula = y ~ x, colour = "firebrick") +
  labs(
    x = "Incidence (proportion of leaves damaged)",
    y = "Severity (mean % area lost, as proportion)",
    title = "Incidence and severity are different quantities",
    subtitle = "A plant can be attacked often but lightly, or rarely but severely"
  )
p_two_responses
# save_figure(p_two_responses, "A1_incidence_vs_severity.png")

## A2. Shape of each response -------------------------------------------------
p_dist_inc <- ggplot(dat_inc, aes(x = incidence)) +
  geom_histogram(bins = 30, fill = "steelblue", colour = "white") +
  labs(x = "Incidence", y = "Number of plants",
       title = "Incidence: bounded at 0 and 1")
p_dist_inc
# save_figure(p_dist_inc, "A2a_incidence_distribution.png")

p_dist_sev <- ggplot(dat_sev, aes(x = severity_prop)) +
  geom_histogram(bins = 30, fill = "darkorange", colour = "white") +
  labs(x = "Severity (proportion)", y = "Number of plants",
       title = "Severity: piled up near zero, hard floor at 0",
       subtitle = "Gaussian assumes symmetry and constant variance — neither holds")
p_dist_sev
# save_figure(p_dist_sev, "A2b_severity_distribution.png")


## A3. Variation between species ----------------------------------------------
p_species <- dat_inc |>
  ggplot(aes(x = reorder(species_acronym, incidence, FUN = median),
             y = incidence)) +
  geom_boxplot(outlier.alpha = 0.2, fill = "grey90") +
  coord_flip() +
  labs(x = "Species", y = "Incidence",
       title = "Herbivory incidence by tree species")
p_species
# save_figure(p_species, "A3_incidence_by_species.png")


## A4. Variation between plots — the reason random effects are needed ---------
p_plots <- dat_inc |>
  ggplot(aes(x = reorder(plot_id, incidence, FUN = median), y = incidence)) +
  geom_boxplot(outlier.alpha = 0.2, fill = "grey90") +
  coord_flip() +
  labs(x = "Plot", y = "Incidence",
       title = "Between-plot variation in incidence",
       subtitle = "Structure like this is why (1 | plot_id/block/subblock) exists")
p_plots
# save_figure(p_plots, "A4_incidence_by_plot.png")


## A5. Trait coverage — sets up Exercise 2 ------------------------------------
p_coverage <- dat |>
  group_by(species_acronym) |>
  summarise(n_plants = n(),
            n_with_traits = sum(!is.na(toughness)), .groups = "drop") |>
  tidyr::pivot_longer(c(n_plants, n_with_traits)) |>
  ggplot(aes(x = species_acronym, y = value, fill = name)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(n_plants = "grey70", n_with_traits = "firebrick"),
                    labels = c("all plants", "with trait data"), name = NULL) +
  labs(x = "Species", y = "Number of plants",
       title = "Trait data exist for only some species",
       subtitle = "The effective n for a species-level trait effect is the number of RED bars")
p_coverage
# save_figure(p_coverage, "A5_trait_coverage.png")


## =============================================================================
## EXERCISE 1 — Does your random-effects structure match the design?
##
## T-REX treatments are applied to whole SUBBLOCKS, not individual trees.
## 16 trees in a subblock are not 16 independent replicates.
##
## Thanks, Julie!
## =============================================================================

## WRONG: every plant treated as an independent replicate
m_inc_naive <- glm(
  cbind(n_leaves_damaged, n_leaves_undamaged) ~ species_acronym,
  data = dat_inc, family = binomial()
)
summary(m_inc_naive)

m_inc_mixed <- glmer(
  cbind(n_leaves_damaged, n_leaves_undamaged) ~ species_acronym +
    (1 | plot_id/block/subblock),
  data = dat_inc, family = binomial(),
  control = glmerControl(optimizer = "bobyqa")
)
summary(m_inc_mixed)

se_naive <- summary(m_inc_naive)$coefficients[, "Std. Error"]
se_mixed <- summary(m_inc_mixed)$coefficients[, "Std. Error"]

se_compare <- data.frame(
  term      = names(se_naive),
  se_naive  = round(se_naive, 4),
  se_mixed  = round(se_mixed[names(se_naive)], 4),
  inflation = round(se_mixed[names(se_naive)] / se_naive, 2),
  row.names = NULL
)
print(se_compare)
# save_result(se_compare, "ex1_se_inflation.csv")

## CAVEAT: glm() and glmer() coefficients for a binomial/logit model
## are not on exactly the same scale (glmer = conditional
## on the random effects; glm = population-averaged). This ratio mixes the
## pseudoreplication effect with a small scale change — a 
## illustration of direction and rough magnitude, not a formal statistic.
## The design-effect cross-check just below sidesteps that issue entirely.

p_inflation <- se_compare |>
  filter(term != "(Intercept)") |>
  ggplot(aes(x = reorder(term, inflation), y = inflation)) +
  geom_col(fill = "firebrick") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  labs(x = NULL, y = "SE(mixed) / SE(naive)",
       title = "How much the naive model understates uncertainty",
       subtitle = "Dashed line = no inflation. Above it = the naive CI was too narrow.")
p_inflation
# save_figure(p_inflation, "EX1_se_inflation.png")

print(check_overdispersion(m_inc_mixed))
# save_result(summary(m_inc_mixed), "ex1_incidence_mixed.txt")


## =============================================================================
## EXERCISE 1B — Does the herbivory-species relationship change with age?
##
## In an ideal world we would have environmental variables that change with
## time, rather than time only. E.g. canopy cover and soil moisture.
## =============================================================================

## Naive (no random effects yet) — but WITH the interaction, to test it first
m_inc_noint_age <- glm(
  cbind(n_leaves_damaged, n_leaves_undamaged) ~ species_acronym + age_2026,
  data = dat_inc, family = binomial()
)
m_inc_naive_age <- glm(
  cbind(n_leaves_damaged, n_leaves_undamaged) ~ species_acronym * age_2026,
  data = dat_inc, family = binomial()
)

## Test the interaction AS A WHOLE before reading 15 individual interaction
## p-values — scanning that many at once will find "significant" ones by
## chance alone even if nothing is really going on.
message("\nEXERCISE 1B: does species x age improve on species + age?")
print(anova(m_inc_noint_age, m_inc_naive_age, test = "LRT")) # Yes!

## Rescale age before fitting the mixed version. age_2026 runs 0-50 while
## every other predictor is a 0/1 dummy; multiplying a 0-50 variable against
## dummies to build interaction terms creates a poorly-scaled optimization
## surface. This is why the unscaled version below hit a convergence warning
## during development of this script.
## Please be aware that this is not a property of the data, but a property of the scale.
dat_inc <- dat_inc |> mutate(age_scaled = as.numeric(scale(age_2026)))

## optimizer = "nloptwrap": NOT glmer's default (that is bobyqa then
## Nelder_Mead). This was chosen here because allFit() comparison across five
## optimizers on this exact model showed nloptwrap.NLOPT_LN_BOBYQA converged
## to the same estimates (log-likelihoods agreed to 3 decimal places) roughly
## 10x faster than plain bobyqa. Again, an empirical choice for this model, not a
## rule to apply everywhere.
m_inc_mixed_age <- glmer(
  cbind(n_leaves_damaged, n_leaves_undamaged) ~ species_acronym * age_scaled +
    (1 | plot_id/block/subblock),
  data = dat_inc, family = binomial(),
  control = glmerControl(optimizer = "nloptwrap")
)

## Confirming convergence is trustworthy

run_allfit_crosscheck <- TRUE
if (run_allfit_crosscheck) {
  message("\nCross-checking convergence across optimizers (this takes a while)...")
  fits <- allFit(m_inc_mixed_age)
  print(summary(fits)$which.OK)
  print(range(unlist(summary(fits)$llik)))  # should agree to several decimals
}

## Check dispersion properly: do not assume the binomial variance is correct
## just because the model converged.
inc_age_res <- simulateResiduals(m_inc_mixed_age)
print(testDispersion(inc_age_res))
print(check_overdispersion(m_inc_mixed_age))
# png(here::here("output", "figures", "EX1B_dharma_dispersion.png"),
#     width = 700, height = 500)
plot(inc_age_res, quantreg = FALSE)
# dev.off()
testOutliers(inc_age_res, type = "bootstrap")
plotResiduals(inc_age_res, form = dat_inc$n_leaves, xlab = "n_leaves (trial size)")

mf <- model.frame(m_inc_mixed_age)
resp <- model.response(mf)                 # matrix: columns = damaged, undamaged
n_leaves_used <- rowSums(resp)              # reconstructs n_leaves, guaranteed aligned

plotResiduals(inc_age_res, form = n_leaves_used, xlab = "n_leaves (trial size)")

## If the ratio above is meaningfully > 1 (not just "significant" — with
## ~3000+ rows even a trivial ratio can hit p < 0.001, see EX1B_dharma note
## below): the beta-binomial is the direct fix for a cbind(success, failure)
## response, preferred over an observation-level random effect for binomial
## data specifically (Harrison 2015, PeerJ 3:e1114, compares the two head to
## head). Same formula, but with one extra dispersion parameter.
m_inc_mixed_age_bb <- glmmTMB(
  cbind(n_leaves_damaged, n_leaves_undamaged) ~ species_acronym * age_scaled +
    (1 | plot_id/block/subblock),
  data = dat_inc, family = betabinomial(link = "logit")
)
message("\nAIC comparison — binomial vs beta-binomial:")
print(AIC(m_inc_mixed_age, m_inc_mixed_age_bb))
message(">>> Lower AIC wins. If they're close, prefer the simpler (binomial) model.")

# save_result(summary(m_inc_mixed_age), "ex1b_incidence_age_mixed.txt")

## Confirm the species pattern directly
dat_inc |>
  group_by(species_acronym) |>
  summarise(median_n_leaves = median(n_leaves), max_n_leaves = max(n_leaves), .groups = "drop") |>
  arrange(desc(median_n_leaves))

## -----------------------------------------------------------------------------
## The beta-binomial fits dramatically better (AIC drop ~4551) — this is now
## the adopted model, not a comparison candidate. Re-run the same diagnostic
## sequence on IT, not on the discarded binomial model, to confirm the fix
## actually worked rather than assuming it did because AIC improved.
## -----------------------------------------------------------------------------
bb_res <- simulateResiduals(m_inc_mixed_age_bb)
print(testDispersion(bb_res))
print(testOutliers(bb_res, type = "bootstrap"))

# png(here::here("output", "figures", "EX1B_dharma_betabinomial.png"),
#     width = 700, height = 500)
plot(bb_res, quantreg = FALSE)
# dev.off()

mf_bb <- model.frame(m_inc_mixed_age_bb)
n_leaves_bb <- rowSums(model.response(mf_bb))
# png(here::here("output", "figures", "EX1B_dharma_betabinomial_vs_nleaves.png"),
#     width = 700, height = 500)
plotResiduals(bb_res, form = n_leaves_bb, xlab = "n_leaves (trial size)")
# dev.off()

## This is now the model to report and save.
# save_result(summary(m_inc_mixed_age_bb), "ex1b_incidence_age_betabinomial.txt")

## Results plot: predicted species-specific trajectories with age -------------
## This is the actual answer to "does it depend on species and age,
## considering their interaction" — far more direct than 29 coefficients.
## re.form = NA: population-level prediction, ignoring plot/block/subblock
## random effects (we want the average trajectory, not one specific plot's).

age_range <- range(dat_inc$age_2026, na.rm = TRUE)
pred_grid <- expand.grid(
  species_acronym = levels(dat_inc$species_acronym),
  age_2026 = seq(age_range[1], age_range[2], length.out = 100)
)
age_mean <- mean(dat_inc$age_2026, na.rm = TRUE)
age_sd   <- sd(dat_inc$age_2026, na.rm = TRUE)
pred_grid$age_scaled <- (pred_grid$age_2026 - age_mean) / age_sd
pred_grid$fit <- predict(m_inc_mixed_age_bb, newdata = pred_grid,
                         type = "response", re.form = NA)   # <- swapped model

p_age_curves <- ggplot(pred_grid, aes(x = age_2026, y = fit, colour = species_acronym)) +
  geom_line(linewidth = 0.8) +
  labs(x = "Plot age (years)", y = "Predicted incidence",
       colour = "Species",
       title = "Species-specific herbivory trajectories with forest age",
       subtitle = "Converging lines = species differences shrink with recovery (H3b); parallel/diverging = they don't")
p_age_curves
# save_figure(p_age_curves, "EX1B_age_by_species_curves.png", width = 9, height = 5.5)



## =============================================================================
## SEVERITY AS A BETA GLMM
## =============================================================================

m_sev <- glmmTMB(
  severity_prop_sv ~ species_acronym + (1 | plot_id/block/subblock),
  data = dat_sev, family = beta_family(link = "logit")
)
summary(m_sev)
# save_result(summary(m_sev), "severity_beta_glmm.txt")

co_sev <- summary(m_sev)$coefficients$cond
sev_est <- data.frame(
  term = rownames(co_sev),
  est  = co_sev[, "Estimate"],
  se   = co_sev[, "Std. Error"],
  row.names = NULL
) |>
  filter(term != "(Intercept)") |>
  mutate(
    lwr = est - 1.96 * se,
    upr = est + 1.96 * se,
    term = gsub("species_acronym", "", term)
  )

p_sev_est <- ggplot(sev_est, aes(x = reorder(term, est), y = est)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(ymin = lwr, ymax = upr)) +
  coord_flip() +
  labs(x = "Species (vs. reference level)", y = "Effect on severity (logit scale)",
       title = "Severity beta-GLMM: species effects",
       subtitle = "Points crossing the dashed line are not distinguishable from the reference")
save_figure(p_sev_est, "RES_severity_species_effects.png")

sev_res <- simulateResiduals(m_sev)
# png(here::here("output", "figures", "RES_severity_dharma.png"),
#     width = 900, height = 450)
plot(sev_res)
# dev.off()
# message("  saved figure: RES_severity_dharma.png")

## Now what would you do to improve this model fit?

## =============================================================================
## EXERCISE 2 — Effective n
##
## Thanks, Hanna!
## =============================================================================

dat_tr <- dat_sev |> filter(!is.na(toughness), !is.na(chl_ug_cm2))

message("\nEXERCISE 2: ", nrow(dat_tr), " plants, but only ",
        n_distinct(dat_tr$species_acronym), " species carry trait values.")
message(">>> The effective n for a species-level trait effect is that second number.")

## WRONG: toughness borrows precision from every leaf-damage observation
m_tr_naive <- glmmTMB(
  severity_prop_sv ~ toughness + chl_ug_cm2 + (1 | plot_id/block/subblock),
  data = dat_tr, family = beta_family(link = "logit")
)

m_tr_mixed <- glmmTMB(
  severity_prop_sv ~ toughness + chl_ug_cm2 +
    (1 | species_acronym) + (1 | plot_id/block/subblock),
  data = dat_tr, family = beta_family(link = "logit")
)

co_naive <- summary(m_tr_naive)$coefficients$cond
co_mixed <- summary(m_tr_mixed)$coefficients$cond

trait_compare <- data.frame(
  term      = rownames(co_naive),
  est_naive = round(co_naive[, "Estimate"], 3),
  se_naive  = round(co_naive[, "Std. Error"], 3),
  p_naive   = signif(co_naive[, "Pr(>|z|)"], 3),
  est_mixed = round(co_mixed[rownames(co_naive), "Estimate"], 3),
  se_mixed  = round(co_mixed[rownames(co_naive), "Std. Error"], 3),
  p_mixed   = signif(co_mixed[rownames(co_naive), "Pr(>|z|)"], 3),
  row.names = NULL
)
print(trait_compare)
# save_result(trait_compare, "ex2_trait_effective_n.csv")

p_trait <- trait_compare |>
  filter(term != "(Intercept)") |>
  tidyr::pivot_longer(
    cols = c(est_naive, se_naive, est_mixed, se_mixed),
    names_to = c(".value", "model"), names_sep = "_"
  ) |>
  mutate(lwr = est - 1.96 * se, upr = est + 1.96 * se,
         model = factor(model, levels = c("naive", "mixed"),
                        labels = c("without (1|species)", "with (1|species)"))) |>
  ggplot(aes(x = term, y = est, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(ymin = lwr, ymax = upr),
                  position = position_dodge(width = 0.4)) +
  coord_flip() +
  scale_colour_manual(values = c("firebrick", "steelblue"), name = NULL) +
  labs(x = NULL, y = "Effect on severity (logit scale)",
       title = "Exercise 2: what the species random effect does to trait effects",
       subtitle = "Same data and predictors, but the replication level differs")
p_trait
# save_figure(p_trait, "EX2_trait_effective_n.png")

message("\n>>> Discussion point: whether or not the species random effect looks ",
        "'significant', it belongs in the model. Its job is to stop the trait ",
        "term from claiming statistical power it does not have.")