## =============================================================================
## 03_glmm_mortality.R
##
## IN:  data/processed/trex_analysis.csv
## OUT: output/figures/, output/results/
##
## T-REX plants 16 species across 25 plots, with individuals of each species
## drawn from DIFFERENT MOTHER TREES and mixed within each plot, specifically
## so genotype effects can be tested.
## `(1 | mother_tree)` would be that test.
##
## Run via: source(here::here("utils", "workflow.R"))
## =============================================================================

library(dplyr)
library(ggplot2)
library(lme4)
library(DHARMa)
library(performance)

source(here::here("utils", "functions.R"))

theme_set(theme_bw(base_size = 12))


## -----------------------------------------------------------------------------
## Build the mortality-complete dataset
## -----------------------------------------------------------------------------

censo_clean <- readr::read_csv(here::here("data", "clean", "censo_clean.csv"),
                               show_col_types = FALSE)
land_clean  <- readr::read_csv(here::here("data", "clean", "landcover_clean.csv"),
                               show_col_types = FALSE)

dat_mort <- censo_clean |>
  left_join(land_clean, by = "plot_id") |>
  filter(!is.na(died)) |>
  mutate(age_scaled = as.numeric(scale(age_2026)),
         plot_id = factor(plot_id), block = factor(block),
         subblock = factor(subblock), species_acronym = factor(species_acronym),
         mother_tree = factor(mother_tree))

message(sum(dat_mort$species_id_flag), " plant(s) excluded: unreliable species ID.")
message(sum(dat_mort$risk_time_flag), " plant(s) excluded: individual negative risk time.")
message(sum(dat_mort$plot_control_date_flag), " plant(s) excluded: whole-plot control-date issue.")

dat_mort_offset <- dat_mort |>
  filter(!species_id_flag, !risk_time_flag, !plot_control_date_flag,
         !is.na(days_at_risk)) |>
  droplevels()

## =============================================================================
## PART A — EXPLORATORY PLOTS
## =============================================================================

mort_by_sp <- dat_mort_offset |>
  group_by(species_acronym) |>
  summarise(n = n(), n_died = sum(died), prop_died = n_died / n, .groups = "drop") |>
  mutate(se = sqrt(prop_died * (1 - prop_died) / n))

p_sp <- ggplot(mort_by_sp, aes(x = reorder(species_acronym, prop_died), y = prop_died)) +
  geom_col(fill = "grey70") +
  geom_errorbar(aes(ymin = pmax(prop_died - se, 0), ymax = pmin(prop_died + se, 1)),
                width = 0.25) +
  geom_text(aes(label = paste0("n=", n)), hjust = -0.2, size = 3) +
  coord_flip(clip = "off") +
  labs(x = "Species", y = "Proportion died",
       title = "Raw mortality by tree species",
       subtitle = "Error bars are naive binomial SEs — ignore plot structure and follow-up time")
p_sp
# save_figure(p_sp, "A1_mortality_by_species.png")

p_age_raw <- dat_mort_offset |>
  group_by(plot_id, age_2026) |>
  summarise(prop_died = mean(died), n = n(), .groups = "drop") |>
  ggplot(aes(x = age_2026, y = prop_died)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, colour = "firebrick") +
  labs(x = "Plot age (years)", y = "Proportion died",
       title = "Raw mortality vs. plot age",
       subtitle = "Motivates testing species x age below — not yet accounting for follow-up time")
p_age_raw
# save_figure(p_age_raw, "A2_mortality_by_age_raw.png")

mort_by_plot <- dat_mort_offset |>
  group_by(plot_id, treatment) |>
  summarise(n = n(), prop_died = mean(died), .groups = "drop")

p_plot <- ggplot(mort_by_plot, aes(x = reorder(plot_id, prop_died), y = prop_died,
                                   fill = treatment)) +
  geom_col() +
  coord_flip() +
  labs(x = "Plot", y = "Proportion died", fill = NULL,
       title = "Mortality varies substantially between plots",
       subtitle = "This spread is what (1 | plot_id/block/subblock) absorbs")
p_plot
# save_figure(p_plot, "A3_mortality_by_plot.png")

p_mt <- dat_mort_offset |>
  group_by(species_acronym, mother_tree) |>
  summarise(n = n(), prop_died = mean(died), .groups = "drop") |>
  filter(n >= 5) |>
  ggplot(aes(x = mother_tree, y = prop_died)) +
  geom_col(fill = "darkolivegreen4") +
  facet_wrap(~ species_acronym, scales = "free_x") +
  labs(x = "Mother tree", y = "Proportion died",
       title = "Mortality by maternal line, within species",
       subtitle = "Groups with n < 5 dropped. Variation here = the genotype signal")
p_mt
# save_figure(p_mt, "A4_mortality_by_mother_tree.png", width = 10, height = 7)


## =============================================================================
## PART B — THE MODEL
## =============================================================================
## - species x age (LRT-tested first, as in 02_glmm_herbivory.R)
## - offset(log(days_at_risk)): plants had different follow-up time before the
##   census; without this, species/plots planted earlier look falsely riskier
## - (1 | plot_id/block/subblock): same design structure as the herbivory model
## - (1 | species_acronym:mother_tree)

m_mort_noint <- glmer(
  died ~ species_acronym + age_scaled +
    (1 | plot_id/block/subblock) + (1 | species_acronym:mother_tree),
  data = dat_mort_offset, family = binomial(),
  control = glmerControl(optimizer = "nloptwrap")
)

m_mort <- glmer(
  died ~ species_acronym * age_scaled +
    (1 | plot_id/block/subblock) + (1 | species_acronym:mother_tree),
  data = dat_mort_offset, family = binomial(),
  control = glmerControl(optimizer = "nloptwrap")
)

message(">>> Follow-up time (days_at_risk) varies a lot across plants (1-381 ",
        "days) and was tested as an offset, but it destabilized model fitting ",
        "when combined with several near-zero-mortality genotype groups. ",
        "Dropped here, but a real limitation worth paying attention to. ")

message("\nDoes species x age improve on species + age for mortality?")
print(anova(m_mort_noint, m_mort, test = "LRT"))

print(summary(m_mort))
save_result(summary(m_mort), "mortality_glmm.txt")


## -----------------------------------------------------------------------------
## Model checking
## -----------------------------------------------------------------------------
## NOTE: `died` is one Bernoulli trial per plant, not an aggregated
## cbind(success, failure) with varying trial size like herbivory incidence.
## There is no beta-binomial-style fix here — for genuinely binary data the
## variance is pinned to the mean, so "overdispersion" isn't the same kind of
## identifiable, fixable thing. DHARMa's checks below still catch OTHER
## misspecification (missing predictors, wrong functional form) — that part
## still matters, just don't expect a family-switch fix if something's off.

mort_res <- simulateResiduals(m_mort)
print(testDispersion(mort_res))
print(testOutliers(mort_res, type = "bootstrap"))

# png(here::here("output", "figures", "RES_mortality_dharma.png"),
#     width = 900, height = 450)
plot(mort_res, quantreg = FALSE)
# dev.off()


## -----------------------------------------------------------------------------
## Fixed effects plot
## -----------------------------------------------------------------------------

co <- summary(m_mort)$coefficients
fix_est <- data.frame(term = rownames(co), est = co[, "Estimate"], se = co[, "Std. Error"],
                      row.names = NULL) |>
  filter(term != "(Intercept)") |>
  mutate(lwr = est - 1.96 * se, upr = est + 1.96 * se,
         term = gsub("species_acronym", "", term))

p_fix <- ggplot(fix_est, aes(x = reorder(term, est), y = est)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(ymin = lwr, ymax = upr)) +
  coord_flip() +
  labs(x = NULL, y = "Effect on mortality (log-odds)",
       title = "Mortality GLMM: species, age, and their interaction")
p_fix
# save_figure(p_fix, "RES_mortality_species_effects.png")


## -----------------------------------------------------------------------------
## Variance components — plot/block/subblock vs. genotype
## -----------------------------------------------------------------------------

vc <- as.data.frame(VarCorr(m_mort))
vc_out <- data.frame(group = vc$grp, variance = round(vc$vcov, 4), sd = round(vc$sdcor, 4))
print(vc_out)
# save_result(vc_out, "mortality_variance_components.csv")

p_vc <- ggplot(vc_out, aes(x = reorder(group, sd), y = sd)) +
  geom_col(fill = "steelblue", width = 0.6) +
  coord_flip() +
  labs(x = NULL, y = "Standard deviation (log-odds scale)",
       title = "Where does the variation in mortality sit?",
       subtitle = "species:mother_tree vs. plot/block/subblock — genotype vs. environment")
p_vc
# save_figure(p_vc, "RES_mortality_variance_components.png")

## -----------------------------------------------------------------------------
## Predicted trajectories: species-specific mortality risk with age
## -----------------------------------------------------------------------------
## No offset in the final model (see note above — days_at_risk destabilized
## fitting and was dropped). Predictions are simply "probability of dying by
## the census," not standardised to a fixed exposure window.
age_range <- range(dat_mort_offset$age_2026, na.rm = TRUE)
pred_grid <- expand.grid(
  species_acronym = levels(dat_mort_offset$species_acronym),
  age_2026 = seq(age_range[1], age_range[2], length.out = 100)
)
age_mean <- mean(dat_mort_offset$age_2026, na.rm = TRUE)
age_sd   <- sd(dat_mort_offset$age_2026, na.rm = TRUE)
pred_grid$age_scaled <- (pred_grid$age_2026 - age_mean) / age_sd
pred_grid$fit <- predict(m_mort, newdata = pred_grid, type = "response", re.form = NA)

p_age_curves <- ggplot(pred_grid, aes(x = age_2026, y = fit, colour = species_acronym)) +
  geom_line(linewidth = 0.8) +
  labs(x = "Plot age (years)", y = "Predicted mortality probability",
       colour = "Species",
       title = "Species-specific mortality risk across the chronosequence")
p_age_curves
# save_figure(p_age_curves, "RES_mortality_age_by_species_curves.png", width = 9, height = 5.5)