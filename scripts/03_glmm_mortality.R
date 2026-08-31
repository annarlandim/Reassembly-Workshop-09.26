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

source(here::here("utils", "functions.R"))

dat <- readr::read_csv(
  here::here("data", "processed", "trex_analysis.csv"),
  show_col_types = FALSE
) |>
  mutate(
    plot_id         = factor(plot_id),
    species_acronym = factor(species_acronym),
    mother_tree     = factor(mother_tree)
  )

theme_set(theme_bw(base_size = 12))


## -----------------------------------------------------------------------------
## The response variable
## -----------------------------------------------------------------------------
## `died` came from the censo column NAMED `Survival`, which is actually a
## mortality flag. Renamed in 01_clean_data.R so the model reads the way it behaves.

dat_mort <- dat |> filter(!is.na(died)) |>
  mutate(age_scaled = as.numeric(scale(age_2026)))

message("mortality model: ", nrow(dat_mort), " plants")
print(dat_mort |> count(died))
message("  (1 = died, 0 = alive at census)")


## =============================================================================
## PART A — EXPLORATORY PLOTS
## =============================================================================

## A1. Mortality by species ---------------------------------------------------

mort_by_sp <- dat_mort |>
  group_by(species_acronym) |>
  summarise(n = n(), n_died = sum(died),
            prop_died = n_died / n, .groups = "drop") |>
  ## binomial SE for a proportion, purely for the error bars in this plot
  mutate(se = sqrt(prop_died * (1 - prop_died) / n))

p_sp <- ggplot(mort_by_sp,
               aes(x = reorder(species_acronym, prop_died), y = prop_died)) +
  geom_col(fill = "grey70") +
  geom_errorbar(aes(ymin = pmax(prop_died - se, 0),
                    ymax = pmin(prop_died + se, 1)), width = 0.25) +
  geom_text(aes(label = paste0("n=", n)), hjust = -0.2, size = 3) +
  coord_flip(clip = "off") +
  labs(x = "Species", y = "Proportion died",
       title = "Raw mortality by tree species",
       subtitle = "Error bars are naive binomial SEs — they ignore plot structure")
p_sp
# save_figure(p_sp, "A1_mortality_by_species.png")


## A2. Mortality by plot ------------------------------------------------------
## Between-plot spread is the visual case for (1 | plot_id).

mort_by_plot <- dat_mort |>
  group_by(plot_id, treatment) |>
  summarise(n = n(), prop_died = mean(died), .groups = "drop")

p_plot <- ggplot(mort_by_plot,
                 aes(x = reorder(plot_id, prop_died), y = prop_died,
                     fill = treatment)) +
  geom_col() +
  coord_flip() +
  labs(x = "Plot", y = "Proportion died", fill = NULL,
       title = "Mortality varies substantially between plots",
       subtitle = "This spread is what (1 | plot_id) absorbs")
p_plot
# save_figure(p_plot, "A2_mortality_by_plot.png")


## A3. Mortality by mother tree, within species -------------------------------
## The genotype question, before any model. Each panel is one species; bars
## within a panel are different maternal lines of that species.

p_mt <- dat_mort |>
  group_by(species_acronym, mother_tree) |>
  summarise(n = n(), prop_died = mean(died), .groups = "drop") |>
  filter(n >= 5) |>
  ggplot(aes(x = mother_tree, y = prop_died)) +
  geom_col(fill = "darkolivegreen4") +
  facet_wrap(~ species_acronym, scales = "free_x") +
  labs(x = "Mother tree", y = "Proportion died",
       title = "Mortality by maternal line, within species",
       subtitle = "Groups with n < 5 plants dropped. Variation here = the genotype signal")
p_mt
# save_figure(p_mt, "A3_mortality_by_mother_tree.png", width = 10, height = 7)


## =============================================================================
## PART B — THE MODEL
## =============================================================================
## Species FIXED: 16 known, deliberately chosen levels; the contrasts are the
##   question.
## Plot RANDOM: 25 sampled sites standing in for a wider population.
## Mother tree RANDOM: the genotype question. Random not because we don't care,
##   but because we want its VARIANCE, not nine individual contrasts.

m_mort <- glmer(
  died ~ species_acronym * age_scaled + offset(log(days_at_risk)) +
    (1 | plot_id/block/subblock) + (1 | species_acronym:mother_tree),
  data = dat_mort, family = binomial(),
  control = glmerControl(optimizer = "nloptwrap")
)

print(summary(m_mort))
# save_result(summary(m_mort), "mortality_glmm.txt")


## B1. Fixed effects with CIs -------------------------------------------------

co <- summary(m_mort)$coefficients
fix_est <- data.frame(
  term = rownames(co), est = co[, "Estimate"], se = co[, "Std. Error"],
  row.names = NULL
) |>
  filter(term != "(Intercept)") |>
  mutate(lwr = est - 1.96 * se, upr = est + 1.96 * se,
         term = gsub("species_acronym", "", term))

p_fix <- ggplot(fix_est, aes(x = reorder(term, est), y = est)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(ymin = lwr, ymax = upr)) +
  coord_flip() +
  labs(x = "Species (vs. reference level)", y = "Effect on mortality (log-odds)",
       title = "Mortality GLMM: species effects",
       subtitle = "Above 0 = higher mortality than the reference species")
p_fix
# save_figure(p_fix, "RES_mortality_species_effects.png")


## B2. Variance components — the actual result of interest --------------------
## How much variation sits between mother trees vs. between plots? That
## comparison IS the genotype-versus-environment answer, and it comes from the
## variance estimates, not from a p-value.

vc <- as.data.frame(VarCorr(m_mort))
vc_out <- data.frame(
  group = vc$grp, variance = round(vc$vcov, 4), sd = round(vc$sdcor, 4)
)
print(vc_out)
save_result(vc_out, "mortality_variance_components.csv")

p_vc <- ggplot(vc_out, aes(x = reorder(group, sd), y = sd)) +
  geom_col(fill = "steelblue", width = 0.6) +
  coord_flip() +
  labs(x = NULL, y = "Standard deviation (log-odds scale)",
       title = "Where does the variation in mortality sit?",
       subtitle = "Taller bar = that grouping explains more variation")
save_figure(p_vc, "RES_mortality_variance_components.png")

message("\n>>> Compare the two SDs. A larger mother-tree SD than plot SD would ",
        "mean genotype matters more than site for early mortality — exactly ",
        "what T-REX's mixed-maternal-line design was built to detect.")


## B3. Random effect estimates (caterpillar plot) -----------------------------
## Which individual plots and maternal lines depart most from the average.

re <- ranef(m_mort, condVar = TRUE)
re_df <- bind_rows(lapply(names(re), function(g) {
  x  <- re[[g]]
  sd <- sqrt(as.numeric(attr(x, "postVar")))
  data.frame(grouping = g, level = rownames(x),
             est = x[, 1], se = sd, row.names = NULL)
}))

p_re <- re_df |>
  mutate(lwr = est - 1.96 * se, upr = est + 1.96 * se) |>
  group_by(grouping) |>
  arrange(est) |>
  mutate(rank = row_number()) |>
  ggplot(aes(x = rank, y = est)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(ymin = lwr, ymax = upr), size = 0.3) +
  facet_wrap(~ grouping, scales = "free") +
  labs(x = "Ranked level", y = "Random effect (log-odds)",
       title = "Caterpillar plot: departures from the average",
       subtitle = "Intervals crossing 0 are indistinguishable from the overall mean")
save_figure(p_re, "RES_mortality_random_effects.png", width = 9, height = 4.5)


## B4. Model checking ---------------------------------------------------------

mort_res <- simulateResiduals(m_mort)
png(here::here("output", "figures", "RES_mortality_dharma.png"),
    width = 900, height = 450)
plot(mort_res)
dev.off()
message("  saved figure: RES_mortality_dharma.png")


## -----------------------------------------------------------------------------
## CAVEAT to state when presenting
## -----------------------------------------------------------------------------
message("\n>>> Mother-tree IDs run MT1-MT9 but are only meaningful WITHIN a ",
        "species: 'MT4' of Cecropia and 'MT4' of Piper are unrelated trees. ",
        "As written, (1 | mother_tree) pools them — fine as a teaching example, ",
        "WRONG for a real genotype analysis. The correct term is ",
        "(1 | species_acronym:mother_tree). Good exercise: change it and see ",
        "what happens to the variance estimate.")