## ----setup, include=FALSE-----------------------------------------------------
knitr::knit_hooks$set(purl = knitr::hook_purl)
knitr::opts_chunk$set(echo = TRUE)

## ----Load packages & data-----------------------------------------------------
library(dplyr)
library(robvis)
library(metafor)
library(readxl)
library(meta)


## ----Univariate analysis------------------------------------------------------

# ---------------------------------------------------------------
## 1. Define assumed pre–post correlation and helper function
# ---------------------------------------------------------------

r_within <- 0.50   

sd_change <- function(sd_pre, sd_post, r = r_within) {
  sqrt(sd_pre^2 + sd_post^2 - 2 * r * sd_pre * sd_post)
}

# ---------------------------------------------------------------
## 2. Compute pre–post change scores for each group and outcome
# ---------------------------------------------------------------

dat <- dat %>%
  mutate(
    ## PTSD
    PTSD_mean_int1_change = PTSD_mean_int1_post - PTSD_mean_int1_pre,
    PTSD_sd_int1_change   = sd_change(PTSD_sd_int1_pre, PTSD_sd_int1_post, r_within),
    PTSD_mean_ctrl_change = PTSD_mean_ctrl_post - PTSD_mean_ctrl_pre,
    PTSD_sd_ctrl_change   = sd_change(PTSD_sd_ctrl_pre, PTSD_sd_ctrl_post, r_within),
    
    ## Depression
    Depression_mean_int1_change = Depression_mean_int1_post - Depression_mean_int1_pre,
    Depression_sd_int1_change   = sd_change(Depression_sd_int1_pre, Depression_sd_int1_post, r_within),
    Depression_mean_ctrl_change = Depression_mean_ctrl_post - Depression_mean_ctrl_pre,
    Depression_sd_ctrl_change   = sd_change(Depression_sd_ctrl_pre, Depression_sd_ctrl_post, r_within),
    
    ## Anxiety
    Anxiety_mean_int1_change = Anxiety_mean_int1_post - Anxiety_mean_int1_pre,
    Anxiety_sd_int1_change   = sd_change(Anxiety_sd_int1_pre, Anxiety_sd_int1_post, r_within),
    Anxiety_mean_ctrl_change = Anxiety_mean_ctrl_post - Anxiety_mean_ctrl_pre,
    Anxiety_sd_ctrl_change   = sd_change(Anxiety_sd_ctrl_pre, Anxiety_sd_ctrl_post, r_within)
  )

# ---------------------------------------------------------------
## 3. Compute effect sizes based on change scores
# ---------------------------------------------------------------

# PTSD
c.ptsd <- metacont(
  n.e = n_int1_post,
  mean.e = PTSD_mean_int1_change,
  sd.e = PTSD_sd_int1_change,
  n.c = n_ctrl_post,
  mean.c = PTSD_mean_ctrl_change,
  sd.c = PTSD_sd_ctrl_change,
  studlab = study_id,
  data = dat,
  sm = "SMD", method.smd = "Hedges",
  random = TRUE, common = FALSE,
  method.tau = "REML",
  method.random.ci = "HK",
  title = "PTSD (Change Scores)"
)
summary(c.ptsd)
meta::forest(c.ptsd, prediction = TRUE)

# Depression
c.dep <- metacont(
  n.e = n_int1_post,
  mean.e = Depression_mean_int1_change,
  sd.e = Depression_sd_int1_change,
  n.c = n_ctrl_post,
  mean.c = Depression_mean_ctrl_change,
  sd.c = Depression_sd_ctrl_change,
  studlab = study_id,
  data = dat,
  sm = "SMD", method.smd = "Hedges",
  random = TRUE, common = FALSE,
  method.tau = "REML",
  method.random.ci = "HK",
  title = "Depression (Change Scores)"
)
summary(c.dep)
meta::forest(c.dep, prediction = TRUE)

# Anxiety
c.anx <- metacont(
  n.e = n_int1_post,
  mean.e = Anxiety_mean_int1_change,
  sd.e = Anxiety_sd_int1_change,
  n.c = n_ctrl_post,
  mean.c = Anxiety_mean_ctrl_change,
  sd.c = Anxiety_sd_ctrl_change,
  studlab = study_id,
  data = dat,
  sm = "SMD", method.smd = "Hedges",
  random = TRUE, common = FALSE,
  method.tau = "REML",
  method.random.ci = "HK",
  title = "Anxiety (Change Scores)"
)
summary(c.anx)
meta::forest(c.anx, prediction = TRUE)




## ----Multivariate Analysis----------------------------------------------------
# ---------------------------------------------------------------
# 1) Pull study-level effects from each meta object
# ---------------------------------------------------------------

es_ptsd <- data.frame(
  study_id = c.ptsd$studlab,
  outcome  = "PTSD",
  yi       = as.numeric(c.ptsd$TE),          # Hedges' g per study
  vi       = as.numeric(c.ptsd$seTE)^2       # sampling variance
)

es_dep <- data.frame(
  study_id = c.dep$studlab,
  outcome  = "Depression",
  yi       = as.numeric(c.dep$TE),
  vi       = as.numeric(c.dep$seTE)^2
)

es_anx <- data.frame(
  study_id = c.anx$studlab,
  outcome  = "Anxiety",
  yi       = as.numeric(c.anx$TE),
  vi       = as.numeric(c.anx$seTE)^2
)

dat_es <- bind_rows(es_ptsd, es_dep, es_anx) %>%
  mutate(eff_id = paste(study_id, outcome, sep = "_"))

# ---------------------------------------------------------------
# 2) Multilevel meta-analysis (effects nested within studies)
# ---------------------------------------------------------------

fit_all <- rma.mv(yi, vi,
                  random = ~ 1 | study_id/eff_id,   # within-study dependence
                  data   = dat_es,
                  method = "REML")
summary(fit_all)

#Forest plot
forest(fit_all,
       slab = paste(dat_es$study_id, dat_es$outcome, sep = ", "),
       xlab = "Hedges' g (negative = greater pre–post improvement for intervention)",
       alim = c(-2, 2),
       at = seq(-2, 2, 0.5),
       header = "Study, Outcome",
       mlab = "Overall pooled effect (multilevel model)")
abline(v = 0, lty = 2, col = "gray50")

# ---------------------------------------------------------------
# 3) Publication bias diagnostics (with study labels)
# ---------------------------------------------------------------

# a) Funnel plot with study IDs
funnel(fit_all,
       xlab = "Hedges' g",
       ylab = "Standard Error",
       main = "Funnel Plot (All Outcomes)",
       pch = 21, bg = "lightblue", cex = 1.2)

# Add study labels
text(x = dat_es$yi,
     y = sqrt(dat_es$vi),           # SE = sqrt(vi)
     labels = dat_es$study_id,
     cex = 0.7, pos = 4, col = "darkblue")

abline(v = fit_all$b, lty = 2, col = "gray40")

# b) Univariate random-effects model (for bias tests)
fit_uni <- rma(yi, vi, data = dat_es, method = "REML")

# c) Egger’s regression test (now on the rma model)
egger_test <- regtest(fit_uni, predictor = "sei")
egger_test

# d) Trim and fill
tf <- trimfill(fit_uni)
summary(tf)

# e) Funnel plot after trim and fill with labels
funnel(tf,
       xlab = "Hedges' g",
       main = "Trim and Fill Funnel Plot",
       pch = 21, bg = "lightgreen", cex = 1.2)

text(x = dat_es$yi,
     y = sqrt(dat_es$vi),
     labels = dat_es$study_id,
     cex = 0.7, pos = 4, col = "darkgreen")

abline(v = fit_uni$b, lty = 2, col = "gray40")

# f) Compare pooled estimates before vs after trim & fill
cat("\n--- Publication Bias Sensitivity ---\n")
cat("Original pooled g:", round(fit_uni$b, 3),
    "\nTrim & fill pooled g:", round(tf$b, 3),
    "\nNumber of filled studies:", tf$k0, "\n")

#ROB plot 
rob_data <- read_excel("trauma_data_set.xlsx", sheet = "ROB")
rob_summary(data = rob_data, 
            tool = "ROB2",
            weighted = FALSE)




## ----Moderator Analyses-------------------------------------------------------

# ---------------------------------------------------------------
## 1. Store study info
# ---------------------------------------------------------------

study_info <- dat %>%
  distinct(study_id,
           age_mean_study,
           int1_support,            # e.g., "Self-Guided/Human Supported"
           control_type,        # e.g., "Active"/"Passive"
           app_type,            # e.g., "CBT","Mindfulness","Other"
           app_name,            # e.g., "PTSD Coach","Wysa", etc.
           PTSD_mean_int1_pre,       # baseline PTSD mean (if available)
           Depression_mean_int1_pre, # baseline Depression mean (if available)
           Anxiety_mean_int1_pre)    # baseline Anxiety mean (if available)

dat_es <- bind_rows(es_ptsd, es_dep, es_anx) %>%
  mutate(eff_id = paste(study_id, outcome, sep = "_")) %>%
  left_join(study_info, by = "study_id")

dat_es <- dat_es %>%
  # baseline aligned to the outcome on that row
  mutate(
    baseline_symptom_mean = dplyr::case_when(
      outcome == "PTSD"      ~ PTSD_mean_int1_pre,
      outcome == "Depression"~ Depression_mean_int1_pre,
      outcome == "Anxiety"   ~ Anxiety_mean_int1_pre,
      TRUE ~ NA_real_
    ),
    # moderators: center continuous (age, baseline intervention mean)
    age_c  = scale(age_mean_study, center = TRUE, scale = FALSE),
    base_c = scale(baseline_symptom_mean, center = TRUE, scale = FALSE),
    # support moderator
    support = factor(int1_support, levels = c("Self-guided","Human-supported")),
    
    # control subgroup
    control_type = factor(control_type, levels = c("Passive","Active")),
    
    # app subgroups
    app_cbt = factor(ifelse(tolower(app_type) == "cbt", "CBT", "Other"),
                     levels = c("Other","CBT")),
    app_ptsdcoach = factor(ifelse(tolower(app_name) %in% c("ptsd coach","ptsdcoach","ptsd_coach"),
                                  "PTSDCoach","Other"),
                           levels = c("Other","PTSDCoach"))
  )
fit_all <- rma.mv(yi, vi,
                  random = ~ 1 | study_id/eff_id,
                  data   = dat_es, method = "REML")
summary(fit_all)

# ---------------------------------------------------------------
## 2. Fit moderators
# ---------------------------------------------------------------

fit_mod <- rma.mv(yi, vi,
                  mods   = ~ 1 + age_c + base_c + support,
                  random = ~ 1 | study_id/eff_id,
                  data   = dat_es, method = "REML")
summary(fit_mod)

# Age only
fit_age <- rma.mv(yi, vi,
                  mods = ~ 1 + age_c,
                  random = ~ 1 | study_id/eff_id,
                  data = dat_es, method = "REML")
summary(fit_age)

# Baseline only
fit_base <- rma.mv(yi, vi,
                   mods = ~ 1 + base_c,
                   random = ~ 1 | study_id/eff_id,
                   data = dat_es, method = "REML")
summary(fit_base)

# Guidance only
fit_guid <- rma.mv(yi, vi,
                   mods = ~ 1 + support,
                   random = ~ 1 | study_id/eff_id,
                   data = dat_es, method = "REML")
summary(fit_guid)


#Cbt vs other apps 
fit_sub_cbt <- rma.mv(yi, vi,
                      mods   = ~ 0 + app_cbt,   # gives separate intercepts per level
                      random = ~ 1 | study_id/eff_id,
                      data   = dat_es, method = "REML")
summary(fit_sub_cbt)

# Wald contrast (CBT - Other)
anova(fit_sub_cbt, X = rbind(c(-1, 1)))

#PTSD coach vs others
fit_sub_ptsdcoach <- rma.mv(yi, vi,
                            mods   = ~ 0 + app_ptsdcoach,
                            random = ~ 1 | study_id/eff_id,
                            data   = dat_es, method = "REML")
summary(fit_sub_ptsdcoach)

anova(fit_sub_ptsdcoach, X = rbind(c(-1, 1)))  # (PTSDCoach - Other)

#control type 
fit_sub_ctrl <- rma.mv(yi, vi,
                       mods   = ~ 0 + control_type,
                       random = ~ 1 | study_id/eff_id,
                       data   = dat_es, method = "REML")
summary(fit_sub_ctrl)

anova(fit_sub_ctrl, X = rbind(c(-1, 1)))  # (Active - Waitlist)





