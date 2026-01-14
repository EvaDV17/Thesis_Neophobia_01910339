##Prepare data

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(readr)
  library(glmmTMB)
  library(lme4)
  library(emmeans)
  library(performance)
  library(DHARMa)
  library(broom.mixed)
  library(car)
  library(slider)
})

safe_AIC <- function(...) tryCatch(AIC(...), error = function(e) NA_real_)

library(dplyr)
library(tidyr)

#Load data
df <- read_csv("~/Downloads/rank_and_nn_outputLAST.csv") %>%
  mutate(
    BlockID   = str_extract(Group, "^B\\d+_\\d+_[A-Z]"),
    Condition = str_extract(Group, "(control|novel)") |>
      str_to_title() |>
      factor(c("Control","Novel"))
  )

Chick_data <- read_csv2("~/Downloads/ALL_CHICK_ID_FINAL.csv")
keep <- vapply(Chick_data, function(x) !all(is.na(x) | (is.character(x) & x == "")), logical(1))
Chick_data <- Chick_data[, keep, drop = FALSE]

#Merge data
merged_data <- left_join(Chick_data, df, by = "Bird_ID")

#Check assumptions
m1 <- lmer(nearest_neighbor_distance ~ Condition + (1|Chick_ID), data = merged_data)
summary(m1)

car::Anova(m1, type="II")

#Response: nearest_neighbor_distance
#Chisq Df Pr(>Chisq)
#Condition 0.8214  1     0.3648

#Assumptions are not good; log transformation seems to be a lot better; use this in the analysis

check_model(m1)

merged_data <- merged_data %>%
  rename(rank_from_food = rank_from_center) %>%  # just renaming the column
  mutate(
    logNND = log(nearest_neighbor_distance + 1)
  )

#Load fps (frames per second) per BlockID
fps_by_block <- tribble(
  ~BlockID,  ~fps,
  "B2_1_A",30, "B2_1_B",30,
  "B2_3_A",30, "B2_3_B",30,
  "B2_5_A",30, "B2_5_B",30,
  "B2_7_A",30, "B2_7_B",30,
  "B3_2_A",30, "B3_2_B",30,
  "B3_4_A",30, "B3_4_B",30,
  "B3_6_A",30, "B3_6_B",30,
  "B3_8_A",30, "B3_8_B",30,
  "B4_1_B",24,
  "B4_3_A",24, "B4_3_B",25,
  "B4_5_A",30, "B4_5_B",30,
  "B4_7_A",25, "B4_7_B",25,
  "B5_2_A",30, "B5_2_B",30,
  "B5_4_A",30, "B5_4_B",30,
  "B5_6_A",25, "B5_6_B",25,
  "B5_8_A",25, "B5_8_B",25,
  "B6_1_A",30, "B6_1_B",30,
  "B6_3_A",30, "B6_3_B",30,
  "B6_5_A",25, "B6_5_B",25,
  "B6_7_A",25, "B6_7_B",25,
  "B7_2_A",24, "B7_2_B",24,
  "B7_4_A",30, "B7_4_B",30,
  "B7_6_A",30, "B7_6_B",30,
  "B7_8_A",30, "B7_8_B",30,
  "B8_1_A",30, "B8_1_B",30,
  "B8_3_A",30, "B8_3_B",30,
  "B8_5_A",24, "B8_5_B",30,
  "B8_7_A",30, "B8_7_B",30,
  "B9_2_A",30, "B9_4_A",30,
  "B9_6_A",30, "B9_8_A",30
)

#Merge fps with dataset
merged_data <- merged_data %>%
  left_join(fps_by_block, by = "BlockID")

#Load groupsize per BlockID
groupsize_by_block <- tribble(
  ~BlockID,   ~groupsize,
  "B2_1_A",4, "B2_1_B",4, "B2_3_A",4, "B2_3_B",4,
  "B2_5_A",4, "B2_5_B",4, "B2_7_A",4, "B2_7_B",4,
  "B3_2_A",4, "B3_2_B",4, "B3_4_A",4, "B3_4_B",4,
  "B3_6_A",4, "B3_6_B",4, "B3_8_A",4, "B3_8_B",4,
  "B4_1_B",4, "B4_3_A",4, "B4_3_B",4, "B4_5_A",4, "B4_5_B",4,
  "B4_7_A",4, "B4_7_B",4,
  "B7_2_A",4, "B7_4_A",4, "B7_6_A",4, "B7_8_A",4,
  "B8_1_A",4, "B8_1_B",4, "B8_3_B",4, "B8_5_B",4, "B8_7_B",4,
  "B5_2_A",5, "B5_2_B",5, "B5_4_A",5, "B5_4_B",5,
  "B5_6_A",5, "B5_6_B",5, "B5_8_A",5, "B5_8_B",5,
  "B6_1_A",5, "B6_1_B",5, "B6_3_A",5, "B6_3_B",5,
  "B6_5_A",5, "B6_5_B",5, "B6_7_A",5, "B6_7_B",5,
  "B7_2_B",5, "B7_4_B",5, "B7_6_B",5, "B7_8_B",5,
  "B8_3_A",5, "B8_5_A",5, "B8_7_A",5,
  "B9_2_A",5, "B9_4_A",5, "B9_6_A",5, "B9_8_A",5
)

#Merge data
merged_data <- merged_data %>%
  left_join(groupsize_by_block, by = "BlockID") %>%
  mutate(
    groupsize = factor(groupsize),
    Condition = factor(Condition, levels = c("Control","Novel"))
  )

#Sanity checks
stopifnot(all(!is.na(merged_data$fps)))
stopifnot(all(!is.na(merged_data$groupsize)))

##Research questions
-----------
#R1:The first question explores whether the presence of a novel stimulus influences the group cohesion of herring gulls. 
-----------
  
  #Trial-level mean logNND with groupsize
  trial_data_rq1 <- merged_data %>%
  group_by(BlockID, Condition, groupsize) %>%
  summarise(
    mean_logNND = mean(logNND, na.rm = TRUE),
    n_frames    = n_distinct(frame),
    .groups     = "drop"
  ) %>%
  drop_na(mean_logNND, n_frames, groupsize, Condition)

#Possible models for AIC comparison
m0_rq1 <- lm(mean_logNND ~ 1, data = trial_data_rq1)
m1_rq1 <- lm(mean_logNND ~ Condition, data = trial_data_rq1)
m3_rq1 <- lm(mean_logNND ~ groupsize, data = trial_data_rq1)
m4_rq1 <- lm(mean_logNND ~ Condition + groupsize, data = trial_data_rq1)
m5_rq1 <- lm(mean_logNND ~ Condition * groupsize, data = trial_data_rq1)

aic_rq1 <- tibble(
  model = c("m0","m1","m3","m4","m5"),
  AIC   = c(safe_AIC(m0_rq1), safe_AIC(m1_rq1), safe_AIC(m3_rq1),
            safe_AIC(m4_rq1), safe_AIC(m5_rq1))
) %>% arrange(AIC)

cat("\n-- RQ1 AIC (cohesion) --\n"); print(aic_rq1)
#1 m0    -2.96  
#2 m1    -2.02  
#3 m3    -1.01  
#4 m4    -0.0832
#5 m5     1.58 

#Summary of the best model
summary(m0_rq1)

#The presence of a novel stimulus did not significantly change the spatial cohesion at the trial level, regardless of group size.

#Plot the group cohesion across conditions
ggplot(trial_data_rq1,
       aes(x = Condition, y = mean_logNND, fill = Condition)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.6, color = "grey30") +
  geom_jitter(width = 0.1, alpha = 0.5, size = 2, color = "black") +
  labs(
    title = "Group cohesion across conditions",
    x = "Condition",
    y = "Mean log(NND)"
  ) +
  scale_fill_manual(values = c("Control" = "#a6cee3", "Novel" = "#1f78b4")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

-------
#Preparation for R2 and R4
-------
min_hold_sec <- 0.5

# Detect overtakes (stable ≥ 0.5 s per trial-specific FPS)
ov <- merged_data %>%
  arrange(BlockID, Chick_ID, frame) %>%
  group_by(BlockID) %>%
  group_modify(function(d, key) {
    k <- max(1, round(dplyr::first(d$fps) * min_hold_sec))
    d %>%
      group_by(Chick_ID, .add = TRUE) %>%
      mutate(
        prev_rank  = dplyr::lag(rank_from_food),
        decreased  = rank_from_food < prev_rank,
        fut_max_k  = slider::slide_dbl(rank_from_food, max, .before = 0, .after = k - 1),
        stable     = fut_max_k <= prev_rank,
        overtook   = as.integer(decreased & stable)
      ) %>%
      ungroup()
  }) %>%
  ungroup()

#Summarise overtakes per trial
overtake_trial <- ov %>%
  group_by(BlockID, Condition) %>%
  summarise(
    n_overtakes = sum(overtook, na.rm = TRUE),
    n_frames    = dplyr::n_distinct(frame),
    fps         = dplyr::first(fps),
    n_seconds   = n_frames / fps,
    .groups     = "drop"
  ) %>%
  mutate(log_exposure = log(pmax(n_seconds, 1e-6))) %>%
  left_join(
    merged_data %>% distinct(BlockID, Condition, groupsize),
    by = c("BlockID","Condition")
  ) %>%
  drop_na(groupsize) %>%
  mutate(groupsize = factor(groupsize))

--------
  #R2: The second research question explores whether there are fewer or more overtaking events that take place in the presence of a novel stimulus.
--------
  
  #Possible models for AIC comparison
m0_rq2 <- glmmTMB(n_overtakes ~ 1 + (1|BlockID),
                    family = nbinom2, offset = log_exposure, data = overtake_trial)
m1_rq2 <- glmmTMB(n_overtakes ~ Condition + (1|BlockID),
                  family = nbinom2, offset = log_exposure, data = overtake_trial)
m2_rq2 <- glmmTMB(n_overtakes ~ groupsize + (1|BlockID),
                  family = nbinom2, offset = log_exposure, data = overtake_trial)
m3_rq2 <- glmmTMB(n_overtakes ~ Condition + groupsize + (1|BlockID),
                  family = nbinom2, offset = log_exposure, data = overtake_trial)
m4_rq2 <- glmmTMB(n_overtakes ~ Condition * groupsize + (1|BlockID),
                  family = nbinom2, offset = log_exposure, data = overtake_trial)

aic_rq2_nb <- tibble(
  model = c("m0","m1","m2","m3","m4"),
  AIC   = c(safe_AIC(m0_rq2), safe_AIC(m1_rq2), safe_AIC(m2_rq2),
            safe_AIC(m3_rq2), safe_AIC(m4_rq2))
) %>% arrange(AIC)

cat("\n-- RQ2 AIC NB --\n");  print(aic_rq2_nb)
#1 m4     415.
#2 m3     416.
#3 m2     418.
#4 m0     425.
#5 m1     425.

#m4 seems to be the best; suggestive that the difference between Control and Novel changes with group size? 

best_nb_rq2_name <- aic_rq2_nb$model[1]
best_nb_rq2 <- list(m0=m0_rq2,m1=m1_rq2,m2=m2_rq2,m3=m3_rq2,m4=m4_rq2)[[best_nb_rq2_name]]

if ("Condition" %in% all.vars(formula(best_nb_rq2))) {
  emm_rq2 <- emmeans(best_nb_rq2, ~ Condition | groupsize, type = "response", offset = 0)
  print(emm_rq2)
  print(pairs(emm_rq2, by = "groupsize"))
}
---------
  
  #Per group (4 or 5):
  
  #groupsize = 4:
  #Condition response    SE  df asymp.LCL asymp.UCL
  #Control       1.41 0.225 Inf      1.03      1.93
  #Novel         2.48 0.357 Inf      1.87      3.29
  
  #groupsize = 5:
  #Condition response    SE  df asymp.LCL asymp.UCL
  #Control       3.19 0.480 Inf      2.38      4.28
  #Novel         3.27 0.507 Inf      2.41      4.43
  
  #Confidence level used: 0.95 
  #Intervals are back-transformed from the log scale 
  
  --------------------
  #groupsize = 4:
  #contrast        ratio    SE  df null z.ratio p.value
  #Control / Novel 0.568 0.122 Inf    1  -2.638  0.0084
  
  #groupsize = 5:
  #contrast        ratio    SE  df null z.ratio p.value
  #Control / Novel 0.976 0.211 Inf    1  -0.111  0.9115
  
  #Tests are performed on the log scale 
  
  #Control: response ≈ 1.41 overtakes/s
  #Novel: response ≈ 2.48 overtakes/s
  #Contrast Control/Novel: ratio ≈ 0.568, p ≈ 0.0084
  #--> The novel condition has significantly ~1.76 (1/0.568) more overtakes per second than control in 4-bird groups. 
  #Thus, in groups of 4 birds, overtaking actually increased signifcantly under the novel object. 
  #This effect was not found for the 5-bird groups. 
  
  summary(m4_rq2)$coefficients

#                           Estimate  Std. Error  z value
#(Intercept)                0.3432579  0.1592935  2.154878
#ConditionNovel             0.5660170  0.2145970  2.637581
#groupsize5                 0.8169257  0.2190705  3.729055
#ConditionNovel:groupsize5 -0.5419923  0.3045621 -1.779579

#Pr(>|z|)
#(Intercept)               0.0311714107
#ConditionNovel            0.0083499694
#groupsize5                0.0001921996
#ConditionNovel:groupsize5 0.0751449037

  
#Visual of overtakes per second by condition and group size 
  
library(emmeans)

emm <- emmeans(m4_rq2, ~ Condition | groupsize, type = "response", offset = 0)
emm_df <- as.data.frame(emm)

ggplot(overtake_trial, aes(x = groupsize, y = overtakes_per_sec, fill = Condition)) +
  geom_boxplot(alpha = 0.45, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
    alpha = 0.35, size = 1
  ) +
  geom_pointrange(
    data = emm_df,
    aes(x = groupsize, y = response, ymin = asymp.LCL, ymax = asymp.UCL, color = Condition),
    position = position_dodge(width = 0.8),
    inherit.aes = FALSE
  ) +
  labs(
    title = "Overtakes per second by condition and group size",
    x = "Group size",
    y = "Overtakes per second"
  ) +
  theme_minimal(base_size = 14)

-------
  #R3: The third research question explores whether the same individuals are consistently at the front
-------
  
  #Per-bird summaries by condition
  leadership_condition <- merged_data %>%
  group_by(Chick_ID, Condition) %>%
  summarise(
    total_frames  = n(),
    front_frames  = sum(rank_from_food == 1, na.rm = TRUE),
    prop_front    = front_frames / total_frames,
    trials_present= n_distinct(BlockID),
    .groups = "drop"
  ) %>%
  arrange(Chick_ID, desc(prop_front))

#Summary by condition
leadership_condition %>%
  group_by(Condition) %>%
  summarise(
    mean_prop_front = mean(prop_front, na.rm = TRUE),
    sd_prop_front   = sd(prop_front,   na.rm = TRUE),
    .groups = "drop"
  )

#Wide table: led_Novel vs led_Control
leadership_wide <- leadership_condition %>%
  select(Chick_ID, Condition, prop_front) %>%
  pivot_wider(
    names_from  = Condition,
    values_from = prop_front,
    names_prefix = "led_"
  )

print(leadership_wide)

# A tibble: 67 × 3
#Chick_ID led_Novel led_Control
#<chr>        <dbl>       <dbl>
#1 BB_BB       0.499       0     
#2 BB_BP       0           0     
#3 BB_BY       0.484       0.451 
#4 BB_GR       0.110       0     
#5 BB_GY       0.368       0     
#6 BB_PR       0.55        0.0251
#7 BB_RR       0.0539      0.04  
#8 BG_BP       0.423       0.163 
#9 BG_BR       0.120       0.305 
#10 BG_BY       0           0 
#.....

#0 = never led in any trial
#0.5 = led in half of the trials (1/2)
#1 = led in all of the trials (2/2)

leadership_wide %>%
  mutate(type = case_when(
    led_Control == 0 & led_Novel == 0 ~ "Never leads",
    led_Control > 0 & led_Novel == 0  ~ "Leads only in control",
    led_Control == 0 & led_Novel > 0  ~ "Leads only in novel",
    TRUE                              ~ "Leads in both"
  )) %>%
  count(type)

# A tibble: 4 × 2
#type                        n
#<chr>                     <int>
#1 Leads in both            48
#2 Leads only in control     5
#3 Leads only in novel      11
#4 Never leads               3

#Of the 67 individuals, 48 led at least once in both contexts, 5 led only under control conditions, 
#11 led only under the novel condition, and 3 never led at all. 
#This indicates that leadership was not restricted to a stable smaller subset of individuals.

#Visual showing leadership consistency across conditions
ggplot(leadership_wide,
       aes(x = led_Control, y = led_Novel)) +
  geom_jitter(width = 0.03, height = 0.03, size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  labs(
    title = "Leadership Consistency Across Conditions",
    x = "Proportion leading (Control)",
    y = "Proportion leading (Novel)",
    subtitle = "Dashed line = identical leadership across conditions"
  ) +
  theme_minimal(base_size = 14)

#--> (1,1) in the graph are consistent leaders, regardless of condition
# Above diagonal: lead more in novel
#Below diagonal: lead more in control
# (0,0): little to no leadership

#Correlation and paired t-test
leadership_stats <- leadership_wide %>%
  filter(!is.na(led_Control) & !is.na(led_Novel))

cor_test <- cor.test(leadership_stats$led_Control,
                     leadership_stats$led_Novel, method = "pearson")
cat("\n=== Leadership Consistency — Correlation ===\n")
cat(sprintf("Pearson r = %.3f | t = %.3f | df = %d | p = %.4g\n",
            cor_test$estimate, cor_test$statistic, cor_test$parameter, cor_test$p.value))

#Pearson r = 0.137 | t = 1.118 | df = 65 | p = 0.2678

ttest_lead <- t.test(leadership_stats$led_Novel,
                     leadership_stats$led_Control, paired = TRUE)
cat("\n=== Leadership Consistency — Paired t-test ===\n")
cat(sprintf("Mean (Novel) = %.3f | Mean (Control) = %.3f | Mean diff = %.3f\n",
            mean(leadership_stats$led_Novel, na.rm=TRUE),
            mean(leadership_stats$led_Control, na.rm=TRUE),
            mean(leadership_stats$led_Novel - leadership_stats$led_Control, na.rm=TRUE)))
cat(sprintf("t = %.3f | df = %d | p = %.4g\n",
            ttest_lead$statistic, ttest_lead$parameter, ttest_lead$p.value))

#t = 0.017 | df = 66 | p = 0.9862

## 5.4 ICC per condition (binary GLMM)
merged_data <- merged_data %>%
  mutate(at_front = as.integer(rank_from_food == 1))

compute_icc_and_var <- function(df, cond_name) {
  mod <- glmmTMB(
    at_front ~ 1 + (1 | Chick_ID) + (1 | BlockID),
    family = binomial,
    data = df
  )
  vc <- VarCorr(mod)$cond
  var_chick <- as.numeric(vc$Chick_ID[1])^2
  var_block <- as.numeric(vc$BlockID[1])^2
  v_logit  <- (pi^2) / 3
  icc_val  <- var_chick / (var_chick + var_block + v_logit)
  
  tibble(
    Condition   = cond_name,
    ICC         = icc_val,
    Var_ChickID = var_chick,
    Var_BlockID = var_block
  )
}

icc_per_condition <- merged_data %>%
  split(.$Condition) %>%
  purrr::map_dfr(~ compute_icc_and_var(.x, unique(.x$Condition)))

print(icc_per_condition)

# A tibble: 2 × 4
#Condition   ICC Var_ChickID Var_BlockID
#<fct>       <dbl>       <dbl>       <dbl>
#1 Control   0.950        61.9    4.42e-18
#2 Novel     0.863        20.7    1.74e-17

#Possible models for AIC comparison
m0_rq3_noZI <- glmmTMB(at_front ~ 1 + (1|Chick_ID) + (1|BlockID),
                       family = binomial, data = merged_data)
m1_rq3_noZI <- glmmTMB(at_front ~ Condition + (1|Chick_ID) + (1|BlockID),
                       family = binomial, data = merged_data)
m2_rq3_noZI <- glmmTMB(at_front ~ groupsize + (1|Chick_ID) + (1|BlockID),
                       family = binomial, data = merged_data)
m3_rq3_noZI <- glmmTMB(at_front ~ Condition + groupsize + (1|Chick_ID) + (1|BlockID),
                       family = binomial, data = merged_data)
m4_rq3_noZI <- glmmTMB(at_front ~ Condition * groupsize + (1|Chick_ID) + (1|BlockID),
                       family = binomial, data = merged_data)

m0_rq3_ZI <- glmmTMB(at_front ~ 1 + (1|Chick_ID) + (1|BlockID),
                     ziformula = ~1, family = binomial, data = merged_data)
m1_rq3_ZI <- glmmTMB(at_front ~ Condition + (1|Chick_ID) + (1|BlockID),
                     ziformula = ~1, family = binomial, data = merged_data)
m2_rq3_ZI <- glmmTMB(at_front ~ groupsize + (1|Chick_ID) + (1|BlockID),
                     ziformula = ~1, family = binomial, data = merged_data)
m3_rq3_ZI <- glmmTMB(at_front ~ Condition + groupsize + (1|Chick_ID) + (1|BlockID),
                     ziformula = ~1, family = binomial, data = merged_data)
m4_rq3_ZI <- glmmTMB(at_front ~ Condition * groupsize + (1|Chick_ID) + (1|BlockID),
                     ziformula = ~1, family = binomial, data = merged_data)

aic_rq3_full <- tibble(
  model = c("m0_noZI","m1_noZI","m2_noZI","m3_noZI","m4_noZI",
            "m0_ZI","m1_ZI","m2_ZI","m3_ZI","m4_ZI"),
  AIC   = c(
    safe_AIC(m0_rq3_noZI), safe_AIC(m1_rq3_noZI),
    safe_AIC(m2_rq3_noZI), safe_AIC(m3_rq3_noZI),
    safe_AIC(m4_rq3_noZI),
    safe_AIC(m0_rq3_ZI),  safe_AIC(m1_rq3_ZI),
    safe_AIC(m2_rq3_ZI),  safe_AIC(m3_rq3_ZI),
    safe_AIC(m4_rq3_ZI)
  )
) %>% arrange(AIC)

cat("\n-- Full AIC comparison (RQ3, zero-inflated vs not) --\n")
print(aic_rq3_full)

# A tibble: 10 × 2
#model       AIC
#<chr>      <dbl>
#1 m0_ZI    37444.
#2 m2_ZI    37445.
#3 m1_ZI    37446.
#4 m3_ZI    37447.
#5 m4_ZI    37448.

#6 m0_noZI  40734.
#7 m2_noZI  40735.
#8 m1_noZI  40736.
#9 m3_noZI  40737.
#10 m4_noZI 40738.

#The response at_front was coded as 0 or 1 if a gull occupied the front position in a frame or not. 
#Because many individuals rarely or never led,a zero-inflated binomial GLMM is used to account for excess zeros. 
#Zero-inflated models seem to be way better; choosing m1_rq3_ZI here

summary(m0_rq3_ZI)

#Conditional model:
#.            Estimate Std. Error  z value Pr(>|z|)  
#(Intercept)    2.195      1.165   1.885   0.0595 .
---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Zero-inflation model:
#             Estimate Std. Error z value Pr(>|z|)    
#(Intercept)   0.2004     0.0182   11.01   <2e-16 ***
  ---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  #--> In a zero-inflated model, two processes are modelled at the same time:
  #The conditional model describes the probability that a gull is at the front in a given frame,
  #considering effects such as condition or random variation between individuals and trials.
  #The zero-inflation model captures the excess of zeros that cannot be explained by the main model alone. In this case, gulls that never lead at all.
  
#After accounting for repeated measures and zero inflation, the probability that a gull was leading did not differ between the Control and Novel conditions (p = 0.98).
#The weak correlation indicates that a bird that leads a lot in the control condition does not necessarily lead a lot in novel condition. 
#However, the ICC's per condition are high, so some birds consistently lead more than others within a condition; there strong individual differences within each condition. 
#--> Leadership thus seems to be consistent within conditions, but not across conditions. 

#Visual of per-bird mean rank across trials
rank_summary <- merged_data %>%
  group_by(Chick_ID, Condition, BlockID) %>%
  summarise(mean_rank = mean(rank_from_food, na.rm = TRUE), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(TrialNum = dense_rank(BlockID)) %>%
  ungroup() %>%
  left_join(leadership_condition %>% select(Chick_ID, Condition, prop_front),
            by = c("Chick_ID", "Condition"))

ggplot(rank_summary, aes(x = TrialNum, y = mean_rank, group = Chick_ID, color = prop_front)) +
  geom_line(linewidth = 1) +
  geom_point(aes(size = prop_front), alpha = 0.7) +
  scale_y_reverse() +
  scale_color_viridis_c(option = "plasma") +
  scale_size_continuous(range = c(0.5, 2)) +
  facet_wrap(~Condition, scales = "free_x") +
  labs(title = "Per-bird mean rank across trials",
       x = "Trial number",
       y = "Mean rank (lower = closer to front)",
       color = "Prop. of time at front",
       size = "Prop. of time at front") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45),
        panel.grid.minor = element_blank())

----------
#R4:The fourth research question explores whether less cohesive groups show more overtaking events, or whether more cohesive groups are associated with more synchronized movement and fewer position changes.
----------
  
  # Trial-level cohesion 
  trial_cohesion_rq4 <- merged_data %>%
  group_by(BlockID, Condition) %>%
  summarise(mean_logNND = mean(logNND, na.rm = TRUE), .groups = "drop")

# Merge data so cohesion, overtakes and groupsize are available
trial_summary_rq4 <- trial_cohesion_rq4 %>%
  left_join(
    overtake_trial %>%
      select(BlockID, Condition, n_overtakes, log_exposure, groupsize),
    by = c("BlockID","Condition")
  ) %>%
  drop_na(n_overtakes, log_exposure, mean_logNND, groupsize) %>%
  mutate(
    mean_logNND_scaled = as.numeric(scale(mean_logNND)),
    groupsize = factor(groupsize)
  )

stopifnot(all(is.finite(trial_summary_rq4$log_exposure)))

#Possible models for AIC comparison
m0_rq4 <- glmmTMB(n_overtakes ~ 1 + (1|BlockID),
                  family = nbinom2, offset = log_exposure, data = trial_summary_rq4)
m1_rq4 <- glmmTMB(n_overtakes ~ Condition + (1|BlockID),
                  family = nbinom2, offset = log_exposure, data = trial_summary_rq4)
m2_rq4 <- glmmTMB(n_overtakes ~ mean_logNND_scaled + (1|BlockID),
                  family = nbinom2, offset = log_exposure, data = trial_summary_rq4)
m3_rq4 <- glmmTMB(n_overtakes ~ Condition + mean_logNND_scaled + (1|BlockID),
                  family = nbinom2, offset = log_exposure, data = trial_summary_rq4)
m4_rq4 <- glmmTMB(n_overtakes ~ Condition * mean_logNND_scaled + (1|BlockID),
                  family = nbinom2, offset = log_exposure, data = trial_summary_rq4)
m5_rq4 <- glmmTMB(n_overtakes ~ Condition + groupsize + (1|BlockID),
                  family = nbinom2, offset = log_exposure, data = trial_summary_rq4)
m6_rq4 <- glmmTMB(n_overtakes ~ Condition * groupsize + (1|BlockID),
                  family = nbinom2, offset = log_exposure, data = trial_summary_rq4)
m7_rq4 <- glmmTMB(n_overtakes ~ Condition + mean_logNND_scaled + groupsize + (1|BlockID),
                  family = nbinom2, offset = log_exposure, data = trial_summary_rq4)
m8_rq4 <- glmmTMB(n_overtakes ~ Condition * mean_logNND_scaled + groupsize + (1|BlockID),
                  family = nbinom2, offset = log_exposure, data = trial_summary_rq4)

models_nb_rq4 <- list(
  m0 = m0_rq4, m1 = m1_rq4, m2 = m2_rq4, m3 = m3_rq4, m4 = m4_rq4,
  m5 = m5_rq4, m6 = m6_rq4, m7 = m7_rq4, m8 = m8_rq4
)

aic_rq4_nb <- tibble(
  model = names(models_nb_rq4),
  AIC   = vapply(models_nb_rq4, AIC, numeric(1))
) %>% arrange(AIC)

print(aic_rq4_nb)

#1 m6     415.
#2 m5     416.
#3 m7     418.
#4 m8     419.
#5 m0     425.
#6 m1     425.
#7 m2     426.
#8 m3     427.
#9 m4     428.

summary(m6_rq4)

#Conditional model:
#                           Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                 0.3433     0.1593   2.155 0.031171 *  
#ConditionNovel              0.5660     0.2146   2.638 0.008350 ** 
#groupsize5                  0.8169     0.2191   3.729 0.000192 ***
#ConditionNovel:groupsize5  -0.5420     0.3046  -1.780 0.075145 .

suppressWarnings({
  print(performance::r2(m6_rq4))
  print(performance::check_overdispersion(m6_rq4))
})

#dispersion ratio = 0.907
#p-value =  0.96

#No overdispersion detected

rhs_terms <- all.vars(formula(m6_rq4))  
has_cond   <- "Condition"   %in% rhs_terms
has_gsize  <- "groupsize"   %in% rhs_terms

#groupsize = 4:
#Condition response    SE  df asymp.LCL asymp.UCL
#Control       1.41 0.225 Inf      1.03      1.93
#Novel         2.48 0.357 Inf      1.87      3.29

#groupsize = 5:
#Condition response    SE  df asymp.LCL asymp.UCL
#Control       3.19 0.480 Inf      2.38      4.28
#Novel         3.27 0.507 Inf      2.41      4.43

#Confidence level used: 0.95 
#Intervals are back-transformed from the log scale 

#groupsize = 4:
#contrast        ratio    SE  df null z.ratio p.value
#Control / Novel 0.568 0.122 Inf    1  -2.638  0.0084

#groupsize = 5:
#contrast        ratio    SE  df null z.ratio p.value
#Control / Novel 0.976 0.211 Inf    1  -0.111  0.9115

#Tests are performed on the log scale 

rhs_terms <- all.vars(formula(m7_rq4))  
has_cond   <- "Condition"   %in% rhs_terms
has_gsize  <- "groupsize"   %in% rhs_terms

  suppressWarnings({
    print(performance::r2(m7_rq4))
    print(performance::check_overdispersion(m7_rq4))

#No overdispersion detected.

#Overtake dynamics are better explained by condition × group size than by cohesion. 
#The best-fitting model does not include mean_logNND at all, meaning cohesion ( captured by the NND) 
#is not strongly predictive of overtake counts once condition and group size are in the model.

#Visual of overtakes / second by condition across group sizes
if (has_cond && has_gsize) {
  emm_rq4 <- emmeans(m6_rq4, ~ Condition | groupsize, type = "response", offset = 0)
  print(emm_rq4)
  print(pairs(emm_rq4, by = "groupsize"))
  
  emmip(m6_rq4,
        Condition ~ groupsize,
        type   = "response",
        offset = 0) +
    labs(title = "Overtakes / second by Condition across group sizes",
         x = "Group size (birds)", y = "Overtakes / second (NB2 + offset)") +
    theme_minimal()
}






