# ==============================================================================
# FECAL INDICATOR BACTERIA & MICROBIAL SOURCE TRACKING (FIB-MST) PIPELINE
# ==============================================================================
# Purpose: Analyze fecal indicator bacteria (E. coli) concentrations and microbial source tracking (MST)
# marker amplifications across Austin-area creeks to identify pollution sources and high-risk sites
#
# MST Markers:
# - HF183: Human-associated Bacteroides marker
# - DogBac: Dog-associated Bacteroides marker
# High confidence detection defined as >= 3 amplifications out of 4
#
# Data: final_data_lexi_FIBMST_paper_2025_GOOD.csv
# EPA Benchmark: Geometric mean of 126 cfu/100mL; single-sample max 406 cfu/100mL


# ------------------------------------------------------------------------------
# 1. LOAD PACKAGES
# ------------------------------------------------------------------------------
library(tidyverse)
library(FSA)

# ------------------------------------------------------------------------------
# 2. CONSTANTS
# ------------------------------------------------------------------------------
EPA_GM_THRESHOLD <- 126 

# Colorblind-friendly palette used throughout
MST_COLORS <- c(
  "HF183"             = "#0072B2",
  "HF183_copynumber"  = "#0072B2",
  "DogBac"            = "#E69F00",
  "DogBac_copynumber" = "#E69F00",
  "Both"              = "#CC79A7",
  "Ecoli_count"       = "#CC79A7"
)

CREEK_COLORS <- c(
  "#4A4A4A", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

# ------------------------------------------------------------------------------
# 3. LOAD AND CLEAN MAIN DATA
# ------------------------------------------------------------------------------
# - Trim whitespace from Creek names
# - Remove Waller 2017 samples (unreliable)
# - Strip "/4" denominators from amplification columns and coerce to numeric
# - Recode "N/A" and "N" entries in amplifications/CQ columns to 0

setwd("~/Team Stuart Data Analysis/FIB-MST")

fibmst_raw <- read.csv("data/final data (lexi) - FIBMST paper 2025 - GOOD.csv")

clean_fibmst_data <- fibmst_raw %>%
  mutate(Date = mdy(Date),
         Creek = trimws(Creek)) %>% 
  filter(!(Creek == "Waller" & year(Date) == 2017)) %>% 
  mutate(X..of.HF183.amplifications..out.of.4. = str_replace(X..of.HF183.amplifications..out.of.4., "/.*", ""),
         X..of.DogBac.amplifications..out.of.4. = str_replace(X..of.DogBac.amplifications..out.of.4., "/.*", ""),
         HF183.CQ.mean = str_replace(HF183.CQ.mean, "N/A", "0"),
         DogBac.CQ.mean = str_replace(DogBac.CQ.mean, "N/A", "0")) %>% 
  filter(X..of.HF183.amplifications..out.of.4. != "" & X..of.DogBac.amplifications..out.of.4. != "" &
           E..coli.count..avg.cfu.or.mpn.100mL. != "")

# ------------------------------------------------------------------------------
# 4. MODEL ASSUMPTION CHECKS
# ------------------------------------------------------------------------------
# Verify whether standard linear regression assumptions (linearity, homoscedasticity,
# normality) are met. If not, assess whether log-transformation resolves violations

fib_model_data <- clean_fibmst_data %>% 
  select(Creek,
         Ecoli_count = E..coli.count..avg.cfu.or.mpn.100mL.,
         HF183_amps = X..of.HF183.amplifications..out.of.4.,
         DogBac_amps = X..of.DogBac.amplifications..out.of.4.) %>% 
  filter(!is.na(as.numeric(Ecoli_count)),
         !is.na(as.numeric(HF183_amps)),
         !is.na(as.numeric(DogBac_amps))) %>% 
  mutate(Ecoli_count = as.numeric(Ecoli_count),
         HF183_amps = as.numeric(HF183_amps),
         DogBac_amps = as.numeric(DogBac_amps),
         log_Ecoli = log(Ecoli_count + 1))
  
# --- Linearity: raw E. coli vs. amplifications ---
fib_model_data %>%
  ggplot(aes(x = HF183_amps, y = Ecoli_count)) +
  geom_point() +
  labs(title = "E. coli Count vs. HF183 Amplifications (Raw)") +
  theme_minimal()
# Linearity not met, right-skewed

fib_model_data %>% 
  ggplot(aes(x = DogBac_amps, y = Ecoli_count)) +
  geom_point() +
  labs(title = "E. coli Count vs. DogBac Amplifications (Raw)") +
  theme_minimal()
# Linearity not met

# --- Homoscedasticity and normality: raw vs. log-transformed ---
plot(lm(Ecoli_count ~ factor(HF183_amps), data = fib_model_data))
# Homoscedasticity not met, normality not met
plot(lm(Ecoli_count ~ factor(DogBac_amps), data = fib_model_data))
# Homoscedasticity not met, normality not met

plot(lm(log_Ecoli ~ factor(HF183_amps), data = fib_model_data))
plot(lm(log_Ecoli ~ factor(DogBac_amps), data = fib_model_data))

# Conclusion: log(E.coli + 1) transformation is appropriate for all models

# ------------------------------------------------------------------------------
# 5. AMPLIFICATION TOTALS SUMMARY
# ------------------------------------------------------------------------------
# Summarize how many samples met the high-confidence threshold (>= 3 amps)
sample_total <- nrow(fib_model_data)
human_total <- sum(fib_model_data$HF183_amps >= 3)
dog_total <- sum(fib_model_data$DogBac_amps >= 3)

totals_df <- data.frame(
  Category = c("Total Samples Analyzed", "Human Source (>=3 Amps)", "Dog Sources (>=3 Amps)"),
  Count = c(sample_total, human_total, dog_total)
)

ggplot(totals_df, aes(x = Category, y = Count)) + 
  geom_bar(stat = "identity", width = 0.75) +
    geom_text(aes(label = Count), vjust = -0.5, size = 3.5) +
      labs(title = "Fecal Source Tracking Results: High Confidence Amplifications",
           subtitle = "High confidence defined as >=3 amplifications out of 4",
           x = NULL, y = "Count of Samples") +
        theme_minimal(base_size = 13) +
          theme(plot.title = element_text(face = "bold"),
                axis.text.x = element_text(angle = 15, hjust = 1)
          )

# ------------------------------------------------------------------------------
# 6. MULTIPLE LINEAR REGRESSION
# ------------------------------------------------------------------------------
# Test whether human and dog MST markers jointly predict log E. coli
fib_mlr <- lm(log_Ecoli ~ HF183_amps + DogBac_amps + HF183_amps*DogBac_amps, data = fib_model_data) 
summary(fib_mlr)
# Only HF183 is a significant predictor of log E. coli, DogBac and the interaction term aren't

# ------------------------------------------------------------------------------
# 7. ONE-WAY ANOVA
# ------------------------------------------------------------------------------
# Treats HF183 amplification count as a categorical grouping variable to test whether mean log E. coli
# differs across amplification levels
human_anova <- aov(lm(log_Ecoli~ factor(HF183_amps), data = fib_model_data))
summary(human_anova) 
# There's a significant difference in means across groups

# --- Tukey's HSD post-hoc ---
# To see which groups differ
TukeyHSD(human_anova)
# Groups 2 vs. 0 and 4 vs. 0 are significantly different

# --- Boxplot: log E. coli by HF183 amplification group ---
fib_model_data %>% 
  ggplot(aes(x = factor(HF183_amps), y = log_Ecoli)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(alpha = 0.5) + 
  labs(title = "Log E. coli Concentrations by HF183 Amplification Group",
       subtitle = "Analysis on log(E.coli + 1) transformed data",
       x = "Number of Amplifications (0-4)",
       y = "Log(E. coli + 1) (cfu/100mL)") +
  theme_minimal()
  
# ------------------------------------------------------------------------------
# 8. HELPER FUNCTION: CREEK PROFILE PLOT
# ------------------------------------------------------------------------------
# All per-creek plots follow an identical structure. This function eliminates the repetition:
# filter to one creek and date, pivot to long format, compute a dual-axis scale factor, and plot
make_creek_profile <- function(data, creek_name, sample_date, title) {
  
  creek_long <- data %>%
    filter(Creek == creek_name, Date == as.Date(sample_date)) %>%
    select(Distance_frommouth, Ecoli_count, HF183_copynumber, DogBac_copynumber) %>%
    pivot_longer(cols      = -Distance_frommouth,
                 names_to  = "Type",
                 values_to = "Value") %>%
    arrange(desc(Distance_frommouth))
  
  max_ecoli    <- max(creek_long$Value[creek_long$Type == "Ecoli_count"], na.rm = TRUE)
  max_mst      <- max(creek_long$Value[creek_long$Type %in%
                                         c("HF183_copynumber", "DogBac_copynumber")], na.rm = TRUE)
  scale_factor <- max_mst / max_ecoli
  
  creek_long <- creek_long %>%
    mutate(
      Value_scaled = ifelse(Type == "Ecoli_count", Value * scale_factor, Value),
      Type = factor(Type, levels = c("HF183_copynumber", "DogBac_copynumber", "Ecoli_count"))
    )
  
  ggplot(creek_long, aes(as.factor(Distance_frommouth), Value_scaled, fill = Type)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.7) +
    scale_y_continuous(
      name     = "MST (Copy Number/100mL)",
      sec.axis = sec_axis(~ . / scale_factor, name = "FIB (MPN/100mL)")
    ) +
    scale_fill_manual(
      values = MST_COLORS[c("HF183_copynumber", "DogBac_copynumber", "Ecoli_count")],
      labels = c("HF183_copynumber"  = "HF-183",
                 "DogBac_copynumber" = "DogBac",
                 "Ecoli_count"       = "FIB")
    ) +
    labs(title = title,
         x     = "Distance from Mouth of Creek (m)") +
    theme_minimal() +
    theme(plot.title         = element_text(face = "bold", hjust = 0.5),
          legend.position    = "bottom",
          legend.title       = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
}

# ------------------------------------------------------------------------------
# 9. PIE CHART: FECAL SOURCE TYPE AMONG POSITIVE SAMPLES
# ------------------------------------------------------------------------------
# A sample is considered positive if it meets BOTH criteria:
#   1. >= 3 amplifications out of 4
#   2. CQ mean < 36.9

piechart_data <- clean_fibmst_data %>%
  select(X..of.HF183.amplifications..out.of.4., HF183.CQ.mean,
         X..of.DogBac.amplifications..out.of.4., DogBac.CQ.mean) %>%
  rename(HF183_amps  = X..of.HF183.amplifications..out.of.4.,
         DogBac_amps = X..of.DogBac.amplifications..out.of.4.) %>%
  mutate(HF183_amps  = as.numeric(HF183_amps),
         DogBac_amps = as.numeric(DogBac_amps)) %>% 
  filter(HF183_amps >= 3 | DogBac_amps >= 3,
         HF183.CQ.mean < 36.9, DogBac.CQ.mean < 36.9,
         HF183.CQ.mean != "", DogBac.CQ.mean != "") %>% 
  drop_na(HF183_amps, DogBac_amps)

n_positive <- nrow(piechart_data)

piechart_df <- data.frame(
  Type = c("HF-183", "DogBac", "Both"),
  Perc = c(
    mean(piechart_data$HF183_amps >= 3 & piechart_data$DogBac_amps == 0),
    mean(piechart_data$DogBac_amps == 0 & piechart_data$HF183_amps >= 3),
    mean(piechart_data$DogBac_amps >= 3 & piechart_data$HF183_amps >= 3)
  )
) %>%
  arrange(Perc) %>%
  mutate(Amps = scales::percent(Perc))

ggplot(piechart_df, aes("", Perc, fill = Type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste(Type, "=", Amps)),
            position    = position_stack(vjust = 0.5),
            show.legend = FALSE, size = 4.5) +
  scale_fill_manual(
    values = c("HF-183" = "#0072B2", "DogBac" = "#E69F00", "Both" = "#CC79A7"),
    labels = c("HF-183" = "Only HF-183 Amplification",
               "DogBac" = "Only DogBac Amplification",
               "Both"   = "Both HF-183 and DogBac")
  ) +
  labs(title = paste0("Distribution of Source Types Among Positive Samples (n=", n_positive, ")")) +
  theme_void() +
  theme(plot.title      = element_text(size = 15, face = "bold", hjust = 0.5),
        legend.title    = element_blank(),
        legend.position = "bottom")


