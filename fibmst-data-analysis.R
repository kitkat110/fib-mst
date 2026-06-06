library(tidyverse)

setwd("~/Desktop/Team Stuart Data Analysis/FIB-MST")

fib_data_raw <- Copy_of_UrbEco_FIB_MST_Data_2021_Present_Good_correct_data

# See column names
names(fib_data_raw)

# CLEANING
fib_data_clean <- fib_data_raw %>%
  select(Date, # Select relevant columns
         Creek,
         Ecoli_count = `E. coli count (avg cfu or mpn/100mL)`,
         HF183_amps_raw = `# of HF183 amplifications (out of 4)`,
         DogBac_amps_raw = `# of DogBac amplifications (out of 4)`) %>% 
  mutate(HF183_amps = as.numeric(str_replace(HF183_amps_raw, "/.*", "")), # Clean and convert amplification columns
         DogBac_amps = recode(DogBac_amps_raw,
                              "N/A" = NA_character_,
                              "0/4*" = "0"),
         DogBac_amps = as.numeric(str_replace(DogBac_amps_raw, "/.*", "")))  %>% 
  mutate(Date = mdy(Date), 
         Month = month(Date, label = TRUE, abbr = TRUE),
         Season = case_when( # Create season variable
           Month %in% c("Dec", "Jan", "Feb") ~ "Winter",
           Month%in% c("Mar", "Apr", "May") ~ "Spring",
           Month %in% c("Jun", "Jul", "Aug") ~ "Summer",
           Month %in% c("Sep", "Oct", "Nov") ~ "Fall")
         ) %>% 
  drop_na(Ecoli_count, HF183_amps, DogBac_amps) %>% 
  mutate(Ecoli_count = as.numeric(Ecoli_count)) %>% 
  select(Date, Month, Season, Creek, Ecoli_count, HF183_amps, DogBac_amps)

# See column names
names(fib_data_clean)

# TOTALS
sample_total <- nrow(fib_data_clean)

dog_total <- fib_data_clean %>% 
  filter(DogBac_amps >= 3) %>% 
  nrow()

human_total <- fib_data_clean%>% 
  filter(HF183_amps >= 3) %>% 
  nrow()

totals_df <- data.frame(
  Category = c("Total Samples Analyzed", "Human Source (>=3 Amps)", "Dog Sources (>= Amps)"),
  Count = c(sample_total, human_total, dog_total)
)

# PLOT AMPLIFICATION TOTALS
ggplot(totals_df, aes(x=Category, y=Count)) + 
  geom_bar(stat = "identity", width=0.75, fill="#047857") +
    geom_text(aes(label=Count), vjust=-0.5, size=3.5) +
      labs(title = "Fecal Source Tracking Results: High Confidence Amplifications",
           subtitle = "High confidence defined as 3 or more amplifications (out of 4) in the relevant dataset",
           x = NULL,
           y = "Count of Samples") +
        theme_minimal(base_size = 13) +
          theme(
            plot.title = element_text(face = "bold"),
            axis.text.x = element_text(angle = 15, hjust = 1)
          )

# REGRESSION MODEL OF HF183, DOGBAC, & E. COLI
# Assumptions tests
fib_data_clean %>% # Linearity not met
  ggplot(aes(x = HF183_amps, y = Ecoli_count)) +
  geom_point() +
  labs(title = "E. coli Count vs. HF183 Amplifications",
       x = "HF183 Amplifications (0-4)",
       y = "E. coli Concentration (avg cfu/100mL)")

fib_data_clean %>% # Linearity not met
  ggplot(aes(x = DogBac_amps, y = Ecoli_count)) +
  geom_point() +
  labs(title = "E. coli Count vs. DogBac Amplifications",
       x = "DogBac Amplifications (0-4)",
       y = "E. coli Concentration (avg cfu/100mL)")

plot(lm(Ecoli_count ~ factor(HF183_amps), data = fib_data_clean)) # Homoscedasticity not met (better at predicting low E. coli concentrations), normality not met
plot(lm(Ecoli_count ~ factor(DogBac_amps), data = fib_data_clean)) # Homoscedasticity not met, normality not met
# Standard linear regression not appropriate, GLM or ANOVA better

# E.coli data has to be transformed to compress extreme values
fib_data_clean$log_Ecoli_count <- log(fib_data_clean$Ecoli_count + 1) 

plot(lm(log_Ecoli_count ~ factor(HF183_amps), data = fib_data_clean)) # Variance mostly good, few outliers still, assumptions met
plot(lm(log_Ecoli_count ~ factor(DogBac_amps), data = fib_data_clean)) # Few outliers but assumptions met

# MODEL & PLOTS
# Multiple Linear Regression with Interaction
fib_glm_log <- fib_data_clean %>% 
  lm(log_Ecoli_count ~ HF183_amps + DogBac_amps + HF183_amps*DogBac_amps, data = .) 
summary(fib_glm_log)
# Only Human FIB is a significant predictor of E. coli concentrations

# One-way ANOVA through Linear Model 
human_anv_log <- aov(lm(log_Ecoli_count ~ factor(HF183_amps), data = fib_data_clean))
summary(human_anv_log) 
# There's a difference in means in the groups

# Tukey's HSD to see which amplification groups are different from each other
human_tukey_test <- TukeyHSD(human_anv_log)
# 2-0 and 4-0 are significantly dif.

# Boxplot of transformed E. coli data by HF183 Amplifications
human_graph_log <- fib_data_clean %>% 
  ggplot(aes(x = factor(HF183_amps), y = log_Ecoli_count)) +
  geom_boxplot(fill = "steelblue", alpha = 0.8) +
  labs(title = "Log E. coli Concentrations by HF183 Amplification",
       subtitle = "Analysis on log(E.coli + 1) transformed, unfiltered data",
       x = "Number of Amplifications (HF183)",
       y = "Log(E. coli Concentrations + 1)") +
  theme_minimal()
  
human_graph_log
# Jump from 1 to 2 amplifications, but not much change after that

# Non-parametric test - Kruskal-Wallis Rank Sum Test
kruskal.test(log_Ecoli_count ~ factor(HF183_amps), data = fib_data_clean)
# There are significant differences in log E.coli counts between the amplification groups

# Non-parametric - Dunn's Test to see which amplification groups are different from each other
install.packages("FSA")
library(FSA)

dunnTest(log_Ecoli_count ~ factor(HF183_amps), data = fib_data_clean, method = "bonferroni")
# 0-2 and 0-4 are significantly dif.



# BAR GRAPH OF AVERAGE E.COLI COUNT BY CREEK 
# Locations with chronically high FIB lvls & chronically low FIB lvls
# EPA recommendation - Geometric mean (GM) of 126 cfu/100mL, shouldn't be more than 406 cfu/100mL in any 1 sample

# Average FIB lvls
shoal_avg <- relevant_data %>% 
  filter(Creek == "Shoal")  %>% 
  summarize(geo_mean = exp(mean(log(`E. coli count (avg cfu or mpn/100mL)`)))) %>% 
  pull(geo_mean) %>% 
  print() # 260.80 cfu/100mL

boggy_avg <- relevant_data %>% 
  filter(Creek == "Boggy", `E. coli count (avg cfu or mpn/100mL)` > 0)  %>% 
  summarize(geo_mean = exp(mean(log(`E. coli count (avg cfu or mpn/100mL)`)))) %>% 
  pull(geo_mean) %>% 
  print() # 260.25 cfu/100mL

tannehill_avg <- relevant_data %>% 
  filter(Creek == "Tannehill")  %>% 
  summarize(geo_mean = exp(mean(log(`E. coli count (avg cfu or mpn/100mL)`)))) %>% 
  pull(geo_mean) %>% 
  print() # 329.35 cfu/100mL

blunn_avg <- relevant_data %>% 
  filter(Creek == "Blunn")  %>% 
  summarize(geo_mean = exp(mean(log(`E. coli count (avg cfu or mpn/100mL)`)))) %>% 
  pull(geo_mean) %>% 
  print() # 313.49 cfu/100mL

# Chronically high FIB lvls - consistently had >406 cfu/100mL in any single sample
waller_avg <- relevant_data %>% 
  filter(Creek == "Waller", `E. coli count (avg cfu or mpn/100mL)` > 0)  %>% 
  summarize(geo_mean = exp(mean(log(`E. coli count (avg cfu or mpn/100mL)`)))) %>% 
  pull(geo_mean) %>% 
  print() # 1148.33 cfu/100mL

# Chronically low FIB lvls - consistently had <406 cfu/100mL in any single sample
onion_avg <- relevant_data %>% 
  filter(Creek == "Onion")  %>% 
  summarize(geo_mean = exp(mean(log(`E. coli count (avg cfu or mpn/100mL)`)))) %>% 
  pull(geo_mean) %>% 
  print() # 27.00 cfu/100mL

bull_avg <- relevant_data %>% 
  filter(Creek == "Bull", `E. coli count (avg cfu or mpn/100mL)` > 0)  %>% 
  summarize(geo_mean = exp(mean(log(`E. coli count (avg cfu or mpn/100mL)`)))) %>% 
  pull(geo_mean) %>% 
  print() # 26.38 cfu/100 mL

creek_avg_df <- data.frame(
  Creek = c("Shoal", "Boggy", "Tannehill", "Blunn", "Waller", "Onion", "Bull"),
  Ecoli_avg = c(shoal_avg, boggy_avg, tannehill_avg, blunn_avg, waller_avg, onion_avg, bull_avg)
)
print(creek_avg_df)

# Graph of avg. E.coli counts of creeks
creek_avg_df %>% 
  ggplot( aes(x = Creek, y = Ecoli_avg, fill = Creek)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 126, color = "black", linetype = "solid", size = 0.5) +
  annotate("text", x = 3.5, y = 145, label = "EPA recommendation - 126 cfu/100mL", color = "black", size = 4) +
  annotate("text", x = 3.5, y = 45, label = "Chronically low FIB lvls.", color = "black", size = 3) +
  annotate("text", x = 7, y = 1170, label = "Chronically high FIB lvls.", color = "black", size = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = expression("Average " * italic("E.coli") * " Count by Creek"),
       subtitle = "Compared to EPA guidelines",
       x = "Creek",
       y = expression("Geometric Mean of " * italic("E.coli") * " Count (cfu/100mL)")) +
  theme_minimal()

relevant_data %>% 
  ggplot(aes(x = Creek, y = `E. coli count (avg cfu or mpn/100mL)`, fill = Creek)) +
  geom_boxplot() +
  labs(title = expression("Average " * italic("E.coli") * " Count by Creek"),
       x = "Creek",
       y = expression("Geometric Mean of " * italic("E.coli") * " Count (cfu/100mL)")) +
  theme_minimal()






