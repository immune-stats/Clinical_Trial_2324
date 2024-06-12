#Loading the data

library(dplyr)
file_path <- '......'

data <- read.csv(file_path, header = TRUE)

#Analysis of the data
summary(data)
head(data)
colSums(is.na(data))
sapply(data[, c("treatment_arm", "sex", "CYP2D6_phenotype")], function(x) length(unique(x)))

table(data$gametocytes_day_7)

data <- data %>%
  mutate(
    treatment = factor(treatment_arm, levels = c(0.00, 0.25, 0.40), 
                       labels = c("Placebo", "PQ_0.25", "PQ_0.40")),
    safe = case_when(
      age <= 5 & hb_day_7 >= 11.0 ~ 1,
      age >= 6 & age <= 12 & hb_day_7 >= 11.5 ~ 1,
      age >= 13 & age <= 14 & hb_day_7 >= 12.0 ~ 1,
      TRUE ~ 0
    ),
    efficacy = ifelse(gametocytes_day_7 == 0, 1, 0)
  )


data$Gender <- as.numeric(ifelse(data$sex == 'F',1,0))
#Treating lost in the follow-up as non-responders to the treatments
data$efficacy[is.na(data$efficacy)] <- 0

library(ggplot2)

# Hemoglobin Levels on Day 7 by Treatment
ggplot(data, aes(x = treatment, y = hb_day_7, fill = treatment)) +
  geom_boxplot(color = "black") +
  labs(title = "Hemoglobin Levels on Day 7 by Treatment", x = "Treatment", y = "Hemoglobin (g/dL)") +
  scale_fill_brewer(palette = "Set2")

# Gametocyte Presence on Day 7 by Treatment
ggplot(data, aes(x = treatment, fill = factor(gametocytes_day_7))) +
  geom_bar(position = "dodge", color = "black") +
  labs(title = "Gametocyte Presence on Day 7 by Treatment", x = "Treatment", fill = "Gametocytes on Day 7")+
  scale_fill_brewer(palette = "Set2")

# Treatment Efficacy on Day 7
ggplot(data, aes(x = treatment, fill = factor(efficacy))) +
  geom_bar(position = "dodge", color = "black") +
  labs(title = "Treatment Efficacy on Day 7", x = "Treatment", fill = "Efficacy (1 = No Gametocytes)")+
  scale_fill_brewer(palette = "Set2")

# Safety of Treatment on Day 7
ggplot(data, aes(x = treatment, fill = factor(safe))) +
  geom_bar(position = "dodge", color = "black") +
  labs(title = "Safety of Treatment on Day 7", x = "Treatment", fill = "Safety (1 = Safe)")+
  scale_fill_brewer(palette = "Set2")

# Descriptive analysis
summary_stats <- data %>%
  group_by(treatment) %>%
  summarise(
    count = n(),
    mean_age = mean(age, na.rm = TRUE),
    mean_hb_day_7 = mean(hb_day_7, na.rm = TRUE),
    efficacy_rate = mean(efficacy, na.rm = TRUE),
    safety_rate = mean(safe, na.rm = TRUE)
  )

print(summary_stats)

## Statistical analysis

#Efficacy endpoint
table_efficacy <- table(data$treatment_arm, data$efficacy)
print(table_efficacy)

#Age
table_efficacy_age <- table(data$age, data$efficacy)
print(table_efficacy_age)

data_eff_age <- data %>%
  group_by(age, efficacy) %>%
  summarise(count = n()) %>%
  ungroup()

ggplot(data_eff_age, aes(x = age, y = count, fill = factor(efficacy))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Efficacy by Age", x = "Age", y = "Count", fill = "Efficacy") +
  scale_fill_brewer(palette = "Set2")+
  scale_x_continuous(breaks = seq(min(data$age), max(data$age), by = 1)) 

#Gender
table_efficacy_sex <- table(data$Gender, data$efficacy) #1=female
print(table_efficacy_sex)

chisq_test <- chisq.test(table_efficacy)
print(chisq_test)

data_eff_sex <- data %>%
  mutate(Gender = factor(Gender, levels = c(0, 1), labels = c("Male", "Female")))

data_eff_sex <- data_eff_sex %>%
  group_by(Gender, efficacy) %>%
  summarise(count = n()) %>%
  ungroup()

ggplot(data_eff_sex, aes(x = Gender, y = count, fill = factor(efficacy))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Efficacy by Gender", x = "Gender", y = "Count", fill = "Efficacy") +
  scale_fill_brewer(palette = "Set2")


#Pairwise chi-squared
chi_square_pairwise <- function(data, group1, group2) {
  sub_data <- data %>% filter(treatment_arm %in% c(group1, group2))
  table_sub <- table(sub_data$treatment_arm, sub_data$efficacy)
  test_result <- chisq.test(table_sub)
  return(test_result$p.value)
}

groups <- unique(data$treatment_arm)

p_values <- matrix(nrow = length(groups), ncol = length(groups))
rownames(p_values) <- colnames(p_values) <- groups


for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    p_values[i, j] <- chi_square_pairwise(data, groups[i], groups[j])
    p_values[j, i] <- p_values[i, j]  
  }
}

# Bonferroni correction
p_adjusted <- p.adjust(p_values[upper.tri(p_values)], method = "bonferroni")

p_values_adj <- matrix(NA, nrow = length(groups), ncol = length(groups))
rownames(p_values_adj) <- colnames(p_values_adj) <- groups
p_values_adj[upper.tri(p_values_adj)] <- p_adjusted


print("P-value matrix before the correction:")
print(p_values)
print("P-value matrix after the correction:")
print(p_values_adj)


# Safety endpoint

table_safe <- table(data$treatment_arm, data$safe)
print(table_safe)

#Age
table_safety_age <- table(data$age, data$safe)
print(table_safety_age)

data_safe_age <- data %>%
  group_by(age, safe) %>%
  summarise(count = n()) %>%
  ungroup()

ggplot(data_safe_age, aes(x = age, y = count, fill = factor(safe))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Safety by Age", x = "Age", y = "Count", fill = "Safety") +
  scale_fill_brewer(palette = "Set2")+
  scale_x_continuous(breaks = seq(min(data$age), max(data$age), by = 1)) 

#Gender
table_safety_sex <- table(data$Gender, data$safe) #1=female
print(table_safety_sex)

data_safe_sex <- data %>%
  mutate(Gender = factor(Gender, levels = c(0, 1), labels = c("Male", "Female")))

data_safe_sex <- data_safe_sex %>%
  group_by(Gender, safe) %>%
  summarise(count = n()) %>%
  ungroup()

ggplot(data_safe_sex, aes(x = Gender, y = count, fill = factor(safe))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Safety by Gender", x = "Gender", y = "Count", fill = "Safety") +
  scale_fill_brewer(palette = "Set2")


bartlett.test(hb_day_7 ~ treatment_arm, data = data)
model_hb_day_7_treatment <- aov(hb_day_7 ~ treatment_arm, data = data)
table_hb_day_7_treatment <- anova(model_hb_day_7_treatment)
table_hb_day_7_treatment$"Pr(>F)"[1]
model_hb_day_7 <- kruskal.test(hb_day_7 ~ treatment_arm, data = data)
shapiro.test(rstandard(model_hb_day_7_treatment))

chisq_test <- chisq.test(table_safe)
print(chisq_test)

## AGE AND GENDER

glm_efficacy <- glm(efficacy ~ treatment_arm + age + sex, data = data, family = binomial(link=
                                                                                           "logit"))
summary(glm_efficacy)

# Fit GLM for Safety
glm_safety <- glm(safe ~ treatment_arm + age + sex, data = data, family = binomial(link="logit"))
summary(glm_safety)

ggplot(data, aes(x = age)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black", alpha = 0.7) +
  facet_wrap(~ treatment, scales = "free_y") +
  labs(title = "Age Distribution by Treatment Arm",
       x = "Age",
       y = "Frequency") +
  theme_minimal()


### CYP2D6 phenotype
#Efficacy
table_efficacy_cyp2d6 <- table(addNA(data$CYP2D6_phenotype), addNA(data$efficacy))
print(table_efficacy_cyp2d6)

data_efficacy_cyp2d6 <- data %>%
  mutate(CYP2D6_phenotype = factor(CYP2D6_phenotype, 
                                   levels = c(0, 1), 
                                   labels = c("Slow/Intermediate Metabolizer", "Extensive/Ultrarapid Metabolizer")))

data_efficacy_cyp2d6 <- data_efficacy_cyp2d6 %>%
  group_by(CYP2D6_phenotype, efficacy) %>%
  summarise(count = n()) %>%
  ungroup()

ggplot(data_efficacy_cyp2d6, aes(x = CYP2D6_phenotype, y = count, fill = factor(efficacy))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Efficacy by CYP2D6 Phenotype", x = "CYP2D6 Phenotype", y = "Count", fill = "Efficacy") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal()

#Model
glm_efficacy_CYP <- glm(efficacy ~ treatment_arm + age + sex + CYP2D6_phenotype, data = data, family = binomial(link = 'logit'))
summary(glm_efficacy_CYP)

table_CYP_efficacy <- table(data$CYP2D6_phenotype, data$efficacy)
print(table_CYP_efficacy)

chisq.test(table_CYP_efficacy)


#Safety
table_safety_cyp2d6 <- table(addNA(data$CYP2D6_phenotype), addNA(data$safe))
print(table_safety_cyp2d6)

data_safety_cyp2d6 <- data %>%
  mutate(CYP2D6_phenotype = factor(CYP2D6_phenotype, 
                                   levels = c(0, 1), 
                                   labels = c("Slow/Intermediate Metabolizer", "Extensive/Ultrarapid Metabolizer")))

data_safety_cyp2d6 <- data_safety_cyp2d6 %>%
  group_by(CYP2D6_phenotype, safe) %>%
  summarise(count = n()) %>%
  ungroup()

ggplot(data_safety_cyp2d6, aes(x = CYP2D6_phenotype, y = count, fill = factor(safe))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Safety by CYP2D6 Phenotype", x = "CYP2D6 Phenotype", y = "Count", fill = "Safety") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal()

glm_safety_CYP <- glm(safe ~ treatment_arm + age + sex+CYP2D6_phenotype , data = data, family = binomial(link = 'logit'))
summary(glm_safety)

table_CYP_safe <- table(data$CYP2D6_phenotype, data$safe)
print(table_CYP_safe)

chisq.test(table_CYP_safe)

# Stratified analysis
#Without phenotype 

#Efficacy
stratified_analysis <- by(data, data$sex, function(subdata) {
  lm_model <- glm(efficacy ~ treatment_arm, data = subdata, family = binomial(link="logit"))
  return(summary(lm_model))
})
stratified_analysis

stratified_analysis_2 <- by(data, data$age, function(subdata) {
  lm_model <- glm(efficacy ~ treatment_arm, data = subdata, family = binomial(link="logit"))
  return(summary(lm_model))
}) #! treatment_arm age: 3, 4, 11
stratified_analysis_2

#Safety
stratified_analysis <- by(data, data$sex, function(subdata) {
  lm_model <- glm(safe ~ treatment_arm, data = subdata, family = binomial(link="logit"))
  return(summary(lm_model))
})
stratified_analysis

stratified_analysis_2 <- by(data, data$age, function(subdata) {
  lm_model <- glm(safe ~ treatment_arm, data = subdata, family = binomial(link="logit"))
  return(summary(lm_model))
}) #treatment_arm, age: 9
stratified_analysis_2

#With phenotype
#Efficacy
stratified_analysis <- by(data, data$sex, function(subdata) {
  lm_model <- glm(efficacy ~ treatment_arm + CYP2D6_phenotype, data = subdata, family = binomial(link="logit"))
  return(summary(lm_model))
})#! treatment_arm, phenotype, sex: F
stratified_analysis

stratified_analysis_2 <- by(data, data$age, function(subdata) {
  lm_model <- glm(efficacy ~ treatment_arm + CYP2D6_phenotype, data = subdata, family = binomial(link="logit"))
  return(summary(lm_model))
})
stratified_analysis_2

#Safety
stratified_analysis <- by(data, data$sex, function(subdata) {
  lm_model <- glm(safe ~ treatment_arm + CYP2D6_phenotype, data = subdata, family = binomial(link="logit"))
  return(summary(lm_model))
}) #! phenotype, sex: M
stratified_analysis

stratified_analysis_2 <- by(data, data$age, function(subdata) {
  lm_model <- glm(safe ~ treatment_arm + CYP2D6_phenotype, data = subdata, family = binomial(link="logit"))
  return(summary(lm_model))
})
stratified_analysis_2
#! phenotype, age: 2