# install.packages("dplyr")
library(dplyr)
# install.packages("skimr")
library(skimr)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("GGally")
library(GGally)
# install.packages("VIM")
library(VIM)
# install.packages("car")
library(car)
# install.packages("glmnet")
library(glmnet)


set.seed(123)
df <- read.csv("Data_Clinical_Trial_Fluge_2015_for_CT_course.csv")

#We exclude variables that are more of notes rather than explanatory 
#variables for the model.
df <- df %>% select(-ID, -Exclusion_Treatment, -Reasons_Exclusion, -Reasons_Discontinued)

df <- df %>% mutate(across(where(is.character), as.factor)) %>% 
  mutate(across(where(~!is.numeric(.)), as.factor))

num_cols <- sapply(df, is.numeric)

#We define a binary variable that carries infomation whether there was a response
# in the treatment (at least 6 consecutive weeks of 4.5 fatigue score) and prepare data
# to analyze the response duration
Response <- ifelse(df$Response_duration >= 6, 1, 0)
X <- df

Response_dur <- df$Response_duration
X2 <- subset(df, select=-c(Response_duration))

df$Response <- Response

#-------------------------------------------------------------------------------------
# install.packages("table1")
library(table1)
todrop <- c("Inclusion", "Received_Treatment", "Analyzed")

?table1
table1(~. | factor(Response), data=df[-which(is.na(df$Response)), which(!(names(df) %in% todrop))])
table1(~Age+Gender | factor(Response), data=df[-which(is.na(df$Response)), which(!(names(df) %in% todrop))])
fisher.test(matrix(c(7, 2, 13, 6), nrow=2))


#-------------------------------------------------------------------------------------
#Histograms of numeric variables
for (col in names(df)[num_cols]) {
  print(ggplot(df, aes_string(x = col)) +
          geom_histogram(aes(y = ..count.., fill = is.na(get(col))), binwidth = 10, color = "black") +
          scale_fill_manual(values = c("royalblue", "black"), labels = c("Non-NA", "NA")) +
          labs(title = paste("Histogram of", col), x = col, y = "Frequency") +
          theme(legend.title = element_blank(),
                plot.title = element_text(hjust = 0.5)))
}
cat_cols <- sapply(df, is.factor)

# Bar charts for each categorical column
for (col in names(df)[cat_cols]) {
  print(ggplot(df, aes_string(x = col)) +
          geom_bar(fill = "royalblue", color = "black") +
          labs(title = paste("Bar Chart of", col), x = col, y = "Count") +
          theme(plot.title = element_text(hjust = 0.5)))
}

#Numerical variables vs response in treatment
for (col in names(df)[num_cols]) {
  print(ggplot(df, aes_string(x = "as.factor(Response)", y = col)) +
          geom_boxplot(aes(fill = as.factor(is.na(df[[col]])))) +
          scale_fill_manual(values = c("royalblue", "black"), labels = c("Non-NA", "NA")) +
          labs(title = paste("Boxplot of", col, "by Response"), x = "Response", y = col) +
          theme(legend.title = element_blank(),
                plot.title = element_text(hjust = 0.5)))  # Centering the title
}

ggplot(df, aes(x=Age, fill=factor(Response)))+
  geom_histogram(color='black',binwidth=10) +
  scale_fill_manual(values=c('cornflowerblue', 'royalblue'), labels = c("No", "Yes"))+
  labs(title = paste("Histogram of", "Age", "by Response"), x = "Age", y = "count")+
  guides(fill=guide_legend(title="Response"))+ 
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))


ggplot(df, aes(y=Age, fill=factor(Response)))+
  geom_boxplot() +
  scale_fill_manual(values=c('cornflowerblue', 'royalblue', 'black'), labels = c("No", "Yes", "NA"))+
  labs(title = paste("Boxplot of", "Age", "by Response"), y = "Age")+
  guides(fill=guide_legend(title="Response"))+ 
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

ggplot(df, aes(x=Gender, fill=factor(Response))) +
  geom_bar(col="black",stat="count", position="dodge", width=0.7)+
  scale_fill_manual(values=c('cornflowerblue', 'royalblue', 'black'), labels = c("No", "Yes", "NA"))+
  labs(title = "Response count by Gender", y = "count")+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))


ggplot(df, aes(x=Response, fill=Family_AD	)) +
  geom_bar(col="black",stat="count", width=0.7)+
  scale_fill_manual(values=c('cornflowerblue', 'royalblue', 'black'))+
  labs(title = paste("Boxplot of", "Family Autoimmune Diseases", "by Response"), y = "count")+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

ggplot(df, aes(x=Response, fill=Discontinued_Treatment	)) +
  geom_bar(col="black",stat="count", width=0.7)+
  scale_fill_manual(values=c('cornflowerblue', 'royalblue', 'black'))+
  labs(title = paste("Boxplot of", "Discontinued Treatment", "by Response"), y = "count")+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

ggplot(df, aes(x=Severity, fill=factor(Response))) +
  geom_bar(col="black",stat="count", position="dodge", width=0.7)+
  scale_fill_manual(values=c('cornflowerblue', 'royalblue', 'black'), labels = c("No", "Yes", "NA"))+
  labs(title = "Response count by Severity", y = "count")+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

#We drop all subjects who didn't begin treatment.
rows_to_remove <- which(is.na(df$Response_duration) | df$Analyzed != "Yes")

# Remove rows from matrix X
X <- X[-rows_to_remove, ]
X2 <- X2[-rows_to_remove,]
# Remove corresponding elements from vector Response and Response_dur
Response <- Response[-rows_to_remove]
Response_dur <- Response_dur[-rows_to_remove]


#We exclude the variables [Inclusion, Received_Treatment, Analyzed]
#since they all have the same value for patients who
#have undergone the treatment
X <- subset(X, select = -c(Inclusion, Received_Treatment, Analyzed, Response_duration, Response_end, Longest_period))
X2 <- subset(X2, select = -c(Inclusion, Received_Treatment, Analyzed))
colnames(X)
colnames(X2)
#Models
# checking influence on response with and without secondary endpoints
X_scaled <- X %>%
  mutate(across(where(is.numeric), scale))

X_scaled <- makeX(X_scaled, na.impute =TRUE)

cv_fit <- cv.glmnet(X_scaled, Response, alpha = 1)
lambda_min <- cv_fit$lambda.min

lasso_model_scaled <- glmnet(X_scaled, Response, family="binomial", alpha=1, lambda=lambda_min)
coefficients_scaled <- coef(lasso_model_scaled)
print(coefficients_scaled)

# without secondary endpoints

X_without_sec <- X[,1:8]

X_scaled_ws <- X_without_sec %>%
  mutate(across(where(is.numeric), scale))

X_scaled_ws <- makeX(X_scaled_ws, na.impute =TRUE)

cv_fit_ws <- cv.glmnet(X_scaled_ws, Response, alpha = 1)
lambda_min <- cv_fit_ws$lambda.min

lasso_model_scaled_ws <- glmnet(X_scaled_ws, Response, family="binomial", alpha=1, lambda=lambda_min)
coefficients_scaled_ws <- coef(lasso_model_scaled_ws)
print(coefficients_scaled_ws)

# checking influence on response duration with or without secondary endpoints
X2_scaled <- X2 %>%
  mutate(across(where(is.numeric), scale))

X2_scaled <- makeX(X2_scaled, na.impute =TRUE)

cv2_fit <- cv.glmnet(X2_scaled, Response, alpha = 1)
lambda_min <- cv2_fit$lambda.min

lasso_model2_scaled <- glmnet(X2_scaled, Response_dur, alpha=1, lambda=lambda_min)
coefficients2_scaled <- coef(lasso_model2_scaled)
print(coefficients2_scaled)
lasso_model2_scaled
# without secondary endpoints

X2_without_sec <- X2[,1:8]

X2_scaled_ws <- X2_without_sec %>%
  mutate(across(where(is.numeric), scale))

X2_scaled_ws <- makeX(X2_scaled_ws, na.impute =TRUE)

cv2_fit_ws <- cv.glmnet(X2_scaled_ws, Response_dur, alpha = 1)
lambda_min <- cv2_fit_ws$lambda.min

lasso_model2_scaled_ws <- glmnet(X2_scaled_ws, Response_dur, alpha=1, lambda=lambda_min)
coefficients2_scaled_ws <- coef(lasso_model2_scaled_ws)
print(coefficients2_scaled_ws)




#----------------------------------------------------------------------------------
#Confidence intervals + SE
se <- function (p,n) {
  se_p <- sqrt(p*(1-p)/n)
  return(se_p)
}
#Followed protocol
n <- 29
conf_level <- 0.95

prop_result <- prop.test(19, n, conf.level = conf_level)
print(prop_result)
prop_ci <- prop_result$conf.int
print(paste("Confidence interval using prop.test:", prop_ci[1], "-", prop_ci[2]))
p <- prop_result$estimate
print(paste("SE:", se(p,n)))
#Adverse events

prop_result <- prop.test(8, n, conf.level = conf_level)
print(prop_result)
prop_ci <- prop_result$conf.int
print(paste("Confidence interval using prop.test:", prop_ci[1], "-", prop_ci[2]))
p <- prop_result$estimate
print(paste("SE:", se(p,n)))
