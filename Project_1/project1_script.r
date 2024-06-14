library(table1)
library(pwr)
library(dplyr)
library(tidyr)
library(Biobase)
library(GEOquery)


### calculation of the methylation data

# data_txt <- getGEO(filename='GSE149747_series_matrix.txt.gz', getGPL=FALSE)
# pheno_txt_table <- as.data.frame(data_txt@phenoData@data)

# sample annotation file
# sample_annotation <- pheno_txt_table[, c('title', 'geo_accession')]
# sample_annotation$Tissue = 'Saliva'
# sample_annotation$Age <- substr(pheno_txt_table$`age:ch1` , 1, 2)
# sample_annotation$Sex <- ifelse(pheno_txt_table$`gender:ch1` == "male", 0, 1)
# sample_annotation$Group <- pheno_txt_table$`individual:ch1`
# sample_annotation$Timepoint <- pheno_txt_table$`timepoint (intervention):ch1`
# names(sample_annotation)[names(sample_annotation) == 'title'] <- 'Sample_Name'

# methylation data
#signals <- read.csv('GSE149747_MDL_Matrix_Signal_intensities.csv.gz')
#betas <- data.frame("ID_REF" = signals$ID_REF)

# sample names
#list_of_names <- sample_annotation$Sample_Name
#list_of_names2 <- gsub("-", ".", list_of_names)

# beta values
# for (i in list_of_names2) {
#   
#   unmethylated_name <- paste(i, ".Unmethylated.signal", sep="")
#   methylated_name <- paste(i, ".Methylated.signal", sep="")
#   betas[, c(i)] <- signals[, c(methylated_name)]/(signals[, c(unmethylated_name)] + signals[, c(methylated_name)])
#   
# }

#DNAmAge calculation and saving the results to csv files

#install.packages("devtools")
#devtools::install_github("yiluyucheng/dnaMethyAge")
#library('dnaMethyAge')
#clock_name <- "HorvathS2013"
#clock_name <- "HannumG2013"
#horvath_age <- methyAge(betas, clock=clock_name)
#data <- read.csv("methylation_data2.csv")
#
# data <- data[!duplicated(data$ID_REF), ]
# 
# row_names <- data[,1]
# data <- data[,-1]
# rownames(data) <- row_names
# 
# horvath_age <- methyAge(data, clock=clock_name)
# write.csv(horvath_age, "hannum_age.csv", row.names=TRUE)




# data preparation

horvath <- read.csv("horvath_age.csv")
hannum <- read.csv("hannum_age.csv")
sample_data <- read.csv("sample_annotation.csv")

sample_names <- sample_data$Sample_Name

sample_data <- sample_data[, c("Sample_Name", "Age", "Group", "Timepoint")]
sample_data$patient_id <- substr(sample_data$Sample_Name, 1, 8)

hannum <- hannum[, c("Sample", "mAge")]
colnames(hannum) <- c("Sample_Name", "mAgeHannum")
hannum$Sample_Name <- sample_names

horvath <- horvath[, c("Sample", "mAge")]
colnames(horvath) <- c("Sample_Name", "mAgeHorvath")
horvath$Sample_Name <- sample_names

data <- merge(x=sample_data, y=horvath, by="Sample_Name")
data <- merge(data, hannum, by="Sample_Name")
data$Timepoint <-  factor(data$Timepoint, levels = c("Before", "Midpoint", "After"))


# number of observations by group 

table(data$Group, data$Timepoint)


# define age differences
data$diff_hannum <- data$Age - data$mAgeHannum
data$diff_horvath <- data$Age - data$mAgeHorvath



#hist(data[data$Group == "Control", "diff_hannum"])
#hist(data[data$Group == "Control", "diff_horvath"])

#hist(data[data$Group == "Subject", "diff_hannum"])
#hist(data[data$Group == "Subject", "diff_horvath"])


# compare the mean and median ages 

table1(~ Age 
       + mAgeHorvath
       + diff_horvath
       + mAgeHannum
       + diff_hannum| Timepoint, data=data[data$Group == "Control",])
table1(~ Age
       + mAgeHorvath
       + diff_horvath
       + mAgeHannum
       + diff_hannum| Timepoint, data=data[data$Group == "Subject",])


boxplot(data$Age ~ data$Group)
boxplot(data$mAgeHorvath ~ data$Group)
boxplot(data$mAgeHannum ~ data$Group)

boxplot(data$diff_horvath ~ data$Group)

boxplot(data$Age ~ data$Group)


## testing the differences: Horvath

data_pivot <- data[, c("patient_id", "Group", "Timepoint", "mAgeHorvath")] %>% 
  pivot_wider(names_from = Timepoint, values_from = mAgeHorvath)

data_pivot$change <- data_pivot$After - data_pivot$Before


change_control <- unlist(as.vector(data_pivot %>% filter(Group=="Control") %>% select(change)))
change_treatment <- unlist(as.vector(data_pivot %>% filter(Group=="Subject") %>% select(change)))
change_control <- change_control[!is.na(change_control)]
change_treatment <- change_treatment[!is.na(change_treatment)]


#authors are using 2-tailed t-test
t.test(change_treatment, change_control, alternative="two.sided")
#p-value = 0.0834 => difference is not significant

#one-tailed t-test seems more reasonable as we want to check whether treatment group has better results (i.e. lower change)
t.test(change_treatment, change_control, alternative="less")
#p-value = 0.0417 => treatment group has significantly lower change of age

mean(change_treatment) # -1.41
mean(change_control) # 0.46



## testing the differences: Hannum

data_pivot <- data[, c("patient_id", "Group", "Timepoint", "mAgeHannum")] %>% 
  pivot_wider(names_from = Timepoint, values_from = mAgeHannum)

data_pivot$change <- data_pivot$After - data_pivot$Before

change_control <- unlist(as.vector(data_pivot %>% filter(Group=="Control") %>% select(change)))
change_treatment <- unlist(as.vector(data_pivot %>% filter(Group=="Subject") %>% select(change)))
change_control <- change_control[!is.na(change_control)]
change_treatment <- change_treatment[!is.na(change_treatment)]


#authors are using 2-tailed t-test
t.test(change_treatment, change_control, alternative="two.sided")
#p-value = 0.1686 => difference is not significant

#one-tailed t-test seems more reasonable as we want to check whether treatment group has better results (i.e. lower change)
t.test(change_treatment, change_control, alternative="less")
#p-value = 0.9157 => difference is not significant

mean(change_treatment) # 0.2790348
mean(change_control) # -1.072345


### power calculations

# first we check the power for different effect sizes

effect_sizes = seq(0, 2, 0.05)

power_two_sided <- pwr.2p2n.test(h = effect_sizes, 
                                 n1=length(change_control), 
                                 n2=length(change_treatment),
                                 alternative = 'two.sided')$power

power_one_sided <- pwr.2p2n.test(h = effect_sizes, 
                                 n1=length(change_control), 
                                 n2=length(change_treatment),
                                 alternative = 'greater')$power

plot(effect_sizes, power_two_sided,
     ylab='power',
     xlab='effect size', 
     col='blue')
points(effect_sizes, power_one_sided,
       col='green')
legend(x = 'topright',
       legend=c('two sided', 'one sided'),
       fill=c('blue', 'green'))


# we calculate the effect size for our means by Cohens method 

effect_size = cohensD(change_treatment, change_control)

# for the two-sided test we get power = 0.425
pwr.2p2n.test(h = effect_size, 
              n1=length(change_control), 
              n2=length(change_treatment),
              alternative = 'two.sided')$power


# for the one-sided test power = 0.55
pwr.2p2n.test(h = effect_size, 
              n1=length(change_control), 
              n2=length(change_treatment),
              alternative = 'greater')$power
