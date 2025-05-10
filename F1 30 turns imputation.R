install.packages("Metrics")  # 仅需安装一次
library(Metrics)
library(ggplot2)
library(tibble) 
library(mice)
library(dplyr)
library(stats)    

# read data
data <- read.csv("G:/项目/All data Cyert for paper.csv")

# 30turns round imputation
#Data containing missing values were multiply interpolated 30 times with up to 10 rounds of iterations using predictive mean matching, and random seeds were set to ensure reproducible results
imputed_data <- mice(data, method = "pmm", m = 30, maxit = 10, seed = 123)

# Inspect the imputation process
summary(imputed_data)


# Converts the mice() interpolation result object, imputed_data, into a long form dataframe (containing all interpolated versions) into completed_data.
completed_data <- complete(imputed_data, action = "long")
#missing check
sum(is.na(completed_data))  


# View the imputed dataset
head(completed_data)
view(completed_data)



# Calculate the mean value of each round of interpolated data according to .imp
mean_by_round <- aggregate(. ~ .imp, 
                           data = completed_data[sapply(completed_data, is.numeric)], 
                           FUN = function(x) mean(x, na.rm = TRUE))
view(mean_by_round)


# Rank each column (in descending order)
ranked_data <- apply(mean_by_round[, -1], 2, function(x) rank(-x, ties.method = "first"))
view(ranked_data)

ranked_data <- as.data.frame(ranked_data)
# Calculate the sum of columns 1-24 in each row
rank_sums <- rowSums(ranked_data[, 1:24])

# Add the calculated rank_sums to ranked_data as a new column
ranked_data <- cbind(ranked_data, rank_sum = rank_sums)

# Sort by rank_sum in descending order
ranked_data <- ranked_data[order(-ranked_data$rank_sum), ]

# View Results
View(ranked_data)

#Delete the rounds ranked in the top 10% and the bottom 10%
data_80 <- subset(completed_data, !(.imp %in% c(2, 26, 28, 10, 19, 23)))
view(data_80)

#Group by round and take the average value of the numerical variables in each group
imputed_means <- data_80 %>%
  group_by(.id) %>% 
  summarise(
    YORF = first(YORF),   
    NAME = first(NAME),   
    across(where(is.numeric), mean) 
  )
view(imputed_means)

#Screen the sample column
data_80 <- imputed_means[, 2:27]
view(data_80)

#save the file
write.csv(data_80, "G:/项目/data_80.csv", row.names = TRUE)

#calculate the mean value
mean_values4 <- apply(data_80[, 3:26], 2, mean, na.rm = TRUE)
print(mean_values4)
#calculate the standard deviation
imputed_sd3 <- apply(data_80[, 3:26], 2, sd, na.rm = TRUE)
print(imputed_sd3)
#summary
sum_mean3 <- apply(data_80[, 3:26], 2, summary, na.rm = TRUE)
print(sum_mean3)

#calculate the original data mean value
mean_values <- apply(data[, 3:26], 2, mean, na.rm = TRUE)
print(mean_values)

#Convert the matrix or the data frame with row names to the standard data frame and extract the row names as the "variable" column.
mean_values <- as.data.frame(mean_values) %>% rownames_to_column(var = "variable") 
mean_values4 <- as.data.frame(mean_values4) %>% rownames_to_column(var = "variable") 

#Converting variable columns to characters
mean_values4$variable <- as.character(mean_values4$variable)
mean_values$variable <- as.character(mean_values$variable)

#Merge the original and interpolated mean tables, aligned by variable name (“variable”).
m_data1 <- merge(mean_values, mean_values4, by = "variable", all = TRUE)
str(m_data1)

#Plot a side-by-side bar graph of the original and interpolated means.
ggplot(m_data1, aes(x = variable)) +
  geom_bar(aes(y = mean_values, fill = "Original"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = mean_values4, fill = "Impute"), stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Comparison of mean values of Original and Imputed Values",
       x = "Variable",
       y = "Mean values") +
  scale_fill_manual(name = "Data Type", values = c("Original" = "blue", "Impute" = "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  scale_y_continuous(expand = expansion(mult = c(0, 0.001)))  

#Calculate the standard deviation for each column
sd_values <- apply(data[, 3:26], 2, sd, na.rm = TRUE)
print(sd_values)

#Convert the matrix or the data frame with row names to the standard data frame and extract the row names as the "variable" column.
sd_values <- as.data.frame(sd_values) %>% rownames_to_column(var = "variable") 
imputed_sd3  <- as.data.frame(imputed_sd3 ) %>% rownames_to_column(var = "variable") 

#Converting variable columns to characters
sd_values$variable <- as.character(sd_values$variable)
imputed_sd3$variable <- as.character(imputed_sd3$variable)

#Merge the original and interpolated standard deviation tables, aligned by variable name (“variable”).
s_data1 <- merge(sd_values, imputed_sd3, by = "variable", all = TRUE)
str(s_data1)

#Plot a side-by-side bar graph of the original and interpolated standard deviation.
ggplot(s_data1, aes(x = variable)) +
  geom_bar(aes(y = sd_values, fill = "Original"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = imputed_sd3, fill = "Impute"), stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Comparison of standard deviation of Original and Imputed Values",
       x = "Variable",
       y = "SD values") +
  scale_fill_manual(name = "Data Type", values = c("Original" = "blue", "Impute" = "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # 旋转x轴标签
  scale_y_continuous(expand = expansion(mult = c(0, 0.001)))  # 调整Y轴对比度

#Calculate the MAD for each column
mad_values <- apply(data[, 3:26], 2, function(x) mad(x, center = median(x, na.rm = TRUE), constant = 1.4826, na.rm = TRUE))
print(mad_values)
mad_values1 <- apply(imputed_data1[, 3:26], 2, function(x) mad(x, center = median(x, na.rm = TRUE), constant = 1.4826, na.rm = TRUE))
print(mad_values1)
mad_values2 <- apply(data_80[, 3:26], 2, function(x) mad(x, center = median(x, na.rm = TRUE), constant = 1.4826, na.rm = TRUE))
print(mad_values2)

#Convert the matrix or the data frame with row names to the MAD frame and extract the row names as the "variable" column.
mad_values <- as.data.frame(mad_values) %>% rownames_to_column(var = "variable") 
print(mad_values)
mad_values2 <- as.data.frame(mad_values2) %>% rownames_to_column(var = "variable") 
print(mad_values2)

#Side-by-side MAD histogram of raw and interpolated data
par(mfrow=c(1,2))
barplot(mad_values$mad_values, 
        names.arg = mad_values$variable, 
        col = "skyblue", 
        las = 2, 
        main = "Original MAD", 
        xlab = "Variable", 
        ylab = "MAD Value")
barplot(mad_values2$mad_values2, 
        names.arg = mad_values2$variable, 
        col = "red", 
        las = 2, 
        main = "30 turn imputed MAD", 
        xlab = "Variable", 
        ylab = "MAD Value")




