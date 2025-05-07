install.packages("mice")
install.packages("VIM")
install.packages("naniar")
install.packages("Hmisc")
install.packages("Metrics")

library(mice)
library(VIM)
library(stats)
library(naniar)
library(zoo)
library(Hmisc)
library(Metrics)

data <- read.csv("G:/项目/All data Cyert for paper.csv")
# Determine whether each value is a missing value
is.na(data)
sum(is.na(data))

# Calculate the missing percentage of each column
Rate <- function(x){sum(is.na(x))/length(x)*100}
apply(data, 2, Rate)

# Calculate whether the sample is complete
compelete.cases(data)
sum(!complete.cases(data))

# Count the missing of each line
apply(data, 1, Rate)
# Look for features with missing values
md.pattern(data)
# Visualize the above results
aggr_plot <- aggr(data, col=c('navyblue', 'red'),
                  numbers = TRUE,
                  SORTVARS = TRUE,
                  labels = names(data),
                  cex.axis = 0.7,
                  gap = 3,
                  ylab = c("Missing data count", "Histogram of missing data"))


# Explore the relationship between missing variables
x<-as.data.frame(abs(is.na(data)))
head(data, n=5)

head(x, n=5)
# Calculate the relevant coefficient matrix
y <- x[, apply(x, MARGIN=2, sd) > 0]
# Calculate the relationship between missing variables
cor(y)

numeric_data <- data[, sapply(data, is.numeric)]
cor(numeric_data, y, use = "pairwise.complete.obs")

#Perform missing data MCAR test
mcar_test(data)

#Based on a p-value of less than 0.05, the missing value is not MCAR
#Visualize the distribution of missing data
vis_miss(data)
gg_miss_upset(data)