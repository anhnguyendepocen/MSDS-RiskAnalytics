library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(GGally)
library(ggthemes)
library(scales)
library(reshape2)
library(skimr)
library(gridExtra)
library(Ecdat)
library(faraway)
library(fGarch)

#####################################################################
######################### Ford Case Study ###########################
#####################################################################

theme_set(theme_sjplot())

# Theme Overrides
theme_update(plot.title = element_text(hjust = 0.5),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 10),
             axis.title = element_text(face = "bold", size = 12, colour = "steelblue4"),
             legend.position = "top", legend.title = element_blank())

path.data <- "D:/Projects/MSDS-RiskAnalytics/datasets"
setwd(path.data)

# 4.11
# Ford Case Study

data.ford <- as.data.table(read.csv("ford.csv", header = T))

colnames(data.ford) <- c("n", "date", "return")

head(data.ford)

mean(data.ford$return)
sd(data.ford$return)
median(data.ford$return)

ggplot(data.ford, aes(sample = return)) +
  geom_qq() +
  geom_qq_line()

shapiro.test(data.ford$return)

n <- nrow(data.ford)
df_canidates <- c(1, 2, 4, 6, 10, 20)
q_range <- (1:n) / (n+1)

plots <- rep(0, 6)

i <- 1

for(i in range(1, 6)) {
  df <- df_canidates[i]
  
  data <- data.table(ret = data.ford$return, theoretical = qt(q_range, df))
  
  print(head(data))
  
  ggplot(data, aes(sample = ret)) +
    geom_qq(line.p = data$theoretical) +
    geom_qq_line()
}

grid.arrange(plots, nrow = 3)



df <- df_canidates[i]

data <- data.table(ret = data.ford$return, theoretical = qt(q_range, df))
data$theoretical <- sort(data$theoretical)

model <- lm(qt(c(0.25,0.75), df = df) ~ quantile(data$ret,c(0.25,0.75)))

ggplot(data, aes(x = sort(ret), y = theoretical)) +
  geom_point() +
  geom_abline(slope = model$coefficients[2], intercept = model$coefficients[1])
  labs(title = paste("QQ-Plot: Returns vs t-distribution, df = ", df), 
       xlab = "Theoretical", 
       ylab = "Returns")
  
qqplot(data$ret, data$theoretical)
