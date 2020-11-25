# Teodora Tockovska
# Nov 24, 2020

# I will use this script to perform the chi-squared tests and the logistic regression

# --- Loading functions ----

library("tidyverse")

# function to read the files:
read_files <- function(file1, file2){
  return_lst <- list()
  sample_names <- read.table(file1) %>% unlist() %>% as.character()
  significance_dataframe <- read_csv(file2) %>% data.frame() %>% `rownames<-`(sample_names)
  return_lst[[1]] <- sample_names
  return_lst[[2]] <- significance_dataframe
  return(return_lst)
}

# For the chi-squared tests, I will need to create a table. I created 2 helper functions and 1 main function to build the confusion matrix. 
make_table <- function(){
  conf_matrix <- matrix(0, nrow = 2, ncol = 2)
  rownames(conf_matrix) <- c("Env", "Host")
  colnames(conf_matrix) <- c("Significant", "InSignificant")
  return(conf_matrix)
}

fill_table <- function(sample_df){
  lst_info <- list()
  lst_info[[1]] <- sum(sample_df$Significance)
  lst_info[[2]] <- nrow(sample_df) - sum(sample_df$Significance)
  return(lst_info)
}

create_confusion_matrix <- function(data1, data2){
  confusion_matrix <- make_table()
  
  # Filling in the environmental information:
  env_info <- fill_table(data2)
  confusion_matrix[1,1] <- env_info[[1]]
  confusion_matrix[1,2] <- env_info[[2]]
  
  # Filling in the host information:
  host_info <- fill_table(data1)
  confusion_matrix[2,1] <- host_info[[1]]
  confusion_matrix[2,2] <- host_info[[2]]
  
  return(confusion_matrix)
}

plot_significance <- function(sample_df){
  plot(sample_df$Maximum.Distance..Km., 
       sample_df$Significance, 
       pch = 16, 
       xlab = "Maximum distance (Km)", 
       ylab = "Significance",
       main = "Plotting Significance Against Maximum\nDistance (Km) between Samples")
}

# --- Reading Data and Performing Statistical Models ----

# Read the host-associated data files:
file_contents_hosts <- read_files("combined_host_names_23Nov2020.txt", "host_significance_dataframe_23Nov2020.csv")
host_names <- file_contents_hosts[[1]]
significance_dataframe_hosts <- file_contents_hosts[[2]]

model_host <- glm(Significance ~ Maximum.Distance..Km., data = significance_dataframe_hosts, family="binomial")
summary(model_host)

plot_significance(significance_dataframe_hosts)

# Read the environmental data files:
env_file_contents <- read_files("combined_study_ids_ENV_23Nov2020.txt", "ENV_significance_dataframe_23Nov2020.csv")
study_ids <- env_file_contents[[1]]
significance_dataframe_ENV <- env_file_contents[[2]]

table(gsub(pattern = "StudyID_[0-9.]+_", x = rownames(significance_dataframe_ENV), replacement = ""))

model_env <- glm(Significance ~ Maximum.Distance..Km., data = significance_dataframe_ENV, family="binomial")
summary(model_env)

plot_significance(significance_dataframe_ENV)

# Creating the confusion matrix:
significance_matrix <- create_confusion_matrix(significance_dataframe_hosts, significance_dataframe_ENV)

# computing the chi-square test: I was using the function below to calculate the p-values but it was giving an error message "Chi-squared approximation may be incorrect". I found this thread that said for low counts, I should calculate the Fisher's Exact test: https://stats.stackexchange.com/questions/81483/warning-in-r-chi-squared-approximation-may-be-incorrect 
chisq.test(significance_matrix)
fisher.test(significance_matrix)

# Anyway, the p-value is insignificant.

# Creating the Summary
df <- data.frame(significance_matrix)
df$Type <- rownames(df)
df$Total <- 0
df$Total[1] <- sum(df$Significant[1], df$InSignificant[1])
df$Total[2] <- sum(df$Significant[2], df$InSignificant[2])

tiff(filename = "Significance_dotplot.tiff", units="in", width = 5, height = 4, res = 300)
ggplot(data=df, aes(x=Type, y=Significant, group=1, colour = Type, shape = Type)) +
  geom_line(size = 0.75, colour = "black") +
  geom_point(size = 4) +
  expand_limits(y = 8) + 
  expand_limits(y = 0) +
  geom_text(aes(label=paste(Significant, "/", Total, sep = "")), hjust=1.6, colour="black", size=3.5)+
  theme_minimal() +
  xlab("Microbial Data Type") +
  ylab("Number of Significant Samples") +
  ggtitle("Significant Samples from Environmental and Host-associated\nMicrobial Data") +
  theme(legend.position = "none") + 
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
dev.off()
