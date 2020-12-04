# Teodora Tockovska
# Dec 1, 2020

# This script contains the code to calculate the Fisher exact test and logistic regression models. Paneled figures are also created in this script.

# Loading my script which contains functions and libraries
source("Functions.R")

# --- Reading Data and Performing Statistical Models ----

# --- Read the host-associated data files and save the information into variables ---
file_contents_hosts <- read_files("combined_host_names_23Nov2020.txt", "host_significance_dataframe_23Nov2020.csv")
host_names <- file_contents_hosts[[1]]
significance_dataframe_hosts <- file_contents_hosts[[2]]

# Calculate the logistic regression model by regressing Significance on the Maximum Distance (km) and passing "binomial" in the family parameter. Then, view the summary of the model. The maximum distance predictor is insignficant.
model_host <- glm(Significance ~ Maximum.Distance..Km., data = significance_dataframe_hosts, family="binomial")
summary(model_host)

# --- Read the environmental data files ---
env_file_contents <- read_files("combined_study_ids_ENV_23Nov2020.txt", "ENV_significance_dataframe_23Nov2020.csv")
study_ids <- env_file_contents[[1]]
significance_dataframe_ENV <- env_file_contents[[2]]

# Different types of environmental materials
table(gsub(pattern = "StudyID_[0-9.]+_", x = rownames(significance_dataframe_ENV), replacement = ""))

# Calculate the logistic regression model by regressing Significance on the Maximum Distance (km) and passing "binomial" in the family parameter. Then, view the summary of the model. The maximum distance predictor is insignficant.
model_env <- glm(Significance ~ Maximum.Distance..Km., data = significance_dataframe_ENV, family="binomial")
summary(model_env)

# Remove unnecessary variables.
rm(file_contents_hosts, env_file_contents, model_host, model_env)

# Creating the confusion matrix using the funciton create_confusion_matrix().
significance_matrix <- create_confusion_matrix(significance_dataframe_hosts, significance_dataframe_ENV)

# computing the chi-square test: I was using the function below to calculate the p-values but it was giving an error message "Chi-squared approximation may be incorrect". I found this thread that said for low counts, I should calculate the Fisher's Exact test: https://stats.stackexchange.com/questions/81483/warning-in-r-chi-squared-approximation-may-be-incorrect 
chisq.test(significance_matrix)
fisher.test(significance_matrix)

# Anyway, the p-value is insignificant so that represents there is no association between the maximum distance (km) of sample extraction and microbial diversities within samples.

# Next, I created a dataframe that contains the confusion matrix, and a new column containing the  total amount of samples. This dataframe will be used for plotting in the next section.
df <- data.frame(significance_matrix)
df$Total <- 0
df$Total[1] <- sum(df$Significant[1], df$InSignificant[1])
df$Total[2] <- sum(df$Significant[2], df$InSignificant[2])

# --- Creating Figures ----

# --- Building the paneled world maps for env and host data, respectively ---

coordinate_env_joined <- read_csv("dataframe_Environmental.csv") %>% data.frame()
coordinate_host_joined <- read_csv("dataframe_hostAssociated.csv") %>% data.frame()

worldmap_env <- plot_world_map(coordinate_env_joined, 
                          "#882255", 
                          Longitude, 
                          Latitude, 
                          "Environmental",
                          shape = 21,
                          tiff = F)

worldmap_host <- plot_world_map(coordinate_host_joined, 
                          "#D55E00", 
                          Longitude, 
                          Latitude, 
                          "Host-Associated",
                          shape = 24,
                          tiff = F)

tiff("panelled_world_map.tiff", units="in", width = 8, height = 8, res = 300)
plot_grid(worldmap_env, worldmap_host, nrow = 2, labels = c("A", "B"))
dev.off()

# --- Building the summary plot of env vs host significant samples ---

# tiff(filename = "Significance_dotplot.tiff", units="in", width = 5, height = 4, res = 300)
shapes_plot <- ggplot(data=df, aes(x=rownames(df), y=Significant, group=1, shape = rownames(df), colour = rownames(df))) +
  scale_colour_manual(values=c("Env" = "#882255", "Host" = "#D55E00")) +
  geom_line(size = 0.75, colour = "black") +
  geom_point(size = 4) +
  expand_limits(y = 8) + 
  expand_limits(y = 0) +
  geom_text(aes(label=paste(Significant, "/", Total, sep = "")), hjust=1.6, colour="black", size=3.5)+
  theme_minimal() +
  xlab("Microbial Data Type") +
  ylab("Number of Significant Datasets") +
  ggtitle("Significant Samples from\nMicrobial Data") +
  theme(legend.position = "none") + 
  theme(text = element_text(size = 11),
        plot.title = element_text(hjust = 0.5))
print(shapes_plot)
# dev.off()

# --- Building the paneled figures for the env and host data ---

# Contains the plots of significance regressed on maximum distances (km) for env and host data, respectively.
plot1_logreg <- plot_significance(significance_dataframe_ENV, "Environmental", 21, "#882255")
plot2_logreg <- plot_significance(significance_dataframe_hosts, "Host", 24, "#D55E00")

# putting together in one panel
panel1_plots <- plot_grid(plot1_logreg, plot2_logreg, nrow = 2, labels = c("B", "C"))

# adding the summary figure and labelling them.
tiff("LogReg_Plots.tiff", units="in", width = 8, height = 4.5, res = 300)
plot_grid(shapes_plot, panel1_plots, ncol = 2, labels = c("A", "", ""), align = "v")
dev.off()
