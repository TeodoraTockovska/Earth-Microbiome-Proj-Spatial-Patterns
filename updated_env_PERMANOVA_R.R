# Teodora Tockovska
# Nov 19. 2020

# This script contains the code to perform spatial pattern analyses on free-living microbiomes from the Earth Microbiome Project.

# Loading my script which contains functions 
source("Functions.R")

# Read the host associated metadata into R
env_microbiome_dataset <- read.csv(file = "free_living_EMP.csv", sep = "\t")
nrow(env_microbiome_dataset)

# Checking and Removing any NA values:
env_microbiome_dataset <- check_and_remove_NAs(env_microbiome_dataset)
nrow(env_microbiome_dataset)

# Filter the dataset to extract the environmental material samples with 30 or more records. I chose the environmental materials because I am interested to see what the microbiome diversity is within each sample type. Bacteria are specified to their environments so the materials felt like the most logical thing to use to select data.
env_with_30_records <- as.data.frame(table(env_microbiome_dataset$study_id)) %>% filter(Freq>=30)
env_with_30_records[order(env_with_30_records$Freq, decreasing = T),]

# I created a list of dataframes called "lst_df_env_with_30_records". The dataframes represent the hosts' information that was indexed from the original EMP dataset. The for loop is used to loop through the environmental materials to extract the sample types. I used grep() to match the strings so that I can extract the datasets. The dataframe will be used to extract the coordinates and sample IDs for the hosts.
lst_df_env_with_30_records <- list()
for (i in 1:length(env_with_30_records$Var1)){
  lst_df_env_with_30_records[[i]] <- env_microbiome_dataset[grep(env_with_30_records$Var1[[i]],
                                                                 env_microbiome_dataset$study_id),]
}

# The next goal is to plot the coordinate of each host to visualize their samples and determine outliers. However, there are hosts that were only sampled from 1 unique coordinate (Lat, Long). I call the function chosen_samples() that will check the unique coordinates per host and exclude hosts that were sampled from only 1 region (unique coordinate). The function will also plot the coordinates of the chosen hosts. The function also returns a list of the selected hosts. The function also prints out the number of selected hosts. There are 17 environmental samples that were selected based on the criteria.
greater_than_1_region <- chosen_samples(lst_df_env_with_30_records, env_with_30_records, type = "Study ID")

# Remove variables that are no longer needed:
rm(lst_df_env_with_30_records, env_with_30_records)

# Next, extract the coordinate dataframes per host into a list by using the function get_coordinates(). Using the list of coordinate dataframes, I check which UTM zones the samples were taken from. This is important because samples can be taken from all over the world and I would like to ensure that the samples that I am working with are not too far away from each other. So, I check the zones per coordinate dataframe and I create a dictionary that stores the UTM zone and the associated longitudes for each host. The dictionary is stored into a variable called "dict_UTM_zones".
coordinate_list <- get_coordinates(data = greater_than_1_region, lat = latitude_deg, long = longitude_deg, type = "study_id")

# Checking whether the coordinates fall in the correct hemispheres:
coordinate_list %>% map(~.x %>% check_hemisphere(6))

# There were two samples from Canada (coordinate_list[[21]]) that showed were from NE... that makes no sense. I found that IDs 415 and 416 had strange longitudes, so I will remove these. Additionally, why are there samples for the North Atlantic Ocean with both NE and NW coordinates (coordinate_list[[31]])? Looks like there was a positive longitude. Note: Kenya's coordinates (coordinate_list[[28]]) are correct because the Equator passes through the country so some of the coordinates are in the Northern hemisphere vs some are in the Southern hemisphere.

# First, dealing with the Canadian samples. I will remove the 2 samples that had non-sensical longitudes:
canadian_ids <- which(coordinate_list[[21]]$country == "Canada")
wrong_ids <- which(coordinate_list[[21]][canadian_ids,1] > 0)
coordinate_list[[21]] <- coordinate_list[[21]][-canadian_ids[wrong_ids],]

# Second, dealing with the North Atlantic Ocean's mislabeled coordinates:
mislabeled_NAO <- which(coordinate_list[[31]][,1] > 0)
coordinate_list[[31]][mislabeled_NAO,1] <- -coordinate_list[[31]][mislabeled_NAO,1]

# Validating the changes to the coordinates. Everything is looking good.
coordinate_list %>% map(~.x %>% check_hemisphere(6))

# what env materials are part of the studyids?
info_env_cohorts <- coordinate_list %>% map(~.x %>% select(study_id, country, env_material) %>% table())
# capture.output(info_env_33cohorts, file = "info_env_33cohorts_28Oct2020.txt")

# coordinate_list <- coordinate_list %>% map(~.x %>% unite("Study_Material", study_id:env_material, sep = "-", remove = T))
dict_UTM_zones <- coordinate_list %>% map(~.x %>% group_by(Longitude) %>% check_zones())

# conducting the spatial pattern analysis because I cannot compare microbiomes of the hosts from such large distances. 
lst_uniq_UTM_zones <- return_unique_UTM_zones(dict_UTM_zones)

# How many hosts were sampled within 1 UTM zone? 18 environmental variables, whereas the rest (15) were sampled from multiple UTM zones. Next, I extract the env samples that were sampled from 1 UTM zone and put them into a list by getting their ids and indexing from "coordinate_list".
sum(lst_uniq_UTM_zones %>% map(~.x %>% length(.))==1)
single_zone_ids <- which(lst_uniq_UTM_zones %>% map(~.x %>% length(.))==1)
single_zone_coordinate_datasets <- coordinate_list[single_zone_ids]
lst_uniq_UTM_zones[single_zone_ids]

# what env materials are part of the studyids? I had to run this after 
info_env_17singleUTMcohorts <- single_zone_coordinate_datasets %>% map(~.x %>% select(study_id, country, env_material) %>% table())
# capture.output(info_env_17singleUTMcohorts, file = "info_env_singleUTMcohorts_28Oct2020.txt")

# I extracted the host species' dataframes of coordinates that were sampled from more than 1 UTM zone. 
multiple_UTM_zones <- coordinate_list[-single_zone_ids]
lst_uniq_UTM_zones[-single_zone_ids]

# Checking the names, which generates a tibble
study_id_names <- greater_than_1_region %>% map(~.x %>% group_by(study_id) %>% summarise(n=n()) %>% select(study_id))
# greater_than_1_region %>% map(~.x %>% group_by(host_common_name_provided) %>% summarise(n=n()) )

# Hosts that were sampled from 1 UTM zone are the following: (Numbers in parentheses represent the total names recorded per host. The updated names per host are shown to the left of the parenthesis.)
single_zone_study_id <- extract_names(study_id_names, single_zone_ids, "study_id")
single_zone_study_id <- gsub("^", "StudyID_", single_zone_study_id)
single_zone_study_id

# Adding names to the list:
names(single_zone_coordinate_datasets) <- single_zone_study_id

# Creating the list of dataframes with multiple UTM zone samples
multiple_zones_lst <- lst_uniq_UTM_zones[-single_zone_ids]

dict_multiple_UTM_zones <- dict_UTM_zones[-single_zone_ids]
study_id_longitudes_lst <- return_longitudes_per_UTM(multiple_zones_lst, dict_multiple_UTM_zones, "Multiple")

# Using the function extract_longitudes(), I extracted the hosts with the target longitudes. The function takes 2 things as parameters: the list of coordinate dataframes per host and the list of longitudes within a specific UTM zone per host. The longitudes from the second parameter will be used to select the longitudes from the coordinates dataframes for each host. Hence, I will extract the coordinates of interest for each host.
updated_coordinates <- extract_longitudes(multiple_UTM_zones, study_id_longitudes_lst)
updated_coordinates <- rename_samples_list(updated_coordinates, "StudyID")

output <- rename_inner_lists(updated_coordinates, "StudyID")
updated_coordinates <- output[[1]]
removed_hosts <- output[[2]]
rm(output)
updated_coordinates <- remove_empty(updated_coordinates)

multiple_zones_names <- extract_names(study_id_names, -single_zone_ids, "study_id", "list")
multiple_zones_names <- multiple_zones_names %>% map(~.x %>% gsub("^", "studyid ", .))
output <- return_ids(removed_hosts, updated_coordinates, "StudyID")

# I want to change the names so that they have "study_id_[ID]_num"
multiple_zones_names <- change_all_multiple_samples_names(output, multiple_zones_names)

# There are some samples that have very low number of records and I will need to filter samples<10 records from my dataset.
max_dist_table <- data.frame(max_dist(updated_coordinates, multiple_zones_names))
plot_type(updated_coordinates, multiple_zones_names, "StudyID")

# Split dataset by environmental materials ####

# I combined the lists so that I can do the PERMANOVA test
combined_list_study_ids <- c(single_zone_coordinate_datasets, updated_coordinates)
# combined_study_ids <- c(single_zone_study_id, multiple_zones_names)

# It is clear that 9 studies have multiple environmental materials. I need to split these up.
combined_list_study_ids %>% map(~.x %>% select(env_material) %>% unique() %>% nrow())
which(combined_list_study_ids %>% map(~.x %>% select(env_material) %>% unique() %>% nrow()) > 1)

combined_list_study_ids %>% map(~.x %>% select(env_material) %>% table())

#TODO: separate by ENV MATERIAL

# The idea behind these code chunks: I find the unique environment types and save the vector in a variable. Then, I find the indices for the one type. I will generate n variables with the indices (based on the number of env types found). Then, I generate n number of dataframes that contain the same information but separated by their environmental types. So, in this case, there were 2 different env types. I created 2 dataframes in the end. I also kept the samples with records >= 30

combined_list_study_ids <- subsetting_by_env_materials(combined_list_study_ids) %>% change_element_names()
combined_study_ids <- names(combined_list_study_ids)

max_dist_table <- max_dist(combined_list_study_ids, combined_study_ids) %>% as.data.frame()

# Filtering Outliers ####

# Now, I need to filter based on the max distances and the plots! First, remove the IDs that have records <=15. Then, remove the outliers.
remove_0dist <- which(max_dist_table$`Maximum Distance (Km)` > 0)
combined_list_study_ids <- combined_list_study_ids[remove_0dist]
combined_study_ids <- combined_study_ids[remove_0dist]
max_dist_table <- max_dist(combined_list_study_ids, combined_study_ids) %>% as.data.frame()

# remove samples with records < 15
remove_low_records <- which(max_dist_table$`Number of Records` >= 15)
combined_list_study_ids <- combined_list_study_ids[remove_low_records]
combined_study_ids <- combined_study_ids[remove_low_records]
max_dist_table <- max_dist(combined_list_study_ids, combined_study_ids) %>% as.data.frame()

# Remove the samples with very large distances (1000 (km)+)
remove_700dist <- which(max_dist_table$`Maximum Distance (Km)` < 700)
combined_list_study_ids <- combined_list_study_ids[remove_700dist]
combined_study_ids <- combined_study_ids[remove_700dist]

# checking if they are the same:
max_dist(combined_list_study_ids, combined_study_ids)
plot_type(combined_list_study_ids, combined_study_ids, "StudyID")

# removing outliers in the remaining samples:
combined_list_study_ids[[2]] <- combined_list_study_ids[[2]][which(combined_list_study_ids[[2]]$Longitude > -42.122),]
combined_list_study_ids[[3]] <- combined_list_study_ids[[3]][which(combined_list_study_ids[[3]]$Longitude > -42.122),]
combined_list_study_ids[[4]] <- combined_list_study_ids[[4]][which(combined_list_study_ids[[4]]$Longitude > -106.8),]
combined_list_study_ids[[7]] <- combined_list_study_ids[[7]][which(combined_list_study_ids[[7]]$Longitude < -110.868),]
combined_list_study_ids[[8]] <- combined_list_study_ids[[8]][which(combined_list_study_ids[[8]]$Longitude > 13.0),]
combined_list_study_ids[[13]] <- combined_list_study_ids[[13]][which(combined_list_study_ids[[13]]$Longitude > -110.8),]
combined_list_study_ids[[14]] <- combined_list_study_ids[[14]][which(combined_list_study_ids[[14]]$Latitude > 71.23),]
combined_list_study_ids[[15]] <- combined_list_study_ids[[15]][which(combined_list_study_ids[[15]]$Longitude < -155.74),]
combined_list_study_ids[[16]] <- combined_list_study_ids[[16]][which(combined_list_study_ids[[16]]$Longitude < -73.80),]
combined_list_study_ids[[17]] <- combined_list_study_ids[[17]][which(combined_list_study_ids[[17]]$Longitude < -79.5),]
combined_list_study_ids[[21]] <- combined_list_study_ids[[21]][which(combined_list_study_ids[[21]]$Longitude > -42.124),]
combined_list_study_ids[[22]] <- combined_list_study_ids[[22]][which(combined_list_study_ids[[22]]$Longitude > 8.76),]
combined_list_study_ids[[26]] <- combined_list_study_ids[[26]][which(combined_list_study_ids[[26]]$Longitude > 161.65),]
combined_list_study_ids[[28]] <- combined_list_study_ids[[28]][which(combined_list_study_ids[[28]]$Longitude < -86.0),]
combined_list_study_ids[[32]] <- combined_list_study_ids[[32]][which(combined_list_study_ids[[32]]$Longitude > -106.0),]
combined_list_study_ids[[38]] <- combined_list_study_ids[[38]][which(combined_list_study_ids[[38]]$Longitude > -104.80),]
combined_list_study_ids[[41]] <- combined_list_study_ids[[41]][which(combined_list_study_ids[[41]]$Longitude < 19.4),]
combined_list_study_ids[[43]] <- combined_list_study_ids[[43]][which(combined_list_study_ids[[43]]$Longitude < -147),]

combined_list_study_ids[[45]] <- combined_list_study_ids[[45]][which(combined_list_study_ids[[45]]$Latitude < 68.64),]
combined_list_study_ids[[45]] <- combined_list_study_ids[[45]][which(combined_list_study_ids[[45]]$Longitude < -149.56),]

combined_list_study_ids[[46]] <- combined_list_study_ids[[46]][which(combined_list_study_ids[[46]]$Longitude < -149.56),]
combined_list_study_ids[[47]] <- combined_list_study_ids[[47]][which(combined_list_study_ids[[47]]$Longitude < -149.56),]
combined_list_study_ids[[52]] <- combined_list_study_ids[[52]][which(combined_list_study_ids[[52]]$Longitude > -73.5),]
combined_list_study_ids[[53]] <- combined_list_study_ids[[53]][which(combined_list_study_ids[[53]]$Latitude < 41.04),]

# checking the plots and calculating the distances again:
updated_max_dist_table <- max_dist(combined_list_study_ids, combined_study_ids) %>% as.data.frame()

# remove the low records
updated_removelow_records <- which(updated_max_dist_table$`Number of Records` >= 15)
combined_list_study_ids <- combined_list_study_ids[updated_removelow_records]
combined_study_ids <- combined_study_ids[updated_removelow_records]

# Checking one last time:
max_dist(combined_list_study_ids, combined_study_ids)
plot_type(combined_list_study_ids, combined_study_ids, "StudyID (Filtered)")

# Lastly, making the groupings for the PERMANOVA tests:
combined_list_study_ids <- combined_list_study_ids %>% map(~.x %>% unite("groupings", Longitude:Latitude, sep = "-", remove = F))

# Plotting on world map
joined <- bind_rows(combined_list_study_ids)

plot_world_map(joined, 
               "#882255", 
               Longitude, 
               Latitude, 
               "World Map of the Filtered Environmental Microbiome Data", 
               tiff = T, 
               "Filtered_Env_map")

# Testing out ####

# Here I will plot the coordinates on the world map, but zoomed in to specific countries from which the samples were taken:

# library("rnaturalearth")
# library("rnaturalearthdata")
# 
# world <- ne_countries(scale = "medium", returnclass = "sf")
# class(world)
# 
# ggplot(data = world) +
#   geom_sf() +
#   coord_sf(xlim = c(-124.8, -124.4), ylim = c(48, 48.5), expand = FALSE) +
#   geom_point(data = combined_list_study_ids[[1]], aes(x = Longitude, y = Latitude), size = 5, shape = 21, fill = "blue")
# 
# # not meaningful
# ggplot(data = world) +
#   geom_sf() +
#   coord_sf(xlim = c(-110, -100), ylim = c(43, 55), expand = FALSE) + 
#   geom_point(data = combined_list_study_ids[[4]], aes(x = Longitude, y = Latitude), size = 5, shape = 21, fill = "blue")
# 
# # not meaningful
# ggplot(data = world) +
#   geom_sf() +
#   coord_sf(xlim = c(-120, -105), ylim = c(40, 50), expand = FALSE) + 
#   geom_point(data = combined_list_study_ids[[7]], aes(x = Longitude, y = Latitude), size = 5, shape = 21, fill = "blue")
# 
# ggplot(data = world) +
#   geom_sf() +
#   coord_sf(xlim = c(-90, -80), ylim = c(38, 45), expand = FALSE) + 
#   geom_point(data = combined_list_study_ids[[10]], aes(x = Longitude, y = Latitude), size = 5, shape = 21, fill = "blue")

# Conducting Spatial Pattern Analysis on 49 Env studies ####

# Before I do the PERMANOVA tests, I need to get the BIOM tables for the 11 host species. I will extract the sample IDs.So, I extract the sample IDs and saving them into files. I will transfer those files to Graham so that I can extract the BIOM tables.
for (i in 1:length(combined_list_study_ids)){
  file_name <- paste("studyid", as.character(i), "23Nov2020.tsv", sep = "_")
  sample_ids <- combined_list_study_ids[[i]][4]
  
  # commented this part out since I already did it:
  # write.table(sample_ids, file=file_name, quote=FALSE, sep='\t', row.names = FALSE)
}

# NOTE: The TSV files were transferred to Graham. I extracted the BIOM tables using my custom Bash scripts and transferred the BIOM tables to my desktop to conduct the rest of the analysis.

BIOM_files <- list.files("Data_Files/env/")[grep("*filtered.qza", list.files("Data_Files/env/"))]

# I sorted the names using code from the following link: https://stackoverflow.com/questions/17531403/how-to-sort-a-character-vector-where-elements-contain-letters-and-numbers-in-r
sorting_names <- as.numeric(gsub("[^[:digit:]]", "", BIOM_files))
names(sorting_names) <- seq_along(sorting_names)
BIOM_files <- BIOM_files[as.numeric(names(sort(sorting_names)))]
directory_env_files <- paste("Data_Files/env/", BIOM_files, sep = "")

all_otu_tables <- read_qza_files(directory_env_files)
lst_community_matrices <- assign_otu_m(all_otu_tables)

# Now, conducting the PERMANOVA tests. I made sure that they are all in order
groups <- list()
for (i in 1:length(combined_list_study_ids)){
  groups[[i]] <- combined_list_study_ids[[i]]$groupings
}

# Performing the PERMANOVA tests:
set.seed(5)
PERMANOVA_tests_env_associated <- compute_PERMANOVA_tests(groups, lst_community_matrices)
results <- lapply(PERMANOVA_tests_env_associated, function(x) print(x))

# I save the output into a TXT file using the function capture.output()
capture.output(results, file = "PERMANOVA_results_env_23Nov2020.txt")

# In this next section, I created dataframes of the samples with their associated PERMANOVA test results. I will be using those results along with their maximum distances for a logistic regression model. I printed the results and then created the dataframe
significance_of_tests <- bind_rows(PERMANOVA_tests_env_associated) %>% 
  na.omit() %>% 
  `rownames<-`(combined_study_ids) %>% 
  print() %>% 
  data.frame()

significance_of_tests

# There are 7 environmental samples that shows significance:
significance_of_tests %>% filter(Pr..F. < 0.05)

# the studies were: 1, 4, 7, 10, 22, 26, 33
which(significance_of_tests$Pr..F.<0.05)
combined_list_study_ids[which(significance_of_tests$Pr..F.<0.05)] %>% map(~.x %>% select(country) %>% table())

significance_of_tests <- check_significance(significance_of_tests)
dataframe_to_significance <- cbind(significance_of_tests, max_dist(combined_list_study_ids, combined_study_ids))
dataframe_to_significance$`Maximum Distance (Km)` <- as.numeric(dataframe_to_significance$`Maximum Distance (Km)`)
dataframe_to_significance$`Number of Records` <- as.numeric(dataframe_to_significance$`Number of Records`)

# Saving the data 
write_csv(dataframe_to_significance, "ENV_significance_dataframe_23Nov2020.csv")
write.table(combined_study_ids, "combined_study_ids_ENV_23Nov2020.txt", row.names = F, col.names = F, sep = ",")

# Plotting the significant results

combined_list_study_ids[which(significance_of_tests$Pr..F.<0.05)]
plot_type(combined_list_study_ids[which(significance_of_tests$Pr..F.<0.05)], combined_study_ids[which(significance_of_tests$Pr..F.<0.05)], "StudyID (Filtered)")

# testing ggmap

# Okay, it looks more complicated than I initially thought. I will need to get a google API key (https://developers.google.com/maps/documentation/maps-static/get-api-key). After I create the API key, I can use Google maps to plot where the samples were taken from in the world! So, that's great. I would definitely like to learn to do that. Also, check out https://cran.r-project.org/web/packages/ggmap/readme/README.html. The initial tutorial that I found: https://www.r-bloggers.com/2017/11/how-to-plot-basic-maps-with-ggmap/. 
library(ggmap)
?register_google
map_USA <- get_googlemap("USA", zoom = 12) %>% ggmap()

