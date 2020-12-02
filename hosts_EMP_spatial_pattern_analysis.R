# Teodora Tockovska
# Dec 1, 2020

# This script contains the code to perform spatial pattern analyses on host-associated from the Earth Microbiome Project.

# Loading my script which contains functions 
source("Functions.R")

# Read the host associated metadata
host_associated_microbiome_dataset <- read.csv(file = "host_associated_EMP.csv", sep = "\t")
nrow(host_associated_microbiome_dataset)

# Checking and removing NA values:
host_associated_microbiome_dataset <- check_and_remove_NAs(host_associated_microbiome_dataset)
nrow(host_associated_microbiome_dataset)

# Filter the dataset to extract the hosts with 30 or more records. I saved the dataframe of the host records into a variable. I chose "host_scientific_name" column instead of "host_species" because 882 hosts do not have a species names. Next, I viewed the dataframe and ordered it according to the decreasing records. The hosts have species names except for 4 hosts that only have their genus name. Homo sapiens (Humans) have the most records (3171) and Felis catus (domestic cats) have the fewest records (30).
hosts_with_30_records <- as.data.frame(table(host_associated_microbiome_dataset$host_scientific_name)) %>% filter(Freq>=30)
hosts_with_30_records[order(hosts_with_30_records$Freq, decreasing = T),]

# I created a list of dataframes called "lst_df_hosts_with_30_records". The dataframes represent the hosts' information that was indexed from the original EMP dataset. The for loop is used to loop through the hosts from "hosts_with_30_records" to extract the names. I used grep() to match the names so that I can extract the datasets. The dataframe will be used to extract the coordinates and sample IDs for the hosts.
lst_df_hosts_with_30_records <- list()
for (i in 1:nrow(hosts_with_30_records)){
  lst_df_hosts_with_30_records[[i]] <- host_associated_microbiome_dataset[grep(hosts_with_30_records$Var1[[i]],
                                                                  host_associated_microbiome_dataset$host_scientific_name),]
}

# The next goal is to plot the coordinate of each host to visualize their samples and determine outliers. However, there are hosts that were only sampled from 1 unique coordinate (Lat, Long). I call the function chosen_samples() that will check the unique coordinates per host and exclude hosts that were sampled from only 1 region (unique coordinate). The function will also plot the coordinates of the chosen hosts. The function also returns a list of the selected hosts ad prints out the number of selected hosts. There are 20 hosts that were selected based on the criteria.
greater_than_1_region <- chosen_samples(lst_df_hosts_with_30_records, hosts_with_30_records, type = "Host Genus/Species")

# Remove variables that are no longer needed:
rm(lst_df_hosts_with_30_records, hosts_with_30_records)

# Next, extract the coordinate dataframes per host into a list by using the function get_coordinates(). Using the list of coordinate dataframes, I check which UTM zones the samples were taken from. This is important because samples can be taken from all over the world and I would like to ensure that the samples that I am working with are not too far away from each other. 
coordinate_list <- get_coordinates(data = greater_than_1_region, lat = latitude_deg, long = longitude_deg, type = "Host")

# Next, I check whether the coordinates fall in the correct hemispheres. I use the custom function check_hemisphere() and pass 4 as an argument because that represents the "Country" column.
coordinate_list %>% map(~.x %>% check_hemisphere(4))

# There are errors in the first list. The country is USA but the hemispheres were recorded as NW and NE. USA is in the NW hemisphere, so that means the longitudes need to be negative. However, NE indicates that some longitudes were recorded as positive. I will change those.
coordinate_list[[1]][,1][which(coordinate_list[[1]][,1] > 0)] <- -coordinate_list[[1]][,1][which(coordinate_list[[1]][,1] > 0)]

# After I fixed the incorrect longitudes, I check the zones per coordinate dataframe and I create a dictionary that stores the UTM zone and the associated longitudes for each host. The dictionary is stored into a variable called "dict_UTM_zones".
dict_UTM_zones <- coordinate_list %>% map(~.x %>% group_by(Longitude) %>% check_zones())

# I view the lists of unique UTM zones per host by calling the function return_unique_UTM_zones(). It appears that many of the samples have multiple zones. For example, there are 7 zones from which homo sapiens were sampled (lst_uniq_UTM_zones[[7]]). This causes a problem when conducting the spatial pattern analysis because I cannot compare microbiomes of the hosts from such large distances. 
lst_uniq_UTM_zones <- return_unique_UTM_zones(dict_UTM_zones)

# How many hosts were sampled within 1 UTM zone? 12 hosts, whereas the rest (8) were sampled from multiple UTM zones. Next, I extract the hosts that were sampled from 1 UTM zone and put them into a list by getting their ids and indexing from "coordinate_list".
sum(lst_uniq_UTM_zones %>% map(~.x %>% length(.))==1)
single_zone_ids <- which(lst_uniq_UTM_zones %>% map(~.x %>% length(.))==1)
single_zone_coordinate_datasets <- coordinate_list[single_zone_ids]
lst_uniq_UTM_zones[single_zone_ids]

# I extracted the host species' dataframes of coordinates that were sampled from more than 1 UTM zone. 
multiple_UTM_zones <- coordinate_list[-single_zone_ids]
lst_uniq_UTM_zones[-single_zone_ids]

# Checking the names, which generates a tibble
host_species_names <- greater_than_1_region %>% map(~.x %>% group_by(host_scientific_name) %>% summarise(n=n()) %>% select(host_scientific_name))

# Hosts that were sampled from 1 UTM zone are the following: (Numbers in parentheses represent the total names recorded per host. The updated names per host are shown to the left of the parenthesis.)
single_zone_host_names <- extract_names(host_species_names, single_zone_ids, "host_scientific_name")
single_zone_host_names

# I want to split the hosts according to their UTM zones and then remove outliers. First, I want to extract the UTM zones that have more than 1 longitude for each host. I need to extract those hosts based on their UTM zones. I extracted the longitudes associated with UTM zones for each host using the function return_longitudes_per_UTM() which takes 2 things as parameters: the list of unique UTM zones for hosts and a dictionary of the UTM zones and longitudes. The variables contain the information for hosts sampled from multiple UTM zones. The function returns a list of longitudes. Each element belongs to a host and the longitudes belong to a specific UTM zone.
multiple_zones_lst <- lst_uniq_UTM_zones[-single_zone_ids]
dict_multiple_UTM_zones <- dict_UTM_zones[-single_zone_ids]
host_species_longitudes_lst <- return_longitudes_per_UTM(multiple_zones_lst, dict_multiple_UTM_zones, "Multiple")

# Using the function extract_longitudes(), I extracted the hosts with the target longitudes. The function takes 2 things as parameters: the list of coordinate dataframes per host and the list of longitudes within a specific UTM zone per host. The longitudes from the second parameter will be used to select the longitudes from the coordinates dataframes for each host. Hence, I will extract the coordinates of interest for each host.
updated_coordinates <- extract_longitudes(multiple_UTM_zones, host_species_longitudes_lst)
updated_coordinates <- rename_samples_list(updated_coordinates, "Host")

# I renamed the inner dataframes and removed empty dataframes, updating the original list of dataframes. 
output <- rename_inner_lists(updated_coordinates, "Host")
updated_coordinates <- output[[1]]
removed_hosts <- output[[2]]
rm(output)
updated_coordinates <- remove_empty(updated_coordinates)

# I fixed the names of the hosts using the function extract_names(). The function returns an updated list of names for hosts or environmental samples. This function's main purpose is meant to deal with multiple names for the host or environmental samples.
multiple_zones_names <- extract_names(host_species_names, -single_zone_ids, "host_scientific_name", "list")

# I kept track of the removed dataframes by using the function return_ids() and then I updated the names of the hosts in the variable multiple_zones_names.
output <- return_ids(removed_hosts, updated_coordinates, "Host")
multiple_zones_names <- change_all_multiple_samples_names(output, multiple_zones_names)

# Next, I check the greatest distances per region for the host species that had multiple UTM zones. It is important to note that hosts 5 and 6 in this dataframe didn't get filtered which is why the distances are greater for hosts 5 and 6. But overall, the distances are fairly manageable. I also print the names of the host species to remind myself what they were. The maximum distances were calculated using the function max_dist(). A chart is printed as output, where the host names are the row names, and the maximum distances are the columns. 
max_dist(updated_coordinates, multiple_zones_names)

# Using the chart from above, I can see that I need to filter outliers from the data. To see which hosts contain outliers, I plot the coordinates by using the function plot_type(). There are 3 parameters: the list of coordinate dataframes per host, the host names in order to the list of coordinates, and the data type which is set to "Host". 
plot_type(updated_coordinates, multiple_zones_names, "Host")

# After viewing the plots, I need to remove the outlier in host species 1 (Canis lupus familiaris) and remove species 3 (Felis catus) because their distances are very large and host species 3 doesn't have enough samples after filtering out the American samples. I will filter the 12 samples that had a latitude of -29 (that would leave 32 samples to work with). 
updated_coordinates[[1]] <- updated_coordinates[[1]][updated_coordinates[[1]]$Latitude < -30,]

# Removing the dataframes and the host species
updated_coordinates <- updated_coordinates[-3]
multiple_zones_names <- multiple_zones_names[-3]

# Viewing the changes by printing the maximum distance chart and plotting the coordinates:
max_dist(updated_coordinates, multiple_zones_names)
plot_type(updated_coordinates, multiple_zones_names, "Host")

# Adding the groupings for the PERMANOVA tests:
updated_coordinates <- updated_coordinates %>% map(~.x %>% unite("groupings", Longitude:Latitude, sep = "-", remove = F))

# Likewise, I check the greatest distances among the hosts that were sampled within 1 UTM coordinate by using the function max_dist(). The majority of the hosts were sampled within a close proximity. Note: the first host species shows a distance of 0 because I did not remove the host species after I fixed the longitude. I also plotted the hosts to view the outliers.
max_dist(single_zone_coordinate_datasets, single_zone_host_names)
plot_type(single_zone_coordinate_datasets, single_zone_host_names, "Host")

# Now, I will remove the outliers from the other host species. I need to remove host species 1 (Ateles geoffroyi) and 10 (Turdus merula). I need to filter out the 2 samples from host species 3 (Clamator glandarius). I am unsure if host species 5 has a mis-recording of 27 rabbit samples, but I will filter those out for now as well because the distance between the 27 rabbits and the other ones is about 8.3km.

# Filtering the outliers. After the filtering, I also printed the maximum distances and I plotted the coordinates.
single_zone_coordinate_datasets[[2]] <- single_zone_coordinate_datasets[[2]][single_zone_coordinate_datasets[[2]]$Latitude < 39.15,]
single_zone_coordinate_datasets[[3]] <- single_zone_coordinate_datasets[[3]][single_zone_coordinate_datasets[[3]]$Latitude < 37.6,]
single_zone_coordinate_datasets[[5]] <- single_zone_coordinate_datasets[[5]][single_zone_coordinate_datasets[[5]]$Latitude < -30,]
single_zone_coordinate_datasets[[6]] <- single_zone_coordinate_datasets[[6]][single_zone_coordinate_datasets[[6]]$Latitude < 37.25,]
single_zone_coordinate_datasets[[7]] <- single_zone_coordinate_datasets[[7]][single_zone_coordinate_datasets[[7]]$Latitude > 37.32,]
single_zone_coordinate_datasets[[8]] <- single_zone_coordinate_datasets[[8]][single_zone_coordinate_datasets[[8]]$Longitude > 145.4,]
single_zone_coordinate_datasets[[12]] <- single_zone_coordinate_datasets[[12]][single_zone_coordinate_datasets[[12]]$Longitude > 145.4,]

# Removing the dataframes and the host species
single_zone_coordinate_datasets <- single_zone_coordinate_datasets[-c(1,10)]
single_zone_host_names <- single_zone_host_names[-c(1,10)]

# Viewing the changes by printing the maximum distance chart and plotting the coordinates:
max_dist(single_zone_coordinate_datasets, single_zone_host_names)
plot_type(single_zone_coordinate_datasets, single_zone_host_names,type = "Host")

# Lastly, making the groupings for the PERMANOVA tests:
single_zone_coordinate_datasets <- single_zone_coordinate_datasets %>% map(~.x %>% unite("groupings", Longitude:Latitude, sep = "-", remove = F))

# How many host species do I have left? 17 host species
sum(length(single_zone_coordinate_datasets), length(updated_coordinates))

# I combined the lists so that I can do the PERMANOVA test
combined_list_host_species <- c(single_zone_coordinate_datasets, updated_coordinates)
combined_host_names <- c(single_zone_host_names, multiple_zones_names)

# Viewing the maximum distances
max_dist(combined_list_host_species, combined_host_names)

# Combined the rows to create one dataframe. The idea is so that I could save the dataframe into a CSV file. 
joined_hosts <- bind_rows(combined_list_host_species)

# Next, Saving the dataframe into a CSV file. It will be used in the script "Statistical_Models_and_Plotting.R".
# write_csv(joined_hosts, "dataframe_hostAssociated.csv")

# Conducting Spatial Pattern Analysis on 17 host species ####

# Before I do the PERMANOVA tests, I need to get the BIOM tables for the 11 host species. I will extract the sample IDs.So, I extract the sample IDs and saving them into files. I will transfer those files to Graham so that I can extract the BIOM tables.
for (i in 1:length(combined_list_host_species)){
  file_name <- paste("Host_Sample", as.character(i), "23Nov2020.tsv", sep = "_")
  sample_ids <- combined_list_host_species[[i]][4]
  
  # commented this part out since I already did it:
  # write.table(sample_ids, file=file_name, quote=FALSE, sep='\t', row.names = FALSE)
}

# NOTE: The TSV files were transferred to Graham. I extracted the BIOM tables using my custom Bash scripts and transferred the BIOM tables to my desktop to conduct the rest of the analysis.

# BIOM_files <- list.files()[grep("Data_files/*_filtered.qza",list.files())]
BIOM_files <- list.files("Data_Files/host/")[grep("*filtered.qza", list.files("Data_Files/host/"))]

# I sorted the names using code from the following link: https://stackoverflow.com/questions/17531403/how-to-sort-a-character-vector-where-elements-contain-letters-and-numbers-in-r
sorting_names <- as.numeric(gsub("[^[:digit:]]", "", BIOM_files))
names(sorting_names) <- seq_along(sorting_names)
BIOM_files <- BIOM_files[as.numeric(names(sort(sorting_names)))]
directory_BIOM_files <- paste("Data_Files/host/", BIOM_files, sep = "")

all_otu_tables <- read_qza_files(directory_BIOM_files)
lst_community_matrices <- assign_otu_m(all_otu_tables)

# Now, conducting the PERMANOVA tests. I made sure that they are all in order
groups <- list()
for (i in 1:length(combined_list_host_species)){
  groups[[i]] <- combined_list_host_species[[i]]$groupings
}

# Performing the PERMANOVA tests:
set.seed(5)
PERMANOVA_tests_host_associated <- compute_PERMANOVA_tests(groups, lst_community_matrices)
results <- lapply(PERMANOVA_tests_host_associated, function(x) print(x))

# I save the output into a TXT file using the function capture.output()
# capture.output(results, file = "PERMANOVA_results_host_23Nov2020.txt")

# In this next section, I created dataframes of the samples with their associated PERMANOVA test results. I will be using those results along with their maximum distances for a logistic regression model. I printed the results and then created the dataframe
significance_of_tests <- bind_rows(PERMANOVA_tests_host_associated) %>% 
  na.omit() %>% 
  `rownames<-`(combined_host_names) %>% 
  print() %>% 
  data.frame()

significance_of_tests

# There is only 1 host that shows significance:
significance_of_tests %>% filter(Pr..F. < 0.05)

# What sample types were taken from the significant host?
lapply(greater_than_1_region, function(x) {unique(x$env_material[grep("1696", x$study_id)])})
lapply(greater_than_1_region, function(x) {unique(x$country[grep("1696", x$study_id)])})

# Getting a better idea of the host species:
which(significance_of_tests$Pr..F.<0.05)
combined_list_study_ids[which(significance_of_tests$Pr..F.<0.05)] %>% map(~.x %>% select(country) %>% table())

# Creating a dataframe so that I can save it into a CSV file. The dataframe contains the PERMANOVA test results for each sample and the output from the function max_dist().
significance_of_tests <- check_significance(significance_of_tests)
dataframe_to_significance <- cbind(significance_of_tests, max_dist(combined_list_host_species, combined_host_names))
dataframe_to_significance$`Maximum Distance (Km)` <- as.numeric(dataframe_to_significance$`Maximum Distance (Km)`)
dataframe_to_significance$`Number of Records` <- as.numeric(dataframe_to_significance$`Number of Records`)

# Saving the data because they will be used in the script "Statistical_Models_and_Plotting.R".
# write_csv(dataframe_to_significance, "host_significance_dataframe_23Nov2020.csv")
# write.table(combined_host_names, "combined_host_names_23Nov2020.txt", row.names = F, col.names = F, sep = ",")
