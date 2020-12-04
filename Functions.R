# Teodora Tockovska
# Dec 1, 2020

# This script contains the functions for my main scripts that perform spatial patterns for host-associated microbiomes and free-living microbiomes.

# Loading libraries
library("vegan")
library("phyloseq")
library("qiime2R")
library("tidyverse")
library("sp")
library("rgdal")
library("geosphere")
library("fields")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("cowplot")

#' Remove rows that contain empty latitude and/or longitude coordinates. The parameter "dataset" is a dataframe that contained n samples (rows) and 76 observations (columns). 
#' 
#' @param dataset A dataframe of n samples (rows) and 76 observations (columns).
#' @return Updated dataframe
check_and_remove_NAs <- function(dataset){
  naLONG_ids <- which(is.na(dataset$longitude_deg))
  naLAT_ids <- which(is.na(dataset$latitude_deg))
  
  # Remove empty coordinates (LAT/LONG). 
  if(length(naLONG_ids) > 0){
    
    # The LAT and LONG row IDs are the same (contents and lengths).
    if((naLONG_ids == naLAT_ids) && (length(naLONG_ids) == length(naLAT_ids))){
      dataset <- dataset[-naLONG_ids,]
      print(paste("Removed ", length(naLONG_ids), " rows.", sep =""))
    }
    
    # Differing LAT and LONG lengths and row IDs
    else{
      original_len <- nrow(dataset)
      dataset <- dataset[-naLONG_ids,]
      naLAT_ids <- which(is.na(dataset$latitude_deg))
      dataset <- dataset[-naLAT_ids,]
      total_removed <- sum(original_len, -length(dataset))
      print(paste("Removed ", total_removed, " rows.", sep =""))
    }
    
    return(dataset)
  }
  
  else{
    print(paste(dataset, "contains 0 NA values.", sep = " "))
  }
}

#' Returns the hemisphere where the parameter "longlat" should be located in, and returns the Country that was recorded in the metadata. This function is used to validate the coordinates and remove or correct improperly recorded coordinates in the data. 
#' 
#' @param longlat Dataframe. Contains 1 row for a sample and 3 columns: Latitude, Longitude, and Country. The columns represent the location of sample extraction.
#' @return A list of strings. First element contains the predicted hemisphere in which the coordinate is found, and the second element contains the recorded Country. 
return_hemisphere <- function(longlat){
  
  if(longlat$Latitude > 0){
    if(longlat$Longitude > 0){
      hemisphere <- "NE"
    }
    
    else{
      hemisphere <- "NW"
    }
  }
  
  else{
    if(longlat$Longitude > 0){
      hemisphere <- "SE"
    }
    
    else{
      hemisphere <- "SW"
    }
  }
  
  return(c(hemisphere, longlat$country))
}

#' Checks the validity of the coordinates for each dataframe in the list "df_coordinates" by returning a list of predicted hemispheres in which the coordinates are located. It also returns the recorded Country name for the samples. The results from the list must be analyzed to determine whether the coordinates are located in the correct country and in the correct hemisphere. The "country_column" parameter must be reported because it varies between the host-associated and environmental datasets.
#' 
#' @param df_coordinates A list of dataframes. Each dataframe belongs to a study (environmental data) or host species/class (host-associated data). Each dataframe contains Latitude, Longitude, Country, Study ID, and sample ID columns. The environmental data also contains the environmental material column. 
#' @param country_column Int. A number that represents the country column in df_coordinates.
#' @return
check_hemisphere <- function(df_coordinates, country_column){
  hemisphere_check <- list()
  for(i in 1:nrow(df_coordinates)){
    
    # country_column is 6 for environmental data, and 4 for host-associated data.
    hemisphere_check[[i]] <- return_hemisphere(df_coordinates[i,c(1:2, country_column)])
  }
  return(unique(hemisphere_check))
}

#' Return a list of dataframes. This function loops through the list of host dataframes, host_scientific_name_30, that were extracted from EMP and the list of host names, species_30_records. The function plots the latitude and longitudes (coordinates) of each host, making sure that there are more than 1 latitude and longitude coordinates. This ensures that a host is sampled from more than 1 region. Hosts that were sampled from only 1 region (1 unique latitude and 1 unique longitude) will be ignored and not plotted. The hosts that were sampled from more than 1 region will be plotted and their dataset will be saved into a list. The list, called greater_than_1_region, will be filtered to remove empty elements, and then the list of datasets will be returned.
#' 
#' @param host_scientific_name_30 A list that contains n dataframes. The dataframes are host-specific and were extracted from the original EMP dataset. There are 76 columns in each dataframe. The rows of each dataframe differ because they reflect the number of records per host. 
#' @param species_30_records A dataframe that contains 2 columns and n rows. The columns, Var1 and Freq, refer to the host name (either the full species name or only the genus name) and the number of records, respectively. 
#' @param type A string to label the titles of the graphs.
#' @return A list of dataframes of hosts that were sampled from more than 1 unique coordinate (LAT/LONG).
chosen_samples <- function(host_scientific_name_30, species_30_records, type){
  greater_than_1_region <- list()
  for (i in 1:length(host_scientific_name_30)){
    
    # if more than 1 unique coordinate (LAT/LONG), plot the host and save the dataframe into the list.
    if (length(unique(host_scientific_name_30[[i]]$longitude_deg)) > 1 || length(unique(host_scientific_name_30[[i]]$latitude_deg)) > 1){
      title=paste(type, "#", as.character(i), ":", as.character(species_30_records$Var1[i]), sep = " ")
      plot(host_scientific_name_30[[i]]$longitude_deg, 
           host_scientific_name_30[[i]]$latitude_deg, 
           main = title,
           ylab = "Latitude",
           xlab = "Longitude")
      greater_than_1_region[[i]] <- host_scientific_name_30[[i]]
    }
  }
  # Removing the empty lists:
  greater_than_1_region <- greater_than_1_region[lengths(greater_than_1_region) != 0]
  
  print(paste("Selected", as.character(length(greater_than_1_region)), "cohorts that were sampled from more than 1 unique coordinate.", sep = " "))
  return(greater_than_1_region)
}

#' Return a community matrix from an OTU table. This function will take a phyloseq object (OTU table) and will convert the phyloseq object to a matrix. Then, the otu_matrix will be transformed so that the samples will be the rows and the OTUs will be the columns.
#' 
#' @param phyloseq_OTU_table An OTU table as a phyloseq object. 
#' @return The transformed OTU table of phyloseq_OTU_table.
create_comm_matrix <- function(phyloseq_OTU_table){
  otu_matrix <- as(phyloseq_OTU_table@.Data, "matrix")
  return(t(otu_matrix))
}

#' Returns the list containing the unique UTM zones within dictionary_of_UTM_zones. The parameter, dictionary_of_UTM_zones, contains lists that belong to each species. Within the lists are other lists that contain 2 elements for the UTM zones and for the longitudes. Within the function, the 
#' 
#' @param dictionary_of_UTM_zones 
#' @return 
return_unique_UTM_zones <- function(dictionary_of_UTM_zones){
  unique_zones <- list()
  for (i in 1:length(dictionary_of_UTM_zones)){
    zones_lst <- list()
    for(zone in 1:length(dictionary_of_UTM_zones[[i]])){
      if(!(dictionary_of_UTM_zones[[i]][[zone]] %in% zones_lst)){
        zones_lst[[zone]] <- (dictionary_of_UTM_zones[[i]][[zone]]$Zone)
      }
    }
    unique_zones[[i]] <- unique(unlist(zones_lst))
  }
  return(unique_zones)
}

#' Return a dataframe that contains the coordinates and the collection of sampling.
#' 
#' @param data The dataset or a list of many datasets of interest that comes from EMP.
#' @param lat The latitude column name that is present in data.
#' @param sampleid The sample ID (first column)
#' @param long The longitude column name that is present in data.
#' @param type Str. The description of the data.
#' @return A dataframe of 5 columns and nrows(data). The columns include: latitude, longitude, year of sampling, month of sampling, and day of sampling.
#' @examples 
#' get_coordinates(birds, latitude_deg, longitude_deg, "study_id")
get_coordinates <- function(data, lat, long, type){
  
  if(class(data) == "list"){
    coordinate_lst <- list()
    for(i in 1:length(data)){
      
      if(type == "Host" || type == "host"){
        coords <- data.frame(Longitude=data[[i]]$long, 
                             Latitude=data[[i]]$lat, 
                             sampleid=data[[i]]$X.SampleID,
                             country=data[[i]]$country,
                             study_id=data[[i]]$study_id)
      }
      
      else if(type == "Study_ID" || type == "study_id" || type == "Study_id"){
        coords <- data.frame(Longitude=data[[i]]$long, 
                             Latitude=data[[i]]$lat, 
                             sampleid=data[[i]]$X.SampleID, 
                             study_id=data[[i]]$study_id, 
                             env_material=data[[i]]$env_material, 
                             country=data[[i]]$country)
      }
      
      coordinate_lst[[i]] <- coords
    }
    return(coordinate_lst)
  }
}

#' Returns the unique UTM zones for coord_df. It creates a dictionary (or a list of many nested lists). The first list holds the information on each sample. For example, each sample contains a list of 2 elements, the first element holds the UTM zone integer, and the second element contains the longitude integer belonging to the UTM zone integer. 
#' 
#' @param coord_df Dataframe which contains the Latitude in the first column and Longitude in the second column
#' @return Dictionary; list of lists. It contains the UTM zones and the Longitudes
check_zones <- function(coord_df){
  uniq_coords <- unique(coord_df$Longitude)
  dict <- vector(mode = "list", length = length(uniq_coords))
  
  for (i in 1:length(uniq_coords)){
    names(dict)[[i]] <- paste("UTM",  "Zone", i, sep = "_")
    
    dict[[i]][[1]] <- floor((uniq_coords[[i]] + 180) / 6) + 1
    names(dict[[i]])[[1]] <- "Zone"
    
    dict[[i]][[2]] <- uniq_coords[[i]]
    names(dict[[i]])[[2]] <- "Longitude"
  }
  
  return(dict)
}

#' Reads the QZA QIIME2 artifacts (QZA files) into R and converts them into OTU tables. 
#' 
#' @param list_of_BIOM A list that contains the names of the QZA files that are located in the working directory. The files must be QZA QIIME2 artifact files. The names can contain the path to the files.
#' @return A list of OTU tables.
read_qza_files <- function(list_of_BIOM){
  otu_tables_lst <- list()
  for (i in 1:length(list_of_BIOM)){
    otu_tables_lst[[i]] <- qza_to_phyloseq(features = list_of_BIOM[[i]])
  }
  return(otu_tables_lst)
}

#' Return a list of community matrices. This function takes in a list of OTU tables. The OTU tables list object is generated when reading the QZA artifacts into R and saving the variables into a list. Then, the function will loop through the list to convert each OTU table into a OTU matrix. It uses the function create_comm_matrix(). 
#' 
#' @param lst_of_otu_tables A list object of OTU tables.
#' @return A list of community matrices.
assign_otu_m <- function(lst_of_otu_tables){
  lst_of_community_matrices <- list()
  for (i in 1:length(lst_of_otu_tables)){
    lst_of_community_matrices[[i]] <- create_comm_matrix(lst_of_otu_tables[[i]])
  }
  return(lst_of_community_matrices)
}

#' Conduct PERMANOVA tests for the df_groupings using comm_matrix, and return the test.
#' 
#' @param df_groupings A list of dataframes that contain the groupings, latitude, longitude, and sample ID of a host or environmental sample.
#' @param comm_matrix A list of community matrices that are in the same orders as the hosts or environmental samples as in df_groupings list.
#' @return A list of PERMANOVA tests
compute_PERMANOVA_tests <- function(df_groupings, comm_matrix){
  Ptests <- list()
  for(i in 1:length(df_groupings)){
    print(paste("Started PERMANOVA for Sample", i, sep = " "))
    Ptests[[i]] <- adonis2(decostand(comm_matrix[[i]], "hel") ~ ., data = as.data.frame(df_groupings[[i]]))
    print(paste("Finished PERMANOVA test for Sample", i, sep = " "))
  }
  
  print("Printing Results")
  return(Ptests)
}

#' HELPER FUNCTION. Return a list of longitudes that are specific to a UTM zone, zone_int. The function loops through the UTM_zone_lst and matches the zone in the list to zone_int, and will add the longitudes that belong to the matched UTM zone, and that list is returned. This function is called in longitudes_uniq_dict().
#' 
#' @param zone_int An integer that represents the UTM zone.
#' @param UTM_zone_lst A list that belonging to a sample (host-associated or environmental dataset) that contain a series of lists of UTM zones in the first elements and longitudes in the second elements.
#' @return A list containing the longitudes that are associated to the specified UTM zone.
get_longitudes <- function(zone_int, UTM_zone_lst){
  longitudes_per_UTM <- list()
  for(zone_lst in 1:length(UTM_zone_lst)){
    if(UTM_zone_lst[[zone_lst]]$Zone == zone_int){
      longitudes_per_UTM[[zone_lst]] <- UTM_zone_lst[[zone_lst]]$Longitude
    }
  }
  return(unique(unlist(longitudes_per_UTM)))
}

#' HELPER FUNCTION. Returns the vector of longitudes that are associated with a list of UTM zones. It is called in return_longitudes_per_UTM().
#' 
#' @param lst_of_UTM_zones A list of all UTM zones belonging to a host or environmental sample from the dictionary, dict_UTM_zones.
#' @param dict_UTM_zones A dicitonary that contains the UTM zones and longitudes for a host or environmental sample. The dictionary can contain several lists that belong to multiple cohorts, or the dictionary can just be a list of UTM zones and longitudes for one cohort.
#' @return A list of longitudes that are associated with a list of UTM zones.
longitudes_uniq_dict <- function(lst_of_UTM_zones, dict_UTM_zones){
  new_lst <- list()
  for(i in 1:length(lst_of_UTM_zones)){
    
    # this function will be called, where I pass on the integer -- specific UTM zone, and the list of zones to get the longitudes
    longitudes <- get_longitudes(lst_of_UTM_zones[[i]], dict_UTM_zones)
    if(length(longitudes) == 1){
      next
    }
    new_lst[[i]] <- longitudes
    
  }
  # return(unlist(new_lst))
  return(new_lst)
}

#' Return a list of several lists that contain longitudes per UTM zone. 
#' 
#' @param lst_of_UTM_zones A list of all UTM zones belonging to a host or environmental sample from the dictionary, dict_UTM_zones.
#' @param dict_UTM_zones A dicitonary that contains the UTM zones and longitudes for a host or environmental sample. The dictionary can contain several lists that belong to multiple cohorts, or the dictionary can just be a list of UTM zones and longitudes for one cohort.
#' @param type Str or missing value. It indicates whether the dictionary contains lists for multiple cohorts or for one cohort. If type is not included as a parameter, then that means the dictionary is species/sample specific (i.e. only contains lists of the UTM zones and longitudes for one host species or one environmental sample). If the type parameter is included when the function is called (it can be anything but should be a string), then that indicates that the dictionary contains nested lists, where the first set of lists represent hosts or environmental samples, and their inner lists contain the UTM zones and longitudes.
#' @return A list containing the longitudes for a specified UTM zone.
return_longitudes_per_UTM <- function(lst_of_UTM_zones, dict_UTM_zones, type){
  
  # The dictionary/list belongs to only 1 species or environmental sample
  if(missing(type)){
    return(longitudes_uniq_dict(lst_of_UTM_zones, dict_UTM_zones))
  }
  
  # The dictionary contains lists that belong to multiple cohorts.
  else{
    nested_list <- list()
    
    # loop through each host species or environmental sample, extracting the species 
    for(i in 1:length(dict_UTM_zones)){
      nested_list[[i]] <- longitudes_uniq_dict(lst_of_UTM_zones[[i]], dict_UTM_zones[[i]])
    }
    return(nested_list)
  }
}

#' Returns a list of dataframes of n rows and 3 columns (Longitude, Latitiude, and sampleid). Each dataframe belongs to a host or environmental sample. This function will index the dataframes using the indices of the longitudes of interest from target_longitudes. It prints out which host or environmental sample was changed within the list of dataframes, multiple_UTM_zones. The samples that were not changed are ultimately not filtered. 
#' 
#' @param multiple_UTM_zones A list that contains dataframes of 3 columns (Longitude, Latitude, and sampleid) per host/environmental sample.
#' @param target_longitudes A list of longitude numerical vectors per host/environmental sample.
#' @return A list of dataframes per host/environmental sample.
extract_longitudes <- function(multiple_UTM_zones, target_longitudes){
  for(host in 1:length(target_longitudes)){
    matched <- list()
    
    # skips samples that don't have longitudes (len=0)
    if(length(target_longitudes[[host]]) == 0){
      next
    }
    
    print(paste("Extracted coordinates for sample/host", as.character(host)))
    
    for(longitudes in 1:length(target_longitudes[[host]])){
      
      matched <- which(multiple_UTM_zones[[host]]$Longitude %in% target_longitudes[[host]][[longitudes]])
      
      # sort the indices that matched:
      matched <- sort(unlist(matched))
      
      if(length(target_longitudes[[host]]) == 1){
        
        # extracting and replacing the coordinates from the dataframe per host species and extracting the rows. 
        target_longitudes[[host]] <- multiple_UTM_zones[[host]][matched,]
      }
      
      else{
        target_longitudes[[host]][[longitudes]] <- multiple_UTM_zones[[host]][matched,]
        
      }
    }
  }
  
  return(target_longitudes)
}

#' Return the list of renamed list elements
#' 
#' @param coordinate_lst_of_dfs A list of dataframes that contain the coordinates of samples
#' @param type Str. A description of the sample types. 
#' @return A list of dataframes 
rename_samples_list <- function(coordinate_lst_of_dfs, type, version = F){
  for(sample in 1:length(coordinate_lst_of_dfs)){
    
    if(version){
      names(coordinate_lst_of_dfs)[[sample]] <- paste(type, sample, sep = ".")
    }
    
    else{
      names(coordinate_lst_of_dfs)[[sample]] <- paste(type, sample, sep = "_")
    }
    
  }
  return(coordinate_lst_of_dfs)
}

#' HELPER FUNCTION. Some times there are hosts that do not have more than 2 listed names so this function is used to deal with the multiple names. A list of names are provided into the function as a parameter and the function iterates through the list. If there is only 1 host name (which means it could be a class or genus name), the function skips that. Once the function sees a name that is > 1 length (meaning it has a species name) it will change the lst_name to the name and it returns the one character vector.
#' 
#' @param lst_names A list of host or environmental sample names. 
#' @return A character vector.
pick_host_name <- function(lst_names){
  for(i in 1:length(lst_names)){
    if(length(unlist(lst_names[[i]])) > 1){
      lst_names <- paste(lst_names[[i]], collapse = " ")
      break
    }
  }
  return(lst_names)
}

#' Returns an updated list of names for hosts or environmental samples. This function's main purpose is meant to deal with multiple names for the host or environmental samples. The function loops through the list of names, extracts the names from each tibble, converts the names into a character vector and calls another function, pick_host_name(), that will choose a name. The chosen name will overwrite the current host or environmental name at index i. The function loops through all host or environmental samples to return a list of updated names. It also saves the updated name followed by total number of names recorded for the sample in parentheses.
#' 
#' @param names_lst A tibble list that contains unique names of hosts or environmental samples. 
#' @param indices An integer vector of indices for samples. These indices are used to extract the names from names_lst by indexing.  
#' @param column A target column to access the tibble. The tibble has columns for each name so using the name of the column provides a way to access (index) each tibble to get the list of names.
#' @param type Empty or a string. Represents what the function should return
#' @return A character vector or a list of the host or environmental names that were updated.
extract_names <- function(names_lst, indices, column, type){
  
  # For cases that indexing must happen to extract names:
  if(!(missing(indices))){
    names <- names_lst[indices]
  }
  
  else{
    names <- names_lst
  }
  
  for(i in 1:length(names)){
    
    # the column_name could either be "host_scientific_name" or "env_material"
    if(length(str_split(names[[i]][[column]], " ")) > 1){
      sample_names <- str_split(names[[i]][[column]], " ")
      
      # saves the updated name followed by total number of names recorded for the sample in parentheses.
      names[[i]] <- paste(pick_host_name(sample_names), " num", as.character(length(sample_names)), "", sep = "")
    }
    colnames(names[[i]]) <- NULL
  }
  
  if(missing(type)){
    return(unlist(names))
  }
  
  else if(type == "list"){
    return(names)
  }
}

#' Change the inner lists' names of "lst_of_dfs". The dataframes names contain the "type" ID ("host" or "StudyID"), followed by the study ID and environmental type (environmental data) or host species names (host-associated data). The names of the lists are changed if they are greater than 0.
#' 
#' @param lst_of_dfs A list that contains dataframes or lists. Each dataframe/list of dataframes belong to a host species or a study (environmental samples). There could be empty lists/dataframes which will be removed.
#' @param type Str. Describes the type of data that is being handled. Either "Host" or "StudyID".
#' @return List. The first element contains the dataframe with updated element names (versions added to the dataframes). The second element contains the removed dataframes.
rename_inner_lists <- function(lst_of_dfs, type){
  new_lst <- list()
  for(host in length(lst_of_dfs):1){
    if(class(lst_of_dfs[[host]]) == "list"){
      new_lst <- append(new_lst, names(lst_of_dfs)[[host]])
      if(length(lst_of_dfs[[host]]) > 0){
        lst_of_dfs[[host]] <- rename_samples_list(lst_of_dfs[[host]], paste(type, as.character(host), sep = "_"), version = T)
        lst_of_dfs <- append(lst_of_dfs, lst_of_dfs[[host]], host)
      }
    }
  }
  new_lst <- unlist(new_lst)
  
  return(list(lst_of_dfs, new_lst))
}

#' Remove all the lists within updated_coordinates and remove any empty dataframes. Return the updated dataframe.
#' 
#' Note: I iterated the list backwards because the ordering of the elements wouldn't change when I updated the list when I removed the empty elements. Iterating backwards perserves the order of the list. I used this link to help me -- https://stackoverflow.com/questions/50403752/r-remove-list-full-of-na-from-a-list-of-lists.
#' 
#' @param lst_of_dfs
#' @return Updated list of dataframes after the empty dataframes were removed.
remove_empty <- function(lst_of_dfs){
  lst_of_dfs <- Filter(function(x) !all(class(x) == "list"), lst_of_dfs)
  
  for(i in length(lst_of_dfs):1){

    if(nrow(lst_of_dfs[[i]]) == 0){
      print(paste("Removed", names(lst_of_dfs[i]), sep = " "))
      lst_of_dfs <- lst_of_dfs[-i]
    }
  }
  
  return(lst_of_dfs)
}

#' HELPER FUNCTION. Adds version numbers to the host_name using num_of_subsets parameter. Returns a list of sample names with their version numbers.
#' 
#' @param host_name Str. The target sample name that needs to be given version numbers.
#' @param num_of_subsets Int. The number of dataframes that belong to host_name.
#' @return A list of sample names that contain version numbers.
extend_names_of_hosts <- function(host_name, num_of_subsets){
  names <- list()
  for(i in 1:num_of_subsets){
    sentence <- paste(host_name, as.character(i), sep = " ")
    names <- append(names, sentence)
  }
  return(names)
}

#' A list of dataframes will be edited to remove a set of elements that were named. Those elements were saved into a vector called samples_vector. This vector will be used along with the lst_of_dfs to determine how many dataframes belong to the non-empty dataframes and also to keep track of the empty dataframes which used to be part of the list of dataframes. The parameter, type, represents the type of samples that are being used. The type parameter will be used to name the new set of elements part of a list. 
#' 
#' A new list will be returned that will contain the removed samples from the list of dataframes. The first element of the list will be renamed using the type parameter, following by the ID, such as "Host_4". The first element of the new list contains the nonzero_IDs, which means that the listed IDs of the samples contained more than 1 dataframe of information. The second element contains the subset_length_IDs, which are the IDs that belong to the list in the 4th element (to be explained). The third element contains the IDs of the empty dataframes. These samples contained no information, hence their removal. Finally, the fourth element in the list contains the subset_length. It is a list of each removed sample along with their lengths of dataframes. The 2nd element (subset_length_IDs) will use this information to extract the lengths per sample.
#' 
#' The new list containing these 4 elements of information will be returned.
#' 
#' @param samples_vector A vector of characters that represent sample names that will be removed from a list of dataframes.
#' @param lst_of_dfs A list of dataframes that will be used to count the number of dataframes per sample in samples_vector.
#' @param type Str. A character that represents the type of data (either "Host" for host-associated or "StudyID" for environmental).
return_ids <- function(samples_vector, lst_of_dfs, type){
  subset_length <- list()
  nonzero_ids <- list()
  subset_length_ids <- list()
  empty_lists <- list()
  
  for(i in 1:length(samples_vector)){
    subset_length[[i]] <- length(grep(samples_vector[[i]], names(lst_of_dfs)))
    names(subset_length)[[i]] <- samples_vector[[i]]
    
    if(subset_length[[i]] > 0){
      nonzero_ids <- append(nonzero_ids, as.integer(gsub(paste(type, "_", sep = ""), "", samples_vector[i])))
      subset_length_ids <- append(subset_length_ids, i)
    }
    
    else{
      empty_lists <- append(empty_lists, as.integer(gsub(paste(type, "_", sep = ""), "", samples_vector[i])))
    }
  }
  
  # Saving variables
  nonzero_ids <- unlist(nonzero_ids)
  subset_length_ids <- unlist(subset_length_ids)
  empty_lists <- sort(unlist(empty_lists), decreasing = F)
  
  # Creating a list and renaming the lists based on their contents.
  lists_to_return <- list(nonzero_ids, subset_length_ids, empty_lists, subset_length)
  names(lists_to_return) <- c("nonzero_ids", "subset_length_ids", "empty_lists", "subset_length")
  
  return(lists_to_return)
}

#' HELPER FUNCTION. Updates the character vector of names by adding the version number to each name. The updated names are returned. 
#' 
#' @param lst_of_lengths A list of lengths. The elements of the list contains the sample names that have some integer that represents the lengths of dataframes.
#' @param lst_of_host_names The associated sample names belonging to lst_of_lengths.
#' @return Character vector of updated sample names with the version numbers.
create_extended_host_names <- function(lst_of_lengths, lst_of_host_names){
  if(length(lst_of_host_names) == length(lst_of_lengths)){
    for(i in 1:length(lst_of_host_names)){
      
      # change the names of the samples by adding the versions to the strings based on the lengths 
      lst_of_host_names[[i]] <- extend_names_of_hosts(lst_of_host_names[[i]], lst_of_lengths[[i]])
    }
  }
  return(lst_of_host_names)
}

#' HELPER FUNCTION. Updates the character vector of names by inserting the updates sample names with their version numbers to the non-zero IDs. The original sample names are removed after the new names are inserted. The updated character vector is returned. 
#' 
#' @param names_vector A character vector of sample names.
#' @param lst_of_host_names A list of updated sample names with version numbers.
#' @param nonzero_ids An int vector that contains the indices of sample names that need to be changed. 
replace_ids_with_changed_names <- function(names_vector, lst_of_host_names, nonzero_ids){
  for(i in 1:length(nonzero_ids)){
    for(j in length(names_vector):1){
      if(nonzero_ids[i] == j){
        names_vector <- append(names_vector, lst_of_host_names[[i]], j)
        names_vector <- names_vector[-j]
      }
    }
  }
  return(names_vector)
}

#' This is the main function that will perform the removal of the sample names that contained empty lists, but also change the names of the list elements for the samples that had multiple dataframes. The parameter output_list contains a list of 4 elements with "nonzero_ids", "subset_length_ids", "empty_lists", and "subset_length". The first element of the list contains the numbers of list elements in the main list of dataframes that contained a list of dataframes (instead of 1) which corresponds to the subsets of the same sample type but based on its UTM zone or location. These lists of dataframes need to be flattened out, and those flattened lists of dataframes need to be renamed. The subset_length_ids element of output_list contains the index positions of the non-zero IDs in the subset_length list. The IDs were provided to easily access those IDs (to find their samples and lengths of dataframes). Finally, the empty_lists element contains the samples that had 0 lengths, which means that those samples need to be removed. 
#' 
#' This function uses many helper functions which would change the names of the samples and returns the changed list of names. Why was this function created? Essentially, there was a list of dataframes and some of them were nested. The nested lists of dataframes belong to 1 type of sample or host. They were nested because the samples were taken from multiple UTM zones, hence different dataframes of coordinates. The nested lists of dataframes were flattened, and the original nested list of dataframes was removed, but the new dataframes that belonged to 1 sample type were renamed with versions to keep track of how many dataframes belong to 1 sample (i.e. studyid_100.2 means that for study ID 100, version 2 contains coordinates). Some dataframes were empty which must be removed. This function was created to handle one side: the character vector of the sample names. The function keeps track of which study IDs or host species had multiple dataframes in order to record their verisons to the host species names or study IDs. 
#' 
#' @param output_list List of 4 elements containing information on the lists of dataframes that had either empty or more than 1 dataframe.
#' @param multiple_zones_names A character vector of sample names that are in the exact order of the list of dataframes (of coordinates) and where each dataframe belongs to one sample.
#' @return Str vector. The updated names of the samples. 
change_all_multiple_samples_names <- function(output_list, multiple_zones_names){
  
  # Retrieve the samples that had mutliple dataframes. The IDs represent the positions of the nested lists of dataframes. Also, retrieving the lengths of dataframes that belong to each nested list of dataframes.
  change_host_names <- multiple_zones_names[output_list$nonzero_ids]
  hosts_info_lengths <- output_list$subset_length[output_list$subset_length_ids]
  
  # Add version numbers fo the names of the samples using the helper function. Then, extract IDs of the names that need to be removed. The IDs reflect the sample IDs for the main list of dataframes (of coordinates) per samples.
  change_host_names <- create_extended_host_names(hosts_info_lengths, change_host_names)
  hosts_to_remove <- multiple_zones_names[output_list$empty_lists]
  
  # Replace the orignal names with the updated names, remove the names that belong to empty dataframes. The changes to the names results into a list, so use unlist() to generate a character vector of the updated names. 
  multiple_zones_names <- replace_ids_with_changed_names(multiple_zones_names, change_host_names, output_list$nonzero_ids)
  multiple_zones_names <- multiple_zones_names[-which(multiple_zones_names %in% hosts_to_remove)]
  multiple_zones_names <- unlist(multiple_zones_names)
  
  return(multiple_zones_names)
}

#' Prints a dataframe of the maximum distance between records for each host or environmental sample. The rows represent the sample. There are 2 columns, where the first column contains the names of the samples and the second column represents the maximum distances.
#' 
#' @param coordinate_df A list of dataframes that contain the longitude, latitude, and sampleid for each host or environmental sample.
#' @param host_names A character vector of host or environmental sample names.
#' @return Null
max_dist <- function(coordinate_df, host_names){
  print(coordinate_df %>% 
          map(~.x %>% select(Longitude, Latitude) %>% distm(fun = distHaversine) %>% max()/1000) %>% 
          data.frame() %>% 
          t() %>% 
          cbind(coordinate_df %>% map(~.x %>% nrow(.))) %>% 
          `rownames<-` (host_names) %>% 
          `colnames<-` (c("Maximum Distance (Km)", "Number of Records"))
  )
}

#' Plots the longitude and latitude coordinates for host or environmental samples.
#' 
#' @param coordinate_df A list of dataframes that contain the longitude, latitude, and sampleid for each host or environmental sample.
#' @param host_names A character vector of host or environmental sample names.
#' @return Null
plot_type <- function(coordinate_df, host_names, type){
  for (i in 1:length(coordinate_df)){
    save_plot <- ggplot(data = coordinate_df[[i]]) +
      geom_point(aes(x = Longitude, y = Latitude)) +
      ggtitle(paste(type, as.character(i), ":", as.character(host_names[i]), sep = " "))
    
    print(save_plot)
  }
}

#' A dataframe of the combined PERMANOVA tests is required as a parameter, called "df_PERMtests". This functions adds a column for the significant values by using the p-values from the column "Pr..F.". 0 is an insignificant p-value and 1 is a significant p-value.
#' 
#' @param df_PERMtests Dataframe that contains n rows (samples) and p columns (PERMANOVA test results). A column called "Pr..F." is required.
#' @return Updated dataframe.
check_significance <- function(df_PERMtests){
  df_PERMtests$Significance <- 0
  df_PERMtests$Significance[which(significance_of_tests$Pr..F.<0.05)] <- 1
  return(df_PERMtests)
}

#' HELPER FUNCTION. Environmental material column is used within each dataframe in the list of dataframes, "studyID_df". All materials are counted and saved in a variable. Materials with records > 15 will be saved. A list is created that contains dataframes for each saved environmental materials and the element names are renamed to have the data type (i.e. "StudyID"), the ID followed by the dataframe number or "version". The version indicates the number of dataframe that belongs to a specific study. This function is called in subsetting_by_env_materials(). 
#' 
#' @param title The name of the current list element (names(studyID_df[x]), where x is the index when this function was called)
#' @param studyID_df A list of dataframes that contain coordinate information per study. The environmental material column is required.
#' @return A list of updated dataframes that contain dataframes per study. The dataframes are environmental material-specific. Each element is named using the title (dataType_num.version).
update_env_materials <- function(title, studyID_df){
  updated_studyID_df <- list()
  table_env_types <- studyID_df$env_material %>% table() %>% as.data.frame()
  env_types <- table_env_types$.[which(table_env_types$Freq >= 15)] %>% as.character()
  
  if(length(env_types) > 1){
    for (i in 1:length(env_types)){
      updated_studyID_df[[i]] <- studyID_df[which(studyID_df$env_material %in% env_types[i]),]
      names(updated_studyID_df)[i] <- paste(title, env_types[i], sep = "_")
    }
  }
  
  return(updated_studyID_df)
}

#' This is the main function that extracts the multiple environmental material rows and creates new dataframes. Hence, the list of dataframes, "studyID_df", grows larger. The helper function, update_env_materials(), is required to perform the extraction of the environmental material rows and compile them into their own dataframes. 
#' 
#' @param studyID_df A list of dataframes that contain coordinate information per study. The environmental material column is required.
#' @return An updated dataframe. Studies with multiple environmental types will be paritioned into several dataframes for each material type.Version numbers will also be added in the names of the list elements. The old dataframes (that were paritioned) will be removed.
subsetting_by_env_materials <- function(studyID_df){
  env_material_lst <- vector(mode = "list", length = length(studyID_df))
  for(i in length(studyID_df):1){
    env_material_lst[[i]] <- update_env_materials(names(studyID_df[i]), studyID_df[[i]])
  }
  
  for(j in length(env_material_lst):1){
    if(length(env_material_lst[[j]]) > 1){
      studyID_df <- append(studyID_df, env_material_lst[[j]], j)
      studyID_df <- studyID_df[-j]
    }
  }
  
  return(studyID_df)
}

#' Renames the list names in "studyID_df" if there are no environmental materials detected. So, first the function detects list names that do not have environmental materials, retrieves the lists' IDs, and then in the for loop, gets those environmental materials and concatenates them to the current list names. The list names are updated. The updated list of dataframes is returned.
#' 
#' @param studyID_df A list of dataframes that contain coordinate information per study. The environmental material column is required.
#' @return An updated dataframe where the list names have been modified.
change_element_names <- function(studyID_df){
  edit_names_IDs <- grep("*[0-9+]$", names(studyID_df))
  get_names <- studyID_df[edit_names_IDs]
  
  for(i in 1:length(get_names)){
    names(get_names)[i] <- paste(names(get_names[i]), get_names[[i]]$env_material[1], sep = "_")
  }
  
  names(studyID_df)[edit_names_IDs] <- names(get_names)
  
  return(studyID_df)
}

#' Creates a world map and plots all longitude and latitude coordinates of samples. ggplot is used to build the figure, as well as a few variables that are required to assemble the world map. The "sample_df" parameter is the necessary dataset which contains n rows (for each sample), and longitude ("long" parameter), latitude ("lat" parameter) columns. The "title" parameter is required to add a title to the figure, "colour_id" is required to input a colour for the points, and "shape" is required for selecting the shape for the points. "tiff" is set to false by default, but it could be set to true if the user wants to save the figure as a tiff figure. Lastly, "file_name" is set to Null by default, but it is required only if tiff is set to true. The parameter "file_name" is for naming the figure to save to desktop. Note that double curly brackets, {{}}, wrapped around the long and lat parameters in order to use the direct column names. This link was used to accomplish that: https://stackoverflow.com/questions/2641653/pass-a-data-frame-column-name-to-a-function
#' 
#' @param sample_df Dataframe that contains longitude and latitudes. 
#' @param colour_id Str. A colour code to colour the points on the world map.
#' @param long Int vector. Longitude column in the dataframe, sample_df.
#' @param lat Int vector. Latitude column in the dataframe, sample_df.
#' @param title Str. The title to the world map.
#' @param shape Int. The number to set the points on the world map.
#' @param tiff Boolean. Default set to False. Set to True to save the figure as a tiff.
#' @param file_name Str. Default set to Null. Required if tiff is set to True. Names the tiff file that you save to desktop.
#' @return Null
plot_world_map <- function(sample_df, colour_id, long, lat, title, shape, tiff = FALSE, file_name = NULL){
  world <- ne_countries(scale = "medium", returnclass = "sf")
  theme_set(theme_bw())
  world_plot <- ggplot(data = world) + 
    geom_sf() + 
    geom_point(data = sample_df, aes(x = {{long}}, y = {{lat}}), size = 2, shape = shape, fill = colour_id) + 
    labs(title = title, x = "Longitude", y = "Latitude") + 
    theme(text = element_text(size = 10), 
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  if(tiff && class(file_name) == "character"){
    tiff(paste(file_name, "tiff", sep = "."), units="in", width = 6, height = 3.5, res = 300)
    print(world_plot)
    dev.off()
  }
  
  print(world_plot)
}

#' This function takes two strings as parameters which represent 2 files: file1 and file2. They will be read and the contents will be saved into a list. The list is returned.
#'
#' @param file1 Str. The file name and prefix (i.e. txt) that contains the names of the host species/study IDs
#' @param file2 Str. The file name and prefix that contains a dataframe that contains the PERMANOVA results in the first 4 columns (including the P-values), a Significance column that contains either 0's or 1's, and 2 final columns for Maximum Distance (km) and Number of Records. The rows are the hosts or studies.
#' @return List. The first element contains the list of strings from file1, and element 2 contains the dataframe from file2.
read_files <- function(file1, file2){
  return_lst <- list()
  sample_names <- read.table(file1) %>% unlist() %>% as.character()
  significance_dataframe <- read_csv(file2) %>% data.frame() %>% `rownames<-`(sample_names)
  return_lst[[1]] <- sample_names
  return_lst[[2]] <- significance_dataframe
  return(return_lst)
}

# For the chi-squared tests, I will need to create a table. I created 2 helper functions and 1 main function to build the confusion matrix.

#' HELPER FUNCTION. Creates and returns a 2x2 confusion matrix. The rows represent the data type, "Env" and "Host". The columns represent the PERMANOVA p-value results, "Significant" and "Insignificant". This function is called in create_confusion_matrix().
#' 
#' @return 2x2 confusion matrix.
make_table <- function(){
  conf_matrix <- matrix(0, nrow = 2, ncol = 2)
  rownames(conf_matrix) <- c("Env", "Host")
  colnames(conf_matrix) <- c("Significant", "InSignificant")
  return(conf_matrix)
}

#' HELPER FUNCTION. The sample_df parameter represents a dataframe that must contain a Significance column with 0's and 1's. The dataframe must have n rows of samples. Samples that have a 0 in the Significance column represent an insignificant p-value result, and samples with 1 in the column represent a significant p-value result (from the PERMANOVA tests). The function will fill return a list of 2 elements. The first element contains the total number of significant results, and the second element contains the number of insignificant results. For example, if there are 15 samples with 3 1's and 12 0's, the list that is returned is new_lst = [3, 12]. This function is called in create_confusion_matrix().
#' 
#' @param sample_df Dataframe. Must contain n samples (rows) and a column called "Significance".
#' @return List of 2 elements of int type.
fill_table <- function(sample_df){
  lst_info <- list()
  lst_info[[1]] <- sum(sample_df$Significance)
  lst_info[[2]] <- nrow(sample_df) - sum(sample_df$Significance)
  return(lst_info)
}

#' Build and return a confusion matrix of the significant PERMANOVA results for both host-associated and environmental microbiome data. The basic confusion matrix is created using a function called make_table(). The host-associated and environmental results are used to complete the table, one at a time, using a function called fill_table(). Two parameters are required: data1 and data2, which are dataframes for the environmental microbiome data and host-associated microbiome data, respectively. 
#' 
#' @param data1 Dataframe for the environmental data. Must contained a "Significance" column.
#' @param data2 Dataframe for the host-associated data. Must contained a "Significance" column.
#' @return 2x2 confusion matrix.
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

#' Creates a figure that represents the confusion matrix. It takes in 4 parameters, 1 is the data, and the rest are required for building the plot. The figure that is created plots the significant samples along the y-axis and the maximum distance between samples (km) along the x-axis. All samples that were significant are plotted at y=1, and all insignificant samples are plotted along y = 0. 
#' 
#' @param sample_df A dataframe that must contain 2 columns: "Maximum.Distance..Km." and "Significance".
#' @param title Str. The desired title for the figure.
#' @param shape Int. Required for the points on the plot. 
#' @param colour_id Str. A "#" followed by numbers and letters that represent a specific colour code. The colour is used to fill in the shapes (the points) on the figure. 
#' @return Figure
plot_significance <- function(sample_df, title, shape, colour_id){
  dot_plot <- ggplot(data = sample_df) + geom_point(aes(x = Maximum.Distance..Km., y = Significance), shape = shape, size = 3, colour = "black", fill = colour_id) + 
    labs(title = title, x = "Maximum Distance (km)", y = "Significance") +
    theme(text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
  return(dot_plot)
}
