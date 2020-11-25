# Teodora Tockovska
# Nov 20, 2020

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
    print(paste(dataset, "contains 0 NA values."))
  }
}

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
#' @return The transformed OTU table of \code{phyloseq_OTU_table}.
create_comm_matrix <- function(phyloseq_OTU_table){
  otu_matrix <- as(phyloseq_OTU_table@.Data, "matrix")
  return(t(otu_matrix))
}

# I made this function to get the list that contains the unique UTM zones within the dictionary: I got the dictionary with the longitudes that were associated with multiple UTM zones. Now, I can count how many coordinates are associated with a UTM zone. If the length of longitudes associated with a UTM zone is 1, ignore and move on to a UTM zone that has more than 1 longitude. That way, I can just remove the host species that were only sampled from 1 region.

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
#' @param lat The latitude column name that is present in \code{data}.
#' @param sampleid The sample ID (first column)
#' @param long The longitude column name that is present in \code{data}.
#' @param type Str. The description of the data.
#' @return A dataframe of 5 columns and nrows(data). The columns include: latitude, longitude, year of sampling, month of sampling, and day of sampling.
#' @examples 
#' get_coordinates(birds, latitude_deg, longitude_deg, collection_timestamp)
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

#' Returns a dataframe that contains the UTM coordinates. I found this function on a Stack Overflow page: https://stackoverflow.com/questions/18639967/converting-latitude-and-longitude-points-to-utm. 
#' 
#' @param x The Longitude values
#' @param y The Latitude values
#' @param zone The geographical zone of interest
#' @return A dataframe with UTM coordinates. There are 2 columns, where the first column contains the latitude, named "X", and the second column contains the longitude coordinates, names as "Y".
#' @examples 
#' x <- c( -94.99729,-94.99726,-94.99457,-94.99458,-94.99729)
#' y <- c( 29.17112, 29.17107, 29.17273, 29.17278, 29.17112)
#' LongLatToUTM(x,y,15)
LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
  return(as.data.frame(res))
}

#' This function returns a list of dataframes that contain the XY UTM coordinates. The function loops through the coordinate dataframes and will call the function LongLatToUTM() to convert the latitudes and longitudes to UTM coordinates. I called LongLatToUTM() using zone 30, which is for Spain (from the link https://mangomap.com/robertyoung/maps/69585/what-utm-zone-am-i-in-#). The first column of each UTM dataframe was excluded because it wasn't necessary. 
#' 
#' @param coord_df A list of dataframes that contain the geographic coordinates.
#' @param zone Int or vector of ints. The zone of interest
#' @return A list of dataframes
get_lst_UTM_matrices <- function(coord_df, zone){
  if(class(zone) == "list"){
    if(length(coord_df) != length(zone)){
      print("Error: Your list of coordinates and your vector of zones do not have the same lengths!")
      break
    }
    for (i in 1:length(coord_df)){
      # It is assumed that the zone is a vector of the same length as the list of coordinates.
      df_XY_UTM <- LongLatToUTM(coord_df[[i]]$Longitude, coord_df[[i]]$Latitude, zone[[i]])[,-1]
      coord_df[[i]] <- cbind(coord_df[[i]], df_XY_UTM)
    }
  }
  
  # the zone is only an integer
  else{
    for (i in 1:length(coord_df)){
      df_XY_UTM <- LongLatToUTM(coord_df[[i]]$Longitude, coord_df[[i]]$Latitude, zone)[,-1]
      coord_df[[i]] <- cbind(coord_df[[i]], df_XY_UTM)
    }
  }
  
  return(coord_df)
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

#' Return a list of longitudes that are specific to a UTM zone, zone_int. The function loops through the UTM_zone_lst and matches the zone in the list to zone_int, and will add the longitudes that belong to the matched UTM zone, and that list is returned. This function is called in return_longitudes_per_UTM().
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

#' Returns the vector of longitudes that are associated with a list of UTM zones. 
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


flatten_lst_of_df <- function(coordinate_lst_of_dfs){
  for(sample in 1:length(coordinate_lst_of_dfs)){
    if(class(coordinate_lst_of_dfs$sample) == "list"){
      coordinate_lst_of_dfs <- unlist(coordinate_lst_of_dfs, recursive = F)
    }
  }
  return(coordinate_lst_of_dfs)
}

#' Some times there are hosts that do not have more than 2 listed names so this function is used to deal with the multiple names. A list of names are provided into the function as a parameter and the function iterates through the list. If there is only 1 host name (which means it could be a class or genus name), the function skips that. Once the function sees a name that is > 1 length (meaning it has a species name) it will change the lst_name to the name and it returns the one character vector.
#' 
#' @param lst_names A list of host or environmental sample names. 
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

#' Change the inner lists' names
rename_inner_lists <- function(lst_of_dfs, type){
  
  # Changing the names of the inner lists within lst_of_dfs:
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

#' Remove all the lists within updated_coordinates and remove any empty dataframes: I used this link to help me -- https://stackoverflow.com/questions/50403752/r-remove-list-full-of-na-from-a-list-of-lists
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

extend_names_of_hosts <- function(host_name, num_of_subsets){
  names <- list()
  for(i in 1:num_of_subsets){
    sentence <- paste(host_name, as.character(i), sep = " ")
    names <- append(names, sentence)
  }
  return(names)
}

#' creating a list of the host lists removed that contained either no lists or had lists of dataframes. The lengths are saved to indicate the number of lists of dataframes for each removed host species.
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
  
  nonzero_ids <- unlist(nonzero_ids)
  subset_length_ids <- unlist(subset_length_ids)
  empty_lists <- sort(unlist(empty_lists), decreasing = F)
  
  lists_to_return <- list(nonzero_ids, subset_length_ids, empty_lists, subset_length)
  names(lists_to_return) <- c("nonzero_ids", "subset_length_ids", "empty_lists", "subset_length")
  
  return(lists_to_return)
}

create_extended_host_names <- function(lst_of_lengths, lst_of_host_names){
  if(length(lst_of_host_names) == length(lst_of_lengths)){
    for(i in 1:length(lst_of_host_names)){
      
      # change the names of the host species
      lst_of_host_names[[i]] <- extend_names_of_hosts(lst_of_host_names[[i]], lst_of_lengths[[i]])
    }
  }
  return(lst_of_host_names)
}

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

change_all_multiple_samples_names <- function(output_list, multiple_zones_names){
  
  change_host_names <- multiple_zones_names[output_list$nonzero_ids]
  hosts_info_lengths <- output_list$subset_length[output_list$subset_length_ids]
  
  change_host_names <- create_extended_host_names(hosts_info_lengths, change_host_names)
  hosts_to_remove <- multiple_zones_names[output_list$empty_lists]
  
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
    plot(coordinate_df[[i]]$Longitude, coordinate_df[[i]]$Latitude, xlab = "Longitude", ylab = "Latitude")
    title(paste(type, as.character(i), ":", as.character(host_names[i]), sep = " "))
  }
}

# This function will add a column for the significant values. 0 is not significant and 1 is significant.
check_significance <- function(df_PERMtests){
  df_PERMtests$Significance <- 0
  df_PERMtests$Significance[which(significance_of_tests$Pr..F.<0.05)] <- 1
  return(df_PERMtests)
}

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

change_element_names <- function(studyID_df){
  edit_names_IDs <- grep("*[0-9+]$", names(studyID_df))
  get_names <- studyID_df[edit_names_IDs]
  
  for(i in 1:length(get_names)){
    names(get_names)[i] <- paste(names(get_names[i]), get_names[[i]]$env_material[1], sep = "_")
  }
  
  names(studyID_df)[edit_names_IDs] <- names(get_names)
  
  return(studyID_df)
}

# https://stackoverflow.com/questions/22309285/how-to-use-a-variable-to-specify-column-name-in-ggplot
plot_world_map <- function(sample_df, colour_id, long, lat, title, tiff = FALSE, file_name = NULL){
  world <- ne_countries(scale = "medium", returnclass = "sf")
  theme_set(theme_bw())
  world_plot <- ggplot(data = world) + 
    geom_sf() + 
    
    # why I used the {{}} -- https://stackoverflow.com/questions/2641653/pass-a-data-frame-column-name-to-a-function
    geom_point(data = sample_df, aes(x = {{long}}, y = {{lat}}), size = 2, shape = 21, fill = colour_id) + 
    labs(title = title, x = "Longitude", y = "Latitude") + 
    theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5))
  
  if(tiff && class(file_name) == "character"){
    tiff(paste(file_name, "tiff", sep = "."), units="in", width = 6, height = 3.5, res = 300)
    print(world_plot)
    dev.off()
  }
  
  print(world_plot)
}
