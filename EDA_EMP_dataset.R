# Teodora Tockovska
# Sep 22, 2020

# This script contains the exploratory data analysis for the Earth Microbiome metadata.

# Reading Data and Loading Script
library("tidyverse")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")

# Reading the EMP dataset file into R.
EMP_metadata <- read.csv(file = "emp_qiime_mapping_qc_filtered.csv", sep = "\t")

# Sample Exploration ####

# Renamed the countries by removing "GAZ:", then I plotted the number of records per country. Nearly half of the samples were taken from the United States of America. I viewed the table of records per country in decreasing order.
EMP_metadata$country <- sub("GAZ:", "", EMP_metadata$country)
pie(sort(table(EMP_metadata$country), decreasing = T))
sort(table(EMP_metadata$country), decreasing = T)

# I changed the Venezuelan bird samples to have the correct longitude sign because the values were incorrectly recorded. To do this, I extracted the row IDs of all the birds and matched the IDs of the birds with the Venezuelan samples. Then, I changed the signs of the longitudes. I used this link to help me use the intersect() function: https://stat.ethz.ch/pipermail/r-help/2013-March/349521.html. Intersect() is from base R and it performs an intersection between 2 vectors.
ven_bird_ids <- intersect(grep("Aves", EMP_metadata$host_class), grep("Venezuela", EMP_metadata$country))
EMP_metadata$longitude_deg[ven_bird_ids] <- -EMP_metadata$longitude_deg[ven_bird_ids]

# I plotted all of the coordinates using plot() just to see the visualization of all coordinates. Recall that latitude is a range from N to S (60deg or 60N, to -60 deg or 60S). The longitude is a range from W to E (-180deg or 180W, to 180deg or 180E). Using these parameters, I plotted the latitude on the Y axis, and the longitude on the X axis. It is hard to tell if there are outliers in the data using this graph.
plot(EMP_metadata$longitude_deg, EMP_metadata$latitude_deg)

# There are countries that have very few env_material samples taken and depending on the type of data and host species, the samples may even decrease further (i.e. the frequency of host species per country could potentially be small). This could pose a problem in my analyses. The first line of code shows the samples taken per country, keeping values that are > 0. I plotted the histogram to view the distribution of samples. The majority of the samples have records < 500. 
samples_per_country <- as.data.frame(table(EMP_metadata$country, EMP_metadata$env_material)) %>% filter(Freq>0)
hist(samples_per_country$Freq)

# I viewed the dataframe of the environmental material frequencies per country that have records > 100. In total, there are 33 rows. I showed the dataframe ordered by the records. USA had the most samples (3578) and the Gulf of Mexico had the fewest samples (103). There are about 26 groups that have <= 1000 records, about 7 groups that have records between 1000 and 2000, and about 2 groups that have records between 3500 and 4000 records. 
samples_per_country <- as.data.frame(table(EMP_metadata$country, EMP_metadata$env_material)) %>% filter(Freq>100)
samples_per_country[order(samples_per_country$Freq, decreasing = T),]

# The histogram plots the frequencies from the dataframe, samples_per_country.
hist(samples_per_country$Freq)

# I checked how many samples were taken from each country per biome. USA has about 5018 samples from the urban biome.
biome_country_EMP <- table(EMP_metadata$env_biome, EMP_metadata$country) %>% data.frame() %>% filter(Freq>0)
biome_country_EMP[order(biome_country_EMP$Freq, decreasing = T),]

# Plotting samples on world map and Splitting Metadata into 2 files ####

# To begin plotting the data, I have found 2 tutorials online (https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html and https://www.r-spatial.org/r/2018/10/25/ggplot2-sf-2.html), which use ggplot2 and sf packages. 
world <- ne_countries(scale = "medium", returnclass = "sf")

# I plotted the EMP dataset on the world map by using the coordinates. I set the theme of the map to black and white. The blue dots represent from where the hosts were sampled. There seem to be samples taken from the Atlantic ocean and a sample in the Indian ocean but they look odd and I am not sure if they are outliers.
theme_set(theme_bw())
ggplot(data = world) + geom_sf() + geom_point(data = EMP_metadata, aes(x = longitude_deg, y = latitude_deg), size = 1, shape = 23, fill = "blue") + labs(title = "World Map of EMP Metadata", x = "Longitude", y = "Latitude")

# So, the next step is to get the number or records per sample selected by host species and country. Also, keep in mind that the number of host species is only ideal for samples taken from hosts and not environmental data. 

# Checking number of host classes in the dataset. There are 12,451 empty host classes, which could mean that those samples are either environmental or just not recorded. There are 1991 samples labelled as "c__". I don't know what that could mean. There are 235 Actinopteri host classes (ray-finned fish), 39 Amphibia, 48 Anthozoa (marine invertebrates such as sea anemones and corals), 716 Aves (birds), 6 Bivalia (molluscs), 23 Florideophyceae (red algae), 839 Insecta, 613 Liliopsida (vascular), and 6867 Mammalia.
table(EMP_metadata$host_class)

# Now, analyzing the empty host classes labelled as "c__". In this group, the host_class was not recorded correctly. Some of them do have the species name recorded. The species might have been recorded because they are samples taken from the environment and they are associated with a specific host.
c__hosts <- EMP_metadata[grep("c__$", EMP_metadata$host_class),]
unique(c__hosts$country)
unique(c__hosts$host_scientific_name)
unique(c__hosts$host_common_name_provided)
unique(c__hosts$host_kingdom)
unique(c__hosts$host_common_name)
unique(c__hosts$host_species)

# Interestingly, there are some samples here that are environmental: It seems that all samples here are hosts from the environment. These are the environments in which the samples are involved or partake.
unique(c__hosts$sample_scientific_name)

# Apparently, many countries did not record the host information but this is because these are environmental samples. Now, checking the empty host class data. Ultimately, I will need to combine the two environmental host samples together into one dataframe so that I can run the spatial pattern analyses.
empty_class <- EMP_metadata[grep("^$", EMP_metadata$host_class),]
unique(empty_class$country)
head(unique(empty_class$Description))
unique(empty_class$sample_scientific_name)

# Well, I can keep the host classes that have > 30 records, and split the environmental microbiomes from the host-associated ones. First, I split the data according to the sample types.
as.data.frame(table(EMP_metadata$host_class)) %>% filter(Freq>30)

# Combine the available host classes together into one dataframe. First, get the ids and then extract those host species from the EMP metadata.
host_species_ids <- grep("Actinopteri|Aves|Insecta|Liliopsida|Mammalia|Amphibia|Anthozoa", EMP_metadata$host_class)
host_associated_microbiome_dataset <- EMP_metadata[host_species_ids,]
environmental_microbiome_dataset <- EMP_metadata[-host_species_ids,]

# Plotting the two datasets on the world map. First one will plot the host species. The second world map plots the environmental metadata.
tiff("Hosts_world_map.tiff", units="in", width = 6, height = 3.5, res = 300)
ggplot(data = world) + 
  geom_sf() + 
  geom_point(data = host_associated_microbiome_dataset, aes(x = longitude_deg, y = latitude_deg), size = 2, shape = 21, fill = "#D55E00") + 
  labs(title = "World Map of Host Species in EMP Metadata", x = "Longitude", y = "Latitude") +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
dev.off()

tiff("ENV_world_map.tiff", units="in", width = 6, height = 3.5, res = 300)
ggplot(data = world) + 
  geom_sf() + 
  geom_point(data = environmental_microbiome_dataset, aes(x = longitude_deg, y = latitude_deg), size = 2, shape = 21, fill = "#882255") + 
  labs(title = "World Map of Environmental Microbiome Data", x = "Longitude", y = "Latitude") +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))
dev.off()

# Removing environmental variables that I don't need:
rm(EMP_metadata, empty_class, world, host_species_ids, samples_per_country, biome_country_EMP, c__hosts, ven_bird_ids)

# Saving the new dataframes into separated CSV files so that I can open them in the main analysis script.
# write.table(environmental_microbiome_dataset, "free_living_EMP.csv", quote = F, sep = "\t")
# write.table(host_associated_microbiome_dataset, "host_associated_EMP.csv", quote = F, sep = "\t")
