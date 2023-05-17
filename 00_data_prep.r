############################################################################# #
# Data prep for paper'More than 75% of undescribed plant species are likely
# threatened'
# MJM Brown, SP Bachman, E Nic Lughadha
##############################################################################

# load packages ----
library(tidyverse)
library(rWCVP)
library(sjPlot)
library(sf)
library(lme4)

# get WCVP list of accepted species and distributions ----

#load data
dist <- rWCVPdata::wcvp_distributions
names <- rWCVPdata::wcvp_names
tdwg <- rWCVPdata::wgsrpd3 %>% st_drop_geometry()


# filter to accepted species
df <- names %>%
  filter(taxon_status=="Accepted",
         taxon_rank == "Species",
         is.na(genus_hybrid),
         is.na(species_hybrid),
         is.na(infraspecies)) #should be redundant but better to be safe

df_nobasio <- df %>% select(plant_name_id, taxon_name, family, genus,
                            lifeform_description, first_published,
                            basionym_plant_name_id) %>%
  filter(is.na(basionym_plant_name_id))


df_basio <- df %>% select(plant_name_id, taxon_name, family, genus,
                          lifeform_description, first_published,
                          basionym_plant_name_id) %>%
  filter(!is.na(basionym_plant_name_id)) %>%
  left_join(rWCVPdata::wcvp_names %>%
              select(plant_name_id, basio_year=first_published),
            by=c("basionym_plant_name_id"="plant_name_id"))

df_basio$year <- df_basio$basio_year
df_nobasio$year <- df_nobasio$first_published

df <- bind_rows(df_nobasio, df_basio)


# join ids to distributions
accepted_distributions <-
  dist %>%
  mutate(area_code_l3 = toupper(area_code_l3)) %>%
  filter(extinct + introduced + location_doubtful == 0) %>%
  select(area_code_l3, plant_name_id) %>%
  inner_join(
    df,
    by="plant_name_id"
  )

# clean up code names
accepted_distributions <-
  accepted_distributions %>%
  filter(area_code_l3 %in% tdwg$LEVEL3_COD)


distcentroidlatitude <- function(x){
 suppressWarnings(
   abs(
     mean(
       sf::st_bbox(
         st_combine(
           suppressMessages(sf::st_union(
           rWCVPdata::wgsrpd3[which(rWCVPdata::wgsrpd3$LEVEL3_COD %in% x),]
           )
         ))
         )[4],
     sf::st_bbox(
       st_combine(
         suppressMessages(sf::st_union(
           rWCVPdata::wgsrpd3[which(rWCVPdata::wgsrpd3$LEVEL3_COD %in% x),]
         )
       ))
     )[2]
     )
     )
   )
}

# commented out because it is slooow
# climdist <- accepted_distributions %>%
#   group_by(plant_name_id) %>%
#   summarise(c.lat = distcentroidlatitude(area_code_l3)) %>%
#   mutate(clim_bin = case_when(c.lat<23.3~"tropical",
#                               c.lat>23~"temperate",
#                               TRUE~NA_character_)
# write_csv(climdist, "wcvpv10_climatemapping.csv")

climdist <- read_csv("wcvpv10_climatemapping.csv")



# remove infraspecies ids and deduplicate by species and region
accepted_distributions <-
  accepted_distributions %>%
  select(area_code_l3, taxon_name, plant_name_id) %>%
  group_by(plant_name_id, area_code_l3) %>%
  filter(row_number() == 1) %>%
  ungroup()

endemics <- accepted_distributions %>% group_by(plant_name_id) %>%
  summarise(n=n()) %>%
  filter(n==1)


lifeform_mapping <- read_csv("lifeform_mapping_2023.csv")


df <-
  df %>%
  left_join(lifeform_mapping, by="lifeform_description") %>%
  rename(standard_lifeform=humphreys_lifeform) %>%
  left_join(climdist)

redlist <- read.csv("C:/Users/mbr10kg/OneDrive - The Royal Botanic Gardens, Kew/Desktop/Projects/Name-matching/RedList2023-1/RL-2023-1_WCVP-v10.csv") %>%
  filter(criteriaVersion==3.1)


df <- df %>%
  left_join(redlist, by=c("plant_name_id"="wcvpv10_id"))

threatened_cats <- c("Vulnerable","Endangered","Critically Endangered","Extinct in the Wild", "Extinct")
df$threatened <- as.factor(df$redlistCategory %in% threatened_cats)
df$threatened[which(is.na(df$redlistCategory))] <- NA
df$threatened[which(df$redlistCategory=="Data Deficient")] <- NA
df$threatened_dd <- as.factor(df$redlistCategory %in% c(threatened_cats,"Data Deficient"))
df$threatened_dd[which(is.na(df$redlistCategory))] <- NA

df$threatenedEN <- as.factor(df$redlistCategory %in% threatened_cats[2:5])
df$threatenedEN[which(is.na(df$redlistCategory))] <- NA
df$threatenedEN[which(df$redlistCategory=="Data Deficient")] <- NA
df$threatenedCR <- as.factor(df$redlistCategory %in% threatened_cats[3:5])
df$threatenedCR[which(is.na(df$redlistCategory))] <- NA
df$threatenedCR[which(df$redlistCategory=="Data Deficient")] <- NA

df$endemic <- df$plant_name_id %in% endemics$plant_name_id

df$year_orig <- df$year


df$year[str_detect(df$year,"19553")] <- "1953" # Erycibe subglabra
df$year[str_detect(df$year,"18881")] <- "1881" # Inula racemosa
df$year[str_detect(df$year,"18412")] <- "1841" # Lycopodium ceylanicum
df$year <- gsub("\\(", "",df$year )
df$year <- gsub("\\)", "",df$year )
years <- str_extract_all(df$year, "\\d\\d\\d\\d")
df$year <- unlist(lapply(years, function(x) x[1]))
df$year <- as.numeric(df$year)

df <- df %>% filter(!is.na(year)) %>% filter(year<2022)

write_csv(df, "df_inclNE.csv")


df <- df %>% select(taxon_name, genus, family, threatened, threatenedEN, threatenedCR, redlistCategory, year, standard_lifeform, clim_bin, endemic) %>%
  unique() %>%
  filter(!taxon_name %in% .$taxon_name[duplicated(.$taxon_name)] | threatened==FALSE)

write_csv(df, "df_fit_inclNA.csv")
write_csv(df %>% na.omit(), "df_fit.csv")
