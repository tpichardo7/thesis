#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cleaning data from hilic column
# ESI: negative
# 09/05/2025
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)
library(janitor)

# Functions ----
replacezero <- function(x) "[<-"(x, !x|is.na(x), min(x[x>0], na.rm = T)/2)

# Call in mapping file ----
map <- read_tsv("./data/HILpos_batch_file.txt") # change for HILIC

# Call in raw feature table ----
ft <- read_tsv("./data/ComBat_mzcalibrated_untargeted_mediansummarized_featuretable.txt") # change for HILIC

## filter to retain sample files (remove QC and NIST samples) ----
samples <- map %>% 
    filter(!grepl("NIST", `Sample ID`)) %>% 
    filter(!grepl("QC", `Sample ID`)) %>% 
    separate(`Sample ID`, into = c("ID", "inj"), sep = "_") %>% 
    filter(inj == 1) # for HILIC, change to inj == 1

# Merge sample data with feature table ----
int <- ft %>% 
    select(mz, time, contains("SD")) %>% 
    unite("mz_time", c(mz, time), sep = "_") %>% 
    t() %>% 
    as.data.frame() %>% 
    row_to_names(1) %>% 
    rownames_to_column(var = "File Name")

merged <- samples %>% 
    left_join(., int, by = "File Name")

############### Pre-processing: Filter, replace and normalize #############################

# Detection rate per feature ----
features <- merged %>%
    select(-c(1:4)) %>% 
    map_df(~as.character(.x)) %>%
    map_df(~as.numeric(.x))

# Proportion of samples with 0 intensity per feature ----
sample_zero_prop <- colSums(features==0)/dim(features)[1]
sample_zero_prop_merge <- as.data.frame(sample_zero_prop) %>% 
    rownames_to_column(var= "mz_time")  %>% 
    separate(mz_time, into = c("mz", "time"), sep = "_", remove = FALSE)

# Filter to retain features with detection > 70% of the samples ----
mztime_70 <- ft %>% 
    merge(sample_zero_prop_merge, ., by = c("mz", "time")) %>% 
    filter(sample_zero_prop < 0.3) %>% 
    select(mz_time)

# Select retained features ----
ft_70 <- merged %>% 
    select(ID, `File Name`, all_of(mztime_70$mz_time))

# Replace 0 values with 1/2 minimum intensity and Log2 transform the data ----
ft_70_normalize <- ft_70 %>% 
    select(all_of(mztime_70$mz_time)) %>% 
    map_df(~as.character(.x)) %>%
    map_df(~as.numeric(.x)) %>% 
    map_df(~replacezero(.x)) %>% 
    map_df(~log2(.x)) %>% 
    mutate(ID = ft_70$ID) %>% 
    mutate(`File Name` = ft_70$`File Name`) %>% 
    select(ID, `File Name`, everything())

# Save to disc for merging with meta data and downstream analysis ---- 
ft_70_normalize %>% write_tsv("data/hilic_70detect_log2.csv") # change for HILIC    
mztime_70 %>% write_tsv("data/hilic_70detect_mz_time.csv") # change for HILIC

library(data.table)
hilic_data <- fread("./data/hilic_70detect_log2.csv")

