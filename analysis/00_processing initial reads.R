
# processing initial raw data
library(tidyverse)

# Set the directory where the files are located
folder_path <- "rawdata/rawcounts"

# Get a list of all file names in the directory
file_list <- list.files(path = folder_path, pattern = "\\.tabular$", full.names = TRUE)

# Initialize an empty list to store data frames
data_list <- list()

# Loop through the files and read each file into a data frame
for (file in file_list) {
  # Read the file as a single-column data frame, skipping the first row, using tab as the separator
  temp_df <- read.table(file, row.names = 1, header = FALSE, sep = "\t", skip = 1)
  
  # Use the file name (without path and extension) as the column name
  col_name <- tools::file_path_sans_ext(basename(file))
  
  # Keep only the first 5 characters of the column name
  col_name <- substr(col_name, 1, 5)
  
  # Rename the column of the temporary data frame
  colnames(temp_df) <- col_name
  
  # Add the data frame to the list
  data_list[[col_name]] <- temp_df
}

# Combine all data frames efficiently using bind_cols
final_df <- bind_cols(data_list)

write.csv(final_df, "rawdata/KT_Counts.csv")
