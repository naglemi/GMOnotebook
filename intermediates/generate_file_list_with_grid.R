# Prepare file list with same hyperchannel_to_csv.jpg
# to match with every image instead of testing
# alignment for each individually

setwd("/media/michael/Elements_121/GWAS_Transformation/GTK/wk7")

library(data.table)
file_list <- fread("file_list.csv")
file_list$hyper_img <- paste0(dirname(file_list$rgb_images),
                                      "/hypercube_to_csv.jpg")

#file_list

file_list <- file_list[, c(2, 1)]
#head(file_list)

fwrite(file_list, "rgb_and_hyper_channel_lists.csv",
       quote = FALSE,
       row.names = FALSE)
