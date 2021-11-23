#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-r", "--raw_data"), type="character", default=NULL,
              help="dataset folder containing .jpg, .raw and .hdr files, for one timepoint (only jpg needed)", metavar="character"),
  make_option(c("-R", "--Regression"), type="character", default=12,
              help="dataset folder ", metavar="character"),
  make_option(c("-i", "--y_intercept"), type="numeric", default=20,
              help="0 or 1, need to know so we look in right folder for CLS results", metavar="numeric"),
  make_option(c("-d", "--date"), type="character", default=0,
              help="formatted as 2020-03-15 - this is the date regression was run", metavar="character")
);

# Because datestamp in bash ends up with quotes...

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

arg_out_path <- "prelabel_args.rds"

### IF DEBUGGING IN RSTUDIO, UNCOMMENT THIS LINE INSTEAD OF USING OptParser
#opt <- readRDS("/scratch2/NSF_GWAS/GMOdetectoR/prelabel_args.rds")

saveRDS(opt, file = arg_out_path)

library(data.table)
library(stringr)

print(opt$date)
this_date <- gsub('”', '', opt$date)
CLS_dir <- paste0(opt$Regression, str_split_fixed(opt$raw_data, "/", 5)[5])
print(paste0("Looking for CLS data in: ", CLS_dir))

print(this_date)

# Load all data, organize in data tables to be merged

all_segment_files <- list.files(opt$raw_data,
                                pattern = 'segment_uncropped_processed',
                                full.names = TRUE)

# Before splitting any names, get rid of the _cyan if it is there
#print("Head of segment files")
#print(head(all_segment_files))

segment_data_table <- data.table(cbind(all_segment_files, str_split_fixed(basename(gsub("_cyan", "", all_segment_files)), "_", 10)))[,1:9]

all_rgb_files <- list.files(opt$raw_data,
                            pattern = 'rgb_processed',
                            full.names = TRUE)

#print("Head of RGB files")
#print(head(all_rgb_files))

rgb_data_table <- data.table(cbind(all_rgb_files, str_split_fixed(basename(gsub("_cyan", "", all_rgb_files)), "_", 10)))[,1:9]

print(paste0("Looking for CLS files in directory ", CLS_dir))
all_CLS_files <- list.files(CLS_dir,
                            full.names = TRUE,
                            pattern = 'hdf')

all_CLS_files <- all_CLS_files[!grepl("Broadband", all_CLS_files)]

print(paste0("How many CLS files? ", length(all_CLS_files)))

CLS_data_table <- data.table(cbind(all_CLS_files, str_split_fixed(basename(gsub("_cyan", "", all_CLS_files)), "_", 10)))[,1:9]

# Remove rows with chroma standard
#CLS_data_table[!grepl("hroma", CLS_data_table$all_CLS_files),]


# Format and merge tables

colnames(rgb_data_table) <- colnames(segment_data_table) <- colnames(CLS_data_table) <- c("path",
                                                                                         "Tray",
                                                                                         "IntegrationTime",
                                                                                         "Focus",
                                                                                         "LaserIntensity",
                                                                                         "Timestamp",
                                                                                         "ImageRankofAllImagedonTray",
                                                                                         "Row",
                                                                                         "Col")
CLS_data_table <- CLS_data_table[!grepl("hroma", CLS_data_table$path),]

combined_data_table.1 <- merge(rgb_data_table,
                            segment_data_table,
                            by=c("Tray",
                                 "IntegrationTime",
                                 "Focus",
                                 "LaserIntensity",
                                 "Timestamp",
                                 "Row",
                                 "Col"))

combined_data_table.2 <- merge(combined_data_table.1,
                            CLS_data_table,
                            by=c("Tray",
                                 "IntegrationTime",
                                 "Focus",
                                 "LaserIntensity",
                                 "Timestamp",
                                 "Row",
                                 "Col"),
                             all.x = TRUE,
                             all.y = TRUE)

# Now sort and remove first instances of duplicates (bad images)

combined_data_table.3 <- combined_data_table.2[order(-Timestamp)]

combined_data_table.3$Tray_Row_Col <- paste0(combined_data_table.3$Tray,
                                             "_",
                                            combined_data_table.3$Row,
                                             "_",
                                            combined_data_table.3$Col)

combined_data_table.3$SeenBefore <- numeric(0)

Tray_Row_Col_seen <- c('NA')

for(i in 1:nrow(combined_data_table.3)){
    #print(paste0('Have so far seen ', length(Tray_Row_Col_seen), ' total'))
    this_plate <- combined_data_table.3$Tray_Row_Col[i]
    #print(paste0('Checking if ', this_plate, ' is in list of those seen already'))
    if(sum(grepl(this_plate, Tray_Row_Col_seen))>=1){
        print(paste0('We will only include in analysis the FINAL of multiple images for ', combined_data_table.3$Tray_Row_Col[i]))
        combined_data_table.3$SeenBefore[i] <- 1
    }
    else{
        #print(paste0('Not already seen. Adding to list.'))
        combined_data_table.3$SeenBefore[i] <- 0
        Tray_Row_Col_seen <- append(Tray_Row_Col_seen, combined_data_table.3$Tray_Row_Col[i])
    }
}

combined_data_table.4 <- combined_data_table.3[SeenBefore==0]

# Prepare output

out <- as.data.table(cbind(combined_data_table.4$path.x, combined_data_table.4$path.y, combined_data_table.4$path)) #combined_data_table.4[,9:11]
out$mean_callus_signal <- out$mean_shoot_signal <- out$callus_signal_total <- out$shoot_signal_total <- NA
out$n_pixels_callus_transgenic <- out$n_pixels_callus_escape <- NA
out$n_pixels_shoot_transgenic <- out$n_pixels_shoot_escape <- out$threshold <- NA
colnames(out)[1:3] <- c("rgb", "segment", "CLS_data")
setcolorder(out, c("segment", "CLS_data", "rgb", "mean_callus_signal", "mean_shoot_signal",
                  "callus_signal_total", "shoot_signal_total", "n_pixels_callus_transgenic", "n_pixels_callus_escape",
                  "n_pixels_shoot_transgenic", "n_pixels_shoot_escape", "threshold"))

sample_out_path <- paste0(opt$raw_data, "/samples_pre_labeling.csv")

if( nrow(na.omit(out[,1:3])) < nrow(out[,1:3])) {
  warning("There is an inconsistency. There are some files for which we do not have both hyperspectral and RGB data. These have been excluded.")
  cat("\n")
  print("Rows missing RGB data:")
  print(out[ is.na(out$rgb) , 1:3])
  cat("\n")
  print("Rows missing segment data:")
  print(out[ is.na(out$segment) , 1:3])
  cat("\n")
  print("Rows missing regressopm data:")
  print(out[ is.na(out$CLS_data) , 1:3])
  cat("\n")
}

out <- out[ !is.na(out$rgb) , ]
out <- out[ !is.na(out$segment) , ]
out <- out[ !is.na(out$CLS_data) , ]

print(paste0('Writing ', nrow(out), ' rows to ', sample_out_path))

fwrite(out,
       sample_out_path,
       quote = FALSE,
       row.names = FALSE)
