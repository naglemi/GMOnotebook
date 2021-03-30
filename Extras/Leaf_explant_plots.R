library(data.table)
library(scales)
#library(GMOdetectoR)
library(readxl)
library(ggplot2)
library(randomcoloR)
library(optparse)
library(tools)
library(stringr)


# Read arguments from command line  ---------------------------------------


option_list = list(
  make_option(c("-d", "--datapath1"),
              type="character",
              default=NULL,
              help="data output from GMOlabeleR",
              metavar="character"),
  make_option(c("-r", "--randomization_datasheet_path"),
              type="character",
              default=NULL,
              help="path to randomization datasheet",
              metavar="character"),
  make_option(c("-p", "--pixel_threshold"),
              type="numeric",
              default=5,
              #help="Number of pixels passing intensity threshold for binary classifications of transgenic or not",
              metavar="numeric"),
  make_option(c("-v", "--variable"),
              type="character",
              default="categorical",
              help="categorical or continuous",
              metavar="numeric"),
  make_option(c("-m", "--missing"),
              type="numeric",
              default="FALSE",
              help="TRUE (1) if some explants are missing. Must also specify list with -M",
              metavar="numeric"),
  make_option(c("-M", "--MissingList"),
              type="character",
              default="categorical",
              help="A file output from automated detection of missing explants",
              metavar="numeric"),
  make_option(c("-g", "--grid_type"),
              type="numeric",
              default="12",
              help="Grid type (currently supporting 12 and 20)",
              metavar="numeric"),
  make_option(c("-s", "--sort"),
              type="numeric",
              default="0",
              help="Whether to sort genotypes according to difference in effects of treatment on trangenic callus - numeric 0 (False) or 1 (True)",
              metavar="numeric"),
  make_option(c("-H", "--height"),
              type="numeric",
              default="7",
              help="Plot height (inches)",
              metavar="numeric"),
  make_option(c("-w", "--width"),
              type="numeric",
              default=-9,
              help="Plot width (inches)",
              metavar="numeric"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

setwd("/home/labgroup/code/GMOlabeler/")

### IF DEBUGGING IN RSTUDIO, UNCOMMENT THIS LINE INSTEAD OF USING OptParser
# opt <- readRDS("/home/labgroup/code/GMOlabeler/plots/Elements_10/Transformation/GRF_Poplar/GRF1_Phase_2/wk3/gridplot_args.rds")

# Import and preprocess data ----------------------------------------------

rundir <- getwd()

wd <- paste0(rundir, "/plots/", opt$datapath1)
if(!dir.exists(wd)) dir.create(wd, recursive = TRUE)
setwd(wd)

arg_out_path <- paste0(wd, "gridplot_args.rds")
print(paste0("Saving list of input arguments to : ", arg_out_path))
saveRDS(opt, file = arg_out_path)

#datapath <- paste0("/scratch2/NSF_GWAS/GMOlabeler/output/", opt$datapath1, "stats.csv")
datapath <- paste0(rundir, "/output/", opt$datapath1, "stats_with_sums_over_tissues.csv")
#randomization_datasheet_path <- "/scratch2/NSF_GWAS/macroPhor_Array/T16_DEV_genes/EA_randomized.xlsx"
print(paste0("Reading in output from GMOlabeler at path: ", datapath))

pixel_demographics = data.frame(cbind(c('Shoot', 'Callus', 'Stem', 'All_tissue', 'All_regenerated_tissue', 'Background'),
                                      c('00CC11', '0006CC', 'CC0000', 'All_tissue', 'All_regenerated_tissue', '000000'),
                                      c('green', 'blue', 'red', 'black', NA, NA)))
colnames(pixel_demographics) <- c('Tissue', 'hex_code', 'color')

output <- fread(datapath)
cat("\n")
print(paste0("Rows in output from GMOlabeler: ", nrow(output)))

cat("\n")
print(paste0("Max n_pixels_passing_threshold in output from GMOlabeler: ", max(na.omit(output$n_pixels_passing_threshold))))

cat("\n")
print(paste0("Max total_signal in output from GMOlabeler: ", max(na.omit(output$total_signal))))

output$filename <- gsub("_segment_uncropped_processed", "", output$filename)
output$filename <- gsub("__", "_", output$filename) # Saw this naming error in T19 (e.g. TAP1_I5.0_F1.9__L100_105739_1_0_1)

# Merge in missing explant data -------------------------------------------


if(opt$missing==TRUE){
  library(tidyr)
  missing_explant_datapath <- opt$MissingList
  missing_explant_data <- fread(missing_explant_datapath,
                                header=TRUE,
				colClasses=c("character")) # this should be specified explicitly to stop columns from being read as numeric when missing data comes from score data
  #print("preparing to pivot missing explant data table")
  missing_explant_data_tidy <- pivot_longer(data= missing_explant_data,
                                            cols = colnames(missing_explant_data)[-1],
                                            names_to = "grid_item",
                                            values_to = "present",
                                            values_drop_na = TRUE
  )
  
  colnames(missing_explant_data_tidy)[1] <- "filename"
  missing_explant_data_tidy$filename <- gsub("__", "_", missing_explant_data_tidy$filename) # Saw this naming error in T19 (e.g. TAP1_I5.0_F1.9__L100_105739_1_0_1)

  output$filename <- gsub(".jpg", "", output$filename) # Get rid of file extension to be consistent with outputs from automated missing explant script

  
  cat("\n")
  print("Look at the top of output from GMOlabeler")
  print(head(output$filename))
  
  cat("\n")
  print("Look at the top of missing explant data")
  print(head(missing_explant_data_tidy))
  
  output$filename <- basename(output$filename)
  output$grid_item <- as.character(output$grid_item)
  
  cat("\n")
  print("Look at the head of filenames for missing explant data")
  print(head(missing_explant_data_tidy$filename))
  
  cat("\n")
  print("Look at the colnames of output from GMOlabeler")
  print(colnames(output))
  
  cat("\n")
  print("Look at the colnames of missing explant data")
  print(colnames(missing_explant_data_tidy))
  
  # We merge like this since we only import missing explant data for explant
  # that are missing. The rest is to be imputed as not being missing.
  
  cat("\n")
  print(paste0("Before merging the GMOlabeler and missing output data, they have rows of ",
        nrow(output), " and ", nrow(missing_explant_data_tidy), " respectively."))
  output <- merge(x = output,
                  y = missing_explant_data_tidy,
                  by = c("filename", "grid_item"),
                  all.x = TRUE,
                  all.y = FALSE)
  # Remove explants labeled as missing or contaminated
  cat("\n")
  print(paste0("After merging but before removing missing/contminated explants, merged output has ",
               nrow(output), " rows."))
  output$present[is.na(output$present)] <- "Y"
  
  # cat("\n")
  # print(paste0("All of these ",
  #              nrow(output[which(output$present=="N"),]),
  #              " rows are for missing explant"))
  # print(output[which(output$present=="N"),])
  # 
  # cat("\n")
  # print(paste0("All of these ",
  #              nrow(output[which(output$present=="C"),]),
  #              " rows are for contaminated explant"))
  # output[which(output$present=="C"),]
  
  missing_explant_data[missing_explant_data=='NC'] <- 1
  missing_explant_data[missing_explant_data=='Y'] <- 1
  missing_explant_data[missing_explant_data=='N'] <- 0
  missing_explant_data[missing_explant_data=='C'] <- 0
  missing_explant_data[missing_explant_data=='P'] <- 0
  missing_explant_data[missing_explant_data=='M'] <- 0
  cat("\n")
  print(paste0("After removing the missing/contaminated explants, merged output has ",
               nrow(output)))
}

cat("\n")
print("Afer processing GMOlabeler output to merge with missing explant data, Head of data (w/ first 5 col): ")
output[1:5,1:5]

for(i in 1:nrow(output)){
  #' Given macroPhor Array output filename, parse out tray and plate IDs
  #'
  #' @param filename for macroPhor Array output, with file naming as used in Strauss Lab
  #'
  #' @return A character string with tray ID and plate ID delimited by "_"
  #' @export
  #'
  #' @examples
  parse_trayplateID <- function(name_being_parsed){
    pass_to_dodge_error <- name_being_parsed
    imgpath_stripped <- file_path_sans_ext(basename(pass_to_dodge_error))
    trayID <- str_split_fixed(imgpath_stripped, "_", 2)[1]
    
    assign_ID_index_from_row_column_on_tray <- function(data_to_parse = filename, components_list, mode="table", verbose=FALSE){
      dictionary <- cbind(c(0,0,0,0,0,0,0,
                            1,1,1,1,1,1,1,
                            2,2,2,2,2,2,2),
                          c(0,1,2,3,4,5,6,
                            0,1,2,3,4,5,6,
                            0,1,2,3,4,5,6),
                          c(1:21))
      dictionary <- as.data.table(dictionary)
      colnames(dictionary) <- c("row", "column", "ID")
      # Get the ID of position in tray in according to row and column
      dictionary$row_column <- paste0(dictionary$row, "_", dictionary$column)
      dictionary[,1:2] <- NULL
      if(mode=="table"){
        # Set colnames for spectral components if multiple are same
        #colnames(data_to_parse)[1:length(components_list)] <- components_list
        data_merged <- merge(data_to_parse, dictionary, by="row_column", all.x = TRUE, all.y = TRUE)
        return(data_merged)
      }
      if(mode=="filename"){
        
        # Patch added in v0.19 for compatibility regardless of whether "_cyan" is at end of filename
        if(grepl("cyan", data_to_parse)==1){
          ndelimiters=9
        }else{
          ndelimiters=8
        }
        
        row <- str_split_fixed(basename(file_path_sans_ext(data_to_parse)), "_", 9)[ndelimiters-1]
        # Changed in v0.19 along with patch above
        col <- str_split_fixed(basename(file_path_sans_ext(data_to_parse)), "_", ndelimiters)[ndelimiters]
        row_col <- paste0(row, "_", col)
        ID <- dictionary[which(dictionary$row_column == row_col),]$ID
        # Debugging lines added in v0.19
        if(verbose==TRUE){
          print(paste0("This row is ", row))
          print(paste0("This col is ", col))
          print(paste0("This row_col is ", row_col))
          print(paste0("This filename (stripped) is ", basename(file_path_sans_ext(data_to_parse))))
          print(paste0("This ID about to be returned from assign_ID_index_from_roW_column_on_tray is ", ID))
        }
        
        return(ID)
      }
    }
    
    plateID <- assign_ID_index_from_row_column_on_tray(data_to_parse = imgpath_stripped, mode="filename")
    trayplateID <- paste0(trayID, "_", plateID)
    return(trayplateID)
  }
  output$ID[i] <- parse_trayplateID(name_being_parsed = output$filename[i])
}


# Read randomization datasheet, clean if needed ---------------------------


randomization_datasheet <- read_excel(opt$randomization_datasheet_path)

# If there are leaf explants, note that we can only analyze stem explants in the current version
if (sum(grepl('Tissue type', colnames(randomization_datasheet))) >= 1){
  message("Warning! This dataset contains multiple explant types. We will only analyze LEAF explants.")
  print("Randomization datasheet contains levels of tissue type: ")
  print(levels(factor(randomization_datasheet$`Tissue type`)))
  print(paste0("Nrow before subsetting to leaf only: ", nrow(randomization_datasheet)))
  randomization_datasheet <- randomization_datasheet[which(randomization_datasheet$`Tissue type` == "L"), ]
  print(paste0("Nrow after subsetting to leaf only: ", nrow(randomization_datasheet)))
}

colnames(randomization_datasheet) <- gsub("total_explants", "total_explants_from_master_data", colnames(randomization_datasheet))

## Let's deal with inconsistently named columns here
colnames(randomization_datasheet) <- gsub("Genotype ID", "Genotype_ID", colnames(randomization_datasheet))
# If we don't already have a column named Treatment name, conver the Treatment column to this.
if(sum(grepl("Treatment name", colnames(randomization_datasheet))) == 0){
  colnames(randomization_datasheet) <- gsub("Treatment", "Treatment name", colnames(randomization_datasheet))
}


print(paste0("Upon loading randomization_datasheet, it has how many rows? ",
             nrow(randomization_datasheet)))

colnames(randomization_datasheet)[1:2] <- c("Image #", "Tray ID")
randomization_datasheet$ID <- paste0(randomization_datasheet$`Tray ID`,
                             "_",
                             randomization_datasheet$`Image #`)

output$transgenic <- rep(0, nrow(output))
output$transgenic[which(output$n_pixels_passing_threshold > opt$pixel_threshold)] <- 1
output$segment_present <- 0
output$segment_present[!is.na(output$total_signal)] <- 1
for (i in 1:nrow(pixel_demographics)){
    output$segment_hex <- gsub(pixel_demographics$hex_code[i],
                                      pixel_demographics$Tissue[i],
                                      output$segment_hex)
}

## Calculate total number of missing explants per plate

# # Redundant; done above
# missing_explant_data[missing_explant_data=='Y'] <- 1
# missing_explant_data[missing_explant_data=='N'] <- 0
# missing_explant_data[missing_explant_data=='C'] <- 0
# missing_explant_data[missing_explant_data=='P'] <- 0

cat("\n")
print("Calculating total numbers of missing explants per plate")
missing_explant_data_with_totals <- missing_explant_data
missing_explant_data_matrix <- as.matrix(missing_explant_data[,2:ncol(missing_explant_data)])
class(missing_explant_data_matrix) <- "numeric"
missing_explant_data_matrix[is.na(missing_explant_data_matrix)] <- 0
missing_explant_data_with_totals$total_explants <- rowSums(missing_explant_data_matrix)


print("Let's figure out why contaminated/missing plates are being said to have 0 instead of NA")
cat("\n")
print("Dimensions of missing explant data with totals")
print(dim(missing_explant_data_with_totals))
print("Minimum number of explants on any plate after removing missing and contaminated explants is...")
print(min(na.omit(missing_explant_data_with_totals$total_explants)))

print(head(missing_explant_data_with_totals))

# Format to get ready to merge with randomization datasheet (which contains totals values per plate)
colnames(missing_explant_data_with_totals)[1] <- "filename"
missing_explant_data_with_totals <- cbind(missing_explant_data_with_totals, rep(NA, nrow(missing_explant_data_with_totals)))
colnames(missing_explant_data_with_totals)[length(missing_explant_data_with_totals)] <- "ID"
missing_explant_data_with_totals$filename <- gsub("_rgb", "", missing_explant_data_with_totals$filename)
missing_explant_data_with_totals$filename <- gsub(".jpg",
                                           "",
                                           missing_explant_data_with_totals$filename)

cat("\n")
print("About to add IDs to missing explant data based on filenames. Here is head of filenames.")
print(head(missing_explant_data_with_totals$filename))
cat("\n")
for(i in 1:nrow(missing_explant_data_with_totals)){
  missing_explant_data_with_totals$ID[i] <- parse_trayplateID(name_being_parsed = missing_explant_data_with_totals$filename[i])
}

print("We have added IDs to missing explant data based on filenames. Here is head of IDs.")
print(head(missing_explant_data_with_totals$ID))
cat("\n")

print("N rows of randomization datasheet before merging with missing explant data: ")
print(nrow(randomization_datasheet))
print("Head of IDs in missing explant data")
print(head(missing_explant_data_with_totals$ID))
print("Head of IDs in randomization datasheet")
print(head(randomization_datasheet$ID))

cat("\n")
print("Now about to merge randomization_datasheet and missing explant data.")
print(paste0("Rows in randomization_datasheet: ", nrow(randomization_datasheet)))
print(paste0("Rows in missing explant data: ", nrow(missing_explant_data_with_totals)))

cat("\n")
print("Now about to merge randomization datasheet and missing explant data, by ID. Look at head of IDs for both.")
print(head(randomization_datasheet$ID))
print(head(missing_explant_data_with_totals$ID))

randomization_datasheet <- merge(x = randomization_datasheet,
                                 y = missing_explant_data_with_totals,
                                 by = "ID",
                                 all.x = TRUE,
                                 all.y = FALSE)

print(head(randomization_datasheet))
print(randomization_datasheet$total_explants)

print("Still figuring, now in merged rd, out why contaminated/missing plates are being said to have 0 instead of NA")
cat("\n")
print("Dimensions of missing explant data with totals")
print(dim(missing_explant_data_with_totals))
print("Levels of total explants on any plate after removing missing and contaminated explants is...")
print(levels(factor(randomization_datasheet$total_explants)))
print("Minimum number of explants on any plate after removing missing and contaminated explants is...")
print(min(na.omit(randomization_datasheet$total_explants)))
cat("\n")

print(paste0("Rows in merged output: ", nrow(randomization_datasheet)))
print(randomization_datasheet$total_explants)

cat("\n")
randomization_datasheet$total_explants[is.na(randomization_datasheet$total_explants)] <- opt$grid_type
print("Assuming for all explants for which we have no missing explant data that there are:")
opt$grid_type
cat("\n")

print("N rows of randomization datasheet after merging with missing explant data: ")
print(nrow(randomization_datasheet))

print("Column names of randomization datasheet afer merging with missing explant data: ")
print(colnames(randomization_datasheet))

cat("\n")
print(paste0("In output, Maximum observations for a segment in a grid item in a filename (should be 1): ", 
             max(table(output$filename, output$segment_hex, output$grid_item))))
print(paste0("In randomization datasheet, Maximum observations for an ID (should be 1): ", 
             max(table(randomization_datasheet$ID))))


# Debugging before calculations -------------------------------------------

cat("\n")
print("Before calculating totals, let's look at the column total_explants in randomization datasheet")
print(head(randomization_datasheet$total_explants))
print(levels(factor(randomization_datasheet$total_explants)))
cat("\n")


# Remove missing grid positions from output -------------------------------

## This specific step is needed for the situation where a grid item is removed from total explants
## But there's still something there (due to overlap) so total % phenotype can be 
## incorrectly counted as 1.

missing_explant_data_tidy$filename <- gsub("_rgb", "", missing_explant_data_tidy$filename)
missing_explant_data_tidy$ID <- rep(NA, nrow(missing_explant_data_tidy))
for(i in 1:nrow(missing_explant_data_tidy)){
  missing_explant_data_tidy$ID[i] <- parse_trayplateID(name_being_parsed = missing_explant_data_tidy$filename[i])
}

IDs_to_drop <- missing_explant_data_tidy[which(missing_explant_data_tidy$present!="NC"),]
#IDs_to_drop <- missing_explant_data_tidy[which(missing_explant_data_tidy$present!="P"),]
IDs_to_drop$ID_exp <- paste0(IDs_to_drop$ID, "_", IDs_to_drop$grid_item)
output$ID_exp <- paste0(output$ID, "_", output$grid_item)

'%!in%' <- function(x,y)!('%in%'(x,y))
cat("\n")
print(paste0("Nrow of output before dropping missing explants is: ", nrow(output)))
print("Head of ID_exp columns from both output data table and IDs_to_drop data tbales, respectively:")
print(head(output$ID_exp))
print(head(IDs_to_drop$ID_exp))
output <- output[output$ID_exp %!in% IDs_to_drop$ID_exp,]
print(paste0("Nrow of output AFTER dropping missing explants is: ", nrow(output)))
print("Head of output")
print(head(output, 1))
cat("\n")

print(head(IDs_to_drop))
print(colnames(output))
#stop()

# Calculate transgenic ----------------------------------------------------

## First calculate all transgenic
# Initialize new columns
for (j in 1:nrow(pixel_demographics)){
    new_column_name <- paste0('n_transgenic_', pixel_demographics$Tissue[j])
    randomization_datasheet <- cbind(randomization_datasheet, rep(NA, nrow(randomization_datasheet)))
    colnames(randomization_datasheet)[ncol(randomization_datasheet)] <- new_column_name
    
    new_column_name <- paste0('portion_transgenic_', pixel_demographics$Tissue[j])
    randomization_datasheet <- cbind(randomization_datasheet, rep(NA, nrow(randomization_datasheet)))
    colnames(randomization_datasheet)[ncol(randomization_datasheet)] <- new_column_name
}

for (i in 1:nrow(randomization_datasheet)){
    for (j in 1:nrow(pixel_demographics)){
        column_of_interest <- paste0('n_transgenic_', pixel_demographics$Tissue[j])
        total_transgenic <- sum(na.omit(output$transgenic[which(output$ID==randomization_datasheet$ID[i] & output$segment_hex == pixel_demographics$Tissue[j])]))
        # print(paste0('Total transgenic explants on plate ',
        #             randomization_datasheet$ID[i],
        #             ' and tissue ',
        #             pixel_demographics$Tissue[j],
        #             ' is ',
        #             total_transgenic))
        randomization_datasheet[i, eval(column_of_interest)] <- total_transgenic 
        
        column_of_interest <- paste0('portion_transgenic_', pixel_demographics$Tissue[j])
        total_transgenic <- sum(na.omit(output$transgenic[which(output$ID==randomization_datasheet$ID[i] & output$segment_hex == pixel_demographics$Tissue[j])]))
        # print(paste0('Total transgenic explants on plate ',
        #             randomization_datasheet$ID[i],
        #             ' and tissue ',
        #             pixel_demographics$Tissue[j],
        #             ' is ',
        #             total_transgenic))
        
        # Adding this patch because otherwise we get all NaN for min and max of each trait in debugging stdout (debugging needs only...?)
        if(randomization_datasheet$total_explants[i]==0){
          randomization_datasheet[i, eval(column_of_interest)] <- NA
        }
        if(randomization_datasheet$total_explants[i]>0){
          randomization_datasheet[i, eval(column_of_interest)] <- total_transgenic / randomization_datasheet$total_explants[i]
        }
        
        #randomization_datasheet[i, eval(column_of_interest)] <- total_transgenic / randomization_datasheet$total_explants[i]
    }
}
print(paste0("Maximum # grid positions with transgenic stem in any plate: ", max(na.omit(randomization_datasheet$n_transgenic_Stem))))
print(paste0("Maximum # grid positions with transgenic callus in any plate: ", max(na.omit(randomization_datasheet$n_transgenic_Callus))))
print(paste0("Maximum # grid positions with transgenic shoot in any plate: ", max(na.omit(randomization_datasheet$n_transgenic_Shoot))))





# Calculate total ---------------------------------------------------------
## Now calculate all of each tissue whether transgenic or not
# Initialize new columns
for (j in 1:nrow(pixel_demographics)){
    new_column_name <- paste0('n_', pixel_demographics$Tissue[j])
    randomization_datasheet <- cbind(randomization_datasheet, rep(NA, nrow(randomization_datasheet)))
    colnames(randomization_datasheet)[ncol(randomization_datasheet)] <- new_column_name
    
    new_column_name <- paste0('portion_', pixel_demographics$Tissue[j])
    randomization_datasheet <- cbind(randomization_datasheet, rep(NA, nrow(randomization_datasheet)))
    colnames(randomization_datasheet)[ncol(randomization_datasheet)] <- new_column_name
}

for (i in 1:nrow(randomization_datasheet)){
    for (j in 1:nrow(pixel_demographics)){
        column_of_interest <- paste0('n_', pixel_demographics$Tissue[j])
        total_tissue <- sum(na.omit(output$segment_present[which(output$ID==randomization_datasheet$ID[i] & output$segment_hex == pixel_demographics$Tissue[j])]))
        randomization_datasheet[i, eval(column_of_interest)] <- total_tissue
        if(total_tissue>20){
          cat("\n")
          cat("\n")
          print("There are more than 20 grid positions with this tissue on this plate? This can't be right.")
          print(paste0("This ID is : ", randomization_datasheet$ID[i]))
          print(paste0("This segment_hex is : ", pixel_demographics$Tissue[j]))
          print(paste0("The rows we are summing over from output to calculate this obviously wrong statistic:"))
          print(output$segment_present[which(output$ID==randomization_datasheet$ID[i] & output$segment_hex == pixel_demographics$Tissue[j])])
          
        }
        
        column_of_interest <- paste0('portion_', pixel_demographics$Tissue[j])
        total_tissue <- sum(na.omit(output$segment_present[which(output$ID==randomization_datasheet$ID[i] & output$segment_hex == pixel_demographics$Tissue[j])]))
        randomization_datasheet[i, eval(column_of_interest)] <- total_tissue / randomization_datasheet$total_explants[i]
        

    }
}

print(paste0("Minimum # grid positions with explants in any plate: ", min(na.omit(randomization_datasheet$total_explants))))
print(paste0("Maximum # grid positions with explants in any plate: ", max(na.omit(randomization_datasheet$total_explants))))
print(paste0("Maximum # grid positions with callus in any plate: ", max(na.omit(randomization_datasheet$n_Callus))))
print(paste0("Maximum # grid positions with shoot in any plate: ", max(na.omit(randomization_datasheet$n_Shoot))))
if(max(na.omit(randomization_datasheet$n_Callus))>20){
  stop("Error! Why are there more than 20 calli on a plate? Are there two images per sample?")
}

# We need to replace NaN with NA if we want to be able to look at minimum and maximum rates in our debugging output (needed for anything else too...?)
#setDT(randomization_datasheet)[, lapply(.SD, function(x) ifelse(is.nan(x), NA, x))] # https://stackoverflow.com/questions/25013446/how-to-turn-nans-in-a-data-frame-into-nas
#randomization_datasheet[randomization_datasheet=="NaN"] <- NA
#randomization_datasheet[sapply(randomization_datasheet, is.nan)] <- NA

for (j in 1:nrow(pixel_demographics)){
  column_of_interest <- paste0('portion_transgenic_', pixel_demographics$Tissue[j])
  print(paste0("Column of interest is ", column_of_interest))
  print(paste0("For transgenic tissue, What is the range of values for portion for ", pixel_demographics$Tissue[j]))
  print(min(randomization_datasheet[, eval(column_of_interest)]),
        NA.rm = TRUE)
  print(max(randomization_datasheet[, eval(column_of_interest)]),
        NA.rm = TRUE)
  
  print(paste0("What are factor levels?"))
  print(levels(factor(randomization_datasheet[, eval(column_of_interest)])))
  cat("\n")
  
  
  column_of_interest <- paste0('portion_', pixel_demographics$Tissue[j])
  print(paste0("Column of interest is ", column_of_interest))
  print(paste0("For total tissue, What is the range of values for portion for ", pixel_demographics$Tissue[j]))
  print(min(randomization_datasheet[, eval(column_of_interest)]),
        NA.rm = TRUE)
  print(max(randomization_datasheet[, eval(column_of_interest)]),
        NA.rm = TRUE)
  
  print(paste0("What are factor levels?"))
  print(levels(factor(randomization_datasheet[, eval(column_of_interest)])))
  cat("\n")
}

# Remove the "Treatment" Row that is in some data
#randomization_datasheet <- randomization_datasheet[!which(randomization_datasheet$Treatment=="Treatment")]

## Make plots
print("N rows of randomization datasheet: ")
print(nrow(randomization_datasheet))

randomization_datasheet <- randomization_datasheet[!is.na(randomization_datasheet$Genotype_ID),]

print("N rows of randomization datasheet: ")
print(nrow(randomization_datasheet))

cat("\n")
print("Plots will be saved to: ")
print(wd)
cat("\n")


# Calculate size of plot output -------------------------------------------

## Before we make any plots, let's make sure their size will be ok
## We want a plot size that gives 1/4in per genotype*treatment

if(opt$width == -9){
  n_genotypes <- length(levels(factor(randomization_datasheet$Genotype_ID)))
  n_treatments <- length(levels(factor(randomization_datasheet$Treatment)))
  plot_horz_in <- max(c(n_genotypes*n_treatments/3),
                      7) # We want the plots to be no smaller than 7 in wide to allow for titles 
}
if(opt$width != -9){
  plot_horz_in <- opt$width
}

if(plot_horz_in > 10){
  angle <- 75
}
if(plot_horz_in <= 10){
  angle <- 0
}

# Sort genotypes for plots ------------------------------------------------

if(opt$sort==1){
  colnames(randomization_datasheet) <- gsub("Treatment name", "Treatment", colnames(randomization_datasheet))
  library(tidyr)
  
  rd_to_melt <- data.frame(cbind(randomization_datasheet$Genotype_ID,
                                 randomization_datasheet$Treatment,
                                 randomization_datasheet$portion_transgenic_Callus))
  colnames(rd_to_melt) <- c("Genotype_ID", "Treatment", "portion_transgenic_Callus")
  rd_to_melt$portion_transgenic_Callus <- as.numeric(as.character(rd_to_melt$portion_transgenic_Callus))
  
  rd_to_melt_aggregated <- aggregate(portion_transgenic_Callus ~ Genotype_ID + Treatment,
                                     data = rd_to_melt,
                                     FUN = mean,
                                     drop = TRUE)
  
  rd_to_melt_aggregated_spread <- spread(rd_to_melt_aggregated,
                                         key = "Treatment",
                                         value = "portion_transgenic_Callus")
  
  rd_to_melt_aggregated_spread$difference <- rd_to_melt_aggregated_spread[,ncol(rd_to_melt_aggregated_spread)] - rd_to_melt_aggregated_spread[,2]
  
  ordered_genotypes <- rd_to_melt_aggregated_spread[order(rd_to_melt_aggregated_spread$difference),]$Genotype_ID
  
  randomization_datasheet$Genotype_ID <- factor(randomization_datasheet$Genotype_ID,
                                                levels = ordered_genotypes)
  colnames(randomization_datasheet) <- gsub("Treatment", "Treatment name", colnames(randomization_datasheet))
}


# All_tissue portion plots (2D1 and 2D2) ----------------------------------------------------
ggplot(randomization_datasheet, aes(x=`Treatment name`, y=portion_All_tissue, group=`Treatment name`)) + 
  geom_boxplot(outlier.shape=NA,
               fill=randomColor(luminosity=c("bright")),
               color="white",
               alpha=0.2) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  ylab("Grid positions with any tissue") +
  theme_dark() +
  geom_jitter(width=0.50, height=0.001, size = 1) +
  facet_grid(~`Genotype_ID`) +
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  + 
  stat_summary(fun.y=mean, geom="point") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = rel(1.4)),
        axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(size = rel(1.5), vjust=-1),
        axis.title.y = element_text(size = rel(1.3)),
        strip.text = element_text(angle = angle, size=rel(1.3)),
        plot.title = element_text(size=rel(1.7))) +
  ggtitle("Rates of any tissue")#+
#scale_y_continuous(labels = scales::scientific)

ggsave(
  "./LEAF_2D1_portion_All_tissue_total.png",
    plot = last_plot(),   width = plot_horz_in,   height = opt$height,   units = "in")

ggplot(randomization_datasheet, aes(x=`Treatment name`, y=portion_transgenic_All_tissue, group=`Treatment name`)) + 
  geom_boxplot(outlier.shape=NA,
               fill=randomColor(luminosity=c("bright")),
               color="white",
               alpha=0.2) + 
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  theme_dark() + 
  geom_jitter(width=0.50, height=0.001, size = 1) +
  ylab("Grid positions with any transgenic tissue") +
  facet_grid(~`Genotype_ID`) +
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  + 
  stat_summary(fun.y=mean, geom="point") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = rel(1.4)),
        axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(size = rel(1.5), vjust=-1),
        axis.title.y = element_text(size = rel(1.3)),
        strip.text = element_text(angle = angle, size=rel(1.3)),
        plot.title = element_text(size=rel(1.7))) +
  ggtitle("Rates of any transgenic tissue")#+
#scale_y_continuous(labels = scales::scientific)library(randomcoloR)

ggsave(
  "./LEAF_2D2_portion_All_tissue_transgenic.png",
    plot = last_plot(),   width = plot_horz_in,   height = opt$height,   units = "in")







# All_tissue portion plots - GENOTYPES COMBINE (2D1 and 2D2) ----------------------------------------------------
ggplot(randomization_datasheet, aes(x=`Treatment name`, y=portion_All_tissue, group=`Treatment name`)) + 
  geom_violin(outlier.shape=NA,
               fill=randomColor(luminosity=c("bright")),
               color="white",
               alpha=0.2) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  ylab("Grid positions with any tissue") +
  theme_dark() +
  #geom_jitter(width=0.10, height=0.001, size = 1) +
  #facet_grid(~`Genotype_ID`) +
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  + 
  stat_summary(fun.y=mean, geom="point") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = rel(1.4)),
        axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(size = rel(1.5), vjust=-1),
        axis.title.y = element_text(size = rel(1.3)),
        strip.text = element_text(angle = angle, size=rel(1.3)),
        plot.title = element_text(size=rel(1.7))) +
  ggtitle("Rates of any tissue (in leaf explant grid positions)")#+
#scale_y_continuous(labels = scales::scientific)

ggsave(
  "./LEAF_2D1c_portion_All_tissue_total.png",
  plot = last_plot(),   width = opt$width,   height = opt$height,   units = "in")

ggplot(randomization_datasheet, aes(x=`Treatment name`, y=portion_transgenic_All_tissue, group=`Treatment name`)) + 
  geom_violin(outlier.shape=NA,
               fill=randomColor(luminosity=c("bright")),
               color="white",
               alpha=0.2) + 
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  theme_dark() + 
  #geom_jitter(width=0.10, height=0.001, size = 1) +
  ylab("Grid positions with any transgenic tissue") +
  #facet_grid(~`Genotype_ID`) +
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  + 
  stat_summary(fun.y=mean, geom="point") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = rel(1.4)),
        axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(size = rel(1.5), vjust=-1),
        axis.title.y = element_text(size = rel(1.3)),
        strip.text = element_text(angle = angle, size=rel(1.3)),
        plot.title = element_text(size=rel(1.7))) +
  ggtitle("LEAF_Rates of transgenic tissue (in leaf explants)")#+
#scale_y_continuous(labels = scales::scientific)library(randomcoloR)

ggsave(
  "./2D2c_portion_All_tissue_transgenic.png",
  plot = last_plot(),   width = opt$width,   height = opt$height,   units = "in")


# All_tissue sum plots (3D1 and 3D2)  ------------------------------------------------------

ggplot(randomization_datasheet, aes(x=`Treatment name`, y=n_All_tissue, group=`Treatment name`)) + 
  geom_boxplot(outlier.shape=NA,
               fill=randomColor(luminosity=c("bright")),
               color="white",
               alpha=0.2) +
  ylab("# Grid positions with any tissue") +
  theme_dark() +
  geom_jitter(width=0.05, height=0.01, size = 1) +
  facet_grid(~`Genotype_ID`) +
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  + 
  stat_summary(fun.y=mean, geom="point") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = rel(1.4)),
        axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(size = rel(1.5), vjust=-1),
        axis.title.y = element_text(size = rel(1.3)),
        strip.text = element_text(angle = angle, size=rel(1.3)),
        plot.title = element_text(size=rel(1.7))) +
  ggtitle("Rates of any tissue (in leaf explant grid positions)")#+
#scale_y_continuous(labels = scales::scientific)

ggsave(
  "./LEAF_3D1_All_tissue_total.png",
    plot = last_plot(),   width = plot_horz_in,   height = opt$height,   units = "in")

ggplot(randomization_datasheet, aes(x=`Treatment name`, y=n_transgenic_All_tissue, group=`Treatment name`)) + 
  geom_boxplot(outlier.shape=NA,
               fill=randomColor(luminosity=c("bright")),
               color="white",
               alpha=0.2) + 
  theme_dark() + 
  ylab("# Grid positions with any transgenic tissue") +
  geom_jitter(width=0.05, height=0.01, size = 1) +
  facet_grid(~`Genotype_ID`) +
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  + 
  stat_summary(fun.y=mean, geom="point") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = rel(1.4)),
        axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(size = rel(1.5), vjust=-1),
        axis.title.y = element_text(size = rel(1.3)),
        strip.text = element_text(angle = angle, size=rel(1.3)),
        plot.title = element_text(size=rel(1.7))) +
  ggtitle("Rates of transgenic tissue (in leaf explants)")#+
#scale_y_continuous(labels = scales::scientific)library(randomcoloR)

ggsave(
  "./LEAF_3D2_All_tissue_transgenic.png",
    plot = last_plot(),   width = plot_horz_in,   height = opt$height,   units = "in")


# Per explant signal plots ------------------------------------------------

cat("\n")
print("List of treatments in this randomization_datasheet: ")
print(levels(factor(randomization_datasheet$`Treatment name`)))
print("Now about to merge randomization_datasheet and GMOlabeler output.")
print(paste0("Rows in randomization_datasheet: ", nrow(randomization_datasheet)))
print(paste0("Rows in GMOlabeler output: ", nrow(output)))

print("Before merging, look at table of ID in both. ")
print(table(randomization_datasheet$ID))
print(table(output$ID))

combined_data <- merge(output, randomization_datasheet, by="ID", all.x=TRUE)

print(paste0("Rows in merged output: ", nrow(combined_data)))
print("List of treatments in this combined data: ")
print(levels(factor(combined_data$`Treatment name`)))
cat("\n")

combined_data$Genotype_ID <- as.factor(combined_data$Genotype_ID)

for(i in 4){ # ONLY take stats for total tissue
  print(paste0("Making final plots for tissue: "))
  print(pixel_demographics[i, ])
  print("Dim before and after subsetting combined data to this tissue only")
  print(dim(combined_data))
  data_subset <- combined_data[which(combined_data$segment_hex == pixel_demographics$Tissue[i]),]
  print(dim(data_subset))
  print("Max mean_signal")
  print(max(na.omit(data_subset$mean_signal)))
  if (is.infinite(max(na.omit(data_subset$mean_signal)))){
    next
}

  print("Max signal plots")
  
  p2 <- ggplot(data_subset, aes(x=`Treatment name`, y=max_signal, group=`Treatment name`)) + 
    geom_boxplot(outlier.shape=NA,
                 fill=randomColor(luminosity=c("bright")),
                 color="black",
                 alpha=0.2) +
    theme_gray() +
    geom_jitter(width=0.05, height=0.01, size = 1) +
    facet_grid(~`Genotype_ID`) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = rel(1.4)),
          axis.title.x = element_text(size = rel(1.3)),
          axis.text.y = element_text(size = rel(1.5), vjust=-1),
          axis.title.y = element_text(size = rel(1.3)),
          strip.text = element_text(angle = angle, size=rel(1.3)),
          plot.title = element_text(size=rel(1.7))) +
    ylab("Max reporter signal") +
    #stat_summary(fun.y=mean, geom="line", aes(group=1))  + 
    stat_summary(fun.y=mean, geom="point") +
    ggtitle(paste0("Max reporter signal in ", pixel_demographics$Tissue[i]))+
    scale_y_continuous(labels = scales::scientific)
  print(p2)
  
  ggsave(
    paste0("LEAF_B",i,"_", pixel_demographics$Tissue[i],"_Max_signal.png"),
      plot = last_plot(),   width = plot_horz_in,   height = opt$height,   units = "in")
  
  print("Total signal plots")

  p3 <- ggplot(data_subset, aes(x=`Treatment name`, y=total_signal, group=`Treatment name`)) + 
    geom_boxplot(outlier.shape=NA,
                 fill=randomColor(luminosity=c("bright")),
                 color="black",
                 alpha=0.2) +
    theme_gray() +
    geom_jitter(width=0.05, height=0.01, size = 1) +
    facet_grid(~`Genotype_ID`) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = rel(1.4)),
          axis.title.x = element_text(size = rel(1.3)),
          axis.text.y = element_text(size = rel(1.5), vjust=-1),
          axis.title.y = element_text(size = rel(1.3)),
          strip.text = element_text(angle = angle, size=rel(1.3)),
          plot.title = element_text(size=rel(1.7))) +
    ylab("Total reporter signal") +
    #stat_summary(fun.y=mean, geom="line", aes(group=1))  + 
    stat_summary(fun.y=mean, geom="point") +
    ggtitle(paste0("Total reporter signal in ", pixel_demographics$Tissue[i])) +
    scale_y_continuous(labels = scales::scientific)
  print(p3)
  
  ggsave(
    paste0("LEAF_A",i,"_", pixel_demographics$Tissue[i],"_Total_signal.png"),
      plot = last_plot(),   width = plot_horz_in,   height = opt$height,   units = "in")
}



# Write out summary statistics of # transgenic per plate ------------------

message(paste("Writing raw statistics out to",
              paste0(rundir, "/output/", opt$datapath1, "LEAF_plants_over_plates.csv")))

fwrite(randomization_datasheet,
       paste0(rundir, "/output/", opt$datapath1, "LEAF_plants_over_plates.csv")
)
