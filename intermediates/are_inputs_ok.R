#!/usr/bin/env Rscript

# IFS=','  # Setting the Internal Field Separator to ',' for array joining
#
# echo "opt\$data <- \"$data\""
# echo "opt\$randomization_datasheet <- \"$randomization_datasheet\""
# echo "opt\$segmentation_mode <- \"$segmentation_mode\""
# echo "opt\$unregenerated_tissues <- \"$unregenerated_tissues\""
# echo "opt\$grid <- \"$grid\""
# echo "opt\$missing_explants <- \"$missing_explants\""
# echo "opt\$fluorophores <- \"${fluorophores[*]}\""  # Joining array elements with IFS
# echo "opt\$desired_wavelength_range <- \"${desired_wavelength_range[*]}\""  # Joining array elements with IFS
# echo "opt\$FalseColor_channels <- \"${FalseColor_channels[*]}\""  # Joining array elements with IFS
# echo "opt\$FalseColor_caps <- \"${FalseColor_caps[*]}\""  # Joining array elements with IFS
# echo "opt\$reporters <- \"${reporters[*]}\""  # Joining array elements with IFS
# echo "opt\$pixel_threshold <- \"$pixel_threshold\""
# echo "opt\$reporter_threshold <- \"$reporter_threshold\""
# echo "opt\$segmentation_model_key <- \"$segmentation_model_key\""
# echo "opt\$segmentation_model_path <- \"$segmentation_model_path\""
# echo "opt\$gmodetector_wd <- \"$gmodetector_wd\""
# echo "opt\$spectral_library_path <- \"$spectral_library_path\""
# echo "opt\$deeplab_path <- \"$deeplab_path\""
# echo "opt\$cubeml_path <- \"$cubeml_path\""
# echo "opt\$alignment_path <- \"$alignment_path\""
# echo "opt\$gmolabeler_path <- \"$gmolabeler_path\""
# echo "opt\$contamination_path <- \"$contamination_path\""
# echo "opt\$data_prefix <- \"$data_prefix\""
# echo "opt\$output_directory_prefix <- \"$output_directory_prefix\""
# echo "opt\$cwd <- \"$cwd\""
#
# unset IFS  # Resetting IFS back to default

# opt <- lapply(opt, function(x) if(is.character(x)) gsub("/mnt/drives/", "/media/gmobot/", x) else x)
# opt <- lapply(opt, function(x) if(is.character(x)) gsub("/home/", "~/", x) else x)
# opt <- lapply(opt, function(x) if(is.character(x)) gsub("cubeglm", "gmodetector_py", x) else x)

library(optparse)

# Define command line arguments
option_list <- list(
  make_option(c("--data"), type="character", default="/mnt/Elements_22/GTNEC_GWAS_poplar_transformation_necrotic_test/day1/", help="Data directory path"),
  make_option(c("--randomization_datasheet"), type="character", default="/mnt/Elements_22/GTNEC_GWAS_poplar_transformation_necrotic_test/GTNEC_labels.xlsx", help="Path to randomization datasheet"),
  make_option(c("--segmentation_mode"), type="character", default="hyperspectral", help="Segmentation mode (rgb or hyperspectral)"),
  make_option(c("--unregenerated_tissues"), type="character", default="Background Stem Necrotic", help="Unregenerated tissues"),
  make_option(c("--grid"), type="integer", default=12, help="Grid size (12 or 20)"),
  make_option(c("--missing_explants"), type="character", default="None", help="Missing explants (None, Automatic, or filepath)"),
  make_option(c("--fluorophores"), type="character", default="GFP,Chl,Noise", help="Fluorophores (space-separated)"),
  make_option(c("--desired_wavelength_range"), type="character", default="500,900", help="Desired wavelength range (first last)"),
  make_option(c("--FalseColor_channels"), type="character", default="Chl,GFP,Noise", help="FalseColor channels (Red Green Blue)"),
  make_option(c("--FalseColor_caps"), type="character", default="200,200,200", help="FalseColor caps (Red Green Blue)"),
  make_option(c("--reporters"), type="character", default="GFP,Chl", help="Reporters (space-separated)"),
  make_option(c("--pixel_threshold"), type="integer", default=3, help="Pixel threshold"),
  make_option(c("--reporter_threshold"), type="integer", default=38, help="Reporter threshold"),
  make_option(c("--segmentation_model_key"), type="character", default="/home/models/poplar_training_a2_v7.key.csv", help="Segmentation model key path"),
  make_option(c("--segmentation_model_path"), type="character", default="/home/models/poplar_model_a2_v7_GBC.pkl", help="Segmentation model path"),
  make_option(c("--gmodetector_wd"), type="character", default="/home/cubeglm/", help="GMODetector working directory"),
  make_option(c("--spectral_library_path"), type="character", default="/home/cubeglm/spectral_library/", help="Spectral library path"),
  make_option(c("--deeplab_path"), type="character", default="/home/gmobot/poplar_model_2_w_contam/", help="DeepLab path"),
  make_option(c("--cubeml_path"), type="character", default="/home/cubeml/", help="CubeML path"),
  make_option(c("--alignment_path"), type="character", default="/home/ImageAlignment/", help="Alignment path"),
  make_option(c("--gmolabeler_path"), type="character", default="/home/GMOlabeler/", help="GMOlabeler path"),
  make_option(c("--contamination_path"), type="character", default="/home/DenseNet", help="Contamination path"),
  make_option(c("--data_prefix"), type="character", default="/mnt/output/", help="Data prefix"),
  make_option(c("--output_directory_prefix"), type="character", default="/mnt/output/gmodetector_out/", help="Output directory prefix"),
  make_option(c("--cwd"), type="character", default="/home/GMOnotebook", help="Current working directory")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Function to check if a file or directory exists
check_existence <- function(path, type="file") {
  if (!file.exists(path)) {
    warning(paste(type, "does not exist:", path))
  }
}

# Check if files and directories exist
check_existence(opt$data, "Directory")
check_existence(opt$randomization_datasheet)
if(opt$segmentation_mode=="hyperspectral"){
  check_existence(opt$segmentation_model_key)
  check_existence(opt$segmentation_model_path)
}
check_existence(opt$gmodetector_wd, "Directory")
check_existence(opt$spectral_library_path, "Directory")
check_existence(opt$deeplab_path, "Directory")
check_existence(opt$cubeml_path, "Directory")
check_existence(opt$alignment_path, "Directory")
check_existence(opt$gmolabeler_path, "Directory")
check_existence(opt$contamination_path, "Directory")
check_existence(opt$data_prefix, "Directory")
check_existence(opt$output_directory_prefix, "Directory")
check_existence(opt$cwd, "Directory")

# Check segmentation_mode
if (!(opt$segmentation_mode %in% c("rgb", "hyperspectral"))) {
  warning("segmentation_mode must be either 'rgb' or 'hyperspectral'")
}

# Check grid
if (!(opt$grid %in% c(12, 20))) {
  warning("grid must be either 12 or 20")
}

# Check missing_explants
if (!(opt$missing_explants %in% c("None", "Automatic")) && !file.exists(opt$missing_explants)) {
  warning("missing_explants must be 'None', 'Automatic', or an existing filepath")
}

# Check fluorophores, FalseColor_channels, and reporters against spectral_library_path
check_spectral_library <- function(names, library_path) {
  for (name in strsplit(names, ",")[[1]]) {
    if (!any(grepl(name, list.files(library_path)))) {
      warning(paste("Name not found in spectral library:", name))
    }
  }
}

check_spectral_library(opt$fluorophores, opt$spectral_library_path)
check_spectral_library(opt$FalseColor_channels, opt$spectral_library_path)
check_spectral_library(opt$reporters, opt$spectral_library_path)

# Check desired_wavelength_range
desired_wavelength_range <- as.integer(strsplit(opt$desired_wavelength_range, ",")[[1]])
if (length(desired_wavelength_range) != 2 || !all(sapply(desired_wavelength_range, is.numeric))) {
  warning("desired_wavelength_range must contain two integers")
}

# Function to check if a string represents a valid integer
is_valid_integer <- function(s) {
  return(grepl("^\\d+$", s))
}

# Parse and check desired_wavelength_range
desired_wavelength_range <- unlist(strsplit(opt$desired_wavelength_range, ","))
if (length(desired_wavelength_range) != 2 || !all(sapply(desired_wavelength_range, is_valid_integer))) {
  warning("desired_wavelength_range must contain two integers")
}

# Parse and check FalseColor_caps
FalseColor_caps <- unlist(strsplit(opt$FalseColor_caps, ","))
if (length(FalseColor_caps) != 3 || !all(sapply(FalseColor_caps, is_valid_integer))) {
  warning("FalseColor_caps must contain three integers")
}

library(readxl)

# Read in the randomization datasheet
randomization_datasheet <- suppressMessages(read_excel(opt$randomization_datasheet))

# Check if the first column is an integer with colname Image#
if (names(randomization_datasheet)[1] != "Image#") {
  warning("The first column must be an integer with column name 'Image#'")
}

# Check if the second column is a string with colname TrayID
if (names(randomization_datasheet)[2] != "Tray_ID") {
  warning("The second column must be a string with column name 'Tray_ID'")
}

# Check for the presence of either "Treatment" or "Treatment name" column
if (!("Treatment" %in% names(randomization_datasheet)) && !("Treatment name" %in% names(randomization_datasheet))) {
  warning("There must be a column named either 'Treatment' or 'Treatment name'")
}

# Check for the presence of "Genotype_ID" column
if (!("Genotype_ID" %in% names(randomization_datasheet))) {
  warning("There must be a column named 'Genotype_ID'")
}

# Check for the presence of "Block" column
if (!("Block" %in% names(randomization_datasheet))) {
  warning("There must be a column named 'Block'")
}

phase_ID <- gsub("\\d+$", "", randomization_datasheet$Tray_ID)
phase_ID <- unique(phase_ID)
if(length(phase_ID) > 1){
  warning("There is more than one phase ID (letter string) in the Tray_ID column.")
  print(phase_ID)
}

files <- list.files(opt$data, pattern = "hdr")
files <- str_split_fixed(files, "_", 2)[,1]
files <- files[!grepl("chroma", files)]
files_phase <- gsub("\\d+$", "", files)
phase_ID.2 <- unique(files_phase)
if(length(phase_ID) > 1){
  warning("There is more than one phase ID (letter string) for files in the folder.")
  print(phase_ID.2)
}

if(phase_ID != phase_ID.2){
  warning("Mismatch between phase ID in randomization datasheet and filenames.")
  print("RD:")
  print(phase_ID)
  print("Files:")
  print(phase_ID.2)
}

message("Check complete!")

