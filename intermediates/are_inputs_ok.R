#!/usr/bin/env Rscript

library(optparse)

# Define command line arguments
option_list <- list(
  make_option(c("-d", "--data"), type="character", default="/mnt/Elements_22/GTNEC_GWAS_poplar_transformation_necrotic_test/day1/", help="Data directory path"),
  make_option(c("-r", "--randomization_datasheet"), type="character", default="/mnt/Elements_22/GTNEC_GWAS_poplar_transformation_necrotic_test/GTNEC_labels.xlsx", help="Path to randomization datasheet"),
  make_option(c("-s", "--segmentation_mode"), type="character", default="hyperspectral", help="Segmentation mode (rgb or hyperspectral)"),
  make_option(c("-u", "--unregenerated_tissues"), type="character", default="Background Stem Necrotic", help="Unregenerated tissues"),
  make_option(c("-g", "--grid"), type="integer", default=12, help="Grid size (12 or 20)"),
  make_option(c("-m", "--missing_explants"), type="character", default="None", help="Missing explants (None, Automatic, or filepath)"),
  make_option(c("-f", "--fluorophores"), type="character", default="GFP Chl Noise", help="Fluorophores (space-separated)"),
  make_option(c("-w", "--desired_wavelength_range"), type="character", default="500 900", help="Desired wavelength range (first last)"),
  make_option(c("-fc", "--FalseColor_channels"), type="character", default="Chl GFP Noise", help="FalseColor channels (Red Green Blue)"),
  make_option(c("-fcap", "--FalseColor_caps"), type="character", default="200 200 200", help="FalseColor caps (Red Green Blue)"),
  make_option(c("-rep", "--reporters"), type="character", default="GFP Chl", help="Reporters (space-separated)"),
  make_option(c("-pt", "--pixel_threshold"), type="integer", default=3, help="Pixel threshold"),
  make_option(c("-rt", "--reporter_threshold"), type="integer", default=38, help="Reporter threshold"),
  make_option(c("-smk", "--segmentation_model_key"), type="character", default="/home/models/poplar_training_a2_v7.key.csv", help="Segmentation model key path"),
  make_option(c("-smp", "--segmentation_model_path"), type="character", default="/home/models/poplar_model_a2_v7_GBC.pkl", help="Segmentation model path"),
  make_option(c("-gwd", "--gmodetector_wd"), type="character", default="/home/cubeglm/", help="GMODetector working directory"),
  make_option(c("-slp", "--spectral_library_path"), type="character", default="/home/cubeglm/spectral_library/", help="Spectral library path"),
  make_option(c("-dlp", "--deeplab_path"), type="character", default="/home/gmobot/poplar_model_2_w_contam/", help="DeepLab path"),
  make_option(c("-cmp", "--cubeml_path"), type="character", default="/home/cubeml/", help="CubeML path"),
  make_option(c("-alp", "--alignment_path"), type="character", default="/home/ImageAlignment/", help="Alignment path"),
  make_option(c("-glp", "--gmolabeler_path"), type="character", default="/home/GMOlabeler/", help="GMOlabeler path"),
  make_option(c("-ctp", "--contamination_path"), type="character", default="/home/DenseNet", help="Contamination path"),
  make_option(c("-dp", "--data_prefix"), type="character", default="/mnt/output/", help="Data prefix"),
  make_option(c("-odp", "--output_directory_prefix"), type="character", default="/mnt/output/gmodetector_out/", help="Output directory prefix"),
  make_option(c("-c", "--cwd"), type="character", default="/home/GMOnotebook", help="Current working directory")
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
check_existence(opt$segmentation_model_key)
check_existence(opt$segmentation_model_path)
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
  for (name in strsplit(names, " ")[[1]]) {
    if (!any(grepl(name, list.files(library_path)))) {
      warning(paste("Name not found in spectral library:", name))
    }
  }
}

check_spectral_library(opt$fluorophores, opt$spectral_library_path)
check_spectral_library(opt$FalseColor_channels, opt$spectral_library_path)
check_spectral_library(opt$reporters, opt$spectral_library_path)

# Check desired_wavelength_range
desired_wavelength_range <- as.integer(strsplit(opt$desired_wavelength_range, " ")[[1]])
if (length(desired_wavelength_range) != 2 || !all(sapply(desired_wavelength_range, is.numeric))) {
  warning("desired_wavelength_range must contain two integers")
}

# Check FalseColor_caps
FalseColor_caps <- as.integer(strsplit(opt$FalseColor_caps, " ")[[1]])
if (length(FalseColor_caps) != 3 || !all(sapply(FalseColor_caps, is.numeric))) {
  warning("FalseColor_caps must contain three integers")
}

# Check pixel_threshold and reporter_threshold
if (!is.numeric(opt$pixel_threshold)) {
  warning("pixel_threshold must be an integer")
}
if (!is.numeric(opt$reporter_threshold)) {
  warning("reporter_threshold must be an integer")
}

library(readxl)

# Read in the randomization datasheet
randomization_datasheet <- read_excel(opt$randomization_datasheet)

# Check if the first column is an integer with colname Image#
if (!is.integer(randomization_datasheet[[1]]) || names(randomization_datasheet)[1] != "Image#") {
  warning("The first column must be an integer with column name 'Image#'")
}

# Check if the second column is a string with colname TrayID
if (!is.character(randomization_datasheet[[2]]) || names(randomization_datasheet)[2] != "TrayID") {
  warning("The second column must be a string with column name 'TrayID'")
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

message("Check complete!")

