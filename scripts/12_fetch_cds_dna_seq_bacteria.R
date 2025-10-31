# Print the current working directory
print(getwd())
cat("Starting script execution...\n")

# Clear the working environment
rm(list = ls())
library(Biostrings)

# Set the working directory and output path
setwd("/s3/SHARE/woodman/Prokka2/output/")
output_path <- "/s3/SHARE/woodman/Prokka2/bacteria_cds/12w/"

# Create the output path and log file path
if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}
log_dir <- "/s3/SHARE/woodman/Prokka2/logs"
if (!dir.exists(log_dir)) {
    dir.create(log_dir)
}
error_log_path <- file.path(log_dir, "tsv_error.log")

# Get the list of folders
file_folder <- list.files()

# Process each folder
for (i in 1:length(file_folder)) {
    anno_file_path <- paste0(file_folder[i], "/", file_folder[i], ".tsv")
    out_files <- paste0(output_path, file_folder[i], ".fasta")
    
    # Check if the output file already exists
    if (file.exists(out_files)) {
        cat(file_folder[i], "already processed, skipping.\n")
        next
    }

    # Check if the annotation file exists
    if (file.exists(anno_file_path)) {
        anno_file <- read.delim(anno_file_path, sep = "\t", header = TRUE)
        
        if (dim(anno_file)[1] > 0) {
            cds_anno <- anno_file[which(anno_file$ftype == "CDS"), ]
            if (dim(cds_anno)[1] > 0) {
                cds_id <- paste0(cds_anno$locus_tag, " ", cds_anno$product)
                fa_file_path <- paste0(file_folder[i], "/", file_folder[i], ".ffn")
                
                if (file.exists(fa_file_path)) {
                    fa_file <- readDNAStringSet(fa_file_path)
                    sub_fa_file <- fa_file[cds_id]
                    writeXStringSet(sub_fa_file, out_files, format = "fasta")
                } else {
                    cat(file_folder[i], "\n", file = error_log_path, append = TRUE)
                }
            }
        }
    } else {
        cat(file_folder[i], "\n", file = error_log_path, append = TRUE)
    }
}

cat("Script execution completed.\n")