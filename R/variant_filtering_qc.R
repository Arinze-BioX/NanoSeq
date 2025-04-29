library(UpSetR)
library(ggplot2)

# Get command line arguments, excluding the script name
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments were provided
if (length(args) < 3 || length(args) > 4) {
  stop(paste0("Please provide these 3 positional arguments in this order:",
  " 1. discarded_variant_file, 2. passed_variants_file, 3. output_folder"), call. = FALSE)
}

# Assign arguments to variables (they are initially characters)
discarded_variants <- args[1]  #path to discarded variants CSV
passed_variants <- args[2]  #path to passed variants CSV
out_dir <- args[3]
n_threads <- as.numeric(args[4])

# Check if the output directory exists, if not create it
if (!dir.exists(paste0(out_dir,"/variant_filter_qc"))) {
  dir.create(paste0(out_dir,"/variant_filter_qc"), recursive = TRUE)
}

# Load input csv files
tryCatch({
  discarded_variants <- read.csv(discarded_variants, stringsAsFactors = FALSE)
}, error = function(e) {
  stop("Error reading discarded variants: ", e$message, 
  "\nPlease check the file path and format. Ensure the csv file has headers.")
})

tryCatch({
  passed_variants <- read.csv(passed_variants, stringsAsFactors = FALSE)
}, error = function(e) {
  stop("Error reading discarded variants: ", e$message, 
  "\nPlease check the file path and format. Ensure the csv file has headers.")
})



# Specify Columns of Interest
cols_of_interest <- c("dplx_clip_filter","alignment_score_filter", "mismatch_filter", 
"matched_normal_filter", "duplex_filter", "consensus_base_quality_filter",
"indel_filter", "five_prime_trim_filter", 
"three_prime_trim_filter", "proper_pair_filter", "vaf_filter")


# --- Optional: Validate that columns exist and are binary ---
if (!all(cols_of_interest %in% names(discarded_variants))) {
  stop("One or more QC columns not present in discarded variants csv.")
}

if (!all(cols_of_interest %in% names(discarded_variants))) {
  stop("One or more QC columns not present in discarded variants csv.")
}

# Basic check for 0/1 values (can be adapted if data might have NAs etc.)
# This checks the first 100 rows for simplicity, adjust if needed
# check_data <- data[1:min(nrow(data), 100), cols_of_interest]
# if (!all(sapply(check_data, function(col) all(col %in% c(0, 1, NA))))) {
#    warning("Some specified columns may contain values other than 0 or 1.")
# }


# Function to generate upset plot for discarded variants
plot_discarded_upset <- function(discarded_variants, cols_of_interest) {
  # Body of the function: R code that performs operations
  upset_plot <- upset(
  discarded_variants,
  sets = cols_of_interest,      # The columns (sets) to analyze
  nsets = 11,                   # Number of sets to include (redundant if 'sets' is specified)
  intersections = list(list("dplx_clip_filter","alignment_score_filter", "mismatch_filter", 
  "matched_normal_filter", "duplex_filter", "consensus_base_quality_filter",
  "indel_filter", "five_prime_trim_filter", 
  "three_prime_trim_filter", "proper_pair_filter", "vaf_filter"), 
  list("dplx_clip_filter","alignment_score_filter", "mismatch_filter", 
  "duplex_filter", "consensus_base_quality_filter",
  "indel_filter", "five_prime_trim_filter", 
  "three_prime_trim_filter", "proper_pair_filter"),
  list("dplx_clip_filter","alignment_score_filter", 
  "mismatch_filter", "consensus_base_quality_filter",
  "indel_filter", "five_prime_trim_filter", 
  "three_prime_trim_filter", "proper_pair_filter"),
  list("dplx_clip_filter","alignment_score_filter", "mismatch_filter", 
  "matched_normal_filter", "duplex_filter", "consensus_base_quality_filter",
  "indel_filter", "five_prime_trim_filter", 
  "three_prime_trim_filter", "proper_pair_filter"),
  list("dplx_clip_filter","alignment_score_filter", 
  "matched_normal_filter", "duplex_filter", "consensus_base_quality_filter",
  "indel_filter", "five_prime_trim_filter", 
  "three_prime_trim_filter", "proper_pair_filter", "vaf_filter"),
  list("dplx_clip_filter","mismatch_filter", 
  "matched_normal_filter", "duplex_filter", "consensus_base_quality_filter",
  "indel_filter", "five_prime_trim_filter", 
  "three_prime_trim_filter", "proper_pair_filter", "vaf_filter"),
  list("dplx_clip_filter","mismatch_filter", 
  "matched_normal_filter", "consensus_base_quality_filter",
  "indel_filter", "five_prime_trim_filter", 
  "three_prime_trim_filter", "proper_pair_filter")))

  # Return the result
  return(upset_plot) # Explicit return
}

pdf(paste0(out_dir, "/variant_filter_qc/upset_file_discarded.pdf")) # Save the plot as a PDF
plot_discarded_upset(discarded_variants, cols_of_interest)
dev.off()



#variant qc metrics
normal_coverage_discarded = rowSums(discarded_variants[discarded_variants$commonSNP == 0 & 
discarded_variants$shearwater == 0, 
names(discarded_variants) %in% c("bulkForwardA", "bulkForwardC", "bulkForwardG", 
"bulkForwardT", "bulkForwardIndel", "bulkReverseA", "bulkReverseC", 
"bulkReverseG", "bulkReverseT")])

position_call = discarded_variants[discarded_variants$commonSNP == 0 & 
discarded_variants$shearwater == 0, 
names(discarded_variants) %in% "call"]

call_coverage_normal = rowSums(discarded_variants[discarded_variants$commonSNP == 0 & discarded_variants$shearwater == 0, 
names(discarded_variants) %in% c(paste0("bulkForward",position_call), paste0("bulkReverse",position_call))])

normal_vafs = call_coverage_normal/normal_coverage_discarded
discarded_variants_vaf_df = data.frame(normal_coverage_discarded, 
call_coverage_normal, normal_vafs, 
"Log_coverage"=log10(normal_coverage_discarded))

pdf(paste0(out_dir, "/variant_filter_qc/normal_coverage_plot.pdf"))
ggplot(discarded_variants_vaf_df, aes(x = normal_coverage_discarded)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black", alpha = 0.8) +
  ggtitle("Histogram of Normal Coverage Discarded Vars") +
  xlab("Coverage") +
  ylab("Frequency") +
  theme_minimal() # Using a minimal theme
dev.off()

pdf(paste0(out_dir, "/variant_filter_qc/normal_logcoverage_plot.pdf"))
ggplot(discarded_variants_vaf_df, aes(x = Log_coverage)) +
  geom_histogram(binwidth = 0.2, fill = "skyblue", color = "black", alpha = 0.8) +
  ggtitle("Histogram of Normal Coverage Discarded Vars") +
  xlab("Log Coverage") +
  ylab("Frequency") +
  theme_minimal() # Using a minimal theme
dev.off()

pdf(paste0(out_dir, "/variant_filter_qc/normal_vafs_plot.pdf"))
ggplot(discarded_variants_vaf_df, aes(x = normal_vafs)) +
  geom_histogram(binwidth = 0.02, fill = "skyblue", color = "black", alpha = 0.8) +
  ggtitle("Histogram of total_coverage_snps") +
  xlab("Coverage") +
  ylab("Frequency") +
  theme_minimal() # Using a minimal theme
dev.off()

# Get duplex coverage for genome
total_pos = 0
for (i in 1:n_threads) {
  print(paste0("Processing variants in ", i, ".dsa.bed.gz, going to ", n_threads))
  command_string = paste0("zcat ",out_dir,"/dsa/",i,".dsa.bed.gz | wc -l")
  thread_total_output <- system(command_string, intern = TRUE)
  
  status <- attr(thread_total_output, "status") # Get exit status

  # --- Basic Error Handling ---
  # Check if system command failed
  if (!is.null(status) && status != 0) {
    warning("Command failed for file ", i, ".dsa.bed.gz with status ", status, ": ", command_string)
    next # Skip to the next iteration
  }
   # Check if command returned any output
  if (length(thread_total_output) == 0) {
      warning("Command returned no output for file ", i, ": ", command_string)
      next
  }
  # --- End Error Handling ---


  # Trim whitespace and convert to numeric (take first line of output just in case)
  line_count <- as.numeric(trimws(thread_total_output[1]))

  # Check if conversion to numeric failed
  if (is.na(line_count)) {
    warning("Could not convert line count to numeric for file ", i, 
    ".dsa.bed.gz. Output: '", thread_total_output[1], "'")
    next # Skip to the next iteration
  }

  # Accumulate the count (subtracting 1 as in the original code)
  total_pos <- total_pos + (line_count - 1)
}

#Also get the percent genome coverage
cov_percent = (total_pos / 3300000000) * 100

cov_string = paste0("Total positions in covered by valid duplexes: ", total_pos, "bp. \n", 
"Amounting to: ", cov_percent, "% of the genome\n")
cov_file = paste0(out_dir, "/post/duplex_coverage.txt")
cat(cov_string, file = cov_file, sep = "")