
# Add Genus & Family Data 
library(taxize)

Sys.setenv(ENTREZ_KEY = "bc06c001a1756b20c2809290f211984ee007")
Sys.setenv(ENTREZ_EMAIL = "matt@opendream.co.th")

# Robust batch processing function with proper error handling
get_genus_info_batch <- function(tax_ids, batch_size = 50) {
  # Clean input data
  original_length <- length(tax_ids)
  valid_indices <- which(!is.na(tax_ids) & tax_ids != "" & tax_ids != 0)
  valid_tax_ids <- tax_ids[valid_indices]
  unique_tax_ids <- unique(valid_tax_ids)
  
  cat("Processing", length(unique_tax_ids), "unique tax IDs out of", original_length, "total\n")
  
  if (length(unique_tax_ids) == 0) {
    return(data.frame(
      genus_name = rep(NA_character_, original_length),
      genus_id = rep(NA_character_, original_length),
      stringsAsFactors = FALSE
    ))
  }
  
  # Create lookup table for results
  genus_lookup <- data.frame(
    tax_id = unique_tax_ids,
    genus_name = rep(NA_character_, length(unique_tax_ids)),
    genus_id = rep(NA_character_, length(unique_tax_ids)),
    stringsAsFactors = FALSE
  )
  
  # Split into batches
  batches <- split(unique_tax_ids, ceiling(seq_along(unique_tax_ids) / batch_size))
  
  # Process each batch
  for (i in seq_along(batches)) {
    cat("Processing batch", i, "of", length(batches), "- IDs:", min(batches[[i]]), "to", max(batches[[i]]), "\n")
    
    batch_ids <- batches[[i]]
    
    tryCatch({
      # Get classifications for this batch
      classifications <- classification(batch_ids, db = "ncbi")
      
      # Safely extract genus information
      for (j in seq_along(batch_ids)) {
        current_id <- batch_ids[j]
        
        # Find the corresponding classification
        if (j <= length(classifications) && !is.null(classifications[[j]])) {
          class_data <- classifications[[j]]
          
          # Check if classification data is valid
          if (is.data.frame(class_data) && nrow(class_data) > 0) {
            # Look for genus rank
            genus_rows <- class_data[!is.na(class_data$rank) & class_data$rank == "genus", ]
            
            if (nrow(genus_rows) > 0) {
              # Update lookup table
              lookup_index <- which(genus_lookup$tax_id == current_id)
              if (length(lookup_index) == 1) {
                genus_lookup$genus_name[lookup_index] <- as.character(genus_rows$name[1])
                genus_lookup$genus_id[lookup_index] <- as.character(genus_rows$id[1])
              }
            }
          }
        }
      }
      
    }, error = function(e) {
      warning(paste("Error processing batch", i, ":", e$message))
    })
    
    # Rate limiting - be nice to NCBI
    if (i < length(batches)) {
      Sys.sleep(1)
    }
  }
  
  # Map results back to original order including NAs and duplicates
  result <- data.frame(
    genus_name = rep(NA_character_, original_length),
    genus_id = rep(NA_character_, original_length),
    stringsAsFactors = FALSE
  )
  
  # Fill in results for valid tax_ids
  for (i in seq_along(tax_ids)) {
    if (!is.na(tax_ids[i]) && tax_ids[i] != "" && tax_ids[i] != 0) {
      # Find matching result in lookup table
      match_indices <- which(genus_lookup$tax_id == tax_ids[i])
      if (length(match_indices) > 0) {
        result$genus_name[i] <- genus_lookup$genus_name[match_indices[1]]
        result$genus_id[i] <- genus_lookup$genus_id[match_indices[1]]
      }
    }
  }
  
  # Report results
  successful <- sum(!is.na(result$genus_name))
  cat("Successfully retrieved genus information for", successful, "out of", length(valid_tax_ids), "valid tax IDs\n")
  
  return(result)
}