
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

# Function to read and process a single CSV file
process_csv_file <- function(file_path) {
  # Extract sample name from file name
  sample_name <- gsub(".csv$", "", basename(file_path))
  
  # Read the CSV file
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Add sample name column
  data$sample_name <- sample_name
  
  return(data)
}

# Function to calculate overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  ratio <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq, ratio=ratio, rdf=rdf, p=pval)
}

# ---- DOWNSAMPLING SIMULATION ----
simulate_downsampling <- function(count_mat, depths, n_reps) {
    taxa <- rownames(count_mat)
    samples <- colnames(count_mat)
    
    results <- list()
    
    for (depth in depths) {
      for (rep in seq_len(n_reps)) {
        # Rarefy each sample using vegan::rrarefy
        rarefied_list <- lapply(samples, function(s) {
          x <- count_mat[, s]
          x[is.na(x)] <- 0           # replace NAs
          total_reads <- sum(x, na.rm = TRUE)
          if (total_reads == 0) return(rep(0, length(x)))
          as.numeric(rrarefy(matrix(x, nrow = 1), sample = min(depth, total_reads)))
        })
        
        # Combine into matrix
        rarefied <- do.call(cbind, rarefied_list)
        colnames(rarefied) <- samples
        rownames(rarefied) <- taxa
        
        # Diversity metrics
        richness <- apply(rarefied, 2, function(x) sum(x > 0))
        shannon  <- apply(rarefied, 2, function(x) vegan::diversity(x, index = "shannon"))
        
        sim_df <- tibble(
          sample = samples,
          depth = depth,
          rep = rep,
          richness = as.numeric(richness),
          shannon = as.numeric(shannon)
        )
        
        # Join library prep info to simulated results
        sim_df <- sim_df %>%
          left_join(sample_meta, by = c("sample" = "sample_name"))
        
        results[[length(results) + 1]] <- sim_df
      }
    }
    
    bind_rows(results)
}

# Simulate downsampling for specific targets
simulate_downsampling_targets <- function(count_mat, depths, n_reps, patterns) {
  taxa <- rownames(count_mat)
  samples <- colnames(count_mat)
  
  results <- list()
  
  for (depth in depths) {
    for (rep in seq_len(n_reps)) {
      
      # Rarefy each sample
      rarefied_list <- lapply(samples, function(s) {
        x <- count_mat[, s]
        x[is.na(x)] <- 0           # replace NAs
        total_reads <- sum(x)
        if (total_reads == 0) return(rep(0, length(x)))
        as.numeric(vegan::rrarefy(x, sample = min(depth, total_reads)))
      })
      
      rarefied <- do.call(cbind, rarefied_list)
      colnames(rarefied) <- samples
      rownames(rarefied) <- taxa
      
      # Build base results table
      sim_df <- tibble(
        sample = samples,
        depth = depth,
        rep = rep
      ) %>%
        left_join(sample_meta, by = c("sample" = "sample_name"))
      
      # Optional: calculate total reads per target group for manual filtering later
      for (vname in names(patterns)) {
        match_rows <- grepl(patterns[[vname]], rownames(rarefied), ignore.case = TRUE)
        if (any(match_rows)) {
          sim_df[[vname]] <- colSums(rarefied[match_rows, , drop = FALSE])
        } else {
          sim_df[[vname]] <- 0
        }
      }
      
      results[[length(results) + 1]] <- sim_df
    }
  }
  
  bind_rows(results)
}

# ---- SAMPLE METADATA PARSING ----
# Create sample metadata from sequencing metrics
# Parses sample naming conventions across 3 experiments and extracts:
# - library_prep (MSSPE vs NGS)
# - experiment number (1, 2, or 3)
# - base_sample name
# - sample_type (environmental vs slurry)
create_sample_metadata <- function(seq_metrics, pond_samples = c("s9","s10","s19","s20")) {
  seq_metrics %>% 
    mutate(
      # Library prep: MSSPE if ends with "_M"
      library_prep = if_else(str_ends(sample_name, "_M"), "MSSPE", "NGS"),
      
      # Remove the _M suffix to get the core sample identifier
      core = str_remove(sample_name, "_M$"),
      
      # Detect experiment 3 first (Cxxx_LPx pattern)
      is_exp3 = str_detect(core, "C\\d+.*LP\\d+"),
      
      # Base: first sample ID (e.g. S10, PC1)
      base = str_extract(core, "^[A-Za-z]+\\d+"),
      
      # Append -2 or -3 suffix only if not experiment 3
      suffix = case_when(
        !is_exp3 & str_detect(core, "_2$") ~ "-2",
        !is_exp3 & str_detect(core, "_3$") ~ "-3",
        TRUE ~ ""
      ),
      
      # Lowercase base sample name
      base_sample = tolower(paste0(base, suffix)),
      
      # Mark which experiment number
      experiment = case_when(
        is_exp3 ~ 3,
        str_detect(core, "_3$") ~ 2,
        TRUE ~ 1
      ),
      
      # Group slurry vs environmental samples
      sample_type = if_else(base_sample %in% pond_samples, "environmental", "slurry")
    )
}

# ---- VIRAL NAME GROUPING ----
# Add standardized viral pathogen name groups to virus data
# Maps species names to standardized abbreviations (PBoV, PCV, PKV, etc.)
add_viral_name_groups <- function(virus_data) {
  # Define regex patterns for viral species grouping
  ent_v = "Enterovirus G"
  circo_v = "Circovirus"
  sapo_v = paste(c("Sapovirus","Sapporo"), collapse = "|")
  tes_v = "Teschovirus"
  aich_v = "Aichivirus"
  astr_v = paste(c("Porcine astrovirus","Mamastrovirus 3"), collapse = "|")
  hpv18 = "Alphapapillomavirus"
  sapel_v = "Sapelovirus"
  boc_v = paste(c("Porcine bocavirus","Bocaparvovirus"), collapse = "|")
  
  virus_data %>% 
    mutate(name_group = case_when(
      str_detect(name, regex(hpv18, ignore_case = TRUE)) ~ "HPV18",
      str_detect(name, regex(aich_v, ignore_case = TRUE)) ~ "PKV",
      str_detect(name, regex(tes_v, ignore_case = TRUE)) ~ "PTeV",
      str_detect(name, regex(circo_v, ignore_case = TRUE)) ~ "PCV",
      str_detect(name, regex(sapo_v, ignore_case = TRUE)) ~ "PSapeV",
      str_detect(name, regex(sapel_v, ignore_case = TRUE)) ~ "PSapoV",
      str_detect(name, regex(ent_v, ignore_case = TRUE)) ~ "EV-G",
      str_detect(name, regex(boc_v, ignore_case = TRUE)) ~ "PBoV", 
      str_detect(name, regex(astr_v, ignore_case = TRUE)) ~ "PAstV",
      TRUE ~ name  # Keep original name if no pattern matches
    ))
}

  