###############################################################################
# COMPREHENSIVE DATA ASSESSMENT SCRIPT
# This script examines all input files for structure, content, and potential issues
###############################################################################

# 1) LOAD LIBRARIES
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
library(ggplot2)
library(dlnm)

# Helper functions for assessment
print_separator <- function() {
  cat("\n", paste(rep("=", 80), collapse=""), "\n\n")
}

print_header <- function(text) {
  print_separator()
  cat(toupper(text), "\n")
  print_separator()
}

assess_vector <- function(x, name) {
  cat("* ", name, ":\n")
  cat("  - Class: ", class(x), "\n")
  cat("  - Length: ", length(x), "\n")
  cat("  - NA Count: ", sum(is.na(x)), " (", round(sum(is.na(x))/length(x)*100, 2), "%)\n", sep="")
  
  if(is.numeric(x)) {
    cat("  - Range: [", min(x, na.rm=TRUE), ", ", max(x, na.rm=TRUE), "]\n", sep="")
    cat("  - Mean: ", mean(x, na.rm=TRUE), "\n", sep="")
    cat("  - Median: ", median(x, na.rm=TRUE), "\n", sep="")
    cat("  - SD: ", sd(x, na.rm=TRUE), "\n", sep="")
  }
  
  if(is.factor(x) || is.character(x)) {
    n_levels <- if(is.factor(x)) length(levels(x)) else length(unique(x))
    cat("  - Unique values: ", n_levels, "\n", sep="")
    
    if(n_levels <= 10) {
      if(is.factor(x)) {
        cat("  - Levels: ", paste(levels(x), collapse=", "), "\n", sep="")
      } else {
        cat("  - Unique values: ", paste(unique(x), collapse=", "), "\n", sep="")
      }
    } else {
      # Show counts of top 5 values
      val_counts <- if(is.factor(x)) {
        sort(table(x), decreasing=TRUE)[1:min(5, n_levels)]
      } else {
        sort(table(x), decreasing=TRUE)[1:min(5, n_levels)]
      }
      cat("  - Top 5 values: \n")
      for(i in 1:length(val_counts)) {
        cat("    * ", names(val_counts)[i], ": ", val_counts[i], " occurrences\n", sep="")
      }
    }
  }
  
  if(is.Date(x) || inherits(x, "POSIXt")) {
    cat("  - Date range: [", min(x, na.rm=TRUE), ", ", max(x, na.rm=TRUE), "]\n", sep="")
    cat("  - Time span: ", difftime(max(x, na.rm=TRUE), min(x, na.rm=TRUE), units="days"), " days\n", sep="")
  }
  
  cat("\n")
}

assess_matrix <- function(x, name) {
  cat("* ", name, ":\n")
  cat("  - Class: ", class(x), "\n")
  cat("  - Dimensions: ", paste(dim(x), collapse=" x "), "\n")
  cat("  - NA Count: ", sum(is.na(x)), " (", round(sum(is.na(x))/prod(dim(x))*100, 2), "%)\n", sep="")
  
  if(is.numeric(x)) {
    cat("  - Range: [", min(x, na.rm=TRUE), ", ", max(x, na.rm=TRUE), "]\n", sep="")
    cat("  - Mean: ", mean(x, na.rm=TRUE), "\n", sep="")
  }
  
  cat("  - Column names: ", if(is.null(colnames(x))) "None" else 
    paste(head(colnames(x), 10), collapse=", "), if(ncol(x) > 10) "..." else "", "\n", sep="")
  
  cat("\n")
}

assess_dataframe <- function(df, name) {
  cat("* ", name, ":\n")
  cat("  - Dimensions: ", nrow(df), " rows x ", ncol(df), " columns\n", sep="")
  cat("  - Total NA Count: ", sum(is.na(df)), " (", round(sum(is.na(df))/(nrow(df)*ncol(df))*100, 2), "%)\n", sep="")
  
  # Column information
  cat("  - Column types:\n")
  col_classes <- sapply(df, class)
  class_counts <- table(unlist(lapply(col_classes, `[`, 1)))
  for(cls in names(class_counts)) {
    cat("    * ", cls, ": ", class_counts[cls], " columns\n", sep="")
  }
  
  # Key date and identifier variables
  date_cols <- names(df)[sapply(df, function(x) inherits(x, "Date") || inherits(x, "POSIXt"))]
  id_cols <- names(df)[grepl("id$|^id|code|region", names(df), ignore.case = TRUE)]
  
  if(length(date_cols) > 0) {
    cat("  - Date columns: ", paste(date_cols, collapse=", "), "\n", sep="")
    
    # Check for temporal range
    for(date_col in date_cols) {
      if(inherits(df[[date_col]], "Date") || inherits(df[[date_col]], "POSIXt")) {
        date_range <- range(df[[date_col]], na.rm=TRUE)
        cat("    * ", date_col, " range: [", date_range[1], ", ", date_range[2], "] (", 
            difftime(date_range[2], date_range[1], units="days"), " days)\n", sep="")
      }
    }
  }
  
  if(length(id_cols) > 0) {
    cat("  - ID/Region columns: ", paste(id_cols, collapse=", "), "\n", sep="")
    
    # Check for unique values in ID columns
    for(id_col in id_cols) {
      n_unique <- length(unique(df[[id_col]]))
      cat("    * ", id_col, ": ", n_unique, " unique values\n", sep="")
    }
  }
  
  # Check for specific columns we know are important from the model
  key_cols <- c("Salmonella", "Tmean", "Rainfall", "Floodwave", "Humidity", "LHD_code")
  present_key_cols <- key_cols[key_cols %in% names(df)]
  
  if(length(present_key_cols) > 0) {
    cat("  - Key variables found:\n")
    
    for(col in present_key_cols) {
      cat("    * ", col, ":\n", sep="")
      col_dat <- df[[col]]
      
      # Basic statistics
      if(is.numeric(col_dat)) {
        cat("      - Range: [", min(col_dat, na.rm=TRUE), ", ", max(col_dat, na.rm=TRUE), "]\n", sep="")
        cat("      - Mean: ", mean(col_dat, na.rm=TRUE), "\n", sep="")
        cat("      - Median: ", median(col_dat, na.rm=TRUE), "\n", sep="")
        cat("      - SD: ", sd(col_dat, na.rm=TRUE), "\n", sep="")
        cat("      - NA count: ", sum(is.na(col_dat)), " (", round(sum(is.na(col_dat))/length(col_dat)*100, 2), "%)\n", sep="")
        
        # Check for zero count
        if(is.integer(col_dat) || all(col_dat == round(col_dat), na.rm=TRUE)) {
          zero_count <- sum(col_dat == 0, na.rm=TRUE)
          cat("      - Zero count: ", zero_count, " (", round(zero_count/length(col_dat)*100, 2), "%)\n", sep="")
        }
      } else if(is.factor(col_dat) || is.character(col_dat)) {
        n_unique <- length(unique(col_dat))
        cat("      - Unique values: ", n_unique, "\n", sep="")
        
        if(n_unique <= 10) {
          val_counts <- sort(table(col_dat), decreasing=TRUE)
          val_props <- round(100 * val_counts / sum(val_counts), 2)
          
          for(i in 1:length(val_counts)) {
            cat("        + ", names(val_counts)[i], ": ", val_counts[i], 
                " (", val_props[i], "%)\n", sep="")
          }
        } else {
          val_counts <- sort(table(col_dat), decreasing=TRUE)[1:5]
          val_props <- round(100 * val_counts / sum(table(col_dat)), 2)
          
          cat("      - Top 5 values:\n")
          for(i in 1:length(val_counts)) {
            cat("        + ", names(val_counts)[i], ": ", val_counts[i], 
                " (", val_props[i], "%)\n", sep="")
          }
        }
      }
    }
  }
  
  # Check for columns with high NA percentage
  na_percentage <- colMeans(is.na(df)) * 100
  high_na_cols <- names(na_percentage[na_percentage > 10])
  
  if(length(high_na_cols) > 0) {
    cat("  - Columns with high NA percentage (>10%):\n")
    for(col in high_na_cols) {
      cat("    * ", col, ": ", round(na_percentage[col], 2), "% NA\n", sep="")
    }
  }
  
  # Check for potential date issues
  if("Date" %in% names(df) && inherits(df$Date, "Date")) {
    date_diffs <- diff(sort(unique(df$Date)))
    if(any(date_diffs > 1)) {
      gap_count <- sum(date_diffs > 1)
      max_gap <- max(date_diffs)
      cat("  - Date gaps found: ", gap_count, " gaps, maximum gap: ", max_gap, " days\n", sep="")
    }
  }
  
  cat("\n")
}

###############################################################################
# 2) ASSESS CASE-CROSSOVER DATA
###############################################################################
print_header("CASE-CROSSOVER DATA ASSESSMENT")

# Try-catch to handle potential loading issues
tryCatch({
  # Load data
  cat("Loading data_casecrossover.RData...\n")
  load("model_input/data_casecrossover.RData")
  
  # Check if the expected object exists
  if(!exists("data_cco")) {
    cat("ERROR: Expected object 'data_cco' not found in the RData file.\n")
    cat("Available objects:", paste(ls(), collapse=", "), "\n")
  } else {
    cat("Successfully loaded 'data_cco' with dimensions", nrow(data_cco), "rows x", ncol(data_cco), "columns.\n\n")
    
    # Assess the dataframe
    assess_dataframe(data_cco, "data_cco")
    
    # Check for duplicated column names
    if(anyDuplicated(names(data_cco))) {
      cat("WARNING: Duplicated column names found in data_cco!\n")
      dup_cols <- names(data_cco)[duplicated(names(data_cco))]
      cat("Duplicated columns:", paste(dup_cols, collapse=", "), "\n\n")
    }
    
    # Check for duplicated rows
    if(anyDuplicated(data_cco)) {
      cat("WARNING: Duplicated rows found in data_cco!\n")
      cat("Number of duplicated rows:", sum(duplicated(data_cco)), "\n\n")
    }
    
    # Check for exposure variables
    exposure_cols <- c("Tmean", "Rainfall", "Floodwave", "Humidity")
    missing_cols <- exposure_cols[!exposure_cols %in% names(data_cco)]
    if(length(missing_cols) > 0) {
      cat("WARNING: Missing expected exposure columns:", paste(missing_cols, collapse=", "), "\n\n")
    } else {
      # Check for NAs in exposure variables
      na_counts <- sapply(data_cco[exposure_cols], function(x) sum(is.na(x)))
      cat("NA counts in exposure variables:\n")
      for(col in exposure_cols) {
        cat("  -", col, ":", na_counts[col], "NAs\n")
      }
      
      # Total rows affected by NA
      na_any_exposure <- sum(is.na(data_cco$Tmean) | is.na(data_cco$Rainfall) | 
                               is.na(data_cco$Floodwave) | is.na(data_cco$Humidity))
      cat("Total rows with any exposure NA:", na_any_exposure, 
          "(", round(na_any_exposure/nrow(data_cco)*100, 2), "%)\n\n")
    }
    
    # Check for LHD_code
    if("LHD_code" %in% names(data_cco)) {
      lhd_counts <- table(data_cco$LHD_code)
      cat("LHD_code distribution:\n")
      for(lhd in names(lhd_counts)) {
        cat("  -", lhd, ":", lhd_counts[lhd], "observations\n")
      }
      cat("\n")
    } else {
      cat("WARNING: Expected column 'LHD_code' not found in data_cco!\n\n")
    }
    
    # Check for strata (or create if not exists)
    if("strata" %in% names(data_cco)) {
      cat("Strata already defined in data_cco.\n")
      cat("Number of unique strata:", length(unique(data_cco$strata)), "\n\n")
    } else {
      # Try to create strata if we have LHD_code and year
      if("LHD_code" %in% names(data_cco)) {
        if("year" %in% names(data_cco)) {
          # Year is available
          strata <- paste(data_cco$LHD_code, data_cco$year, sep=":")
        } else if("Date" %in% names(data_cco) && inherits(data_cco$Date, "Date")) {
          # Extract year from Date
          year <- as.numeric(format(data_cco$Date, "%Y"))
          strata <- paste(data_cco$LHD_code, year, sep=":")
        } else {
          strata <- NULL
          cat("WARNING: Cannot create strata, missing year or Date column.\n\n")
        }
        
        if(!is.null(strata)) {
          cat("Created strata from LHD_code and year.\n")
          cat("Number of unique strata:", length(unique(strata)), "\n\n")
        }
      } else {
        cat("WARNING: Cannot create strata, missing LHD_code column.\n\n")
      }
    }
    
    # Visualize Salmonella counts if available
    if("Salmonella" %in% names(data_cco) && "Date" %in% names(data_cco)) {
      cat("Salmonella statistics:\n")
      cat("  - Range:", range(data_cco$Salmonella, na.rm=TRUE), "\n")
      cat("  - Mean:", mean(data_cco$Salmonella, na.rm=TRUE), "\n")
      cat("  - Zero count:", sum(data_cco$Salmonella == 0, na.rm=TRUE), 
          "(", round(sum(data_cco$Salmonella == 0, na.rm=TRUE)/nrow(data_cco)*100, 2), "%)\n")
      cat("  - NA count:", sum(is.na(data_cco$Salmonella)), "\n\n")
      
      # Try to plot histogram
      if(require(ggplot2)) {
        cat("Creating histograms for Salmonella counts...\n")
        hist_plot <- ggplot(data_cco, aes(x=Salmonella)) + 
          geom_histogram(bins=30) +
          theme_minimal() +
          labs(title="Histogram of Salmonella Counts", x="Count", y="Frequency")
        print(hist_plot)
        cat("\n")
      }
    }
  }
}, error = function(e) {
  cat("ERROR loading or analyzing data_casecrossover.RData:\n")
  cat(conditionMessage(e), "\n")
})

###############################################################################
# 3) ASSESS TIME-SERIES DATA
###############################################################################
print_header("TIME-SERIES DATA ASSESSMENT")

# Try-catch to handle potential loading issues
tryCatch({
  # Load data
  cat("Loading data_timeseries.RData...\n")
  load("model_input/data_timeseries.RData")
  
  # Check if the expected object exists
  if(!exists("data_ts")) {
    cat("ERROR: Expected object 'data_ts' not found in the RData file.\n")
    cat("Available objects:", paste(ls(), collapse=", "), "\n")
  } else {
    cat("Successfully loaded 'data_ts' with dimensions", nrow(data_ts), "rows x", ncol(data_ts), "columns.\n\n")
    
    # Assess the dataframe
    assess_dataframe(data_ts, "data_ts")
    
    # Check for duplicated column names
    if(anyDuplicated(names(data_ts))) {
      cat("WARNING: Duplicated column names found in data_ts!\n")
      dup_cols <- names(data_ts)[duplicated(names(data_ts))]
      cat("Duplicated columns:", paste(dup_cols, collapse=", "), "\n\n")
    }
    
    # Check for exposure variables
    exposure_cols <- c("Tmean", "Rainfall", "Floodwave", "Humidity")
    missing_cols <- exposure_cols[!exposure_cols %in% names(data_ts)]
    if(length(missing_cols) > 0) {
      cat("WARNING: Missing expected exposure columns:", paste(missing_cols, collapse=", "), "\n\n")
    } else {
      # Check for NAs in exposure variables
      na_counts <- sapply(data_ts[exposure_cols], function(x) sum(is.na(x)))
      cat("NA counts in exposure variables:\n")
      for(col in exposure_cols) {
        cat("  -", col, ":", na_counts[col], "NAs\n")
      }
      
      # Total rows affected by NA
      na_any_exposure <- sum(is.na(data_ts$Tmean) | is.na(data_ts$Rainfall) | 
                               is.na(data_ts$Floodwave) | is.na(data_ts$Humidity))
      cat("Total rows with any exposure NA:", na_any_exposure, 
          "(", round(na_any_exposure/nrow(data_ts)*100, 2), "%)\n\n")
    }
    
    # Check for LHD_code
    if("LHD_code" %in% names(data_ts)) {
      lhd_counts <- table(data_ts$LHD_code)
      cat("LHD_code distribution:\n")
      for(lhd in names(lhd_counts)) {
        cat("  -", lhd, ":", lhd_counts[lhd], "observations\n")
      }
      cat("\n")
    } else {
      cat("WARNING: Expected column 'LHD_code' not found in data_ts!\n\n")
    }
    
    # Check for month/year columns
    if("month" %in% names(data_ts) && "year" %in% names(data_ts)) {
      cat("Month and Year columns found in data_ts.\n")
      cat("Year range:", range(data_ts$year), "\n")
      cat("Months present:", paste(sort(unique(data_ts$month)), collapse=", "), "\n\n")
    } else if("Date" %in% names(data_ts) && inherits(data_ts$Date, "Date")) {
      cat("Date column found but month/year not explicitly provided.\n")
      cat("Date range:", as.character(range(data_ts$Date)), "\n")
      
      # Extract month and year
      month <- as.numeric(format(data_ts$Date, "%m"))
      year <- as.numeric(format(data_ts$Date, "%Y"))
      cat("Year range:", range(year), "\n")
      cat("Months present:", paste(sort(unique(month)), collapse=", "), "\n\n")
    } else {
      cat("WARNING: Missing month/year information for time series modeling.\n\n")
    }
    
    # Check temporal coverage
    if("Date" %in% names(data_ts) && inherits(data_ts$Date, "Date")) {
      dates <- sort(unique(data_ts$Date))
      date_diffs <- diff(dates)
      
      if(all(date_diffs == 1)) {
        cat("Continuous daily time series with no gaps.\n\n")
      } else {
        gap_count <- sum(date_diffs > 1)
        max_gap <- max(date_diffs)
        cat("Non-continuous time series with", gap_count, "gaps. Maximum gap:", max_gap, "days.\n\n")
        
        if(gap_count <= 10) {
          gap_indices <- which(date_diffs > 1)
          cat("Gaps found between:\n")
          for(i in gap_indices) {
            cat("  -", dates[i], "and", dates[i+1], "(", date_diffs[i], "days)\n")
          }
          cat("\n")
        }
      }
    }
    
    # Visualize Salmonella counts over time if available
    if("Salmonella" %in% names(data_ts) && "Date" %in% names(data_ts)) {
      cat("Salmonella statistics:\n")
      cat("  - Range:", range(data_ts$Salmonella, na.rm=TRUE), "\n")
      cat("  - Mean:", mean(data_ts$Salmonella, na.rm=TRUE), "\n")
      cat("  - Zero count:", sum(data_ts$Salmonella == 0, na.rm=TRUE), 
          "(", round(sum(data_ts$Salmonella == 0, na.rm=TRUE)/nrow(data_ts)*100, 2), "%)\n")
      cat("  - NA count:", sum(is.na(data_ts$Salmonella)), "\n\n")
      
      # Try to plot time series
      if(require(ggplot2)) {
        if("LHD_code" %in% names(data_ts)) {
          # If multiple regions, just plot a few
          regions <- unique(data_ts$LHD_code)
          if(length(regions) > 6) regions <- regions[1:6]
          
          cat("Creating time series plot for first", length(regions), "regions...\n")
          subset_data <- data_ts %>% filter(LHD_code %in% regions)
          ts_plot <- ggplot(subset_data, aes(x=Date, y=Salmonella, color=factor(LHD_code))) + 
            geom_line() +
            theme_minimal() +
            labs(title="Salmonella Counts Over Time", x="Date", y="Count", color="Region")
        } else {
          cat("Creating overall time series plot...\n")
          ts_plot <- ggplot(data_ts, aes(x=Date, y=Salmonella)) + 
            geom_line() +
            theme_minimal() +
            labs(title="Salmonella Counts Over Time", x="Date", y="Count")
        }
        
        print(ts_plot)
        cat("\n")
      }
    }
  }
}, error = function(e) {
  cat("ERROR loading or analyzing data_timeseries.RData:\n")
  cat(conditionMessage(e), "\n")
})

###############################################################################
# 4) ASSESS SPATIAL ADJACENCY DATA
###############################################################################
print_header("SPATIAL ADJACENCY DATA ASSESSMENT")

# Try-catch to handle potential loading issues
tryCatch({
  # Load data
  cat("Loading spatial_adjacency.RData...\n")
  load("model_input/spatial_adjacency.RData")
  
  # Check for expected objects
  expected_objects <- c("adj_vec", "wts_vec", "num_vec")
  missing_objects <- expected_objects[!sapply(expected_objects, exists)]
  
  if(length(missing_objects) > 0) {
    cat("ERROR: Some expected objects not found in the RData file:\n")
    cat("Missing objects:", paste(missing_objects, collapse=", "), "\n")
    cat("Available objects:", paste(ls(), collapse=", "), "\n\n")
  } else {
    cat("Successfully loaded spatial adjacency data.\n\n")
    
    # Assess each component
    assess_vector(adj_vec, "adj_vec")
    assess_vector(wts_vec, "wts_vec")
    assess_vector(num_vec, "num_vec")
    
    # Check consistency
    if(length(wts_vec) != length(adj_vec)) {
      cat("WARNING: Length mismatch between adj_vec and wts_vec!\n")
      cat("adj_vec length:", length(adj_vec), "\n")
      cat("wts_vec length:", length(wts_vec), "\n\n")
    }
    
    # Count adjacencies per region
    if(all(c("adj_vec", "num_vec") %in% ls())) {
      cat("Adjacency summary:\n")
      cat("  - Total regions:", length(num_vec), "\n")
      cat("  - Total adjacency pairs:", length(adj_vec), "\n")
      cat("  - Range of neighbors per region:", range(num_vec), "\n")
      cat("  - Average neighbors per region:", mean(num_vec), "\n\n")
      
      # Check for regions with few neighbors
      few_neighbors <- which(num_vec <= 1)
      if(length(few_neighbors) > 0) {
        cat("Regions with ≤1 neighbor:", length(few_neighbors), "\n")
        cat("Indices:", paste(few_neighbors, collapse=", "), "\n\n")
      }
      
      # Check if adjacency structure is symmetric
      # This is a bit complex to check from the CAR format, but we can do a basic check
      # by counting how many times each region appears in adj_vec
      region_counts <- table(adj_vec)
      expected_counts <- num_vec
      
      # If symmetric, each region should appear in adj_vec exactly as many times as specified in num_vec
      if(length(region_counts) != length(expected_counts) || 
         !all(sort(as.numeric(names(region_counts))) == sort(seq_along(expected_counts)))) {
        cat("WARNING: Possible asymmetry in adjacency structure. Not all regions appear as neighbors.\n\n")
      }
    }
  }
}, error = function(e) {
  cat("ERROR loading or analyzing spatial_adjacency.RData:\n")
  cat(conditionMessage(e), "\n")
})

###############################################################################
# 5) ASSESS MODEL 3 RESULTS IF AVAILABLE
###############################################################################
print_header("MODEL 3 RESULTS ASSESSMENT")

# Try-catch to handle potential loading issues
tryCatch({
  # Check if file exists
  if(!file.exists("model_output/model3_results_nimble.RData")) {
    cat("model3_results_nimble.RData file not found. Skipping assessment.\n")
  } else {
    # Load data
    cat("Loading model3_results_nimble.RData...\n")
    load("model_output/model3_results_nimble.RData")
    
    # Check if the expected object exists
    if(!exists("model3_samples")) {
      cat("ERROR: Expected object 'model3_samples' not found in the RData file.\n")
      cat("Available objects:", paste(ls(), collapse=", "), "\n")
    } else {
      cat("Successfully loaded 'model3_samples' with dimensions", 
          nrow(model3_samples), "MCMC iterations x", ncol(model3_samples), "parameters.\n\n")
      
      # Check basic MCMC structure
      if(!"D" %in% colnames(model3_samples)) {
        cat("WARNING: 'D' (deviance) not found in model3_samples.\n\n")
      } else {
        # Check if D has NA values
        na_count_d <- sum(is.na(model3_samples[, "D"]))
        if(na_count_d > 0) {
          cat("WARNING:", na_count_d, "NA values found in deviance (", 
              round(na_count_d/nrow(model3_samples)*100, 2), "%).\n\n")
        }
      }
      
      # Check for lambda parameters
      lambda_cols <- grep("^lambda\\[", colnames(model3_samples), value=TRUE)
      if(length(lambda_cols) == 0) {
        cat("WARNING: No lambda parameters found in model3_samples.\n")
        cat("This will prevent manual deviance calculation.\n\n")
      } else {
        cat("Found", length(lambda_cols), "lambda parameters.\n")
        
        # Check for NA values in lambda
        lambda_na_count <- sum(is.na(model3_samples[, lambda_cols[1]]))
        if(lambda_na_count > 0) {
          cat("WARNING:", lambda_na_count, "NA values found in first lambda parameter (", 
              round(lambda_na_count/nrow(model3_samples)*100, 2), "%).\n\n")
        }
        
        # Check lambda ranges (just first iteration for efficiency)
        first_iter_lambdas <- model3_samples[1, lambda_cols]
        lambda_range <- range(first_iter_lambdas, na.rm=TRUE)
        cat("Lambda range in first iteration:", lambda_range, "\n\n")
        
        # Look for potential numerical issues
        very_small_lambdas <- sum(first_iter_lambdas < 1e-10, na.rm=TRUE)
        very_large_lambdas <- sum(first_iter_lambdas > 1e10, na.rm=TRUE)
        
        if(very_small_lambdas > 0) {
          cat("WARNING:", very_small_lambdas, "very small lambda values (<1e-10) in first iteration.\n")
        }
        if(very_large_lambdas > 0) {
          cat("WARNING:", very_large_lambdas, "very large lambda values (>1e10) in first iteration.\n")
        }
        cat("\n")
      }
      
      # Check for beta parameters
      beta_temp_cols <- grep("^beta_temp\\[", colnames(model3_samples), value=TRUE)
      if(length(beta_temp_cols) > 0) {
        cat("Found", length(beta_temp_cols), "beta_temp parameters.\n")
        
        # Check for NA values
        beta_na_count <- sum(is.na(model3_samples[, beta_temp_cols[1]]))
        if(beta_na_count > 0) {
          cat("WARNING:", beta_na_count, "NA values found in first beta_temp parameter (", 
              round(beta_na_count/nrow(model3_samples)*100, 2), "%).\n")
        }
        cat("\n")
      }
      
      # Perform basic convergence diagnostics
      if(require(coda)) {
        cat("Performing basic convergence diagnostics...\n")
        
        # Convert to mcmc object for diagnostics
        mcmc_obj <- as.mcmc(model3_samples)
        
        # Effective sample size for key parameters
        key_params <- c("D", 
                        if("D" %in% colnames(model3_samples)) "D" else NULL,
                        if(length(beta_temp_cols) > 0) beta_temp_cols[1] else NULL)
        
        if(length(key_params) > 0) {
          cat("Effective sample sizes for key parameters:\n")
          ess <- effectiveSize(mcmc_obj[, key_params, drop=FALSE])
          for(param in names(ess)) {
            cat("  -", param, ":", round(ess[param]), "\n")
          }
          cat("\n")
        }
      }
    }
  }
}, error = function(e) {
  cat("ERROR loading or analyzing model3_results_nimble.RData:\n")
  cat(conditionMessage(e), "\n")
})

###############################################################################
# 6) ASSESS MODEL 4 RESULTS IF AVAILABLE
###############################################################################
print_header("MODEL 4 RESULTS ASSESSMENT")

# Try-catch to handle potential loading issues
tryCatch({
  # Check if file exists
  if(!file.exists("model_output/model4_results_nimble.RData")) {
    cat("model4_results_nimble.RData file not found. Skipping assessment.\n")
  } else {
    # Load data
    cat("Loading model4_results_nimble.RData...\n")
    load("model_output/model4_results_nimble.RData")
    
    # Check if the expected object exists
    if(!exists("model4_samples")) {
      cat("ERROR: Expected object 'model4_samples' not found in the RData file.\n")
      cat("Available objects:", paste(ls(), collapse=", "), "\n")
    } else {
      cat("Successfully loaded 'model4_samples' with dimensions", 
          nrow(model4_samples), "MCMC iterations x", ncol(model4_samples), "parameters.\n\n")
      
      # Check basic MCMC structure
      if(!"D" %in% colnames(model4_samples)) {
        cat("WARNING: 'D' (deviance) not found in model4_samples.\n\n")
      } else {
        # Check if D has NA values
        na_count_d <- sum(is.na(model4_samples[, "D"]))
        if(na_count_d > 0) {
          cat("WARNING:", na_count_d, "NA values found in deviance (", 
              round(na_count_d/nrow(model4_samples)*100, 2), "%).\n\n")
        }
      }
      
      # Check for lambda parameters
      lambda_cols <- grep("^lambda\\[", colnames(model4_samples), value=TRUE)
      if(length(lambda_cols) == 0) {
        cat("WARNING: No lambda parameters found in model4_samples.\n")
        cat("This will prevent manual deviance calculation.\n\n")
      } else {
        cat("Found", length(lambda_cols), "lambda parameters.\n")
        
        # Check for NA values in lambda
        lambda_na_count <- sum(is.na(model4_samples[, lambda_cols[1]]))
        if(lambda_na_count > 0) {
          cat("WARNING:", lambda_na_count, "NA values found in first lambda parameter (", 
              round(lambda_na_count/nrow(model4_samples)*100, 2), "%).\n\n")
        }
        
        # Check lambda ranges (just first iteration for efficiency)
        first_iter_lambdas <- model4_samples[1, lambda_cols]
        lambda_range <- range(first_iter_lambdas, na.rm=TRUE)
        cat("Lambda range in first iteration:", lambda_range, "\n\n")
        
        # Look for potential numerical issues
        very_small_lambdas <- sum(first_iter_lambdas < 1e-10, na.rm=TRUE)
        very_large_lambdas <- sum(first_iter_lambdas > 1e10, na.rm=TRUE)
        
        if(very_small_lambdas > 0) {
          cat("WARNING:", very_small_lambdas, "very small lambda values (<1e-10) in first iteration.\n")
        }
        if(very_large_lambdas > 0) {
          cat("WARNING:", very_large_lambdas, "very large lambda values (>1e10) in first iteration.\n")
        }
        cat("\n")
      }
      
      # Check for beta parameters
      beta_temp_cols <- grep("^beta_temp\\[", colnames(model4_samples), value=TRUE)
      if(length(beta_temp_cols) > 0) {
        cat("Found", length(beta_temp_cols), "beta_temp parameters.\n")
        
        # Check for NA values
        beta_na_count <- sum(is.na(model4_samples[, beta_temp_cols[1]]))
        if(beta_na_count > 0) {
          cat("WARNING:", beta_na_count, "NA values found in first beta_temp parameter (", 
              round(beta_na_count/nrow(model4_samples)*100, 2), "%).\n")
        }
        cat("\n")
      }
      
      # Check for specific time series parameters
      gamma_cols <- grep("^gamma\\[", colnames(model4_samples), value=TRUE)
      delta_cols <- grep("^delta\\[", colnames(model4_samples), value=TRUE)
      
      if(length(gamma_cols) > 0) {
        cat("Found", length(gamma_cols), "gamma (trend) parameters.\n")
      }
      
      if(length(delta_cols) > 0) {
        cat("Found", length(delta_cols), "delta (seasonal) parameters.\n")
      }
      cat("\n")
      
      # Perform basic convergence diagnostics
      if(require(coda)) {
        cat("Performing basic convergence diagnostics...\n")
        
        # Convert to mcmc object for diagnostics
        mcmc_obj <- as.mcmc(model4_samples)
        
        # Effective sample size for key parameters
        key_params <- c(
          if("D" %in% colnames(model4_samples)) "D" else NULL,
          if(length(beta_temp_cols) > 0) beta_temp_cols[1] else NULL,
          if(length(gamma_cols) > 0) gamma_cols[1] else NULL
        )
        
        if(length(key_params) > 0) {
          cat("Effective sample sizes for key parameters:\n")
          ess <- effectiveSize(mcmc_obj[, key_params, drop=FALSE])
          for(param in names(ess)) {
            cat("  -", param, ":", round(ess[param]), "\n")
          }
          cat("\n")
        }
      }
    }
  }
}, error = function(e) {
  cat("ERROR loading or analyzing model4_results_nimble.RData:\n")
  cat(conditionMessage(e), "\n")
})

###############################################################################
# 7) ASSESS DIC RESULTS IF AVAILABLE
###############################################################################
print_header("DIC RESULTS ASSESSMENT")

# Try-catch to handle potential loading issues
tryCatch({
  # Check if file exists
  if(!file.exists("model_output/dic_results.RData")) {
    cat("dic_results.RData file not found. Skipping assessment.\n")
  } else {
    # Load data
    cat("Loading dic_results.RData...\n")
    load("model_output/dic_results.RData")
    
    # Check if the expected object exists
    if(!exists("dic_results")) {
      cat("ERROR: Expected object 'dic_results' not found in the RData file.\n")
      cat("Available objects:", paste(ls(), collapse=", "), "\n")
    } else {
      cat("Successfully loaded 'dic_results'.\n\n")
      
      # Check structure
      if(!is.list(dic_results)) {
        cat("ERROR: 'dic_results' is not a list as expected.\n")
        cat("Class:", class(dic_results), "\n\n")
      } else {
        # Check for model components
        if(!"model3" %in% names(dic_results)) {
          cat("WARNING: 'model3' component not found in dic_results.\n")
        } else if(is.null(dic_results$model3)) {
          cat("WARNING: 'model3' component is NULL.\n")
        } else {
          cat("Model 3 DIC results:\n")
          for(metric in names(dic_results$model3)) {
            value <- dic_results$model3[[metric]]
            if(is.na(value)) {
              cat("  -", metric, ": NA\n")
            } else {
              cat("  -", metric, ":", round(value, 2), "\n")
            }
          }
          cat("\n")
        }
        
        if(!"model4" %in% names(dic_results)) {
          cat("WARNING: 'model4' component not found in dic_results.\n")
        } else if(is.null(dic_results$model4)) {
          cat("WARNING: 'model4' component is NULL.\n")
        } else {
          cat("Model 4 DIC results:\n")
          for(metric in names(dic_results$model4)) {
            value <- dic_results$model4[[metric]]
            if(is.na(value)) {
              cat("  -", metric, ": NA\n")
            } else {
              cat("  -", metric, ":", round(value, 2), "\n")
            }
          }
          cat("\n")
        }
        
        # Compare models if both results are available
        if(!is.null(dic_results$model3) && !is.null(dic_results$model4)) {
          cat("Model comparison:\n")
          
          if(!is.na(dic_results$model3$DIC) && !is.na(dic_results$model4$DIC)) {
            dic_diff <- dic_results$model3$DIC - dic_results$model4$DIC
            cat("  - DIC difference (Model 3 - Model 4):", round(dic_diff, 2), "\n")
            cat("  - Preferred model (by DIC):", ifelse(dic_diff < 0, "Model 3", "Model 4"), "\n")
          } else {
            cat("  - DIC comparison not possible (NA values)\n")
          }
          
          if(!is.na(dic_results$model3$BIC) && !is.na(dic_results$model4$BIC)) {
            bic_diff <- dic_results$model3$BIC - dic_results$model4$BIC
            cat("  - BIC difference (Model 3 - Model 4):", round(bic_diff, 2), "\n")
            cat("  - Preferred model (by BIC):", ifelse(bic_diff < 0, "Model 3", "Model 4"), "\n")
          } else {
            cat("  - BIC comparison not possible (NA values)\n")
          }
          cat("\n")
        }
      }
    }
  }
}, error = function(e) {
  cat("ERROR loading or analyzing dic_results.RData:\n")
  cat(conditionMessage(e), "\n")
})

###############################################################################
# 8) SUMMARY OF FINDINGS
###############################################################################
print_header("OVERALL ASSESSMENT SUMMARY")

cat("This assessment has examined the structure and content of the following files:\n")
cat("1. model_input/data_casecrossover.RData\n")
cat("2. model_input/data_timeseries.RData\n")
cat("3. model_input/spatial_adjacency.RData\n")
cat("4. model_output/model3_results_nimble.RData (if available)\n")
cat("5. model_output/model4_results_nimble.RData (if available)\n")
cat("6. model_output/dic_results.RData (if available)\n\n")

cat("Key findings and recommendations:\n")

# List potential issues found
issues_found <- FALSE

# Data issues
if(exists("data_cco") && anyDuplicated(names(data_cco))) {
  cat("- ISSUE: data_cco has duplicated column names which should be fixed\n")
  issues_found <- TRUE
}

if(exists("data_ts") && anyDuplicated(names(data_ts))) {
  cat("- ISSUE: data_ts has duplicated column names which should be fixed\n")
  issues_found <- TRUE
}

# Missing exposure columns
if(exists("data_cco") && !all(c("Tmean", "Rainfall", "Floodwave", "Humidity") %in% names(data_cco))) {
  cat("- ISSUE: data_cco is missing some expected exposure columns\n")
  issues_found <- TRUE
}

if(exists("data_ts") && !all(c("Tmean", "Rainfall", "Floodwave", "Humidity") %in% names(data_ts))) {
  cat("- ISSUE: data_ts is missing some expected exposure columns\n")
  issues_found <- TRUE
}

# NA issues in data
if(exists("data_cco")) {
  na_any_exposure_cco <- sum(is.na(data_cco$Tmean) | is.na(data_cco$Rainfall) | 
                               is.na(data_cco$Floodwave) | is.na(data_cco$Humidity))
  if(na_any_exposure_cco > 0) {
    cat("- ISSUE: data_cco has", na_any_exposure_cco, "rows with NA in exposure variables\n")
    issues_found <- TRUE
  }
}

if(exists("data_ts")) {
  na_any_exposure_ts <- sum(is.na(data_ts$Tmean) | is.na(data_ts$Rainfall) | 
                              is.na(data_ts$Floodwave) | is.na(data_ts$Humidity))
  if(na_any_exposure_ts > 0) {
    cat("- ISSUE: data_ts has", na_any_exposure_ts, "rows with NA in exposure variables\n")
    issues_found <- TRUE
  }
}

# Model issues
if(exists("model3_samples") && "D" %in% colnames(model3_samples) && 
   sum(is.na(model3_samples[, "D"])) > 0) {
  cat("- ISSUE: model3_samples has NA values in deviance calculation\n")
  issues_found <- TRUE
}

if(exists("model4_samples") && "D" %in% colnames(model4_samples) && 
   sum(is.na(model4_samples[, "D"])) > 0) {
  cat("- ISSUE: model4_samples has NA values in deviance calculation\n")
  issues_found <- TRUE
}

# Lambda monitoring issues
if(exists("model3_samples") && 
   length(grep("^lambda\\[", colnames(model3_samples), value=TRUE)) == 0) {
  cat("- ISSUE: model3_samples does not include lambda parameters, which are needed for DIC calculation\n")
  issues_found <- TRUE
}

if(exists("model4_samples") && 
   length(grep("^lambda\\[", colnames(model4_samples), value=TRUE)) == 0) {
  cat("- ISSUE: model4_samples does not include lambda parameters, which are needed for DIC calculation\n")
  issues_found <- TRUE
}

# DIC calculation issues
if(exists("dic_results") && !is.null(dic_results$model3) && 
   (is.na(dic_results$model3$DIC) || is.na(dic_results$model3$BIC))) {
  cat("- ISSUE: DIC/BIC calculation for Model 3 resulted in NA values\n")
  issues_found <- TRUE
}

if(exists("dic_results") && !is.null(dic_results$model4) && 
   (is.na(dic_results$model4$DIC) || is.na(dic_results$model4$BIC))) {
  cat("- ISSUE: DIC/BIC calculation for Model 4 resulted in NA values\n")
  issues_found <- TRUE
}

if(!issues_found) {
  cat("- No major issues detected in the examined files.\n")
}

cat("\nRecommendations:\n")
cat("1. Ensure all NA values in exposure variables are properly handled before modeling\n")
cat("2. Include 'lambda' parameters in MCMC monitors for proper deviance calculation\n")
cat("3. Use 'dpois(y[i], lambda[i], log=TRUE)' for deviance calculation in NIMBLE\n")
cat("4. Consider implementing a post-hoc deviance calculation if MCMC-monitored deviance has issues\n")
cat("5. Check MCMC convergence carefully, especially for complex spatial models\n")

print_separator()
cat("Assessment completed.\n")