create_list_of_fits <- function(...,
                                list) {
  # Allow for either sequence of fits or list of fits, convert both to list
  fits <- if(missing(list)) {
    base::list(...)
  } else {
    c(base::list(...), list)
  }
  
  # Subset to objects which have the right class, discarding others
  which_gibbsm <- sapply(fits, function(fit) inherits(fit, "gibbsm"))
  if(!all(which_gibbsm)) {
    warning("Provided some arguments that were not fit objects, these were discarded.")
  }
  if(length(which_gibbsm) == 0) {
    stop("No valid fits were provided.")
  }
  fits <- fits[which_gibbsm]
  
  # If user did not give names to the fits, add some default names
  default_names <- paste0("Fit ", seq_len(length(fits)))
  if(is.null(names(fits))) {
    names(fits) <- default_names
  }
  names(fits) <- ifelse(names(fits) == "", default_names, names(fits))
  
  fits
}


access_coefficient <- function(fits,
                               coefficient) {
  
    is_alpha <- gsub("^alpha([0-9]*)", "\\1", coefficient)
    if("" == is_alpha[1] & length(is_alpha) == 1) { # it starts with alpha, but has no integer afterwards
      index <- seq_len(max(sapply(fits, function(fit) {
        length(fit$coefficients$alpha)
      }))) # all the alphas are of interest
      access <- function(obj) obj[["alpha"]]
      identification <- "alpha"
    } else if(length(coefficient) == 1 & is_alpha[1] != coefficient[1]) { # it starts with alpha, with an integer afterwards
      identification <- "alpha"
      access <- function(obj) obj[["alpha"]][as.numeric(is_alpha)]
      index <- is_alpha
    } else {
      is_beta <- gsub("^beta([0-9]*)", "\\1", coefficient)
      identification <- "beta"
      list_covariates <- lapply(fits, function(fit) { # Contains all the covariates
        colnames(as.matrix(fit$coefficients[["beta"]]))
      })
      list_covariates <- as.character(unique(Reduce(c, list_covariates)))
      if("" == is_beta[1] & length(is_beta) == 1) { # it starts with beta, but has no integer afterwards
        index <- list_covariates
      } else if(!all(is.na(suppressWarnings(as.numeric(is_beta))))) { # is_beta is an integer, telling us which beta to access
        index <- as.numeric(is_beta)
      } else if(all(coefficient %in% list_covariates)) {
        index <- Reduce(c, coefficient)
      } else {
        stop("Unrecognised format for the coefficient.")
      }
      
      # If it starts with beta, then this works. If not, assume that coefficient is a covariate and access in this way too.
      access <- function(obj) {
        beta <- as.matrix(obj[["beta"]])
        if(is.character(index)) {
          i <- intersect(index, colnames(beta))
          z <- as.matrix(beta[, i])
          colnames(z) <- i
          rownames(z) <- rownames(beta)
          base::list(z)
        } else {
          base::list(beta)
        }
      }
    }
  
  list(access = access,
       identification = identification,
       index = index)
}




convert_names <- function(x,
                          named_list,
                          ...) {
  possible_names <- list(...)
  if(!is.list(named_list)) {
    stop("named_list should be a named list, with names corresponding to the abbreviated names.")
  }
  sapply(as.character(x), function(ty) {
    if(length(named_list[[ty]]) > 0) {
      return(named_list[[ty]])
    }
    for(possible_names in possible_names) {
      if(!is.list(possible_names)) {
        stop("The extra supplied objects should all be (named) lists.")
      }
      if(length(possible_names[[ty]]) > 0) {
        if(length(named_list[[possible_names[[ty]]]]) > 0) {
          return(named_list[[possible_names[[ty]]]])
        }
      }
    }
    ty
  })
}




#just adding defaults and get rid of how argument because confusing
make_sum_df <- function(fits, #Must be a list of fits
                        coefficient = "alpha",
                        summ, # Must be a list of sums 
                        only_statistically_significant = FALSE,
                        which = "all",
                        full_names = NULL,
                        compute_confidence_intervals = TRUE,
                        classes = NULL,
                        involving = NULL) {
  # Take care of the classes argument
  if(!is.null(classes)) {
    classes <- as.list(classes)
    if(length(fits) > 1) {
      stop("If classes is supplied, then there should be a single fit. Otherwise, colours are used to distinguish fits, not classes.")
    }
  }
  
  # Take care of the full_names argument
  if(!is.null(full_names)) {
    full_names <- as.list(full_names)
    if(length(unique(full_names)) != length(full_names)) {
      stop("full_names does not contain unique names: this will cause some issues later on, please supply distinct names.")
    }
  }
  
  read_coefficient <- access_coefficient(fits = fits, coefficient = coefficient)
  access <- read_coefficient$access
  identification <- read_coefficient$identification
  index <- read_coefficient$index
  
  # If coefficient is not an interaction parameter, need to do extra work
  if(identification == "beta" & which != "all") {
    warning("Plotting regression coefficient but parameter \"which\" was set to something other than \"all\". Assuming this is a typo and setting \"which\" to \"all\".")
    which <- "all"
  }
  
  # Do not compute CIs if average_between
  compute_confidence_intervals <- compute_confidence_intervals & which != "average_between"
  
  if(compute_confidence_intervals) {
    if(!missing(summ)) { # Make sure summaries and fits are compatible
      if(!is(summ, "list")) {
        summ <- base::list(summ)
      }
      stopifnot(length(fits) == length(summ))
    } else { # Construct the summaries
      summ <- lapply(fits, function(f) summary(f))
    }
  }
  
  # Extract list of coefficients, each one corresponding to one of the fits
  estimates <- lapply(fits, function(f) {
    tryCatch(access(f$coefficients), error = function(err) NA)
  })
  
  # Since our estimates are a list of fits, convert df to a list of dataframes, each one corresponding to one of the fits.
  dfs <- lapply(seq_len(length(estimates)), function(k) {
    dfs_by_potential <- lapply(seq_len(length(estimates[[k]])), function(n) {
      # Get the types of the current fit
      types <- rownames(estimates[[k]][[n]])

      # Types involving the focal ones
      types_subset <- if(!is.null(involving)) {
        involving <- Reduce(c, as.list(involving))
        intersect(types, involving)
      } else {
        types
      }

      
      if(identification == "alpha") { # Create dataframe of two columns with all possible pairs of types
        d <- as.data.frame(expand.grid(from = types_subset, to = types_subset, stringsAsFactors = FALSE))
        d <- d[!duplicated(t(apply(d, 1, sort))), ]
      
        
        if(which == "within") { # in this case, remove columns that are the same, so shows only between interactions
          d <- d[d$to == d$from, ]
        } else if(which == "between" | which == "average_between") {
          d <- d[d$to != d$from, ]
        }
      } else if(is.numeric(index)) {
        d <- as.data.frame(Reduce(rbind, lapply(index, function(i) data.frame(from = types_subset, to = i))))
      } else {
        d <- as.data.frame(Reduce(rbind, lapply(intersect(index, colnames(fits[[k]]$coefficients$beta)), function(i) data.frame(from = types_subset, to = i))))
      }
      
      d$E <- sapply(seq_len(nrow(d)), function(i) { #these and below parts specify the values of point and bars
        tryCatch(estimates[[k]][[n]][d$from[i], d$to[i]], error = function(err) NA)
      })
      
      if(compute_confidence_intervals) {
        d$lo <- sapply(seq_len(nrow(d)), function(i) { # Get the lower-endpoint of the CIs, apply function to sequence of numbers 1 to nrow.df
          tryCatch(access(summ[[k]]$lo)[[n]][d$from[i], d$to[i]], error = function(err) NA) # retrieve numbers from summ and apply to the rows
        })
        
        d$hi <- sapply(seq_len(nrow(d)), function(i) { # Get the upper-endpoint of the CIs
          tryCatch(access(summ[[k]]$hi)[[n]][d$from[i], d$to[i]], error = function(err) NA)
        })
        
        d$lo_numerical <- sapply(seq_len(nrow(d)), function(i) { # Lower-endpoint of the numerical CI (relating to the numerical error due to dummy points)
          tryCatch(access(summ[[k]]$lo_numerical)[[n]][d$from[i], d$to[i]], error = function(err) NA)
        })
        
        d$hi_numerical <- sapply(seq_len(nrow(d)), function(i) { # Upper-endpoint of the numerical CI (relating to the numerical error due to dummy points)
          tryCatch(access(summ[[k]]$hi_numerical)[[n]][d$from[i], d$to[i]], error = function(err) NA)
        })
      }
      
      # Convert names to full names
      if(!is.null(full_names)) {
        # We do not yet convert the names, because to decide class we want to look at original names first
        d$new_from <- convert_names(d$from, full_names)
        if(identification != "beta") {
          d$new_to <- convert_names(d$to, full_names)
        }
      }
      
      # If classes was supplied, fill in class
      if(!is.null(classes)) {
        if(!is.null(full_names)) {
          d$class_from <- convert_names(d$from, classes, full_names)
          d$class_to <- convert_names(d$to, classes, full_names)
        } else {
          d$class_from <- convert_names(d$from, classes)
          d$class_to <- convert_names(d$to, classes)
        }
      }
      
      # At this point, we can discard the temporary column names chosen above
      if("new_from" %in% colnames(d)) {
        d$from <- d$new_from
        d$new_from <- NULL
      }
      if("new_to" %in% colnames(d)) {
        d$to <- d$new_to
        d$new_to <- NULL
      }
      
      # Set name of the fit
      if(nrow(d) > 0) {
        d$Fit <- names(fits)[k]
        d$Potential <- paste0("Potential ", index[n])
      }
      
      # Compute average of the coefficient for each type
      if(which == "average_between") {
        d <- as.data.frame(Reduce(rbind, lapply(union(unique(d$from), unique(d$to)), function(ty) {
          g <- d[d$from == ty | d$to == ty, ][1, ]
          g$from <- ty
          g$to <- ty
          g$E <- mean(d$E[d$from == ty | d$to == ty], na.rm = TRUE)
          g
        })))
      }
      
      if(only_statistically_significant) { # in this case, remove non-stat significant
        d <- as.data.frame(d[d$lo > 0 | d$hi < 0, ]) # If non-statistically significant, remove row,. i.e. if low is > 0 or if high < 0, the row is kept
      }
      
      d
    })
    
    Reduce(rbind, dfs_by_potential)
  })
  
  # Flatten dfs
  df <- as.data.frame(Reduce(rbind, dfs))
  
  # Make into factor and order them correctly
  if(!is.null(df$Fit)) {
    df$Fit <- factor(df$Fit, levels = names(fits))
  }
  if(!is.null(df$Potential)) {
    df$Potential <- factor(df$Potential)
  }
  if(!is.null(df$class_from)) {
    df$class_from <- factor(df$class_from)
  }
  if(!is.null(df$class_to)) {
    df$class_to <- factor(df$class_to)
  }
  
  colnames(df)[colnames(df) == "E"] <- identification
  
  df
}


# 
# # 
# # ## Test function
# df1 <- make_sum_df(fits = list("ANU101" = fit_anu101_2, "ANU363" = fit_anu363),
#             summ = list(sum_anu101_2, sum_anu363),
#             coefficient = "alpha")
# 
# 
# fits <- anu101$fit
# sum <- anu101$sum
# 
# df2 <- make_sum_df(fits = list(fits), 
#                    summ = list(sum))
# 
# 
