# This 'module' contains functions to process data from guard cell starch dynamics stats

# nolint start
library('tidyverse') # a collection of R packages designed for data science
library('ggthemes')  # includes colorblind safe color palette
library('rstatix')   # for effect size calculations

load_starch_data <- function(file_path) {
    # Load the starch dynamics results
    # Args:
    #    file_path (chr): xlx file path
    # Returns:
    #    a tibble
    readxl::read_xlsx(file_path, na = 'NA', skip = 0) %>% 
        dplyr::mutate(
            dplyr::across(
                c(`Species`, `Common name`, `Treatment`, `image name`),
                \(x) factor(x)
            )
        )  %>% 
        dplyr::select(
            dplyr::where(is.factor), 
            `starch guard cell area ratio (%)`
        ) %>% 
        dplyr::mutate_if(is.character, as.numeric) %>% 
        tidyr::drop_na() %>% 
        dplyr::filter(
          `Treatment` != 'EON',
          `Species` %in% c('Arabidopsis thaliana', 'Lactuca sativa', 'Hordeum vulgare', 'Brachypodium distachyon')          
          )
}

test_treatment_differences <- function(starch_data) {
  # Performs non-parametric statistical tests between treatments 'RL' and 'RLBL'
  # for each species in the dataset, including effect size
  # 
  # Args:
  #    starch_data (tibble): The loaded starch data from load_starch_data()
  # Returns:
  #    A tibble with species, p-values, and effect sizes
  
  # Function to perform Wilcoxon test for a single species
  test_single_species <- function(species_name) {
    # Filter data for the current species and treatments
    sp_data <- starch_data %>% 
      dplyr::filter(Species == species_name, Treatment %in% c('RL', 'RLBL'))
    
    # Get data for each treatment
    rl_data <- sp_data %>% 
      dplyr::filter(Treatment == 'RL') %>% 
      pull(`starch guard cell area ratio (%)`)
    
    rlbl_data <- sp_data %>% 
      dplyr::filter(Treatment == 'RLBL') %>% 
      pull(`starch guard cell area ratio (%)`)
    
    if (length(rl_data) >= 3 && length(rlbl_data) >= 3) {  # Ensure minimum sample size
      # Perform Wilcoxon test
      test_result <- wilcox.test(
        rl_data, rlbl_data,
        exact = FALSE
      )
      
      # Calculate effect size using rstatix - fixing the issue with NA results
      # Directly calculate effect size without using the formula syntax
      starch_var <- 'starch guard cell area ratio (%)'
      eff_size_result <- tryCatch({
        # Create a simplified dataset for effect size calculation
        test_data <- data.frame(
          value = c(rl_data, rlbl_data),
          group = factor(c(rep('RL', length(rl_data)), rep('RLBL', length(rlbl_data))))
        )
        
        # Use the simplified dataset for effect size calculation
        res <- rstatix::wilcox_effsize(
          data = test_data,
          value ~ group,
          ci = TRUE
        )
        res
      }, error = function(e) {
        # Print the error to help with debugging
        message('Error in effect size calculation: ', e$message)
        return(data.frame(effsize = NA_real_, conf.low = NA_real_, conf.high = NA_real_))
      })
      
      # Safely extract effect size values
      effect_size <- if(is.data.frame(eff_size_result) && 'effsize' %in% names(eff_size_result)) {
        eff_size_result$effsize[1]
      } else {
        NA_real_
      }
      
      effect_size_lower <- if(is.data.frame(eff_size_result) && 'conf.low' %in% names(eff_size_result)) {
        eff_size_result$conf.low[1]
      } else {
        NA_real_
      }
      
      effect_size_upper <- if(is.data.frame(eff_size_result) && 'conf.high' %in% names(eff_size_result)) {
        eff_size_result$conf.high[1]
      } else {
        NA_real_
      }
      
      # Calculate medians
      median_rl <- median(rl_data)
      median_rlbl <- median(rlbl_data)
      
      # Create result row with interpretation
      tibble(
        Species = species_name,
        p_value = test_result$p.value,
        effect_size = effect_size,
        effect_size_lower = effect_size_lower,
        effect_size_upper = effect_size_upper,
        effect_magnitude = dplyr::case_when(
          is.na(effect_size) ~ 'Unknown',
          abs(effect_size) < 0.1 ~ 'Negligible',
          abs(effect_size) < 0.3 ~ 'Small',
          abs(effect_size) < 0.5 ~ 'Medium',
          TRUE ~ 'Large'
        ),
        median_RL = median_rl,
        median_RLBL = median_rlbl,
        n_RL = length(rl_data),
        n_RLBL = length(rlbl_data)
      )
    } else {
      # Return NA if not enough data
      tibble(
        Species = species_name,
        p_value = NA_real_,
        effect_size = NA_real_,
        effect_size_lower = NA_real_,
        effect_size_upper = NA_real_,
        effect_magnitude = 'Insufficient data',
        median_RL = if(length(rl_data) >= 1) median(rl_data) else NA_real_,
        median_RLBL = if(length(rlbl_data) >= 1) median(rlbl_data) else NA_real_,
        n_RL = length(rl_data),
        n_RLBL = length(rlbl_data)
      )
    }
  }
  
  # Get unique species and apply the test to each
  species_levels <- levels(starch_data$Species)
  results <- purrr::map_dfr(species_levels, test_single_species)
  
  # Add multiple testing correction
  results <- results %>%
    # Allow selection of adjustment method
    dplyr::mutate(
      p_adjusted = p.adjust(p_value, method = 'BH'),
      # Add option for different correction methods
      p_bonferroni = p.adjust(p_value, method = 'bonferroni'),
      # Add explanation of adjustment methods
      adjustment_note = paste0(
        'BH controls false discovery rate, recommended for exploratory analyses. ',
        'Bonferroni is more conservative and controls family-wise error rate.'
      ),
      # Add significance based on RAW p-value (not adjusted)
      significance = dplyr::case_when(
        p_value < 0.001 ~ '***',
        p_value < 0.01 ~ '**',
        p_value < 0.05 ~ '*',
        p_value < 0.1 ~ '.',
        TRUE ~ 'ns'
      ),
      # Keep adjusted significance for reference
      significance_adjusted = dplyr::case_when(
        p_adjusted < 0.001 ~ '***',
        p_adjusted < 0.01 ~ '**',
        p_adjusted < 0.05 ~ '*',
        p_adjusted < 0.1 ~ '.',
        TRUE ~ 'ns'
      )
    )
  
  return(results)
}

# Add a function to visualize the results
plot_species_differences <- function(starch_data, test_results) {
  # Create plots for visualizing differences between treatments for each species
  # 
  # Args:
  #   starch_data: The original starch data
  #   test_results: Results from test_treatment_differences function
  # Returns:
  #   A ggplot object showing significance based on raw p-values
  
  # Join significance information with the data
  plot_data <- starch_data %>%
    filter(Treatment %in% c('RL', 'RLBL')) %>%
    left_join(
      test_results %>% 
        select(Species, p_value, effect_magnitude, significance) %>% 
        mutate(
          label = paste0(significance, ' (', effect_magnitude, ')')
        ),
      by = 'Species'
    ) %>%
    # Add stomata shape classification
    mutate(
      stomata_shape = case_when(
        Species %in% c('Arabidopsis thaliana', 'Lactuca sativa') ~ 'Kidney-stomata',
        Species %in% c('Hordeum vulgare', 'Brachypodium distachyon') ~ 'Dumbbell-stomata',
        TRUE ~ 'Kidney-stomata'
      ),
      # Create a combined facet label with species name and stomata shape
      species_with_shape = paste0(Species, '\n(', stomata_shape, ')')
    )
  
  # Explicitly set the order of stomata shapes with kidney first
  plot_data$stomata_shape <- factor(
    plot_data$stomata_shape,
    levels = c('Kidney-stomata', 'Dumbbell-stomata')
  )
  
  # Order species by stomata shape (with kidney first) and then by species name
  plot_data$species_with_shape <- factor(
    plot_data$species_with_shape,
    levels = unique(plot_data$species_with_shape[order(plot_data$stomata_shape, plot_data$Species)])
  )
  
  ggplot(plot_data, aes(x = Treatment, y = `starch guard cell area ratio (%)`, color = Treatment)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    # Use facet_wrap with the combined label
    facet_wrap(~ species_with_shape, scales = 'free_y', ncol = 4) +
    geom_text(
      data = plot_data %>% group_by(species_with_shape) %>% summarize(
        y_pos = max(`starch guard cell area ratio (%)`, na.rm = TRUE) * 1.1,
        label = first(label),
        Treatment = 'RL',
        .groups = 'drop'
      ),
      aes(y = y_pos, label = label),
      color = 'black', size = 3.5
    ) +
    # ylim(0, 100) +
    labs(
      # subtitle = 'Significance: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1, ns p≥0.1 (unadjusted p-values)',
      y = 'Starch area / guard cell area (%)\n'
    ) +
    theme_bw() +
    ggthemes::scale_color_colorblind(
      name = 'Treatment',
      labels = c(
        'RL' = 'Red light 2 h',
        'RLBL' = 'Red + Blue light 1 h'
      )
    ) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = 'black', fill = NA, linewidth = 1),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 20),
      legend.key = ggplot2::element_blank(),
      legend.position = 'bottom',
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(3, 'line'),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 25),
      axis.text.x = ggplot2::element_text(size = 18),
      axis.text.y = ggplot2::element_text(size = 18),
      strip.text = ggplot2::element_text(size = 16, face = 'italic'),
      strip.background = ggplot2::element_rect(fill = 'grey98'),
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 10, r = 20, b = 40, l = 10, unit = 'pt')  # Increased bottom margin
    )
}
# nolint end