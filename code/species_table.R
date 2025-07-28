# This 'module' contains functions to process data from Licor LI-6800

# nolint start
library('tidyverse') # a collection of R packages designed for data science
library('ggthemes')  # includes colorblind safe color palette
library('ggtext')  # for styled text elements
library('randomForest') # for random forest
library('xgboost')  # for xgboost
library(ggpubr)
library(janitor)
library('ggpattern')  # Add this



load_species_table <- function(file_path) {
    # Load the species table xlsx file and filter out the parameters of no importance
    # Args:
    #    file_path (chr): xlx file path
    # Returns:
    #    a tibble
    readxl::read_xlsx(file_path, na = 'NA', skip = 1) %>% 
        dplyr::mutate(
            dplyr::across(
                c(
                `Reference`, `Scientific name`, `Common name`, `Crop`, `Category`,
                 `Photosynthetic pathway`, `Climate`
                ),
                \(x) factor(x)
            ),
            `Category` = factor(`Category`,
                levels = c('Angiosperm (Eudicot)', 'Angiosperm (Magnoliids)', 'Angiosperm (Monocot)', 'Fern (Leptosporangiopsida)', 'Gymnosperm')
                ),
            `Photosynthetic pathway` = factor(`Photosynthetic pathway`,
                levels = c('C3', 'C2', 'C4', 'CAM')
                ),
            `Climate` = factor(`Climate`,
                levels = c('Temperate sun', 'Temperate shade', 'Tropical sun', 'Tropical shade')
                ),

        )  %>% 
        dplyr::select(
            dplyr::where(is.factor), 
            `Amax/gswmax (µmol mol-1)`, `A_max (µmol m-2 s-1)`, `gsw_ini (mmol m-2 s-1)`, `gsw_max  (mmol m-2 s-1)`, `gsw_ini / gsw_max`,
            `tA90 (min)`, `tgs90 (min)`, `Low light intensity`, `High light intensity`, `Exp_CO2 (ppm)`, `Exp_Temperature`
        ) %>%
        dplyr::mutate_if(is.character, as.numeric) %>%
        dplyr::filter(
            !is.na(`Reference`),
            `Exp_CO2 (ppm)` < 600
        )
}


t90_plot_crop <- function(species_data) {
    # Plot crop and Arabidopsis tA90 (min) and tgs90 (min)
    # Args:
    #   species_data (tibble): dataset containing species info
    # Returns:
    #   a ggplot object
    ranked <- species_data %>%
      dplyr::filter(Crop == 'Yes' | `Scientific name` == 'Arabidopsis thaliana') %>%
      dplyr::group_by(`Common name`) %>%
      dplyr::summarise(avgTA = mean(`tA90 (min)`, na.rm = TRUE)) %>%
      dplyr::arrange(avgTA)

    species_data %>%
      dplyr::filter(Crop == 'Yes' | `Scientific name` == 'Arabidopsis thaliana') %>%
      dplyr::mutate(
        `Common name` = factor(`Common name`, levels = ranked$`Common name`)
      ) %>%
      tidyr::pivot_longer(
        cols = c('tA90 (min)', 'tgs90 (min)'),
        names_to = 'Parameter',
        values_to = 'Value'
      ) %>%
      dplyr::mutate(
        Parameter = factor(Parameter, 
                          levels = c('tA90 (min)', 'tgs90 (min)'),
                          labels = c("italic(t[A90])~'(min)'", "italic(t[gs90])~'(min)'"))
      ) %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = `Common name`,
          y = Value,
          fill = Parameter
        )
      ) +
      ggplot2::geom_boxplot(
        position = ggplot2::position_dodge(width = 0.8),
        alpha = 0.7
      ) +
      ggplot2::geom_jitter(
        position = ggplot2::position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8),
        alpha = 0.4,
        size = 1.5
      ) +
      ggthemes::scale_fill_colorblind(labels = function(x) parse(text = x)) +
      ggplot2::labs(x = 'Common name', y = 'Time (min)', fill = 'Parameter') +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = 'black', fill = NA, linewidth = 1),
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 20),
        legend.key = ggplot2::element_blank(),
        legend.background = ggplot2::element_blank(),
        legend.box.background = ggplot2::element_blank(),
        legend.key.size = ggplot2::unit(3, 'line'),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggtext::element_markdown(size = 25),
        axis.text.x = ggplot2::element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = ggplot2::element_text(size = 20),
        strip.text = ggplot2::element_text(size = 25),
        strip.background = ggplot2::element_rect(fill = 'grey98'),
        plot.margin = ggplot2::margin(t = 10, r = 10, b = 40, l = 10, unit = 'pt')
      )
}


t90_plot_category <- function(species_data) {
    # Plot plant Category: tA90 (min) and tgs90 (min)
    # Args:
    #   species_data (tibble): dataset containing species info
    # Returns:
    #   a ggplot object
  
  species_data$Category <- as.character(species_data$Category)
  species_data$Category[species_data$Category == 'Fern (Leptosporangiopsida)'] <- 'Modern Fern'  # Replace the specific long name with 'Modern Fern'

  # Apply formatting to each category name
  species_data$Category <- sapply(species_data$Category, function(category) {
    if (grepl('\\(', category)) {
      gsub(' \\(', '\n(', category, fixed = FALSE)  # Replace space before '(' with newline
    } else {
      category
    }
  })

  rankedCat <- species_data %>%
    dplyr::group_by(Category) %>%
    dplyr::summarise(avgTA = mean(`tA90 (min)`, na.rm = TRUE)) %>%
    dplyr::arrange(avgTA)

  species_data %>%
    dplyr::mutate(Category = factor(Category, levels = rankedCat$Category)) %>%
    tidyr::pivot_longer(
      cols = c('tA90 (min)', 'tgs90 (min)'),
      names_to = 'Parameter',
      values_to = 'Value'
    ) %>%
    dplyr::mutate(
      Parameter = factor(Parameter, 
                        levels = c('tA90 (min)', 'tgs90 (min)'),
                        labels = c("italic(t[A90])~'(min)'", "italic(t[gs90])~'(min)'"))
    ) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = Category,
        y = Value,
        fill = Parameter
      )
    ) +
    ggplot2::geom_boxplot(
      position = ggplot2::position_dodge(width = 0.8),
      alpha = 0.7
    ) +
    ggplot2::geom_jitter(
      position = ggplot2::position_jitterdodge(jitter.width = 0.25, dodge.width = 0.8),
      alpha = 0.4,
      size = 1.5
    ) +
    ggthemes::scale_fill_colorblind(labels = function(x) parse(text = x)) +
    ggplot2::labs(x = NULL, y = 'Time (min)', fill = 'Parameter') +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = 'black', fill = NA, linewidth = 1),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 20),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(3, 'line'),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggtext::element_markdown(size = 25),
      axis.text.x = ggplot2::element_text(size = 18, angle = 45, hjust = 0.5, vjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 20),
      strip.text = ggplot2::element_text(size = 25),
      strip.background = ggplot2::element_rect(fill = 'grey98'),
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 40, l = 10, unit = 'pt')  # Increased bottom margin
    )
}


compute_feature_importances <- function(species_data, response_var = 'tgs90_min', top_n = 15) {
    # Compute the feature importances for a given response variable 
    # Args:
    #   species_data (tibble): dataset containing species info
    #   response_var (chr): name of the response variable
    #   top_n (num): number of top features to include
    # Returns:
    #   a list containing data, models, importances, and a ggplot object
    
    # ------------------------------------------------
    # 1. Data Cleaning
    # ------------------------------------------------
    data_clean <- species_data %>%
      dplyr::select(-`Reference`, -`Scientific name`, -`Common name`, -Crop) %>%
      tidyr::drop_na() %>%
      janitor::clean_names(replace = c('\u00b5' = 'u')) %>%
      dplyr::mutate(dplyr::across(where(is.factor), droplevels)) %>%
      dplyr::mutate(dplyr::across(where(is.character), as.factor)) %>%
      dplyr::rename_with(~ str_remove(., '_umol_mol_1$|_umol_m_2_s_1$|_mmol_m_2_s_1$|_ppm$')) %>%
      dplyr::rename(amax_gswmax_ratio = amax_gswmax, gsw_ini_gsw_max_ratio = gsw_ini_gsw_max)
    
    # ------------------------------------------------
    # 2. Prepare data so that every factor is expanded
    #    into dummy variables, and numeric columns are scaled.
    # ------------------------------------------------
    form <- stats::as.formula(paste(response_var, '~ .'))
    mm <- stats::model.matrix(form, data = data_clean)[, -1]
    mm_df <- as.data.frame(mm)
    mm_df[] <- lapply(mm_df, as.numeric)
    
    numeric_cols <- sapply(mm_df, is.numeric)
    mm_df[numeric_cols] <- scale(mm_df[numeric_cols])
    y_scaled <- scale(data_clean[[response_var]])
    
    feature_names <- colnames(mm_df)
    rf_data_list <- list(
      X = mm_df,
      y = as.numeric(y_scaled),
      feature_names = feature_names
    )
    
    # ------------------------------------------------
    # 3. Train Random Forest on the dummy matrix
    # ------------------------------------------------
    set.seed(123)
    rf_model <- randomForest::randomForest(
      x = rf_data_list$X,
      y = rf_data_list$y,
      importance = TRUE,
      ntree = 500
    )
    
    rf_importance_raw <- randomForest::importance(rf_model, type = 1)
    rf_importance_df <- data.frame(
      Feature   = rownames(rf_importance_raw),
      IncMSE    = rf_importance_raw[, 1],
      stringsAsFactors = FALSE
    ) %>%
      dplyr::mutate(Importance = IncMSE / sum(IncMSE)) %>%
      dplyr::arrange(dplyr::desc(Importance))
    
    # ------------------------------------------------
    # 4. Train XGBoost on the same dummy matrix
    # ------------------------------------------------
    xgb_dtrain <- xgboost::xgb.DMatrix(
      data = as.matrix(rf_data_list$X),
      label = rf_data_list$y
    )
    
    set.seed(123)
    xgb_model <- xgboost::xgb.train(
      params = list(
        objective = 'reg:squarederror',
        eta = 0.03,
        max_depth = 4,
        min_child_weight = 3,
        subsample = 0.8,
        colsample_bytree = 0.8,
        gamma = 1
      ),
      data = xgb_dtrain,
      nrounds = 200,
      verbose = 0
    )
    
    xgb_importance <- xgboost::xgb.importance(
      model = xgb_model,
      feature_names = rf_data_list$feature_names
    ) %>%
      dplyr::mutate(Gain = Gain / sum(Gain)) %>%
      dplyr::rename(Importance = Gain) %>%
      dplyr::arrange(dplyr::desc(Importance))
    
    # ------------------------------------------------
    # 5. Compute correlation-based directions
    # ------------------------------------------------
    feature_cor <- sapply(rf_data_list$X, function(col) stats::cor(col, rf_data_list$y))
    direction_df <- data.frame(
      Feature   = names(feature_cor),
      Direction = ifelse(feature_cor > 0, 'Positive', 'Negative'),
      stringsAsFactors = FALSE
    )
    
    # ------------------------------------------------
    # 6. Combine, Sign the Importance, and Plot
    # ------------------------------------------------
    plot_importance_comparison <- function(rf_importance, xgb_importance, direction_df, top_n) {
      rf_data <- rf_importance %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::mutate(Method = 'Random Forest') %>%
        dplyr::left_join(direction_df, by = 'Feature') %>%
        dplyr::mutate(Signed_Importance = Importance * ifelse(Direction == 'Positive', 1, -1))
      
      xgb_data <- xgb_importance %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::mutate(Method = 'XGBoost') %>%
        dplyr::left_join(direction_df, by = 'Feature') %>%
        dplyr::mutate(Signed_Importance = Importance * ifelse(Direction == 'Positive', 1, -1))
      
      # Define intrinsic features - list of original feature names
      intrinsic_features <- c(
        't_a90_min', 'gsw_ini_gsw_max_ratio', 'a_max', 'amax_gswmax_ratio', 
        'gsw_ini', 'gsw_max', 'photosynthetic_pathwayC4', 'photosynthetic_pathwayC2',
        'categoryAngiosperm (Monocot)', 'categoryAngiosperm (Magnoliids)', 
        'categoryFern (Leptosporangiopsida)'
      )
      
      combined_data <- dplyr::bind_rows(rf_data, xgb_data) %>%
        dplyr::group_by(Feature) %>%
        dplyr::mutate(MaxAbsImp = max(abs(Signed_Importance))) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(dplyr::desc(MaxAbsImp)) %>%
        # Add FeatureType column before other operations
        dplyr::mutate(
          # Flag to identify intrinsic vs environmental features
          FeatureType = ifelse(Feature %in% intrinsic_features, 'Intrinsic', 'Environmental'),
          Feature = factor(Feature, levels = unique(Feature))
        ) %>%
        mutate(
          Feature = recode_factor(Feature,
            't_a90_min'             = 't[A90]',  # intrinsic
            'gsw_ini_gsw_max_ratio' = 'g[swini]/g[swmax]', # intrinsic
            'a_max'                 = 'A[max]', # intrinsic
            'amax_gswmax_ratio'     = 'A[max]/g[swmax]', # intrinsic
            'gsw_ini'               = 'g[swini]', # intrinsic
            'gsw_max'               = 'g[swmax]', # intrinsic
            'photosynthetic_pathwayC4' = 'C[4]~photosynthesis', # intrinsic
            'photosynthetic_pathwayC2' = 'C[2]~photosynthesis', # intrinsic
            'exp_co2'                  = 'CO[2]~levels',
            'categoryAngiosperm (Monocot)'       = 'Monocot', # intrinsic
            'climateTropical sun'                = 'Tropical~sun',
            'low_light_intensity'                = 'PPFD[shade]',
            'high_light_intensity'               = 'PPFD[sun]',
            'exp_temperature'                    = 'Temperature',
            'categoryAngiosperm (Magnoliids)'    = 'Magnoliids', # intrinsic
            'categoryFern (Leptosporangiopsida)' = 'Fern') # intrinsic
            )
  
      
      # Create the plot with fixed pattern types
      ggplot2::ggplot(
        combined_data, 
        ggplot2::aes(
          x = Signed_Importance,
          y = forcats::fct_reorder(Feature, abs(Signed_Importance))
        )
      ) +
        ggpattern::geom_col_pattern(
          ggplot2::aes(
            fill = Method,
            pattern = FeatureType
          ),
          position = ggplot2::position_dodge(width = 0.8),
          width = 0.55,
          pattern_density = 0.4,
          pattern_spacing = 0.04,
          pattern_alpha = 0.5
        ) +
        ggplot2::geom_vline(xintercept = 0, linetype = 'dashed', color = 'gray40') +
        ggplot2::scale_y_discrete(labels = function(x) parse(text = x)) +
        ggplot2::scale_x_continuous(labels = abs) +
        ggthemes::scale_fill_colorblind(
          name = 'Methods', 
          breaks = c('XGBoost', 'Random Forest')
        ) +
        ggpattern::scale_pattern_manual(
          name = 'Feature type',
          values = c('Intrinsic' = 'weave', 'Environmental' = 'none'),
          breaks = c('Intrinsic', 'Environmental')
        ) +
        # Try a different approach to separate legends
        guides(
          fill = guide_legend(
            title = 'Methods', 
            override.aes = list(pattern = 'none'),
            keywidth = unit(1.5, 'cm'),
            keyheight = unit(1, 'cm'),
            order = 1  # Make Methods appear first
          ),
          pattern = guide_legend(
            title = 'Feature type', 
            override.aes = list(
              fill = c('white', 'white'),
              pattern_density = 1.0,
              pattern_spacing = 0.02,
              pattern_alpha = 1.0,
              pattern_fill = 'white',
              color = 'black',
              size = 1
            ),
            keywidth = unit(1.5, 'cm'),
            keyheight = unit(1, 'cm'),
            order = 2  # Make Feature type appear second
          )
        ) +
        ggplot2::labs(
            title = 'Drivers of speedy stomata during induction:\nRanked by their feature importances',
            x = expression(atop('', paste('Feature Importance Score on Predicting ', italic(t[gs90])))),
            y = 'Features\n',
            fill = 'Methods'
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_blank(), 
          panel.background = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = 'black', fill = NA, linewidth = 1),
          legend.position = c(0.78, 0.35),
          legend.box = 'vertical',
          legend.spacing.y = unit(1.5, "cm"), # This controls spacing between legends
          legend.title = ggplot2::element_text(size = 23),
          legend.text = ggplot2::element_text(size = 21),
          legend.key = ggplot2::element_blank(),
          legend.background = ggplot2::element_blank(),
          legend.box.background = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_text(size = 25),
          axis.title.y = ggplot2::element_text(size = 25),
          axis.text.x = ggplot2::element_text(size = 20),
          axis.text.y = ggplot2::element_text(size = 20),
          strip.text = ggplot2::element_text(size = 23),
          strip.background = ggplot2::element_rect(fill = 'grey98'),
          plot.margin = ggplot2::margin(t = 10, r = 10, b = 40, l = 10, unit = 'pt')
        )
    }
    
    comparison_plot <- plot_importance_comparison(
      rf_importance_df,
      xgb_importance,
      direction_df,
      top_n = top_n
    )
    
    print(comparison_plot)

    return(list(
      data_clean          = data_clean,
      rf_data_list        = rf_data_list,
      rf_model            = rf_model,
      rf_importance_df    = rf_importance_df,
      xgb_model           = xgb_model,
      xgb_importance      = xgb_importance,
      direction_df        = direction_df,
      comparison_plot     = comparison_plot
    ))
  }


# nolint end