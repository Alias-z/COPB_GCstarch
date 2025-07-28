# This "module" contains functions to process data from Licor LI-6800

# nolint start
library('tidyverse') # a collection of R packages designed for data science
library('ggthemes')  # includes colorblind safe color palette
library('plantecophys')  # to fit the FVcB model
library('lubridate')  # to work with dates and times
library('ggtext')  # for styled text elements

load_licor_log <- function(file_path) {
    # Load an xlsx file from Licor and filter out the parameters of no importance
    # Args:
    #    file_path (chr): xlx file path
    # Returns:
    #    a tibble
    readxl::read_xlsx(file_path, sheet = 2, na = 'NA', skip = 14) %>% 
        dplyr::slice(-1) %>%
        dplyr::mutate(
            `time` = as.POSIXct(`hhmmss...5`, format = '%H:%M:%S'), 
            dplyr::across(
                c(`Age`, `Species`, `Genotype`, `Replicates`, `Treatment`),
                \(x) factor(x)
            )
        ) %>% 
        dplyr::select(
            dplyr::where(is.factor), 
            `time`, `Fv/Fm`, `Q`, `Q_blue`, `Q_red`, `Ca`, `Ci`, `A`, `gsw`, `E`, `Tleaf`
        ) %>% 
        dplyr::mutate_if(is.character, as.numeric) %>% 
        dplyr::group_by(`Replicates`) %>% 
        dplyr::mutate(
            dplyr::across(
                c(`Q`, `Q_blue`, `Q_red`),
                \(x) round(x, 0)
            ),
            `Genotype_Replicate` = factor(interaction(`Genotype`, `Replicates`, sep = '_')),
            `Time_points` = dplyr::case_when(
                `Treatment` != 'Dark adapted Fv/Fm' ~ {
                    start_time <- min(`time`[`Treatment` == 'RL100 (1)'], na.rm = TRUE)
                    round(as.numeric(difftime(`time`, start_time, units = 'mins')), 0)
                }
            ),
            `A_normalized` = dplyr::case_when(
                `Treatment` != 'Dark adapted Fv/Fm' ~ {
                    A_max <- max(`A`[`Treatment` == 'RLBL1000'], na.rm = TRUE)
                    `A` / A_max
                },
                TRUE ~ NA_real_
            ),
            `gsw_normalized` = dplyr::case_when(
                `Treatment` != 'Dark adapted Fv/Fm' ~ {
                    gsw_max <- max(`gsw`[`Treatment` == 'RLBL1000'], na.rm = TRUE)
                    `gsw` / gsw_max
                },
                TRUE ~ NA_real_
            )
        ) %>%
        # Calculate changes every 30 minutes
        dplyr::mutate(
            time_bin = cut(Time_points, breaks = seq(-30, max(Time_points, na.rm = TRUE) + 30, by = 30)),
            A_change_30min = dplyr::case_when(
                `Treatment` != 'Dark adapted Fv/Fm' ~ {
                    ave(A_normalized, time_bin, FUN = function(x) {
                        if(length(x) >= 2) max(x) - min(x) else NA_real_
                    })
                }
            ),
            gsw_change_30min = dplyr::case_when(
                `Treatment` != 'Dark adapted Fv/Fm' ~ {
                    ave(gsw_normalized, time_bin, FUN = function(x) {
                        if(length(x) >= 2) max(x) - min(x) else NA_real_
                    })
                }
            )
        ) %>%
        dplyr::ungroup()
}

batch_licor_logs <- function(root) {
    # Combine data from all xlsx files
    # Args:
    #    root (chr): the root directory
    # Returns:
    #    a tibble with all xlsx files data combined
    logs <- list.files(
        path = root, 
        pattern = '\\.xlsx$', 
        full.names = TRUE)  # find all xlsx files
    
    logs %>% 
        lapply(load_licor_log) %>% 
        dplyr::bind_rows()  # load all and bind them together
}


load_starch_data <- function(file_path) {
    # Load starch data from Excel and calculate changes between time points
    # Args:
    #    file_path (chr): Excel file path
    # Returns:
    #    a tibble with calculated changes for each time window
    
    data <- readxl::read_xlsx(file_path) %>%
        dplyr::mutate(
            Time_points = factor(Time_points, 
                               levels = c("0", "30", "60", "90", "120", "150", "180", "210"))
        ) %>%
        dplyr::arrange(Species, Time_points) %>%
        dplyr::group_by(Species, Time_points) %>%
        dplyr::summarise(
            starch_ratio = mean(`starch guard cell area ratio (%)`, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        dplyr::group_by(Species) %>%
        dplyr::mutate(
            time_bin = paste0("(", 
                as.numeric(as.character(Time_points)) - 30, ",",
                as.numeric(as.character(Time_points)), "]"),
            starch_change_30min = starch_ratio - dplyr::lag(starch_ratio)
        ) %>%
        dplyr::ungroup()
    
    return(data)
}


batch_load_starch_data <- function(root) {
    # Combine data from all xlsx files with starch measurements
    # Args:
    #    root (chr): the root directory
    # Returns:
    #    a tibble with all xlsx files data combined
    logs <- list.files(
        path = root, 
        pattern = 'starch.*\\.xlsx$',  # find xlsx files with "starch" in name
        full.names = TRUE
    )
    
    logs %>%
        lapply(load_starch_data) %>%
        dplyr::bind_rows()  # load all and bind them together
}



rlbl_gs_curve <- function(`data_rlbl`, `starch_data`, plot_species) {
    # Create data frame for bin centers and changes
    time_bins <- data.frame(
        bin_center = seq(-15, 195, by = 30)
    )
    
    changes_data <- data_rlbl %>%
        dplyr::group_by(time_bin) %>%
        dplyr::summarise(
            A_max_diff = max(A_normalized) - min(A_normalized),
            A_change = ifelse(A_normalized[which.max(Time_points)] > A_normalized[which.min(Time_points)], 
                            A_max_diff, -A_max_diff),
            gsw_max_diff = max(gsw_normalized) - min(gsw_normalized),
            gsw_change = ifelse(gsw_normalized[which.max(Time_points)] > gsw_normalized[which.min(Time_points)], 
                                gsw_max_diff, -gsw_max_diff),
            bin_center = mean(as.numeric(Time_points)),
            Time_points = mean(Time_points)
        ) %>%
        dplyr::filter(!is.na(A_change), Time_points > -30 & Time_points <= 210) %>%
        dplyr::left_join(
            starch_data %>%
                dplyr::filter(Species == plot_species) %>%
                dplyr::select(time_bin, starch_change_30min),
            by = "time_bin"
        )

    treatment_ranges <- data.frame(
        `Treatment` = factor(
            c('EON', 'RL100 (1)', 'RL1000 (1)', 'RL100 (2)', 'RLBL1000', 'RL100 (3)', 'RL1000 (2)'),
            levels = c('EON', 'RL100 (1)', 'RL1000 (1)', 'RL100 (2)', 'RLBL1000', 'RL100 (3)', 'RL1000 (2)')
        ),
        `start` = c(0, 30, 60, 90, 120, 180, 210),
        `end`   = c(30, 60, 90, 120, 180, 210, 240),
        `colors` = c('#d3d3d3', '#ffd6d6', '#ffa6a6', '#ffd6d6', '#cbb2ea', '#ffd6d6', '#ffa6a6')
    )

    main_plot <- ggplot2::ggplot(`data_rlbl`, ggplot2::aes(x = `Time_points`)) +
        ggplot2::stat_summary(
            ggplot2::aes(x = `Time_points`, y = `gsw_normalized`, color = 'gsw'),
            fun = mean,
            geom = 'point',
            alpha = 0.4,
            size = 3
        ) +
        ggplot2::stat_summary(
            ggplot2::aes(x = `Time_points`, y = `gsw_normalized`, color = 'gsw'),
            fun.data = mean_se,
            geom = 'errorbar',
            width = 1
        ) +
        ggplot2::stat_summary(
            ggplot2::aes(x = `Time_points`, y = `A_normalized`, color = 'A'),
            fun = mean,
            geom = 'point',
            alpha = 0.4,
            size = 3
        ) +
        ggplot2::scale_color_manual(
            name = 'Parameter',
            values = c('A' = '#D55E00', 'gsw' = '#0072B2')
        ) +
        ggplot2::ylab('Normalized Value') +
        ggplot2::scale_y_continuous(
            name = paste0(
                "Normalized <span style='color:#D55E00;'>A</span> ",
                "and <span style='color:#0072B2;'>gsw</span>"
            ),
            breaks = seq(-0.2, 1.0, 0.2)
        ) +
        ggplot2::scale_x_continuous(breaks = seq(0, 210, 30)) +
        ggplot2::xlab(expression(paste('Time since EON (min)'))) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            legend.position = 'none',
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_rect(color = 'black', fill = NA, linewidth = 1),
            legend.title = ggplot2::element_text(size = 25),
            legend.text = ggplot2::element_text(size = 23),
            legend.key = ggplot2::element_blank(),
            legend.background = ggplot2::element_blank(),
            legend.box.background = ggplot2::element_blank(),
            legend.key.size = ggplot2::unit(3, 'line'),
            axis.title.x = ggplot2::element_text(size = 25),
            axis.title.y = ggtext::element_markdown(size = 25),
            axis.text.x = ggplot2::element_text(size = 23),
            axis.text.y = ggplot2::element_text(size = 23),
            strip.text = ggplot2::element_text(size = 25),
            strip.background = ggplot2::element_rect(fill = 'grey98')
        ) +
        # Add change values at bin centers
        ggplot2::geom_text(
            data = changes_data,
            mapping = ggplot2::aes(
                x = bin_center,
                y = 1.3,
                label = sprintf("%.2f", A_change)
            ),
            color = "#D55E00",
            size = 8
        ) +
        ggplot2::geom_text(
            data = changes_data,
            mapping = ggplot2::aes(
                x = bin_center,
                y = 1.2,
                label = sprintf("%.2f", gsw_change)
            ),
            color = "#0072B2",
            size = 8
        ) +
        ggplot2::geom_text(
            data = changes_data %>% dplyr::filter(!is.na(starch_change_30min)),
            mapping = ggplot2::aes(
                x = bin_center,
                y = 1.1,
                label = sprintf("%.2f", starch_change_30min)
            ),
            color = "#228B22",  # Forest green color for starch
            size = 8
        ) +
        ggplot2::annotate(
            "text", x = 160, y = 0.1, hjust = 0, vjust = 1,
            label = plot_species, size = 10, color = "black"
        )

    treatment_bar <- ggplot2::ggplot() +
        ggplot2::geom_rect(
            data = treatment_ranges,
            mapping = ggplot2::aes(
                xmin = `start`,
                xmax = `end`,
                ymin = 0,
                ymax = 200,
                fill = `Treatment`
            )
        ) +
        ggplot2::geom_text(
            data = treatment_ranges,
            mapping = ggplot2::aes(
                x = (`start` + `end`) / 2,
                y = 100,
                label = stringr::str_wrap(
                    stringr::str_remove(`Treatment`, '\\(.*\\)'),
                    width = 4
                )
            ),
            color = 'black',
            size = 4.5,
            angle = 0,
            hjust = 0.5
        ) +
        ggplot2::scale_fill_manual(values = treatment_ranges$colors) +
        ggplot2::scale_y_continuous(name = '30 min', limits = c(0, 200)) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            legend.position = 'none'
        )

    patchwork::wrap_plots(
        treatment_bar, main_plot,
        ncol = 1,
        heights = c(0.05, 0.95)
    )
}
# nolint end