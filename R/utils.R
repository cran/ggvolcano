#' Publication-ready volcano plots
#'
#' Creates a volcano plot to visualize differential expression results.
#' This function is highly configurable to suit publication standards.
#'
#' @param data A data frame containing test statistics. Requires at least columns for variable names, log2 fold changes, and p-values.
#' @param labels Column name or row names for variable names.
#' @param logFC_col Column name for log2 fold changes.
#' @param pval_col Column name for nominal or adjusted p-values.
#' @param x_limits Limits of the x-axis (default auto-calculated).
#' @param y_limits Limits of the y-axis (default auto-calculated).
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param caption Plot caption.
#' @param pval_cutoff P-value cutoff for significance.
#' @param logFC_cutoff Log2 fold-change cutoff for significance.
#' @param cutoff_line List of options for cutoff lines (`type`, `color`, `width`).
#' @param point_aes List of aesthetic options for points (`size`, `shape`, `color`, `alpha`).
#' @param label_aes List of aesthetic options for labels (`size`, `color`, `face`, `parse`).
#' @param legend_aes List of aesthetic options for legend (`labels`, `position`, `label_size`, `icon_size`).
#' @param shade_options List of options for shading regions in the plot.
#' @param connector_aes List of aesthetic options for connectors (`line_width`, `arrow_type`, `arrow_ends`, `arrow_length`, `line_color`, `direction`, `draw_arrowheads`).
#' @param gridlines List with logical values indicating whether to draw gridlines (`major`, `minor`).
#' @param plot_border Add a border for plot axes (`"partial"` or `"full"`).
#' @param border_width Width of the border.
#' @param border_color Color of the border.
#' @param horizontal_line Numeric value(s) for drawing horizontal line(s).
#' @param horizontal_line_aes List of aesthetic options for the horizontal line(s) (`type`, `color`, `width`).
#'
#' @return A \code{ggplot2} object representing the volcano plot.
#'
#' @examples
#' data <- read.csv(system.file("extdata", "example.csv", package = "ggvolcano"))
#'
#' ggvolcano(data,
#'   logFC_col = "log2FoldChange",
#'   pval_col = "pvalue",
#'   pval_cutoff = 10e-4,
#'   logFC_cutoff = 1.5,
#'   x_limits = c(-5.5, 5.5),
#'   y_limits = c(0, -log10(10e-12)),
#'   title = "Example Volcano plot",
#'   caption = "FC cutoff, 1.5; p-value cutoff, 10e-4",
#'   gridlines = list(major = TRUE, minor = TRUE),
#'   horizontal_line = 10e-8, # Draw horizontal line for p-value cutoff
#'   horizontal_line_aes = list(type = "dashed", color = "red", width = 0.5)
#' )
#'
#' @import ggplot2
#' @import ggrepel
#' @export
ggvolcano <- function(data,
                      labels = "",
                      logFC_col,
                      pval_col,
                      x_limits = c(min(data[[logFC_col]], na.rm = TRUE) - 1.5, max(data[[logFC_col]], na.rm = TRUE) + 1.5),
                      y_limits = c(0, max(-log10(data[[pval_col]]), na.rm = TRUE) + 5),
                      xlab = bquote( ~ Log[2] ~ "fold change"),
                      ylab = bquote( ~ -Log[10] ~ italic(P)),
                      title = "Volcano plot",
                      subtitle = "",
                      caption = paste0("total = ", nrow(data), " variables"),
                      pval_cutoff = 1e-6,
                      logFC_cutoff = 1.0,
                      cutoff_line = list(type = "longdash",
                                         color = "black",
                                         width = 0.4),
                      point_aes = list(
                        size = 1.5,
                        shape = c(19, 19, 19, 19),
                        color = c("grey30", "#00CD6C", "#009ADE", "#FF1F5B"),
                        alpha = 0.9
                      ),
                      label_aes = list(
                        size = 2.5,
                        color = "black",
                        face = "plain",
                        parse = FALSE
                      ),
                      legend_aes = list(
                        labels = c(
                          "NS",
                          expression(Log[2] ~ FC),
                          "p-value",
                          expression(p - value ~ and ~ log[2] ~ FC)
                        ),
                        position = "right",
                        label_size = 14,
                        icon_size = 5
                      ),
                      shade_options = NULL,
                      connector_aes = list(
                        line_width = 0.5,
                        arrow_type = "closed",
                        arrow_ends = "first",
                        arrow_length = unit(0.01, "npc"),
                        line_color = "grey10",
                        direction = "both",
                        draw_arrowheads = TRUE
                      ),
                      gridlines = list(major = TRUE, minor = TRUE),
                      plot_border = "partial",
                      border_width = 0.8,
                      border_color = "black",
                      horizontal_line = NULL,
                      horizontal_line_aes = list(type = "longdash",
                                                 color = "black",
                                                 width = 0.4))
{

  # Check if the fold change column is numeric
  if (!is.numeric(data[[logFC_col]])) {
    stop(paste(logFC_col, " is not numeric!", sep = ""))
  }

  # Check if the p-value column is numeric
  if (!is.numeric(data[[pval_col]])) {
    stop(paste(pval_col, " is not numeric!", sep = ""))
  }

  # Prepare data table by converting to a data frame
  data <- as.data.frame(data)

  # Define significance groups based on p-value and logFC cutoffs
  data$significance <- "NS" # Not significant by default
  data$significance[(abs(data[[logFC_col]]) > logFC_cutoff)] <- "FC" # Significant fold change
  data$significance[(data[[pval_col]] < pval_cutoff)] <- "P" # Significant p-value
  data$significance[(data[[pval_col]] < pval_cutoff) & (abs(data[[logFC_col]]) > logFC_cutoff)] <- "FC_P" # Significant in both
  data$significance <- factor(data$significance, levels = c("NS", "FC", "P", "FC_P"))

  # Handle cases where p-value is zero
  if (min(data[[pval_col]], na.rm = TRUE) == 0) {
    warning(
      paste(
        "One or more p-values is 0.",
        "Converting to 10^-1 * current",
        "lowest non-zero p-value..."
      ),
      call. = FALSE
    )
    data[which(data[[pval_col]] == 0), pval_col] <- min(data[which(data[[pval_col]] != 0), pval_col], na.rm = TRUE) * 10^-1
  }

  data$labels <- labels
  data$xvals <- data[[logFC_col]]
  data$yvals <- data[[pval_col]]

  base_theme <- theme_bw(base_size = 24) +
    theme(
      legend.background = element_rect(),
      plot.title = element_text(angle = 0, size = 18, face = "bold", vjust = 1),
      plot.subtitle = element_text(angle = 0, size = 14, face = "plain", vjust = 1),
      plot.caption = element_text(angle = 0, size = 14, face = "plain", vjust = 1),
      axis.text.x = element_text(angle = 0, size = 18, vjust = 1),
      axis.text.y = element_text(angle = 0, size = 18, vjust = 0.5),
      axis.title = element_text(size = 18),
      legend.position = legend_aes$position,
      legend.key = element_blank(),
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_text(size = legend_aes$label_size),
      title = element_text(size = legend_aes$label_size),
      legend.title = element_blank()
    )

  # Plot setup depending on color gradient option
  if (is.null(point_aes$color_gradient)) {
    plot <- ggplot(data, aes(x = xvals, y = -log10(yvals))) +
      base_theme +
      guides(colour = guide_legend(
        order = 1,
        override.aes = list(
          shape = c(NS = point_aes$shape[1], FC = point_aes$shape[2], P = point_aes$shape[3], FC_P = point_aes$shape[4]),
          size = legend_aes$icon_size
        )
      )) +
      geom_point(aes(color = significance, shape = significance), alpha = point_aes$alpha, size = point_aes$size, na.rm = TRUE) +
      scale_color_manual(
        values = c(NS = point_aes$color[1], FC = point_aes$color[2], P = point_aes$color[3], FC_P = point_aes$color[4]),
        labels = legend_aes$labels,
        drop = TRUE
      ) +
      scale_shape_manual(
        values = c(NS = point_aes$shape[1], FC = point_aes$shape[2], P = point_aes$shape[3], FC_P = point_aes$shape[4]),
        guide = "none",
        drop = TRUE
      )
  } else {
    # If using a color gradient
    plot <- ggplot(data, aes(x = xvals, y = -log10(yvals))) +
      base_theme +
      geom_point(aes(color = yvals, shape = significance), alpha = point_aes$alpha, size = point_aes$size, na.rm = TRUE) +
      scale_colour_gradient(
        low = point_aes$color_gradient[1],
        high = point_aes$color_gradient[2],
        limits = point_aes$color_gradient_limits,
        breaks = point_aes$color_gradient_breaks,
        labels = point_aes$color_gradient_labels
      ) +
      scale_shape_manual(
        values = c(NS = point_aes$shape[1], FC = point_aes$shape[2], P = point_aes$shape[3], FC_P = point_aes$shape[4]),
        guide = "none",
        drop = TRUE
      )
  }

  # Add cutoff lines
  plot <- plot +
    xlab(xlab) +
    ylab(ylab) +
    xlim(x_limits[1], x_limits[2]) +
    ylim(y_limits[1], y_limits[2]) +
    geom_vline(
      xintercept = c(-logFC_cutoff, logFC_cutoff),
      linetype = cutoff_line$type,
      colour = cutoff_line$color,
      linewidth = cutoff_line$width
    ) +
    geom_hline(
      yintercept = -log10(pval_cutoff),
      linetype = cutoff_line$type,
      colour = cutoff_line$color,
      linewidth = cutoff_line$width
    )

  # Add title, subtitle, caption
  plot <- plot + labs(
    title = title,
    subtitle = subtitle,
    caption = caption
  )

  # Add vertical and horizontal lines if specified
  if (!is.null(horizontal_line)) {
    plot <- plot + geom_hline(
      yintercept = -log10(horizontal_line),
      linetype = horizontal_line_aes$type,
      colour = horizontal_line_aes$color,
      linewidth = horizontal_line_aes$width
    )
  }

  # Add plot borders
  if (plot_border == "full") {
    plot <- plot + theme(panel.border = element_rect(colour = border_color, fill = NA, size = border_width))
  } else if (plot_border == "partial") {
    plot <- plot + theme(
      axis.line = element_line(linewidth = border_width, colour = border_color),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
  }

  # Add gridlines
  if (gridlines$major) {
    plot <- plot + theme(panel.grid.major = element_line())
  } else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines$minor) {
    plot <- plot + theme(panel.grid.minor = element_line())
  } else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }

  # Handle text labels
  plot <- plot + geom_text_repel(
    data = subset(data, data[[pval_col]] < pval_cutoff & abs(data[[logFC_col]]) > logFC_cutoff),
    aes(label = subset(data, data[[pval_col]] < pval_cutoff & abs(data[[logFC_col]]) > logFC_cutoff)[["labels"]]),
    size = label_aes$size,
    segment.color = connector_aes$line_color,
    segment.size = connector_aes$line_width,
    arrow = if (connector_aes$draw_arrowheads) {
      arrow(
        length = connector_aes$arrow_length,
        type = connector_aes$arrow_type,
        ends = connector_aes$arrow_ends
      )
    } else {
      NULL
    },
    colour = label_aes$color,
    fontface = label_aes$face,
    parse = label_aes$parse,
    na.rm = TRUE,
    direction = connector_aes$direction,
    max.overlaps = 15,
    min.segment.length = 0
  )

  # Optional shading if specified
  if (!is.null(shade_options)) {
    plot <- plot +
      stat_density2d(
        data = subset(data, rownames(data) %in% shade_options$variables),
        fill = shade_options$fill,
        alpha = shade_options$alpha,
        geom = "polygon",
        contour = TRUE,
        size = shade_options$size,
        bins = shade_options$bins,
        show.legend = FALSE,
        na.rm = TRUE
      )
  }

  # Ensure the plot does not clip any labels
  plot <- plot + coord_cartesian(clip = "off")

  # Return the plot object
  return(plot)
}


utils::globalVariables(c("xvals"))
utils::globalVariables(c("yvals"))
utils::globalVariables(c("significance"))
