## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(ggplot2)
library(grid)
library(ggvolcano)


## ----example-1, fig.width=7.2, fig.height=5-----------------------------------
data <- read.csv(system.file("extdata", "example.csv", package = "ggvolcano"))

ggvolcano(data,
          logFC_col = "log2FoldChange",
          pval_col = "pvalue",
          pval_cutoff = 1e-6,
          logFC_cutoff = 1.0,
          title = "Basic Volcano Plot",
          caption = paste("Total variables:", nrow(data)))

## ----example-2, fig.width=7.2, fig.height=5-----------------------------------
ggvolcano(data,
          logFC_col = "log2FoldChange",
          pval_col = "pvalue",
          pval_cutoff = 1e-4,
          logFC_cutoff = 1.5,
          title = "Custom Points Volcano Plot",
          caption = "Customized point shapes and colors",
          point_aes = list(
            size = 2,
            shape = c(16, 17, 18, 19),  # Different shapes for NS, FC, P, and FC_P
            color = c("lightgrey", "blue", "green", "red"),
            alpha = 0.8
          ))


## ----example-3, fig.width=7.2, fig.height=5-----------------------------------
ggvolcano(data,
          logFC_col = "log2FoldChange",
          pval_col = "pvalue",
          pval_cutoff = 1e-6,
          logFC_cutoff = 1.0,
          title = "Volcano Plot with Color Gradient",
          caption = "Color gradient represents the range of p-values",
          point_aes = list(
            size = 2,
            shape = c(16, 17, 18, 19),
            alpha = 0.85,
            # Activate the color gradient by specifying two colors and associated options:
            color_gradient = c("yellow", "red"),
            color_gradient_limits = c(min(data$pvalue, na.rm = TRUE), max(data$pvalue, na.rm = TRUE)),
            color_gradient_breaks = seq(min(data$pvalue, na.rm = TRUE), max(data$pvalue, na.rm = TRUE), length.out = 5),
            color_gradient_labels = round(seq(min(data$pvalue, na.rm = TRUE), max(data$pvalue, na.rm = TRUE), length.out = 5), 3)
          ))



## ----example-4, fig.width=7.2, fig.height=5-----------------------------------
ggvolcano(data,
          labels = data$X,            # Use the column 'X' for labels
          logFC_col = "log2FoldChange",
          pval_col = "pvalue",
          pval_cutoff = 1e-6,
          logFC_cutoff = 2.0,
          xlab = "Log2 Fold Change",   # Custom x-axis label
          ylab = "-Log10(p-value)",    # Custom y-axis label
          title = "Enhanced Volcano Plot with Labels and Connectors",
          subtitle = "Differential Expression Analysis",
          caption = "p-value gradient: yellow (low) to red (high)",
          point_aes = list(
            size = 1.5,
            shape = c(16, 17, 18, 19),
            alpha = 0.9,
            color_gradient = c("blue", "red"),
            color_gradient_limits = c(min(data$pvalue, na.rm = TRUE),
                                      max(data$pvalue, na.rm = TRUE)),
            color_gradient_breaks = seq(min(data$pvalue, na.rm = TRUE),
                                        max(data$pvalue, na.rm = TRUE),
                                        length.out = 5),
            color_gradient_labels = round(seq(min(data$pvalue, na.rm = TRUE),
                                              max(data$pvalue, na.rm = TRUE),
                                              length.out = 5), 3)
          ))



## ----example-5, fig.width=7.2, fig.height=5-----------------------------------
ggvolcano(data,
          logFC_col = "log2FoldChange",
          pval_col = "pvalue",
          pval_cutoff = 1e-4,
          logFC_cutoff = 1.5,
          title = "Volcano Plot with Custom Gridlines & Border",
          caption = "Major gridlines only; full border in blue",
          gridlines = list(major = TRUE, minor = FALSE),
          plot_border = "full",
          border_width = 1.0,
          border_color = "blue")



## ----example-6, fig.width=7.2, fig.height=5-----------------------------------

ggvolcano(data,
          logFC_col = "log2FoldChange",
          pval_col = "pvalue",
          pval_cutoff = 0.05,
          logFC_cutoff = 1.0,
          title = "Volcano Plot with Reference Line",
          caption = "Extra horizontal line",
          horizontal_line = 1e-9,
          horizontal_line_aes = list(type = "dotted", color = "purple", width = 1))


## ----example-7, fig.width=7.2, fig.height=5-----------------------------------
ggvolcano(data,
          logFC_col = "log2FoldChange",
          pval_col = "pvalue",
          pval_cutoff = 1e-6,
          logFC_cutoff = 1.0,
          xlab = "Log2 Fold Change",
          ylab = "-Log10(p-value)",
          title = "Volcano Plot with Custom Axes",
          subtitle = "Differential Expression Analysis",
          caption = "Customized axis labels and subtitle")



## ----example-8, fig.width=7.2, fig.height=6-----------------------------------

ggvolcano(data,
          logFC_col = "log2FoldChange",
          pval_col = "pvalue",
          pval_cutoff = 1e-4,
          logFC_cutoff = 1.5,
          title = "Volcano Plot with Custom Legend",
          caption = "Legend labels and position modified",
          legend_aes = list(
            labels = c("Non-sig", "FC Sig", "P-value Sig", "Both Sig"),
            position = "bottom",
            label_size = 12,
            icon_size = 6
          ))

