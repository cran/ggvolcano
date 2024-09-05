
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ggvolcano

<!-- badges: start -->
<!-- badges: end -->

The goal of **ggvolcano** is to provide a flexible and customizable
solution for generating publication-ready volcano plots in R. It
simplifies the process of visualizing differential expression results
from analyses like RNA-seq, making it easier to communicate key
findings.

## Installation

You can install the development version of **ggvolcano** from GitHub
using the following commands:

``` r
# Install devtools if you haven't already
install.packages("devtools")
```

## Example

This is a basic example showing how to create a volcano plot using the
ggvolcano package:

``` r
# Load the package
library(ggvolcano)

# Load example data
data <- read.csv(system.file("extdata", "example.csv", package = "ggvolcano"))

# Create a volcano plot
ggvolcano(data,
  logFC_col = "log2FoldChange", # Column containing log2 fold changes
  pval_col = "pvalue",          # Column containing p-values
  pval_cutoff = 10e-4,          # Cutoff for significance in p-values
  logFC_cutoff = 1.5,           # Cutoff for significance in fold change
  x_limits = c(-5.5, 5.5),      # X-axis limits
  y_limits = c(0, -log10(10e-12)), # Y-axis limits
  title = "Example Volcano plot",   # Title of the plot
  caption = "FC cutoff, 1.5; p-value cutoff, 10e-4", # Caption of the plot
  legend_aes = list(position = "right"),  # Position legend on the right
  gridlines = list(major = TRUE, minor = TRUE), # Show major and minor gridlines
  horizontal_line = 10e-8, # Draw a horizontal line at this p-value cutoff
  horizontal_line_aes = list(type = "dashed", color = "red", width = 0.5)  # Aesthetics for horizontal line
)
```
