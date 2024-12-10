#Plot without real anomaly indices. 
library(ggplot2)

plot_density_ratios <- function(density_ratios, threshold, is_outlier, title) {
  
  # Create a data frame with a logical column for real anomaly indices
  df <- data.frame(
    Index = seq_along(density_ratios),
    DensityRatio = density_ratios,
    Outlier = is_outlier
  )
  
  # Plot with ggplot2
  ggplot(df, aes(x = Index, y = DensityRatio)) +
    geom_point(aes(color = Outlier)) +  # Color points based on 'Outlier'
    geom_hline(yintercept = threshold, color = "red", linetype = "dashed") +
    labs(title = title, x = "Test Data Index", y = "Density Ratio") +  # Dynamically set the title
    scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "blue"),
      labels = c("TRUE" = "Outlier", "FALSE" = "Inlier")
    ) +
    theme_minimal()
}
