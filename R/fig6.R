# 
# Plot two individual proteoform matching percentages. 
# Each proteoform is alterd and matched against proteoforms in the database that share the same accession.
# From all possible proteoforms, we calculate the percentage of them that got matched using each matching.
#
start_time <- Sys.time()

# Libraries

library(ggplot2)
library(cowplot)
library(plyr)
source("R/plotPercentages.R")

# Parameters
percentagesFile1 <- "resources/sensitivity/matchesFile1.csv"
percentagesFile2 <- "resources/sensitivity/matchesFile2.csv"

# Main script

print(paste("Loading data from", percentagesFile1))
percentages1 <- read.csv(percentagesFile1, header = T)
print(paste("Loading data from", percentagesFile2))
percentages2 <- read.csv(percentagesFile2, header = T)


print(paste("Plotting percentages"))
plot1 <- plotPercentagesSeparated(percentages1)
plot1
plot2 <- plotPercentagesSeparated(percentages2)
plot2

print("Making grid")
grid <- plot_grid(plot1 + theme(legend.position = "none"), 
                  plot2 + theme(legend.position = "none"),
                  labels = c("A", "B"),
                  align = "vh",
                  hjust = -1, nrow = 1)
legend <- get_legend(plot1)
grid <- plot_grid(grid, legend, rel_widths = c(2, 0.3))
grid
print("Saving plot")
save_plot("docs/fig6.png", grid,
          ncol = 2, nrow =1,
          base_aspect_ratio = 2:1
)

print("Finished")
print(Sys.time() - start_time)

