# 
# Plot percentage of modified proteoforms matched with the proteoforms of proteins with multiple proteoforms.
# Each proteoform is altered and matched against proteoforms in the database that share the same accession.
# From all possible proteoforms, we calculate the percentage of them that got matched using each matching.
#
start_time <- Sys.time()

# Libraries

library(ggplot2)
library(cowplot)
library(plyr)
source("R/plotPercentages.R")

# Parameters
percentagesFile <- "resources/sensitivity/percentagesFileMultiproteoforms.csv"

# Main script

print(paste("Loading data from", percentagesFile))
percentages <- read.csv(percentagesFile, header = T)

print(paste("Plotting percentages"))
percentagesPlot <- plotPercentagesSeparated(percentages)
percentagesPlot

png("docs/fig7.png", height = 12, width = 12, units = "cm", res = 600)
plot(percentagesPlot)
dummy <- dev.off()

print("Finished")
print(Sys.time() - start_time)
